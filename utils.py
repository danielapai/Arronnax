# Shared modules
from astropy import constants as ct
from astropy import units as u
from matplotlib import pyplot as plt
import numpy as np
import tqdm

# Directory containing all models
models_dir = 'models/'

def dimensions_check(x,accepted_dims=None):
	# Check that the array has units
	if type(x) is not u.quantity.Quantity:
		print("Error: array must have units!")
		return None
	
	# Get the physical type and dimensions
	ptype = x.unit.physical_type
	unit = x.unit
	
	# Check that the dimensions are consistent with 'accepted' (optional)
	if accepted_dims is not None:
		# Ensure that this is a list
		accepted_dims = np.array(accepted_dims)
		ptype = x.unit.physical_type
		if ptype not in accepted_dims:
			print("Error: array must have dimensions in: {0}".format(accepted_dims))
			return None
	
	return unit,ptype

def interpolate_model(model_file,x):
	# Interpolates a model onto wavelength grid x
	# The first line in the model file should be the unit, otherwise assumed to be the same as x
	
	# Load the numerical values
	xm,ym = np.loadtxt(model_file,unpack='True')
	
	# Sort by x
	idx = np.argsort(xm)
	xm, ym = xm[idx], ym[idx]
	
	# The unit should be listed on the first line e.g., '# nm'
	with open(model_file,'r') as f:
		lines = f.readlines()
		if '#' in lines[0]:
			unit = lines[0].strip('#').split()[0]
		else:
			unit = x.unit
	xm = xm*u.Unit(unit)
	
	# Convert xm to the units of x
	xm = xm.to(x.unit)

	# Interpolate the y values onto the x grid
	ym = np.interp(x,xm,ym,left=0.,right=0.)
	
	return ym

def read_params(params_file):
	# Read the params file into a dictionary
	with open(params_file,'r') as f:
		lines = f.readlines()
		
	# Loop through the lines; look for ! flags to start a new level in the dictionary
	flag = None
	p = {}
	for line in lines:
		if '!' in line:
			flag = line.replace('(','').replace(')','').split('!')[1].strip()
			p[flag] = {}
			continue
		elif flag is None or line.strip() == '' or line.strip()[0] == '#':
			continue
			
		key,unit,val = line.split('#')[0].strip().split()
		
		# Determine if val is a float, int, string, or None
		if val.lower()=='none':
			p[flag][key] = None
		elif '_' in val: # Odd exception for some numeric strings with underscores in them
			p[flag][key] = val
		else:	
			try:
				p[flag][key] = int(val) # Is it an integer?
			except ValueError:
				try:
					p[flag][key] = float(val) # Is it a float?
				except ValueError: 
					p[flag][key] = val # Otherwise it's a string
		
		# If a unit is specified, apply the unit
		if unit.lower() not in ['-','none']:
			p[flag][key] = p[flag][key]*u.Unit(unit)
		
	return p

def rebin_array(x,factor,axis):
	# Strip the units for now
	if type(x) is u.quantity.Quantity:
		x = x.decompose()
		unit = x.unit
		x = x.value
	else:
		unit = 1

	# Rebins array x along axis based on the integer factor; discards any remainder bins
	shape = np.array(x.shape)
	shape[axis] = int(shape[axis]/factor)
	
	# Loop through the new bins
	left = 0
	x_out = []
	for i in range(shape[axis]):
		# Indices along axis corresponding to this bin
		idxes = np.arange(left,left+factor,1,dtype=int)
		
		# Retrieve an array of values in this bin
		bin_vals = np.take(x,idxes,axis=axis)
		
		# Sum along the axis and save to the rebinned array
		bin_sum = np.sum(bin_vals,axis=axis)
		x_out.append(bin_sum)
		
		left += factor
	
	# Resize x_out to the proper shape
	x_out = np.array(x_out)
	x_out = np.resize(x_out,shape)
	
	return x_out*unit
