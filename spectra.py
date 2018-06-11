from utils import *
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
from astropy.convolution import Gaussian1DKernel
from scipy.signal import convolve

def convolve_spectrum(x,dx,y,R):
	# Calculate the median line width in this range in the wavelength units
	sigma = np.median(x)/R
	
	# Divide by the bin size of x to get the dimensionless sigma (i.e. in units of x bins)
	sigma /= dx*x.unit
	
	# Calculate the 1D Gaussian kernel (in dimensionless x) and convolve y
	kernel = Gaussian1DKernel(sigma.decompose().value)
	yc = convolve(y.value,kernel,mode='same') * y.unit
		
	return yc

def generate_planck_blackbody(x,T_eff,frequency=None):
	# Generates a Planck SED
	ptype,unit = dimensions_check(T_eff,accepted_dims=[u'temperature'])
	ptype,unit = dimensions_check(x,accepted_dims=[u'length',u'frequency'])
	if ptype is None or unit is None:
		return None
	
	if (frequency is None and ptype==u'frequency') or frequency:
		y = blackbody_nu(x,T_eff)
	else:
		y = blackbody_lambda(x,T_eff)
		
	return y

def read_PHOENIX_spectrum(filename,x):
	# Read in a PHOENIX atmosphere spectrum

	# Interpolate the photosphere model onto x
	ym = interpolate_model(models_dir+'/star/'+filename,x)
	
	# Add units (for PHOENIX models, it's ergs/s/cm2/cm)
	ym *= u.erg/u.cm**2/u.s/u.sr/u.cm
	
	return ym

def SED_to_observed(y,dist,R_st):
	# Converts an SED (e.g. Planck's law or PHOENIX model) to an observed spectrum
	# taking into account the size of and distance to the star
	
	ptype,unit = dimensions_check(dist,accepted_dims=[u'length'])
	ptype,unit = dimensions_check(R_st,accepted_dims=[u'length'])
	
	y *= (R_st**2/dist**2)*u.sr
	
	return y.to(u.erg/u.s/u.cm**2/u.nm)
