from utils import *
import batman


def get_atmosphere_model(filename,x,R_st):
	# Loads an atmospheric model from file and interpolates onto x
	# Then modulates y to match the stellar radius

	# Interpolate the atmosphere model onto x
	ym = interpolate_model(models_dir+'/planet/'+filename,x)
	
	# Divide the ym values by (R_st/R_sol)
	ym /= (R_st/ct.R_sun).decompose().value
	
	return ym
	
def get_contamination_model(filename,x):
	# Loads a contamination model (wavelength vs. (Rp/Rp0))
	if filename is None:
		return np.ones_like(x,dtype=float)

def base_transit_params(M_st=1.0*ct.M_sun,R_st=1.0*ct.R_sun,RpRs=None,t0=0.,per=None,aRs=None,inc=90.,ecc=0.,w=90.,uld=[0.1,0.3],limb_dark="quadratic"):
	# Calculate the period from the semi-major axis, or semi-major axis from period
	if per is None and aRs is None:
		print("Must specify either period or a/Rs!")
		return None
	if per is None:
		per = ((4*np.pi**2/(ct.G*M_st))*(aRs*R_st)**3)**0.5
	if aRs is None:
		aRs = (((ct.G*M_st)/(4*np.pi**2))*(per)**2)**(1./3.)/R_st
		aRs = aRs.decompose().value
		
	# Convert the period to days and strip the units
	per = per.to(u.day).value
	
	# Convert t0 to days
	t0 = t0.to(u.day).value
	
	params = batman.TransitParams()
	params.t0 = t0                      #time of inferior conjunction
	params.per = per                    #orbital period
	params.rp = RpRs                    #planet radius (in units of stellar radii)
	params.a = aRs                      #semi-major axis (in units of stellar radii)
	params.inc = inc                    #orbital inclination (in degrees)
	params.ecc = ecc                    #eccentricity
	params.w = w                       	#longitude of periastron (in degrees)
	params.u = uld              		#limb darkening coefficients
	params.limb_dark = limb_dark      	#limb darkening model
	return params

def generate_transit_models(RpRs,params,t):
	# Convert t to days and strip the unit
	t = t.to(u.day).value
	
	# Loop through the RpRs values and calculate the transit model
	RpRs = np.array(RpRs)
	model_LC = np.zeros((len(t),len(RpRs)),dtype=float)
	for i in range(len(RpRs)):
		params.rp = RpRs[i]
		mod = batman.TransitModel(params,t)
		model_LC[:,i] = mod.light_curve(params)

	return model_LC








