from utils import *
from spectra import *
from transit import *

# Wavelength range parameters
xmin,xmax,dx = 300,3000,1.0
wl_unit = u.nm
t_unit = u.day

# Initialize

# Read in the parameters file
p = read_params('params')

# Unpack and apply units
planet_file = p['planet_name']+'.planet'
star_file = None if p['star_name'] is None else p['star_name']+'.star'
det_file = None if p['det_name'] is None else p['det_name']+'.detector'
a = p['a']*u.AU
R_st = p['R_st']*u.Rsun
M_st = p['M_st']*u.Msun
T_eff = p['T_eff']*u.K
d = p['d']*u.pc 
D = p['D']*u.m
N_units = p['N_units']
px_scale = p['px_scale']*u.arcsecond/u.pixel
gain = p['gain']*u.electron/u.photon
read_noise = p['read_noise']*u.electron
throughput = p['throughput']
R = p['R']
N_exp = p['N_exp']
t_exp = p['t_exp']*u.s
binning = p['binning']*u.nm

# Generate the wavelength and time grids
x = np.arange(xmin,xmax,dx,dtype=float)*wl_unit
te = t_exp.to(t_unit).value
t = np.arange(0,N_exp*te,te)*t_unit

# If a detector model is specified, interpolate onto the wavelength grid
if det_file is not None:
	throughput = interpolate_model(models_dir+'/detector/'+det_file,x)
	# Remove all wavelengths where throughput < 1%
	mask = throughput > 0.01
	x = x[mask]
	throughput = throughput[mask]

# Part 1: Compute the apparent photon flux of the source vs. wavelength and time in each "micro" bin, in units of phot/s/cm^2
def compute_apparent_flux(x,t):
	# STEP 1: Generate a spectrum for the target (erg/s/cm^2/sr/nm)
	if star_file is None:
		y = generate_planck_blackbody(x,T_eff)
	else:
		y = read_PHOENIX_spectrum(star_file,x)
	
	# STEP 2: Read in the atmosphere+stellar contamination models (RpRs vs wavelength)
	RpRs_m = get_atmosphere_model(planet_file,x,R_st)
	eps = get_contamination_model(None,x)
	RpRs_m *= eps.value

	# STEP 3: Build the normalized transit light curves for each wavelength; multiply the stellar spectrum by this model
	params = base_transit_params(t0=np.ptp(t)/2.,M_st=M_st,R_st=R_st,per=365.25*u.day)
	LC_m = generate_transit_models(RpRs_m,params,t)
	y = y[None,:]*LC_m
	
	# STEP 4: Integrate over the solid angle of the star and calculate the photon flux, then the photon flux in each bin
	y *= (np.pi*R_st**2/d**2)*u.sr
	e_phot = (ct.h*ct.c/x).to(u.erg)/u.photon
	y = (y[None,:]/e_phot)[0]
	y *= dx*wl_unit

	return y.to(u.photon/u.s/u.cm**2),LC_m
	
y_physical,LC_m = compute_apparent_flux(x,t)

# Part 2: Compute the ADU spectrum recorded by a single telescope unit in a single exposure, in units of e-
def expose_single_unit(t,x,y_in):
	# STEP 0: Make a copy of y_in so as to avoid affecting the original
	y = np.copy(y_in)*y_in.unit

	# STEP 1: Multiply the incident flux by the objective area and detector QE
	y *= np.pi*(D/2.)**2
	y *= throughput
	
	# STEP 2: Multiply by exposure time
	y *= t_exp
	
	# STEP 3: Apply photon noise (carry over the units)
	y = np.random.normal(y,scale=y**0.5)*y.unit

	# STEP 4: Convolve the spectrum by the spectral resolution for each time index
	for i in range(len(t)):
		y[i,:] = convolve_spectrum(x,dx,y[i,:],R)

	# STEP 5: Convert photons to electrons via the gain
	y *= gain

	return y.to(u.electron)

y_single = expose_single_unit(t,x,y_physical)

"""
# STEP 5: Calculate the pixel size in wavelength units and apply read noise to each pixel
def apply_read_noise(y_in):
	y = np.copy(y_in)
	fwhm = (np.median(x)/D).decompose()*206265*u.arcsecond/px_scale # FWHM in pixels
	px = (x-x[0])/(dx*u.nm)*fwhm
	px = np.array(px.value,dtype=int)
	unique_px = np.unique(px)
	for i in range(len(unique_px)):
		mask = px == unique_px[i]
		px_sum = np.sum(y[:,mask],axis=1)
		rn_delta = np.random.normal(px_sum,scale=read_noise/gain)-px_sum
		y[:,mask] += rn_delta[:,None]/np.sum(mask)
	return y

# STEP 6: Combine the results of each unit, calculating the read noise each time
y_combined = np.zeros_like(y,dtype=float)
for i in range(N_units):
	y_combined += apply_read_noise(y)

summed = np.sum(y_combined,axis=1)
plt.plot(t,summed/np.median(summed[:50])*1e6)
plt.show()
"""

# Part 3: Get the binned light curves, summed from all telescope units
def combined_binned_lightcurve(t,x,y_in,y_mod):
	# STEP 1: Calculate the light curve for each unit, then sum the results
	y_obs = np.zeros((len(t),len(x),N_units),dtype=float)*u.electron
	bar = tqdm.tqdm
	for i in bar(range(N_units)):
		y_obs[:,:,i] = expose_single_unit(t,x,y_in)
	y_sum = np.sum(y_obs,axis=2)

	# STEP 2: Bin down to match the desired binning
	bin_factor = int(binning/(dx*x.unit))
	y_binned = np.zeros((len(t),int(len(x)/bin_factor)),dtype=float)
	y_mod_binned = np.zeros((len(t),int(len(x)/bin_factor)),dtype=float)
	for i in range(len(t)):
		# axis = 1 mode isn't working...
		y_binned[i,:] = rebin_array(y_sum[i,:],bin_factor,axis=0)
		y_mod_binned[i,:] = rebin_array(y_mod[i,:],bin_factor,axis=0)
	
	# and get the average value of x in each bin
	x_binned = rebin_array(x,bin_factor,axis=0)/bin_factor
	x_binned = x_binned.to(wl_unit)

	return x_binned,y_binned,y_mod_binned

xb,yb,ymb = combined_binned_lightcurve(t,x,y_physical,LC_m)

# Part 4: Retrieve RpRs for each wavelength bin
def retrieve_transmission_spectrum_simple(t,x,y,ymod):
	# STEP 1: Mask out the points in the baseline
	mask = np.ones((len(t),len(x)),dtype=bool)
	mask[ymod==np.amax(ymod)] = False
	
	# STEP 2: Use the masked baseline points to scale the model to the data
	baseline = np.zeros(len(x),dtype=float)
	for i in range(len(x)):
		baseline[i] = np.median(y[mask[:,i],i])/np.median(ymod[mask[:,i],i])
	ymod = ymod*baseline[None,:]
	
	# STEP 3: Compute the residuals, the ppm scatter, and the ppm scatter for the entire transit width
	res = y-ymod
	sigma_ppm = np.std(res,axis=0)/np.median(y,axis=0)*1e6
	sigma_ppm_RpRs2 = sigma_ppm/np.sum(mask,axis=0)**0.5

	# STEP 4: Randomly re-draw RpRs based on this sigma value
	RpRs = get_atmosphere_model(planet_file,x,R_st)
	dRpRs = (sigma_ppm_RpRs2/1e6)/(2*RpRs**0.5)
	RpRs = np.random.normal(RpRs,scale=dRpRs)
	
	return RpRs,dRpRs
	
RpRs,dRpRs = retrieve_transmission_spectrum_simple(t,xb,yb,ymb)

# Plot the transmission spectrum
RpRs_fine = get_atmosphere_model(planet_file,x,R_st)

plt.errorbar(xb.value,RpRs,xerr=np.ptp(xb.value)/len(xb)/2,yerr=dRpRs,linestyle='None',marker='s')
plt.plot(x.value,RpRs_fine,c='green')
plt.show()












