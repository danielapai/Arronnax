from matplotlib import pyplot as plt
import numpy as np

# Simulate yields for a several years long mission based on Ehrenreich et al. (2006) scaling relations

# PARAMETERS
# "Effective" aperture size of telescope (in m)
eps_D = 0.9 * 50

# Maximum duration of mission (in days)
t_dur = 365.25*10

# Density of G,K,M stars in pc^-3 (e.g., http://www.pas.rochester.edu/~emamajek/memo_star_dens.html)
n_G = 0.004
n_K = 0.013
n_F = 0.002

# Cutoff distance for the sample (in pc)
d_max = 1000.

# SNR threshold for a significant detection
SNR_thresh = 5

# Period of a typical HZ planet for G,K,M stars (in days)
P_G = 400.
P_K = 80.
P_F = 550.

# Reference V-band absolute magnitudes for each spectral type (X5V stars from http://www.pas.rochester.edu/~emamajek/EEM_dwarf_UBVIJHK_colors_Teff.txt)
MV_G = 4.98
MV_K = 7.36
MV_F = 3.40

# Likelihood of a transit event
p_transit = 0.005

# Number of Earth-twins (e.g., size & in habitable zone) per star
p_Earth = 0.33

# Abundances being probed
species = ['H2O','CO2','O3','O2']

# Reference values from Ehrenreich et al. (2006) for epsilonD and V band magnitude
eps_D_ref = 10. # in meters
V_ref = 8.
# placeholders:
t = 1.
tau = 1.

# Reference significance levels for each species during a single transit observation
SNR_G_ref = [0.8,0.5,1.2,0.2]
SNR_K_ref = [1.7,1.1,1.9,0.2]
SNR_F_ref = [0.5,0.3,0.9,0.1]

# SIMULATION
# Generate the G,K,M stars isotropically for d < d_max (i.e., P(d) proportional to d^2)
# Also modulate by the likelihood of a transiting & habitable zone Earth-like planet
vol = (4.*np.pi/3.)*(d_max/2.)**3
N_G,N_K,N_F = int(n_G*vol*p_transit*p_Earth),int(n_K*vol*p_transit*p_Earth),int(n_F*vol*p_transit*p_Earth)

d0 = np.linspace(0.1,d_max,10000)
p0 = d0**2 / np.sum(d0**2)

d_G,d_K,d_F = np.random.choice(d0,p=p0,size=N_G),np.random.choice(d0,p=p0,size=N_K),np.random.choice(d0,p=p0,size=N_F)

# Calculate the apparent V magnitudes for each
V_G,V_K,V_F = 5*np.log10(d_G)-5+MV_G,5*np.log10(d_K)-5+MV_K,5*np.log10(d_F)-5+MV_F

# Calculate the SNR ratio for each (relative to the values in the Ehr. 2006 table)
r_G,r_K,r_F = (eps_D/eps_D_ref) * (t/tau)**0.5 * 10**(-0.2*(np.array((V_G,V_K,V_F))-8))

# Calculate the significances for each in a single detection
SNR_G = np.array([SNR*r_G for SNR in SNR_G_ref])
SNR_K = np.array([SNR*r_K for SNR in SNR_K_ref])
SNR_F = np.array([SNR*r_F for SNR in SNR_F_ref])

SNR_G1, SNR_K1, SNR_F1 = np.copy(SNR_G), np.copy(SNR_K), np.copy(SNR_F)

# Generate a timeline for the entire mission, in days
t = np.linspace(0,t_dur,20) 

# Calculate how many times each planet has been observed as a function of time, given a random starting phase between -1 and 0
N_tr_G, N_tr_K, N_tr_F = np.random.uniform(-1,0,size=N_G),np.random.uniform(-1,0,size=N_K),np.random.uniform(-1,0,size=N_F)
N_tr_G = np.array([(t0/P_G)+N_tr_G for t0 in t])
N_tr_F = np.array([(t0/P_F)+N_tr_F for t0 in t])
N_tr_K = np.array([(t0/P_K)+N_tr_K for t0 in t])

# Mask out values where N_tr < 0
N_tr_G[N_tr_G<0] = 0.
N_tr_F[N_tr_F<0] = 0.
N_tr_K[N_tr_K<0] = 0.

# To get the detection significance achieved after N_tr transits, multiple SNR by sqrt(N_tr)
SNR_G = np.array([(N_tr_G[i,None,:]**0.5)*SNR_G for i in range(N_tr_G.shape[0])])
SNR_F = np.array([(N_tr_F[i,None,:]**0.5)*SNR_F for i in range(N_tr_F.shape[0])])
SNR_K = np.array([(N_tr_K[i,None,:]**0.5)*SNR_K for i in range(N_tr_K.shape[0])])

# Make a threshold of the number of significant detections vs time in each spectral type
N_det_G = np.array([np.sum(SNR_G[i]>SNR_thresh,axis=1) for i in range(len(t))])
N_det_F = np.array([np.sum(SNR_F[i]>SNR_thresh,axis=1) for i in range(len(t))])
N_det_K = np.array([np.sum(SNR_K[i]>SNR_thresh,axis=1) for i in range(len(t))])

# PLOTS

# Plot parameters
plt.rcParams['font.size']=20.
labelfontsize=34.

# First convert t from days to years
t /= 365.25

# Plot the number of significant detections for each species
fig, ax = plt.subplots(1,3)

Ns = [N_det_F,N_det_G,N_det_K]
for i in range(3):
	for j in range(len(species)):
		ax[i].plot(t,Ns[i][:,j],label=species[j])
	ax[i].set_title('FGK'[i],fontsize=labelfontsize)
	ax[i].set_xlim([0,np.amax(t)])
	ax[i].set_ylim([0,np.amax(Ns[i])+1])
	if i == 0:
		ax[i].legend(loc='upper left')

fig.text(0.5,0.04,"time elapsed (yr)", ha="center", va="center",fontsize=labelfontsize)
fig.text(0.05,0.5,"total SNR > {0} detections".format(int(SNR_thresh)), ha="center", va="center", rotation=90,fontsize=labelfontsize)
plt.show()

# Plot a histogram of the SNR achieved for each species in a single observation (WIP)
fig, ax = plt.subplots(1,3)
sigs = [SNR_F1,SNR_G1,SNR_K1]

ymax = [0,0,0]
for i in range(3):
	for j in range(len(species)):
		y,x=np.histogram(sigs[i][j])
		x = (x[1:]+x[:-1])/2.
		ax[i].plot(x,y,label=species[j])
		ymax[i] = max(ymax[i],y[np.argmin(np.abs(x-1.))])
	ax[i].set_ylim([0,ymax[i]])
plt.show()









