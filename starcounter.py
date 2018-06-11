# Calculate number of stars, planets, transiting planets as a function of
# spectral type and distance, based on stellar number density
# Then calculate the apparent magnitudes of each star, which can be used to
# scale the Phoenix magnitudes
# Daniel Apai, April 2018
# ==========================================================================
import numpy as np
import math
import matplotlib.pyplot as plt
from astropy import constants as const
import pdb
 # Calculates the inner and outer boundaries of the HZ based on Kopparapu et al. 2013

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def distmod(distance):
    return 5*np.log10(distance)-5;


def hzboundaries(teff,lum,msgs):
 if teff > 7200:
      return np.NaN,np.NaN;
 elif(teff < 7200):
     seff = [0,0,0,0,0,0]
     seffsun  = [1.776,1.107, 0.356, 0.320, 1.188, 0.99]
     a = [2.136e-4, 1.332e-4, 6.171e-5, 5.547e-5, 1.433e-4, 1.209e-4]
     b = [2.533e-8, 1.580e-8, 1.698e-9, 1.526e-9, 1.707e-8, 1.404e-8]
     c = [-1.332e-11, -8.308e-12, -3.198e-12, -2.874e-12, -8.968e-12, -7.418e-12]
     d = [-3.097e-15, -1.931e-15, -5.575e-16, -5.011e-16, -2.084e-15, -1.713e-15]

     tstar = teff - 5780.0
     for i in range(len(a)):
        seff[i] = seffsun[i] + a[i]*tstar + b[i]*tstar**2 + c[i]*tstar**3 + d[i]*tstar**4

     if msgs == 1:
       print 'Conservative HZ Limits'
       print 'Inner HZ - Recent Venus Limit:', (lum/seff[1])**0.5, 'au'
       print 'Outer HZ - Early Mars Limit:', (lum/seff[2])**0.5, 'au'
       print 'Optimistic HZ Limits'
       print 'Inner HZ - Runaway Greenhouse Limit:', (lum/seff[0])**0.5, 'au'
       print 'Outer HZ - Maximum Greenhouse Limit:', (lum/seff[3])**0.5, 'au'
# This script and results were compared against the online calculator at VPL:
#https://depts.washington.edu/naivpl/sites/default/files/hz.shtml


     rin=(lum/seff[0])**0.5
     rout=(lum/seff[3])**0.5

     return rin, rout;


# Stellar densities based on the 10 pc RECONS sample: http://www.recons.org/census.posted.htm
#
#WDs               18        20        20       white dwarfs
#Os                 0         0         0
#Bs                 0         0         0
#As                 4         4         4
#Fs                 6         6         6
#Gs                20        20        20
#Ks                44        44        44
#Ms               198       246       248       25% increase since 2000
#Ls                 0         4         5
#Ts                 1         9        10

TargetSampleSize=1000 # Set the size of the desired target sample (number of transiting, Earth-sized planets)

recons_sptype=['O','B','A','F','G','K','M','L']
recons_numbers=[0,0,4,6,20,44,248,5]
recons_density=np.array([0.,0.,0.,0.,0.,0.,0.,0.]) # Number density of stars by sp. type: star pc^-3
tenpcvol=(4./3.*np.pi*10**3)
print tenpcvol

for i in range(len(recons_sptype)):
   recons_density[i]=recons_numbers[i]/tenpcvol

#recons_density=recons_numbers/tenpcvol
#recons_density=np.array(recons_numbers/tenpcvol)

# Luminosity and effective temperature values from Mamajek's table [Pecaut & Mamajek (2013, ApJS, 208, 9;]
loglum=[4.99,2.49,1.06,0.36,-0.12,-0.98,-3.09,-4.24] # The typical luminosity of the stars as a function of sp type
teff=[34500.,14000.,7800.,6240.,5530.,4070.,2710,1590.] # The typical luminosity of the stars as a function of sp type
Mstar=[22.9,7.9,1.76,1.21,0.96,0.65,0.10,0.050]
Rstar=[8.75,2.99,1.86,1.30,0.949,0.654,0.118,0.099]
Mv=[-4.8,-0.4,2.07,3.87,5.18,8.15,17.07,float('nan')]
VIc=[0.,-0.133,0.242,0.579,0.766,1.529,4.31,float('nan')]
sptype=['O8V','B7V','A7V','F7V','G7V','K7V','M6.5','L5']
#         0      1    2      3    4     5     6     7
rin=np.array([0.,0.,0.,0.,0.,0.,0,0.]) # Array to store the inner radius of the habitable zone
rout=np.array([0.,0.,0.,0.,0.,0.,0,0.]) # Array to store the outer radius of the habitable zone
period_hz_in=np.array([0.,0.,0.,0.,0.,0.,0,0.]) # Array stores the orbital periods of the inner edge of the HZ
period_hz_out=np.array([0.,0.,0.,0.,0.,0.,0,0.])

# Assumptions about planet occurrence rates: only considering earth-sized planets here:
eta_earth=[0.1,0.2,0.3,0.3,0.3,0.6,1.5,0.2]

transit_prob_out=np.array([0.,0.,0.,0.,0.,0.,0,0.])
transit_prob_in = np.array([0.,0.,0.,0.,0.,0.,0,0.])

msgs=0 # If this is on, detailed messages will be shown and test cases will be evaluated, too
if msgs == 1: #Detailed Messages/test cases?
    print '----------------------------'
    print 'Test case: Sun'
    msgs=1
    rintest, routtest = hzboundaries(5780.,1.0,msgs)
    print 'Another Test case: K dwarf'
    rintest, routtest = hzboundaries(4000.,0.6,msgs)
    print "Optimistic Habitable Zone:"

print "Cycling through the stellar types to calculate optimistic habitable zone boundaries"
i=0
msgs=0
while i <= 7: # Calculate HZ boundaries and periods for sp types
      print 'Caclculating', sptype[i]
      rin[i], rout[i] = hzboundaries(teff[i],10**loglum[i],msgs)
      #print rin[i], rout[i]
      #print math.isnan(rin[i]), math.isnan(rout[i])
      if not math.isnan(rin[i]) and not math.isnan(rout[i]):
          period_hz_in[i] = np.sqrt( (rin[i]* const.au.value)**3 * 4*(np.pi**2)/const.G.value /(Mstar[i]*const.M_sun.value)) /3600./24. #in days
          period_hz_out[i] = np.sqrt( (rout[i]* const.au.value)**3 * 4*(np.pi**2)/const.G.value /(Mstar[i]*const.M_sun.value)) /3600./24. #in days
          transit_prob_in[i] =  ( Rstar[i]*const.R_sun.value)/(rin[i]* const.au.value)
          transit_prob_out[i] = ( Rstar[i]*const.R_sun.value)/(rout[i]* const.au.value)

          print "-----------"
          print sptype[i],period_hz_in[i], period_hz_out[i], 'days'
          print 'Transit probabilities', transit_prob_in[i], transit_prob_out[i]
          if msgs ==1:
              print "Teff", teff[i], "Lum", 10**loglum[i]
              print "Rin", rin[i], "Rout", rout[i]
      else:
          print 'No valid HZ boundaries for this spectral type'
      print '...'
      i=i+1

print 'Habitable Zone boundaries'
print 'Inner radii:', rin
print 'Outer radii:', rout



#rin.replace(-1,np.NaN)

# Create a plot showing the number of transiting HZ planets as a function of distance around stars with different spectral types

dstep = 1.0 # Distance step
maxdist=400.0 #Maximum distance to consider
distance=np.arange(2,maxdist,dstep) # The range variable for plots, distance measured in pc



# Calculate the number of transiting planets:
# numtr = np.array((len(rin),maxdist/dstep))  # A 2-D array to store number of transiting planets as a function of sp type and distance
#numtr = [[0 for x in range(rin)] for y in range(maxdist/dstep)]
numtr = np.zeros((rin.size,distance.size))
apparentmag=np.zeros((rin.size,distance.size))
absdepth=np.zeros((rin.size,distance.size))

Ic_ZP=9.68e-6 # ergs cm-2 s-1  , from Rydgren, Schmelz, Zak, & Vrba, U.S. Naval Observatory PublicaBons XXV, Pt. 1 (1984)


numtr.shape
for sp in range(len(rin)):
    for di in range(len(distance)):
        numtr[sp,di]= 4*np.pi*(distance[di])**3 /3. * recons_density[sp] * eta_earth[sp] * (transit_prob_in[sp]+transit_prob_out[sp])/2.
        apparentmag[sp,di] = (Mv[sp]-VIc[sp]) + distmod(distance[di])
        absdepth[sp,di] = 10**(-0.4*apparentmag[sp,di]+Ic_ZP) * (const.R_earth.value/(Rstar[sp]*const.R_sun.value))**2

totalplanets = numtr[3,]+numtr[4,]+numtr[5,]+numtr[6,]
fgkplanets = numtr[3,]+numtr[4,]+numtr[5,]

# The index of the distance step where there are enough M dwarfs for 50% of the sample
fifty_mswitch_idx=find_nearest(numtr[6,:],TargetSampleSize/2.)

mplanets = numtr[6,:].copy()
mplanets[mplanets >TargetSampleSize/2.]=TargetSampleSize/2.
fiftyfifty = fgkplanets + mplanets

# Find distances that encompass different samples
print find_nearest(numtr[6,:],TargetSampleSize)
mplanets_dist=distance[find_nearest(numtr[6,:],TargetSampleSize)]
fiftyfifty_dist=distance[find_nearest(fiftyfifty,TargetSampleSize)]
fgkplanets_dist=distance[find_nearest(fgkplanets,TargetSampleSize)]
sunlike_dist=distance[find_nearest(numtr[4,],TargetSampleSize)]

print 'Distances:', mplanets_dist, fiftyfifty_dist, fgkplanets_dist,sunlike_dist


#= Create plots
plt.figure(1)
plt.subplot(2,1,1)

# red dashes, blue squares and green triangles
#plt.plot(t, t, 'r--', t, t**2, 'bs', t, t**3, 'g^')
plt.axis([0, 400, 0, 1200.])
#plt.yscale('log')
#plt.xscale('log')
plt.plot(distance,numtr[3,],'blue',marker='o')
plt.plot(distance,numtr[3,],'blue')
plt.plot(distance,numtr[4,],'yellow',marker='o')
plt.plot(distance,numtr[4,],'yellow')
plt.plot(distance,numtr[5,],'orange',marker='o')
plt.plot(distance,numtr[5,],'orange')
plt.plot(distance,numtr[6,],'red',marker='o')
plt.plot(distance,numtr[6,],'red')
plt.plot(distance,totalplanets,'bo',distance,totalplanets,'b')
plt.plot(distance,fgkplanets,'black',distance,fgkplanets,'black')
plt.plot(distance,fiftyfifty,'black',distance,fiftyfifty,'black')


plt.plot([mplanets_dist,mplanets_dist],[0,TargetSampleSize],'--')
plt.plot([fiftyfifty_dist,fiftyfifty_dist],[0,TargetSampleSize],'--')
plt.plot([fgkplanets_dist,fgkplanets_dist],[0,TargetSampleSize],'--')
plt.plot([sunlike_dist,sunlike_dist],[0,TargetSampleSize],'--')

plt.text(distance[20]-15, numtr[6,20], sptype[6])
plt.text(distance[120]+8, numtr[5,120], sptype[5])
plt.text(distance[140]+8, numtr[4,140]+70, sptype[4])
plt.text(distance[350], numtr[3,350]+70, sptype[3])
plt.text(distance[155]-25, fgkplanets[155], 'FGK Combined')
plt.text(distance[55]+8, fiftyfifty[55]-30, '50% M 50% FGK')

plt.plot([0,400],[TargetSampleSize/2.,TargetSampleSize/2.],'g-')
plt.plot([0,400],[TargetSampleSize,TargetSampleSize],'g-')

#plt.text(sp_t, sp_m, sp_name)
plt.title('Cumulative Number of Transiting HZ Earth-sized Planets')
plt.ylabel('Number of Planets')
plt.xlabel('Distance from Sun [pc]')

#Calculate the signal to be detected for an Earth-sized planet as a function of stellar spectral type and Distance
#
plt.subplot(2,1,2)
plt.yscale('log')
plt.axis([0, 400, 1e-10, 1e-6])
#plt.xscale('log')
plt.title('Absolute Transit Signal vs. Host Star Spectral Type and Distance')
plt.ylabel('Intensity Drop During Earth Transit [ergs cm-2 s-1]]')
plt.xlabel('Distance from Sun [pc]')
plt.plot(distance,absdepth[3,],'blue',distance,absdepth[3,],'blue')
plt.plot(distance,absdepth[4,],'y',distance,absdepth[4,],'y')
plt.plot(distance,absdepth[5,],'orange',distance,absdepth[5,],'orange')  #
plt.plot(distance,absdepth[6,],'r',distance,absdepth[6,],'r') # M dwarfs

# Draw vertical lines to mark the mamimum distances in each sample
plt.plot([mplanets_dist,mplanets_dist],[10*absdepth[0,1],0.01*absdepth[3,1]],'--')
plt.plot([fiftyfifty_dist,fiftyfifty_dist],[10*absdepth[0,1],0.005*absdepth[3,1]],'--')
plt.plot([fgkplanets_dist,fgkplanets_dist],[10*absdepth[0,1],0.005*absdepth[3,1]],'--')
plt.plot([sunlike_dist,sunlike_dist],[10*absdepth[0,1],0.005*absdepth[3,1]],'--')

# Draw horizontal lines to mark the smallest absolute signal in each sample
idx_mstars=find_nearest(numtr[6,:],TargetSampleSize)
idx_fifty=find_nearest(fiftyfifty,TargetSampleSize)
idx_fgk=find_nearest(fgkplanets,TargetSampleSize)
idx_sunlike=find_nearest(numtr[4,],TargetSampleSize)
plt.plot([0,mplanets_dist],[absdepth[6,idx_mstars],absdepth[6,idx_mstars]],'--')
plt.plot([0,fiftyfifty_dist],[absdepth[6,fifty_mswitch_idx],absdepth[6,fifty_mswitch_idx]],'--')
plt.plot([0,fgkplanets_dist],[absdepth[5,idx_fgk],absdepth[5,idx_fgk]],'--')
plt.plot([0,sunlike_dist],[absdepth[4,idx_sunlike],absdepth[4,idx_sunlike]],'--')

plt.text(distance[20]-15, absdepth[6,20], sptype[6])
plt.text(distance[120]+8, absdepth[5,120], sptype[5])
plt.text(distance[170]+2, absdepth[4,170], sptype[4])
plt.text(distance[340]+10, absdepth[3,340], sptype[3])

plt.show()


# Write out the compositions of the four samples for the concept paper:
#   \caption{Properties of possible target samples. \label{T:TargetSample}}
#   \tablehead{\colhead{Sample}&\colhead{\#M}&\colhead{\#K}&\colhead{\#G} &
#   \colhead{\#F}  & \colhead{Max. Dist. [pc]} & \colhead{V-mag range} & \colhead{I-mag range}& \colhead{HZE Period [day]} }
# ======================

#Sample 1
print 'Sample 1 - M star planets only'
print 'Max Distance:', mplanets_dist, 'pc'
print 'M dwarf planets:', TargetSampleSize
print 'K dwarf planets:', 0.
print 'G dwarf planets:', 0.
print 'F dwarf planets:', 0.
print 'V-band apparent magnitude of faintest target'
print 'I-band apparent magnitude of faintest target', apparentmag[6,int(mplanets_dist-1)]
print 'HZ Periods:', sptype[6],period_hz_in[6], period_hz_out[6], 'days'

print
print

#Sample 2
print 'Sample 2 - 50% M star planets and 50% FGK star planets'
print 'Max Distance:', fiftyfifty_dist, 'pc'
print 'M dwarf planets:', TargetSampleSize/2.
print 'K dwarf planets:', round(numtr[5,int(fiftyfifty_dist-1)])
print 'G dwarf planets:', round(numtr[4,int(fiftyfifty_dist-1)])
print 'F dwarf planets:', round(numtr[3,int(fiftyfifty_dist-1)])
print 'V-band apparent magnitude of faintest target'
print 'I-band apparent magnitude of faintest M dwarf', apparentmag[6,fifty_mswitch_idx-1]
print 'I-band apparent magnitude of faintest K dwarf', apparentmag[5,int(fiftyfifty_dist-1)]
print 'I-band apparent magnitude of faintest G dwarf', apparentmag[4,int(fiftyfifty_dist-1)]
print 'I-band apparent magnitude of faintest F dwarf', apparentmag[3,int(fiftyfifty_dist-1)]
print
print
print 'HZ Periods:', sptype[6],period_hz_in[6], period_hz_out[6], 'days'
print 'HZ Periods:', sptype[5],period_hz_in[5], period_hz_out[5], 'days'
print 'HZ Periods:', sptype[4],period_hz_in[4], period_hz_out[4], 'days'
print 'HZ Periods:', sptype[3],period_hz_in[3], period_hz_out[3], 'days'

#Sample 3
print 'Sample 3 - Closest FGK star planets'
print 'Max Distance:', fgkplanets_dist, ' pc'
print 'M dwarf planets:', 0
print 'K dwarf planets:', round(numtr[5,int(fgkplanets_dist-1)])
print 'G dwarf planets:', round(numtr[4,int(fgkplanets_dist-1)])
print 'F dwarf planets:', round(numtr[3,int(fgkplanets_dist-1)])
print 'V-band apparent magnitude of faintest target', apparentmag[sp,di]
print 'I-band apparent magnitude of faintest target', apparentmag[sp,di]
print 'I-band apparent magnitude of faintest K dwarf', apparentmag[5,int(fgkplanets_dist-1)]
print 'I-band apparent magnitude of faintest G dwarf', apparentmag[4,int(fgkplanets_dist-1)]
print 'I-band apparent magnitude of faintest F dwarf', apparentmag[3,int(fgkplanets_dist-1)]
print 'HZ Periods:', sptype[5],period_hz_in[5], period_hz_out[5], 'days'
print 'HZ Periods:', sptype[4],period_hz_in[4], period_hz_out[4], 'days'
print 'HZ Periods:', sptype[3],period_hz_in[3], period_hz_out[3], 'days'
print
print

#Sample 4
print 'Sample 4 - G star planets only'
print 'Max Distance:', sunlike_dist, ' pc'
print 'M dwarf planets:',0
print 'K dwarf planets:',0
print 'G dwarf planets:', TargetSampleSize
print 'F dwarf planets:',0
print 'V-band apparent magnitude of faintest target'
print 'I-band apparent magnitude of faintest target'
print 'I-band apparent magnitude of faintest G dwarf', apparentmag[4,int(sunlike_dist-1)]
print 'HZ Periods:', sptype[4],period_hz_in[4], period_hz_out[4], 'days'
