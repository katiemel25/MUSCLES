
# coding: utf-8

# In[107]:

from __future__ import print_function
from scipy.optimize import curve_fit

import matplotlib.pyplot as plt
import numpy as np

import astropy.io.fits as pyfits

get_ipython().magic(u'matplotlib inline')


# In[108]:

hdu = pyfits.open("Spectra/GJ176/ADP.2016-10-28T01:01:09.313.fits")
hdu.info()
hdu[1].header


# In[109]:

wave = hdu[1].data[0][0]
flux = hdu[1].data[0][1]
err = np.nan_to_num(hdu[1].data[0][2])
print("HARPS Error Array: ", err)


# In[110]:

#speed of light, target star RV, H-alpha in meters per second
c = 299792458
RV = 26410
Ha = 6562.8e-10

#accounting for target star RV conversion in angstroms
#RV given in SIMBAD (note: HARPS already corrects for movement of earth)
RV_conv = (Ha * RV / c) * (10 ** 10)


# In[111]:

#fnding mean continuum flux
cont_6500 = 6500 + RV_conv
cont_6550 = 6550 + RV_conv
cont_6575 = 6575 + RV_conv
cont_6625 = 6625 + RV_conv
idx_cont = np.where(((wave > cont_6500) & (wave < cont_6550))|((wave > cont_6575) & (wave < cont_6625)))
sub_cont = flux[idx_cont]
err_cont = np.std(sub_cont)
mean_cont = np.mean(sub_cont)
print("Mean Continuum Flux (Fc): ", mean_cont)
print("Mean Continuum Flux STD: ", err_cont)


# In[112]:

#integrating wavelength to get equivalent width (via eqn 1. in Newton et al. 2017)
idx_wave = np.array(np.where((wave > 6558.8 + RV_conv) & (wave < 6566.8 + RV_conv)))
err_flux = err[idx_wave]
delta_wave = wave[idx_wave[0][400]]-wave[idx_wave[0][399]] #pixel width, assumed to be same for all 
sub_flux = (1-(flux[idx_wave]/mean_cont))*delta_wave
ew = np.sum(sub_flux)
print("Equivalent Width (Angstroms): ", ew)


# In[119]:

#error analysis of equivalent width -- note: assuming pixel width error is negligible
x = flux[idx_wave]
y = mean_cont
z = delta_wave
func_var = np.array([(err_flux[i]**2)*((-y*z)**2)+(err_cont**2)*((x[i]*z/y**2)**2) for i in range(len(sub_flux))])
err_ew = np.sum(func_var)
print("Equivalent Width VAR: ", err_ew)


# In[120]:

plt.plot(wave,flux)
plt.xlim(6558.8+RV_conv,6566.8+RV_conv)
plt.axvline(6562.8+RV_conv) 


# In[ ]:



