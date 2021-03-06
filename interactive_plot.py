from __future__ import print_function

import matplotlib.pyplot as plt
import numpy as np

import astropy.io.fits as pyfits

hdu2 = pyfits.open("ADP.2014-09-25T15_35_39.157.fits")
hdu2.info()
hdu2[1].header

wave2 = hdu2[1].data[0][0]
flux2 = hdu2[1].data[0][1]

plt.plot(wave2,flux2)
plt.xlim(3930,3938)
plt.ylim(-50,400)
print("Please click")
x = plt.ginput(2)
print("Wavelength bounds selected: ", x)
plt.show()