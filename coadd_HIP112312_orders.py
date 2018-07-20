import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import numpy as np

plt.ion()

input_filename = '/Users/aayoungb/Data/HST/FUMES/HIP112312/HIP112312_FUVe.fits'

spec_hdu = pyfits.open(input_filename)
spec = spec_hdu[1].data
spec_header = spec_hdu[1].header

plt.figure()

min_wave = np.min(spec['WAVELENGTH'][36])
max_wave = np.max(spec['WAVELENGTH'][33])

wave_array = np.arange(min_wave,max_wave,0.013)

order_33 = np.interp(wave_array,spec['WAVELENGTH'][33],spec['FLUX'][33])
order_34 = np.interp(wave_array,spec['WAVELENGTH'][34],spec['FLUX'][34])
order_34 = np.interp(wave_array,spec['WAVELENGTH'][35],spec['FLUX'][35])
order_34 = np.interp(wave_array,spec['WAVELENGTH'][34],spec['FLUX'][34])
order_35 = np.interp(wave_array,spec['WAVELENGTH'][35],spec['FLUX'][35])
order_36 = np.interp(wave_array,spec['WAVELENGTH'][36],spec['FLUX'][36])
error_33 = np.interp(wave_array,spec['WAVELENGTH'][33],spec['ERROR'][33])
error_34 = np.interp(wave_array,spec['WAVELENGTH'][34],spec['ERROR'][34])
error_35 = np.interp(wave_array,spec['WAVELENGTH'][35],spec['ERROR'][35])
error_36 = np.interp(wave_array,spec['WAVELENGTH'][36],spec['ERROR'][36])

plt.plot(wave_array,order_33)
plt.plot(wave_array,order_34)
plt.plot(wave_array,order_35)
plt.plot(wave_array,order_36)

flux_array = np.zeros(len(wave_array))
error_array = np.zeros(len(wave_array))

overlap_33_34 = (wave_array >= np.min(spec['WAVELENGTH'][33])) & (wave_array <= np.max(spec['WAVELENGTH'][34]))
flux_array[overlap_33_34] = (order_33[overlap_33_34] + order_34[overlap_33_34])/2.

overlap_34_35 = (wave_array >= np.min(spec['WAVELENGTH'][34])) & (wave_array <= np.max(spec['WAVELENGTH'][35]))
flux_array[overlap_34_35] = (order_34[overlap_34_35] + order_35[overlap_34_35])/2.

overlap_35_36 = (wave_array >= np.min(spec['WAVELENGTH'][35])) & (wave_array <= np.max(spec['WAVELENGTH'][36]))
flux_array[overlap_35_36] = (order_35[overlap_35_36] + order_36[overlap_35_36])/2.

order_33_mask = wave_array >= np.max(spec['WAVELENGTH'][34])
flux_array[order_33_mask] = order_33[order_33_mask]

order_34_mask = (wave_array >= np.max(spec['WAVELENGTH'][35])) & (wave_array <= np.min(spec['WAVELENGTH'][33]))
flux_array[order_34_mask] = order_34[order_34_mask]

order_35_mask = (wave_array >= np.max(spec['WAVELENGTH'][36])) & (wave_array <= np.min(spec['WAVELENGTH'][34]))
flux_array[order_35_mask] = order_35[order_35_mask]

order_36_mask = wave_array <= np.min(spec['WAVELENGTH'][35])
flux_array[order_36_mask] = order_36[order_36_mask]

plt.plot(wave_array,flux_array)

error_array[order_36_mask] = error_36[order_36_mask]
error_array[order_35_mask] = error_35[order_35_mask]
error_array[order_34_mask] = error_34[order_34_mask]
error_array[order_33_mask] = error_33[order_33_mask]

error_array[overlap_33_34] = np.sqrt(error_33[overlap_33_34]**2 / 4. + error_34[overlap_33_34]**2 / 4.)
error_array[overlap_34_35] = np.sqrt(error_34[overlap_34_35]**2 / 4. + error_35[overlap_34_35]**2 / 4.)
error_array[overlap_35_36] = np.sqrt(error_35[overlap_35_36]**2 / 4. + error_35[overlap_35_36]**2 / 4.)

np.savetxt('HIP112312_FUV_coadd.txt',np.transpose(np.array([wave_array,flux_array,error_array])))
