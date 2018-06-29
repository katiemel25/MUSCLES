from astroquery.simbad import Simbad
import numpy as np

#star_names = np.genfromtxt("star_names.txt", dtype='str', invalid_raise=False)
customSimbad = Simbad()
customSimbad.add_votable_fields("fluxdata(I)", "fluxdata(K)")
#print customSimbad.get_votable_fields()
table = customSimbad.query_object("GJ 876")
#print table.dtype.names
i_flux = table["FLUX_I"][0]
k_flux = table["FLUX_K"][0]

print "I Flux: ", i_flux, "K Flux: ", k_flux