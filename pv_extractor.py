import pylab as pl
import numpy as np
from astropy.visualization import quantity_support
from astropy.coordinates import Galactic
from astropy import units as u
from astropy import wcs
from spectral_cube import SpectralCube
from pvextractor import extract_pv_slice, Path

# set so that these display properly on black backgrounds
# pl.rcParams['figure.facecolor']='w'
cube = SpectralCube.read('./Output/Program1_13.fits')
g = Galactic([156, 148] * u.deg, [7.7, 1.7] * u.deg)
path = Path(g, width=2 * u.arcsec)
slice1 = extract_pv_slice(cube, path)
slice1.writeto('./my_slice.fits', overwrite=True)
# cube = SpectralCube.read('my_cube.fits')
# slice2 = extract_pv_slice(cube, path)

# slice3 = extract_pv_slice('my_cube.fits', path)

