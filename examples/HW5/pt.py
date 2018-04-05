import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

mu = np.array(f['mu'])
phi = np.array(f['scalar_flux'])

pt = 
