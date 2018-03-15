import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

k = f['ksearch/k_cycle'];
plt.plot(k,'x');
plt.xlabel("cycle#");
plt.ylabel("k");
plt.ylim(0)
plt.show(); 
