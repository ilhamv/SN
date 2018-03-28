import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

phi = f['scalar_flux_time'];
phi = np.array(phi);
z = f['z'];

f = h5py.File('steady.h5', 'r');
phi0 = f['scalar_flux']

for k in [0,1,2,3,4,8,20,40,100,300,500]:
    plt.plot(z,phi[k], 'o')
    plt.plot(z,phi0)
    plt.title("k = %s"%(k));
    plt.xlabel("z, cm");
    plt.ylabel("Scalar flux, /cm^2.s");
    plt.show();
