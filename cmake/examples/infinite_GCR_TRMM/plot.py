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

flux = f['MG/flux/mean'];
energy = f['MG/energy'];

energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

plt.loglog(energy,flux[0]);
#plt.ylim(1,1E9);
#plt.xlim(1E-9,2E1);
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.grid();
plt.show(); 

x = f['TRMM/alpha/real'];
y = f['TRMM/alpha/imag'];

plt.plot(x,y,'o');
plt.xlabel("Re");
plt.ylabel("Im");
plt.grid();
plt.show(); 

flux = f['TRMM/flux'];
energy = f['MG/energy'];

energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

plt.loglog(energy,flux[0]);
plt.loglog(energy,flux[1]);
plt.loglog(energy,flux[2]);
plt.loglog(energy,flux[3]);
plt.loglog(energy,flux[4]);
plt.ylim(1,);
plt.xlim(1E-9,2E1);
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.grid();
plt.show();
