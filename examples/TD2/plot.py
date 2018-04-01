import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

f = h5py.File('output.h5', 'r');

phi = np.array(f['scalar_flux_time']);
t = np.array(f['time']);
z = np.array(f['z']);

f = h5py.File('steady.h5', 'r');
phi0 = np.array(f['scalar_flux'])
z_s = np.array(f['z'])

fig = plt.figure()
ax = plt.axes(xlim=(-0.5, 32), ylim=(0.0, 0.3))
line, = ax.plot([], [], 'o', lw=2)
line2, = ax.plot([], [], lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
plt.xlabel("z, cm");
plt.ylabel("Scalar flux, /cm^2.s");

def init():
    line.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return time_text, line, line2

def animate(i):
    line.set_data(z, phi[i])
    line2.set_data(z_s, phi0)
    time_text.set_text('time = %.4f s' %t[i])
    return time_text, line, line2

inter = 5000 / len(phi)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(phi), interval=inter, blit=True)

plt.show()
