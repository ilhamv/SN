import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

f = h5py.File('output.h5', 'r');

phi = np.array(f['scalar_flux_time']);
t = np.array(f['time']);
z = np.array(f['z']);

fig = plt.figure()
ax = plt.axes(xlim=(-0.5, 10.5), ylim=(-2, 10.0))
line, = ax.plot([], [], 'o', lw=2)
line2, = ax.plot([], [], lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
plt.xlabel("z, cm");
plt.ylabel("Scalar flux, /cm^2.s");

Q = 0.1
SigmaA = 0.4
alpha = 1.0
beta = 3.0
A = 2*alpha+2./3*beta-Q/SigmaA
v = 5000
def f_phi(tt):
    return A * np.exp(-v*SigmaA*tt) + Q/SigmaA

def init():
    line.set_data([], [])
    line2.set_data([], [])
    time_text.set_text('')
    return time_text, line, line2

def animate(i):
    line.set_data(z, phi[i])
    line2.set_data(z, np.ones(len(z))*f_phi(t[i]))
    time_text.set_text('time = %.4f s' %t[i])
    return time_text, line, line2

inter = 5000 / len(phi)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(phi), interval=inter, blit=True)

plt.show()

error = abs(f_phi(t[-1])-phi[-1][4])/f_phi(t[-1])

print(error)
