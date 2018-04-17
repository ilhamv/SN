import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Problem 1a
#==============================================================================

# Quadrature sets
N = 16
mu, w = np.polynomial.legendre.leggauss(N)

# Function of omega
def omega_func(tau,SigmaT_h):
    if SigmaT_h == 0: tau = 0
    omega = 0.0
    for n in range(int(N/2),N):
        num = ( np.cos(tau)**2 - 3.0*mu[n]**2 ) * w[n]
        if SigmaT_h == 0:
            denom = 1.0
        else:
            denom = ( np.cos(tau)**2 
                      + ( 2.0 / SigmaT_h * np.sin(tau) )**2 * mu[n]**2 )
        omega = omega + num / denom
    return omega

# List of tau
N_tau = 100
tau_list = np.linspace(0.0,np.pi/2.0,N_tau)

# Cases run
I = 16
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# Loop over cases
for i in range(I):
    SigmaT_h[i] = 0.2 * i

    # Loop over tau_list and find the maximum rho
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)

    # Printout
    print(SigmaT_h[i],rho[i])

# Plot
plt.plot( SigmaT_h, rho )
plt.plot( SigmaT_h, np.ones(I) )
plt.show()
