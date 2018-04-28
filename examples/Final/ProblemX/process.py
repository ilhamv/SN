from subprocess import call
import h5py
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Problem 1a
#==============================================================================

# Quadrature sets
N = 16
mu, w = np.polynomial.legendre.leggauss(N)

# Function of omega (h=0)
def omega_func_zero(lamb,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        num = ( 1.0 - 3.0 * mu[n]**2 ) * w[n]
        denom = 1.0 + lamb**2 * mu[n]**2
        omega = omega + num / denom
    return omega

# Function of omega
def omega_func(tau,SigmaT_h):
    if SigmaT_h == 0: tau = 0
    omega = 0.0
    for n in range(int(N/2),N):
        num = ( np.cos(tau)**2 * ( 1.0 - 3.0*mu[n]**2 ) ) * w[n]
        if SigmaT_h == 0:
            denom = 1.0
        else:
            denom = ( np.cos(tau)**2 
                      + ( 2.0 / SigmaT_h * np.sin(tau) )**2 * mu[n]**2 )
        omega = omega + num / denom
    return omega

# List of lambda (h=0)
N_lamb = 100
lamb_list = np.linspace(0.0,10.0,N_lamb)

# List of tau
N_tau = 100
tau_list = np.linspace(0.0,np.pi/2.0,N_tau)

# Cases run
I = 16
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# For h=0
SigmaT_h[0] = 0.0
om = []
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    om.append(abs(omega))
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)
# Printout
print(SigmaT_h[0],rho[0])
plt.plot(tau_list,om)
plt.title(r'$\Sigma_th$ = %.3f'%SigmaT_h[0])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$|\omega|$')
plt.grid()
plt.show()

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 3.0 / (I-i)
    om = []

    # Loop over tau_list and find the maximum rho
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        om.append(abs(omega))
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)

    # Printout
    print(SigmaT_h[i],rho[i])
    if i == 1 or i == 12 or i == 15:
        plt.plot(tau_list,om)
        plt.title(r'$\Sigma_th$ = %.3f'%SigmaT_h[i])
        plt.xlabel(r'$\tau$')
        plt.ylabel(r'$|\omega|$')
        plt.grid()
        plt.show()

# Plot
plt.plot( SigmaT_h, rho, '*-', label="Theory" )


#==============================================================================
# Problem 1b
#==============================================================================

# Cases run
I = 15
rho = np.zeros(I)
SigmaT_h = np.zeros(I)
rho = np.zeros_like(SigmaT_h)
mesh = np.zeros(2)

args = ["./../../../SN.exe","."]

for i in range(len(SigmaT_h)):
    SigmaT_h[i] = 3.0 / (I-i)
    mesh[0] = 18.0 / SigmaT_h[i]
    mesh[1] = 6.0 / SigmaT_h[i]
    
    with open('input.xml', 'r') as file:
        data = file.readlines()
    data[18] = "    <mesh> %i %i %i </mesh>\n"%(mesh[0],mesh[1],mesh[1]);
    with open('input.xml', 'w') as file:
        file.writelines( data )
    call(args)
    
    f = h5py.File('output.h5', 'r');
    rho[i] = np.array(f['spectral_radius'])[-1]
    f.close()


print(SigmaT_h,rho)
plt.plot(SigmaT_h,rho,"o",label="Code")
plt.xlabel(r'$\Sigma_th$')
plt.ylabel(r'$\rho$')
plt.grid()
plt.legend()
plt.show()
