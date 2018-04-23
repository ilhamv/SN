from subprocess import call
import h5py
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Problem 3a
#==============================================================================

# Quadrature sets
N = 16    
mu, w = np.polynomial.legendre.leggauss(N)

# Function of omega (h=0)
def omega_func_zero(lamb,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        num = w[n]
        denom = ( lamb * mu[n] )**-2 + 1.0
        omega = omega + num / denom
    return 1.0 - ( 1.0 + 3.0 / lamb**2 ) * omega

# Function of omega
def omega_func(tau,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        coth = np.tanh( SigmaT_h/2/mu[n] )**-1
        num =  2.0*mu[n]/SigmaT_h * coth * w[n]
        denom = np.tan(tau)**-2 + coth**2
        omega = omega + num / denom
    return 1.0 - ( 1.0 + 3.0 * ( SigmaT_h/2.0/np.sin(tau) )**2 ) * omega

# List of lambda (h=0)
N_lamb = 100
lamb_list = np.linspace(0.0,np.pi/2.0,N_lamb)

# List of tau
N_tau = 100
tau_list = np.linspace(0.0,np.pi/2.0,N_tau)

# Cases run
I = 16
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# For h=0
SigmaT_h[0] = 0.0
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)

# Printout
print(SigmaT_h[0],rho[0])

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 0.2 * i

    # Loop over tau_list and find the maximum rho
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)

    # Printout
    print(SigmaT_h[i],rho[i])

# Plot
plt.plot( SigmaT_h, rho, '*-', label="Theory" )
plt.plot( SigmaT_h, np.ones(I), label=r"$\rho$ = c" )


#==============================================================================
# Problem 3b
#==============================================================================

SigmaT_h = [ 0.2, 0.3, 0.5, 0.6, 1.0, 1.2, 1.5, 2.0, 3.0 ]
mesh = np.zeros(2)
rho = np.zeros_like(SigmaT_h)

args = ["./../../../SN.exe","."]

for i in range(len(SigmaT_h)):
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
