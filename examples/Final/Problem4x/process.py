from subprocess import call
import h5py
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Problem 4 Theoretical
#==============================================================================

# Quadrature sets
N = 16    
mu, w = np.polynomial.legendre.leggauss(N)

# The beta function
def beta_f(SigmaT_h):
    return 4.0/3.0 / ( SigmaT_h**2 + 4./3.)

# Cases run
I = 16
beta = np.zeros(I)
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# Beta plot
SigmaT_h[0] = 0.0
beta[0] = beta_f(SigmaT_h[0])

for i in range(1,I):
    SigmaT_h[i] = 0.2 * i
    beta[i] = beta_f(SigmaT_h[i])
plt.plot(SigmaT_h,beta,"*-")
plt.ylabel(r'$\beta$')
plt.xlabel(r'$\Sigma_th$')
plt.show()

# Function of omega (h=0)
def omega_func_zero(lamb,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        num = w[n]
        denom = ( lamb * mu[n] )**-2 + 1.0
        omega = omega + num / denom
    return 1.0 - ( 1.0 + 3.0 / lamb**2 * beta_f(SigmaT_h) ) * omega

# Function of omega
def omega_func(tau,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        coth = np.tanh( SigmaT_h/2/mu[n] )**-1
        num =  2.0*mu[n]/SigmaT_h * coth * w[n]
        denom = np.tan(tau)**-2 + coth**2
        omega = omega + num / denom
    return 1.0 - ( 1.0 + 3.0 * ( SigmaT_h/2.0/np.sin(tau) )**2
                         * beta_f(SigmaT_h) ) * omega

# List of lambda (h=0)
N_lamb = 100
lamb_list = np.linspace(0.0,10.0,N_lamb)

# List of tau
N_tau = 100
tau_list = np.linspace(0.0,np.pi/2.0,N_tau)

# Cases run
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# For h=0
SigmaT_h[0] = 0.0
om=[]
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    om.append(omega)
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)
om = np.array(om)
print(SigmaT_h[0],abs(om[1]),abs(om[-1]),rho[0])
plt.plot(tau_list,om)
plt.show()

# Printout
print(SigmaT_h[0],rho[0])

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 0.2 * i

    # Loop over tau_list and find the maximum rho
    om=[]
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        om.append(omega)
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)
    om = np.array(om)
    print(SigmaT_h[i],abs(om[1]),abs(om[-1]),rho[i])
    plt.plot(tau_list,abs(om))
    plt.show()
    # Printout

# Plot
plt.plot( SigmaT_h, rho, '*-', label="Theory - Relaxed" )

#==============================================================================
# Problem 4 Experimental
#==============================================================================

SigmaT_h = [ 0.2, 0.3, 0.5, 0.6, 1.0, 1.2, 1.5, 2.0, 3.0 ]
mesh = np.zeros(2)
rho = np.zeros_like(SigmaT_h)

args = ["./../../../SN.exe","."]

for i in range(len(SigmaT_h)):
    mesh[0] = 18.0 / SigmaT_h[i]
    mesh[1] = 6.0 / SigmaT_h[i]
    beta_code = beta_f(SigmaT_h[i])
    
    with open('input.xml', 'r') as file:
        data = file.readlines()
    data[5] = "<Accelerator type=\"IDSA\" beta=\"%f\"/>\n"%beta_code
    data[18] = "    <mesh> %i %i %i </mesh>\n"%(mesh[0],mesh[1],mesh[1])
    with open('input.xml', 'w') as file:
        file.writelines( data )
    call(args)
    
    f = h5py.File('output.h5', 'r');
    rho[i] = np.array(f['spectral_radius'])[-1]
    f.close()

plt.plot(SigmaT_h,rho,"o",label="Code - Relaxed")
plt.xlabel(r'$\Sigma_th$')
plt.ylabel(r'$\rho$')
plt.grid()


#==============================================================================
# Problem 3 omega
#==============================================================================

SigmaT_h = np.zeros(I)
rho = np.zeros(I)

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

# For h=0
SigmaT_h[0] = 0.0
om=[]
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    om.append(omega)
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)
#plt.plot(tau_list,om)
#plt.show()

# Printout
print(SigmaT_h[0],rho[0])

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 0.2 * i

    # Loop over tau_list and find the maximum rho
    om=[]
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        om.append(omega)
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)
#plt.plot(tau_list,om)
#plt.show()
    # Printout
    print(SigmaT_h[i],rho[i])

# Plot
plt.plot( SigmaT_h, rho, '-*', label="Theory - Original" )
plt.plot( SigmaT_h, np.ones(I), label=r"$\rho$ = c" )
plt.legend()
plt.show()
