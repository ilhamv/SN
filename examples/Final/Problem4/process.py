from subprocess import call
import h5py
import numpy as np
import matplotlib.pyplot as plt

#==============================================================================
# Problem 4 Theoretical
#==============================================================================

c = 1.0

# Quadrature sets
N = 16    
mu, w = np.polynomial.legendre.leggauss(N)

# The beta function
def beta_f(SigmaT_h):
    tanh = 0
    coth = 0
    for n in range(int(N/2),N):
        tanh = tanh + mu[n] * np.tanh(SigmaT_h/2/mu[n]) * w[n]
        coth = coth + mu[n] * np.tanh(SigmaT_h/2/mu[n])**-1 * w[n]
    return ( 2.0 - 2./SigmaT_h*tanh ) / ( 1 - ((c-1-c*2/SigmaT_h*tanh)/(1-c+1./3*(2./SigmaT_h)**2)) )
#    return ( 2.0 - 2./SigmaT_h*tanh ) / (3.0*SigmaT_h/2)/( coth + tanh )

# Cases run
I = 16
beta = np.zeros(I)
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# Beta plot
SigmaT_h[0] = 0.0
beta[0] = 1.0

for i in range(1,I):
    SigmaT_h[i] = 3.0 / (I-i)
    beta[i] = beta_f(SigmaT_h[i])
plt.plot(SigmaT_h,beta,"*-")
plt.ylabel(r'$\beta$')
plt.xlabel(r'$\Sigma_th$')
plt.show()
SigmaTh_beta = SigmaT_h

# Function of omega (h=0)
def omega_func_zero(lamb,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        num = w[n]
        denom = ( lamb * mu[n] )**-2 + 1.0
        omega = omega + num / denom
    return c / (1-c+1./3*lamb**2) * ( 1./3*lamb**2 - ( 1./3 * lamb**2 + 1.0 ) * omega )

# Function of omega
def omega_func(tau,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        coth = np.tanh( SigmaT_h/2/mu[n] )**-1
        num =  2.0*mu[n]/SigmaT_h * coth * w[n]
        denom = np.tan(tau)**-2 + coth**2
        omega = omega + num / denom
    return c / ( 1-c+1./3*(2./SigmaT_h*np.sin(tau))**2 ) * ((1+(beta_f(SigmaT_h)-1)*c+1./3*(2./SigmaT_h*np.sin(tau))**2)*(1-omega)-beta_f(SigmaT_h))

# Function of omega (h=0)
def omega_func_zero_ori(lamb,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        num = w[n]
        denom = ( lamb * mu[n] )**-2 + 1.0
        omega = omega + num / denom
    return c * ( 1.0 - ( 1.0 + 3.0 / lamb**2 ) * omega )

# Function of omega
def omega_func_ori(tau,SigmaT_h):
    omega = 0.0
    for n in range(int(N/2),N):
        coth = np.tanh( SigmaT_h/2/mu[n] )**-1
        num =  2.0*mu[n]/SigmaT_h * coth * w[n]
        denom = np.tan(tau)**-2 + coth**2
        omega = omega + num / denom
    return c * ( 1.0 - ( 1.0 + 3.0 * ( SigmaT_h/2.0/np.sin(tau) )**2 ) * omega )

# List of lambda (h=0)
N_lamb = 1000
lamb_list = np.linspace(0.0,100.0,N_lamb)

# List of tau
N_tau = 100
tau_list = np.linspace(0.0,np.pi/2.0,N_tau)

# Cases run
I = 16
rho = np.zeros(I)
SigmaT_h = np.zeros(I)

# For h=0
SigmaT_h[0] = 0.0
om=[]
om_ori=[]
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    omega_ori = omega_func_zero_ori(lamb,SigmaT_h[0])
    om.append(abs(omega))
    om_ori.append(abs(omega_ori))
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)
# Printout
print(SigmaT_h[0],rho[0])
plt.plot(lamb_list,om,label="Relaxed")
plt.plot(lamb_list,om_ori,'--',label="Original")
plt.title(r'$\Sigma_th$ = %.3f'%SigmaT_h[0])
plt.xlabel(r'$\lambda$')
plt.ylabel(r'$|\omega|$')
plt.grid()
plt.legend()
plt.show()

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 3.0 / (I-i)
    om = []
    om_ori = []

    # Loop over tau_list and find the maximum rho
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        omega_ori = omega_func_ori(tau,SigmaT_h[i])
        om.append(abs(omega))
        om_ori.append(abs(omega_ori))
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)
    # Printout
    print(SigmaT_h[i],rho[i])
    
    plt.plot(tau_list,om,label="Relaxed")
    plt.plot(tau_list,om_ori,'--',label="Original")
    plt.title(r'$\Sigma_th$ = %.3f'%SigmaT_h[i])
    plt.xlabel(r'$\tau$')
    plt.ylabel(r'$|\omega|$')
    plt.grid()
    plt.legend()
    plt.show()
    

# Plot
plt.plot( SigmaT_h, rho, '*-', label="Theory - Relaxed" )

#==============================================================================
# Problem 4 Experimental
#==============================================================================

# Cases run
I = 15
rho = np.zeros(I)
SigmaT_h = np.zeros(I)
rho = np.zeros_like(SigmaT_h)
mesh = np.zeros(2)

args = ["./../../../SN.exe","."]

with open('input.xml', 'r') as file:
    data = file.readlines()
data[11] = "        <scatter xs=\"%f\"/>\n"%c
with open('input.xml', 'w') as file:
    file.writelines( data )

for i in range(len(SigmaT_h)):
    SigmaT_h[i] = 3.0 / (I-i)
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
plt.ylabel(r'$\rho$ or $\beta$')
plt.grid()


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
for lamb in lamb_list:
    omega = omega_func_zero(lamb,SigmaT_h[0])
    if abs(omega) > rho[0]:
        rho[0] = abs(omega)

# Loop over cases
for i in range(1,I):
    SigmaT_h[i] = 3.0 / (I-i)

    # Loop over tau_list and find the maximum rho
    for tau in tau_list:
        omega = omega_func(tau,SigmaT_h[i])
        if abs(omega) > rho[i]:
            rho[i] = abs(omega)

# Plot
plt.plot( SigmaT_h, rho, '*-', label="Theory - Original" )
plt.plot( SigmaT_h, np.ones(I)*c, label=r"$\rho$ = c" )
plt.plot(SigmaTh_beta,beta,"*-",label=r'$\beta$')
plt.legend(loc=0)
plt.show()
