# This is a code to solve the SSA. 

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy.signal as signal

##########
# Inputs #
##########

# Constants
rho_i = 917.0 # ice density (kg/m3)
rho_sw = 1028.0 # sea water density (kg/m3)
g = 9.8 # gravity (m/s2)
yearinsec = 365.25*24*60*60

# Load info about glacier geometry
flowline = np.loadtxt("flowline.dat",skiprows=1)
velocity = np.loadtxt("velocity.dat",skiprows=1)
glacierwidth = np.loadtxt("width.dat",skiprows=1)

# Choose time and timestep
dt = 0.2*yearinsec # Years
t = 0.4*yearinsec # Years
Nt = int(np.ceil(t/dt)) # number of timesteps

# Choose node spacing
dx_0 = 100.0 # desired grid size in meters
L_0 = flowline[-1,0]-flowline[0,0] # length of flowline

# Inflow boundary conditions
u_inflow = 200.0/yearinsec

# Flow law parameter
A = 1.7e-24 # Pa-3 s-1

# Cutoff for iterations
cutoff = 1.0 # maximum m/yr

# Distance in front of ice front
pad = 1.0e3
extent = dist_0+pad

##############
# Mode setup #
##############

# Set up initial nodes and values for velocity and height. We're using a staggered grid for
# number stability, so the height nodes will be located halfway between the 
# the velocity nodes.
norm_nodes_0 = np.linspace(0,extent,extent/dx_0+1) # Nodes
stag_nodes_0 = norm_nodes_0 + (norm_nodes_0[1]-norm_nodes_0[0])/2 # Staggered nodes
Nnodes = len(norm_nodes_0)

# Find glacier width, surface elevation, and bed at nodes
filt_len = 5.0e3
filt_cutoff=(1/filt_len)/(1/(np.diff(flowline[1:3,0])*2))
b,a=signal.butter(4,filt_cutoff,btype='low')
zs_filtered=signal.filtfilt(b,a,flowline[:,4])
zb_filtered=signal.filtfilt(b,a,flowline[:,3])
width_filtered=signal.filtfilt(b,a,glacierwidth[:,1])
u_filtered=signal.filtfilt(b,a,velocity[:,1])
del filt_len,filt_cutoff,b,a

# Interpolate ice thickness and velocity to nodes
u_0 = np.interp(norm_nodes_0,velocity[:,0],u_filtered)/yearinsec # Initial velocities
H_0 = np.interp(norm_nodes_0,flowline[:,0],zs_filtered-zb_filtered) # Initial heights

# Decide on a friction coefficient
x = np.linspace(0,1,Nnodes)
#Beta2 = (2.0e9)*np.linspace(1,0,Nnodes)
Beta2 = 2.0e10*(np.exp(-2.0*x)-np.exp(-2.0))

zb = np.interp(norm_nodes_0,flowline[:,0],zb_filtered)
W = np.interp(norm_nodes_0,glacierwidth[:,0],width_filtered)

# Set up matrices for results
norm_nodes = np.zeros([len(norm_nodes_0),Nt])
stag_nodes = np.zeros([len(stag_nodes_0),Nt])

u = np.zeros([len(norm_nodes_0),Nt])
u_new = np.zeros(Nnodes)
H_stag = np.zeros([len(stag_nodes_0),Nt])
H = np.zeros([len(norm_nodes_0),Nt])
zs = np.zeros([len(norm_nodes_0),Nt])

norm_nodes[0:Nnodes,0] = norm_nodes_0
stag_nodes[0:Nnodes,0] = stag_nodes_0
#H_stag[0:Nnodes-1,0] = H_0[0:-1]
H[:,0] = H_0
H[0,:] = H_0[0]
u_new[0:Nnodes] = u_0

zs[0:Nnodes,0] = np.interp(norm_nodes_0,flowline[:,0],zs_filtered)
# zs[0:Nnodes,0] = np.linspace(zs_filtered[0],zs_filtered[-1],Nnodes)
# zs[0:Nnodes,0] = np.linspace(flowline[0,4],flowline[-1,4],Nnodes)


##################
# Run time steps #
##################


CF_ind = len(norm_nodes_0)-1
for i in range(0,Nt-1):
  
  # Set up specifics for this timestep
  dx = norm_nodes[1,i]-norm_nodes[0,i]
  
  # Set up variables for velocity iterations
  if i > 0:
    u_new = np.interp(norm_nodes[:,i],norm_nodes[:,i-1],u_new)
  else:
    u_new = u_0
  u_last = np.zeros(Nnodes)
  
  # Count the number of iterations
  n=0

  
  while np.max(abs(u_last-u_new))*yearinsec > cutoff:
  
    # Effective viscosity
    nu = np.zeros([len(norm_nodes)])
    nu[0:-1] = A**(-1.0/3)*abs((u_new[1:]-u_new[0:-1])/dx)**(-2.0/3.0)
    nu[-1] = nu[-2]
  
    # Set up coefficients for matrix problem
    Diagonal = np.zeros(Nnodes) # diagonal terms
    Left = np.zeros(Nnodes) # terms on left of diagonal
    Right = np.zeros(Nnodes) # terms on right side of diagonal
    Vector = np.zeros(Nnodes) # terms on right side of matrix problem
  
    # Find values at staggered nodes
    H_stag[0:Nnodes-1,i] = (H[1:,i]+H[0:-1,i])/2
    nu_stag = (nu[1:]+nu[0:-1])/2
  
    # Interior nodes behind calving front
    Left[0:CF_ind-1] = (2.0/dx**2.0)*(H_stag[0:CF_ind-1,i]*nu_stag[0:CF_ind-1]) 
  
    gamma = np.sign(u_new[1:CF_ind])*(abs(u_new[1:CF_ind]))**(-2.0/3.0)
    Diagonal[1:CF_ind] = -(2.0/dx**2.0)*(H_stag[1:CF_ind,i]*nu_stag[1:CF_ind]+H_stag[0:CF_ind-1,i]*nu_stag[0:CF_ind-1]) - \
      Beta2[1:CF_ind]-gamma*(H[1:CF_ind,i]/W[1:CF_ind])*(5.0/(2.0*A*W[1:CF_ind]))**(1.0/3.0)     
  
    Right[2:CF_ind+1] = (2.0/dx**2.0)*(H_stag[1:CF_ind,i]*nu_stag[1:CF_ind])
  
    Vector[1:CF_ind] = rho_i*g*H[1:CF_ind,i]*(zs[2:CF_ind+1,i]-zs[1:CF_ind,i])/dx
  
    # Interior nodes in front of calving front
    Left[CF_ind-1:]=-1.0
    Diagonal[CF_ind:] = 1.0
    Vector[CF_ind:] = A*dx*((rho_i*g*H[CF_ind:,i])/4*(1-rho_i/rho_sw))**(3.0)
  
    # Inflow boundary condition
    Diagonal[0] = 1.0 
    Vector[0] = u_inflow
  
    # Create sparse matrix
    Matrix = sparse.spdiags([Left,Diagonal,Right],[-1,0,1],Nnodes,Nnodes)

    # Save iterated velocities so that we can check to see if the solution has converged
    u_last = np.array(u_new)
    u_new = linalg.spsolve(Matrix,Vector)
    
    # Number of iterations
    n=n+1
  
  # The velocities have converged, so we save the converged solution.
  u[0:Nnodes,i] = u_new
  
  # Update ice thickness
  dH = np.zeros(Nnodes)
  dH = -(1.0/W[0:-1])*(((u[1:,i]*W[1:]*H[1:,i])-(u[0:-1,i]*W[0:-1]*H[0:-1,i]))/dx)*dt
  
  H[1:,i+1] = dH+H[1:Nnodes,i]
  
  # Find new calving front