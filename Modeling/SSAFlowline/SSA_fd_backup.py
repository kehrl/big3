# This is a code to solve the SSA. 

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy.signal as signal

##########
# Inputs #
##########

# Load info about CFacier geometry
flowline = np.loadtxt("flowline.dat",skiprows=1)
velocity = np.loadtxt("velocity.dat",skiprows=1)
glacierwidth = np.loadtxt("width.dat",skiprows=1)

# Choose time and timestep
yearinsec = 365.25*24*60*60
dt = 0.5*yearinsec # Years
t = 1.0*yearinsec # Years
Nt = int(np.ceil(t/dt)) # number of timesteps

# Choose node spacing
dx_0 = 2000.0 # desired grid size in meters
dist_0 = flowline[-1,0]-flowline[0,0] # length of flowline

# Inflow boundary conditions
u_inflow = 200.0/yearinsec

# Constants
rho_i = 917.0 # ice density (kg/m3)
rho_sw = 1028.0 # sea water density (kg/m3)
g = 9.8 # gravity (m/s2)

# Flow law parameter
A = 1.7e-24 # Pa-3 s-1

# Cutoff for iterations
cutoff = 1.0e-1 # maximum m/yr

##############
# Mode setup #
##############

# Set up initial nodes and values for velocity and height. We're using a staggered grid for
# number stability, so the height nodes will be located halfway between the 
# the velocity nodes.
norm_nodes_0 = np.linspace(0,dist_0,dist_0/dx_0+1) # Nodes
stag_nodes_0 = norm_nodes_0 + (norm_nodes_0[1]-norm_nodes_0[0])/2 # Staggered nodes
Nnodes = len(norm_nodes_0)

u_0 = np.interp(norm_nodes_0,velocity[:,0],velocity[:,1])/yearinsec # Initial velocities
H_0 = np.interp(norm_nodes_0,flowline[:,0],flowline[:,4]-flowline[:,3]) # Initial heights

# Find glacier width and bed at nodes
zb = np.interp(norm_nodes_0,flowline[:,0],flowline[:,3])
W = np.interp(norm_nodes_0,glacierwidth[:,0],glacierwidth[:,1])

# Decide on a friction coefficient
Beta2 = 2.0*1.0e5*np.interp(norm_nodes_0,glacierwidth[:,0],glacierwidth[:,1]) # (Pa s m-1)

# Set up matrices for results
norm_nodes = np.zeros([len(norm_nodes_0),Nt])
stag_nodes = np.zeros([len(stag_nodes_0),Nt])

u = np.zeros([len(norm_nodes_0),Nt])
u_new = np.zeros([len(norm_nodes_0)])
H = np.zeros([len(norm_nodes_0),Nt])
zs = np.zeros([len(norm_nodes_0),Nt])

norm_nodes[0:Nnodes,0] = norm_nodes_0
stag_nodes[0:Nnodes,0] = stag_nodes_0
u_new[0:Nnodes] = u_0
H[0:Nnodes,0] = H_0

filt_len = 5.0e3
cutoff=(1/filt_len)/(1/(np.diff(flowline[1:3,0])*2))
b,a=signal.butter(4,cutoff,btype='low')
zs_filtered=signal.filtfilt(b,a,flowline[:,4])
zs[0:Nnodes,0] = np.interp(norm_nodes_0,flowline[:,0],zs_filtered)
#zs[0:Nnodes,0] = np.linspace(flowline[0,4],flowline[-1,4],Nnodes)

##################
# Run time steps #
##################

CF_ind = len(norm_nodes_0)-1 # Assume calving front is at terminus
for i in range(0,1):
  
  # Set up specifics for this timestep
  dx = norm_nodes[1,i]-norm_nodes[0,i]
  
  # Define u_last
  u_last = np.zeros(Nnodes)
  
  n=0
  while np.max(abs(u_new-u_last))*yearinsec > cutoff:
    
    # Effective viscosity
    nu = np.zeros([len(norm_nodes)])
    nu[0:-1] = A**(-1.0/3.0)*(abs(u_new[1:]-u_new[0:-1])/dx)**(-2.0/3.0)
    nu[-1] = nu[-2]
  
    # Set up coefficients for matrix problem
    Diagonal = np.zeros(Nnodes) # diagonal terms
    Left = np.zeros(Nnodes) # terms on left of diagonal
    Right = np.zeros(Nnodes) # terms on right side of diagonal
    Vector = np.zeros(Nnodes) # terms on right side of matrix problem
  
    # Find values at staggered nodes
    H_stag = (H[1:,i]+H[0:-1,i])/2.0
    nu_stag = (nu[1:]+nu[0:-1])/2.0
    zs_stag = (zs[1:,i]+zs[0:-1,i])/2.0
  
    # Interior nodes behind calving front
    Left[0:CF_ind-1] = (2.0/dx**2.0)*(H_stag[0:CF_ind-1]*nu_stag[0:CF_ind-1]) 
  
    gamma = np.sign(u_new[1:CF_ind])*(abs(u_new[1:CF_ind]))**(-2.0/3.0)
    Diagonal[1:CF_ind] = -(2.0/dx**2.0)*(H_stag[1:CF_ind]*nu_stag[1:CF_ind]+H_stag[0:CF_ind-1]*nu_stag[0:CF_ind-1]) - \
      Beta2[1:CF_ind]-gamma*(2*H[1:CF_ind,i]/W[1:CF_ind])*(5/(A*W[1:CF_ind]))**(1.0/3.0)     
  
    Right[2:CF_ind+1] = (2.0/dx**2.0)*(H_stag[1:CF_ind]*nu_stag[1:CF_ind])
  
    Vector[1:CF_ind] = rho_i*g*H[1:CF_ind,i]*(zs_stag[1:CF_ind]-zs_stag[0:CF_ind-1])/dx
    # Vector[1:CF_ind] = rho_i*g*H[1:CF_ind,i]*(zs[2:CF_ind+1,i]-zs[1:CF_ind,i])/dx
  
    # Interior nodes in front of calving front
    Left[CF_ind-1:]= -1.0
    Diagonal[CF_ind:] = 1.0
    Vector[CF_ind:] = A*dx*((rho_i*g*H[CF_ind:,i])/4*(1-rho_i/rho_sw))**(3.0)
  
    # Inflow boundary condition
    Diagonal[0] = 1.0 
    Vector[0] = u_inflow
  
    # Create sparse matrix
    Matrix = sparse.dia_matrix(([Left,Diagonal,Right],[-1,0,1]),shape=(Nnodes,Nnodes))

    u_last = np.array(u_new) # From last iteration
    u_new = linalg.spsolve(Matrix,Vector) # New iteration
    
    n=n+1
  
  u[:,i] = u_new