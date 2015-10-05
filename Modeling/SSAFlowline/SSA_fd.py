# This is a code to solve the SSA. 

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy.signal as signal
import sys
import os


###########
# Options #
###########
CF_move = 0

##########
# Inputs #
##########

# Constants
rho_i = 917.0 # ice density (kg/m3)
rho_sw = 1028.0 # sea water density (kg/m3)
rho_w = 1000.0
g = 9.8 # gravity (m/s2)
yearinsec = 365.25*24*60*60

# Load info about glacier geometry
DIR="/Users/kehrl/Models/Helheim/Meshes/Flowline/MorlighemNew_SmoothedVelocity/Inputs/"
flowline = np.loadtxt(DIR+"flowline.dat",skiprows=1)
velocity = np.loadtxt(DIR+"velocity.dat",skiprows=1)
glacierwidth = np.loadtxt(DIR+"width.dat",skiprows=1)

# Choose time step and length of run
dt = (1.0/6.0)/(365.0)*yearinsec # Years
t = 5.0*yearinsec # Years
Nt = int(np.ceil(t/dt)) # number of timesteps

# Choose node spacing
dx_0 = 200.0 # desired grid size in meters
L_0 = 85115.0 # initial length of flowline 
L_0_ind = np.argmin(abs(flowline[:,0]-L_0)) # index for this ice front position in flowline

# Set inflow boundary condition (Dirichlet BC)
u_inflow = 200.0/yearinsec

# Flow law parameter
A = 1.7e-24 # Pa-3 s-1

# Calving law parameters
d_w = 10.0 # water depth in crevasses (m)

# Cutoff for iterations
cutoff = 10.0 # change between iterations in m/yr
nmax = 1e3 # maximum number of iteration, need to add to code

# Padding
# In case the glacier advances, we need to add some extra same in our matrices so that 
# we don't run out of room.
pad = np.ceil((5.0e4*t/yearinsec)/dx_0) # number of nodes

###############
# Model setup #
###############

# Set up initial nodes and values for velocity and height. We're using a staggered grid 
# for numerical stability. We calculate height on the normal grid and velocity/flux on the
# staggered grid.
Nnodes = np.ceil((L_0)/dx_0+1)
norm_nodes_0 = np.linspace(0,L_0,Nnodes) # Nodes
stag_nodes_0 = norm_nodes_0 + (norm_nodes_0[1]-norm_nodes_0[0])/2 # Staggered nodes
Nnodes = len(norm_nodes_0) # Length of nodes

# Find glacier width, surface elevation, and bed at nodes. I'm filtering these values so 
# that they have the same spatial scale.
filt_len = 2.0e3 # filter length in m
filt_cutoff=(1/filt_len)/(1/(np.diff(flowline[1:3,0])*2))
b,a=signal.butter(4,filt_cutoff,btype='low')
zs_filtered=signal.filtfilt(b,a,flowline[:,4]) # filtered surface elevation
zb_filtered=signal.filtfilt(b,a,flowline[:,3]) # filtered bed elevation
width_filtered=signal.filtfilt(b,a,glacierwidth[:,1]) # filtered width
u_filtered=signal.filtfilt(b,a,velocity[:,1]) #filtered velocity
del filt_len,filt_cutoff,b,a

# Interpolate ice thickness and velocity to nodes
u_0 = np.interp(norm_nodes_0,velocity[:,0],u_filtered)/yearinsec # Initial velocities
H_0 = np.interp(norm_nodes_0,flowline[:,0],zs_filtered-zb_filtered) # Initial heights

# Decide on a friction coefficient
# Beta2_0 is the friction coefficient for the first time step. We update the friction 
# coefficient (Beta2) on each time step to account for floating vs. grounding ice and 
# advance/retreat of the ice front.
temp = np.linspace(0,1,Nnodes)
Beta2_0 = 5.0e3
#Beta2_0 = 5.0e3*(np.exp(-0.5*temp)-np.exp(-0.5)) 
#Beta2_0 = (8.0e9)*np.linspace(1,0,Nnodes)
del temp

# Set up matrices for results
# Since we are using a moving grid, we need to keep track of where the nodes are located
# for each time step.
norm_nodes = np.zeros([len(norm_nodes_0)+pad,Nt]) # node locations for velocity, flux
stag_nodes = np.zeros([len(stag_nodes_0)+pad,Nt]) # node locations for ice thickness change

u = np.zeros([len(norm_nodes_0)+pad,Nt]) # velocity
u_new = np.zeros(Nnodes) # variable for velocity iterations
H_stag = np.zeros([len(stag_nodes_0)+pad,Nt]) # ice thickness on staggered grid
H = np.zeros([len(norm_nodes_0)+pad,Nt]) # ice thickness on normal nodes
zs = np.zeros([len(norm_nodes_0)+pad,Nt]) # surface elevation
L = np.zeros(Nt) # glacier length through time

# Set all variables to NaN
H[:,:] = 'NaN'
H_stag[:,:] = 'NaN'
zs[:,:] = 'NaN'
u[:,:] = 'NaN'

# Set up values for first time step
H[0:Nnodes,0] = H_0
H[0,:] = H_0[0]
u_new[0:Nnodes] = u_0
L[0] = L_0
zs[0:Nnodes,0] = np.interp(norm_nodes_0,flowline[:,0],zs_filtered)

##################
# Run time steps #
##################

CF_ind = np.argmin(abs(flowline[:,0]-L_0)) # Index of calving front
for i in range(0,Nt-1):
  if i % 100 == 0:
    print "Timestep", i, "out of", Nt 
  
  #########
  # Nodes #
  #########
  
  # Set up nodes for this time step
  Nnodes = np.ceil((L[i])/dx_0+1)
  norm_nodes[0:Nnodes,i] = np.linspace(0,L[i],Nnodes) # Nodes
  stag_nodes[0:Nnodes,i] = norm_nodes[0:Nnodes,i] + (norm_nodes[1,i]-norm_nodes[0,i])/2 # Staggered nodes
  dx = norm_nodes[1,i]-norm_nodes[0,i]
  CF_ind = Nnodes-1
  
  # Update ice thickness if we're past the first time step.
  if i > 0:
    H[1:Nnodes,i] = np.interp(norm_nodes[1:Nnodes,i],norm_nodes[1:len(u_new),i-1],dH+H[1:len(u_new),i-1])
  
  # Interpolate width and surface elevation
  W = np.interp(norm_nodes[0:Nnodes,i],glacierwidth[:,0],width_filtered)
  zb = np.interp(norm_nodes[0:Nnodes,i],flowline[:,0],zb_filtered)
  zs[0:Nnodes,i] = H[0:Nnodes,i]+zb
  
  #######################
  # Flotation criterion #
  #######################
  
  # Use flotation criterion to determine if glacier is floating
  # Get water depth at nodes
  D = -np.interp(norm_nodes[0:Nnodes,i],flowline[:,0],flowline[:,3])
  D[D<0] = 0
  
  # Bed float has values of -1 (floating ice), +1 (grounded ice), and 0 (grounding line)
  Bed_Float = np.zeros(len(norm_nodes[0:Nnodes,i]))
  Bed_Float[:] = 1 # Grounded ice
  ind = np.where(D > (rho_i/rho_sw)*H[0:Nnodes,i])
  Bed_Float[ind[0]] = -1 # Floating
  
  ##################
  # Basal friction #
  ##################
  
  # Set up boundary condition at bed. We want zero basal friction if the ice is floating.
  Beta2 = Beta2_0*np.ones(Nnodes)#np.interp(norm_nodes[0:Nnodes,i],norm_nodes_0,Beta2_0)
  Beta2[ind[0]]=0
  del ind
  
  # Set up variables for velocity iterations
  if i > 0:
    # Set velocity for this iteration to the velocity from the last time step.
    u_new = np.interp(norm_nodes[0:Nnodes,i],norm_nodes[0:len(u_new),i-1],u_new)
  else:
    # Set velocity equal to the measured velocity for the first time step.
    u_new = u_0
  u_last = np.zeros(Nnodes) # velocity on last iteration
  
  # Count the number of iterations
  n=0
  
  while (np.max(abs(u_last-u_new))*yearinsec > cutoff) and (n<100):
  
    #######################
    # Effective viscosity #
    #######################
    
    # Calculate effective viscosity
    nu = np.zeros(Nnodes)
    nu[0:-1] = A**(-1.0/3)*abs((u_new[1:]-u_new[0:-1])/dx)**(-2.0/3.0)
    nu[-1]=nu[-2] 

    # Set up coefficients for matrix problem
    Diagonal = np.zeros(Nnodes) # diagonal terms
    Left = np.zeros(Nnodes) # terms on left of diagonal
    Right = np.zeros(Nnodes) # terms on right side of diagonal
    Vector = np.zeros(Nnodes) # terms on right side of matrix problem
  
    # Find values at staggered nodes
    H_stag[0:Nnodes-1,i] = (H[1:Nnodes,i]+H[0:Nnodes-1,i])/2
    nu_stag = (nu[1:]+nu[0:-1])/2

    #######################################
    # Interior nodes behind calving front #
    #######################################
    
    Left[0:CF_ind-1] = (2.0/dx**2.0)*(H_stag[0:CF_ind-1,i]*nu_stag[0:CF_ind-1]) 
  
    # Gamma uses the velocities (u^(-2/3)) on the last time step to linearize the lateral
    # drag term.
    gamma = np.sign(u_new[1:CF_ind])*(abs(u_new[1:CF_ind]))**(-2.0/3.0)
    Diagonal[1:CF_ind] = -(2.0/dx**2.0)*(H_stag[1:CF_ind,i]*nu_stag[1:CF_ind]+H_stag[0:CF_ind-1,i]*nu_stag[0:CF_ind-1]) - \
      Beta2[1:CF_ind]*gamma*(H[1:CF_ind,i]-rho_w/rho_i*D[1:CF_ind])-gamma*(H[1:CF_ind,i]/W[1:CF_ind])*(5.0/(2.0*A*W[1:CF_ind]))**(1.0/3.0)     
  
    Right[2:CF_ind+1] = (2.0/dx**2.0)*(H_stag[1:CF_ind,i]*nu_stag[1:CF_ind])
  
    # Driving stress on right side of matrix problem
    Vector[1:CF_ind] = rho_i*g*H[1:CF_ind,i]*(zs[2:CF_ind+1,i]-zs[1:CF_ind,i])/dx
  
    ###################################
    # Nodes in front of calving front #
    ###################################
    
    # Nodes in front of calving front
    Left[CF_ind-1:]=-1.0
    Diagonal[CF_ind:] = 1.0
    Vector[CF_ind:] = A*dx*((rho_i*g*H[CF_ind:Nnodes,i])/4*(1-rho_i/rho_sw))**(3.0)
  
    ###################
    # Inflow boundary #
    ###################
    
    # Inflow boundary condition
    Diagonal[0] = 1.0 
    Vector[0] = u_inflow
    
    ##################
    # Matrix problem #
    ##################
    
    # Create sparse matrix
    Matrix = sparse.spdiags([Left,Diagonal,Right],[-1,0,1],Nnodes,Nnodes)

    # Save iterated velocities so that we can check to see if the solution has converged
    u_last = np.array(u_new)
    u_new = linalg.spsolve(Matrix,Vector)
    
    # Number of iterations
    n=n+1
  
  # The velocities have converged, so we save the converged solution.
  if (np.isnan(u_new[0])):
    print 
    sys.exit("Velocities became nan")
  
  u[0:Nnodes,i] = u_new
  
  ##########################
  # Find new calving front #
  ##########################
  
  R_xx = 2*((1/A)*(u[1:,i]-u[0:-1,i])/dx)**(1.0/3.0)
  
  d_crev = R_xx/(rho_i*g)+(rho_w/rho_i)*d_w
  
  # Find the glacier advance based on the velocity at the terminus.
  dL=u[Nnodes-1,i]*dt 
  L[i+1]=L[i]+dL
  
  ###########################
  # Change in ice thickness #
  ###########################
  
  # Find change in ice thickness
  dH = np.zeros(Nnodes)
  dH = -(1.0/W[0:-1])*(((u[1:Nnodes,i]*W[1:]*H[1:Nnodes,i])-(u[0:Nnodes-1,i]*W[0:Nnodes-1]*H[0:Nnodes-1,i]))/dx)*dt

  
