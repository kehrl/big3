# This is a code to solve the SSA. 

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as linalg
import scipy.signal as signal
import sys
import matplotlib.pyplot as plt

###########
# Options #
###########

CF_move = "False"

Helheim = "True"

# Numerical method
# Options: Substitution, Newton
#Method = "Substitution"
Method = "Newton"

#############
# Constants #
#############

# Constants
rho_i = 917.0 # ice density (kg/m3)
rho_sw = 1028.0 # sea water density (kg/m3)
rho_w = 1000.0
g = 9.8 # gravity (m/s2)
yearinsec = 365.25*24*60*60

# Flow law exponent
n = 3.0

# Flow law parameter
A = 1.7e-24 # Pa-3 s-1

# Sliding law exponent
m = 3.0 # velocity exponent
p = 1.0 # water pressure / height above buoyancy exponent

# Calving law parameters
d_w = 10.0 # water depth in crevasses (m)

# Accumulation rate
ELA = 600 # equilibrium line altitude in m
mdot = 0/yearinsec # m across the entire width
adot = 0/yearinsec # m across the entire width

##################
# Model Numerics #
##################

# Cutoff for nonlinear iterations
cutoff = 1.0 # change between iterations in m/yr
iter_max = 1.0e3 # maximum number of iteration, need to add to code

# Relaxation parameter for method Substitution
Relax = 0.5

##################
# Model Geometry #
##################

# Choose node spacing
dx_0 = 100.0 # desired grid size in meters

# Set initial length of flowband
if Helheim=="True":
  L_0 = 85115.0 
else:
  L_0 = 10000

# Number of nodes
N = np.ceil((L_0)/dx_0+1) 

# Load info about glacier geometry
DIR="/Users/kehrl/Models/Helheim/Meshes/Flowline/MorlighemNew_SmoothedVelocity/Inputs/"
flowline = np.loadtxt(DIR+"flowline.dat",skiprows=1)
velocity = np.loadtxt(DIR+"velocity.dat",skiprows=1)
glacierwidth = np.loadtxt(DIR+"width.dat",skiprows=1)

if Helheim=="True":
  input_x = flowline[:,0]
  input_u = velocity[:,1]
  input_zs = flowline[:,4]
  input_zb = flowline[:,3]
  input_W = np.interp(input_x,glacierwidth[:,0],glacierwidth[:,1])
  input_H = input_zs - input_zb
else:
  input_x = np.linspace(0,L_0+1.0e4,N)
  input_u = np.linspace(500.0,1000,len(input_x))
  input_zs = 75.0*np.ones(len(input_x))
  input_zb = -2000*np.ones(len(input_x))
  input_H = 700.0*np.ones(len(input_x))
  
del flowline, velocity, glacierwidth

#############
# Time step #
#############

# Choose time step and length of run
dt = 1.0/(365.0)*yearinsec # Years
t = 1.0*yearinsec # Years
Nt = int(np.ceil(t/dt)) # number of timesteps

#############################
# Inflow Boundary Condition #
#############################

# Set inflow boundary condition (Dirichlet BC)
u_inflow = input_u[0]/yearinsec

# In case the glacier advances, we need to add some extra same in our matrices so that 
# we don't run out of room.
pad = np.ceil((2.0e4*t/yearinsec)/dx_0) # number of nodes

###############
# Model setup #
###############

# Set up initial nodes and values for velocity and height. We're using a staggered grid 
# for numerical stability. We calculate height on the normal grid and velocity/flux on the
# staggered grid. The grounding line is located on the staggered grid.

# Nodes on normal grid
norm_nodes_0 = np.linspace(0,L_0,N) 

# Nodes on staggered grid
stag_nodes_0 = norm_nodes_0 + (norm_nodes_0[1]-norm_nodes_0[0])/2 

# Interpolate ice thickness, velocity, and surface geometry to nodes
u_0_stag = np.interp(stag_nodes_0,input_x,input_u)/yearinsec # Initial velocities
H_0 = np.interp(norm_nodes_0,input_x,input_H) # Initial heights

# Decide on a friction coefficient
# Beta2_0 is the friction coefficient for the first time step. We update the friction 
# coefficient (Beta2) on each time step to account for floating vs. grounding ice and 
# advance/retreat of the ice front.
#temp = np.linspace(0,1,N)
Beta2_0 = 8.0e3
#Beta2_0 = 5.0e3*np.interp(norm_nodes_0,input_x,input_W)/np.max(input_W)
#Beta2_0 = 5.0e3*(np.exp(-0.5*temp)-np.exp(-0.5)) 
#Beta2_0 = (5.0e3)*np.linspace(1,0,N)
#del temp

# Set up matrices for results
# Since we are using a moving grid, we need to keep track of where the nodes are located
# for each time step.
norm_nodes = np.zeros([len(norm_nodes_0)+pad,Nt]) # node locations for velocity, flux
stag_nodes = np.zeros([len(stag_nodes_0)+pad,Nt]) # node locations for ice thickness change

u_stag = np.zeros([len(stag_nodes_0)+pad,Nt]) # velocity
u_new = np.zeros(N) # variable for velocity iterations
H = np.zeros([len(norm_nodes_0)+pad,Nt]) # ice thickness on normal nodes
zs = np.zeros([len(norm_nodes_0)+pad,Nt]) # surface elevation
L = np.zeros(Nt) # glacier length through time
GL = np.zeros(Nt)

# Set all variables to NaN
H[:,:] = 'NaN'
zs[:,:] = 'NaN'
u_stag[:,:] = 'NaN'

# Set up values for first time step
H[0:N,0] = H_0
H[0,:] = H_0[0]
u_new[0:N] = u_0_stag
L[0] = L_0
zs[0:N,0] = np.interp(norm_nodes_0,input_x,input_zs)

##################
# Run time steps #
##################

for i in range(0,Nt-1):
  if i % 100 == 0:
    print "Timestep", i, "out of", Nt 
  
  #########
  # Nodes #
  #########
  
  # Set up nodes for this time step
  N = np.ceil((L[i])/dx_0+1)
  norm_nodes[0:N,i] = np.linspace(0,L[i],N) # Nodes
  stag_nodes[0:N,i] = norm_nodes[0:N,i] + (norm_nodes[1,i]-norm_nodes[0,i])/2 # Staggered nodes
  dx = norm_nodes[1,i]-norm_nodes[0,i]
  CF_ind = N-1
  
  # Update ice thickness if we're past the first time step.
  if i > 0:
    H[1:N,i] = np.interp(norm_nodes[1:N,i],norm_nodes[1:len(u_new),i-1],dH+H[1:len(u_new),i-1])
  
  # Interpolate width and surface elevation
  W = np.interp(norm_nodes[0:N,i],input_x,input_W)
  W_stag = np.zeros(N)
  W_stag[0:-1] = (W[0:-1]+W[1:])/2
  W_stag[-1] = np.interp(stag_nodes[N-1,i],input_x,input_W)
  zb = np.interp(norm_nodes[0:N,i],input_x,input_zb)
  
  #######################
  # Flotation criterion #
  #######################
  
  # Use flotation criterion to determine if glacier is floating
  # Get water depth at nodes
  D = -np.interp(norm_nodes[0:N,i],input_x,input_zb)
  D[D<0] = 0
  
  # Find values at staggered nodes
  H_stag = (H[1:N,i]+H[0:N-1,i])/2
  D_stag = (D[1:N]+D[0:N-1])/2
  
  # Bed float has values of -1 (floating ice), +1 (grounded ice), and 0 (grounding line)
  Bed_Float_norm = np.zeros(len(norm_nodes[0:N,i]))
  Bed_Float_stag = np.zeros(len(norm_nodes[0:N,i]))
  Bed_Float_norm[:] = 1 # Grounded ice
  Bed_Float_stag[:] = 1
  ind_norm = np.where(rho_i*H[0:N,i] < rho_sw*D)
  ind_stag = ind_norm[0]-1
  
  if np.any(ind_norm[0]):
    Bed_Float_norm[ind_norm[0][0]:] = -1 # Floating
    Bed_Float_stag[ind_stag[0]:] = -1 # Floating
    FL_ind_norm = ind_norm[0][0] # First floating index on normal grid
    if FL_ind_norm == -1:
      FL_ind_norm = 0
    GL[i] = stag_nodes[FL_ind_norm-1,i]
    zs[0:FL_ind_norm,i] = H[0:FL_ind_norm,i]+zb[0:FL_ind_norm]
    zs[FL_ind_norm:,i] = H[FL_ind_norm:,i]*(1-rho_i/rho_w)
    if (i == 0): 
      print "Grounding line at", GL[i]
    elif abs(GL[i]-GL[i-1]) > 20:
      print "Grounding Line moved to", GL[i]
  else:
    GL[i] = L[i]
    zs[0:N,i] = H[0:N,i]+zb

  
  ##################
  # Basal friction #
  ##################
  
  # Set up boundary condition at bed. We want zero basal friction if the ice is floating.
  Beta2 = Beta2_0*np.ones(N)#np.interp(norm_nodes[0:N,i],norm_nodes_0,Beta2_0)
  Beta2[Bed_Float_stag==-1]=0
  
  # Set up variables for velocity iterations
  if i > 0:
    # Set velocity for this iteration to the velocity from the last time step.
    u_new = np.interp(norm_nodes[0:N,i],norm_nodes[0:len(u_new),i-1],u_new)
  else:
    # Set velocity equal to the measured velocity for the first time step.
    u_new = u_0_stag
  u_last = np.zeros(N) # velocity on last iteration
  
  # Count the number of iterations
  iter = 0
  
  while np.max(abs(u_last-u_new))*yearinsec > cutoff and (iter < iter_max):
  
    #######################
    # Effective viscosity #
    #######################
    
    # Calculate effective viscosity
    nu = np.zeros(N)
    nu[0:-1] = A**(-1.0/n)*abs((u_new[1:]-u_new[0:-1])/dx)**((1-n)/n)
    nu[-1]=nu[-2] 
    
    if Method == "Substitution":
    
      # Set up coefficients for matrix problem
      Diagonal = np.zeros(N) # diagonal terms
      Left = np.zeros(N) # terms on left of diagonal
      Right = np.zeros(N) # terms on right side of diagonal
      Vector = np.zeros(N) # terms on right side of matrix problem
    
      #######################################
      # Interior nodes behind calving front #
      #######################################
    
      Left[0:CF_ind-1] = (2.0/dx**2.0)*(H[1:CF_ind,i]*nu[1:CF_ind]) 
  
      # Gamma uses the velocities (u^((1-n)/n)) on the last time step to linearize the lateral
      # drag term.
    
      gamma_n = (u_new[1:CF_ind])**((1-n)/n)
      gamma_m = (u_new[1:CF_ind])**((1-m)/m)
      
      Diagonal[1:CF_ind] = -(2.0/dx**2.0)*(H[1:CF_ind,i]*nu[1:CF_ind]+H[2:CF_ind+1,i]*nu[2:CF_ind+1]) - \
          gamma_m*Beta2[1:CF_ind]*(abs(H_stag[1:CF_ind]-(rho_sw/rho_i)*D_stag[1:CF_ind]))**(1.0/p)- \
          gamma_n*(H_stag[1:CF_ind]/W_stag[1:CF_ind])*(5.0/(2.0*A*W_stag[1:CF_ind]))**(1.0/n)     
      
      Right[2:CF_ind+1] = (2.0/dx**2.0)*(H[2:CF_ind+1,i]*nu[2:CF_ind+1])
  
      # Driving stress on right side of matrix problem
      Vector[1:CF_ind] = rho_i*g*H_stag[1:CF_ind]*(zs[2:CF_ind+1,i]-zs[1:CF_ind,i])/dx
  
      ###################################
      # Nodes in front of calving front #
      ###################################
    
      # Nodes in front of calving front
      Left[CF_ind-1:]=-1.0
      Diagonal[CF_ind:] = 1.0
      Vector[CF_ind:] = A*dx*((rho_i*g*H[CF_ind:N,i])/4*(1-rho_i/rho_sw))**(n)
  
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
      Matrix = sparse.spdiags([Left,Diagonal,Right],[-1,0,1],N,N)
      
      # Save iterated velocities so that we can check to see if the solution has converged
      u_last = np.array(u_new)
      u_new = Relax*linalg.spsolve(Matrix,Vector)+(1-Relax)*u_last

    elif Method == "Newton":
      
      # Set up variables
      F = np.zeros(N) # SSA equation = 0
      dFdu_m1 = np.zeros(N) # u_(i-1/2) part of Jacobian
      dFdu_0  = np.zeros(N) # u_(i+1/2) part of Jacobian
      dFdu_p1 = np.zeros(N) # u_(i+3/2) part of Jacobian
      
      ###################
      # Inflow boundary #
      ###################
      
      F[0] = u_new[0]-u_inflow
      dFdu_0[0] = 1
      
      ##################
      # Interior Nodes #
      ##################
      
      # Left hand side of equation
      F[1:-1] = (2.0/dx**2.0)*(H[2:N,i]*nu[2:N]*(u_new[2:N]-u_new[1:N-1])-H[1:N-1,i]*nu[1:N-1]*(u_new[1:N-1]-u_new[0:N-2])) - \
          Beta2[1:N-1]*((H_stag[1:N-1]-(rho_sw/rho_i)*D_stag[1:N-1]))**(1.0/p)*u_new[1:N-1]**(1.0/m) - \
          (H_stag[1:N-1]/W_stag[1:N-1])*((5*u_new[1:N-1])/(2*A*W_stag[1:N-1]))**(1.0/n)- \
          rho_i*g*H_stag[1:N-1]*(zs[2:N,i]-zs[1:N-1,i])/dx
    
      # For u_(i-1/2)
      dFdu_m1[0:CF_ind-1] = (2.0/dx**2.0)*(H[1:CF_ind,i]*nu[1:CF_ind]) 
      
      # For u_(i+1/2)
      gamma_n = (u_new[1:CF_ind])**((1-n)/n)
      gamma_m = (u_new[1:CF_ind])**((1-m)/m)
      dFdu_0[1:CF_ind] = -(2.0/dx**2.0)*(H[1:CF_ind,i]*nu[1:CF_ind]+H[2:CF_ind+1,i]*nu[2:CF_ind+1]) - \
          Beta2[1:N-1]*((H_stag[1:N-1]-rho_sw/rho_i*D_stag[1:N-1]))**(1.0/p)*(1.0/m)*gamma_m - \
          H_stag[1:N-1]/W_stag[1:N-1]*((5/(2*A*W_stag[1:N-1]))**(1.0/n))*(1.0/n)*gamma_n
      
      # For u_(i+3/2)
      dFdu_p1[2:CF_ind+1] = (2.0/dx**2.0)*(H[2:CF_ind+1,i]*nu[2:CF_ind+1])
    
      #################
      # Calving front #
      #################
      
      F[-1] = (u_new[-1]-u_new[-2])/dx-A*(rho_i*g/4*H[N-1,i]*(1-rho_i/rho_sw))**n
      
      dFdu_m1[-2] = -1/dx
      dFdu_0[-1] = 1/dx
      
      ################
      # Solve for du #
      ################
      
      Matrix = sparse.spdiags([dFdu_m1,dFdu_0,dFdu_p1],[-1,0,1],N,N)
      
      du = linalg.spsolve(Matrix,-F)
      u_last = np.array(u_new)
      u_new = u_last + Relax*du
      
    
    # Number of iterations
    iter = iter+1
    
    # Check to see if we are on last iteration and if so print an error message stating 
    # that the solution never converged.
    if iter == iter_max:
      print "Nonlinear iteration did not converge for timestep",i,np.max(abs(u_last-u_new))*yearinsec
  
  # If for some reason the velocities become NaN, we want to exit the model.
  if (np.isnan(u_new[0])):
    sys.exit("Velocities became nan")
  
  # The velocities have converged, so we save the converged solution.
  u_stag[0:N,i] = u_new
  
  #if i == 0:
  #  sys.exit("blah")
  
  ##########################
  # Find new calving front #
  ##########################
  
  # Only move the terminus if the simulation permits a moving terminus
  if CF_move=="True":
  
    # To determine the position of the new calving front, we use a crevasse depth calving law.
  
    # Find the longitudinal resistive stress
    R_xx = 2*((1/A)*(abs(u_stag[1:N,i]-u_stag[0:N-1,i])/dx))**(1.0/3.0)
  
    # Find the crevasse depth
    d_crev = R_xx/(rho_i*g)+(rho_w/rho_i)*d_w
  
    # Check to see if the calving condition is met
    ind = np.where(d_crev[0:N-1] > zs[0:N-1,i])
    indices = np.where(Bed_Float[ind[0]]==-1)
    if np.any(indices): 
      CF_ind = np.min(ind[0][indices[0]])
      L[i+1] = norm_nodes[CF_ind,i]
    
    else:
      # Find the glacier advance based on the velocity at the terminus.
      dL=u_stag[N-1,i]*dt 
      L[i+1]=L[i]+dL
      
  # Otherwise maintain the same terminus position  
  else: 
    L[i+1]=L_0
  
  ###########################
  # Change in ice thickness #
  ###########################
  
  # Interpolate mass balance
  Bdot = np.zeros(N)
  ELA_ind = np.argmin(abs(zs[0:N,i]-ELA))
  Bdot[0:ELA_ind]=adot
  Bdot[ELA_ind+1:]=mdot
  
  # Find change in ice thickness
  dH = np.zeros(N-1)
  dH = -(1.0/W[1:N])*(((u_stag[1:N,i]*W_stag[1:N]*H[1:N,i])- \
       (u_stag[0:N-1,i]*W_stag[0:N-1]*H[0:N-1,i]))/dx)*dt+Bdot[1:N]*dt

  
