!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculates the velocity at the inflow boundary in the flowline model, using 
!! the SIA approximation. Currently assumes a basal velocity of 60 m/a at inflow (may
!! want to change when everything else is fixed).
!!
!! LMK, UW, 09/10/2014
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Function Inflow (Model, nodenumber, dumy) RESULT(vel)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber
   	REAL(KIND=dp) :: Mask
   	INTEGER :: DIM, R, i, ind
   	REAL(KIND=dp) :: x, y, vel, thick, mindist, dist, dumy
   	REAL(KIND=dp) :: xf, yf, df, zb, zs, dv, us, ub
   	
   	LOGICAL :: found
   	LOGICAL :: FirsttimeInflow=.True.
     
    ! Save data for future use    
   	SAVE FirsttimeInflow,dv,us,df,zb,zs
	
	! Load data
   	if (FirsttimeInflow) then

   		FirsttimeInflow=.False.
   	   	
   	   	! Read velocity file
   		OPEN(10,file="Inputs/velocity.dat")
 		Read(10,*) R
 		! Only need to read the first row for inflow
   		READ(10,*) dv, us
   		CLOSE(10)

   		! Read flowline so we can figure out thickness
   		OPEN(10,file="Inputs/flowline.dat")
   		Read(10,*) R
   		! Only need to read the first row for inflow
   		READ(10,*) df, xf, yf, zb, zs 
   		CLOSE(10)
   		
   	End if
   	
    ! Get current position
   	x=Model % Nodes % x (nodenumber)
   	y=Model % Nodes % y (nodenumber)
   
    ! Compute velocity according to an approximation to the SIA, assuming n=3
    thick=zs-zb
    ub=200.0_dp
    vel = ub + (1.0_dp - ((zs - y) / (thick))**4_dp) * (us-ub)
   	
   	Return 
End