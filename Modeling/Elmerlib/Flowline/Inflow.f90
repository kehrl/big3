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
   	LOGICAL :: Firsttime1=.true.
   	LOGICAL :: Firsttime2=.true.
        
   	SAVE Firsttime1,Firsttime2,dv,us,df,zb,zs
	
   	if (Firsttime1) then
   	! Read velocity file
   		Firsttime1=.False.
   	
   		OPEN(10,file="Inputs/velocity.dat")
 		Read(10,*) R
   		READ(10,*) dv, us
   		CLOSE(10)
   	End if

   	if (Firsttime2) then
   		! Read velocity file
   		! Read flowline so that we can find thickness
   		Firsttime2=.False.
   		
   		OPEN(10,file="Inputs/flowline.dat")
   		Read(10,*) R
   		READ(10,*) df, xf, yf, zb, zs 
   		CLOSE(10)
   	End if
   
   	x=Model % Nodes % x (nodenumber)
   	y=Model % Nodes % y (nodenumber)
   

    thick=zs-zb
    ub=60
    vel = ub + (1.0_dp - ((zs - y) / (thick))**4) * (us-ub)

   	Return 
End