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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION LateralConvergence (Model, nodenumber, dumy) RESULT(MassLoad)
   	USE Types
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	TYPE(Nodes_t) :: ElementNodes
   	TYPE(Variable_t), POINTER :: WeightVariable,VelocityVariable

   	INTEGER :: nodenumber, i, R, minind
   	INTEGER, POINTER :: WeightPerm(:), VelocityPerm(:)
   	
   	REAL(KIND=dp) :: x, y, velx
	REAL(KIND=dp), ALLOCATABLE :: xhw(:), hw(:), dwdx(:)
   	REAL(KIND=dp) :: dumy, MassLoad, weight, ratio, dist, mindist, width, dw
   	REAL(KIND=dp), POINTER :: WeightValues(:), VelocityValues(:)
   	
   	LOGICAL :: FirstTimeLC = .TRUE.
   	LOGICAL :: found = .FALSE.
	
	! Save inputs from file for future use
	SAVE xhw, hw, dwdx, FirstTimeLC
	
	if (FirstTimeLC) then

    	FirsttimeLC=.False.

        ! open file
        Open(10,file='Inputs/width.dat')
        Read(10,*) R
        Allocate(xhw(R),hw(R),dwdx(R))
        Do i=1,R
        	Read(10,*) xhw(i),hw(i),dwdx(i)
		End do
		Close(10)
    End if
	
	! Get current location
	x = Model % Nodes % x (nodenumber) 
	y = Model % Nodes % y (nodenumber)
	
	! Get element areas
	WeightVariable => VariableGet( Model % Variables, 'Flow Solution Weights' )
	IF (ASSOCIATED(WeightVariable)) THEN
    	WeightPerm    => WeightVariable % Perm
    	WeightValues  => WeightVariable % Values
    ELSE
        CALL FATAL('LateralConvergence','Could not find variable Weights. You need to add the Calculate Weights option.')
	END IF
	
	! Get velocity
	VelocityVariable => VariableGet( Model % Variables, 'Velocity 1' )
	IF (ASSOCIATED(VelocityVariable)) THEN
    	VelocityPerm    => VelocityVariable % Perm
    	VelocityValues  => VelocityVariable % Values
    ELSE
        CALL FATAL('LateralConvergence','Could not find x-velocity.')
	END IF
	
	! Get width at that location
	mindist=10000
   	DO i=1,SIZE(xhw)
      	dist=dabs(xhw(i)-x)
      	IF (dist<=mindist) THEN
        mindist=dist
        minind=i
        found = .true.
      	END IF
   	END DO

	! Interpolate width to current position
   	if (.not.found) then
      	print *, 'LateralConvergence: Could not find a suitable width to interpolate ',x
   	else
   		if (xhw(minind) >= x) then 
      		ratio=(x-xhw(minind))/(xhw(minind+1)-xhw(minind))
      		width=(hw(minind)+ratio*(hw(minind+1)-hw(minind)))
      		dw=(dwdx(minind)+ratio*(dwdx(minind+1)-dwdx(minind)))
      	else 
      		ratio=(x-xhw(minind))/(xhw(minind)-(xhw(minind-1)))
      		width=(hw(minind)+ratio*(hw(minind)-hw(minind-1)))
      		dw=(dwdx(minind)+ratio*(dwdx(minind)-dwdx(minind-1)))
      	endif
   	endif
	
	weight = WeightValues(WeightPerm(nodenumber))
	velx = VelocityValues(VelocityPerm(nodenumber))
	
	MassLoad = (-dw/width)*weight*velx	
	RETURN
END
FUNCTION LateralFrictionCoefficient( Model, nodenumber, dumy) RESULT(kcoeff) !
    USE types
	Use DefUtils
    implicit none
	TYPE(Model_t) :: Model
    Real(kind=dp) :: kcoeff, halfwidth, eta, dumy
    Real(kind=dp),allocatable :: xwidths(:),widths(:),junk(:)
    Real(kind=dp) :: x,y,mindist,dist,ratio,n,yearinsec,rhoi
    INTEGER :: nodenumber,minind,R,i
    TYPE(Variable_t), POINTER :: ViscosityVariable
    INTEGER, POINTER :: ViscosityPerm(:)
  	REAL(KIND=dp), POINTER :: ViscosityValues(:)
		
    LOGICAL :: FirsttimeLFC=.true.
    LOGICAL :: found

    SAVE xwidths, widths
    SAVE FirsttimeLFC
	
	if (FirsttimeLFC) then

    	FirsttimeLFC=.False.

        ! open file
        Open(10,file='Inputs/width.dat')
        Read(10,*) R
        Allocate(xwidths(R),widths(R),junk(R))
        Do i=1,R
        	Read(10,*) xwidths(i),widths(i),junk(i)
		End do
		Close(10)
    End if
    
    
    
    x = Model % Nodes % x (nodenumber)  
    
    ! Now find data closest to current position
   	found = .FALSE.		
	
   	mindist=100000
   	DO i=1,SIZE(xwidths)
      	dist=dabs(xwidths(i)-x)
      	IF (dist<=mindist) THEN
        mindist=dist
        minind=i
        found = .true.
      	END IF
   	END DO

	! Interpolate width to current position, and then find halfwidth for friction coefficient calculation
   	if (.not.found) then
      	print *, 'LateralFrictionCoefficient: Could not find a suitable width to interpolate ',x
   	else
   		if (xwidths(minind) >= x) then 
      		ratio=(x-xwidths(minind))/(xwidths(minind+1)-xwidths(minind))
      		halfwidth=(widths(minind)+ratio*(widths(minind+1)-widths(minind)))/2
      	else 
      		ratio=(x-xwidths(minind))/(xwidths(minind)-(xwidths(minind-1)))
      		halfwidth=(widths(minind)+ratio*(widths(minind)-widths(minind-1)))/2
      	endif
   	endif
    
	
	! Get viscosity at location
	ViscosityVariable => VariableGet( Model % Variables, 'Viscosity' )
   	IF (ASSOCIATED(ViscosityVariable)) THEN
    	ViscosityPerm    => ViscosityVariable % Perm
    	ViscosityValues  => ViscosityVariable % Values
    ELSE
        CALL FATAL('LateralFrictionCoefficient','Could not find variable Viscosity')
  	END IF
  	
  	eta = ViscosityValues(ViscosityPerm(nodenumber))
	n = 3
	yearinsec = 365.25*24*60*60
	rhoi = 917.0/(1.0e6*yearinsec**2)
	
	kcoeff = eta * (n+1)**(1/n) / (halfwidth**(1+(1/n)))
	
    Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION USF_Init (Model, nodenumber) RESULT(vel)
   	USE types
    USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber,ind, Rs, i
   	REAL(KIND=dp) :: x, dist, mindist, output, ratio, vel
   	REAL(KIND=dp), ALLOCATABLE :: xs(:), vs(:)
   	LOGICAL :: found
   
    LOGICAL :: Firsttime4=.true.
        
   	SAVE Firsttime4,xs,vs,Rs
   
   	if (Firsttime4) then
   		! Read file
   		Firsttime4=.False.
   		
   		OPEN(10,file="Inputs/velocity.dat")
   		Read(10,*) Rs
   		ALLOCATE(xs(Rs), vs(Rs))
   		READ(10,*)(xs(i), vs(i), i=1,Rs)
   		CLOSE(10)
   	End if
   
	x=Model % Nodes % x (nodenumber)
   
   	found = .false.
   	mindist=dabs(2*(xs(1)-xs(2)))
   
   	do 20, i=1,Rs
      	dist=dabs(x-xs(i))
      	if (dist<=mindist .and. xs(i)<=x) then
        	mindist=dist
        	ind=i
        	found = .true.
      	endif 
   	20 enddo

   	if (.not.found) then
      	print *, 'Could not find a suitable velocity to interpolate ',x
   	else
      	ratio=(x-xs(ind))/(xs(ind+1)-xs(ind))
      	vel=vs(ind)+ratio*(vs(ind+1)-vs(ind))
   	endif
   
   	Return
End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Viscosity( Model, nodenumber) RESULT(eta) !
    USE types
	Use DefUtils
    implicit none
	TYPE(Model_t) :: Model
    Real(kind=dp) :: T,eta,Aval,yearinsec,mindist,dist
    Real(kind=dp),allocatable :: xA(:),yA(:),A(:)
    Real(kind=dp) :: x,y,z,xnew
    INTEGER :: nodenumber,minind,RA,i
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'Viscosity'
		
    LOGICAL :: Firsttime5=.true.
    LOGICAL :: found

    SAVE xA,yA,A,RA
    SAVE Firsttime5

    if (Firsttime5) then

    	Firsttime5=.False.

        ! open file
        Open(10,file='Inputs/flowparameters.dat')
        Read(10,*) RA
        Allocate(xA(RA),yA(RA),A(RA))
        Do i=1,RA
        	Read(10,*) xA(i),yA(i),A(i)
		End do
		Close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

	! Now find data closest to current position
   	found = .FALSE.		
	
	IF (x > xA(RA)) THEN
		xnew = XA(RA)
		mindist=10000
   	ELSE
   		xnew = x
   		mindist=5000
   	END IF 
   	
   	DO i=1,RA
      	dist=sqrt((xnew-xA(i))**2+(y-yA(i))**2)
      	IF (dist<=mindist) THEN
        mindist=dist
        minind=i
        found = .true.
      	END IF
   	END DO
	
	IF (found) THEN
		yearinsec=365.25*24*60*60
    	Aval=A(minind)
		eta=(2.0*Aval*yearinsec)**(-1.0_dp/3.0_dp)*1.0e-6
		!print *,'',x,y,eta
	else
		print *,'No viscosity at',x,y
		CALL FATAL(SolverName,'No viscosity found for above coordinates')
	end if	
    !print *,'Viscosity is',x,y,eta
    Return
End

