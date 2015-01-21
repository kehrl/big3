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