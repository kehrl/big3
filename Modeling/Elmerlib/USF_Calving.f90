! ******************************************************************************
! *
! *  Laura Kehrl
! *  University of Washington
! *  Dec. 1, 2014
! * 
! *****************************************************************************

! Initialize the calving front
FUNCTION XCalvingFrontIni ( Model, nodenumber, x) RESULT(CalvingFront)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   CalvingFront      
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True.

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'CalvingFront')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsCalvingIni','Could not find variable >Zs Calving<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % x (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % x (i)
          END IF
        END DO
   END IF

       CalvingFront = Zs0(ZsPerm(nodenumber)) 

END FUNCTION XCalvingFrontIni

! Find the calving front
FUNCTION Advance ( Model, nodenumber, dumy) RESULT(XTotDelta)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	Type(Mesh_t), TARGET :: OldMesh
   	TYPE(Mesh_t), POINTER :: Mesh
   	TYPE(Variable_t), POINTER :: XIniVariable, XVariable, XVelVariable, TimestepVariable, TimestepSizeVariable
   	REAL(KIND=dp), POINTER :: XVelValues(:), XIniValues(:), XValues(:)
   	INTEGER, POINTER :: BoundaryPerm(:), XVelPerm(:), XIniPerm(:), XPerm(:)
   	INTEGER, ALLOCATABLE :: boundIndex(:)
   	REAL(KIND=dp) :: minvelcalve, XTotDelta, XDelta, dt, CalvingFront, dumy, Retreated        
	INTEGER :: N, BoundaryNodes, minvelcalveNode, nodenumber,  i, dim, Timestep, LastTimestep
	LOGICAL :: FirstTime = .TRUE.

   	SAVE FirstTime 
	SAVE BoundaryPerm, BoundaryNodes, LastTimestep
	SAVE XVelPerm, XVelValues, XIniPerm, XIniValues, minvelcalveNode 
	SAVE XValues, XPerm
	SAVE Retreated

! Set up variables for first time
IF (FirstTime) then
   LastTimestep=0
   dim = CoordinateSystemDimension()
   
   Mesh => Model % Mesh
   N = Mesh % NumberOfNodes
   ALLOCATE(BoundaryPerm(N))
   BoundaryPerm = 0
   BoundaryNodes = 0
   CALL MakePermUsingMask(Model,Model%Solver,Model%Mesh,'CalvingMask',.FALSE.,BoundaryPerm,BoundaryNodes)

  ! Read in variables
  XVelVariable => VariableGet( Model % Variables, 'Velocity 1' )
  XIniVariable => VariableGet( Model % Variables, 'CalvingFrontIni')
  XVariable => VariableGet( Model % Variables, 'CalvingFront')
  IF ( ASSOCIATED( XIniVariable ) ) THEN
    XIniPerm    => XIniVariable % Perm
    XIniValues  => XIniVariable % Values
  END IF
  IF ( ASSOCIATED( XVariable ) ) THEN
    XPerm    => XVariable % Perm
    XValues  => XVariable % Values
  END IF
  IF ( ASSOCIATED( XVelVariable ) ) THEN
    XVelPerm    => XVelVariable % Perm
    XVelValues  => XVelVariable % Values	
  END IF

  FirstTime = .FALSE.
END IF

! Get time info
TimestepVariable => VariableGet( Model % Variables,'Timestep')
Timestep=TimestepVariable % Values(1)
TimestepSizeVariable => VariableGet( Model % Variables,'Timestep Size')
dt=TimestepSizeVariable % Values(1)

! Find new calving front offset based on minimum velocity at the calving front, 
! only do this once per timestep
IF( Timestep > LastTimestep) THEN
	! Find the minimum velocity at the calving front
  	minvelcalve = 2.0e4
  	DO i = 1,size(BoundaryPerm)
  	IF (BoundaryPerm(i) /=0 ) THEN  		
  		minvelcalve = min(minvelcalve, XVelValues(XVelPerm(i)))
    END IF    
  END DO

	! Find advance in terminus position
  	XDelta=dt*minvelcalve

  	! Update current terminus position
  	DO i = 1,size(BoundaryPerm)
  		IF (BoundaryPerm(i) /=0 ) THEN  		
  			XValues(XPerm(i)) = XDelta + XValues(XPerm(i))
    	END IF    
  	END DO
	
	! Compute mesh update
	XTotDelta= XValues(XPerm(nodenumber)) -  XIniValues(XIniPerm(nodenumber))
	
	LastTimestep = Timestep
ELSE
	! Compute mesh update
	XTotDelta = XValues(XPerm(nodenumber)) -  XIniValues(XIniPerm(nodenumber))
END IF

END FUNCTION Advance
