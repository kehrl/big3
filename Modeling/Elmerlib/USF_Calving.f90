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

FUNCTION XCalvingFront ( Model, nodenumber, dumy) RESULT(XTotDelta)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	TYPE(Mesh_t), POINTER :: Mesh
   	TYPE(Variable_t), POINTER :: XIniVariable, XVariable, XVelVariable, TimestepVariable, TimestepSizeVariable
   	REAL(KIND=dp), POINTER :: XVelValues(:), XIniValues(:), XValues(:)
   	INTEGER, POINTER :: BoundaryPerm(:), XVelPerm(:), XIniPerm(:), XPerm(:)
   	INTEGER, ALLOCATABLE :: boundIndex(:)
   	REAL(KIND=dp) :: Bed, XTotDelta, XDelta, dt, CalvingFront, dumy        
	INTEGER :: N, BoundaryNodes, BedNode, nodenumber,  i, dim, Timestep, LastTimestep
	LOGICAL :: FirstTime = .TRUE.
	INTEGER :: CoupledIter,PrevCoupledIter=-1

   	SAVE FirstTime 
	SAVE BoundaryPerm, BoundaryNodes, LastTimestep
	SAVE XVelPerm, XVelValues, XIniPerm, XIniValues, BedNode 
	SAVE XValues, XPerm

!********************************************************************************!
! ACCESS BOUNDARY ONLY ELEMENTS USING MakePermUsingMask
IF (FirstTime) then
   LastTimestep=0
   dim = CoordinateSystemDimension()
   
   Mesh => Model % Mesh
   N = Mesh % NumberOfNodes
   ALLOCATE(BoundaryPerm(N))
   BoundaryPerm = 0
   BoundaryNodes = 0
   CALL MakePermUsingMask(Model,Model%Solver,Model%Mesh,'CalvingMask',.FALSE.,BoundaryPerm,BoundaryNodes)

  ! READ IN ALL VALUES OF RELEVANT VARIABLES 
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

  !Find node at bed
  Bed=1.0e4
  DO i = 1,size(BoundaryPerm)
  	IF (BoundaryPerm(i) /=0 ) THEN  		
  		Bed = min(bed, Model % Nodes % y (i))
  		BedNode=i
    END IF    
  END DO
  FirstTime = .FALSE.
END IF

TimestepVariable => VariableGet( Model % Variables,'Timestep')
Timestep=TimestepVariable % Values(1)
TimestepSizeVariable => VariableGet( Model % Variables,'Timestep Size')
dt=TimestepSizeVariable % Values(1)

! Get velocity at bottom of calving front
CoupledIter = GetCoupledIter()
IF( Timestep > LastTimestep) THEN
  	XDelta=dt*XVelValues(XVelPerm(BedNode))
  	
  	DO i = 1,size(BoundaryPerm)
  		IF (BoundaryPerm(i) /=0 ) THEN  		
  			XValues(XPerm(i)) = XDelta + XValues(XPerm(i))
    	END IF    
  	END DO
  	!CalvingFront=CalvingFront+XDelta 	

    LastTimestep = Timestep
END IF

XTotDelta= XValues(XPerm(nodenumber)) -  XIniValues(XIniPerm(nodenumber))

END FUNCTION XCalvingFront
