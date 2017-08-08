SUBROUTINE InterpolateAccumulation(Model,Solver,dt,Transient )
  !------------------------------------------------------------------------------
  USE CRSMatrix
  USE GeneralUtils
  USE ElementDescription 
  USE MainUtils
  USE SolverUtils

  IMPLICIT NONE
  !------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: Transient

  !------------------------------------------------------------------------------
  ! Local variables
  !------------------------------------------------------------------------------
   

  TYPE(Variable_t), POINTER :: PointerToVariable, bdotVar, TimeStepVar
  TYPE(ValueList_t), POINTER :: Params
  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='InterpolateAccumulation', FileName,TopMaskName
  
  LOGICAL :: FirstTime=.TRUE., UnfoundFatal=.TRUE.,Found, Steady

  INTEGER :: i, j, nx, ny, dummyint, TimeStep, TimeStepLast, restartposition
  INTEGER, POINTER :: Permutation(:), TopPerm(:)=>NULL()

  REAL(KIND=dp) :: x,y
  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp), ALLOCATABLE :: smbgrid(:,:), xx(:), yy(:)
  REAL(kind=dp) :: LinearInterp
  
  SAVE smbgrid, xx, yy, nx, ny
  SAVE FirstTime, TimeStepLast, restartposition, Steady
  
  Params => GetSolverParams()  

  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values
  
  ! Find out if this run is a restart of a previous run, and if so add that to the timestep
  ! to get the right accumulation rate
  IF (FirstTime) THEN
      restartposition = ListGetInteger(Params,'Restart Timestep',Found)
      IF (.NOT. Found) THEN
        CALL FATAL(SolverName, 'No restart position given.')
      END IF
      Steady = ListGetLogical(Params,'Steady State',Found)
      IF (.NOT. Found) THEN
        CALL INFO(SolverName, 'Steady State not given, assuming it is false.')
      END IF
      print *,'STEADY',Steady
  
      TopMaskName = "Top Surface Mask"
  END IF
  
  IF (Steady) THEN
    
    IF (FirstTime) THEN
      FirstTime = .FALSE.

      ! open file
      OPEN(10,file="inputs/smb.xy")
      READ(10,*) nx
      READ(10,*) ny
      ALLOCATE(xx(nx),yy(ny))
      ALLOCATE(smbgrid(nx,ny))
      DO i=1,nx
        DO j=1,ny
          READ(10,*) xx(i),yy(j),smbgrid(i,j)
        END DO
		  END DO
		  CLOSE(10)
    END IF
  ELSE

    TimestepVar => VariableGet( Model % Variables,'Timestep')
    Timestep = TimestepVar % Values(1) + restartposition

    IF ((FirstTime) .OR. (TimeStep /= TimeStepLast)) THEN
			
			  IF (.NOT. FirstTime) THEN
			    DEALLOCATE(xx,yy,smbgrid)
			  END IF
    	  FirstTime=.False.

        WRITE (FileName, "(A10,I4.4,A3)") "inputs/smb", Timestep,".xy"

        ! open file
        OPEN(10,file=FileName)
        READ(10,*) nx
        READ(10,*) ny
        ALLOCATE(xx(nx),yy(ny))
        ALLOCATE(smbgrid(nx,ny))
        DO i=1,nx
          DO j=1,ny
        	  READ(10,*) xx(i),yy(j),smbgrid(i,j)
          END DO
			  END DO
			  CLOSE(10)
    
        TimestepLast = Timestep
    
    END IF
    
  END IF

  CALL INFO(SolverName, 'Interpolating bdot onto mesh',level=3)

  ALLOCATE(TopPerm( Model % Mesh % NumberofNodes))
  TopPerm = 0
  CALL MakePermUsingMask( Model, Solver, Model % Mesh, TopMaskName, &
         .FALSE., TopPerm, dummyint)

  DO i=1,Model % Mesh % NumberofNodes
    IF ( TopPerm(i) <= 0) CYCLE
    x = Model % Mesh % Nodes % x(i)
    y = Model % Mesh % Nodes % y(i)
    VariableValues(Permutation(i)) = LinearInterp(smbgrid,xx,yy,nx,ny,x,y)
  END DO
  
  CALL INFO(SolverName, 'End', level=3)

!------------------------------------------------------------------------------
END SUBROUTINE InterpolateAccumulation 
!------------------------------------------------------------------------------


!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!
