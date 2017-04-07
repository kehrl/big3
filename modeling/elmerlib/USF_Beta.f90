SUBROUTINE InterpolateBeta(Model,Solver,dt,Transient )
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
   

  TYPE(Variable_t), POINTER :: PointerToVariable, BetaVar
  TYPE(ValueList_t), POINTER :: Params
  
  CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='InterpolateBeta', FileName,BotMaskName
  
  LOGICAL :: FirstTime=.TRUE., UnfoundFatal=.TRUE.,Found

  INTEGER :: i, j, nx, ny, dummyint
  INTEGER, POINTER :: Permutation(:), BotPerm(:)=>NULL()

  REAL(KIND=dp) :: x,y
  REAL(KIND=dp), POINTER :: VariableValues(:)
  REAL(KIND=dp), ALLOCATABLE :: betas(:,:), xb(:), yb(:)
  REAL(kind=dp) :: LinearInterp
  
  SAVE betas, xb, yb, nx, ny
  SAVE FirstTime
  
  Params => GetSolverParams()  
  
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

  IF (FirstTime) THEN
    BotMaskName = "Bottom Surface Mask"
    
    CALL INFO(SolverName, 'Loading beta into memory for the first time')

    FirstTime=.False.

    FileName = ListGetString( Params, 'Data File', Found )

    ! open file
    OPEN(10,file=FileName)
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xb(nx),yb(ny))
    ALLOCATE(betas(nx,ny))
    DO i=1,nx
    	DO j=1,ny
      	READ(10,*) xb(i),yb(j),betas(i,j)
      END DO
		END DO
		CLOSE(10)
		
  END IF

  CALL INFO(SolverName, 'Interpolating beta onto mesh',level=3)

  ALLOCATE(BotPerm( Model % Mesh % NumberofNodes))
  BotPerm = 0
  CALL MakePermUsingMask( Model, Solver, Model % Mesh, BotMaskName, &
         .FALSE., BotPerm, dummyint)

  DO i=1,Model % Mesh % NumberofNodes
    IF ( BotPerm(i) <= 0) CYCLE
    x = Model % Mesh % Nodes % x(i)
    y = Model % Mesh % Nodes % y(i)
    VariableValues(Permutation(i)) = LinearInterp(betas,xb,yb,nx,ny,x,y)
  END DO

!------------------------------------------------------------------------------
END SUBROUTINE InterpolateBeta 
!------------------------------------------------------------------------------


!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!
