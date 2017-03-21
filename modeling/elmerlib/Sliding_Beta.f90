!Get slip coefficient
!--------------------------------------------------------------------------!
FUNCTION Linear( Model, nodenumber, dumy) RESULT(coeff) !
    !------------------------------------------------------------------!
	USE types
	USE CoordinateSystems
	USE SolverUtils
	USE ElementDescription
	USE DefUtils
	IMPLICIT NONE
	TYPE(Model_t) :: Model
	TYPE(Solver_t), TARGET :: Solver
	INTEGER :: nodenumber, nx, ny, i, j
	REAL(KIND=dp) :: x, y, dumy, coeff 
	REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), betas(:,:)
	LOGICAL :: FirstTime=.True.
	CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_linear.xy'
	REAL(kind=dp) :: LinearInterp

  SAVE betas, xb, yb, nx, ny
  SAVE FirstTime

  IF (FirstTime) THEN

    FirstTime=.False.

    ! open file
    OPEN(10,file='inputs/beta_linear.xy')
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

  ! position current point
  x = Model % Nodes % x (nodenumber)
  y = Model % Nodes % y (nodenumber)

  coeff = LinearInterp(betas,xb,yb,nx,ny,x,y)

  Return
End

FUNCTION Weertman( Model, nodenumber, dumy) RESULT(coeff) !
    !------------------------------------------------------------------!
	USE types
	USE CoordinateSystems
	USE SolverUtils
	USE ElementDescription
	USE DefUtils
	IMPLICIT NONE
	TYPE(Model_t) :: Model
	TYPE(Solver_t), TARGET :: Solver
	INTEGER :: nodenumber, i, j, nx, ny
	REAL(KIND=dp) :: x, y, dumy, coeff
	REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), betas(:,:)
	LOGICAL :: FirstTime=.True.
	CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_weertman.xy'
	REAL(kind=dp) :: LinearInterp
	
  SAVE betas, xb, yb, nx, ny
  SAVE FirstTime

  IF (FirstTime) THEN

    FirstTime=.False.

    ! open file
    OPEN(10,file='inputs/beta_weertman.xy')
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

  ! position current point
  x = Model % Nodes % x (nodenumber)
  y = Model % Nodes % y (nodenumber)

  coeff = LinearInterp(betas,xb,yb,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!