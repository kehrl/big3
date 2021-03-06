!Get slip coefficient
!--------------------------------------------------------------------------!
FUNCTION Linear_NearestNeighbor( Model, nodenumber, dumy) RESULT(coeff) !
!------------------------------------------------------------------!
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber, nx, i, ind(1)
   REAL(KIND=dp) :: x, y, dumy, coeff
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), betas(:)
   LOGICAL :: FirstTime=.True.
   CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_linear.dat'
   REAL(kind=dp) :: LinearInterp

   SAVE betas, xb, yb, nx
   SAVE FirstTime

   IF (FirstTime) THEN
      
      FirstTime=.False.

      ! open file
      OPEN(10,file='inputs/beta_linear.dat')
      READ(10,*) nx
      ALLOCATE(xb(nx),yb(nx),betas(nx))
      DO i=1,nx
          READ(10,*) xb(i),yb(i),betas(i)
      END DO
      CLOSE(10)
  
   END IF

   ! position current point
   x = Model % Mesh % Nodes % x (nodenumber)
   y = Model % Mesh % Nodes % y (nodenumber)

   ! Get nearest node
   ind = minloc(sqrt((x-xb)**2+(y-yb)**2),1)
   coeff = betas(ind(1))

   Return
End


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
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

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
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)
  
  coeff = LinearInterp(betas,xb,yb,nx,ny,x,y)

  Return
End

!Get slip coefficient
!--------------------------------------------------------------------------!
FUNCTION Weertman_NearestNeighbor( Model, nodenumber, dumy) RESULT(coeff) !
!------------------------------------------------------------------!
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber, nx, i, ind(1)
   REAL(KIND=dp) :: x, y, dumy, coeff
   REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), betas(:)
   LOGICAL :: FirstTime=.True.
   CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_weertman.dat'
   REAL(kind=dp) :: LinearInterp

   SAVE betas, xb, yb, nx
   SAVE FirstTime

   IF (FirstTime) THEN
      
      FirstTime=.False.

      ! open file
      OPEN(10,file='inputs/beta_weertman.dat')
      READ(10,*) nx
      ALLOCATE(xb(nx),yb(nx),betas(nx))
      DO i=1,nx
          READ(10,*) xb(i),yb(i),betas(i)
      END DO
      CLOSE(10)
  
   END IF

   ! position current point
   x = Model % Mesh % Nodes % x (nodenumber)
   y = Model % Mesh % Nodes % y (nodenumber)

   ! Get nearest node
   ind = minloc(sqrt((x-xb)**2+(y-yb)**2),1)
   coeff = betas(ind(1))

   Return
End

!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!
