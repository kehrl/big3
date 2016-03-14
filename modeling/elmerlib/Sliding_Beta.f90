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
	INTEGER :: nodenumber, Nb, i, ind(1)
	REAL(KIND=dp) :: x, y, z, dumy, coeff, beta 
	REAL(KIND=dp) :: Zbed, znode, dist(1)
	REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:), betas(:)
	LOGICAL :: FirstTime=.True.
	CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_linear.xyz'

	SAVE xb, yb, zb, betas
	SAVE Firsttime

	if (Firsttime) then

    	Firsttime=.False.

       	! open file
       	open(10,file='inputs/beta_linear.xyz')
    	Read(10,*) Nb
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb), betas(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), betas(i), i=1,Nb)
		close(10)
				
	End if

	! position current point
    x=Model % Nodes % x (nodenumber)
    y=Model % Nodes % y (nodenumber)
    ! Find the height of the current point
    z=Model % Nodes % z (nodenumber)
		
	ind=minloc(sqrt((x-xb)**2+(y-yb)**2),1)
    dist=minval(sqrt((x-xb)**2+(y-yb)**2),1)
    beta=betas(ind(1))
    coeff=beta

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
	INTEGER :: nodenumber, Nb, i, ind(1)
	REAL(KIND=dp) :: x, y, z, dumy, coeff, beta 
	REAL(KIND=dp) :: Zbed, znode, dist(1)
	REAL(KIND=dp), ALLOCATABLE :: xb(:), yb(:), zb(:), betas(:)
	LOGICAL :: FirstTime=.True.
	CHARACTER(len=MAX_NAME_LEN) :: filin='inputs/beta_weertman.xyz'

	SAVE xb, yb, zb, betas
	SAVE Firsttime

	if (Firsttime) then

    	Firsttime=.False.

       	! open file
       	open(10,file='inputs/beta_weertman.xyz')
    	Read(10,*) Nb
        ALLOCATE(xb(Nb), yb(Nb), zb(Nb), betas(Nb))
        READ(10,*)(xb(i), yb(i), zb(i), betas(i), i=1,Nb)
		close(10)
				
	End if

	! position current point
    x=Model % Nodes % x (nodenumber)
    y=Model % Nodes % y (nodenumber)
    ! Find the height of the current point
    z=Model % Nodes % z (nodenumber)
		
	ind=minloc(sqrt((x-xb)**2+(y-yb)**2),1)
    dist=minval(sqrt((x-xb)**2+(y-yb)**2),1)
    
    beta=betas(ind(1))
    coeff=beta

    Return
End

!------------------------------------------------------------------!
FUNCTION GuessBeta( Model, nodenumber, dumy) RESULT(coeff) !
!------------------------------------------------------------------!
	USE types
	USE DefUtils
    IMPLICIT NONE
	TYPE(Model_t) :: Model
    REAL(kind=dp) :: dumy,coeff
    INTEGER :: nodenumber
    REAL(kind=dp) :: LinearInterp

    REAL(kind=dp),allocatable :: xx(:),yy(:),beta0(:,:)
    REAL(kind=dp) :: x,y,z
    
    INTEGER :: nx,ny
    INTEGER :: i,j
		
    LOGICAL :: FirstTimeBeta=.true.

    SAVE xx,yy,beta0,nx,ny
    SAVE FirstTimeBeta

    if (FirstTimeBeta) then

    	FirstTimeBeta=.False.


        ! open file
        open(10,file='inputs/beta0.xy')
        Read(10,*) nx
        Read(10,*) ny
        ALLOCATE(xx(nx),yy(ny))
        ALLOCATE(beta0(nx,ny))
        Do i=1,nx
        	Do j=1,ny
                read(10,*) xx(i),yy(j),beta0(i,j)
            End Do
		End do
		close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

    coeff = LinearInterp(beta0,xx,yy,nx,ny,x,y)
		
    Return
End

!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!    
