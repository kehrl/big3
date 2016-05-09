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