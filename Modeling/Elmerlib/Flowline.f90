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
   	
   	logical :: found
   	logical :: Firsttime1=.true.
   	logical :: Firsttime2=.true.
        
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION ShapeFactor (Model, nodenumber) RESULT(output)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber,ind, Rw, i
   	REAL(KIND=dp) :: x, dist, mindist, output, ratio
   	REAL(KIND=dp), ALLOCATABLE :: xw(:), yw(:)
   	logical :: found
   
    logical :: Firsttime3=.true.
        
   	SAVE Firsttime3,xw,yw,Rw
   

   	if (Firsttime3) then
   		! Read file
   		Firsttime3=.False.
   		
   		OPEN(10,file="Inputs/shapefactor.dat")
   		Read(10,*) Rw
   		ALLOCATE(xw(Rw), yw(Rw))
   		READ(10,*)(xw(i), yw(i), i=1,Rw)
   		CLOSE(10)
   	End if

   	x=Model % Nodes % x (nodenumber)   
   
   	found = .false.
   	mindist=dabs(2*(xw(2)-xw(3)))

   	do 20, i=1,Rw
      	dist=dabs(x-xw(i))
      	if (dist<=mindist .and. xw(i)<=x) then
        	mindist=dist
        	ind=i
        	found = .true.
      	endif 
   	20 enddo

   	if (.not.found) then
      	print *, 'Could not find a suitable shapefactor to interpolate ',x
   	else
      	ratio=(x-xw(ind))/(xw(ind+1)-xw(ind))
      	output=yw(ind)+ratio*(yw(ind+1)-yw(ind))
      	!
   	endif
    !print *,'',x,output
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
   	logical :: found
   
    logical :: Firsttime4=.true.
        
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
    Real(kind=dp) :: x,y,z
    Integer :: nodenumber,minind,RA,i
		
    Logical :: Firsttime5=.true.
    logical :: found

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
    x=Model % Nodes % x (nodenumber)
    y=Model % Nodes % y (nodenumber)

   	found = .false.		
	
	mindist=1000
   	do 20, i=1,RA
      	dist=sqrt((x-xA(i))**2+(y-yA(i))**2)
      	if (dist<=mindist) then
        	mindist=dist
        	minind=i
        	found = .true.
      	endif 
   	20 enddo
	
	yearinsec=365.25*24*60*60
    Aval=A(minind)
	eta=(2.0*Aval*yearinsec)**(-1.0_dp/3.0_dp)*1.0e-6
	
    Return
End

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Slip_Linear (Model, nodenumber) RESULT(coefficient1)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber,ind, R6, i
   	REAL(KIND=dp) :: x, dist, mindist, ratio, coefficient1
   	REAL(KIND=dp), ALLOCATABLE :: xb(:), beta(:), yb(:)
   	logical :: found

    Logical :: Firsttime6=.true.

    SAVE xb,yb,beta
    SAVE Firsttime6,R6

    if (Firsttime6) then

    	Firsttime6=.False.

   		! Read file
   		OPEN(10,file="Inputs/beta_linear.xy")
   		Read(10,*) R6
   		ALLOCATE(xb(R6), yb(R6), beta(R6))
   		READ(10,*)(xb(i), yb(i),beta(i), i=1,R6)
   		CLOSE(10)
   		
    End if
   
	x=Model % Nodes % x (nodenumber)
   
   	found = .false.
   	mindist=dabs(2*(xb(1)-xb(2)))
   	do 20, i=1,R6
      	dist=dabs(x-xb(i))
      	if (dist<=mindist .and. xb(i)<=x) then
        	mindist=dist
        	ind=i
        	found = .true.
      	endif 
   	20 enddo

   	if (.not.found) then
      	print *, 'Could not find a suitable beta to interpolate ',x
   	else
      	ratio=(x-xb(ind))/(xb(ind+1)-xb(ind))
      	coefficient1=beta(ind)+ratio*(beta(ind+1)-beta(ind))
   	endif
   
   	Return
End


FUNCTION Slip_Weertman (Model, nodenumber) RESULT(coefficient1)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber,ind, R6, i
   	REAL(KIND=dp) :: x, dist, mindist, ratio, coefficient1
   	REAL(KIND=dp), ALLOCATABLE :: xb(:), beta(:), yb(:)
   	logical :: found

    Logical :: Firsttime6=.true.

    SAVE xb,yb,beta
    SAVE Firsttime6,R6

    if (Firsttime6) then

    	Firsttime6=.False.

   		! Read file
   		OPEN(10,file="Inputs/beta_weertman.xy")
   		Read(10,*) R6
   		ALLOCATE(xb(R6), yb(R6), beta(R6))
   		READ(10,*)(xb(i), yb(i),beta(i), i=1,R6)
   		CLOSE(10)
   		
    End if
   
	x=Model % Nodes % x (nodenumber)
   
   	found = .false.
   	mindist=dabs(2*(xb(1)-xb(2)))
   	do 20, i=1,R6
      	dist=dabs(x-xb(i))
      	if (dist<=mindist) then
        	mindist=dist
        	ind=i
        	found = .true.
      	endif 
   	20 enddo

   	if (.not.found) then
      	print *, 'Could not find a suitable beta to interpolate ',x
   	else
   		if (xb(ind) >= x) then 
      		ratio=(x-xb(ind))/(xb(ind+1)-xb(ind))
      		coefficient1=beta(ind)+ratio*(beta(ind+1)-beta(ind))
      	else 
      		ratio=(x-xb(ind))/(xb(ind)-(xb(ind-1)))
      		coefficient1=beta(ind)+ratio*(beta(ind)-beta(ind-1))
      	endif
   	endif
   
   	Return
End

Function Bedrock (Model, nodenumber, znode) RESULT(BED)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	INTEGER :: nodenumber, ind
   	REAL(KIND=dp) :: Mask
   	INTEGER :: NMAX, i, R7, DIM
   	REAL(KIND=dp) :: x, y, z, znode
   	REAL(KIND=dp) :: BED,ratio,mindist,dist
   	REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), zbed(:)
   	Logical :: Firsttime7=.true.
	Logical :: found

    SAVE xbed,ybed
    SAVE Firsttime7,R7

    if (Firsttime7) then

    	Firsttime7=.False.

   		! Load bedrock
   		OPEN(10,file="Inputs/roughbed.dat")
   		Read(10,*) R7
   		ALLOCATE(xbed(R7), ybed(R7))
   		READ(10,*)(xbed(i), ybed(i), i=1,R7)
   		CLOSE(10)
   		
    End if

	x=Model % Nodes % x (nodenumber)

  	found = .false.
   	mindist=dabs(2*(xbed(1)-xbed(2)))
   	do 20, i=1,R7
      	dist=dabs(x-xbed(i))
      	if (dist<=mindist .and. xbed(i)<=x) then
        	mindist=dist
        	ind=i
        	found = .true.
      	endif 
   	20 enddo

   	if (.not.found) then
      	print *, 'Could not find a suitable bed to interpolate ',x
   	else
      	ratio=(x-xbed(ind))/(xbed(ind+1)-xbed(ind))
      	BED=ybed(ind)+ratio*(ybed(ind+1)-ybed(ind))
      	print *,'bed is',x,BED
   	endif

   Return
End