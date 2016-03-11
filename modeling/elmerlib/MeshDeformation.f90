FUNCTION DisplacYX (Model, nodenumber, BasalVelocity) RESULT(test)
	USE DefUtils
	USE types
	IMPLICIT NONE
	TYPE(Model_t) :: Model
  	TYPE(Variable_t), POINTER :: StressSol, MeshSol
  	REAL(KIND=dp), POINTER :: MeshUpdate(:), MeshVelocity(:)
  	INTEGER, POINTER :: TPerm(:), MeshPerm(:), StressPerm(:)
    INTEGER :: STDOFs, nodenumber
    REAL(KIND=dp) :: test, BasalVelocity, dt
    LOGICAL :: Found

  	StressSol => VariableGet(CurrentModel % Solver % Mesh % Variables, 'Mesh Update' )
  	If( ListGetLogical(CurrentModel % Solver % Values,'Ignore Displacement',Found) ) then
    	print *,' NO DISPLACEMENT---------------'
 	End If


  	If ( ASSOCIATED( StressSol ) )  then
    	StressPerm   => StressSol % Perm
    	STDOFs       =  StressSol % DOFs
    	MeshUpdate => StressSol % Values	
	End If
	
	test=MeshUpdate(nodenumber) + BasalVelocity*(1/365.25)
	print *,'displacement', test
	RETURN
END


Function DisplaceX (Model, nodenumber, BasalVelocity) RESULT(Xdel)
   	USE types
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(Solver_t), TARGET :: Solver
   	TYPE(Nodes_t) :: ElementNodes
   	TYPE(Variable_t), POINTER :: TimeVar, TimeStepSizeVar

   	INTEGER :: nodenumber, i, NMAX
   	
   	REAL(KIND=dp) :: x, y, Xdel, dt 
	REAL(KIND=dp), ALLOCATABLE :: XIni(:),Xdeltot(:), t0(:), t(:)
   	REAL(KIND=dp) :: BasalVelocity, time
   	
   	LOGICAL :: FirstTimeX = .TRUE.

	SAVE XIni, FirstTimeX, Xdeltot, t0, t
	
	TimeVar => VariableGet( Model % Variables,'Time')
	
	! Allocate and set initial values
	if (FirstTimeX) then
		FirstTimeX = .FALSE.
		
		NMAX = Model % NumberOfNodes 
	  	ALLOCATE(XIni(NMAX),Xdeltot(NMAX),t0(NMAX),t(NMAX))
        DO i = 1, NMAX
        	XIni(i) = Model % Nodes % x (i)
        	Xdeltot(i) = 0
        	t0(i)=0
        	t(i)=TimeVar % Values(1)
        END DO
	end if
	
	! Find current model parameters
	x = Model % Nodes % x (nodenumber)
	t(nodenumber)=TimeVar % Values(1) 

	TimeStepSizeVar => VariableGet( Model % Variables,'Timestep Size')
	dt=TimeStepSizeVar % Values(1) 
	
	! Update displacement if we are on a new coordinate or new timestep
	if (t(nodenumber) > t0(nodenumber)) then
		Xdel = Xdeltot(nodenumber)+dt*BasalVelocity
		Xdeltot(nodenumber)=Xdel
		t0(nodenumber)=t(nodenumber)
	else 
		Xdel = Xdeltot(nodenumber)
	end if
	
   Return
End

Function DisplaceY (Model, nodenumber, BasalVelocity) RESULT(dely)
   	USE types
   	USE CoordinateSystems
   	USE SolverUtils
   	USE ElementDescription
   	USE DefUtils
   	IMPLICIT NONE
   	TYPE(Model_t) :: Model
   	TYPE(variable_t), POINTER :: TimeVar, TimeStepSizeVar
   	TYPE(Solver_t), TARGET :: Solver
   	REAL(KIND=dp) :: Mask
   	INTEGER :: nodenumber, i, j, DIM, NMAX, R, ind
   	REAL(KIND=dp) :: x, y, Xdel, dt, dely 
	REAL(KIND=dp), ALLOCATABLE :: XIni(:),y0(:),Xdeltot(:), totaldely(:), t0(:), t(:)
	REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), zbed(:)
   	REAL(KIND=dp) :: BasalVelocity, ratio, dist, mindist,test
   	
   	LOGICAL :: FirstTimeA2 = .TRUE.

	SAVE y0, FirstTimeA2, Xdeltot, totaldely, t0, t
	SAVE xbed, ybed, R
	
	TimeVar => VariableGet( Model % Variables,'Time')
 

    if (FirsttimeA2) then

    	FirsttimeA2=.False.

  		! Load bedrock
   		OPEN(10,file="Inputs/roughbed.dat")
   			Read(10,*) R
   			ALLOCATE(xbed(R), ybed(R), zbed(R))
   			READ(10,*)(xbed(i), ybed(i), zbed(i), i=1,R)
   		CLOSE(10)
		
		NMAX = Model % NumberOfNodes 
	  	ALLOCATE(XIni(NMAX),y0(NMAX),Xdeltot(NMAX),totaldely(NMAX),t0(NMAX),t(NMAX))
        DO j = 1, NMAX
        	XIni(j) = Model % Nodes % x (j)
        	y0(j) = Model % Nodes % y (j)
        	Xdeltot(j) = 0
        	totaldely(j) = 0
        	t0(j)=0
        	t(j)=TimeVar % Values(1)
        END DO
	end if
	
	t(nodenumber)=TimeVar % Values(1) 
	TimeStepSizeVar => VariableGet( Model % Variables,'Timestep Size')
	dt=TimeStepSizeVar % Values(1) 
	
	if (t(nodenumber) > t0(nodenumber)) then
		Xdel = Xdeltot(nodenumber)+dt*BasalVelocity
		Xdeltot(nodenumber)=Xdel
		x = Model % Nodes % x (nodenumber) + Xdel

   		mindist=dabs(2*(xbed(1)-xbed(2)))
   		do 20, i=1,R
      		dist=dabs(x-xbed(i))
      		if (dist<=mindist .and. xbed(i)<=x) then
        		mindist=dist
        		ind=i
      		endif 
   		20 enddo
	
      	ratio=(x-xbed(ind))/(xbed(ind+1)-xbed(ind))
      	y=ybed(ind)+ratio*(ybed(ind+1)-ybed(ind))
		dely = y - y0(nodenumber)
		totaldely(nodenumber)=dely
		t0(nodenumber)=t(nodenumber)
	else 
		Xdel = Xdeltot(nodenumber)
		dely = totaldely(nodenumber)
	end if
	
   Return
End