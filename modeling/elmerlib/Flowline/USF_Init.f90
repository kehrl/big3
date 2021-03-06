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
   	LOGICAL :: found
   
    LOGICAL :: FirsttimeInit=.True.
    
    ! Save variables for future use    
   	SAVE FirsttimeInit,xs,vs,Rs
   
    ! Load velocity data
   	if (FirsttimeInit) then
   		! Read file
   		FirsttimeInit=.False.
   		
   		OPEN(10,file="Inputs/velocity.dat")
   		Read(10,*) Rs
   		ALLOCATE(xs(Rs), vs(Rs))
   		READ(10,*)(xs(i), vs(i), i=1,Rs)
   		CLOSE(10)
   	End if
   
    ! Find current position
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

    ! Linearly interpolate to current position
   	if (.not.found) then
      	print *, 'Could not find a suitable velocity to interpolate ',x
   	else
      	ratio=(x-xs(ind))/(xs(ind+1)-xs(ind))
      	vel=vs(ind)+ratio*(vs(ind+1)-vs(ind))
   	endif
   
   	Return
End