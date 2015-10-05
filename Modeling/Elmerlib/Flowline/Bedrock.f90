!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION BedFromFile( Model, nodenumber) RESULT(zb) !
    USE types
	Use DefUtils
    implicit none
	TYPE(Model_t) :: Model
    Real(kind=dp) :: mindist,dist,zb,yearinsec,ratio
    Real(kind=dp),allocatable :: xbed(:),zbed(:)
    Real(kind=dp) :: x,y,z,xnew
    INTEGER :: nodenumber,minind,RB,i
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'BedFromFile'
		
    LOGICAL :: FirstTimeBed=.true.
    LOGICAL :: found

    SAVE xbed,zbed,RB
    SAVE FirstTimeBed

    if (FirstTimeBed) then

    	FirstTimeBed=.False.

        ! open file
        Open(10,file='Inputs/roughbed.dat')
        Read(10,*) RB
        Allocate(xbed(RB),zbed(RB))
        Do i=1,RB
        	Read(10,*) xbed(i),zbed(i)
		End do
		Close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

	! Now find data closest to current position
   	found = .FALSE.		
	
	IF (x > xbed(RB)) THEN
		xnew = xbed(RB)
		mindist=10000
   	ELSE
   		xnew = x
   		mindist=10000
   	END IF 
   	
   	DO i=1,RB
      	dist=abs(xnew-xbed(i))
      	IF (dist<=mindist) THEN
            mindist=dist
            minind=i
            found = .true.
      	END IF
   	END DO
	
	IF (found) THEN
   		if (xbed(minind) >= xnew) then 
      		ratio=(xnew-xbed(minind))/(xbed(minind+1)-xbed(minind))
      		zb=(zbed(minind)+ratio*(zbed(minind+1)-zbed(minind)))
      	else 
      		ratio=(xnew-xbed(minind))/(xbed(minind)-(xbed(minind-1)))
      		zb=(zbed(minind)+ratio*(zbed(minind)-zbed(minind-1)))
      	endif
      	!print *,'Zb at',xnew,zb
	else
	    print *,'No bed at',xnew
		CALL FATAL(SolverName,'No bed found for above coordinates')
	end if	
    Return
End
