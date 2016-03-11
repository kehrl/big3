!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION Viscosity( Model, nodenumber) RESULT(eta) !
    USE types
	Use DefUtils
    implicit none
	TYPE(Model_t) :: Model
    Real(kind=dp) :: T,eta,Aval,yearinsec,mindist,dist
    Real(kind=dp),allocatable :: xA(:),yA(:),A(:),TA(:)
    Real(kind=dp) :: x,y,z,xnew
    INTEGER :: nodenumber,minind,RA,i
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName = 'Viscosity'
		
    LOGICAL :: Firsttime5=.true.
    LOGICAL :: found
     
    ! Save values for future use
    SAVE xA,yA,A,RA,TA
    SAVE Firsttime5

    if (Firsttime5) then

    	Firsttime5=.False.

        ! open file
        Open(10,file='Inputs/flowparameters.dat')
        Read(10,*) RA
        Allocate(xA(RA),yA(RA),TA(RA),A(RA))
        Do i=1,RA
        	Read(10,*) xA(i),yA(i),TA(i),A(i)
		End do
		Close(10)
    End if

    ! Get position
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

	! Now find data closest to current position
   	found = .FALSE.		
	
	IF (x > xA(RA)) THEN
		xnew = XA(RA)
		mindist=10000
   	ELSE
   		xnew = x
   		mindist=5000
   	END IF 
   	
   	DO i=1,RA
      	dist=sqrt((xnew-xA(i))**2+(y-yA(i))**2)
      	IF (dist<=mindist) THEN
        mindist=dist
        minind=i
        found = .true.
      	END IF
   	END DO
	
	! Linearly interpolate to the current position
	IF (found) THEN
		yearinsec=365.25*24*60*60
    	Aval=A(minind)
		eta=(2.0*Aval*yearinsec)**(-1.0_dp/3.0_dp)*1.0e-6
	else
		print *,'No viscosity at',x,y
		CALL FATAL(SolverName,'No viscosity found for above coordinates')
	end if	

    Return
End
