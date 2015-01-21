FUNCTION LateralFrictionCoefficient( Model, nodenumber, dumy) RESULT(kcoeff) !
    USE types
	Use DefUtils
    implicit none
	TYPE(Model_t) :: Model
    Real(kind=dp) :: kcoeff, halfwidth, eta, dumy
    Real(kind=dp),allocatable :: xwidths(:),widths(:),junk(:)
    Real(kind=dp) :: x,y,mindist,dist,ratio,n,yearinsec,rhoi
    INTEGER :: nodenumber,minind,R,i
    TYPE(Variable_t), POINTER :: ViscosityVariable
    INTEGER, POINTER :: ViscosityPerm(:)
  	REAL(KIND=dp), POINTER :: ViscosityValues(:)
		
    LOGICAL :: FirsttimeLFC=.true.
    LOGICAL :: found

    SAVE xwidths, widths
    SAVE FirsttimeLFC
	
	if (FirsttimeLFC) then

    	FirsttimeLFC=.False.

        ! open file
        Open(10,file='Inputs/width.dat')
        Read(10,*) R
        Allocate(xwidths(R),widths(R),junk(R))
        Do i=1,R
        	Read(10,*) xwidths(i),widths(i),junk(i)
		End do
		Close(10)
    End if
    
    
    
    x = Model % Nodes % x (nodenumber)  
    
    ! Now find data closest to current position
   	found = .FALSE.		
	
   	mindist=100000
   	DO i=1,SIZE(xwidths)
      	dist=dabs(xwidths(i)-x)
      	IF (dist<=mindist) THEN
        mindist=dist
        minind=i
        found = .true.
      	END IF
   	END DO

	! Interpolate width to current position, and then find halfwidth for friction coefficient calculation
   	if (.not.found) then
      	print *, 'LateralFrictionCoefficient: Could not find a suitable width to interpolate ',x
   	else
   		if (xwidths(minind) >= x) then 
      		ratio=(x-xwidths(minind))/(xwidths(minind+1)-xwidths(minind))
      		halfwidth=(widths(minind)+ratio*(widths(minind+1)-widths(minind)))/2
      	else 
      		ratio=(x-xwidths(minind))/(xwidths(minind)-(xwidths(minind-1)))
      		halfwidth=(widths(minind)+ratio*(widths(minind)-widths(minind-1)))/2
      	endif
   	endif
    
	
	! Get viscosity at location
	ViscosityVariable => VariableGet( Model % Variables, 'Viscosity' )
   	IF (ASSOCIATED(ViscosityVariable)) THEN
    	ViscosityPerm    => ViscosityVariable % Perm
    	ViscosityValues  => ViscosityVariable % Values
    ELSE
        CALL FATAL('LateralFrictionCoefficient','Could not find variable Viscosity')
  	END IF
  	
  	eta = ViscosityValues(ViscosityPerm(nodenumber))
	n = 3
	yearinsec = 365.25*24*60*60
	rhoi = 917.0/(1.0e6*yearinsec**2)
	
	kcoeff = eta * (n+1)**(1/n) / (halfwidth**(1+(1/n)))
	
    Return
End
