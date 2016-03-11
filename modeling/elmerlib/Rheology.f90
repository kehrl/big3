!------------------------------------------------------------------!
FUNCTION Viscosity( Model, nodenumber, dumy) RESULT(eta) !
        !------------------------------------------------------------------!
        USE types
		Use DefUtils
        implicit none
		TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,T,eta,A,yearinsec
        INTEGER :: nodenumber,minind(1)

        Real(kind=dp),allocatable :: inds(:),xx(:),yy(:),zz(:),flowparameters(:)
        Real(kind=dp) :: x,y,z

        integer :: n
        integer :: i
		
        logical :: Firsttime=.true.

        SAVE inds,xx,yy,zz,flowparameters,n
        SAVE Firsttime

        if (Firsttime) then

                Firsttime=.False.

                ! open file
                open(10,file='Inputs/flowparameters.dat')
                Read(10,*) n
                allocate(inds(n),xx(n),yy(n),zz(n),flowparameters(n))
                Do i=1,n
                      read(10,*) inds(i),xx(i),yy(i),zz(i),flowparameters(i)
				End do
				close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)
        z=Model % Nodes % z (nodenumber)
		
		yearinsec=365.25*24*60*60
		minind=minloc((x-xx)**2+(y-yy)**2+(z-zz)**2)
    	A=flowparameters(minind(1))
		!print *,'T A', T,A
		eta=(2.0_dp*A*yearinsec)**(-1.0_dp/3.0_dp)*1.0e-6 
		!print *,'eta is',x,xx(minind(1)),A
        Return
End
