!Compute SIA
!--------------------------------------------------------------------------!

FUNCTION USIA( Model, nodenumber, dumy) RESULT(U) !
        !------------------------------------------------------------------!
        USE types

        implicit none
		TYPE(Model_t) :: Model
		Real(kind=dp) :: dumy
        Real(kind=dp) :: x,y,z
        Real(kind=dp),allocatable :: xx(:),yy(:),surfs(:),beds(:),usurf(:),ubed(:)
        Real(kind=dp) :: U
        INTEGER :: nodenumber
        INTEGER :: nx,i,minind(1)

		Real(kind=dp) :: zb,zs,ub,us
		
		logical :: Firsttime3=.true.
    
        
        SAVE Firsttime3,xx,yy,surfs,beds,usurf,ubed

        if (Firsttime3) then

                Firsttime3=.False.

                ! open file
                open(10,file='Inputs/ubdem.xy')
                Read(10,*) nx
                Allocate(xx(nx),yy(nx),surfs(nx),beds(nx),usurf(nx),ubed(nx))
                Do i=1,nx
                	read(10,*) xx(i),yy(i),surfs(i),beds(i),usurf(i),ubed(i)
				End do
				close(10)
        End if
        
        ! get location
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)
        ! Find the height of the current point
        z=Model % Nodes % z (nodenumber)
        
        ! Get index in file
        minind=minloc((x-xx)**2+(y-yy)**2)
        
        zs=surfs(minind(1))
        zb=beds(minind(1))
        ub=ubed(minind(1))
        us=usurf(minind(1))
        
        U = ub + (1.0_dp - ((zs - z) / (zs - zb))**4) * (us - ub)
        
    
        Return
End

FUNCTION VSIA( Model, nodenumber, dumy) RESULT(V) !
        !------------------------------------------------------------------!
        USE types

        implicit none
		TYPE(Model_t) :: Model
		Real(kind=dp) :: dumy
        Real(kind=dp) :: x,y,z,V
        Real(kind=dp),allocatable :: xx(:),yy(:),surfs(:),beds(:),vbed(:),vsurf(:)
        INTEGER :: nodenumber
        INTEGER :: nx,i,minind(1)

        Real(kind=dp) :: zs,zb,vb,vs

		logical :: Firsttime3=.true.
         
        
        SAVE Firsttime3,xx,yy,surfs,beds,vbed,vsurf

        if (Firsttime3) then

                Firsttime3=.False.

                ! open file
                open(10,file='Inputs/vbdem.xy')
                Read(10,*) nx
                Allocate(xx(nx),yy(nx),surfs(nx),beds(nx),vsurf(nx),vbed(nx))
                Do i=1,nx
                	read(10,*) xx(i),yy(i),surfs(i),beds(i),vsurf(i),vbed(i)
				End do
				close(10)
        End if
        
        ! get location
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)
        ! Find the height of the current point
        z=Model % Nodes % z (nodenumber)
        
        ! Get index in file
        minind=minloc((x-xx)**2+(y-yy)**2)
        
        zs=surfs(minind(1))
        zb=beds(minind(1))
        vs=vsurf(minind(1))
        vb=vbed(minind(1))
        
        V = vb + (1.0_dp - ((zs -z) / (zs - zb))**4) * (vs - vb)
        
        Return
End
    