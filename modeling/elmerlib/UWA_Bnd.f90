!------------------------------------------------------------------!
        FUNCTION UBnd( Model, nodenumber, dumy) RESULT(U) !
        !------------------------------------------------------------------!
        USE types

        implicit none
TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y,z
        Real(kind=dp) :: LinearInterp
        Real(kind=dp) :: MaxU

        integer :: nx,ny
        integer :: i,j

        character(len=16) :: glacier
        character(len=MAX_NAME_LEN) :: filin='Inputs/UDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny, MaxU
        SAVE Firsttime

        if (Firsttime) then
MaxU = 0.d0

                Firsttime=.False.

                ! open file
                open(10,file='Inputs/UDEM.xy')
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
End do
close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)
        ! Get the height of the current point too
        z=Model % Nodes % z (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return
End


        !------------------------------------------------------------------!
        FUNCTION VBnd( Model, nodenumber, dumy) RESULT(U) !
        !------------------------------------------------------------------!
        USE types

        implicit none
TYPE(Model_t) :: Model
        Real(kind=dp) :: dumy,U
        INTEGER :: nodenumber

        Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
        Real(kind=dp) :: x,y,z
        Real(kind=dp) :: LinearInterp
        Real(kind=dp) :: MaxV

        integer :: nx,ny
        integer :: i,j

        character(len=16) :: glacier
        character(len=MAX_NAME_LEN) :: filin='Inputs/VDEM.xy'

        logical :: Firsttime=.true.

        SAVE dem,xx,yy,nx,ny, MaxV
        SAVE Firsttime

        if (Firsttime) then
MaxV = 0.d0

                Firsttime=.False.

                ! open file
                open(10,file='Inputs/VDEM.xy')
                Read(10,*) nx
                Read(10,*) ny
                allocate(xx(nx),yy(ny))
                Allocate(dem(nx,ny))
                Do i=1,nx
                   Do j=1,ny
                      read(10,*) xx(i),yy(j),dem(i,j)
                   End Do
End do
close(10)
        End if

        ! position current point
        x=Model % Nodes % x (nodenumber)
        y=Model % Nodes % y (nodenumber)
        ! Find the height of the current point
        z=Model % Nodes % z (nodenumber)

        U=LinearInterp(dem,xx,yy,nx,ny,x,y)

        Return
End


        !------------------------------------------------------------------!
        include 'Interp.f90' !
        !------------------------------------------------------------------!