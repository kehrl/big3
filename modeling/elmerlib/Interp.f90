 !------------------------------------------------------------------!
        Function LinearInterp(dem,xx,yy,nx,ny,x,y) Result(InterP1) !
        !------------------------------------------------------------------!
        USE TYPES
        implicit none
        REAL(KIND=dp) :: dem(nx,ny),xx(nx),yy(ny)
        REAL(KIND=dp) :: Dx,Dy,DxDy
        Real(kind=dp) :: x,y,x_1,y_1,dist,B(4)
        Real(kind=dp) :: InterP1
        integer :: nx,ny,i,j
        integer :: nx_1,ny_1

        logical :: found

        Dx=(xx(nx)-xx(1))/(nx-1)
        Dy=(yy(ny)-yy(1))/(ny-1)
        DxDy=Dx*Dy

        ! lower left point in DEM
        nx_1=floor((x-xx(1))/Dx) + 1
        ny_1=floor((y-yy(1))/Dy) + 1
        nx_1=min(nx_1,nx-1)
        ny_1=min(ny_1,ny-1)

        x_1=xx(nx_1)
        y_1=yy(ny_1)


        ! DEM Value in surroundings points
        ! 4 ----- 3
        ! | |
        ! 1 ----- 2
        B(1)=dem(nx_1,ny_1)
        B(2)=dem(nx_1+1,ny_1)
        B(3)=dem(nx_1+1,ny_1+1)
        B(4)=dem(nx_1,ny_1+1)

        found = .false.

        if (minval(B)/=-2e+9) then
            ! Linear Interpolation at Point x,y
            InterP1=(x-x_1)*(y-y_1)*(B(3)+B(1)-B(2)-B(4))/DxDy
            InterP1=InterP1+(x-x_1)*(B(2)-B(1))/Dx+(y-y_1)*(B(4)-B(1))/Dy+B(1)
            found = .true.
        else
        do i=0,1
                do j=0,1
                    dist = max( dabs(x-xx(nx_1+i)),dabs(y-yy(ny_1+j)) )
                    if (dist<=0.5*dx .and. dem(nx_1+i,ny_1+j)/=-2e+9) then
                        InterP1 = dem(nx_1+i,ny_1+j)
                        found = .true.
                    endif
enddo
            enddo
        endif

if (.not.found) then
print *, 'Could not find a suitable point to interpolate ',x,y
        endif

Return
End