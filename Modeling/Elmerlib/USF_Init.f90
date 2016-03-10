!------------------------------------------------------------------!
FUNCTION UIni( Model, nodenumber, dumy) RESULT(U) !
!------------------------------------------------------------------!
	USE types

  implicit none
	TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,U
  INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp
  Real(kind=dp) :: MaxU

  Integer :: nx,ny
  Integer :: i,j

  Character(len=16) :: glacier
  Character(len=MAX_NAME_LEN) :: filin

  Logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny,MaxU
  SAVE Firsttime

  if (Firsttime) then
		MaxU = 0.d0

    Firsttime=.False.

    ! open file
    Open(10,file='inputs/udem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		Close(10)
  End if

  ! position current point
  x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  U=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION VIni( Model, nodenumber, dumy) RESULT(U) !
!------------------------------------------------------------------!
	USE types

  Implicit none
  TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,U
  INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp
  Real(kind=dp) :: MaxV

  Integer :: nx,ny
  Integer :: i,j

  Character(len=16) :: glacier
  Character(len=MAX_NAME_LEN) :: filin

  Logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny,MaxV
  SAVE Firsttime

  if (Firsttime) then
  	MaxV = 0.d0

    Firsttime=.False.

    ! open file
    Open(10,file='inputs/vdem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		close(10)
  End if

  ! position current point
  x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  U=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION zsIni( Model, nodenumber, dumy) RESULT(zs) !
!------------------------------------------------------------------!
	USE types

  Implicit none
	TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,zs
  INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp

  integer :: nx,ny
  integer :: i,j

  Character(len=16) :: glacier
  Character(len=MAX_NAME_LEN) :: filin='inputs/zsdem.xy'

  Logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny
  SAVE Firsttime

  if (Firsttime) then
		Firsttime=.False.

    ! open file
    open(10,file='inputs/zsdem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		close(10)
  End if

  ! position current point
  x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  zs=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION zbIni( Model, nodenumber, dumy) RESULT(zb) !
!------------------------------------------------------------------!
	USE types
	implicit none
	TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,zb
	INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp

  integer :: nx,ny
  integer :: i,j

	character(len=16) :: glacier
  character(len=MAX_NAME_LEN) :: filin='inputs/zbdem.xy'

  logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny
  SAVE Firsttime

  if (Firsttime) then
		Firsttime=.False.

    ! open file
    open(10,file='inputs/zbdem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		close(10)
  End if

  ! position current point
	x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  zb=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION UbIni( Model, nodenumber, dumy) RESULT(ub) !
!------------------------------------------------------------------!
	USE types
	implicit none
	TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,ub
	INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp

  integer :: nx,ny
  integer :: i,j

	character(len=16) :: glacier
  character(len=MAX_NAME_LEN) :: filin='inputs/ubdem.xy'

  logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny
  SAVE Firsttime

  if (Firsttime) then
		Firsttime=.False.

    ! open file
    open(10,file='inputs/ubdem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		close(10)
  End if

  ! position current point
	x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  ub=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION VbIni( Model, nodenumber, dumy) RESULT(vb) !
!------------------------------------------------------------------!
	USE types
	implicit none
	TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,vb
	INTEGER :: nodenumber

  Real(kind=dp),allocatable :: dem(:,:),xx(:),yy(:)
  Real(kind=dp) :: x,y
  Real(kind=dp) :: LinearInterp

  integer :: nx,ny
  integer :: i,j

	character(len=16) :: glacier
  character(len=MAX_NAME_LEN) :: filin='inputs/vbdem.xy'

  logical :: Firsttime=.true.

  SAVE dem,xx,yy,nx,ny
  SAVE Firsttime

  if (Firsttime) then
		Firsttime=.False.

    ! open file
    open(10,file='inputs/vbdem.xy')
    Read(10,*) nx
    Read(10,*) ny
    Allocate(xx(nx),yy(ny))
    Allocate(dem(nx,ny))
    Do i=1,nx
    	Do j=1,ny
      	Read(10,*) xx(i),yy(j),dem(i,j)
      End Do
		End do
		close(10)
  End if

  ! position current point
	x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)

  vb=LinearInterp(dem,xx,yy,nx,ny,x,y)

  Return
End

!------------------------------------------------------------------!
FUNCTION UWa( Model, nodenumber, dumy) RESULT(U)                   !
!------------------------------------------------------------------!
	USE types

  implicit none
  TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,U
  INTEGER :: nodenumber

  Real(kind=dp) :: x, y, z, zs, zb, us, ub
  Real(kind=dp), external :: UIni, UbIni, zsIni, zbIni

  ! position current point
  x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)
  z=Model % Nodes % z (nodenumber)

  zs = zsIni( Model, nodenumber, dumy )
  zb = zbIni( Model, nodenumber, dumy )

  us = UIni( Model, nodenumber, dumy )
  ub = UbIni( Model, nodenumber, dumy )

  U = ub + (1.0_dp - ((zs - z) / (zs - zb))**4) * (us - ub)

	Return 
End


!------------------------------------------------------------------!
FUNCTION VWa( Model, nodenumber, dumy) RESULT(U)                   !
!------------------------------------------------------------------!
	USE types

  implicit none
  TYPE(Model_t) :: Model
  Real(kind=dp) :: dumy,U
  INTEGER :: nodenumber

  Real(kind=dp) :: x, y, z, zs, zb, vs, vb
  Real(kind=dp), external :: VIni, VbIni, zsIni, zbIni

  ! position current point
  x=Model % Nodes % x (nodenumber)
  y=Model % Nodes % y (nodenumber)
  z=Model % Nodes % z (nodenumber)

  zs = zsIni( Model, nodenumber, dumy )
  zb = zbIni( Model, nodenumber, dumy )

  vs = VIni( Model, nodenumber, dumy )
  vb = VbIni( Model, nodenumber, dumy )

  U = vb + (1.0_dp - ((zs - z) / (zs - zb))**4) * (vs - vb)

  Return 
End



!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!