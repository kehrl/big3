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
FUNCTION GuessBeta( Model, nodenumber, dumy) RESULT(coeff) !
!------------------------------------------------------------------!
		USE types
		USE DefUtils
  	IMPLICIT NONE
		TYPE(Model_t) :: Model
  	REAL(kind=dp) :: dumy,coeff
  	INTEGER :: nodenumber
  	REAL(kind=dp) :: LinearInterp

  	REAL(kind=dp),allocatable :: xx(:),yy(:),beta0(:,:)
    REAL(kind=dp) :: x,y,z
    
    INTEGER :: nx,ny
    INTEGER :: i,j
		
    LOGICAL :: FirstTimeBeta=.true.

    SAVE xx,yy,beta0,nx,ny
    SAVE FirstTimeBeta

    if (FirstTimeBeta) then

    	FirstTimeBeta=.False.


        ! open file
        open(10,file='inputs/beta0.xy')
        Read(10,*) nx
        Read(10,*) ny
        ALLOCATE(xx(nx),yy(ny))
        ALLOCATE(beta0(nx,ny))
        Do i=1,nx
        	Do j=1,ny
                read(10,*) xx(i),yy(j),beta0(i,j)
            End Do
		End do
		close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

    coeff = LinearInterp(beta0,xx,yy,nx,ny,x,y)
		
    Return
End

!------------------------------------------------------------------!
FUNCTION Viscosity( Model, nodenumber, dumy) RESULT(eta) !
!------------------------------------------------------------------!
    USE types
		Use DefUtils
    implicit none
		TYPE(Model_t) :: Model
    Real(kind=dp) :: dumy,eta
    INTEGER :: nodenumber

		Real(kind=dp),allocatable :: dem(:,:,:), xx(:), yy(:)
		Real(kind=dp) :: x, y, z, zs , zb, dz
		Real(kind=dp) :: yearinsec, E, alpha
		integer :: nx, ny, nz, k, i, j
		REAL(kind=dp) :: LinearInterp, zsIni, zbIni
		
		TYPE(Variable_t), POINTER :: dSVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)

		
    logical :: Firsttime=.true.

    SAVE dem,xx,yy,nx,ny,nz
    SAVE Firsttime

    if (Firsttime) then

    	Firsttime=.False.

    	! open file
      open(10,file='inputs/flowA.xyz')
      Read(10,*) nx
      Read(10,*) ny
      Read(10,*) nz
      
      allocate(xx(nx), yy(ny))
      allocate(dem(nx, ny, nz))

      do i = 1, nx
      	do j = 1, ny
        	read(10, *) xx(i), yy(j), dem(i, j, :)
        End do
      End do
      close(10)
      
		End if

    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)
    
    dSVariable => VariableGet( Model % Variables, 'dS' )
    IF (ASSOCIATED(dSVariable)) THEN
    	dSPerm    => dSVariable % Perm
    	dSValues  => dSVariable % Values
    ELSE
      CALL FATAL('USF_Init, Viscosity','Could not find variable >dS<')
    END IF
    z = dSValues(dSPerm(nodenumber))

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )		
		
		yearinsec=365.25d0*24*60*60
		
		! Enhanced factor
		E = 3.d0
		
		! Find which vertical layer the current point belongs to
		dz = (zs - zb) / (nz - 1)
		k = int( (z-zb) / dz)+1
    IF (k < 0) THEN
      print *,k,z,zb,zs,dz
    END IF
    
    ! Interpolate the value of the temperature from nearby points in
    ! the layers above and below it
    alpha = (z - (zb + (k - 1) * dz)) / dz
    eta = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) + alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
    
    ! Get the viscosity in the correct units
    eta = ((E*eta*yearinsec)**(-1.0/3.0d0))*1.0e-6
    
    Return
End

!------------------------------------------------------------------!
FUNCTION SurfaceTemperature( Model, nodenumber, dumy) RESULT(Ts) !
!------------------------------------------------------------------!
		USE types
		USE DefUtils
  	IMPLICIT NONE
		TYPE(Model_t) :: Model
  	REAL(kind=dp) :: dumy,Ts
  	INTEGER :: nodenumber
  	REAL(kind=dp) :: LinearInterp

  	REAL(kind=dp),allocatable :: xx(:),yy(:),Tgrid(:,:)
    REAL(kind=dp) :: x,y,z
    
    INTEGER :: nx,ny
    INTEGER :: i,j
		
    LOGICAL :: FirstTimeTs=.true.

    SAVE xx,yy,Tgrid,nx,ny
    SAVE FirstTimeTs

    if (FirstTimeTs) then

    	FirstTimeTs=.False.


        ! open file
        open(10,file='inputs/t2m.xy')
        Read(10,*) nx
        Read(10,*) ny
        ALLOCATE(xx(nx),yy(ny))
        ALLOCATE(Tgrid(nx,ny))
        Do i=1,nx
        	Do j=1,ny
                read(10,*) xx(i),yy(j),Tgrid(i,j)
            End Do
		End do
		close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

    Ts = LinearInterp(Tgrid,xx,yy,nx,ny,x,y)
		
    Return
End

!------------------------------------------------------------------!
FUNCTION Accumulation( Model, nodenumber, dumy) RESULT(a) !
!------------------------------------------------------------------!
		USE types
		USE DefUtils
  	IMPLICIT NONE
		TYPE(Model_t) :: Model
  	REAL(kind=dp) :: dumy,a
  	INTEGER :: nodenumber
  	REAL(kind=dp) :: LinearInterp

  	REAL(kind=dp),allocatable :: xx(:),yy(:),smbgrid(:,:)
    REAL(kind=dp) :: x,y,z
    
    INTEGER :: nx,ny
    INTEGER :: i,j
		
    LOGICAL :: FirstTimea=.true.

    SAVE xx,yy,smbgrid,nx,ny
    SAVE FirstTimea

    if (FirstTimea) then

    	FirstTimea=.False.


        ! open file
        open(10,file='inputs/smb.xy')
        Read(10,*) nx
        Read(10,*) ny
        ALLOCATE(xx(nx),yy(ny))
        ALLOCATE(smbgrid(nx,ny))
        Do i=1,nx
        	Do j=1,ny
                read(10,*) xx(i),yy(j),smbgrid(i,j)
            End Do
		End do
		close(10)
    End if

    ! position current point
    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)

    a = LinearInterp(smbgrid,xx,yy,nx,ny,x,y)
		
    Return
End

!------------------------------------------------------------------!
FUNCTION IceDivideTemperature( Model, nodenumber, dumy) RESULT(T) !
!------------------------------------------------------------------!
    USE types
		Use DefUtils
    implicit none
		TYPE(Model_t) :: Model
    Real(kind=dp) :: dumy,T
    INTEGER :: nodenumber

		Real(kind=dp),allocatable :: dem(:,:,:), xx(:), yy(:)
		Real(kind=dp) :: x, y, z, zs , zb, dz
		Real(kind=dp) :: alpha
		integer :: nx, ny, nz, k, i, j
		REAL(kind=dp) :: LinearInterp, zsIni, zbIni
		
		TYPE(Variable_t), POINTER :: dSVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)

		
    logical :: Firsttime=.true.

    SAVE dem,xx,yy,nx,ny,nz
    SAVE Firsttime

    if (Firsttime) then

    	Firsttime=.False.

    	! open file
      open(10,file='inputs/tsteady.xyz')
      Read(10,*) nx
      Read(10,*) ny
      Read(10,*) nz
      
      allocate(xx(nx), yy(ny))
      allocate(dem(nx, ny, nz))

      do i = 1, nx
      	do j = 1, ny
        	read(10, *) xx(i), yy(j), dem(i, j, :)
        End do
      End do
      close(10)
      
		End if

    x = Model % Nodes % x (nodenumber)
    y = Model % Nodes % y (nodenumber)
    z = Model % Nodes % z (nodenumber)
    
    !dSVariable => VariableGet( Model % Variables, 'dS' )
    !IF (ASSOCIATED(dSVariable)) THEN
    !	dSPerm    => dSVariable % Perm
    !	dSValues  => dSVariable % Values
    !ELSE
    !  CALL FATAL('USF_Init, IceDivideTemperature','Could not find variable >dS<')
    !END IF
    !z = dSValues(dSPerm(nodenumber))

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )		
		
		! Find which vertical layer the current point belongs to
		dz = (zs - zb) / (nz - 1)
		k = int( (z-zb) / dz)+1
    IF (k < 1) THEN
      T=263.0d0
    ELSE
      ! Interpolate the value of the temperature from nearby points in
      ! the layers above and below it
      alpha = (z - (zb + (k - 1) * dz)) / dz
      T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) + alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
    END IF
    
    
        
    Return
End


!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!