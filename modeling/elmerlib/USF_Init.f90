!------------------------------------------------------------------!
FUNCTION UIni( Model, nodenumber, dumy) RESULT(U) !
!------------------------------------------------------------------!
  USE types

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, U
  INTEGER :: nodenumber

  REAL(kind=dp), ALLOCATABLE :: dem(:,:),xx(:),yy(:)
  REAL(kind=dp) :: x,y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN

    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/udem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx), yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
   CLOSE(10)
  
  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  U = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION VIni( Model, nodenumber, dumy) RESULT(V) !
!------------------------------------------------------------------!
  USE types

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, V
  INTEGER :: nodenumber

  REAL(kind=dp), ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x, y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN

    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/vdem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx), yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  V = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION zsIni( Model, nodenumber, dumy) RESULT(zs) !
!------------------------------------------------------------------!
  USE types

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy,zs
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x, y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN
    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/zsdem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx),yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Nodes % x (nodenumber)
  y = Model % Nodes % y (nodenumber)

  zs = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION zbIni( Model, nodenumber, dumy) RESULT(zb) !
!------------------------------------------------------------------!
  USE types
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy,zb
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x,y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx,ny
  INTEGER :: i,j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN
    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/zbdem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx), yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i),yy(j),dem(i,j)
      END DO
    END DO
    CLOSE(10)
  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  zb = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION Bedrock( Model, nodenumber, dumy) RESULT(zb) !
!------------------------------------------------------------------!
  USE types
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, zb
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x, y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN
    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/bedrock.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx),yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  zb = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END


!------------------------------------------------------------------!
FUNCTION UbIni( Model, nodenumber, dumy) RESULT(ub) !
!------------------------------------------------------------------!
  USE types
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, ub
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x,y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j
  
  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN
    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/ubdem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE( xx(nx), yy(ny))
    ALLOCATE( dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  ub = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION VbIni( Model, nodenumber, dumy) RESULT(vb) !
!------------------------------------------------------------------!
  USE types
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, vb
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:), xx(:), yy(:)
  REAL(kind=dp) :: x, y
  REAL(kind=dp) :: LinearInterp

  INTEGER :: nx, ny, i, j

  LOGICAL :: Firsttime=.true.

  SAVE dem, xx, yy, nx, ny
  SAVE Firsttime

  IF (Firsttime) THEN
    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/vbdem.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(dem(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), dem(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  vb = LinearInterp(dem, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION UWa( Model, nodenumber, dumy) RESULT(U)                   !
!------------------------------------------------------------------!
  USE types

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, U
  INTEGER :: nodenumber

  REAL(kind=dp) :: x, y, z, zs, zb, us, ub
  REAL(kind=dp), EXTERNAL :: UIni, UbIni, zsIni, zbIni

  ! position current point
  x = Model % Nodes % x (nodenumber)
  y = Model % Nodes % y (nodenumber)
  z = Model % Nodes % z (nodenumber)

  zs = zsIni( Model, nodenumber, dumy )
  zb = zbIni( Model, nodenumber, dumy )

  us = UIni( Model, nodenumber, dumy )
  ub = UbIni( Model, nodenumber, dumy )

  IF (z > zs) THEN
    U = us
  ELSE
    U = ub + (1.0_dp - ((zs - z) / (zs - zb))**4) * (us - ub)
  END IF
  
  RETURN 
END


!------------------------------------------------------------------!
FUNCTION VWa( Model, nodenumber, dumy) RESULT(V)                   !
!------------------------------------------------------------------!
  USE types

  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, V
  INTEGER :: nodenumber

  REAL(kind=dp) :: x, y, z, zs, zb, vs, vb
  REAL(kind=dp), EXTERNAL :: VIni, VbIni, zsIni, zbIni

  ! position current point
  x = Model % Nodes % x (nodenumber)
  y = Model % Nodes % y (nodenumber)
  z = Model % Nodes % z (nodenumber)

  zs = zsIni( Model, nodenumber, dumy )
  zb = zbIni( Model, nodenumber, dumy )

  vs = VIni( Model, nodenumber, dumy )
  vb = VbIni( Model, nodenumber, dumy )

  IF (z > zs) THEN
    V = vs
  ELSE
    V = vb + (1.0_dp - ((zs - z) / (zs - zb))**4) * (vs - vb)
  END IF
  
  RETURN 
END

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

  REAL(kind=dp),ALLOCATABLE :: xx(:), yy(:), beta0(:,:)
  REAL(kind=dp) :: x, y, z
    
  INTEGER :: nx, ny, i, j
  LOGICAL :: FirstTimeBeta=.true.

  SAVE xx, yy, beta0, nx, ny
  SAVE FirstTimeBeta

  IF (FirstTimeBeta) THEN

    FirstTimeBeta=.False.

    ! open file
    OPEN(10,file='inputs/beta0.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(beta0(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i),yy(j),beta0(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)
  
  coeff = LinearInterp(beta0, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION SSAViscosity( Model, nodenumber, dumy) RESULT(eta) !
!------------------------------------------------------------------!
  USE types
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, eta, E, yearinsec
  INTEGER :: nodenumber
  REAL(kind=dp) :: LinearInterp

  REAL(kind=dp),ALLOCATABLE :: xx(:), yy(:), flowA(:,:)
  REAL(kind=dp) :: x, y, z
    
  INTEGER :: nx, ny, i, j

  LOGICAL :: FirstTimeSSAViscosity=.TRUE.

  SAVE xx,yy,flowA,nx,ny
  SAVE FirstTimeSSAViscosity

  IF (FirstTimeSSAViscosity) THEN

    FirstTimeSSAViscosity=.False.

    ! open file
    OPEN(10,file='inputs/ssa_flowA.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(flowA(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i),yy(j),flowA(i,j)
      END DO
    END DO
    CLOSE(10)
    
  END IF

  ! year in seconds for conversion
  yearinsec = 365.25*24*60*60
    
  ! Enhancement factor
  E = 3.d0

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  eta = LinearInterp(flowA, xx, yy, nx, ny, x, y)
  eta = ((2 * E * eta * yearinsec)**(-1.0/3.0d0)) * 1.0e-6    

  RETURN
END


!------------------------------------------------------------------!
FUNCTION ModelViscosity( Model, nodenumber, dumy) RESULT(viscosity) !
!------------------------------------------------------------------!
  USE types
  Use DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy, viscosity
  INTEGER :: nodenumber

  REAL(kind=dp),ALLOCATABLE :: dem(:,:,:), xx(:), yy(:)
  REAL(kind=dp) :: x, y, z, zs , zb, dz
  REAL(kind=dp) :: yearinsec, E, alpha, eta
  INTEGER :: nx, ny, nz, k, i, j, ExtrudeLevels
  REAL(kind=dp) :: LinearInterp, zsIni, zbIni

  TYPE(Variable_t), POINTER :: dSVariable
  INTEGER, POINTER :: dSPerm(:) 
  REAL(KIND=dp), POINTER :: dSValues(:)

  LOGICAL :: Firsttime=.true., Found

  SAVE dem, xx, yy, nx, ny, nz
  SAVE Firsttime

  IF (Firsttime) THEN

    Firsttime=.False.

    ! open file
    OPEN(10,file='inputs/modelA.xyz')
    READ(10,*) nx
    READ(10,*) ny
    READ(10,*) nz
      
    ALLOCATE( xx(nx), yy(ny))
    ALLOCATE( dem(nx,ny,nz))

    DO i=1,nx
      DO j=1,ny
        READ(10, *) xx(i), yy(j), dem(i,j,:)
      END DO
    END DO
    CLOSE(10)
      
  END IF

  ! Get number of extruded levels
  ExtrudeLevels = GetInteger(Model % Simulation,'Extruded Mesh Levels',Found)

  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)
  z = Model % Mesh % Nodes % z (nodenumber)

  zs = zsIni( Model, nodenumber, dumy )
  zb = zbIni( Model, nodenumber, dumy )

  ! year in seconds for conversion		
  yearinsec=365.25d0*24*60*60

  ! Enhanced factor
  E = 3.d0

  ! Find which vertical layer the current point belongs to
  dz = (zs - zb) / (nz - 1)
  k = int( (z-zb) / dz)+1
  IF (k < 0) THEN
    PRINT *,'k',k,z,zb,zs,dz
  END IF
    
  ! Interpolate the value of the temperature from nearby points in
  ! the layers above and below it
  alpha = (z - (zb + (k - 1) * dz)) / dz

  IF (k >= ExtrudeLevels) THEN
    k = ExtrudeLevels
    eta = LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)
  ELSE
    IF (k <= 1) THEN
      k = 1
      eta = LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)
    ELSE
      eta = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) &
          + alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
    END IF  
  END IF  

  ! Get the viscosity in the correct units
  viscosity = ((2*E*eta*yearinsec)**(-1.0/3.0d0))*1.0e-6

  RETURN
END


!------------------------------------------------------------------!
FUNCTION ModelTemperature( Model, nodenumber, dumy) RESULT(T) !
!------------------------------------------------------------------!
    USE types
    Use DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    REAL(kind=dp) :: dumy,T
    INTEGER :: nodenumber

    REAL(kind=dp),ALLOCATABLE :: dem(:,:,:), xx(:), yy(:)
    REAL(kind=dp) :: x, y, z, zs , zb, dz
    REAL(kind=dp) :: alpha
    INTEGER :: nx, ny, nz, k, i, j, Timestep, TimestepInit, ExtrudeLevels
    REAL(kind=dp) :: LinearInterp, zsIni, zbIni
    TYPE(Variable_t), POINTER :: dSVariable, TimestepVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)
    LOGICAL :: Found

    LOGICAL :: Firsttime=.TRUE.
    LOGICAL :: NotMapped=.FALSE.
    LOGICAL :: Debug=.FALSE.

    SAVE dem, xx, yy, nx, ny, nz, TimestepInit
    SAVE Firsttime

    IF (Firsttime) THEN

      Firsttime=.False.

      ! open file
      IF (Debug) THEN
        PRINT *,'MODELT: Reading data file'
      END IF
      OPEN(10,file='inputs/modelT.xyz')
      READ(10,*) nx
      READ(10,*) ny
      READ(10,*) nz
      
      IF (Debug) THEN
        PRINT *,'MODELT: Allocating file'
      END IF
      ALLOCATE( xx(nx), yy(ny))
      ALLOCATE( dem(nx,ny,nz))

      DO i = 1, nx
        DO j = 1, ny
          READ(10, *) xx(i), yy(j), dem(i,j,:)
        END DO
      END DO
      CLOSE(10)
      
      TimestepVariable => VariableGet( Model % Variables,'Timestep')
      TimestepInit=TimestepVariable % Values(1)
      
    END IF

    ! Get number of extruded levels
    ExtrudeLevels = GetInteger(Model % Simulation,'Extruded Mesh Levels',Found)

    ! Get coordinates
    x = Model % Mesh % Nodes % x (nodenumber)
    y = Model % Mesh % Nodes % y (nodenumber)
    z = Model % Mesh % Nodes % z (nodenumber)

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )

    ! On the first iteration, we still have z mapped from 0 to 1, so we need to 
    ! check to make sure that it isn't the first iteration. If it is, we just 
    ! set the temperature to a default of -10 deg C.
    TimestepVariable => VariableGet( Model % Variables,'Timestep')
    Timestep = TimestepVariable % Values(1)
    IF (Timestep == TimestepInit) THEN
      IF (z <= 1.0) THEN
        IF (z >= 0.0) THEN          
          NotMapped = .true.
        END IF
      END IF
    END IF
    IF (NotMapped) THEN
      T = -10.0d0
    ELSE
      ! Find which vertical layer the current point belongs to
      dz = (zs - zb) / (nz - 1)
      k = int( (z-zb) / dz)+1
    
      ! Interpolate the value of the temperature from nearby points in
      ! the layers above and below it
      alpha = (z - (zb + (k - 1) * dz)) / dz
      IF (alpha < 0) THEN
        alpha = 0.0d0
      END IF

      IF (Debug) THEN
        PRINT *,'MODELT: k,alpha=',k,alpha
      END IF

      IF (k >= ExtrudeLevels) THEN
        k = ExtrudeLevels
        T = LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) 
      ELSE 
        IF (k <= 1) THEN
          k = 1
          T = LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)
        ELSE
          T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)+ alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
        END IF
      END IF
      
      ! In case we have restarted the file, we don't want to later end up 
      ! with this timestep
      TimestepInit = 0
    
    END IF
    
    RETURN
END

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

  REAL(kind=dp),ALLOCATABLE :: xx(:),yy(:),Tgrid(:,:)
  REAL(kind=dp) :: x,y,z
    
  INTEGER :: nx,ny
  INTEGER :: i,j
  LOGICAL :: FirstTimeTs=.TRUE.

  SAVE xx,yy,Tgrid,nx,ny
  SAVE FirstTimeTs

  IF (FirstTimeTs) THEN

    FirstTimeTs=.False.
 
    ! open file
    OPEN(10,file='inputs/t2m.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(Tgrid(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i), yy(j), Tgrid(i,j)
      END DO
    END DO
    CLOSE(10)

  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  Ts = LinearInterp(Tgrid, xx, yy, nx, ny, x, y)

  RETURN
END

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

  REAL(kind=dp),ALLOCATABLE :: xx(:),yy(:),smbgrid(:,:)
  REAL(kind=dp) :: x,y,z
    
  INTEGER :: nx,ny
  INTEGER :: i,j
  
  LOGICAL :: FirstTimea=.true.

  SAVE xx,yy,smbgrid,nx,ny
  SAVE FirstTimea

  IF (FirstTimea) THEN

    FirstTimea=.False.

    ! open file
    OPEN(10,file='inputs/smb.xy')
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(smbgrid(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i),yy(j),smbgrid(i,j)
      END DO
    END DO
    CLOSE(10)
    
  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  a = LinearInterp(smbgrid, xx, yy, nx, ny, x, y)

  RETURN
END

FUNCTION TimeVaryingAccumulation( Model, nodenumber, dumy) RESULT(a) !
!------------------------------------------------------------------!
  USE types
  USE DefUtils
  IMPLICIT NONE
  TYPE(Model_t) :: Model
  REAL(kind=dp) :: dumy,a
  TYPE(Variable_t), POINTER :: TimestepVar
  INTEGER :: nodenumber, Timestep, TimestepLast=0
  REAL(kind=dp) :: LinearInterp
  CHARACTER(LEN=MAX_NAME_LEN) :: FileName
  
  REAL(kind=dp),ALLOCATABLE :: xx(:), yy(:), smbgrid(:,:)
  REAL(kind=dp) :: x,y
    
  INTEGER :: nx, ny, i, j

  LOGICAL :: FirstTime=.TRUE.

  SAVE xx,yy,smbgrid,nx,ny
  SAVE FirstTime, TimestepLast

  TimestepVar => VariableGet( Model % Variables,'Timestep')
  Timestep = TimestepVar % Values(1)

  IF ((FirstTime) .OR. (TimeStep /= TimeStepLast)) THEN
    !IF (.NOT. FirstTime) THEN
    DEALLOCATE(xx,yy,smbgrid)
    !END IF
    FirstTime=.False.

    WRITE (FileName, "(A10,I4.4,A3)") "inputs/smb", Timestep,".xy"
    PRINT *,FileName

    ! open file
    OPEN(10,file=FileName)
    READ(10,*) nx
    READ(10,*) ny
    ALLOCATE(xx(nx),yy(ny))
    ALLOCATE(smbgrid(nx,ny))
    DO i=1,nx
      DO j=1,ny
        READ(10,*) xx(i),yy(j),smbgrid(i,j)
      END DO
    END DO
    CLOSE(10)
    
    TimestepLast = Timestep
    
  END IF

  ! position current point
  x = Model % Mesh % Nodes % x (nodenumber)
  y = Model % Mesh % Nodes % y (nodenumber)

  a = LinearInterp(smbgrid, xx, yy, nx, ny, x, y)

  RETURN
END

!------------------------------------------------------------------!
FUNCTION IceDivideTemperature( Model, nodenumber, dumy) RESULT(T) !
!------------------------------------------------------------------!
    USE types
    Use DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    REAL(kind=dp) :: dumy,T
    INTEGER :: nodenumber

    REAL(kind=dp),ALLOCATABLE :: dem(:,:,:), xx(:), yy(:)
    REAL(kind=dp) :: x, y, z, zs , zb, dz
    REAL(kind=dp) :: alpha
    INTEGER :: nx, ny, nz, k, i, j, Timestep, TimestepInit
    REAL(kind=dp) :: LinearInterp, zsIni, zbIni

    TYPE(Variable_t), POINTER :: dSVariable, TimestepVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)

    LOGICAL :: Firsttime=.true.
    LOGICAL :: NotMapped=.false.

    SAVE dem,xx,yy,nx,ny,nz,TimestepInit
    SAVE Firsttime

    IF (Firsttime) THEN

      Firsttime=.False.

      ! open file
      OPEN(10,file='inputs/tsteady_icedivide.xyz')
      READ(10,*) nx
      READ(10,*) ny
      READ(10,*) nz
      
      ALLOCATE( xx(nx), yy(ny))
      ALLOCATE( dem(nx,ny,nz))

      DO i = 1, nx
        DO j = 1, ny
          READ(10,*) xx(i), yy(j), dem(i, j, :)
        END DO
      END DO
      CLOSE(10)
      
      TimestepVariable => VariableGet( Model % Variables,'Timestep')
      TimestepInit=TimestepVariable % Values(1)
      
    END IF

    ! Get coordinates
    x = Model % Mesh % Nodes % x (nodenumber)
    y = Model % Mesh % Nodes % y (nodenumber)
    z = Model % Mesh % Nodes % z (nodenumber)

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )

    ! On the first iteration, we still have z mapped from 0 to 1, so we need to check to make
    ! sure that it isn't the first iteration. If it is, we just set the temperature to a 
    ! default of -20 deg C.
    TimestepVariable => VariableGet( Model % Variables,'Timestep')
    Timestep = TimestepVariable % Values(1)
    IF (Timestep == TimestepInit) THEN
      IF (z <= 1.0) THEN
        IF (z >= 0.0) THEN
          NotMapped=.true.         
        END IF
      END IF
    ELSE
      NotMapped=.false.
    END IF

    IF (NotMapped) THEN
      T = 273.15-20.0d0
    ELSE
      ! Find which vertical layer the current point belongs to
      dz = (zs - zb) / (nz - 1)
      k = int( (z-zb) / dz)+1
    
      ! Interpolate the value of the temperature from nearby points in
      ! the layers above and below it
      alpha = (z - (zb + (k - 1) * dz)) / dz
      T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) + alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)

      ! In case we have restarted the file, we DOn't want to later end up with this timestep
      TimestepInit = 0
    END IF
    
    RETURN
END


!------------------------------------------------------------------!
FUNCTION UModel( Model, nodenumber, dumy) RESULT(T) !
!------------------------------------------------------------------!
    USE types
    Use DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    REAL(kind=dp) :: dumy,T
    INTEGER :: nodenumber

    REAL(kind=dp),ALLOCATABLE :: dem(:,:,:), xx(:), yy(:)
    REAL(kind=dp) :: x, y, z, zs , zb, dz
    REAL(kind=dp) :: alpha
    INTEGER :: nx, ny, nz, k, i, j, Timestep, TimestepInit
    REAL(kind=dp) :: LinearInterp, zsIni, zbIni

    TYPE(Variable_t), POINTER :: dSVariable, TimestepVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)

    LOGICAL :: Firsttime=.TRUE.
    LOGICAL :: NotMapped=.FALSE.

    SAVE dem, xx, yy, nx, ny, nz, TimestepInit
    SAVE Firsttime

    IF (Firsttime) THEN

      Firsttime=.FALSE.

      ! open file
      OPEN(10,file='inputs/modelU.xyz')
      READ(10,*) nx
      READ(10,*) ny
      READ(10,*) nz
      
      ALLOCATE( xx(nx), yy(ny))
      ALLOCATE( dem(nx,ny,nz))

      DO i = 1, nx
        DO j = 1, ny
          READ(10, *) xx(i), yy(j), dem(i, j, :)
        END DO
      END DO
      CLOSE(10)
      
      TimestepVariable => VariableGet( Model % Variables,'Timestep')
      TimestepInit = TimestepVariable % Values(1)
      
    END IF

    ! Get coordinates
    x = Model % Mesh % Nodes % x (nodenumber)
    y = Model % Mesh % Nodes % y (nodenumber)
    z = Model % Mesh % Nodes % z (nodenumber)

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )

    ! On the first iteration, we still have z mapped from 0 to 1, so we need to 
    ! check to make sure that it isn't the first iteration. If it is, we just 
    ! set the temperature to a default of -10 deg C.
    TimestepVariable => VariableGet( Model % Variables,'Timestep')
    Timestep = TimestepVariable % Values(1)
    IF (Timestep == TimestepInit) THEN
      IF (z <= 1.0) THEN
        IF (z >= 0.0) THEN          
          NotMapped = .true.
        END IF
      END IF
    END IF
    IF (NotMapped) THEN
      T = 0
    ELSE
      ! Find which vertical layer the current point belongs to
      dz = (zs - zb) / (nz - 1)
      k = int( (z-zb) / dz)+1
    
      ! Interpolate the value of the temperature from nearby points in
      ! the layers above and below it
      alpha = (z - (zb + (k - 1) * dz)) / dz
      IF (alpha < 0) THEN
        alpha = 0.0d0
      END IF
      IF (k == 10) THEN
        T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) 
      ELSE
        T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)+ alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
      END IF
      
      ! In case we have restarted the file, we DOn't want to later end up 
      ! with this timestep
      TimestepInit = 0
    
    END IF
    
    RETURN
END

!------------------------------------------------------------------!
FUNCTION VModel( Model, nodenumber, dumy) RESULT(T) !
!------------------------------------------------------------------!
    USE types
    Use DefUtils
    IMPLICIT NONE
    TYPE(Model_t) :: Model
    REAL(kind=dp) :: dumy,T
    INTEGER :: nodenumber

    REAL(kind=dp),ALLOCATABLE :: dem(:,:,:), xx(:), yy(:)
    REAL(kind=dp) :: x, y, z, zs , zb, dz
    REAL(kind=dp) :: alpha
    INTEGER :: nx, ny, nz, k, i, j, Timestep, TimestepInit
    REAL(kind=dp) :: LinearInterp, zsIni, zbIni

    TYPE(Variable_t), POINTER :: dSVariable, TimestepVariable
    INTEGER, POINTER :: dSPerm(:) 
    REAL(KIND=dp), POINTER :: dSValues(:)


    LOGICAL :: Firsttime=.TRUE.
    LOGICAL :: NotMapped=.FALSE.

    SAVE dem,xx,yy,nx,ny,nz,TimestepInit
    SAVE Firsttime

    IF (Firsttime) THEN

      Firsttime=.FALSE.

      ! open file
      OPEN(10,file='inputs/modelV.xyz')
      READ(10,*) nx
      READ(10,*) ny
      READ(10,*) nz
      
      ALLOCATE(xx(nx), yy(ny))
      ALLOCATE(dem(nx, ny, nz))

      DO i = 1, nx
        DO j = 1, ny
          READ(10, *) xx(i), yy(j), dem(i, j, :)
        END DO
      END DO
      CLOSE(10)
      
      TimestepVariable => VariableGet( Model % Variables,'Timestep')
      TimestepInit = TimestepVariable % Values(1)
      
    END IF

    ! Get coordinates
    x = Model % Mesh % Nodes % x (nodenumber)
    y = Model % Mesh % Nodes % y (nodenumber)
    z = Model % Mesh % Nodes % z (nodenumber)

    zs = zsIni( Model, nodenumber, dumy )
    zb = zbIni( Model, nodenumber, dumy )

    ! On the first iteration, we still have z mapped from 0 to 1, so we need to 
    ! check to make sure that it isn't the first iteration. If it is, we just 
    ! set the temperature to a default of -10 deg C.
    TimestepVariable => VariableGet( Model % Variables,'Timestep')
    Timestep = TimestepVariable % Values(1)
    IF (Timestep == TimestepInit) THEN
      IF (z <= 1.0) THEN
        IF (z >= 0.0) THEN          
          NotMapped = .true.
        END IF
      END IF
    END IF
    IF (NotMapped) THEN
      T = 0
      print *,'First iteration, z is not mapped ',z
    ELSE
      ! Find which vertical layer the current point belongs to
      dz = (zs - zb) / (nz - 1)
      k = int( (z-zb) / dz)+1
    
      ! Interpolate the value of the temperature from nearby points in
      ! the layers above and below it
      alpha = (z - (zb + (k - 1) * dz)) / dz
      IF (alpha < 0) THEN
        alpha = 0.0d0
      END IF
      IF (k == 10) THEN
        T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y) 
      ELSE
        T = (1 - alpha) * LinearInterp(dem(:,:,k), xx, yy, nx, ny, x, y)+ alpha * LinearInterp(dem(:,:,k+1), xx, yy, nx, ny, x, y)
      END IF
      
      ! In case we have restarted the file, we don't want to later end up 
      ! with this timestep
      TimestepInit = 0
    
    END IF
    
    RETURN
END




!------------------------------------------------------------------!
include 'Interp.f90' !
!------------------------------------------------------------------!
