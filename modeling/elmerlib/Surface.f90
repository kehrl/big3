SUBROUTINE Surfrock( Model,Solver,dt,TransientSimulation )

!*************************************************************************
!
! Create the Surfrock 
!
!*************************************************************************

  USE DefUtils
  IMPLICIT NONE

!-----------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------
  TYPE(Element_t),POINTER :: Element
  TYPE(Variable_t), POINTER :: PointerToVariable, MeshVariable
  TYPE(Nodes_t), SAVE :: Nodes
  TYPE(ValueList_t), POINTER :: SolverParams

  INTEGER :: ii, tt, nn, jj, DIM, R, i
  INTEGER, POINTER :: Permutation(:), MeshPerm(:)

  REAL(KIND=dp), POINTER :: VariableValues(:), MeshValues(:)
  REAL(KIND=dp) :: x, y, z, xdelta
  
  REAL(KIND=dp)  :: SurfFromFile

  Logical :: First=.true.

!-----------------------------------------------------------------------------
  ! Surfrock variable
  PointerToVariable => Solver % Variable
  Permutation  => PointerToVariable % Perm
  VariableValues => PointerToVariable % Values

   ! Current value of Mesh Update 1
   MeshVariable => VariableGet( Model % Variables, 'Mesh Update 1' )
   IF (ASSOCIATED(MeshVariable)) THEN
    	MeshPerm    => MeshVariable % Perm
    	MeshValues  => MeshVariable % Values
   !ELSE
        !CALL FATAL('ZsBottomMZsIni','Could not find variable >Mesh Update 1<')
   END IF

  DIM = CoordinateSystemDimension()


  DO ii = 1, Model % NumberofNodes
      !IF ( Permutation(ii) == 0 ) CYCLE
    
    	x = Model % Nodes % x (ii) 
    	xdelta = MeshValues(MeshPerm(ii))
		x = x + xdelta
		
      IF (DIM==3) THEN
            !VariableValues(Permutation(Element % NodeIndexes(ii))) = zSurf(ii) 
      ELSE IF (DIM==2) THEN        
            VariableValues(Permutation(ii)) = SurfFromFile(x)
      		!print *,'cycling2'
      		!if (x > 80000) then
      		!  print *,'x,Surf',x,xdelta,SurfFromFile(x)
      		!end if        
      END IF

  END DO

END SUBROUTINE Surfrock

Function SurfFromFile(x)
   	USE types
   	IMPLICIT NONE
   	INTEGER :: i, R7, ind
   	REAL(KIND=dp) :: x, dist, mindist
   	REAL(KIND=dp) :: SurfFromFile
   	REAL(KIND=dp), ALLOCATABLE :: xSurf(:), ySurf(:), zSurf(:)

	Logical :: Firsttime7=.true.
	logical :: found

    SAVE xSurf,ySurf
    SAVE Firsttime7,R7

    if (Firsttime7) then

    	Firsttime7=.False.

   		! Load Surfrock
   		OPEN(10,file="Inputs/surf.dat")
   		Read(10,*) R7
   		ALLOCATE(xSurf(R7), ySurf(R7))
   		READ(10,*)(xSurf(i), ySurf(i), i=1,R7)
   		CLOSE(10)
   		
    End if

   	found = .false.
   	mindist=dabs(2*(xSurf(1)-xSurf(2)))
   	IF (x < xSurf(1)) THEN
    	ind=i
    	found = .true.
    ELSE
   		DO 20, i=1,R7
      	dist=dabs(x-xSurf(i))
      		IF (dist<=mindist .and. xSurf(i)<=x) THEN
        		mindist=dist
        		ind=i
        		found = .true.
      		Endif 
   		20 Enddo
   	END IF

   	IF (.not.found) THEN
      	print *, 'Could not find a suitable Surf to interpolate ',x
   	ELSE
      	SurfFromFile=((x-xSurf(ind))/(xSurf(ind+1)-xSurf(ind)))*(ySurf(ind+1)-ySurf(ind))+ySurf(ind)
   	END IF

END FUNCTION SurfFromFile