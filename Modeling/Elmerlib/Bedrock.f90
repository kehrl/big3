SUBROUTINE Bedrock( Model,Solver,dt,TransientSimulation )

!*************************************************************************
!
! Create the bedrock 
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
  
  REAL(KIND=dp)  :: BedFromFile

  Logical :: First=.true.

!-----------------------------------------------------------------------------
  ! Bedrock variable
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
            !VariableValues(Permutation(Element % NodeIndexes(ii))) = zbed(ii) 
      ELSE IF (DIM==2) THEN        
            VariableValues(Permutation(ii)) = BedFromFile(x)
      		!print *,'cycling2'
      		!if (x > 80000) then
      		!  print *,'x,bed',x,xdelta,BedFromFile(x)
      		!end if        
      END IF

  END DO

END SUBROUTINE Bedrock

Function BedFromFile(x)
   	USE types
   	IMPLICIT NONE
   	INTEGER :: i, R7, ind
   	REAL(KIND=dp) :: x, dist, mindist
   	REAL(KIND=dp) :: BedFromFile
   	REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), zbed(:)

	Logical :: Firsttime7=.true.
	logical :: found

    SAVE xbed,ybed
    SAVE Firsttime7,R7

    if (Firsttime7) then

    	Firsttime7=.False.

   		! Load bedrock
   		OPEN(10,file="Inputs/roughbed.dat")
   		Read(10,*) R7
   		ALLOCATE(xbed(R7), ybed(R7))
   		READ(10,*)(xbed(i), ybed(i), i=1,R7)
   		CLOSE(10)
   		
    End if

   	found = .false.
   	mindist=dabs(2*(xbed(1)-xbed(2)))
   	IF (x < xbed(1)) THEN
    	ind=i
    	found = .true.
    ELSE
   		DO 20, i=1,R7
      	dist=dabs(x-xbed(i))
      		IF (dist<=mindist .and. xbed(i)<=x) THEN
        		mindist=dist
        		ind=i
        		found = .true.
      		Endif 
   		20 Enddo
   	END IF

   	IF (.not.found) THEN
      	print *, 'Could not find a suitable bed to interpolate ',x
   	ELSE
      	BedFromFile=((x-xbed(ind))/(xbed(ind+1)-xbed(ind)))*(ybed(ind+1)-ybed(ind))+ybed(ind)
   	END IF

END FUNCTION BedFromFile
