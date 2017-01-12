SUBROUTINE ChangeMesh(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  USE DefUtils
  USE CRSMatrix
  USE GeneralUtils
  USE ElementDescription
  USE MeshUtils  
  USE InterpVarToVar
  USE SolverUtils
  USE ModelDescription

  USE LoadMod

  IMPLICIT NONE
!------------------------------------------------------------------------------
  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model
  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: OldMesh => NULL(), NewMesh => NULL(), Mesh, sMesh
  TYPE(Variable_t), POINTER :: Var,NewVar,Var2
  TYPE(Nodes_t), POINTER :: OldNodes, NewNodes
  !TYPE(Element_t), POINTER :: NewElements,OldElements
  TYPE( Matrix_t ), POINTER :: NewMatrix
  INTEGER :: i, j, k, n
  INTEGER :: Def_Dofs(10,6)
  INTEGER, POINTER :: Permutation(:)
  REAL(KIND=dp), POINTER :: Field(:), Work(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, OldMeshName, NewMeshName, SolverName, VarName
  LOGICAL :: Found
  LOGICAL :: BandwidthOptimized, GlobalBubbles
  REAL(kind=dp) :: t

  INTERFACE
     SUBROUTINE InterpolateMeshToMesh( OldMesh, NewMesh, OldVariables, &
          NewVariables, UseQuadrantTree, Projector, MaskName, UnfoundNodes )
       !------------------------------------------------------------------------------
       USE Lists
       USE SParIterComm
       USE Interpolation
       USE CoordinateSystems
       !-------------------------------------------------------------------------------
       TYPE(Mesh_t), TARGET  :: OldMesh, NewMesh
       TYPE(Variable_t), POINTER, OPTIONAL :: OldVariables, NewVariables
       LOGICAL, OPTIONAL :: UseQuadrantTree
       LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:)
       TYPE(Projector_t), POINTER, OPTIONAL :: Projector
       CHARACTER(LEN=*),OPTIONAL :: MaskName
     END SUBROUTINE InterpolateMeshToMesh
  END INTERFACE
  
  SolverName = "ChangeMesh"
  
  ! Get current mesh
  OldMesh => Solver % Mesh
  OldMeshName = OldMesh % Name
  OldNodes => OldMesh % Nodes
  !Var => VariableGet( Mesh % Variables, GetString(GetSolverParams(),'Interpolant'))

  ! Get new mesh
  NewMeshName = "mesh1"
  NewMesh => LoadMesh2( Model, NewMeshName, NewMeshName, &
       .FALSE., Parenv % PEs, ParEnv % myPE) ! May need to adjust parameters to account for parallel mesh
  NewMesh % Name = NewMeshName
  NewNodes => NewMesh % Nodes
  !NewElements => NewMesh % Elements
  
  ! Create variables in new mesh
  !NULLIFY( NewMesh % Variables )
  !DO i=1,99
  !  WRITE( Name, '(A,I0)') 'Variable ',i
  !  VarName = ListGetString( Solver % Values,TRIM(Name),Found)
  !  IF(.NOT. Found ) EXIT

  !   Var => VariableGet( OldMesh % Variables, TRIM( VarName), ThisOnly = .TRUE.) 
  !    IF( ASSOCIATED( Var) ) THEN
  !      NULLIFY( Field, FieldPerm )
  !      ALLOCATE( Field(NewNodes % NumberofNodes), FieldPerm(NewNodes % NumberofNodes) )  

        !TODO InterpolateMeshToMesh should be made to report missing points in 
        !interpolation, as InterpVarToVar currently does.
  !      Field = -HUGE(Field)
  !      FieldPerm = 0

  !      CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, TRIM( VarName), 1, Field, FieldPerm)

  !      IF(ASSOCIATED(Var % PrevValues)) THEN
  !        VarNew => VariableGet( NewMesh % Variables, TRIM( VarName), ThisOnly = .TRUE. )
          !ALLOCATE( VarNew % PrevValues ( SIZE( Var % PrevValues,1), &
          !                                    SIZE( Var % PrevValues,2) ) )
  !      END IF
  !      CALL Info(SolverName,'Created variable: '//TRIM( VarName ) )
  !    ELSE
  !      WRITE(Message,'(a,a,a)') "Requested variable: ", TRIM(VarName), " but it wasn't found!"
  !      CALL Warn(SolverName, Message)
  !    END IF
  !END DO
  
  ! Interpolate variables between meshes
  !CALL InterpolateMeshToMesh( OldMesh, NewMesh, OldMesh % Variables,NewMesh % Variables, .TRUE. )
  !CALL Info('ChangeMesh','Interpolation done')

  !OldNodes % x = NewNodes % x
  !OldNodes % y = NewNodes % y
  !OldNodes % z = NewNodes % z
  !OldNodes % NumberofNodes = NewNodes % NumberofNodes
  !Model % Mesh % Elements => NewMesh % Elements
  !Model % Mesh % Faces => NewMesh % Faces
  !Model % Mesh % Edges => NewMesh % Edges
  !OldMesh % ParallelInfo = NewMesh % ParallelInfo
  !Model % Mesh % RootQuadrant => NewMesh % RootQuadrant

  !Var => NewMesh % Variables
  !DO WHILE(ASSOCIATED(Var))
  !   Var2 => VariableGet(OldMesh % Variables, TRIM(Var % Name), ThisOnly = .TRUE.)
  !   IF(.NOT. ASSOCIATED(Var2)) CALL Fatal(SolverName, "Error while copying back variables")
  !   CALL Info(SolverName,'Copying back variable: '//TRIM( Var % Name ) )
  !   Var2 % Values = Var % Values
  !   PRINT *,'Range',MINVAL( Var2 % Values), MAXVAL( Var2 % Values )

  !   Var => Var % Next
  !END DO

  !DO i=1,99
  !   WRITE( Name, '(A,I0)') 'Nullify ',i
  !   VarName = ListGetString( Solver % Values,TRIM(Name),Found)
  !   IF(.NOT. Found ) EXIT
  !   Var2 => VariableGet( OldMesh % Variables, TRIM( VarName), ThisOnly = .TRUE. ) 
  !   IF( ASSOCIATED( Var2) ) THEN
  !      CALL Info(SolverName,'Zeroing variable: '//TRIM( VarName ) )
  !      Var2 % Values = 0.0_dp
  !   ELSE
  !      WRITE(Message,'(a,a,a)') "Requested nullified variable: ", &
  !           TRIM(VarName), " but it wasn't found!"
  !      CALL Warn(SolverName, Message)
  !   END IF
  !END DO

  !OldMesh % Changed = .TRUE.
  
  !   Add the new mesh to the global list of meshes:
  !   ----------------------------------------------
  NewMesh % Next   => Model % Meshes 
  Model % Meshes   => NewMesh
  OldMesh % Child  => NewMesh
  NewMesh % Parent => OldMesh
  NewMesh % Child => NULL()

  NewMesh % MaxBDOFs = OldMesh % MaxBDOFs


  !   Initialize local variables for the new mesh:
  !   --------------------------------------------
  NULLIFY( NewMesh % Variables )

  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 1', 1, NewMesh % Nodes % x )

  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 2', 1, NewMesh % Nodes % y )

  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
         'Coordinate 3', 1, NewMesh % Nodes % z )

  !   Time must always be there:
  !   --------------------------
  Var => VariableGet( OldMesh % Variables,'Time',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Time', 1, Var % Values )

  Var => VariableGet( OldMesh % Variables,'Timestep',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Timestep', 1, Var % Values )

  Var => VariableGet( OldMesh % Variables,'Timestep size',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                     'Timestep size', 1, Var % Values )

  NewVar => VariableGet( NewMesh % Variables,'Timestep size',ThisOnly=.TRUE. )
  NewVar % PrevValues => Var % PrevValues

  Var => VariableGet( OldMesh % Variables,'Timestep interval',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Timestep interval', 1, Var % Values )

  Var => VariableGet( OldMesh % Variables,'Coupled iter',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Coupled iter', 1, Var % Values )

  Var => VariableGet( OldMesh % Variables,'Nonlin iter',ThisOnly=.TRUE. )
  CALL VariableAdd( NewMesh % Variables, NewMesh, Solver, &
                 'Nonlin iter', 1, Var % Values )

  ! Initialize the field variables for the new mesh. These are
  ! interpolated from the old meshes variables. Vector variables
  ! are in the variable lists in two ways: as vectors and as
  ! vector components. We MUST update the vectors (i.e. DOFs>1)
  ! first!!!!!
  ! -----------------------------------------------------------
  CALL SetCurrentMesh( Model, NewMesh )
  Var => OldMesh % Variables
  DO WHILE( ASSOCIATED( Var ) )
    IF ( Var % DOFs > 1 ) THEN
      NewVar => VariableGet( NewMesh % Variables,Var % Name,.FALSE. )
      k = SIZE( NewVar % Values )
      IF ( ASSOCIATED( NewVar % Perm ) ) THEN
        k = COUNT( NewVar % Perm > 0 )
      END IF
      IF ( GetVarName( NewVar ) == 'flow solution' ) THEN
        NewVar % Norm = 0.0d0
        DO i=1,NewMesh % NumberOfNodes
          DO j=1,NewVar % DOFs-1
            NewVar % Norm = NewVar % Norm + &
            NewVar % Values( NewVar % DOFs*(i-1)+j )**2
          END DO
        END DO
        NewVar % Norm = SQRT( NewVar % Norm / k )
      ELSE
        NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
      END IF
    END IF
    Var => Var % Next
  END DO

!   Second time around, update scalar variables and
!   vector components:
!   -----------------------------------------------
  Var => OldMesh % Variables
  DO WHILE( ASSOCIATED( Var ) )

    IF( SIZE( Var % Values ) == Var % DOFs ) THEN
      Var => Var % Next
      CYCLE
    END IF

    SELECT CASE( Var % Name )
    CASE( 'coordinate 1', 'coordinate 2', 'coordinate 3', 'time' , &
             'timestep', 'timestep size', 'timestep interval', &
             'coupled iter', 'nonlin iter' )
    CASE DEFAULT
      IF ( Var % DOFs == 1 ) THEN
        Found = .FALSE.
        Found = Found .OR. INDEX( Var % Name, '.error'  ) > 0
        Found = Found .OR. INDEX( Var % Name, '.eOld'   ) > 0
        Found = Found .OR. INDEX( Var % Name, '.perror' ) > 0
        IF ( Found ) THEN
          k = Solver % Variable % NameLen
          IF ( Var % Name(1:k) /= Solver % Variable % Name(1:k) ) THEN
            Var => Var % Next
            CYCLE
          END IF
        END IF

        NewVar => VariableGet( NewMesh % Variables, Var % Name, .FALSE. )
        k = SIZE( NewVar % Values )
        IF ( ASSOCIATED( NewVar % Perm ) ) THEN
          k = COUNT( NewVar % Perm > 0 )
        END IF
          NewVar % Norm = SQRT( SUM(NewVar % Values**2) / k )
        END IF
      END SELECT
      Var => Var % Next
    END DO

!-------------------------------------------------------------------    
    WRITE( Message, * )  &
     'Mesh variable update time (cpu-secs):            ',CPUTime()-t
    CALL Info( 'OldineMesh', Message, Level = 6 )
!-------------------------------------------------------------------    

!
!   Update Solver structure to use the new mesh:
!   ---------------------------------------------    
    Solver % Mesh => NewMesh
    CALL MeshStabParams( NewMesh )
!
!   Nothing computed on this mesh yet:
!   ----------------------------------
    NewMesh % SavesDone    = 0  ! start new output file
    NewMesh % OutputActive = .FALSE.
    
    NewMesh % Changed   = .TRUE.

!
!   Update the solvers variable pointer:
!   ------------------------------------    
    Solver % Variable => VariableGet( Solver % Mesh % Variables, &
            Solver % Variable % Name, ThisOnly=.TRUE. )
    !Solver % Variable % PrimaryMesh => NewMesh

!
!   Create matrix structures for the new mesh:
!   ------------------------------------------    
    t = CPUTime()

!
!   Try to account for the reordering of DOFs
!   due to bandwidth optimization:
!   -----------------------------------------
    !GlobalBubbles = ListGetLogical( Solver % Values, &
    !      'Bubbles in Global System', Found )
    !IF ( .NOT. Found ) GlobalBubbles = .FALSE.

    !BandwidthOptimized = ListGetLogical( Solver % Values, &
    !     'Optimize Bandwidth', Found )
    !IF ( .NOT. Found ) BandwidthOptimized = .FALSE.

    !IF ( BandwidthOptimized ) THEN
    !   n = NewMesh % NumberOfNodes
    !   IF ( GlobalBubbles ) &
    !      n = n + NewMesh % MaxBDOFs*NewMesh % NumberOFBulkElements
    !   CALL AllocateVector( Permutation,  n )
    !ELSE
    Permutation => Solver % Variable % Perm
    !END IF

    ! Create the CRS format matrix tables for solving the
    ! current equation on the new mesh. Also do bandwidth
    ! optimization, if requested:
    ! ----------------------------------------------------
    !NewMatrix => CreateMatrix( Model, Solver, Solver % Mesh,  &
    !     Permutation, Solver % Variable % DOFs, MATRIX_CRS, &
    !     .FALSE., ListGetString( Solver % Values, 'Equation' ), &
    !     GlobalBubbles=.FALSE. )

!    IF ( ASSOCIATED( Solver % Matrix ) ) THEN
!       CALL FreeMatrix( Solver % Matrix )
!       NULLIFY( Solver % Matrix )
!    END IF

!   Solver % Matrix % Child  => NewMatrix
!    NewMatrix % Parent => Solver % Matrix
!    NULLIFY( NewMatrix % Child )
!    Solver % Matrix => NewMatrix

!
!   Reorder the primary variable for bandwidth optimization:
!   --------------------------------------------------------
!    IF ( BandwidthOptimized ) THEN
!       n = Solver % Variable % DOFs
!       ALLOCATE( Work(SIZE(Permutation)) )
!       DO i=0,n-1
!#if 0
!          WHERE( Permutation > 0 )
!            Solver % Variable % Values(n*Permutation-i) = &
!                  Solver % Variable % Values(n*Solver % Variable % Perm-i)
!          END WHERE
 
!          IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
!             DO j=1,SIZE( Solver % Variable % PrevValues,2)
!               WHERE( Permutation > 0 )
!                 Solver % Variable % PrevValues(n*Permutation-i,j) = &
!                  Solver % Variable % PrevValues(n*Solver % Variable % Perm-i,j)
!               END WHERE
!             END DO
!          END IF
!#else
!          Work = Solver % Variable % Values(i+1::n)
!          DO j=1,SIZE(Permutation)
!             IF ( Permutation(j) > 0 ) THEN
!                Solver % Variable % Values(n*Permutation(j)-i) = &
!                    Work(Solver % Variable % Perm(j))
!             END IF
!          END DO
!          IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
!             DO j=1,SIZE( Solver % Variable % PrevValues,2)
!               Work = Solver % Variable % PrevValues(i+1::n,j)
!               DO k=1,SIZE(Permutation)
!                  IF ( Permutation(k) > 0 ) THEN
!                     Solver % Variable % PrevValues(n*Permutation(k)-i,j) = &
!                         Work(Solver % Variable % Perm(k))
!                  END IF
!               END DO
!             END DO
!          END IF
!#endif
!       END DO
!       k = SIZE(Permutation)
!       Solver % Variable % Perm(1:k) = Permutation
!       DEALLOCATE( Permutation, Work )
!    END IF

!   TODO: CreateMatrix should do these
!   -----------------------------------

    !Solver % Matrix % Lumped = ListGetLogical( Solver % Values, &
    !        'Lumped Mass Matrix', Found )

    !Solver % Matrix % Symmetric = ListGetLogical( Solver % Values, &
    !        'Linear System Symmetric', Found )

    !CALL AllocateVector( Solver % Matrix % RHS, SIZE(Solver % Variable % Values) )
    !Solver % Matrix % RHS = 0.0d0

!   Transient case additional allocations:
!   --------------------------------------
    !IF (ListGetString(Model % Simulation,'Simulation Type')=='transient') THEN
    !   n = SIZE( Solver % Variable % Values )

    !   CALL AllocateArray( Solver % Matrix % Force, n, Solver % TimeOrder+1 )
    !   Solver % Matrix % Force = 0.0d0
    !END IF

    CALL ParallelInitMatrix( Solver, Solver % Matrix )

    WRITE( Message, * ) 'Matrix structures update time (cpu-secs):        ',CPUTime()-t
    CALL Info( 'OldineMesh', Message, Level=6 )

!
!   Release previous meshes. Keep only the original mesh, and
!   the last two meshes:
!   ---------------------------------------------------------
    !Mesh => OldMesh % Parent
    !DO WHILE( ASSOCIATED(Mesh) )
    !   sMesh => Mesh % Parent
    !   IF ( Mesh % AdaptiveDepth /= 0 ) THEN
    !      IF ( ASSOCIATED( Mesh % Parent ) ) THEN
    !         Mesh % Parent % Child => Mesh % Child
    !      END IF

    !      IF ( ASSOCIATED(Mesh % Child) ) THEN
    !         Mesh % Child % Parent => Mesh % Parent
    !      END IF

    !      CALL ReleaseMesh( Mesh )

    !      Mesh % Child  => NULL()
    !      Mesh % Parent => NULL()
    !   END IF
    !   Mesh => sMesh
    !END DO

  DO i=1,Model % NumberOfSolvers
    IF(ASSOCIATED(Model % Solvers(i) % Mesh,Model % Meshes)) &
      Model % Solvers(i) % Mesh => NewMesh 
  END DO
  NewMesh % Next => Model % Meshes % Next
  Model % Meshes => NewMesh

  
END SUBROUTINE ChangeMesh
