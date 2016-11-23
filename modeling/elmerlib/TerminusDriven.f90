SUBROUTINE ChangeMesh(Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
   USE DefUtils
   USE CRSMatrix
   USE GeneralUtils
   USE ElementDescription
   USE MeshUtils  
   USE InterpVarToVar

   IMPLICIT NONE
!------------------------------------------------------------------------------
   TYPE(Solver_t) :: Solver
   TYPE(Model_t) :: Model
   REAL(KIND=dp) :: dt
   LOGICAL :: TransientSimulation

!------------------------------------------------------------------------------
! Local variables
!------------------------------------------------------------------------------

  TYPE(Mesh_t), POINTER :: OldMesh,MapNewMesh,ExtrudedNewMesh,NewMesh
  TYPE(Variable_t), POINTER :: Var,VarNew
  TYPE(Nodes_t), POINTER :: OldNodes, NewNodes
  INTEGER :: ExtrudeLevels, i, j, k
  INTEGER :: Def_Dofs(10,6)
  INTEGER, POINTER :: FieldPerm(:)
  REAL(KIND=dp), POINTER :: Field(:)
  CHARACTER(LEN=MAX_NAME_LEN) :: Name, OldMeshName, NewMeshName, SolverName, VarName
  LOGICAL :: Found

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
  OldMesh => GetMesh()
  OldMeshName = OldMesh % Name
  OldNodes => OldMesh % Nodes
  !Var => VariableGet( Mesh % Variables, GetString(GetSolverParams(),'Interpolant'))

  ! Get new mesh
  !Def_Dofs = -1;
  !DO k=1,Model % NumberOfSolvers
  !  DO i=1,6
  !    DO j=1,8
  !      Def_Dofs(j,i) = MAXVAL(Model % Solvers(k) % Def_Dofs(j,:,i))
  !      print *,Def_Dofs(j,i)
  !    END DO
  !  END DO
  !END DO

  NewMeshName = "mesh2"
  MapNewMesh => LoadMesh2( Model, NewMeshName, NewMeshName, &
       .FALSE., Parenv % PEs, ParEnv % myPE) ! May need to adjust parameters to account for parallel mesh
  MapNewMesh % Name = NewMeshName
  
  ! Extrude new mesh and map coordinates to new mesh
  ExtrudeLevels = GetInteger(Model % Simulation,'Extruded Mesh Levels',Found)
  IF (Found) THEN
    IF(ExtrudeLevels>1) THEN
      CALL Info(SolverName,'Extruding new mesh')
      ExtrudedNewMesh => MeshExtrude(MapNewMesh, ExtrudeLevels-2)
     DO i=1,Model % NumberOfSolvers
       IF(ASSOCIATED(Solver % Mesh, Model % Meshes)) THEN
         Model % Solvers(i) % Mesh => ExtrudedNewMesh
       END IF  
     END DO
    END IF
  END IF

  IF(ASSOCIATED(Model % Mesh, OldMesh)) THEN
    Model % Mesh => ExtrudedNewMesh
  END IF
  IF(ASSOCIATED(Model % Variables, OldMesh % Variables)) THEN
    Model % Variables => ExtrudedNewMesh % Variables
  END IF

  ExtrudedNewMesh % Next => Model % Meshes % Next
  Model % Meshes => ExtrudedNewMesh

  Var => VariableGet( OldMesh % Variables, GetString(GetSolverParams(),'Interpolant'))


END SUBROUTINE ChangeMesh
