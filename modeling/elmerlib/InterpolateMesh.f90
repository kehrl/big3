! ******************************************************************************
! *
! *  Laura Kehrl
! *  University of Washington
! *  Dec. 1, 2014
! * 
! *****************************************************************************

SUBROUTINE RetreatMesh( Model,Solver,dt,TransientSimulation,Retreated )

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
  REAL(KIND=dp) :: x, y, z, xdelta, Retreated
  

!-----------------------------------------------------------------------------
	!print *,'',Retreated
  
END SUBROUTINE RetreatMesh

SUBROUTINE InterpolateMesh( Model,Solver,dt,TransientSimulation )
  USE DefUtils

  IMPLICIT NONE

  TYPE(Solver_t) :: Solver
  TYPE(Model_t) :: Model

  REAL(KIND=dp) :: dt
  LOGICAL :: TransientSimulation
  TYPE(Variable_t), POINTER :: Var

  TYPE(Mesh_t), POINTER :: Mesh, EvalMesh
  
  !Mesh => GetMesh()
  !EvalMesh => GetMesh()
  !CALL InterpolateMeshToMesh( Mesh, EvalMesh, Mesh % Variables, EvalMesh % Variables, .FALSE. )
  !CALL InterpolateMeshToMesh( Mesh, Mesh, Mesh % Variables, Mesh % Variables, .FALSE. )
  
  Mesh => GetMesh()
  Var => VariableGet( Mesh % Variables, 'ceq residual')
  Var => VariableGet( Mesh % Variables, 'Velocity 1')
  Var => VariableGet( Mesh % Variables, 'Velocity 2')
  Var => VariableGet( Mesh % Variables, 'Stress Vector 1')
  Var => VariableGet( Mesh % Variables, 'Stress Vector 2')
  Var => VariableGet( Mesh % Variables, 'Stress Vector 3')
  Var => VariableGet( Mesh % Variables, 'sxx')
  Var => VariableGet( Mesh % Variables, 'sxy')
  Var => VariableGet( Mesh % Variables, 'szz')
  Var => VariableGet( Mesh % Variables, 'syy')
  Var => VariableGet( Mesh % Variables, 'Pressure')
 
  

END SUBROUTINE InterpolateMesh