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