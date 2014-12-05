!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Olivier Gagliardini, Ga¨el Durand
! *  Modified by Laura Kehrl to account for horizontal mesh movement 
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! * 
! *****************************************************************************
!> ZsMZsIni for variable Zs                   
!> ZsTopMZsIni for variable Zs Top         
!> ZsBottomMZsIni for variable Zs Bottom
!> DyMDyIni for any FS variable name           
FUNCTION ZsIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsIni','Could not find variable >Zs<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsIni

FUNCTION ZsMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsTopMZsIni','Could not find variable >Zs<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs -  Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsMZsIni

!--------------------------------------------------------------------------------

FUNCTION ZsTopIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Top')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsIni','Could not find variable >Zs Top<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsTopIni

FUNCTION ZsTopMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: mu,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Top')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsTopMZsIni','Could not find variable >Zs Top<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

      mu =  Zs -  Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsTopMZsIni

!--------------------------------------------------------------------

FUNCTION ZsBottomIni ( Model, nodenumber, x) RESULT(Zs)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol
   INTEGER, POINTER :: ZsPerm(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim
   REAL(KIND=dp) :: x,   Zs       
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:)       
   LOGICAL :: FirstTime=.True. 

   SAVE FirstTime
   SAVE Zs0 

   ZsSol => VariableGet( Model % Variables, 'Zs Bottom')
   IF (ASSOCIATED(ZsSol)) THEN
        ZsPerm => ZsSol % Perm
   ELSE
        CALL FATAL('ZsBottomIni','Could not find variable >Zs Bottom<')
   END IF

   IF (FirstTime) THEN
        FirstTime = .False.
        dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 ) 
        ALLOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberOfNodes
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF

       Zs = Zs0(ZsPerm(nodenumber)) 

END FUNCTION ZsBottomIni

FUNCTION ZsBottomMzsIni ( Model, nodenumber, Zs) RESULT(mu)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Variable_t), POINTER :: ZsSol, MeshVariable, GroundedMaskVariable 
   INTEGER, POINTER :: ZsPerm(:), MeshPerm(:), GroundedMaskPerm(:) 
   REAL(KIND=dp), POINTER :: MeshValues(:), GroundedMaskValues(:)
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber,  NMAX, i, dim, R, ind
   REAL(KIND=dp) :: mu,   Zs, x , ratio, ZsNew, mindist, dist, junk      
   REAL(KIND=dp), ALLOCATABLE :: Zs0(:), XIni(:), xbed(:), ybed(:)       
   LOGICAL :: FirstTime=.True. 
   LOGICAL :: found=.False.
   Real(kind=dp) :: BedFromFile

   SAVE FirstTime
   SAVE Zs0, ZsPerm, XIni, xbed, ybed 

   IF (FirstTime) THEN
        FirstTime = .False. 
 
    	! Get initial value of Zs Bottom
   		ZsSol => VariableGet( Model % Variables, 'ZsBottomIni')
   		IF (ASSOCIATED(ZsSol)) THEN
        	ZsPerm => ZsSol % Perm
   		ELSE
        	CALL FATAL('ZsBottomMZsIni','Could not find variable >Zs Bottom<')
   		END IF
   		
   		dim = CoordinateSystemDimension()
        NMAX = COUNT( ZsPerm > 0 )
        ALLOCATE(XIni(Model % NumberofNodes))
        AllOCATE(Zs0(NMAX))
        DO i = 1, Model % NumberofNodes
          XIni(i) = Model % Nodes % x (i)
          IF (ZsPerm(i)==0) CYCLE
          IF (dim==2) THEN
             Zs0(ZsPerm(i)) = Model % Nodes % y (i)
          ELSE
             Zs0(ZsPerm(i)) = Model % Nodes % z (i)
          END IF
        END DO
   END IF
   
   ! Get current value of Mesh Update 1
   MeshVariable => VariableGet( Model % Variables, 'Mesh Update 1' )
   IF (ASSOCIATED(MeshVariable)) THEN
    	MeshPerm    => MeshVariable % Perm
    	MeshValues  => MeshVariable % Values
   ELSE
        CALL FATAL('ZsBottomMZsIni','Could not find variable >Mesh Update 1<')
   END IF
   
   ! Get current value of GroundedMask
   GroundedMaskVariable => VariableGet( Model % Variables, 'GroundedMask' )
   IF (ASSOCIATED(GroundedMaskVariable)) THEN
    	GroundedMaskPerm    => GroundedMaskVariable % Perm
    	GroundedMaskValues  => GroundedMaskVariable % Values
   ELSE
        CALL FATAL('ZsBottomMZsIni','Could not find variable >Grounded Mask<')
   END IF

   IF (GroundedMaskValues(GroundedMaskPerm(nodenumber)) > -0.5_dp) THEN
		! Update so that it is equal to the bed
		! Calculate current x-coordinate
		x = XIni(nodenumber) + MeshValues(MeshPerm(nodenumber))
	  	ZsNew = BedFromFile(x)
   ELSE
   		ZsNew = Zs
   END IF
    
    mu =  ZsNew -  Zs0(ZsPerm(nodenumber)) 
      !print *,'Test me', Xini(nodenumber), x, ZsNew

END FUNCTION ZsBottomMzsIni

!---------------------------------------------------------------------------

include 'Bedrock.f90'