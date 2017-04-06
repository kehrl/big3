
!InterpolateVarToVarReduced: Subroutine to interpolate variable values on a reduced
!dimension mesh.  Either once or twice reduced (3D => line search, in the latter case)

! This subroutine is largely based on InterpolateMeshToMesh, with a few modifications
!
! The user specifies a variable 'HeightName' which is to be interpolated. This differs
! from InterpolateMeshToMesh, where all variables are interpolated. This is because
! this subroutine has been principally used for remeshing a glacier, and so only the 
! height variable is required.
! 
! Its also possible to specify a VariableList to get similar behaviour to 
! InterpolateMeshToMesh. If a VariableList is specified which has regular 
! field variables, OldNodeMask and NewNodeMask should be supplied to ensure
! only BC values are interpolated.  Alternatively, one could create a perm
! of an internal layer of an extruded mesh.
! 
! Care should be taken with supplying overly relaxed epsilon values: 
! The subroutine will naively search elements in turn with the given epsilon
! values, and so if a neighbouring element which 'almost' contains the node
! is searched first, this will be accepted by the subroutine. If a more relaxed
! epsilon search is required, one should first call this subroutine with more
! stringent epsilon values, then call again, with more relaxed epsilon, and a mask
! specifying which nodes to look for (i.e. only those not found in the previous run)
!------------------------------------------------------------------------------
MODULE InterpVarToVar

  USE DefUtils
  USE Types
  USE Interpolation

  IMPLICIT NONE

CONTAINS
  SUBROUTINE InterpolateVartoVarReduced( OldMesh, NewMesh, HeightName, HeightDimensions,&
       UnfoundNodes, OldNodeMask, NewNodeMask, OldElemMask, Variables, GlobalEps, LocalEps, NumericalEps)

    TYPE(Mesh_t), TARGET, INTENT(IN)  :: OldMesh
    TYPE(Mesh_t), TARGET, INTENT(INOUT)  :: NewMesh
    TYPE(Variable_t), POINTER, OPTIONAL :: Variables
    CHARACTER(LEN=*) :: HeightName
    INTEGER, POINTER :: HeightDimensions(:)
    REAL(KIND=dp), OPTIONAL :: GlobalEps, LocalEps, NumericalEps
    LOGICAL, POINTER, OPTIONAL :: UnfoundNodes(:), OldNodeMask(:), &
         NewNodeMask(:), OldElemMask(:)
    !------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, nVar, OldVar
    TYPE(Mesh_t), POINTER :: nMesh
    REAL(KIND=dp), POINTER :: OldHeight(:), NewHeight(:), nReal(:)
    INTEGER, POINTER :: OldPerm(:), NewPerm(:), nperm(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i, j, k, l, n, nvars, ierr, npart, nfound, proc, status(MPI_STATUS_SIZE), unfound
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: perm(:), vperm(:)
    REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
    REAL(KIND=dp), POINTER :: ElementValues(:)
    TYPE(Element_t),POINTER :: Element       
    REAL(KIND=dp), ALLOCATABLE, TARGET :: BB(:,:), nodes_x(:),nodes_y(:),nodes_z(:),vstore(:,:),&
         astore(:), xpart(:), ypart(:), zpart(:), PointLocalDistance(:), &
         SendLocalDistance(:), RecvLocalDistance(:), WorkReal(:)
    REAL(KIND=dp) :: detJ, u,v,w,s, dn, NoData = -99999.0_dp
    LOGICAL :: Found, Debug
    REAL(KIND=dp) :: myBB(6), epsBB
    LOGICAL, ALLOCATABLE :: FoundNodes(:), BetterFound(:)

    !------------------------------------------------------------------------------
    TYPE ProcRecv_t
       INTEGER :: n = 0
       REAL(KIND=dp), ALLOCATABLE :: nodes_x(:),nodes_y(:),nodes_z(:)
    END TYPE ProcRecv_t
    TYPE(ProcRecv_t),  ALLOCATABLE, TARGET :: ProcRecv(:)

    TYPE ProcSend_t
       INTEGER :: n = 0
       INTEGER, ALLOCATABLE :: perm(:)
    END TYPE ProcSend_t
    TYPE(ProcSend_t),  ALLOCATABLE :: ProcSend(:)
    !------------------------------------------------------------------------------

    Debug = .FALSE.

    ALLOCATE( FoundNodes(NewMesh % NumberOfNodes),&
         PointLocalDistance(NewMesh % NumberOfNodes))

    FoundNodes=.TRUE.
    PointLocalDistance = 0.0_dp

    IF(PRESENT(UnfoundNodes)) THEN
       IF(ASSOCIATED(UnfoundNodes)) DEALLOCATE(UnfoundNodes)
       ALLOCATE(UnfoundNodes(NewMesh % NumberOfNodes))
    END IF

    IF ( ParEnv % PEs<=1 ) THEN
       CALL InterpolateVarToVarReducedQ( OldMesh, NewMesh, HeightName, HeightDimensions, &
            FoundNodes=FoundNodes, OldNodeMask=OldNodeMask, NewNodeMask=NewNodeMask, &
            OldElemMask=OldElemMask, Variables=Variables, GlobalEps=GlobalEps, &
            LocalEps=LocalEps, NumericalEps=NumericalEps )

       IF(PRESENT(UnfoundNodes)) UnfoundNodes = .NOT. FoundNodes
       RETURN
    END IF

    CALL InterpolateVarToVarReducedQ( OldMesh, NewMesh, HeightName, HeightDimensions, &
         FoundNodes, PointLocalDistance, OldNodeMask, NewNodeMask, &
         OldElemMask, Variables, GlobalEps, LocalEps, NumericalEps )
    CALL MPI_BARRIER(ParEnv % ActiveComm, ierr)

    IF(PRESENT(UnfoundNodes)) UnfoundNodes = .NOT. FoundNodes

    DO i=1,NewMesh % NumberOfNodes
      IF(.NOT. FoundNodes(i)) THEN
        !Mark huge to indicate that this wasn't found
        PointLocalDistance(i) = HUGE(PointLocalDistance(i))
        CYCLE
      END IF

      IF(PointLocalDistance(i) == 0.0_dp) CYCLE
      PRINT *,ParEnv % MyPE,'Debug, point ',i,' found with local dist: ',PointLocalDistance(i)
    END DO

    !Sum up unfound nodes, and those where point wasn't exactly in element
    n = COUNT((.NOT. FoundNodes) .OR. (PointLocalDistance > 0.0_dp))
    dn = n
    IF(Debug) THEN
       PRINT *, 'Partition ',ParEnv % MyPE,' couldnt find ',n,'points!'
    END IF
    CALL SParActiveSUM(dn,2)

    !Special case: all found
    !----------------------
    IF ( dn==0 ) RETURN


    ! Exchange partition bounding boxes:                             
    ! ----------------------------------                              
    myBB(1) = MINVAL(OldMesh % Nodes % x)
    myBB(2) = MINVAL(OldMesh % Nodes % y)
    myBB(3) = MINVAL(OldMesh % Nodes % z)
    myBB(4) = MAXVAL(OldMesh % Nodes % x)
    myBB(5) = MAXVAL(OldMesh % Nodes % y)
    myBB(6) = MAXVAL(OldMesh % Nodes % z)

    myBB(HeightDimensions) = 0.0_dp
    myBB(HeightDimensions + 3) = 0.0_dp

    !Possibly need to adjust this - is it necessary to be extending the 
    !bounding box by a factor of 0.1
    epsBB = 0.1_dp * MAXVAL(myBB(4:6)-myBB(1:3))
    myBB(1:3) = myBB(1:3) - epsBB
    myBB(4:6) = myBB(4:6) + epsBB

    ALLOCATE(BB(6,ParEnv % PEs))
    DO i=1,ParEnv % PEs
       IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
       proc = i-1
       CALL MPI_BSEND( myBB, 6, MPI_DOUBLE_PRECISION, proc, &
            1099, MPI_COMM_WORLD, ierr )
    END DO
    DO i=1,COUNT(ParEnv % Active)-1
       CALL MPI_RECV( myBB, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
            1099, MPI_COMM_WORLD, status, ierr )
       proc = status(MPI_SOURCE)
       BB(:,proc+1) = myBB
    END DO

    Sending:IF ( n==0 ) THEN
       DEALLOCATE(FoundNodes, BB)
       DO i=1,ParEnv % PEs
          IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
          proc = i-1
          CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
               1101, MPI_COMM_WORLD, ierr )
       END DO
    ELSE
       ! Extract nodes that we didnt find from our own partition...   
       ! ------------------------------------------------------------ 
       ALLOCATE( Perm(n), nodes_x(n), nodes_y(n),nodes_z(n) ); Perm=0
       j = 0
       DO i=1,NewMesh % NumberOfNodes
          IF ( FoundNodes(i) .AND. (PointLocalDistance(i) <= 0.0_dp) ) CYCLE
          j = j + 1
          perm(j) = i
          nodes_x(j) = NewMesh % Nodes % x(i)
          nodes_y(j) = NewMesh % Nodes % y(i)
          nodes_z(j) = NewMesh % Nodes % z(i)
       END DO
       IF(ANY(HeightDimensions==1)) nodes_x = 0.0_dp
       IF(ANY(HeightDimensions==2)) nodes_y = 0.0_dp       
       IF(ANY(HeightDimensions==3)) nodes_z = 0.0_dp

       DEALLOCATE(FoundNodes)

       ! ...and ask those from others         
       ! -------------------------------               

       ALLOCATE(ProcSend(ParEnv % PEs))
       DO i=1,ParEnv % PEs
          IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
          proc = i-1
          ! extract those of the missing nodes that are within the other
          ! partions bounding box:                                    
          ! --------------------------------------------------------  
          myBB = BB(:,i) !Actually theirBB, but saves var names...
          npart = 0
          DO j=1,n
             IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                  nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                  nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
             npart = npart+1
          END DO
          ProcSend(proc+1) % n = npart
          IF ( npart>0 ) THEN
             ALLOCATE( xpart(npart),ypart(npart),zpart(npart),ProcSend(proc+1) % perm(npart) )
             npart = 0
             DO j=1,n
                IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                     nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                     nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
                npart=npart+1
                ProcSend(proc+1) % perm(npart)=j
                xpart(npart) = Nodes_x(j)
                ypart(npart) = Nodes_y(j)
                zpart(npart) = Nodes_z(j)
             END DO
          END IF
          ! send count...                                               
          ! -------------                                               
          CALL MPI_BSEND( npart, 1, MPI_INTEGER, proc, &
               1101, MPI_COMM_WORLD, ierr )
          IF ( npart==0 ) CYCLE
          ! ...and points                                               
          ! -------------                                               
          CALL MPI_BSEND( xpart, npart, MPI_DOUBLE_PRECISION, proc, &
               1102, MPI_COMM_WORLD, ierr )
          CALL MPI_BSEND( ypart, npart, MPI_DOUBLE_PRECISION, proc, &
               1103, MPI_COMM_WORLD, ierr )
          CALL MPI_BSEND( zpart, npart, MPI_DOUBLE_PRECISION, proc, &
               1104, MPI_COMM_WORLD, ierr )

          DEALLOCATE(xpart,ypart,zpart)
       END DO
       DEALLOCATE(nodes_x,nodes_y,nodes_z,BB)
    END IF Sending

    ! receive points from others:                                     
    ! ----------------------------                                    
    ALLOCATE(ProcRecv(Parenv % Pes))
    DO i=1,COUNT(ParEnv % Active)-1
       CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
            1101, MPI_COMM_WORLD, status, ierr )

       proc = status(MPI_SOURCE)
       ProcRecv(proc+1) % n = n

       IF ( n<=0 ) CYCLE

       ALLOCATE(ProcRecv(proc+1) % Nodes_x(n), &
            ProcRecv(proc+1) % Nodes_y(n),ProcRecv(proc+1) % Nodes_z(n))

       CALL MPI_RECV( ProcRecv(proc+1) % nodes_x, n, MPI_DOUBLE_PRECISION, proc, &
            1102, MPI_COMM_WORLD, status, ierr )
       CALL MPI_RECV( ProcRecv(proc+1) % nodes_y, n, MPI_DOUBLE_PRECISION, proc, &
            1103, MPI_COMM_WORLD, status, ierr )
       CALL MPI_RECV( ProcRecv(proc+1) % nodes_z, n, MPI_DOUBLE_PRECISION, proc, &
            1104, MPI_COMM_WORLD, status, ierr )
    END DO

    DO i=1,ParEnv % PEs
       IF ( Parenv % mype == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE

       proc = i-1
       n = ProcRecv(i) % n

       IF ( n==0 ) THEN
          CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
               2101, MPI_COMM_WORLD, ierr )
          CYCLE
       END IF

       ! Construct temporary mesh structure for the received points:   
       ! -----------------------------------------------------------   
       Nmesh => AllocateMesh()
       Nmesh % Nodes % x => ProcRecv(i) % nodes_x
       Nmesh % Nodes % y => ProcRecv(i) % nodes_y
       Nmesh % Nodes % z => ProcRecv(i) % nodes_z
       Nmesh % NumberOfNodes = n

       ALLOCATE(nperm(n), nReal(n))
       DO j=1,n
          nPerm(j)=j
       END DO
       nReal = NoData

       CALL VariableAdd( NMesh % Variables, NMesh, CurrentModel % Solver, &
            HeightName, 1, nReal, nPerm )

       NULLIFY(nReal, nPerm)
       nvars = 1 !not 0, because we always have 'height' var

       !-------------------------------------------------------------
       ! If we have a variable list to interpolate, add them to nMesh
       !-------------------------------------------------------------
       IF(PRESENT(Variables)) THEN
          OldVar => Variables
          DO WHILE(ASSOCIATED(OldVar))

             !Is the variable valid?
             IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
                  (OldVar % DOFs > 1) .OR. &                    !-multi-dof
                  (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
                  (OldVar % Name == HeightName) .OR. &          !-already got
                  OldVar % Secondary) THEN                      !-secondary

                OldVar => OldVar % Next
                CYCLE
             END IF

             nvars = nvars + 1
             
             ALLOCATE(nReal(n), nPerm(n))
             DO j=1,n
                nPerm(j)=j
             END DO
             nReal = NoData

             CALL VariableAdd( NMesh % Variables, NMesh, CurrentModel % Solver, & 
                  OldVar % Name, 1, nReal, nPerm)
             NULLIFY(nReal, nPerm)

             nVar => VariableGet(NMesh % Variables, OldVar % Name, .TRUE.)
             IF(.NOT. ASSOCIATED(nVar)) CALL Fatal("InterpolateVarToVarReduced",&
                  "Error trying to add variable to temporary mesh")


             IF ( ASSOCIATED(OldVar % PrevValues) ) THEN
                j = SIZE(OldVar % PrevValues,2)
                nvars = nvars+j
                ALLOCATE(nVar % PrevValues(n,j))
             END IF

             OldVar => OldVar % Next
          END DO
       END IF

       ! try interpolating values for the points:                      
       ! ----------------------------------------                      
       ALLOCATE( FoundNodes(n),&
            SendLocalDistance(n))

       FoundNodes=.FALSE.
       SendLocalDistance = 0.0_dp

       CALL InterpolateVarToVarReducedQ( OldMesh, nMesh, HeightName, HeightDimensions, &
            FoundNodes, SendLocalDistance, OldNodeMask, &
            OldElemMask=OldElemMask, Variables=Variables, GlobalEps=GlobalEps, &
            LocalEps=LocalEps, NumericalEps=NumericalEps )

       nfound = COUNT(FoundNodes)

       CALL MPI_BSEND( nfound, 1, MPI_INTEGER, proc, &
            2101, MPI_COMM_WORLD, ierr )

       unfound = n - nfound
       IF(unfound > 0) THEN
          PRINT *, ParEnv % MyPE, 'InterpVarToVar','Parallel: Found ',nfound,&
               ' nodes but still cant find ',unfound,' nodes!'
       END IF

       ! send interpolated values back to the owner:                   
       ! -------------------------------------------                   
       IF ( nfound>0 ) THEN
          ALLOCATE(vstore(nfound, nvars), vperm(nfound), WorkReal(nfound))
          vstore=0

          k = 0
          DO j=1,n
             IF ( .NOT.FoundNodes(j)) CYCLE
             k = k + 1
             vperm(k) = j
             WorkReal(k) = SendLocalDistance(j)
             Nvar => VariableGet( Nmesh % Variables,HeightName,ThisOnly=.TRUE.)
             nvars = 1
             vstore(k,nvars)=Nvar % Values(j)

             ! send all variables if requested
             !--------------------------------
             IF(PRESENT(Variables)) THEN
                OldVar => Variables
                DO WHILE(ASSOCIATED(OldVar))

                   IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
                        (OldVar % DOFs > 1) .OR. &                    !-multi-dof
                        (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
                        (OldVar % Name == HeightName) .OR. &          !-already got
                        OldVar % Secondary) THEN                      !-secondary

                      OldVar => OldVar % Next
                      CYCLE
                   END IF
                   nVar => VariableGet( Nmesh % Variables,OldVar % Name,ThisOnly=.TRUE.)
                   nvars = nvars+1
                   vstore(k,nvars)=nVar % Values(j)

                   IF ( ASSOCIATED(OldVar % PrevValues) ) THEN
                      DO l=1,SIZE(nVar % PrevValues,2)
                         nvars = nvars+1
                         vstore(k,nvars)=Nvar % PrevValues(j,l)
                      END DO
                   END IF
                   OldVar => OldVar % Next
                END DO
             END IF
             
             IF(Debug) PRINT *,'Partition ',ParEnv % MyPE,' found point with value: ',vstore(k,1)
          END DO

          !Pack up Local Distance from interp for this partition pair
          DEALLOCATE(SendLocalDistance)
          ALLOCATE(SendLocalDistance(nfound))
          SendLocalDistance = WorkReal
          DEALLOCATE(WorkReal)

          CALL MPI_BSEND( SendLocalDistance, nfound,MPI_DOUBLE_PRECISION, proc, &
               2100, MPI_COMM_WORLD,ierr )

          CALL MPI_BSEND( vperm, nfound, MPI_INTEGER, proc, &
               2102, MPI_COMM_WORLD, status, ierr )

          DO j=1,nvars
            CALL MPI_BSEND( vstore(:,j), nfound,MPI_DOUBLE_PRECISION, proc, &
                       2103+j, MPI_COMM_WORLD,ierr )
          END DO

          DEALLOCATE(vstore, vperm)
       END IF

       !These deallocations could be done with fewer lines,
       !but this way is more transparent (Nmesh % nodes point to procrecv)
       DEALLOCATE(ProcRecv(i) % Nodes_x, ProcRecv(i) % Nodes_y,&
            ProcRecv(i) % Nodes_z)
       NULLIFY(Nmesh % Nodes % x, Nmesh % Nodes % y, Nmesh % Nodes % z)
       CALL ReleaseMesh(Nmesh)
       DEALLOCATE(foundnodes, SendLocalDistance, nMesh)
    END DO
    DEALLOCATE(ProcRecv)

    ! Receive interpolated values:                                     
    ! ----------------------------                                     
    DO i=1,COUNT(ParEnv % Active)-1

       ! recv count:                                                   
       ! -----------                                                   
       CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
            2101, MPI_COMM_WORLD, status, ierr )

       proc = status(MPI_SOURCE)
       IF ( n<=0 ) THEN
          IF ( ALLOCATED(ProcSend) ) THEN
             IF ( ALLOCATED(ProcSend(proc+1) % Perm)) &
                  DEALLOCATE(ProcSend(proc+1) % Perm)
          END IF
          CYCLE
       END IF

       ALLOCATE(astore(n),vperm(n), RecvLocalDistance(n))

       ! recv permutation (where in the original array the             
       ! points the partition found are):                              
       ! --------------------------------------------------            
       CALL MPI_RECV( vperm, n, MPI_INTEGER, proc, &
            2102, MPI_COMM_WORLD, status, ierr )

       CALL MPI_RECV( RecvLocalDistance, n, MPI_DOUBLE_PRECISION, proc, &
            2100, MPI_COMM_WORLD, status, ierr )

       ! recv values and store:                                        
       ! ----------------------                                        
       !FIX HERE, INC MPI SEND ID
       CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
            2103+1, MPI_COMM_WORLD, status, ierr )

       Nvar => VariableGet( NewMesh % Variables,HeightName,ThisOnly=.TRUE.)
       IF(.NOT. ASSOCIATED(Nvar)) CALL Fatal("InterpVarToVar","Couldnt get variable from &
            &temporary mesh!")

       !Ddetermine which nodes were found BETTER than previously
       !i.e. compare RecvLocalDistance to PointLocalDistance.
       !In case the node previously wasn't found, BetterFound=.TRUE.,
       !because PointLocalDistance = HUGE()
       ALLOCATE(BetterFound(n))
       BetterFound = .FALSE.
       DO j=1,n 
         k=perm(ProcSend(proc+1) % Perm(vperm(j)))
         IF(RecvLocalDistance(j) < PointLocalDistance(k)) THEN
           BetterFound(j) = .TRUE.

           IF(Debug .AND. (PointLocalDistance(k) /= HUGE(PointLocalDistance(k)))) &
                PRINT *,ParEnv % MyPE,' better find node: ',k,' local dist new, old:',&
                RecvLocalDistance(j), PointLocalDistance(k)

           PointLocalDistance(k) = RecvLocalDistance(j)
         END IF
       END DO

       DO j=1,n
          IF(.NOT. BetterFound(j)) CYCLE
          k=perm(ProcSend(proc+1) % Perm(vperm(j)))
          IF ( Nvar % perm(k)>0 ) THEN
             Nvar % Values(Nvar % Perm(k)) = astore(j)
             IF(PRESENT(UnfoundNodes)) UnfoundNodes(k) = .FALSE.
          END IF
       END DO

       IF(PRESENT(Variables)) THEN
          OldVar => Variables
          nvars = 1
          
          DO WHILE(ASSOCIATED(OldVar))
             IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
                  (OldVar % DOFs > 1) .OR. &                    !-multi-dof
                  (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
                  (OldVar % Name == HeightName) .OR. &          !-already got
                  OldVar % Secondary) THEN                      !-secondary
                
                OldVar => OldVar % Next
                CYCLE
             END IF
             
             nvars = nvars + 1
             CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
                  2103+nvars, MPI_COMM_WORLD, status, ierr )

             nVar => VariableGet( NewMesh % Variables,OldVar % Name,ThisOnly=.TRUE.)

             IF (ASSOCIATED(nVar) ) THEN
                DO j=1,n
                   IF(.NOT. BetterFound(j)) CYCLE
                   k=perm(ProcSend(proc+1) % Perm(vperm(j)))
                   IF(ABS(astore(j) - NoData) < EPSILON(astore(j))) CYCLE !point found but not this specific var
                   IF ( Nvar % perm(k)>0 ) &
                        Nvar % Values(Nvar % Perm(k)) = astore(j)
                END DO
             END IF

             IF ( ASSOCIATED(OldVar % PrevValues) ) THEN
                DO l=1,SIZE(OldVar % PrevValues,2)
                   nvars=nvars+1
                   CALL MPI_RECV( astore, n, MPI_DOUBLE_PRECISION, proc, &
                        2103+nvars, MPI_COMM_WORLD, status, ierr )
                   IF(ASSOCIATED(Nvar)) THEN
                      DO j=1,n
                         IF(.NOT. BetterFound(j)) CYCLE
                         k=perm(ProcSend(proc+1) % Perm(vperm(j)))
                         IF ( Nvar % perm(k)>0 ) &
                              Nvar % PrevValues(Nvar % Perm(k),l) = astore(j)
                      END DO
                   END IF
                END DO
             END IF

             OldVar => OldVar % Next
          END DO
       END IF

       DEALLOCATE(astore,vperm,RecvLocalDistance, BetterFound, ProcSend(proc+1) % perm)

    END DO

    DEALLOCATE(PointLocalDistance)
    IF ( ALLOCATED(Perm) ) DEALLOCATE(Perm,ProcSend)

    !------------------------------------------------------------------------------
  END SUBROUTINE InterpolateVarToVarReduced

  !------------------------------------------------------------------------------
  SUBROUTINE InterpolateVarToVarReducedQ( OldMesh, NewMesh,HeightName,HeightDimensions, &
       FoundNodes, LocalDistances, OldNodeMask, NewNodeMask, OldElemMask, &
       Variables, GlobalEps, LocalEps, NumericalEps)
    !This subroutine takes each boundary node on the specified boundary of the new mesh and finds its height (y coord in 2D) by performing (DIM - 1) interpolaton through boundary elements of the old mesh.

    !-------------------------------------------------------------------------------
    TYPE(Mesh_t), TARGET, INTENT(IN)  :: OldMesh
    TYPE(Mesh_t), TARGET, INTENT(INOUT)  :: NewMesh
    TYPE(Variable_t), POINTER, OPTIONAL :: Variables
    CHARACTER(LEN=*) :: HeightName
    INTEGER, POINTER :: HeightDimensions(:)
    LOGICAL, OPTIONAL :: FoundNodes(:), OldNodeMask(:), &
         NewNodeMask(:), OldElemMask(:)
    REAL(KIND=dp), OPTIONAL :: GlobalEps, LocalEps, NumericalEps
    REAL(KIND=dp), ALLOCATABLE, OPTIONAL :: LocalDistances(:)
    !------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, VarOld, OldVar, NewVar
    TYPE(Element_t),POINTER :: Element       
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp), POINTER :: OldHeight(:), NewHeight(:)
    INTEGER, POINTER :: OldPerm(:), NewPerm(:)
    INTEGER :: i, j, k, l, n, ierr
    REAL(KIND=dp), DIMENSION(3) :: Point
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), DIMENSION(3) :: LocalCoordinates
    REAL(KIND=dp), POINTER :: ElementValues(:)
    REAL(KIND=dp) :: detJ, u,v,w,s, LocalDist
    LOGICAL :: Found, Debug
    REAL(KIND=dp) :: eps_global_limit, eps_local_limit,&
         eps_global_init, eps_local_init, eps_global, eps_local, eps_numeric
    !------------------------------------------------------------------------------

    !========================================
    ! Variables and parameter initializations
    !========================================

    Debug = .FALSE.

    IF(Debug) THEN
       PRINT *, 'Debug, present(OldNodeMask)', PRESENT(OldNodeMask)
       IF(PRESENT(OldNodeMask)) THEN
          PRINT *, 'Size OldNodeMask ', SIZE(OldNodeMask)
          PRINT *, 'count true: ', COUNT(OldNodeMask)
       END IF
    END IF

    ! Get the height variable
    !---------------------------------------    
    VarOld => VariableGet( OldMesh % Variables, HeightName, ThisOnly = .TRUE. )
    OldHeight => VarOld % Values
    OldPerm => VarOld % Perm

    ! If the target variable does not exist, create it
    !----------------------------------------------------------
    Var => VariableGet( NewMesh % Variables, HeightName, ThisOnly = .TRUE. )
    IF( .NOT. ASSOCIATED(Var) ) THEN
       ! This assumes that the mesh connectivity is the same...
       ALLOCATE( NewHeight(SIZE(OldHeight) ) )
       NewHeight = 0.0_dp
       ALLOCATE( NewPerm(SIZE(OldPerm) ) )
       NewPerm = OldPerm
       CALL VariableAdd( NewMesh % Variables, NewMesh, CurrentModel % Solver, &
            HeightName, 1, NewHeight, NewPerm )
       Var => VariableGet( NewMesh % Variables, HeightName, ThisOnly = .TRUE. )
       IF(.NOT. ASSOCIATED(Var)) CALL Fatal("InterpolateVarToVarReduced",&
            "Error adding the height variable to the new mesh.")
    END IF
    NewHeight => Var % Values
    NewPerm => Var % Perm


    ! Get epsilon values if specified
    !------------------------------------------------------------
    eps_global_init = 2.0e-10_dp
    eps_local_init = 1.0e-10_dp
    eps_numeric = 1.0e-10_dp

    IF(PRESENT(GlobalEps)) THEN
       eps_global_limit = GlobalEps
    ELSE
       eps_global_limit = ListGetConstReal( CurrentModel % Simulation,  &
            'Interpolation Global Epsilon', Found)
       IF(.NOT. Found) eps_global_limit = 1.0e-10
    END IF

    IF(PRESENT(LocalEps)) THEN
       eps_local_limit = LocalEps
    ELSE
       eps_local_limit = ListGetConstReal( CurrentModel % Simulation,  &
            'Interpolation Local Epsilon', Found )
       IF(.NOT. Found) eps_local_limit = 1.0e-10
    END IF

    IF(PRESENT(NumericalEps)) THEN
       eps_numeric = NumericalEps
    ELSE
       eps_numeric = ListGetConstReal( CurrentModel % Simulation,  &
            'Interpolation Numerical Epsilon', Found )
       IF(.NOT. Found) eps_numeric = EPSILON(eps_numeric)
    END IF

    !------------------------------------------------------------------------------
    n = OldMesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), &
         ElementNodes % z(n), ElementValues(n) )
    ElementNodes % x = 0.0_dp
    ElementNodes % y = 0.0_dp
    ElementNodes % z = 0.0_dp

    !========================================
    !             Action
    !========================================

    !------------------------------------------------------------------------------
    ! Loop over all nodes in the new mesh
    !------------------------------------------------------------------------------
    DO i=1,NewMesh % NumberOfNodes
       !------------------------------------------------------------------------------

       Found = .FALSE.

       IF( PRESENT(NewNodeMask)) THEN
          IF(NewNodeMask(i)) CYCLE
       END IF
       IF( NewPerm(i) == 0 ) CYCLE

       Point(1) = NewMesh % Nodes % x(i)
       Point(2) = NewMesh % Nodes % y(i)
       Point(3) = NewMesh % Nodes % z(i)
       Point(HeightDimensions) = 0.0_dp

       IF(Debug) PRINT *, 'Debug point no: ',i,' Perm: ',NewPerm(i), Point(1),Point(2),Point(3)
       !------------------------------------------------------------------------------
       ! Go through all old mesh boundary elements
       !------------------------------------------------------------------------------
       eps_global = eps_global_init
       eps_local = eps_local_init

       DO WHILE(.TRUE.)
          
          DO k=OldMesh % NumberOfBulkElements+1,&
               OldMesh % NumberOfBulkElements + OldMesh % NumberOfBoundaryElements

             Element => OldMesh % Elements(k)
             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes

             IF( ANY( OldPerm( NodeIndexes ) == 0 ) ) CYCLE

             IF(PRESENT(OldElemMask)) THEN
                IF(OldElemMask(k)) CYCLE
             END IF
             IF(PRESENT(OldNodeMask)) THEN
                IF( ANY( OldNodeMask( NodeIndexes ))) CYCLE
             END IF

             ElementNodes % x = 0.0_dp
             ElementNodes % y = 0.0_dp
             ElementNodes % z = 0.0_dp

             IF( ALL(HeightDimensions /= 1) ) &
                  ElementNodes % x(1:n) = OldMesh % Nodes % x(NodeIndexes)

             IF( ALL(HeightDimensions /= 2) ) &
                  ElementNodes % y(1:n) = OldMesh % Nodes % y(NodeIndexes)

             IF( ALL(HeightDimensions /= 3) ) &
                  ElementNodes % z(1:n) = OldMesh % Nodes % z(NodeIndexes)

             IF(Debug .AND. .FALSE.) THEN
                PRINT *, 'Debug element coords: '
                DO l=1,n
                   PRINT *, ElementNodes % x(l),ElementNodes % y(l),ElementNodes % z(l)
                END DO
             END IF

             Found = PointInElement( Element, ElementNodes, &
                  Point, LocalCoordinates, eps_global, eps_local, eps_numeric,&
                  LocalDistance=LocalDist)
             IF( Found ) EXIT
          END DO
          IF( Found ) EXIT

          eps_global = eps_global * 10.0_dp
          eps_local = eps_local * 10.0_dp
          IF(eps_global > eps_global_limit) EXIT
          IF(eps_local > eps_local_limit) EXIT

       END DO

       IF (.NOT.Found) THEN
          NULLIFY(Element)
          IF(PRESENT(FoundNodes)) FoundNodes(i) = .FALSE.
          CYCLE
       END IF

       IF(PRESENT(FoundNodes)) FoundNodes(i) = .TRUE.

       ElementValues(1:n) = OldHeight( OldPerm(NodeIndexes) )
       NewHeight(NewPerm(i)) = InterpolateInElement( &
            Element, ElementValues, LocalCoordinates(1), &
            LocalCoordinates(2), LocalCoordinates(3) )

       !Return element distances if requested.
       IF(PRESENT(LocalDistances)) LocalDistances(i) = LocalDist

       
       !-------------------------------------------------------
       ! Interpolate full variable list if requested
       !-------------------------------------------------------
       IF(PRESENT(Variables)) THEN
          OldVar => Variables
          DO WHILE(ASSOCIATED(OldVar))

             IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global   
                  (OldVar % DOFs > 1) .OR. &                    !-multi-dof
                  (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
                  (OldVar % Name == HeightName) .OR. &          !-already got
                  OldVar % Secondary) THEN                      !-secondary

                OldVar => OldVar % Next
                CYCLE
             END IF
             
             NewVar => VariableGet(NewMesh % Variables, OldVar % Name, .TRUE.)
             IF(.NOT. ASSOCIATED(NewVar)) THEN
                OldVar => OldVar % Next
                CYCLE
             END IF

             IF((NewVar % Perm(i) == 0) .OR. ANY(OldVar % Perm(NodeIndexes)==0)) THEN
                ! PRINT *, 'Debug interpvartovar, skipping ',OldVar % Name,' because of zero perm'
                OldVar => OldVar % Next
                CYCLE
             END IF
             
             ElementValues(1:n) = OldVar % Values( OldVar % Perm(NodeIndexes) )
             NewVar % Values(NewVar % Perm(i)) = InterpolateInElement( &
                  Element, ElementValues, LocalCoordinates(1), &
                  LocalCoordinates(2), LocalCoordinates(3) )
             
             IF ( ASSOCIATED( OldVar % PrevValues ) ) THEN
                DO j=1,SIZE(NewVar % PrevValues,2) !NewVar, not OldVar, to prevent error 
                   ElementValues(1:n) = &          !if sizes differ
                        OldVar % PrevValues(OldVar % Perm(NodeIndexes),j)
                   NewVar % PrevValues(NewVar % Perm(i),j) = &
                        InterpolateInElement( Element, ElementValues, &
                        LocalCoordinates(1), &
                        LocalCoordinates(2), LocalCoordinates(3) )
                END DO
             END IF

             OldVar => OldVar % Next
          END DO
       END IF

    END DO

    DEALLOCATE( ElementNodes % x, ElementNodes % y, &
         ElementNodes % z, ElementValues )


    !------------------------------------------------------------------------------
  END SUBROUTINE InterpolateVarToVarReducedQ

  !Subroutine designed to interpolate single missing points which sometimes
  !occur on the base and top interpolation, from surrounding nodes of the same mesh.
  SUBROUTINE InterpolateUnfoundPoint( NodeNumber, Mesh, HeightName, HeightDimensions,&
       ElemMask, NodeMask, Variables )

    TYPE(Mesh_t), TARGET, INTENT(INOUT)  :: Mesh
    TYPE(Variable_t), POINTER, OPTIONAL :: Variables
    CHARACTER(LEN=*) :: HeightName
    INTEGER :: NodeNumber
    INTEGER, POINTER :: HeightDimensions(:)
    LOGICAL, POINTER, OPTIONAL :: ElemMask(:),NodeMask(:)
    !------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: HeightVar, Var
    TYPE(Element_t),POINTER :: Element
    LOGICAL :: Parallel, Debug, HasNeighbours
    LOGICAL, ALLOCATABLE :: ValidNode(:)
    REAL(KIND=dp) :: Point(3), SuppPoint(3), weightsum, weight, Exponent, distance
    REAL(KIND=dp), ALLOCATABLE :: interpedValue(:), PartWeightSums(:), PartInterpedValues(:)
    INTEGER :: i,j,n,idx,NoNeighbours,NoSuppNodes,VarCount,&
         VarNo,proc,status(MPI_STATUS_SIZE), counter, ierr
    INTEGER, ALLOCATABLE :: NeighbourParts(:), WorkInt(:), SuppNodes(:)
    INTEGER, POINTER :: Neighbours(:)

    Debug = .FALSE.
    Parallel = ParEnv % PEs > 1

    HeightVar => VariableGet( Mesh % Variables, HeightName, &
         ThisOnly = .TRUE. )

    !The sought point
    Point(1) = Mesh % Nodes % x(NodeNumber)
    Point(2) = Mesh % Nodes % y(NodeNumber)
    Point(3) = Mesh % Nodes % z(NodeNumber)
    Point(HeightDimensions) = 0.0_dp

    !IDW exponent
    Exponent = 1.0

    !Is another partition also contributing to this
    NoNeighbours = SIZE(Mesh %  ParallelInfo % &
         NeighbourList(NodeNumber) % Neighbours) - 1
    HasNeighbours = NoNeighbours > 0

    !Create list of neighbour partitions (this will almost always be 0 :( )
    IF(HasNeighbours) THEN
      ALLOCATE(NeighbourParts(NoNeighbours))
      counter = 0
      DO i=1,NoNeighbours+1
        IF(Mesh %  ParallelInfo % NeighbourList(NodeNumber) % &
             Neighbours(i) == ParEnv % MyPE) CYCLE
        counter = counter + 1
        NeighbourParts(counter) = Mesh %  ParallelInfo &
             % NeighbourList(NodeNumber) % Neighbours(i)
      END DO
    END IF

    !Count this partition's relevant nodes
    ALLOCATE(ValidNode(Mesh % NumberOfNodes))
    ValidNode = .FALSE.

    !Start by marking .TRUE. based on ElemMask if present
    IF(PRESENT(ElemMask)) THEN
      DO i=1,SIZE(ElemMask)
        IF(ElemMask(i)) CYCLE
        n = Mesh % Elements(i) % TYPE % NumberOfNodes
        ValidNode(Mesh % Elements(i) % NodeIndexes(1:n)) = .TRUE.
      END DO
    ELSE
      ValidNode = .TRUE.
    END IF

    !Knock down by node mask if present
    IF(PRESENT(NodeMask)) THEN
      DO i=1,SIZE(NodeMask)
        IF(NodeMask(i)) ValidNode(i) = .FALSE.
      END DO
    END IF

    !Knock down nodes with 0 perm
    DO i=1,Mesh % NumberOfNodes
      IF(HeightVar % Perm(i) > 0) CYCLE
      ValidNode(i) = .FALSE.
    END DO

    IF(Debug) PRINT *,ParEnv % MyPE,'Debug, seeking nn: ',NodeNumber,' found ',&
         COUNT(ValidNode),' valid nodes.'

    ALLOCATE(WorkInt(100))
    WorkInt = 0

    !Cycle elements containing our node, adding other nodes to list
    NoSuppNodes = 0
    DO i=Mesh % NumberOfBulkElements+1,Mesh % NumberOfBulkElements &
         + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)
      n = Element % TYPE % NumberOfNodes

      !Doesn't contain our point
      IF(.NOT. ANY(Element % NodeIndexes(1:n)==NodeNumber)) CYCLE

      !Cycle element nodes
      DO j=1,n
        idx = Element % NodeIndexes(j)
        IF(idx == NodeNumber) CYCLE !sought node
        IF(ANY(WorkInt == idx)) CYCLE !already got
        IF(.NOT. ValidNode(idx)) CYCLE !invalid
        NoSuppNodes = NoSuppNodes + 1
        WorkInt(NoSuppNodes) = idx
      END DO
    END DO

    !If we aren't the only partition seeking this node, some supporting
    !nodes will also belong to these partitions. Easiest way to remove
    !duplicates is to set priority by partition number. So, if a
    !higher partition number (in NeighbourParts) also has a given supp
    !node, we delete it.
    IF(HasNeighbours) THEN
      DO i=1,NoSuppNodes
        Neighbours => Mesh % ParallelInfo % NeighbourList(WorkInt(i)) % Neighbours

        DO j=1,SIZE(Neighbours)
          IF(Neighbours(j) > ParEnv % MyPE .AND. ANY(NeighbourParts == Neighbours(j))) THEN
            WorkInt(i) = 0
            EXIT
          END IF
        END DO

      END DO

      NoSuppNodes = COUNT(WorkInt > 0)
      IF(Debug) PRINT *,ParEnv % MyPE, ' Debug, seeking ',NodeNumber,&
           ' higher partition has node, so deleting...'
    END IF

    ALLOCATE(SuppNodes(NoSuppNodes))
    SuppNodes = PACK(WorkInt, WorkInt > 0)
    DEALLOCATE(WorkInt)

    IF(Debug) PRINT *,ParEnv % MyPE,'Debug, seeking nn: ',NodeNumber,' found ',&
         NoSuppNodes,' supporting nodes.'

    !count variables if requested
    VarCount = 1 !1 = HeightVar
    IF(PRESENT(Variables)) THEN
      Var => Variables
      DO WHILE(ASSOCIATED(Var))

        !Is the variable valid?
        IF((SIZE(Var % Values) == Var % DOFs) .OR. &    !-global
             (Var % DOFs > 1) .OR. &                    !-multi-dof
             (Var % Name(1:10)=='coordinate') .OR. &    !-coord var
             (Var % Name == HeightName) .OR. &          !-already got
             Var % Secondary) THEN                      !-secondary

          Var => Var % Next
          CYCLE
        END IF
        IF(ANY(Var % Perm(SuppNodes) <= 0) .OR. &
             (Var % Perm(NodeNumber) <= 0)) THEN      !-not fully defined here
          Var => Var % Next
          CYCLE
        END IF

        VarCount = VarCount + 1
        Var => Var % Next
      END DO
    END IF

    ALLOCATE(interpedValue(VarCount))

    !Cycle supporting nodes, gathering weighted contributions
    !to HeightVar, and variables if requested
    weightsum = 0.0_dp
    interpedValue = 0.0_dp
    DO i=1,NoSuppNodes
      SuppPoint(1) = Mesh % Nodes % x(SuppNodes(i))
      SuppPoint(2) = Mesh % Nodes % y(SuppNodes(i))
      SuppPoint(3) = Mesh % Nodes % z(SuppNodes(i))
      SuppPoint(HeightDimensions) = 0.0_dp

      distance = 0.0_dp
      DO j=1,3
        distance = distance + (Point(j) - SuppPoint(j))**2.0_dp
      END DO
      distance = distance**0.5_dp

      weight = distance**(-exponent)
      weightsum = weightsum + weight

      interpedValue(1) = interpedValue(1) + &
           weight * HeightVar % Values(HeightVar % Perm(SuppNodes(i)))

      IF(PRESENT(Variables)) THEN
        VarNo = 1
        Var => Variables
        DO WHILE(ASSOCIATED(Var))

          !Is the variable valid?
          IF((SIZE(Var % Values) == Var % DOFs) .OR. & !-global
               (Var % DOFs > 1) .OR. &                    !-multi-dof
               (Var % Name(1:10)=='coordinate') .OR. &    !-coord var
               (Var % Name == HeightName) .OR. &          !-already got
               Var % Secondary) THEN                      !-secondary
            Var => Var % Next
            CYCLE
          END IF
          IF(ANY(Var % Perm(SuppNodes) <= 0) .OR. &
               (Var % Perm(NodeNumber) <= 0)) THEN      !-not fully defined here
            Var => Var % Next
            CYCLE
          END IF

          VarNo = VarNo + 1

          interpedValue(VarNo) = interpedValue(VarNo) + &
               weight * Var % Values(Var % Perm(SuppNodes(i)))

          Var => Var % Next
        END DO

      END IF
    END DO
    !PRINT *,ParEnv % MyPE,' Supporting nodes: ', NoSuppNodes, &
    !           ' x:', Mesh % Nodes % x(i),&
    !           ' y:', Mesh % Nodes % y(i),&
    !           ' z:', Mesh % Nodes % z(i)

    !PARALLEL STUFF
    IF(HasNeighbours) THEN
      ALLOCATE(PartWeightSums(NoNeighbours+1),&
           PartInterpedValues(VarCount * (NoNeighbours+1)))

      PartWeightSums(1) = WeightSum
      PartInterpedValues(1:VarCount) = interpedValue(1:VarCount)

      DO i=1,NoNeighbours
        proc = NeighbourParts(i)
        CALL MPI_BSEND( interpedValue, VarCount, MPI_DOUBLE_PRECISION, proc, &
             3000, MPI_COMM_WORLD,ierr )
        CALL MPI_BSEND( weightsum, 1, MPI_DOUBLE_PRECISION, proc, &
             3001, MPI_COMM_WORLD,ierr )

        CALL MPI_RECV( PartInterpedValues( (i*VarCount)+1  : (i+1)*VarCount), &
             VarCount, MPI_DOUBLE_PRECISION, proc, 3000, MPI_COMM_WORLD, status, ierr )
        CALL MPI_RECV( PartWeightSums(i+1), 1, MPI_DOUBLE_PRECISION, proc, &
             3001, MPI_COMM_WORLD, status, ierr )
      END DO

      interpedValue = 0.0_dp
      DO i=1,NoNeighbours+1
        DO j=1,VarCount
          interpedValue(j) = interpedValue(j) + PartInterpedValues(((i-1)*VarCount) + j)
        END DO
      END DO
      weightSum = SUM(PartWeightSums)
    END IF

    interpedValue = interpedValue/weightsum

    !Finally, put the interped values in their place
    HeightVar % Values(HeightVar % Perm(NodeNumber)) = interpedValue(1)

    IF(PRESENT(Variables)) THEN
      VarNo = 1
      Var => Variables
      DO WHILE(ASSOCIATED(Var))

        !Is the variable valid?
        IF((SIZE(Var % Values) == Var % DOFs) .OR. & !-global
             (Var % DOFs > 1) .OR. &                    !-multi-dof
             (Var % Name(1:10)=='coordinate') .OR. &    !-coord var
             (Var % Name == HeightName) .OR. &          !-already got
             Var % Secondary) THEN                      !-secondary
          Var => Var % Next
          CYCLE
        END IF
        IF(ANY(Var % Perm(SuppNodes) <= 0) .OR. &
             (Var % Perm(NodeNumber) <= 0)) THEN      !-not fully defined here
          Var => Var % Next
          CYCLE
        END IF

        VarNo = VarNo + 1

        Var % Values(Var % Perm(NodeNumber)) = interpedValue(VarNo)

        Var => Var % Next
      END DO
    END IF

    IF(HasNeighbours) DEALLOCATE(NeighbourParts)

  END SUBROUTINE InterpolateUnfoundPoint

!!!!!!!!!!!!!!!!!!!!!!!
! InterpolateUnfoundPointNearest

  !Subroutine designed to interpolate single missing points from the edges of the mesh 
  ! (using values from the old mesh), which occur due to glacier advance
  SUBROUTINE InterpolateUnfoundPointsNearest(UnfoundNodes, OldMesh, NewMesh, HeightName, & 
      HeightDimensions, Variables, OldNodeMask)

    TYPE(Mesh_t), TARGET, INTENT(IN)  :: OldMesh
    TYPE(Mesh_t), TARGET, INTENT(INOUT)  :: NewMesh
    TYPE(Variable_t), POINTER, OPTIONAL :: Variables
    LOGICAL, POINTER, OPTIONAL :: OldNodeMask(:)
    CHARACTER(LEN=*) :: HeightName
    INTEGER, POINTER :: HeightDimensions(:)
    LOGICAL, POINTER :: UnfoundNodes(:)
    !------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: OldHeightVar, NewHeightVar, OldVar, NewVar
    REAL(KIND=dp), ALLOCATABLE, TARGET :: BB(:,:), nodes_x(:), nodes_y(:), nodes_z(:), &
        xpart(:), ypart(:), zpart(:), mindists(:), minheights(:), partmindists(:),&
        partminheights(:), interped(:,:), partinterped(:,:)
    TYPE(Element_t),POINTER :: Element
    LOGICAL :: Parallel, Debug
    LOGICAL, ALLOCATABLE :: ValidNode(:), NoNearest(:), PartNoNearest(:)
    REAL(KIND=dp) :: Point(3), SuppPoint(3), dist, myBB(6), maxdist, partmindist
    INTEGER :: i, j, k, l, n, m, p, idx, proc, status(MPI_STATUS_SIZE), &
        ierr, npart, nn, VarNo, varn, mpicount
    INTEGER, ALLOCATABLE :: perm(:), mininds(:), partmininds(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName
    !------------------------------------------------------------------------------
    TYPE ProcRecv_t
       INTEGER :: n = 0
       REAL(KIND=dp), ALLOCATABLE :: nodes_x(:), nodes_y(:), nodes_z(:)
    END TYPE ProcRecv_t
    TYPE(ProcRecv_t),  ALLOCATABLE, TARGET :: ProcRecv(:)

    TYPE ProcSend_t
       INTEGER :: n = 0
       INTEGER, ALLOCATABLE :: perm(:)
    END TYPE ProcSend_t
    TYPE(ProcSend_t),  ALLOCATABLE :: ProcSend(:)
    !------------------------------------------------------------------------------

    SolverName='InterpolateUnfoundPointsNearest'
    
    CALL INFO(SolverName, 'Starting to interpolate unfound nodes', level = 3)
    
    Debug = .FALSE.
    Parallel = ParEnv % PEs > 1

    OldHeightVar => VariableGet( OldMesh % Variables, HeightName, &
         ThisOnly = .TRUE. )
    NewHeightVar => VariableGet( NewMesh % Variables, HeightName, &
         ThisOnly = .TRUE. )

    ! Initial minimum distance
    maxdist = 1.0e3

    ! Get this partition's bounding box:                             
    ! ----------------------------------                              
    myBB(1) = MINVAL(OldMesh % Nodes % x)
    myBB(2) = MINVAL(OldMesh % Nodes % y)
    myBB(3) = MINVAL(OldMesh % Nodes % z)
    myBB(4) = MAXVAL(OldMesh % Nodes % x)
    myBB(5) = MAXVAL(OldMesh % Nodes % y)
    myBB(6) = MAXVAL(OldMesh % Nodes % z)
    !TODO: make this related to maxdist, but for some reason not working with large maxdist
    myBB(1:3) = myBB(1:3) - 1.0e3
    myBB(4:6) = myBB(4:6) + 1.0e3

    myBB(HeightDimensions) = 0.0_dp
    myBB(HeightDimensions + 3) = 0.0_dp

    ! Get all partitions' bounding boxes:                             
    ! ---------------------------------- 
    ALLOCATE(BB(6,ParEnv % PEs))
    DO i=1,ParEnv % PEs
       IF ( Parenv % MyPE == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
       proc = i-1
       CALL MPI_BSEND( myBB, 6, MPI_DOUBLE_PRECISION, proc, &
            1399, MPI_COMM_WORLD, ierr )
    END DO
    DO i=1,COUNT(ParEnv % Active)-1
       CALL MPI_RECV( myBB, 6, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, &
            1399, MPI_COMM_WORLD, status, ierr )
       proc = status(MPI_SOURCE)
       BB(:,proc+1) = myBB
    END DO   

    ! Get valid nodes for this partition:                             
    ! ---------------------------------- 
    ALLOCATE(ValidNode(OldMesh % NumberOfNodes))
    ValidNode = .TRUE.
    
    !Knock down by node mask if present
    IF(PRESENT(OldNodeMask)) THEN
      DO i=1,SIZE(OldNodeMask)
        IF(OldNodeMask(i)) ValidNode(i) = .FALSE.
      END DO
    END IF

    !Knock down nodes with 0 perm
    DO i=1,OldMesh % NumberOfNodes
      IF(OldHeightVar % Perm(i) > 0) CYCLE
      ValidNode(i) = .FALSE.
    END DO

    ! Count variables                             
    ! ---------------------------------- 
    IF(PRESENT(Variables)) THEN   
      IF (Debug) PRINT *, ParEnv % MyPE, 'Variables present'
      OldVar => Variables
      VarNo = 1
      DO WHILE(ASSOCIATED(OldVar))
        !Is the variable valid?
        IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
              (OldVar % DOFs > 1) .OR. &                    !-multi-dof
              (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
              (OldVar % Name == HeightName) .OR. &          !-already got
              OldVar % Secondary) THEN                      !-secondary
          OldVar => OldVar % Next
          CYCLE
        END IF
        
        NewVar => VariableGet(NewMesh % Variables, OldVar % Name, .TRUE.)
        IF(.NOT. ASSOCIATED(NewVar)) THEN
          OldVar => OldVar % Next 
          CYCLE
        END IF
          
        VarNo = VarNo+1
        OldVar => OldVar % Next 
      END DO
      IF (Debug) PRINT *,ParEnv % MyPE, 'Counted ',VarNo,'variables'
    END IF
    
    ! Get closest points to unfound nodes for this partition:                             
    ! ---------------------------------- 
    n = COUNT(UnfoundNodes)
    IF (n > 0) THEN
      IF (Debug) PRINT *,ParEnv % MyPE, 'Interpolating ',n,'unfound nodes'
      ALLOCATE(mindists(n),minheights(n),mininds(n),NoNearest(n))
      mindists = maxdist
      minheights = 0
      NoNearest = .TRUE.
      m = 0
      DO nn=1,SIZE(UnfoundNodes)    
        IF (UnfoundNodes(nn)) THEN
          
          m = m+1        
          
          Point(1) = NewMesh % Nodes % x(nn)
          Point(2) = NewMesh % Nodes % y(nn)
          Point(3) = NewMesh % Nodes % z(nn)
          Point(HeightDimensions) = 0.0_dp
          
          IF (Debug) PRINT *,'UnfoundNode',nn,Point(1),Point(2)

          !Cycle through boundary elements of old mesh, find closest node to unfound node in new mesh
          DO i=OldMesh % NumberOfBulkElements+1,OldMesh % NumberOfBulkElements &
            + OldMesh % NumberOfBoundaryElements
            Element => OldMesh % Elements(i)
            n = Element % TYPE % NumberOfNodes
  
            !Cycle element nodes
            DO j=1,n
              idx = Element % NodeIndexes(j)
    
              IF(mininds(m) == idx) CYCLE !already got
              IF(.NOT. ValidNode(idx)) CYCLE !invalid
          
              ! Get coordinates for node
              SuppPoint(1) = OldMesh % Nodes % x(idx)
              SuppPoint(2) = OldMesh % Nodes % y(idx)
              SuppPoint(3) = OldMesh % Nodes % z(idx)
              SuppPoint(HeightDimensions) = 0.0_dp

              ! Get distance from unfound point
              dist = 0.0_dp
              DO k=1,3
                dist = dist + (Point(k) - SuppPoint(k))**2.0_dp
              END DO
              dist = dist**0.5_dp
              
              ! If it's the closest point so far, save it
              IF (dist < mindists(m)) THEN
                mindists(m) = dist
                minheights(m) = OldHeightVar % Values(OldHeightVar % Perm(idx))
                mininds(m) = idx
                !minind = idx
                NoNearest(m) = .FALSE. 
              END IF
            END DO
          END DO
        END IF
      END DO
      
      ! Mark unfound nodes as huge                            
      ! ---------------------------------- 
      IF (Debug) PRINT *,ParEnv % MyPE, 'Marking ',COUNT(NoNearest), 'nodes as still unfound'
      m = 0
      DO i=1,SIZE(UnfoundNodes)
        IF (UnfoundNodes(i)) THEN
          m = m+1
          IF (NoNearest(m)) THEN
            minheights(m) = HUGE(minheights(m))
            mininds(m) = 1
          END IF
        END IF
      END DO     

      ! Save interped values for variables, if minind exists                             
      ! ---------------------------------- 
      IF(PRESENT(Variables)) THEN  
        ALLOCATE(interped(VarNo,COUNT(UnfoundNodes)))
        interped = 0
        OldVar => Variables
        varn = 1
        DO WHILE(ASSOCIATED(OldVar))
          !Is the variable valid?
          IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
              (OldVar % DOFs > 1) .OR. &                    !-multi-dof
              (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
              (OldVar % Name == HeightName) .OR. &          !-already got
              OldVar % Secondary) THEN                      !-secondary
            OldVar => OldVar % Next
            CYCLE
          END IF
        
          NewVar => VariableGet(NewMesh % Variables, OldVar % Name, .TRUE.)
          IF(.NOT. ASSOCIATED(NewVar)) THEN
            OldVar => OldVar % Next 
            CYCLE
          END IF
          
          m = 0
          DO nn=1,SIZE(UnfoundNodes)    
            IF (UnfoundNodes(nn)) THEN
              m = m+1  
              IF (.NOT. NoNearest(m)) THEN
                interped(varn,m) = OldVar % Values(OldVar % Perm(mininds(m)))
              ELSE
                interped(varn,m) = HUGE(interped(varn,m))   
              END IF  
            END IF
          END DO           
          
          varn = varn+1
          OldVar => OldVar % Next 
        END DO
        IF (Debug) PRINT *,ParEnv % MyPE,'Successfully interped ',VarNo,'variables' 
        
      END IF !Variables present
    
    END IF !UnfoundNodes are present in partition

    ! Send out unfound nodes to other processes
    ! ---------------------------------- 
    n = COUNT(UnfoundNodes)
    IF (n == 0) THEN
      DO i=1,ParEnv % PEs
        IF ( (Parenv % MyPE == i-1) .OR. (.NOT. ParEnv % Active(i)) ) CYCLE
        proc = i-1
        CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
            1201, MPI_COMM_WORLD, ierr )
      END DO
    ELSE
      !IF (Debug) PRINT *,ParEnv % MyPE,'Getting unfound nodes ready to send'  
       
      ALLOCATE( perm(n), nodes_x(n), nodes_y(n),nodes_z(n) ); perm=0
      j = 0
      ! Extract nodes that we didn't find
      DO i=1,NewMesh % NumberOfNodes
        IF ( .NOT. UnFoundNodes(i) ) CYCLE
        j = j + 1
        perm(j) = i
        nodes_x(j) = NewMesh % Nodes % x(i)
        nodes_y(j) = NewMesh % Nodes % y(i)
        nodes_z(j) = NewMesh % Nodes % z(i)
      END DO
      IF(ANY(HeightDimensions==1)) nodes_x = 0.0_dp
      IF(ANY(HeightDimensions==2)) nodes_y = 0.0_dp       
      IF(ANY(HeightDimensions==3)) nodes_z = 0.0_dp
    
      ! And ask for nodes from others
      ALLOCATE(ProcSend(ParEnv % PEs))
      DO i=1,ParEnv % PEs
        IF ( Parenv % MyPE == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
        proc = i-1
        ! extract those of the missing nodes that are within the other
        ! partions' bounding boxes                                    
        ! --------------------------------------------------------  
        myBB = BB(:,i) !Actually theirBB, but saves var names...
        npart = 0
        DO j=1,n
          IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                  nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                  nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
          npart = npart+1
        END DO
        ProcSend(proc+1) % n = npart
        IF ( npart>0 ) THEN
          ALLOCATE( xpart(npart),ypart(npart),zpart(npart),ProcSend(proc+1) % perm(npart) )
          npart = 0
          DO j=1,n
            IF ( nodes_x(j)<myBB(1) .OR. nodes_x(j)>myBB(4) .OR. &
                  nodes_y(j)<myBB(2) .OR. nodes_y(j)>myBB(5) .OR. &
                  nodes_z(j)<myBB(3) .OR. nodes_z(j)>myBB(6) ) CYCLE
            npart=npart+1
            ProcSend(proc+1) % perm(npart)=j
            xpart(npart) = nodes_x(j)
            ypart(npart) = nodes_y(j)
            zpart(npart) = nodes_z(j)
          END DO
        END IF
        
        ! send count...                                               
        IF (Debug) PRINT *,ParEnv % MyPE,'Sending count of ',npart,' out of ',n,'to partition',&
            proc                                             
        CALL MPI_BSEND( npart, 1, MPI_INTEGER, proc, &
               1201, MPI_COMM_WORLD, ierr )
        IF ( npart==0 ) CYCLE
        ! ...and points                                                                                             
        CALL MPI_BSEND( xpart, npart, MPI_DOUBLE_PRECISION, proc, &
               1202, MPI_COMM_WORLD, ierr )
        CALL MPI_BSEND( ypart, npart, MPI_DOUBLE_PRECISION, proc, &
               1203, MPI_COMM_WORLD, ierr )
        CALL MPI_BSEND( zpart, npart, MPI_DOUBLE_PRECISION, proc, &
               1204, MPI_COMM_WORLD, ierr )
        DEALLOCATE(xpart,ypart,zpart)
      
      END DO
      
      DEALLOCATE(nodes_x,nodes_y,nodes_z,BB)
 
    END IF

    ! Receive points from others
    ! ----------------------------------
    ALLOCATE(ProcRecv(ParEnv % PEs))                                                                      
    DO i=1,COUNT(ParEnv % Active)-1
      CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
            1201, MPI_COMM_WORLD, status, ierr )

      proc = status(MPI_SOURCE)
      ProcRecv(proc+1) % n = n
      
      IF ( n<=0 ) CYCLE ! No points sent from that partition

      ALLOCATE(ProcRecv(proc+1) % nodes_x(n), ProcRecv(proc+1) % nodes_y(n), &
          ProcRecv(proc+1) % nodes_z(n))

      CALL MPI_RECV( ProcRecv(proc+1) % nodes_x, n, MPI_DOUBLE_PRECISION, proc, &
            1202, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV( ProcRecv(proc+1) % nodes_y, n, MPI_DOUBLE_PRECISION, proc, &
            1203, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV( ProcRecv(proc+1) % nodes_z, n, MPI_DOUBLE_PRECISION, proc, &
            1204, MPI_COMM_WORLD, status, ierr )            
    END DO


    ! Interpolate points sent by other partitions and send back
    ! ----------------------------------
    DO i=1,ParEnv % PEs
      IF ( Parenv % MyPE == i-1 .OR. .NOT. ParEnv % Active(i) ) CYCLE
      
      proc = i-1
      n = ProcRecv(i) % n 
      
      CALL MPI_BSEND( n, 1, MPI_INTEGER, proc, &
            1205, MPI_COMM_WORLD, status, ierr )
     
      IF ( n<=0 ) CYCLE ! No points to interpolate and send
     
      ALLOCATE(partmindists(n),partmininds(n),partminheights(n),PartNoNearest(n))
      partmindists = maxdist
      partmininds = 0
      partminheights = 0
      PartNoNearest = .TRUE.
      DO j=1,n
                  
        Point(1) = ProcRecv(proc+1) % nodes_x(j)
        Point(2) = ProcRecv(proc+1) % nodes_y(j)
        Point(3) = ProcRecv(proc+1) % nodes_z(j)
        Point(HeightDimensions) = 0.0_dp

        !Cycle through boundary elements of old mesh, find closest node to unfound node in new mesh
        DO l=OldMesh % NumberOfBulkElements+1,OldMesh % NumberOfBulkElements &
            + OldMesh % NumberOfBoundaryElements
          Element => OldMesh % Elements(l)
          nn = Element % TYPE % NumberOfNodes
  
          !Cycle element nodes
          DO p=1,nn
            idx = Element % NodeIndexes(p)
    
            IF(partmininds(j) == idx) CYCLE !already got
            IF(.NOT. ValidNode(idx)) CYCLE !invalid
          
            ! Get coordinates for node
            SuppPoint(1) = OldMesh % Nodes % x(idx)
            SuppPoint(2) = OldMesh % Nodes % y(idx)
            SuppPoint(3) = OldMesh % Nodes % z(idx)
            SuppPoint(HeightDimensions) = 0.0_dp

            ! Get distance from unfound point
            dist = 0.0_dp
            DO k=1,3
              dist = dist + (Point(k) - SuppPoint(k))**2.0_dp
            END DO
            dist = dist**0.5_dp
              
            ! If it's the closest point so far, save it
            IF (dist < partmindists(j)) THEN
              partmindists(j) = dist
              partminheights(j) = OldHeightVar % Values(OldHeightVar % Perm(idx))
              partmininds(j) = idx
              PartNoNearest(j) = .FALSE.
            END IF
          END DO
        END DO 
      END DO
     
      ! ... and send back interpolated values
      ! ----------------------------------   
      CALL MPI_BSEND( partmindists, n, MPI_DOUBLE_PRECISION, proc, &
               1206, MPI_COMM_WORLD, ierr )
      CALL MPI_BSEND( partminheights, n, MPI_DOUBLE_PRECISION, proc, &
               1207, MPI_COMM_WORLD, ierr )    

      DEALLOCATE( ProcRecv(i) % nodes_x, ProcRecv(i) % nodes_y, ProcRecv(i) % nodes_z  )
      
      IF(PRESENT(Variables)) THEN  
        ALLOCATE(partinterped(VarNo,n))
        partinterped = 0
        OldVar => Variables
        varn = 1
        DO WHILE(ASSOCIATED(OldVar))
          !Is the variable valid?
          IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
              (OldVar % DOFs > 1) .OR. &                    !-multi-dof
              (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
              (OldVar % Name == HeightName) .OR. &          !-already got
              OldVar % Secondary) THEN                      !-secondary
            OldVar => OldVar % Next
            CYCLE
          END IF
        
          NewVar => VariableGet(NewMesh % Variables, OldVar % Name, .TRUE.)
          IF(.NOT. ASSOCIATED(NewVar)) THEN
            OldVar => OldVar % Next 
            CYCLE
          END IF
          
          DO m=1,n    
            IF (.NOT. PartNoNearest(m)) THEN
              partinterped(varn,m) = OldVar % Values(OldVar % Perm(partmininds(m)))
            ELSE
              partinterped(varn,m) = HUGE(interped(varn,m))   
            END IF  
          END DO           
          
          varn = varn+1
          OldVar => OldVar % Next 
        END DO

        DO j=1,n
          CALL MPI_BSEND(partinterped(:,j),VarNo, MPI_DOUBLE_PRECISION, proc, &
               2200+j, MPI_COMM_WORLD, ierr )
        END DO

        DEALLOCATE(partinterped)        
      END IF !Variables present
     
      DEALLOCATE( partminheights, partmindists, PartNoNearest, partmininds)

    END DO

    ! Collect interpolated values from other partitions and compare to our value
    ! ----------------------------------  
    DO i=1,COUNT(ParEnv % Active)-1

      ! recv count:                                                   
      ! -----------                                                   
      CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
            1205, MPI_COMM_WORLD, status, ierr )
      
      proc = status(MPI_SOURCE)
      If (Debug) PRINT *,ParEnv % MyPE,'Receiving ',n,'from ',proc
      
      IF ( n<=0 ) THEN
         IF ( ALLOCATED(ProcSend) ) THEN
            IF ( ALLOCATED(ProcSend(proc+1) % perm)) &
                 DEALLOCATE(ProcSend(proc+1) % perm)
         END IF
         CYCLE
      END IF  
      
      ALLOCATE( partmindists(n), partminheights(n) )

      CALL MPI_RECV( partmindists, n, MPI_DOUBLE_PRECISION, proc, &
               1206, MPI_COMM_WORLD, status, ierr )
      CALL MPI_RECV( partminheights, n, MPI_DOUBLE_PRECISION, proc, &
               1207, MPI_COMM_WORLD, status, ierr )
            
      DO j=1,n
        k = ProcSend(proc+1) % perm(j)
        IF ( mindists(k) > partmindists(j)) THEN

          mindists(k) = partmindists(j)
          minheights(k) = partminheights(j)
          NoNearest = .FALSE.
            
          IF(PRESENT(Variables)) THEN
            CALL MPI_RECV( interped(:,k), VarNo, MPI_DOUBLE_PRECISION, proc, &
                  2200+j, MPI_COMM_WORLD, status, ierr )
          END IF
        END IF
      END DO
      
      DEALLOCATE( partmindists, partminheights, ProcSend(proc+1) % perm )
    END DO

  IF (ALLOCATED(NoNearest)) THEN
    IF (COUNT(NoNearest) == 0) THEN
      CALL Info(SolverName,'Found all remaining nodes that were not found by InterpVartoVar',&
      level = 3)
    ELSE
      CALL Fatal(SolverName, "Still unable to find nodes ")
    END IF
  END IF
  
  ! Finally, put interped values in their places
  ! ----------------------------------   
  n = COUNT(UnfoundNodes)
  IF (n > 0) THEN
    ! Put interpolated heights in their place
    IF (Debug) PRINT *,ParEnv % MyPE, 'Putting interpolated heights in their places'
    DO i=1,n
      j = perm(i)
      NewHeightVar % Values(NewHeightVar % Perm(j)) = minheights(i)
    END DO

    ! ...and put variable values in their place
    IF(PRESENT(Variables)) THEN
      IF (Debug) PRINT *,ParEnv % MyPE, 'Putting interpolated variables in their places'
      varn = 1
      OldVar => Variables
      DO WHILE(ASSOCIATED(OldVar))

        !Is the variable valid?
        IF((SIZE(OldVar % Values) == OldVar % DOFs) .OR. & !-global
               (OldVar % DOFs > 1) .OR. &                    !-multi-dof
               (OldVar % Name(1:10)=='coordinate') .OR. &    !-coord var
               (OldVar % Name == HeightName) .OR. &          !-already got
               OldVar % Secondary) THEN                      !-secondary
          OldVar => OldVar % Next
          CYCLE
        END IF
        
        NewVar => VariableGet(NewMesh % Variables, OldVar % Name, .TRUE.)
        IF(.NOT. ASSOCIATED(NewVar)) THEN
          OldVar => OldVar % Next
          CYCLE
        END IF

        DO i=1,n
          j = perm(i)
          NewVar % Values(NewVar % Perm(j)) = interped(varn,i)
        END DO
        
        varn = varn + 1        
        OldVar => OldVar % Next
      END DO
      
      DEALLOCATE(interped)
      
    END IF
    
    DEALLOCATE( mindists, minheights, mininds, NoNearest)
    
  END IF !n>0

  IF (Debug) PRINT *,ParEnv % MyPE,'DONE!!!! Running deallocations...'

  ! Deallocations
  ! ---------------------------------- 
  DEALLOCATE(ValidNode)
  DEALLOCATE(ProcRecv)
  IF (ALLOCATED(perm)) DEALLOCATE(perm, ProcSend)

  END SUBROUTINE InterpolateUnfoundPointsNearest


END MODULE InterpVarToVar
