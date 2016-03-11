Function BedMesh ( Model, nodenumber, dumy) RESULT(Zbed)
   USE types
   USE CoordinateSystems
   USE SolverUtils
   USE ElementDescription
   USE DefUtils
   IMPLICIT NONE
   TYPE(Model_t) :: Model
   TYPE(Solver_t), TARGET :: Solver
   INTEGER :: nodenumber
   REAL(KIND=dp) :: Mask
   INTEGER :: NMAX, i, j,Nb, DIM,ind(1)
   REAL(KIND=dp) :: x, y, z
   REAL(KIND=dp) :: Zbed, dumy, dist
   REAL(KIND=dp), ALLOCATABLE :: xbed(:), ybed(:), zbeds(:), indices(:)
   LOGICAL :: FirstTime=.True.
   character(len=60):: cmd

   SAVE FirstTime
   SAVE xbed, ybed, zbeds
   
   DIM = CoordinateSystemDimension()

   IF (FirstTime) THEN
   		Firsttime=.False.
        OPEN(10,file="Inputs/mesh_bed.dat")
		Read(10,*) Nb
        ALLOCATE(xbed(Nb), ybed(Nb), zbeds(Nb), indices(Nb))
        READ(10,*)(indices(i), xbed(i), ybed(i), zbeds(i), i=1,Nb)
        CLOSE(10)

   END IF
! Compute zbed for that point (x,y)
    IF (DIM==3) THEN
        x = Model % Nodes % x (nodenumber)
        y = Model % Nodes % y (nodenumber)
        z = Model % Nodes % z (nodenumber)
        ind=minloc(sqrt((x-xbed)**2+(y-ybed)**2),1)
    	dist=minval(sqrt((x-xbed)**2+(y-ybed)**2),1)
    	Zbed=zbeds(ind(1))
        
    ELSE IF (DIM==2) THEN
        x = Model % Nodes % x (nodenumber)
        Zbed = ybed(nodenumber)
        
    END IF
END FUNCTION