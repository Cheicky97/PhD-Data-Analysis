MODULE mod_tab
 IMPLICIT NONE
    INTEGER                                         :: lenTrajMol, lenTrajX
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: TrajMol
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: trajMolNoSideChain
    DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: TrajX
    CHARACTER(1), dimension(:), ALLOCATABLE         :: symbX
    INTEGER, DIMENSION(:,:), ALLOCATABLE            :: set
    INTEGER, DIMENSION(:), ALLOCATABLE              :: setLastCyle
    INTEGER, DIMENSION(:), ALLOCATABLE              :: lastCycleSet
    INTEGER, DIMENSION(:), ALLOCATABLE              :: beforeLastCycleSet
    CHARACTER(1), dimension(:), ALLOCATABLE         :: symb
    DOUBLE PRECISION, DIMENSION(1,3)                :: normBasePlane
    CHARACTER(:), ALLOCATABLE                       :: inputFile
    !TYPE AtomsOfMonomer
    !	!NbS : nbr de soufre dans un monomer, NbC,NbH pour C et H reps. ###
    !    INTEGER                    	                :: NbS, NbC, NbH
    !END TYPE AtomsOfMonomer

    !TYPE box
    !    INTEGER                    	                :: x, y, z 	!dim x, y et z de la boite
    !END TYPE box
    INTEGER, DIMENSION(3)                           :: monomer = (/6,30,24/)
    DOUBLE PRECISION, DIMENSION(3)                  :: dimBox = (/1.,1.,1./)
    INTEGER                                         :: n_s = 1, n_c = 5 , n_h = 4, nbOfCycleOut = 6
    LOGICAL                                         :: alternate, existe, sidechains = .TRUE.
    DOUBLE PRECISION                                :: dStack = 0, eps = 0, distCH = 1.10, tstep
    INTEGER                                         :: TotStepDiff = 0 ,IntervDiff = 0
    INTEGER                                         :: algoMajTraj = 2, algoMSD = 2, lenLine, step_ini, stack
    INTEGER                                         :: Nb_step = 1, Nb_stack = 1, step_traj, NbX=181, lenFile
END MODULE mod_tab