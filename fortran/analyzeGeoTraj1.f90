MODULE mod_analyzeGeoTraj1
 USE mod_fillset
 USE mod_analyzeTrajPlaneAlgo
 USE mod_indexage
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE planeToGeoProp(trajMol,lenTrajMol,normBasePlane,dimset,Nb_step,Nb_stack,monomer,step_traj,tstep)
        INTEGER, DIMENSION(3), INTENT(IN)                               :: monomer
        INTEGER, INTENT(IN)                                             :: Nb_step,Nb_stack,&
        step_traj,lenTrajMol
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)                    :: normBasePlane
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajMol,3), INTENT(IN)   :: trajMol
        DOUBLE PRECISION, INTENT(IN)                                    :: tstep
        INTEGER, INTENT(IN)                                             :: dimset
        INTEGER, DIMENSION(:,:), ALLOCATABLE                            :: set
        integer, dimension(:), ALLOCATABLE                              :: setB
        integer, dimension(:), ALLOCATABLE                              :: setC  !# 5 puisq'il 1 cyle = 4 C et 1 S
        DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE                 :: posCMA,posCMB,posCMC
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE                   :: normPlane
        DOUBLE PRECISION, DIMENSION(Nb_step)                            :: distPi, tiltAngle, dstack 
        DOUBLE PRECISION                                                :: dPimo, dStamo, tilAnmo, deltaP,&
         deltaA, deltaS
        INTEGER                                                         :: i,j,k,l,m,n,stack
        LOGICAL                                                         :: existe

        ALLOCATE(set(dimset,Nb_stack))
        ALLOCATE(setB(INT(3*monomer(1)/2)));ALLOCATE(setC(INT(5*monomer(1)/2)))
        ALLOCATE(posCMA(Nb_step,Nb_stack,3)); ALLOCATE(posCMB(Nb_step,Nb_stack,3))
        ALLOCATE(posCMC(Nb_step,Nb_stack,3)); ALLOCATE(normPlane(Nb_step,3))
        open(unit = 1, file = 'OUTPUT', status = 'old', position = 'append')
        write(1,*) '________________________________________________________________________________'
        write(1,*) '   Calculation of some structural properties from the trajectory at finite T    '
        write(1,*) '            This calculation uses planes to calculate distances                 '
        write(1,*) '________________________________________________________________________________'
        close(1)
        INQUIRE(FILE='etiquettage.out')
        IF(existe) THEN
            OPEN(UNIT=9,file='etiquettage.out' )
                DO i = 1, dimset
                    READ(1,*) set(i,:)
                ENDDO
            CLOSE(9)
        ELSE
            open(unit = 1, file = 'OUTPUT', status = 'old', position = 'append')
            write(1,*) '[error] : file etiqettage.out not found'
            write(1,*) '>> You can generate it from GEOMETRY.xyz using the &STICKS'
            write(1,*) ' '
            close(1)
        ENDIF

        CALL fillset(set,dimset,setB,SIZE(setB,1),setC,SIZE(setC,1),Nb_stack,monomer)
        distPi = 0.0
        tiltAngle = 0.0
        normPlane = 0.0
        CALL meanPos(trajMol,lenTrajMol,set,dimset,posCMA,Nb_step,Nb_stack)
        CALL meanPos(trajMol,lenTrajMol,setB,SIZE(setB),posCMB,Nb_step,Nb_stack)
        CALL meanPos(trajMol,lenTrajMol,setC,SIZE(setC),posCMC,Nb_step,Nb_stack)
        DO stack = 1, Nb_stack-1
            CALL stackToStack(posCMA,posCMB,posCMC,stack,normPlane,posCMA(:,stack+1,:),distPi,dstack,&
                tiltAngle,normBasePlane,Nb_stack,Nb_step,step_traj,tstep)
        ENDDO
        CALL  shiftSameCycle(trajMol,lenTrajMol,set,dimset,Nb_step,Nb_stack,monomer,step_traj,tstep)
        DEALLOCATE(set)
        DEALLOCATE(setB); DEALLOCATE(setC)
        DEALLOCATE(posCMA); DEALLOCATE(posCMB)
        DEALLOCATE(posCMC); DEALLOCATE(normPlane)
    END SUBROUTINE planeToGeoProp
END MODULE mod_analyzeGeoTraj1