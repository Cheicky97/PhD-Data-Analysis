MODULE mod_exciton
 USE mod_fileGenerator, ONLY: trajFileGenerator
 USE mod_excitonDIFf
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE excitonDiffusion(trajX,lenTraj,Nb_step,TotStepDiff,IntervDiff,algoMajTraj,&
    algoMSD,dimBox,step_traj,tstep)
        DOUBLE PRECISION, DIMENSION(3), INTENT(IN)                  :: dimBox
        INTEGER, INTENT(IN)                                         :: Nb_step,lenTraj
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTraj,3), INTENT(INOUT) :: trajX !traj exciton
        DOUBLE PRECISION, INTENT(IN)                                :: tstep
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTraj,3)              :: newTrajX
        INTEGER, INTENT(IN)                                         :: TotStepDiff,&
        IntervDiff,algoMajTraj, algoMSD, step_traj
        DOUBLE PRECISION                                            :: a
        open(unit=1, file = 'OUTPUT', status = 'old', position = "append")
        write(1,*) '______________________________________________________________________'
        write(1,*) '            Exciton diffusion : Calculation of the MSD                '
        write(1,*) '______________________________________________________________________'
        write(1,*) ' directions {x : backbone  y : sidechains  z : stacking}'
        IF(algoMajTraj == 1) THEN
            write(1,*) '>'
            write(1,*) 'Exciton trajectory unfolded with algorithm 1'
            write(1,*) '(This algo supposes CM of Molecule close to z = 0)'
            a = -dimBox(3)/2.0
            trajX(:,3) = trajX(:,3) + a
            newTrajX = trajX
            CALL majTrajeh1(trajX,lenTraj,newTrajX,1,Nb_step,dimBox)
            CALL majTrajeh1(trajX,lenTraj,newTrajX,2,Nb_step,dimBox)
            CALL majTrajeh1(trajX,lenTraj,newTrajX,3,Nb_step,dimBox)
        ELSE
            write(1,*) '>'
            write(1,*) 'Exciton trajectory unfolded with algorithm 2'
            write(1,*) '            (more general)                  '
            newTrajX = trajX
            CALL majTrajeh2(trajX,newTrajX,lenTraj,Nb_step,1,dimBox(1))
            CALL majTrajeh2(trajX,newTrajX,lenTraj,Nb_step,2,dimBox(2))
            CALL majTrajeh2(trajX,newTrajX,lenTraj,Nb_step,3,dimBox(3))
        ENDIF
        close(1)
        IF (algoMSD == 1) THEN
                CALL msd1(newTrajX,lenTraj,3,Nb_step,step_traj)
                CALL msd1(newTrajX,lenTraj,1,Nb_step,step_traj)
                CALL msd1(newTrajX,lenTraj,0,Nb_step,step_traj)
        ELSE
            CALL MSD2(newTrajX,lenTraj,IntervDIFf,TotStepDIFf,0,"3d",2,Nb_step,step_traj,tstep)
            CALL MSD2(newTrajX,lenTraj,IntervDIFf,TotStepDIFf,3,"z",1,Nb_step,step_traj,tstep)
        ENDIF
    END SUBROUTINE excitonDiffusion
END MODULE mod_exciton