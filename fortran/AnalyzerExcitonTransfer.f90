!##################################"
!# Main program
PROGRAM Analyzer
 USE mod_countNbOfLine
 USE mod_indexage
 USE mod_exciton
 USE mod_analyzeGeoTraj1
 USE mod_analyseGeoalgo2
 !USE mod_fileGenerator
 USE mod_stackMonomer
 USE mod_lenghtening
 USE mod_parseIonsAndWc
 USE mod_tab
 IMPLICIT NONE
    CHARACTER(30)                                   :: control, charac
    LOGICAL                                         :: bool = .TRUE., bool2
    CALL EXECUTE_COMMAND_LINE('date > OUTPUT' )
    OPEN(UNIT=9, FILE = 'OUTPUT',STATUS='old',POSITION='append')
    WRITE(9,*) '***********************************************************************&
            **********************************************'
    WRITE(9,*) '   Cheick Oumar DIARRA, Ph.D working on the MD modeling of the exciton&
            transfer in organic semiconductor  '
    WRITE(9,*) '  This program has been written exclusively for polymers. &
            For instance the P3MT and the P3HT. '
    WRITE(9,*) ' Before using on other molecules, please check the compatibility &
            of the used algorithms with your system (molecule)'
    WRITE(9,*) '                                     Version 1.1  &
       '
    WRITE(9,*) '***********************************************************************&
            **********************************************'
    WRITE(9,*) '                                                       &
            '
    CLOSE(9)
    !lecture des

    DO WHILE (bool)
        READ(*,*) control
        print*, control,' ', bool2
        IF (TRIM(control) == '&MOLECULE') THEN
            bool2 = .TRUE.
            DO WHILE (bool2)
                READ(*,*) control
                IF (TRIM(control) == 'MONOMER') READ(*,*) monomer
                IF (TRIM(control) == 'CYCLATOMS') READ(*,*) n_s, n_c, n_h
                IF (TRIM(control) == 'NBSTACKS') READ(*,*) Nb_stack
                IF (TRIM(control) == '}') bool2 = .FALSE.
            ENDDO
        ENDIF 
        IF (TRIM(control) == '&LENGTHENING') THEN
              bool2 = .TRUE.
            DO WHILE (bool2)
                READ(*,*) control
                IF (TRIM(control) == 'FILENAME') THEN
                    READ(*,*) charac
                    inputFile = TRIM(charac)
                ENDIF
                IF (TRIM(control) == 'SETTWOLASTCYCLES') THEN
                    ALLOCATE(lastCycleSet(n_s+n_c+n_h))
                    ALLOCATE(beforeLastCycleSet(n_s+n_c+n_h))
                    READ(*,*) lastCycleSet  !dernier cycle
                    READ(*,*) beforeLastCycleSet  !avant dernier cycle
                ENDIF
                IF (TRIM(control) == 'CHLENGTH') THEN
                    READ(*,*) distCH, eps
                ENDIF
                IF (TRIM(control) == 'NBCYCLE') READ(*,*) nbOfCycleOut
                IF (TRIM(control) == '}') bool2 = .FALSE.
            ENDDO
            CALL lenghtening(SUM(monomer),SUM(monomer)+2*nbOfCycleOut*(n_s + n_c + n_h),&
                    lastCycleSet,beforeLastCycleSet,TRIM(inputFile),&
                    LEN(TRIM(inputFile)),nbOfCycleOut,monomer,distCH,eps,n_s,n_c,n_h)
            DEALLOCATE(lastCycleSet); DEALLOCATE(beforeLastCycleSet)
        ENDIF
        !bool2 = .TRUE.
        IF (TRIM(control) == '&STACKING') THEN
            bool2 = .TRUE.
            DO WHILE (bool2)
                READ(*,*) control
                IF (TRIM(control) == 'FILENAME') THEN
                    READ(*,*) charac
                    inputFile = TRIM(charac)
                ENDIF
                IF (TRIM(control) == 'NBSTACKS') READ(*,*) Nb_stack
                IF (TRIM(control) == 'CHLENGTH') THEN
                        READ(*,*) distCH, eps
                ENDIF
                IF (TRIM(control) == 'SETCYCLATOMS') THEN
                    ALLOCATE(setLastCyle(n_s+n_c+n_h))
                    READ(*,*) setLastCyle  !dernier cycle
                ENDIF
                IF (TRIM(control) == 'STACKDIST') READ(*,*) dStack
                IF (TRIM(control) == 'ALTERNATE') THEN
                    alternate = .TRUE.
                ENDIF 
                IF (TRIM(control) == 'PERIODICITY') READ(*,*) stack
                IF (TRIM(control) == '}') bool2 = .FALSE.
            ENDDO
            CALL stackMonomer(SUM(monomer),SUM(monomer)*Nb_stack,inputFile,LEN(inputFile),Nb_stack,stack,&
                dStack,monomer,setLastCyle,n_s,n_c,n_h,distCH,eps,alternate)
            DEALLOCATE(setLastCyle)
        ENDIF
        IF (TRIM(control) == '&EXCITON') THEN
            bool2 = .TRUE.
            DO WHILE (bool2)
                READ(*,*) control
                IF (TRIM(control) == 'STEPTRAJ') THEN
                    READ(*,*) step_traj
                ENDIF
                IF (TRIM(control) == 'INITSTEP') READ(*,*) step_ini
                IF (TRIM(control) == 'TSTEP') THEN
                    READ(*,*) tstep !time step in a.u.
                    tstep = tstep*0.024 !in fs
                ENDIF
                IF (TRIM(control) == 'NUMBERBWC') READ(*,*) NbX
                IF (TRIM(control) == 'TOTSTEP') READ(*,*) Nb_step  !dernier cycle
                IF (TRIM(control) == 'BOXDIM') READ(*,*) dimBox
                IF (TRIM(control) == 'MSDALGO') THEN
                    READ(*,*) algoMSD
                    IF (algoMSD == 2) READ(*,*) TotStepDiff, IntervDiff
                ENDIF
                IF (TRIM(control) == 'UNFOLDALGO') READ(*,*) algoMajTraj
                IF (TRIM(control) == '}') bool2 = .FALSE.
            ENDDO
            INQUIRE(FILE='TRAJ_mol.xyz', EXIST = existe)
            IF (existe) THEN
                ALLOCATE(trajX(Nb_step*2,3)); ALLOCATE(symbX(Nb_step*2))
                lenLine = countNbOfLine('TRAJ_mol.xyz',len('TRAJ_mol.xyz'))
                CALL parse(trajX,symbX,2,'TRAJ_mol.xyz',len('TRAJ_mol.xyz'),&
                    lenLine,.TRUE.,Nb_step)
                CALL excitonDiffusion(trajX,2,Nb_step,TotStepDiff,IntervDiff,&
                    algoMajTraj,algoMSD,dimBox,step_traj,tstep)
                DEALLOCATE(trajX); DEALLOCATE(symbX)
            ELSE
                ALLOCATE(trajX(Nb_step*2,3)); ALLOCATE(symbX(Nb_step*2))
                ALLOCATE(trajMol(Nb_step*(monomer(1)+monomer(2)+monomer(3)&
                )*Nb_stack,3))
                ALLOCATE(symb(Nb_step*(monomer(1)+monomer(2)+monomer(3))*Nb_stack))
                lenLine = countNbOfLine('TRAJ_mol.xyz',len('TRAJ_mol.xyz'))
                CALL parseWC(trajMol,symb,SUM(monomer)*Nb_stack,&
                    trajX,symbX,2,Nb_step,Nb_stack,monomer,NbX,step_traj,step_ini)
                CALL excitonDiffusion(trajX,2,Nb_step,TotStepDiff,IntervDiff,&
                    algoMajTraj,algoMSD,dimBox,step_traj,tstep)
                DEALLOCATE(trajMol); DEALLOCATE(symb)
                DEALLOCATE(trajX); DEALLOCATE(symbX)
            ENDIF
        ENDIF
        IF (TRIM(control) == '&TRAJECTORY') THEN
            bool2 = .TRUE.
            DO WHILE (bool2)
                READ(*,*) control
                IF (TRIM(control) == 'STEPTRAJ') THEN
                    READ(*,*) step_traj
                ENDIF
                IF (TRIM(control) == 'INITSTEP') THEN
                    READ(*,*) step_ini
                ENDIF
                IF (TRIM(control) == 'TSTEP') THEN
                    READ(*,*) tstep !time step in a.u.
                    tstep = tstep*0.024 !in fs
                ENDIF
                IF (TRIM(control) == 'TOTSTEP') READ(*,*) Nb_step  !dernier cycle
                IF (TRIM(control) == 'BOXDIM') READ(*,*) dimBox
                IF (TRIM(control) == 'SIDECHAINS') READ(*,*) sidechains
                IF (TRIM(control) == 'NBSTACKS') READ(*,*) Nb_stack
                IF (TRIM(control) == 'NORMALVEC') READ(*,*) normBasePlane
                IF (TRIM(control) == '}') bool2 = .FALSE.
            ENDDO
            INQUIRE(FILE='TRAJEC.xyz', EXIST= existe)
            IF (existe) THEN 
                lenLine = countNbOfLine('TRAJEC.xyz',len('TRAJEC.xyz'))
                ALLOCATE(trajMol(Nb_step*(monomer(1)+monomer(2)+monomer(3))*Nb_stack,3))
                ALLOCATE(symb(Nb_step*(monomer(1)+monomer(2)+monomer(3))*Nb_stack))
                CALL parse(trajMol,symb,(monomer(1)+monomer(2)+monomer(3))*Nb_stack,&
                'TRAJEC.xyz',len('TRAJEC.xyz'),lenLine,.TRUE.,Nb_step)
                IF (sidechains) THEN
                    CALL planeToGeoProp(trajMol,lenTrajMol,normBasePlane,5*monomer(1),&
                    Nb_step,Nb_stack,monomer,step_traj,tstep)                  
                ELSE
                    CALL interatmdist(trajMol,SIZE(trajMol,1),Nb_stack,Nb_step,&
                    monomer,tstep,step_traj)
                ENDIF
                DEALLOCATE(trajMol); DEALLOCATE(symb)
            ELSE
                OPEN(UNIT = 1, FILE='OUTPUT', STATUS='old',POSITION='append')
                    WRITE(1,*) '!! Error: File TRAJEC.xyz not found !!'
                CLOSE(1)
            ENDIF
        ENDIF
        IF (TRIM(control) == '&END') bool = .FALSE. 
    ENDDO
    OPEN(UNIT=9, FILE = 'OUTPUT',STATUS='old',POSITION='append')
    WRITE(9,*) '_______________________________'
    WRITE(9,*) '         Thank you !'
END PROGRAM Analyzer