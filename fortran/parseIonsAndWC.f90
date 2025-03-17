MODULE mod_parseIonsAndWc
 USE mod_countNbOfLine
 USE mod_fileGenerator, ONLY: trajFileGenerator
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE parseWC(trajMol,symb,lenMol,trajX,symbX,lenTrajX,Nb_step,Nb_stack,monomer,Nb_X,step_traj,step_ini)
        !parse TRAJEC.xyz KIND OF FILE /!\ same structure as parseGeo, only difference being that second line of TRAJEC fil is of type integer
        ! len traj  here is the total number of atoms (per step)
!!
        INTEGER, DIMENSION(3), INTENT(IN)                               :: monomer
        INTEGER, INTENT(IN)                                             :: Nb_step, Nb_stack, Nb_X
        INTEGER, INTENT(IN)                                             :: lenMol, lenTrajX , step_traj,step_ini
        ! #lenLine : nb de lines du file   #lenMol: nb tot d'atm/step   #lenFile : nb de caract de filename
        CHARACTER(1), DIMENSION(Nb_step*lenMol), INTENT(INOUT)          :: symb
        CHARACTER(1), DIMENSION(Nb_step*lenTrajX), INTENT(INOUT)        :: symbX
        DOUBLE PRECISION, DIMENSION(Nb_step*lenMol,3), INTENT(INOUT)    :: trajMol 
        !traj de la molecule
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajX,3), INTENT(INOUT)  :: trajX 
        !traj de la molecule
        CHARACTER(1)                                                    :: car
        INTEGER                				                            :: a, lenLine
        INTEGER						                                    :: i,k,j,n
        ! #i : itérant pour parcourir file ligne par ligne, #j : entier pour répérer
        !début où fin d'un step de traj (ex: j=0 <-> 1er step)
        LOGICAL                                                         :: existe
        k=0;j=0;n=1
        lenLine = countNbOfLine('IONS+CENTERS.xyz',LEN('IONS+CENTERS.xyz'))
        INQUIRE( FILE="IONS+CENTERS.xyz", EXIST=existe )
        IF (existe) THEN
            OPEN(UNIT=1, FILE = 'IONS+CENTERS.xyz', ACTION='read')      
                DO i=1, lenLine
                    IF (i == 1 + (lenMol + Nb_X + 2)*(j+1)) j = j+1
                    IF (i <= 2 + (lenMol + Nb_X + 2)*j) THEN
                        READ(1,*) a
                    ENDIF
                    IF (i >= 3 + (lenMol + Nb_X + 2)*j .AND. i <= (lenMol + Nb_X + 2)*(j+1)) THEN
                        k = k+1
                        READ(1,*) symb(k), trajMol(k,:)
                        IF(i >= 3 + (lenMol + Nb_X + 2)*j .AND. i <= 2+(lenMol + Nb_X + 2)*j +&
                         (monomer(1) + monomer(2) + monomer(3))*Nb_stack) THEN
                            k = k+1
                            READ(1,*) symbX(k), trajMol(k,:)
                        ENDIF
                        !__________________extraction data trajectoire de tous les WC____________________________________
                        IF (i > 2 + (lenMol + Nb_X + 2)*j + &
                        (monomer(1) + monomer(2) + monomer(3))*Nb_stack .AND. &
                        i <= 2 + (lenMol + Nb_X + 2)*j + &
                        ((monomer(1) + monomer(2) + monomer(3))*Nb_stack + Nb_X)) THEN
                            READ(1,*) car, a, a, a
                            !__________________extraction data trajectoire exciton___________________________________________
                            IF(i > (lenMol + Nb_X + 2)*j+&
                            ((monomer(1) + monomer(2) + monomer(3))*Nb_stack + Nb_X) .AND.&
                             i <=2+(lenMol + Nb_X + 2)*j+&
                             ((monomer(1) + monomer(2) + monomer(3))*Nb_stack + Nb_X)) THEN
                                READ(1,*) symbX(n), trajX(n,:)
                                n = n+1
                            ENDIF 
                        ENDIF
                    ENDIF
                ENDDO
        CLOSE(1)
        CALL trajFileGenerator(trajMol,SIZE(trajMol,1),symb,lenMol,'mol',3,Nb_step,step_traj,step_ini)
        CALL trajFileGenerator(trajX,SIZE(trajX,1),symb,lenTrajX,'eh',2,Nb_step,step_traj,step_ini)
        open(unit=3, file='OUTPUT',action='write',position='append')
            write(3,*) '>'
            write(3,*) 'files created : TRAJ_mol.xyz, TRAJ_eh.xyz'
        close(3)
        ELSE
            OPEN(unit=2, file='error', action ='write', position='append')
                WRITE(2,*) ' [Error] : file IONS+CENTERS.xyz not found !'
            CLOSE(2)
        ENDIF
    END SUBROUTINE parseWC
END MODULE mod_parseIonsAndWc