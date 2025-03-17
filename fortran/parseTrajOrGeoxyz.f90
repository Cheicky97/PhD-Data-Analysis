!###########################################################################################################
!# module to import geo or traj data in tab
!# if bool == .true. -> input file <-> traj else input file of type geo xyz
MODULE mod_parse
 IMPLICIT NONE
    CONTAINS
    SUBROUTINE parse(trajMol,symb,lenMol,filename,lenFile,lenLine,bool,Nb_step)
        !parse TRAJEC.xyz KIND OF FILE /!\ same structure as parseGeo, only difference
        !being that second line of TRAJEC fil is of type integer
        ! len traj  here is the total number of atoms (per step)
        INTEGER, INTENT(IN)                                                     :: Nb_step
        LOGICAL, INTENT(IN)                                                     :: bool
        INTEGER, INTENT(IN)                                                     :: lenFile,lenMol, lenLine
        ! #lenLine : nb de lines du file   #lenMol: nb tot d'atm/step   #lenFile : nb de caract de filename
        CHARACTER(lenFile), INTENT(IN)                                          :: filename
        ! #filename : name file traj
        CHARACTER(1), DIMENSION(Nb_step*lenMol), INTENT(INOUT)                  :: symb
        DOUBLE PRECISION, DIMENSION(Nb_step*lenMol,3), INTENT(INOUT)            :: trajMol
        !traj de la molecule
        CHARACTER(1)                                                            :: car
        INTEGER                				                                    :: a
        INTEGER						                                            :: i,k,j 
        ! #i : itérant pour parcourir file ligne par ligne, !
        !#j : entier pour répérer début où fin d'un step de traj (ex: j=0 <-> 1er step)
        k=0;j=0
        OPEN(1, FILE = filename)
            DO i=1, lenLine
                IF (i == 1 + (lenMol + 2)*(j + 1)) j = j + 1
                IF (i == 1 + (lenMol + 2)*j) READ(1,*) a
                IF (i == 2 + (lenMol + 2)*j) THEN
                    IF (bool .eqv. .TRUE.) THEN
                        READ(1,*) a
                    ELSE
                        READ(1,*) car   !read only first charater of the second line
                    ENDIF
                ENDIF
                IF(i>=3+(lenMol+2)*j .and. i <= (lenMol+2)*(j+1)) THEN
                    k=k+1
                    READ(1,*) symb(k), trajMol(k,:)
                ENDIF
            ENDDO
        CLOSE(1)
    END SUBROUTINE parse
END MODULE mod_parse