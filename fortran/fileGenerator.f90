!##########################################################################
!# module containing subroutine to generate either geo file or traj file  #
MODULE mod_fileGenerator
 !for geo file
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE geoFileGenerator(geo,lenGeo,symb)
        INTEGER, INTENT(IN)						                :: lenGeo
        ! nb_char : nombre de caractère dans le suffix: 3 pour "mol", par exemple
        CHARACTER(1), DIMENSION(lenGeo),INTENT(IN)		        :: symb
        DOUBLE PRECISION, DIMENSION(lenGeo,3), INTENT(IN)	    :: geo !trajectoire
        INTEGER							                        :: i,j,k,n
        j=0;k=0; n=0
        open(1,file= "GEOMETRY_OUT.xyz")
            DO i=1, 2+lenGeo
                IF(i<=2) THEN
                    IF(i==1) THEN
                        WRITE(1,*) lenGeo
                    ELSE
                        WRITE(1,*) 'GEOMETRY_OUT'
                        n=n+1
                    ENDIF
                ENDIF
                IF(i>=3 .and. i <= 2+lenGeo) THEN
                    k=k+1
                    WRITE(1,*) symb(k), geo(k,:)
                ENDIF
            ENDDO
        close(1)
    END SUBROUTINE geoFileGenerator

    SUBROUTINE trajFileGenerator(traj,nl,symb,Nb_at,suffix,nb_char,Nb_step,step_traj,step_ini)
        ! genère fichier trajectoire du nom TRAJ_suffix.xyz (par exemple suffixe = "mol")
        INTEGER, INTENT(IN)					            :: nl, Nb_at, Nb_step, step_traj, step_ini
        !Nb_at : nb d'atome (pour mol, Nb_S+Nb_H+Nb_C
        INTEGER, INTENT(IN)					            :: nb_char
        !nombre de caractère dans le suffix: 3 pour "mol", par exemple
        DOUBLE PRECISION, dimension(nl,3), INTENT(IN)   :: traj
        !trajectoire
        CHARACTER(nb_char)					            :: suffix
        CHARACTER(1), dimension(nl),INTENT(IN)		    :: symb
        INTEGER						                    :: i,j,k,n
        j=0;k=0; n=0
        open(1,file="TRAJ_"//suffix//".xyz")
            DO i=1, Nb_step*(2+Nb_at)
                IF(i == 1+(Nb_at+2)*(j+1)) THEN
                    j = j+1
                ENDIF
                IF(i<=2+(2+Nb_at)*j) THEN
                    IF(i==1+(2+Nb_at)*j) THEN
                        WRITE(1,*) Nb_at
                    ELSE
                        WRITE(1,*) n*step_traj+step_ini+1
                            n=n+1
                    ENDIF
                ENDIF
                IF(i>=3+(2+Nb_at)*j .and. i <= (2+Nb_at)*(j+1)) THEN
                    k=k+1
                    WRITE(1,*) symb(k) , traj(k,:)
                ENDIF
            ENDDO
        close(1)
    END SUBROUTINE trajFileGenerator

    SUBROUTINE ehFileGenerator(traj,nl,symb,Nb_at,Nb_step,step_traj,step_ini)
        ! genère fichier trajectoire du nom TRAJ_suffix.xyz (par exemple suffixe = "mol")
        INTEGER, INTENT(IN)					                :: nl, Nb_at, Nb_step,step_traj,step_ini
        !Nb_at : nb d'atome (pour mol, Nb_S+Nb_H+Nb_C
        DOUBLE PRECISION, dimension(nl,3), INTENT(IN)	    :: traj !trajectoire
        !CHARACTER(nb_char)					                :: suffix
        CHARACTER(1), dimension(nl),INTENT(IN)		        :: symb
        INTEGER						                        :: i,j,k,n
        j=0;k=0;n=0
        open(1,file="TRAJ_h_CM.xyz")
        open(2,file="TRAJ_e_CM.xyz")
            DO i=1, Nb_step*(2+Nb_at)
                IF(i == 1+(Nb_at+2)*(j+1)) THEN
                    j = j+1
                ENDIF
                IF(i<=2+(2+Nb_at)*j) THEN
                    IF(i==1+(2+Nb_at)*j) THEN
                        WRITE(1,*) 1
                        WRITE(2,*) 1
                    ELSE
                        n=n+1
                        WRITE(1,*) (n-1)*step_traj+step_ini+1
                        WRITE(2,*) (n-1)*step_traj+step_ini+1
                    ENDIF
                ENDIF
                IF(i == 3+(2+Nb_at)*j) THEN
                    k=k+1
                    WRITE(1,*) symb(k), traj(k,:)
                ENDIF
                IF(i==4+(2+Nb_at)*j) THEN
                    k=k+1
                    WRITE(2,*) symb(k), traj(k,:)
                ENDIF
            ENDDO
        close(1)
        close(2)
    END SUBROUTINE ehFileGenerator
END MODULE mod_fileGenerator