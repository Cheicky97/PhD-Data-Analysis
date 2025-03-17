!#####################################################################################
!# Module to sort the file geo (each species separetly) on a given axe
MODULE mod_sort
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE sort(geoOut,lenGeoOut,symbOut,Nb_stack,monomer,axe)
        INTEGER, DIMENSION(3)                                               :: monomer
        INTEGER, INTENT(IN)                                                 :: lenGeoOut, Nb_stack
        DOUBLE PRECISION, DIMENSION(lenGeoOut,3), INTENT(INOUT)		        :: geoOut !geometrie finale
        CHARACTER(1), DIMENSION(lenGeoOut)                                  :: symbOut
        INTEGER, INTENT(IN)							                        :: axe  ! axe de translation 1=x, 2=y, 3=z
        DOUBLE PRECISION, DIMENSION(1,3)                                    :: geoEx !exchange geo
        INTEGER                                                             :: i,j,k,l
        k=0
        DO i=1,Nb_stack
            DO j=1+monomer(1)*(i-1),monomer(1)*i
                geoEx(1,:) = geoOut(j,:)
                DO k = j, monomer(1)*i
                    IF (geoEx(1,axe) >= geoOut(k,axe)) THEN
                        geoOut(j,:) = geoOut(k,:)
                        geoOut(k,:) = geoEx(1,:)
                        geoEx(1,:) = geoOut(j,:)
                    ENDIF
                ENDDO
            ENDDO
            DO j=1+(Nb_stack*monomer(1))+monomer(2)*(i-1), (Nb_stack*monomer(1))+monomer(2)*i
                geoEx(1,:) = geoOut(j,:)
                DO k = j, (Nb_stack*monomer(1))+monomer(2)*i
                    IF (geoEx(1,axe) >= geoOut(k,axe)) THEN
                        geoOut(j,:) = geoOut(k,:)
                        geoOut(k,:) = geoEx(1,:)
                        geoEx(1,:) = geoOut(j,:)
                    ENDIF
                ENDDO
            ENDDO
            DO j=1+Nb_stack*(monomer(1)+monomer(2))+monomer(3)*(i-1), Nb_stack*(monomer(1)+monomer(2))+monomer(3)*i
                geoEx(1,:) = geoOut(j,:)
                DO k = j, Nb_stack*(monomer(1)+monomer(2))+monomer(3)*i
                    IF (geoEx(1,axe) >= geoOut(k,axe)) THEN
                        geoOut(j,:) = geoOut(k,:)
                        geoOut(k,:) = geoEx(1,:)
                        geoEx(1,:) = geoOut(j,:)
                    ENDIF
                ENDDO
            ENDDO  
        ENDDO
    END SUBROUTINE sort
END MODULE mod_sort