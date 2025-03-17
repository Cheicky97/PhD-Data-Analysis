MODULE mod_initGeoForStack
IMPLICIT NONE
CONTAINS
    SUBROUTINE initializeGeo(geo,lenGeo,geoOut,lenGeoOut,symb,symbOut,Nb_stack,d,monomer)
        INTEGER, DIMENSION(3), INTENT(IN)                               :: monomer !contient NbS, NbC, NbH par stack
        INTEGER, INTENT(IN)                                             :: Nb_stack, lenGeo, lenGeoOut
        DOUBLE PRECISION, DIMENSION(lenGeoOut,3), INTENT(INOUT) 	    :: geoOut
        DOUBLE PRECISION, DIMENSION(lenGeo,3), INTENT(IN)		        :: geo
        CHARACTER(1), DIMENSION (lenGeo), INTENT(IN)			        :: symb
        CHARACTER(1), DIMENSION (lenGeoOut), INTENT(INOUT) 		        :: symbOut
        DOUBLE PRECISION, INTENT(IN)					                :: d
        INTEGER							                                :: i,j,k
        DO i=1,Nb_stack
            k=0
            DO j=1+monomer(1)*(i-1),monomer(1)*i
                k = k+1
                symbOut(j) = symb(k)
                geoOut(j,1:2) = geo(k,1:2)
                geoOut(j,3) = geo(k,3) + (i-1)*d
            ENDDO
            k=monomer(1)
            DO j=1+(Nb_stack*monomer(1))+monomer(2)*(i-1), (Nb_stack*monomer(1))+monomer(2)*i
                k = k+1
                symbOut(j) = symb(k)
                geoOut(j,1:2) = geo(k,1:2)
                geoOut(j,3) = geo(k,3) + (i-1)*d
            ENDDO
            k= monomer(1) + monomer(2)
            DO j=1+Nb_stack*(monomer(1)+monomer(2))+monomer(3)*(i-1), Nb_stack*(monomer(1)+monomer(2))+monomer(3)*i
                k = k+1
                symbOut(j) = symb(k)
                geoOut(j,1:2) = geo(k,1:2)
                geoOut(j,3) = geo(k,3) + (i-1)*d
            ENDDO  
        ENDDO
    END SUBROUTINE initializeGeo
END MODULE mod_initGeoForStack