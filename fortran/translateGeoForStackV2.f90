MODULE mod_translate
 IMPLICIT NONE
 CONTAINS
    !###################################################################################
    !# utilisé pour stacker monomer de manière alternée
    SUBROUTINE translation(axe,a,stack,geoOut,lenGeoOut,set,monomer,Nb_stack,n_s,n_c,n_h)
        INTEGER, DIMENSION(3), INTENT(IN)                                   :: monomer
        INTEGER, INTENT(IN)                                                 :: lenGeoOut, &
        Nb_stack, n_s, n_c, n_h
        DOUBLE PRECISION, DIMENSION(lenGeoOut,3), INTENT(INOUT)		        :: geoOut
        !geometrie finale
        DOUBLE PRECISION, INTENT(IN)						                :: a
        !distance de translation
        INTEGER, INTENT(IN)							                        :: axe 
        ! axe de translation 1=x, 2=y, 3=z
        INTEGER, INTENT(IN)							                        :: stack
        !le numéro du stack à transformer 
        ! nombre d'atomes sur lequel opérer une opération si n_s=0 --> tout le stack
        INTEGER, DIMENSION(n_s+n_c+n_h),INTENT(IN)					        :: set 
        ! l'indice du groupe d'atomes à translater
        INTEGER, DIMENSION(INT(AINT(FLOAT(Nb_stack)/FLOAT(stack))),n_s+n_c+n_h)	:: set_stack
        INTEGER								                                :: j,k,i ,l
        l=0
        IF (n_s+n_c+n_h >= 1) THEN
            !operation sur atomes indiqués dans set
            k = stack
            DO while(k <= Nb_stack)
                l = l+1
                DO j=1, n_s+n_c+n_h
                    IF (j<=n_s) THEN
                        set_stack(l,j) = set(j) + (k-1)*monomer(1)
                        geoOut(set_stack(l,j),axe) = geoOut(set_stack(l,j),axe) + a
                    ENDIF
                    IF (j>n_s .and. j<=n_s+n_c) THEN
                        set_stack(l,j) = set(j) + (Nb_stack*monomer(1)) + monomer(2)*(k-1)
                        geoOut(set_stack(l,j),axe) = geoOut(set_stack(l,j),axe) + a    
                    ENDIF
                    IF (j>n_c+n_s .and. j<=n_s+n_c+n_h) THEN
                        set_stack(l,j) = set(j) +&
                         (Nb_stack*(monomer(1) + monomer(2))) + monomer(3)*(k-1)
                        geoOut(set_stack(l,j),axe) = geoOut(set_stack(l,j),axe) + a
                    ENDIF
                ENDDO
                k = k + stack
            ENDDO
        ELSE
            !opération sur l'ensemble des atomes
            k = stack
            DO while (k <= Nb_stack)
                DO j=1 + monomer(1)*(k-1),monomer(1)*k
                    geoOut(j,axe) = geoOut(j,axe) + a
                ENDDO
                DO j=1 + (Nb_stack*monomer(1)) + monomer(2)*(k-1),&
                 (Nb_stack*monomer(1)) + monomer(2)*k
                    geoOut(j,axe) = geoOut(j,axe) + a
                ENDDO
                DO j = 1 + Nb_stack*(monomer(1) + monomer(2)) + monomer(3)*(k-1),&
                 Nb_stack*(monomer(1) + monomer(2)) + monomer(3)*k
                    geoOut(j,axe) = geoOut(j,axe) + a
                ENDDO  
                k = k + stack
            ENDDO
        ENDIF
    end SUBROUTINE translation
END MODULE mod_translate