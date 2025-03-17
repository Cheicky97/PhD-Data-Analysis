!##########################################################################
!# module that contains the geo initializer for the lenghtenin
!# 
MODULE mod_initialize
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE Lenghten(geo,lenGeo,geOut,lenGeoOut,symb,set1,lenSet1,set2,Nb_S,Nb_C,Nb_H,n_s,n_c,n_h,Nb_cyc,t)
        ! Be careful Nb_S here equal the total number of sulfur in the file geo, not the number per stack. Same for the other atoms
        ! n_s,n_c, n_h : resp. the numbers of sulfurs, carbons and hydrogens in one cycle (sidechains included)
        !lenGeout : the total number of atoms in the output geofile
        INTEGER, INTENT(IN)                                         :: lenGeoOut, lenGeo, lenSet1,&
         Nb_S, Nb_H, Nb_C, n_s, n_h,n_c, Nb_cyc
        DOUBLE PRECISION, DIMENSION(lenGeoOut,3), INTENT(INOUT) 	:: geOut    !# geo output data
        DOUBLE PRECISION, DIMENSION(lenGeo,3), INTENT(IN)		    :: geo      !# geo input data 
        INTEGER, DIMENSION(lenSet1), INTENT(IN)				        :: set1, set2
        CHARACTER(1), DIMENSION (lenGeoOut), INTENT(INOUT) 		    :: symb
        DOUBLE PRECISION,INTENT(IN)                                 :: t
        INTEGER							                            :: i,j,k,n,iter,l,m
        k=0;i=0;n=1;iter=0;l=1
        DO j=1,Nb_S !# loop on the sulfur to initialise geo data for them
            k=k+1
            symb(j) = "S"
            geOut(j,:) = geo(k,:)
        ENDDO
        DO j=Nb_S+1, Nb_S+Nb_cyc*n_s*2 !#supplemental sulfur
            symb(j) = "S"
            i=i+1
            IF (i <= n ) THEN
                iter = iter+1
                geOut(j,:) = geo(set1(l),:)
                geOut(j,1) = geo(set1(l),1)+iter*t  !iter count the cycle
            ELSE
                geOut(j,:) = geo(set2(l),:)
                geOut(j,1) = geo(set2(l),1)+iter*t !lengthening 
                n=n+2
            ENDIF
        ENDDO
        k=Nb_S
        ! for the Carbons now
        DO j=Nb_S+Nb_cyc*n_s*2 + 1,Nb_S+Nb_cyc*n_s*2 + Nb_C
            k=k+1
            symb(j) = "C"
            geOut(j,:) = geo(k,:)
        ENDDO
        iter=1;n=1;i=0;l=n_s;m=n_s
        DO j=Nb_S+Nb_cyc*n_s*2 + Nb_C + 1, Nb_S+Nb_cyc*(n_c+n_s)*2 + Nb_C
            symb(j) = "C"
            i=i+1
            IF (i <= n_c*n) THEN
                l=l+1
                geOut(j,:) = geo(Nb_S+set1(l),:)
                geOut(j,1) = geo(Nb_S+set1(l),1)+iter*t
            ELSE
                m=m+1
                geOut(j,:) = geo(Nb_S+set2(m),:)
                geOut(j,1) = geo(Nb_S+set2(m),1)+iter*t
                IF (i== n_c*(n+1)) THEN
                    n=n+2
                    iter = iter+1
                    l=n_s
                    m=n_s
                ENDIF
            ENDIF
        ENDDO
        k= Nb_S + Nb_C
        DO j=Nb_S+Nb_cyc*(n_c+n_s)*2 + Nb_C + 1,Nb_S+Nb_cyc*(n_c+n_s)*2 + Nb_C + Nb_H
            k=k+1
            symb(j) = "H"
            geOut(j,:) = geo(k,:)
        ENDDO
        iter=1;i=0;l=n_s+n_c;m=n_s+n_c;n=1
        DO j=Nb_S+Nb_cyc*(n_c+n_s)*2 + Nb_C + Nb_H + 1, Nb_S+Nb_cyc*(n_c+n_s+n_h)*2 + Nb_C+Nb_H
            symb(j) = "H"
            i=i+1
            IF (i <= n_h*n) THEN
                l=l+1
                geOut(j,:) = geo(Nb_S+Nb_C+set1(l),:)
                geOut(j,1) = geo(Nb_S+Nb_C+set1(l),1)+iter*t
            ELSE
                m=m+1
                geOut(j,:) = geo(Nb_S+Nb_C+set2(m),:)
                geOut(j,1) = geo(Nb_S+Nb_C+set2(m),1)+iter*t
                IF (i== n_h*(n+1)) THEN
                    n=n+2
                    iter = iter+1
                    l=n_s+n_c
                    m=n_s+n_c
                ENDIF
            ENDIF
        ENDDO
    END SUBROUTINE lenghten
END MODULE mod_initialize
