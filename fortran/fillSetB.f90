MODULE mod_fillset
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE fillset(set,dimset, setB,dimsetB,setC,dimsetC,Nb_stack,monomer)
        INTEGER, DIMENSION(3), INTENT(IN)                   :: monomer
        INTEGER, INTENT(IN)                                 :: dimsetB, dimsetC, Nb_stack, dimset
        INTEGER, DIMENSION(dimsetC,Nb_stack), INTENT(INOUT) :: setC
        INTEGER, DIMENSION(dimsetB,Nb_stack), INTENT(INOUT) :: setB
        INTEGER, DIMENSION(dimset,Nb_stack), INTENT(IN)     :: set
        INTEGER                                             ::i,j,k,l,m,n
        DO n=1,Nb_stack
            k=0; l=0;j=1; m = INT(monomer(1)/2) - 1
            DO i = 1, dimset
                IF (i<=monomer(1)) THEN
                    IF(MODULO(i,2) /= 0) THEN
                        k=k+1
                        setB(k,n) = set(i,n)
                    ENDIF
                    IF (i > INT(monomer(1)/2)) THEN
                        l = l+1
                        setC(l,n) =  set(i,n)
                    ENDIF
                ELSE
                    IF (i >= monomer(1) + 4*j+2 .AND.  i <= monomer(1) + 4*j + 3) THEN
                        k = k+1
                        setB(k,n) = set(i,n) !# 4 parcequ'il ya 4 C pour un S dans un sycle 
                    ENDIF
                    j = j + 2
                    IF (i > monomer(1) + 4*m .AND. i<= monomer(1) + 4*(m+1)) THEN
                        l = l+1
                        setC(l,n) = set(i,n)
                    ENDIF
                    m = m+1
                ENDIF
            ENDDO
        ENDDO
    END SUBROUTINE fillset
END MODULE mod_fillset