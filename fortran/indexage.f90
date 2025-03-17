MODULE mod_indexage
 IMPLICIT NONE
 CONTAINS
    !################################################################################################
    !#                 fonction pour calculer distance entre deux positions A et B                   #
    !################################################################################################
    DOUBLE PRECISION FUNCTION distance(A,B)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN) :: A, B
            distance = SQRT(SUM((B-A)**2))
    END FUNCTION distance

    !################################################################################################
    !#          SUBROUTINE etiquette pour avoir le numéro des atomes qui composent les cycles       #
    !#     chaînes alkyls sont donc exclues. cela permettra après de calculer facilement quantités  #
    !#                            tq distance de stacking ou shift des cycles                       #
    !################################################################################################
    SUBROUTINE etiquette(geoOut,lenGeoOut,symb,eps,monomer,Nb_stack,distCH)
!!
        INTEGER, DIMENSION(3), INTENT(IN)                       :: monomer
        ! comtient infos sur atomes et leurs quantités (ici il s'agit de NbS, NbC et NbH)
        INTEGER, INTENT(IN)                                     :: lenGeoOut, Nb_stack
        DOUBLE PRECISION, INTENT(IN)                            :: distCH, eps
        !longeur typique des liasons C-H et incertitudes sur distances liaisons C-H
        !DOUBLE PRECISION                                        :: distance
        ! fonction qui calcule distance, préalablement définie
        CHARACTER(1), DIMENSION(lenGeoOut)                       :: symb 
        ! contient les symboles associés à chaque coord
        DOUBLE PRECISION, DIMENSION(lenGeoOut,3), INTENT(INOUT)  :: geoOut
        !fichier contenant data geométrie en angs
        INTEGER		                                    	    :: i,j,l,n, count_ch
        INTEGER, DIMENSION(5*monomer(1),Nb_stack)              :: set
        DO i=1,Nb_stack
        ! première boucle traite sequentiellement les couches de polymer
            l = 0
            DO j=1+monomer(1)*(i-1),monomer(1)*i
            ! boucle sur les atomes de soufre de la i-ème couche de polymer
                l=l+1
                set(l,i) = j
            ENDDO
            DO j = 1+(Nb_stack*monomer(1))+monomer(2)*(i-1), (Nb_stack*monomer(1))+monomer(2)*i
            !boucle sur les atomes de carbone de la ième couche de polymer
                count_ch = 0
                DO n = Nb_stack*(monomer(1)+monomer(2))+1, lenGeoOut
                    if (symb(n) == 'H') then ! assure qu'il s'agit des coord d'un H
                        if (abs(distance(geoOut(n,:),geoOut(j,:)) - distCH)<= eps) then
                            count_ch = count_ch + 1 !comptage nombre de liaison ch pour carb DOnné
                        ENDif
                    ENDif 
                ENDDO
                if (count_ch <=1) then
                    !****** count_ch <-> nbr de liaison C-H que fait le j-ième carbone si plus d'un -> pas un C de cycle ****
                    ! en effet aucun carbone du cycle ne fait plus d'une liason avec des hydrogènes
                    l=l+1
                    set(l,i) = j
                ENDif 
            ENDDO
        ENDDO
        !*** Ecriture des numéros des atms du cycle dans file ***
        OPEN(1,file="etiquettage.out")
            DO i=1,5*monomer(1)
                write(1,*) set(i,:)
            ENDDO
        close(1)
    END SUBROUTINE etiquette
END MODULE mod_indexage
