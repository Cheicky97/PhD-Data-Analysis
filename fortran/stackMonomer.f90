MODULE mod_stackMonomer
 USE mod_translate
 USE mod_initGeoForStack
 USE mod_fileGenerator, ONLY: geoFileGenerator
 USE mod_parse
 USE mod_sort
 USE mod_indexage
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE stackMonomer(lenGeo,lenGeoOut,inputFile,lenInputFile,Nb_stack,stack,d,monomer,set,n_s,n_c,n_h,distCH,eps,bool)
        INTEGER, DIMENSION(3), INTENT(IN)               :: monomer
        INTEGER, INTENT(IN)                             :: lenInputFile
        CHARACTER(lenInputFile), INTENT(IN)             :: inputFile
        INTEGER, INTENT(IN)                             :: Nb_stack, n_s, n_c, n_h
        DOUBLE PRECISION, INTENT(IN)				    :: d, distCH, eps !distance de translation
        INTEGER, INTENT(IN)							    :: stack !le numéro du stack à transformer
        INTEGER, DIMENSION(n_s+n_c+n_h),INTENT(IN)	    :: set ! l'indice du groupe d'atomes à translater
        LOGICAL, INTENT(IN)                             :: bool ! détermine si oui ou non on veut stack alternés
        INTEGER, INTENT(IN)                             :: lenGeo, lenGeoOut
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: geo !geometrie finale
        CHARACTER(1), DIMENSION(:), ALLOCATABLE		    :: symb !geometrie finale
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: geoOut !geometrie finale
        CHARACTER(1), DIMENSION(:), ALLOCATABLE		    :: symbOut !geometrie finale
        ! axe de translation 1=x, 2=y, 3=z
        LOGICAL                                         :: bool1 = .FALSE.
        DOUBLE PRECISION						        :: a, b !distance de translation
        INTEGER, DIMENSION(0)					        :: set0=0
        ALLOCATE(geo(lenGeo,3)); ALLOCATE(symb(lenGeo))
        CALL parse(Geo,symb,lenGeo,inputFile,lenInputFile,lenGeo+2,bool1,1)
        !!____initialisation des paramètres de translation________
        ALLOCATE(geoOut(lenGeoOut,3)); ALLOCATE(symbOut(lenGeoOut))
        CALL initializeGeo(geo,lenGeo,geoOut,lenGeoOut,symb,symbOut,Nb_stack,d,monomer) !initialize geoOut and stack with inter stack dist d 
        CALL sort(geo,lenGeo,symb,1,monomer,1) !sort file geo; important so to calculate a and b
        IF (bool .eqv. .TRUE.) THEN
            a = abs(geo(monomer(1),1)-geo(monomer(1)-1,1)) !dist betw S of last two cycle (x increasing)
            b = -abs(a+geo(monomer(1),1)-geo(1,1)) !dist betw S of first and last cycle (x increasing)
            CALL translation(1,a,stack,geoOut,lenGeoOut,set0,monomer,Nb_stack,0,0,0) !translate all the stack by a
            CALL translation(1,b,stack,geoOut,lenGeoOut,set,monomer,Nb_stack,n_s,n_c,n_h) !translate only set atoms by b
        ENDIF
        DEALLOCATE(geo); DEALLOCATE(symb)
        CALL sort(geoOut,lenGeoOut,symbOut,Nb_stack,monomer,1)
        CALL etiquette(geoOut,lenGeoOut,symbOut,eps,monomer,Nb_stack,distCH) !calculate and save numbers of each atoms of cycle per stack 
        CALL geoFileGenerator(geoOut,lenGeoOut,symbOut)  !give the ouput geo
        DEALLOCATE(geoOut); DEALLOCATE(symbOut)
        open(unit=1, file = 'OUTPUT', status = 'old', position = 'append')
        write(1,*) '______________________________________________________________________'
        write(1,*) '        Stack an input geometry along stacking direction (z)          '
        write(1,*) '______________________________________________________________________'
        write(1,*) ' geometry in input file has been successfully stacked.                   '
        write(1,*) ' data has been sorted along x (backbone direcrion) axis for each species.'
        write(1,*) ' File created : "etiquette.out" (contains numbers of atoms in cycles)    '
        write(1,*) '                                                                      '
        close(1)
    END SUBROUTINE stackMonomer
END MODULE mod_stackMonomer