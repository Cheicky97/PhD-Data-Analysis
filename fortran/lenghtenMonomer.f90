MODULE mod_lenghtening
 USE mod_parse
 USE mod_initialize
 USE mod_fileGenerator, ONLY: geoFileGenerator
 USE mod_sort
 USE mod_indexage
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE lenghtening(lenGeo,lenGeoOut,set1,set2,inputFile,lenInputFile,nbOfCycleOut,monomer,distCH,eps,n_s,n_c,n_h)
!!
        INTEGER, DIMENSION(3), INTENT(INOUT)                        :: monomer
        INTEGER, INTENT(IN)                                         :: lenInputFile
        INTEGER, INTENT(IN)                                         :: n_s, n_c, n_h  !# resp. the nbr of s, the nbr de C and the nbr de H per aromatic cycle (that include sidechains)
        INTEGER, INTENT(IN)                                         :: nbOfCycleOut   ! number of desired cycle        
        DOUBLE PRECISION, INTENT(IN)						        :: distCH, eps
        CHARACTER(lenInputFile), INTENT(IN)                         :: inputFile
        INTEGER, DIMENSION(n_s+n_c+n_h), INTENT(IN)				    :: set1, set2
        INTEGER, INTENT(IN)          		                        :: lenGeo != monomer(1)+monomer(2)+monomer(3)      !# tot number of atoms in the output file
        INTEGER, INTENT(IN)                                         :: lenGeoOut !lenGeo + INT((nbOfCycleOut-monomer(1)))*(n_s + n_c + n_h)
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE	            :: geoOut
        CHARACTER(1), DIMENSION(:), ALLOCATABLE		                :: symbOut
        CHARACTER(1), DIMENSION(:), ALLOCATABLE                     :: symb
        DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE		        :: geo
        DOUBLE PRECISION					                        :: t1
        LOGICAL                                                     :: bool = .FALSE.
        INTEGER, DIMENSION(3)                                       :: monomerOut
        monomerOut = (/monomer(1)+nbOfCycleOut*2*n_s,monomer(2)+nbOfCycleOut*2*n_c,monomer(3)+nbOfCycleOut*2*n_h/)
        ! importation de data geo input
        ALLOCATE(geo(lenGeo,3)); ALLOCATE(symb(lenGeo))
        CALL parse(geo,symb,lenGeo,inputFile,lenInputFile,lenGeo+2,bool,1)
        !# initialisation des paramètres de translation
        t1 = geo(monomer(1)+18,1)-geo(monomer(1)+9,1)  !!!!!! destiné à automatiser pour rendre le code plus robuste!!!!!!
        !# initialisation de geo_trans avec augmentation nb de cycles selon axe de chaine
        ALLOCATE(geoOut(lenGeoOut,3)); ALLOCATE(symbOut(lenGeoOut))
        CALL Lenghten(geo,lenGeo,geoOut,lenGeoOut,symbOut,set1,&
        size(set1),set2,monomer(1),monomer(2),monomer(3),n_s,n_c,n_h,nbOfCycleOut,t1)
        !trie
        DEALLOCATE(geo); DEALLOCATE(symb)
        CALL sort(geoOut,lenGeoOut,symbOut,1,monomerOut,1)
        !étiquettage
        CALL etiquette(geoOut,lenGeoOut,symbOut,eps,monomerOut,1,distCH)
        !enregistrement
        CALL geoFileGenerator(geoOut,lenGeoOut,symbOut)
        DEALLOCATE(geoOut); DEALLOCATE(symbOut)
        open(unit=1, file = 'OUTPUT', status = 'old', position = 'append')
        write(1,*) '________________________________________________________________________________'
        write(1,*) '       Increase number of aromatic cycle along the backbone direction           '
        write(1,*) '________________________________________________________________________________'
        write(1,*) ' data in file sorted along x coord for each species '
        write(1,*) ' File created : GEOMETRY_OUT.xyz'
        close(1)
    END SUBROUTINE lenghtening
END MODULE mod_lenghtening
