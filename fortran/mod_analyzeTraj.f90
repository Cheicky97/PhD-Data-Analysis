MODULE mod_analyzeTrajPlaneAlgo
 IMPLICIT NONE
  CONTAINS
    !################################################################################################
    !#                 fonction pour calculer distance entre deux points A et B                     #
    DOUBLE PRECISION FUNCTION distance(A,B)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN) :: A, B
        distance = SQRT(SUM((B-A)**2))
    END FUNCTION distance

    !#######################################################
    !# calcul vect AB, sera utile pour rédéfinir  des axes #
    SUBROUTINE vectAB(A,B,AB)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)    :: A,B
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(INOUT) :: AB
        AB = B-A  
    END SUBROUTINE vectAB

    !########################################################
    !#            crosse product AB AC                      #
    SUBROUTINE crossProduct(AB,AC,normPlane)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)     ::AB,AC
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(INOUT)  ::normPlane
            !*************** normPlane= ABxAC ******************
        normPlane(1,1) = AB(1,2)*AC(1,3)-AB(1,3)*AC(1,2)  !ccord x de normPlane
        normPlane(1,2) = AB(1,3)*AC(1,1)-AB(1,1)*AC(1,3)  !ccord y de normPlane
        normPlane(1,3) = AB(1,1)*AC(1,2)-AB(1,2)*AC(1,1)  !ccord z de normPlane
        !PRINT*, '(',normPlane(1,1),' ',normPlane(1,2), normPlane(1,3),')'
    END SUBROUTINE crossProduct

    !#######################################################
    !# Cette fonction calcul le produit scalaire A.B
    DOUBLE PRECISION FUNCTION dotProduct(A,B)
        DOUBLE PRECISION, DIMENSION(1,3),INTENT(IN)  :: A,B
        dotProduct = SUM(A*B)
    END FUNCTION dotProduct

    !########################################################################
    !# Ce SUBROUTINE calcul la position moyenne des atomes dont les numéros 
    !# sont indiqués dans set. dim set est le nombre de ligne dans set
    !# Ce calcul est fait à chaque step d'acquisition de la trajectoire 
    !# pour chacune des stacks
    !# et les positions moynnes sont stockées dans posCM
    SUBROUTINE meanPos(trajMol,lenTrajMol,set,dimset,posCM,Nb_step,Nb_stack)
        INTEGER, INTENT(IN)                                             :: dimset,lenTrajMol,Nb_step,Nb_stack
        DOUBLE PRECISION, DIMENSION(Nb_step,nb_stack,3), INTENT(INOUT)  :: posCM
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajMol,3), INTENT(IN)   :: trajMol
        INTEGER, DIMENSION(dimset,Nb_stack),INTENT(IN)                  :: set
        INTEGER                                                         :: i,j,k
        posCM = 0.0
        DO i=1,Nb_step      !boucle sur step
            DO j=1,Nb_stack   !boucle sur stacks
                DO k=1, dimset  !boucle sur atoms
                    posCM(i,j,:) = posCM(i,j,:)+ trajMol(set(k,j)+(i-1)*lenTrajMol,:)
                ENDDO
                posCM(i,j,:) = posCM(i,j,:)/(dimset)
            ENDDO
        ENDDO
    END SUBROUTINE meanPos

    !###############################################################################
    !# SUBROUTINE pour trouver l'équation du plan contenant trois points A,B et C
    !# retourne coord de la normale au plan, l'équation du plan (normal,d)
    !# d position selon z à en x,y=0
    SUBROUTINE equaPlane(A,B,C,normPlane,equat)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)      :: A,B,C      ! A est la position moyenne du stack
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(INOUT)   :: normPlane  !vecteur normal au plan de coord (a,b,c)
        DOUBLE PRECISION, DIMENSION(1,4), INTENT(INOUT)   :: equat      !équation du plan (a,b,c,d) où d est distance à l'origine
        DOUBLE PRECISION, DIMENSION(1,3)                  :: AB, AC
!        DOUBLE PRECISION                                  :: dotProduct
        INTEGER                                           :: i
        CALL vectAB(A,B,AB) ! calcul de AB vecteur
        CALL vectAB(A,C,AC) ! calcul de AC vecteur
        CALL crossProduct(AB,AC,normPlane)
        DO i=1,3
            equat(1,i) = normPlane(1,i)
        ENDDO
        equat(1,4) = -1*dotProduct(normPlane,A)
    END SUBROUTINE equaPlane

    !########################################################################################
    !Calcul distance entre deux plans parallèles tq second plan contient pointOutOfPlane
    DOUBLE PRECISION FUNCTION distInterPlane(equat,normPlane,pointOutOfPlane)
        DOUBLE PRECISION, DIMENSION(1,4), INTENT(IN)    :: equat
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)    :: normPlane, pointOutOfPlane
        !DOUBLE PRECISION                                :: dotProduct
        distInterPlane = (ABS(dotProduct(normPlane,pointOutOfPlane)+equat(1,4)))/(SQRT(dotProduct(normPlane,normPlane)))
    END FUNCTION distInterPlane

    !#################################################################################################################
    !# fonction calcul et retourn la distance de stacking qui est la distance vertical (selon z) entre deux point
    !#situé sur deux plan parallèle à la même x,y
    DOUBLE PRECISION FUNCTION distStack(equat,pointOutOfPlane)
        DOUBLE PRECISION, DIMENSION(1,4), INTENT(IN)      :: equat
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)      :: pointOutOfPlane
 !       DOUBLE PRECISION                                  :: dotProduct
        DOUBLE PRECISION                                  :: distB ! position en z quand x,y=0 pour le plan parall en pointOutOf..
        distB = -1*(dotProduct(equat(1,1:3),pointOutOfPlane))
        distStack = ABS(distB-equat(1,4))/ABS(equat(1,3))
    END FUNCTION distStack

    !#########################################################################################################################
    !# calcul angle entre deux vecteurs BaseVecteur est vecteur orthogo au plan par rapport auquel on souhaite calculer angle
    DOUBLE PRECISION FUNCTION cosTiltAngle(normPlane,baseVect)
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)    :: normPlane, baseVect
        !DOUBLE PRECISION                                :: dotProduct
        cosTiltAngle = dotProduct(normPlane,baseVect)/(SQRT(dotProduct(normPlane,normPlane))*SQRT(dotProduct(baseVect,baseVect)))
    END FUNCTION cosTiltAngle

    !##############################################################################################
    !# SUBROUTINE pour calculer dist d'écart à la configuration où les souffre de deux stack succ
    !# sont orientés dans le même sens
    !# !!!!!!!!!!!!!!! Algorithme à révoir; (à ne donc pas utiliser pour le moment) !!!!!!!!!!!!!!!!!!
    !# !!!!! algo valable pour un dimer(Nb_stack=2), tetramer p3ht améliorations necessaires !!!!!
    SUBROUTINE shiftSameCycle(trajMol,lenTrajMol,set,dimset,Nb_step,Nb_stack,monomer,step_traj,tstep)
        INTEGER, DIMENSION(3), INTENT(IN)                               :: monomer
        DOUBLE PRECISION, INTENT(IN)                                    :: tstep
        INTEGER, INTENT(IN)                                             :: dimset,Nb_step,Nb_stack,step_traj,lenTrajMol
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajMol,3), INTENT(IN)   :: trajMol
        INTEGER, DIMENSION(dimset,Nb_stack), INTENT(IN)                 :: set
        DOUBLE PRECISION, DIMENSION(Nb_stack,monomer(1),3)              :: z
        DOUBLE PRECISION, DIMENSION(1,monomer(1)-1)                     :: distx, disty
        DOUBLE PRECISION, DIMENSION(1,Nb_stack-1)                       :: distmoyx, distmoyy
        DOUBLE PRECISION                                                :: shiftmoyy =0.0, shiftmoyx=0.0
        INTEGER                                                         :: i,j,k,m,n
        INTEGER                                                         :: stack, cycl
        OPEN(unit=1, file = 'shiftCycleX.out')
        OPEN(unit=2, file = 'shiftCycleY.out')
            DO i=1, Nb_step   ! boucle sur les step de la traj 
                distmoyx = 0; distmoyy=0
                DO j=1,Nb_stack ! boucle sur les stacks au i-ième step  
                    DO k=1,monomer(1)   ! boucle sur les cycles du j-ème stack
                        z(j,k,:) = z(j,k,:) + trajMol(set(k,j)+(i-1)*lenTrajMol,:)
                        DO m=1,4
                            z(j,k,:) = z(j,k,:) + trajMol(set(m+monomer(1)+(k-1)*4,j)+(i-1)*lenTrajMol,:)
                        ENDDO
                        z(j,k,:) = z(j,k,:)/5.0 !moynne des postions des atoms de cycle k, du stack j
                    ENDDO
                ENDDO
                ! à modifier
                DO stack = 1, Nb_stack-1
                    distx = 0; disty = 0; n=0
                    DO cycl =1,monomer(1)-1
                        distx(1,cycl) = ABS(z(stack,cycl,1)-z(stack+1,cycl+1,1)) !calcul des shif x entre cycle de stack et stack+1
                        disty(1,cycl) = ABS(z(stack,cycl,2)-z(stack+1,cycl+1,2)) !calcul des shif y entre cycle de stack et stack+1
                    ENDDO
                    distmoyx(1,stack) = SUM(distx)/(monomer(1)-1) !calcul shift moy entre stack et stack+1 selon axe 1
                    distmoyy(1,stack) = SUM(disty)/(monomer(1)-1) !calcul shift moy entre stack et stack+1 selon axe 1
                ENDDO
                WRITE(1,*) (1+(i-1)*step_traj)*tstep, distmoyx(1,:)
                WRITE(2,*) (1+(i-1)*step_traj)*tstep, distmoyy(1,:)
                shiftmoyx = shiftmoyx + SUM(distmoyx)
                shiftmoyy = shiftmoyy + SUM(distmoyy)
            ENDDO
        CLOSE(1)
        CLOSE(2)
        open(unit=3, file = 'OUTPUT', status = 'old', position = "append")
        write(3,*) 'Shift of equivalent cycle along x and y have been calculated'
        write(3,*) '<Shiftx> (Angstrom) = ', shiftmoyx/Nb_step
        write(3,*) '<Shifty> (Angstrom) = ', shiftmoyy/Nb_step
        write(3,*) 'Files shiftCycleX.out and shiftCycleX.out have been created.'
        write(3,*) 'shift file of structure : time  shiftmoyx(stack1-stack2) ...'
        close(3)
    END SUBROUTINE shiftSameCycle

    !############################################################################################################
    !#  SUBROUTINE pour les différents calculs
    SUBROUTINE stackToStack(posCMA,posCMB,posCMC,stack,normPlane,&
    pointOutOfPlane,distPi,dstack,tiltAngle,normBasePlane,Nb_stack,Nb_step,step_traj,tstep)
        INTEGER, INTENT(IN)                                             :: stack,Nb_step,Nb_stack,step_traj
        DOUBLE PRECISION, DIMENSION(Nb_step,Nb_stack,3), INTENT(INOUT)  :: posCMA,posCMB,posCMC
        DOUBLE PRECISION, DIMENSION(Nb_step,3), INTENT(INOUT)           :: normPlane,pointOutOfPlane
        DOUBLE PRECISION, DIMENSION(Nb_step), INTENT(INOUT)             :: distPi, tiltAngle, dstack 
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)                    :: normBasePlane
        DOUBLE PRECISION, INTENT(IN)                                    :: tstep
!        DOUBLE PRECISION                                                :: distInterPlane, cosTiltAngle, distStack
        DOUBLE PRECISION, DIMENSION(1,4)                                :: equat
        DOUBLE PRECISION                                                :: dPimo, dStamo, tilAnmo, deltaP, deltaA, deltaS
        INTEGER                                                         :: i,j,k
        dPimo=0; dStamo=0; tilAnmo=0; deltaP = 0; deltaS = 0; deltaA=0
        OPEN(unit=1,file='distPi.out')
        OPEN(unit=2,file='distStack.out')
        OPEN(unit=3,file='TiltAngle.out')
            DO i=1,Nb_step
                CALL equaPlane(posCMA(i,stack,:),posCMB(i,stack,:),posCMC(i,stack,:),normPlane(i,:),equat)
                distPi(i) = distInterPlane(equat,normPlane(i,:),pointOutOfPlane(i,:))
                dstack(i) = distStack(equat,pointOutOfPlane(1,:))
                tiltAngle(i) = 180*(1-ACOS(cosTiltAngle(normPlane(i,:),normBasePlane))/ACOS(-1.0))
                WRITE(1,*) (1+(i-1)*step_traj)*tstep, distPi(i)
                WRITE(2,*) (1+(i-1)*step_traj)*tstep, dstack(i)
                WRITE(3,*) (1+(i-1)*step_traj)*tstep, tiltAngle(i)
                dPimo = dPimo + distPi(i)
                deltaP = deltaP + distPi(i)**2
                dStamo = dStamo + dstack(i)
                deltaS = deltaS + dstack(i)**2
                tilAnmo = tilAnmo + tiltAngle(i)
                deltaA = deltaA + tiltAngle(i)**2
            ENDDO
        CLOSE(1)
        CLOSE(2)
        CLOSE(3)
        dPimo = dPimo/Nb_step
        deltaP = SQRT(deltaP/Nb_step - dPimo**2)
        dStamo = dStamo/Nb_step
        deltaS = SQRT(deltaS/Nb_step - dStamo**2)
        tilAnmo = tilAnmo/Nb_step
        deltaA = SQRT(deltaA/Nb_step - tilAnmo**2)
        open(unit=3, file = 'OUTPUT', status = 'old', position = "append")
        write(3,*) 'stacking distance (Angstrom) and pi-pi (Angstrom) distance and the tilt angle (°) have been calculated'
        write(3,*) 'and write in files distPi.out, distStack.out, TiltAngle.out'
        write(3,*) 'stack', stack
        write(3,*) 'distPi moy =',dPimo,'+/-',deltaP
        write(3,*) 'dstack moy =',dStamo,'+/-',deltaS
        write(3,*) 'costiltAn moy =',tilAnmo,'+/-',deltaA
    END SUBROUTINE stackToStack
END MODULE mod_analyzeTrajPlaneAlgo