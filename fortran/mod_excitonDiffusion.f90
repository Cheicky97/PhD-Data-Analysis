MODULE mod_excitonDIFf
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE pasX(traj,lenTrajX,axe,Nb_step)
        ! lenTrajX : toujour = 2
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajX,3), INTENT(IN) :: traj !traj exciton
        INTEGER, INTENT(IN)						                    :: lenTrajX, axe, Nb_step !axe : 1=x, 2=y, 3=z
        DOUBLE PRECISION						                    :: d_h, d_e, meand_h, meand_e
        INTEGER							                            :: i,j,n
        n=1; meand_h = 0.0; meand_e = 0.0
        OPEN(unit = 1, file ='pasElecHol.out')
            DO i=3, Nb_step*lenTrajX-1,2
                j=i+1
                !__pour le h
                d_h = ABS(traj(i,axe)-traj(i-2,axe))
                !__pour l'e
                d_e = ABS(traj(j,axe)-traj(j-2,axe))
                meand_h = meand_h + d_h; meand_e = meand_e + d_e
                WRITE(1,*) n, d_h, d_e
                n=n+1
            ENDDO
        CLOSE(1)
        open(UNIT = 2, FILE = 'OUTPUT', STATUS = 'old', POSITION = "append")
        write(2,*) '>'
        write(2,*) "The distance between two successive steps of the electron&
        (d_e) and the hole (d_h) saved in 'pasElecHol.out'"
        write(2,*) '    d_h (Angstrom) =', meand_h/(Nb_step-1)
        write(2,*) '    d_e (Angstrom) =', meand_e/(Nb_step-1)
        write(2,*) '                                                                  '
        close(2)
    END SUBROUTINE pasX

    !######################################################################################
    !#        Première version de SUBROUTINE pour mettre à jour trajectoire de l'e et h   #
    !#            en s'affranchissant des condition aux limites périodiques               #
    !######################################################################################
    SUBROUTINE majTrajeh1(traj,lenTrajX,majTrajeh,axe,Nb_step,dimBox)
        DOUBLE PRECISION, INTENT(IN), DIMENSION(3)                               :: dimBox
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajX,3), INTENT(IN)	    :: traj !traj exciton
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajX,3), INTENT(INOUT)  :: majTrajeh !traj exciton
        INTEGER, INTENT(IN)						                        :: axe,lenTrajX,Nb_step ! 1=x, 2=y, 3=z
        DOUBLE PRECISION						                        :: d, dim_axe
        INTEGER							                                :: i,j,indx_h,indx_e
        indx_h =0; indx_e=0
        IF(axe==1) THEN
            d = dimBox(1)/2.0
            dim_axe = dimBox(1)
        ENDIF
        IF (axe==2) THEN
            d= dimBox(2)/2.0
            dim_axe = dimBox(2)
        ENDIF
        IF (axe==3) THEN
            d = dimBox(3)/2.0
            dim_axe = dimBox(3)
        ENDIF
        majTrajeh = traj
        !__mise-à-jour deplacement selon axe____ 
        DO i=3, lenTrajX-1,2
            j=i+1
            !__pour le h
            IF (ABS(traj(i,axe)-traj(i-2,axe)) >= d) THEN
                IF (traj(i-2,axe)<traj(i,axe)) THEN
                    indx_h = indx_h-1
                    majTrajeh(i,axe) = traj(i,axe)+indx_h*dim_axe
                ELSE
                    indx_h = indx_h+1
                    majTrajeh(i,axe) = traj(i,axe)+indx_h*dim_axe
                ENDIF
            ELSE
                majTrajeh(i,axe) = traj(i,axe)+indx_h*dim_axe
            ENDIF
            !__pour l'e
            IF (ABS(traj(j,axe)-traj(j-2,axe)) >= d) THEN
                IF (traj(j-2,axe)<traj(j,axe)) THEN
                    indx_e = indx_e-1
                    majTrajeh(j,axe) = traj(j,axe)+indx_e*dim_axe
                ELSE
                    indx_e = indx_e+1
                    majTrajeh(j,axe) = traj(j,axe)+indx_e*dim_axe
                ENDIF
            ELSE
                majTrajeh(j,axe) = traj(j,axe)+indx_e*dim_axe
            ENDIF
        ENDDO
    END SUBROUTINE majTrajeh1

    !#########################################################################################################
    !#           Ce SUBROUTINE permet de s'affranchir des conditions aux limites périodiques                 #
    !#                          met à jour trajectoires de l'exciton                                         # 
    !#  #criterAxe = la distance max entre deux positions successives d'e ou h à partir de laquelle          #
    !#                on considère qu'il y a saut de boîte                                                   #
    !#########################################################################################################
    SUBROUTINE majTrajeh2(oldTrajeh,newTrajeh,lenTraj,Nb_step,axe,dim_axe)
        INTEGER, INTENT(IN)                                                     :: lenTraj
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTraj,3), INTENT(IN)              :: oldTrajeh !traj exciton
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTraj,3),INTENT(INOUT)            :: newTrajeh
        INTEGER, INTENT(IN)                                                     :: axe, Nb_step
        DOUBLE PRECISION, INTENT(IN)                                            :: dim_axe
        DOUBLE PRECISION                                                        :: eps
        INTEGER                                                                 :: i,j,indx_h, indx_e, indx
        indx_h = 0; indx_e = 0
        eps = dim_axe/2.0
        newTrajeh = oldTrajeh
        DO i=2,Nb_step !boucle sur step
            IF (i > 1) THEN
                DO j=1,lenTraj
                    IF ((ABS(oldTrajeh(j+(i-1)*lenTraj,axe)-oldTrajeh(j+(i-2)*lenTraj,axe)) < eps)) THEN
                        IF (oldTrajeh(j+(i-1)*lenTraj,axe) < oldTrajeh(j+(i-2)*lenTraj,axe)) THEN
                            IF (j==1) THEN
                                indx_h = indx_h + 1
                                indx = indx_h
                            ELSE
                                indx_e = indx_e + 1
                                indx = indx_e
                            ENDIF
                        ELSE
                            IF (j==1) THEN
                                indx_h = indx_h - 1
                                indx = indx_h
                            ELSE
                                indx_e = indx_e - 1
                                indx = indx_e
                            ENDIF
                        ENDIF
                    ELSE
                        IF (j==1) THEN
                            indx = indx_h
                        ELSE
                            indx = indx_e
                        ENDIF
                    ENDIF
                    newTrajeh(j+(i-1)*lenTraj,axe) = oldTrajeh(j+(i-1)*lenTraj,axe) + indx*dim_axe
                ENDDO
            ENDIF
        ENDDO
    END SUBROUTINE MajTrajeh2

    !#####################################################################################################
    !#               calcul du produit scalaire de deux point x et x0                                    #
    !#####################################################################################################
    DOUBLE PRECISION FUNCTION msd(x,x0)
        INTEGER						:: i,j
        DOUBLE PRECISION, DIMENSION(1,3), INTENT(IN)		:: x0,x !coord de l'elmnt à un instant DOnné
        msd = SUM((x-x0)**2)
    END FUNCTION msd

    !##################################################################################################
    !# SUBROUTINE calcul de manière traditionelle le msd de l'exciton en fonction du temps            #
    !##################################################################################################
    SUBROUTINE msd1(traj,lenTrajX,axe,Nb_step,step_traj)
        INTEGER, INTENT(IN)						                                :: axe, Nb_step, lenTrajX,step_traj
        DOUBLE PRECISION, DIMENSION(Nb_step*lenTrajX,3), INTENT(IN)	            :: traj !traj exciton
        DOUBLE PRECISION, DIMENSION(Nb_step)				                    :: msd_h, msd_e
        DOUBLE PRECISION, DIMENSION(Nb_step)				                    :: msd_eh, dis
        !DOUBLE PRECISION						                                :: msd
        INTEGER							                                        :: i,j,k,l
        CHARACTER(1)							                                :: car
        j=1;k=0
        msd_h=0; msd_e=0.0
        IF (axe==0) THEN
            OPEN(unit=1,file="msd_eh_3d.res")
                DO i=1, Nb_step
                    k = j+1
                    msd_h(i) = msd(traj(j,:),traj(1,:)) 
                    msd_e(i) = msd(traj(k,:),traj(2,:))
                    msd_eh(i) = (msd_h(i)+msd_e(i))/2.0 
                    j = j+2
                    WRITE(1,*) 1+(i-1)*step_traj, msd_eh(i)
                ENDDO
            CLOSE(1)
            open(UNIT = 2, FILE = 'OUTPUT', STATUS = 'old', POSITION = "append")
            write(2,*) ">"
            write(2,*) "3D MSD of the exciton (X) calculated with MSD(t) = abs(x(t)-x(0))**2"
            write(2,*) "    File created : 'msd_eh_3d.res"
            close(2)
        ELSE
            IF (axe==1) car = "x"
            IF (axe==2) car = "y"
            IF (axe==3) car = "z"
            OPEN(unit=1,file="msd_eh_"//car//".res")
                DO i=1, Nb_step
                    k = j+1
                    msd_eh(i) = ((traj(j,axe)-traj(1,axe))**2 + (traj(k,axe)-traj(2,axe))**2)/2.0
                    j = j+2
                    WRITE(1,*) (i-1)*step_traj, msd_eh(i)
                ENDDO
            CLOSE(1)
            open(UNIT = 2, FILE = 'OUTPUT', STATUS = 'old', POSITION = "append")
            write(2,*) ">"
            write(2,*) "MSD along of the exciton (X) along ",car," calculated with MSD(t) = abs(x(t)-x(0))**2"
            write(2,*) "    File created : msd_eh_3d.res"
            close(2)
        ENDIF
    END SUBROUTINE msd1

    !##################################################################################################
    !#      SUBROUTINE pour calculer MSD, elle dIFfère de la manière classique.                       #
    !#      Dans cette formulation la trajectoire est échantillonée dIFfééremment                     #
    !#  Elle prEND en input la trajectoire #majTrajeh de eh mis à jour, une intervalle              #
    !#  sur laquelle calculer #IntervDIFf, avec un step #TotStepDIFf et selon un #axe DOnné           #
    !#  #nom est le nom du fichier qui sera généré et #len_nom le nombre de carac dans #nom           #
    !#  Propriétaire : Evelyne MARTIN                                                                 #
    !##################################################################################################
    SUBROUTINE MSD2(majTrajeh,lenTrajX,IntervDIFf,TotStepDIFf,axe,nom,len_nom,Nb_step,step_traj,tstep)
        INTEGER :: step,i,n,m,j,k
        INTEGER :: NbCalcMSD,NbCalcMSDPred
        INTEGER :: NorigCalcMSD, norig, nmin
        INTEGER :: Nbeh
        INTEGER :: axe
        INTEGER :: len_nom
        INTEGER, INTENT(IN) :: IntervDIFf,TotStepDIFf,Nb_step,step_traj,lenTrajX
        INTEGER, DIMENSION(:), allocatable :: InitStep
        DOUBLE PRECISION, INTENT(IN)                :: tstep
        DOUBLE PRECISION, DIMENSION(lenTrajX,3),INTENT(IN) :: majTrajeh !traj exciton
        DOUBLE PRECISION, DIMENSION (:,:), allocatable :: UnfoldedPos
        DOUBLE PRECISION, DIMENSION (:,:,:), allocatable :: InitPos
        DOUBLE PRECISION, DIMENSION(:), allocatable :: MSD_R_eh
        DOUBLE PRECISION, DIMENSION(:), allocatable :: DIFf_R_eh
        DOUBLE PRECISION	:: Temps
        CHARACTER(len_nom), INTENT(IN) :: nom
        ! initialize statistics and other:
        k=0
        NbCalcMSD       = 0
    
        IF (IntervDIFf > 0) NbCalcMSDPred   = (Nb_step-TotStepDIFf*IntervDIFf)/IntervDIFf+1
        ! Dynamic allocations 
        allocate( UnfoldedPos(2,3))
        IF (IntervDIFf> 0) THEN
            allocate(InitPos(2,3,TotStepDIFf+1))
            allocate(InitStep(TotStepDIFf+1)) 
            allocate(MSD_R_eh(TotStepDIFf))     ! tableaux pour garder le  MSD du eh
            allocate(DIFf_R_eh(TotStepDIFf))    ! tableaux pour garder le  MSD du eh
            InitPos      = 0.d0
            DIFf_R_eh    = 0.d0
            MSD_R_eh     = 0.d0
            Nbeh         = 2
        ENDIF
    
        DO step=1,Nb_step !time
            UnfoldedPos = majTrajeh(step+k:step+k+1,:)
            k=k+1
            !  Compute mean square displacements for the dIFfusivity
            IF (IntervDIFf > 0 ) THEN
                IF ( mod(step,IntervDIFf) == 0) THEN
                    IF ( (step+(TotStepDIFf-1)*IntervDIFf) <= Nb_step) THEN
                        NbCalcMSD = NbCalcMSD + 1
                        NorigCalcMSD = mod(NbCalcMSD-1,TotStepDIFf)+1
                        InitPos(:,:,NorigCalcMSD) = UnfoldedPos
                        InitStep(NorigCalcMSD) = step
                    ENDIF
                    IF (NbCalcMSD<TotStepDIFf) THEN
                        nmin = 1
                    ELSE
                        nmin = NbCalcMSD-TotStepDIFf+1
                    ENDIF      
                    DO n=nmin,NbCalcMSD
                        norig = mod(n-1,TotStepDIFf)+1            
                        IF ( (InitStep(norig)+(TotStepDIFf-1)*IntervDIFf) >= step) THEN 
                            m = (step - InitStep(norig))/IntervDIFf+1   
                            DO i=1,2
                                IF (axe == 0 ) THEN
                                    DO j =1,3 
                                        DIFf_R_eh(m) = DIFf_R_eh(m)+(UnfoldedPos(i,j)-InitPos(i,j,norig))**2 !(UnfoldedPos(i,j)-InitPos(i,j,norig))
                                    ENDDO
                                ELSE 
                                    DIFf_R_eh(m) = DIFf_R_eh(m)                 &
                                    +(UnfoldedPos(i,axe)-InitPos(i,axe,norig))*  &
                                    (UnfoldedPos(i,axe)-InitPos(i,axe,norig))
                                ENDIF
                            ENDDO   
                        ENDIF
                    ENDDO      
                ENDIF       
            ENDIF
        ENDDO !time
        
        IF ((IntervDIFf > 0) .and. (Nb_step>=TotStepDIFf*IntervDIFf)) THEN
            IF (NbCalcMSD/=NbCalcMSDPred) WRITE(6,*)'DIFf : Pb avec le nombre d origines',NbCalcMSD,NbCalcMSDPred
                                
            OPEN(unit=3,file='msd_'//nom//'.time') 
                DO m=1,TotStepDIFf
                    Temps = (m-1)*step_traj*IntervDIFf*tstep*1e-3
                    MSD_R_eh(m) = DIFf_R_eh(m)/(NbCalcMSDPred)/Nbeh
                    WRITE(3,'(4E15.7)') Temps, MSD_R_eh(m)
                ENDDO
            CLOSE(3)   
            open(UNIT = 2, FILE = 'OUTPUT', STATUS = 'old', POSITION = "append")
            write(2,*) ">"
            write(2,*) "MSD of the exciton (X) calculated with algorithm 2"
            write(2,*) '    File created : msd_'//nom//'.time'
            close(2)        
        ENDIF       
        deallocate(UnfoldedPos)
        IF (IntervDIFf > 0) THEN
            deallocate(InitPos,InitStep)
            deallocate(DIFf_R_eh,MSD_R_eh)
        ENDIF
    END SUBROUTINE MSD2

END MODULE mod_excitonDIFf