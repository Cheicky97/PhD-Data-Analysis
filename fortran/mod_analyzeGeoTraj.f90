MODULE mod_analyseGeoalgo2
 IMPLICIT NONE
 CONTAINS
    SUBROUTINE interatmdist(traj,lenTraj,Nb_stack,Nb_step,monomer,tstep,step_traj)
        INTEGER, DIMENSION(3), INTENT(IN)                       :: monomer
        INTEGER, INTENT(IN)                                     :: lenTraj, Nb_stack, Nb_step, step_traj
        DOUBLE PRECISION, INTENT(IN)                            :: tstep
        DOUBLE PRECISION, DIMENSION(lenTraj,3), INTENT(IN) 	    :: traj
        DOUBLE PRECISION, DIMENSION(1,Nb_stack-1)			    :: dist 
        DOUBLE PRECISION, DIMENSION(Nb_stack)	        	    :: z
        INTEGER							                        :: i,j,k,l
        z = 0.0
        OPEN(1,file="interatmdist.out")
            DO l=0, Nb_step-1
                DO i=1,Nb_stack
                    DO j=1+monomer(1)*(i-1)+l*(lenTraj), monomer(1)*i+l*(lenTraj)
                        z(i) = z(i) + traj(j,3) 
                    ENDDO
                    k=monomer(1)
                    DO j=1+(Nb_stack*monomer(1))+monomer(2)*(i-1)+l*(lenTraj),&
                     (Nb_stack*monomer(1))+monomer(2)*i+l*(lenTraj)
                        z(i) = z(i) + traj(j,3)
                    ENDDO
                ENDDO
                z = z/(monomer(2)+monomer(1))
                DO i =1,Nb_stack-1
                    dist(1,i) = ABS(z(i)-z(i+1))
                ENDDO
                WRITE(1,*) 1+(l*step_traj)*tstep*1e-3, dist(1,:)
            ENDDO
        CLOSE(1)
        open(unit=2, file = 'OUTPUT', status = 'old', position = "append")
        write(2,*) '_________________________________________________________________________________'
        write(2,*) '   calculation of the interstack distance using the CM of the differents stack   '
        write(2,*) '_________________________________________________________________________________'
        write(2,*) '    file created : interatmdist.out (time(ps) | dist(stack1,stack2) | ... |)     '
        write(2,*) '                                                                                 '
        close(2)
    END SUBROUTINE interatmdist
END MODULE mod_analyseGeoalgo2