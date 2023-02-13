!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
 Subroutine LUDCMP(A,N,INDX,D,CODE)
        integer :: N, I, IMAX, J, K
        double precision :: TINY=1.5D-16
        double precision :: AMAX,DUM, SUM, VV(100)
        double precision, dimension(8,8) :: A
        integer CODE, D, INDX(N)

        D=1; CODE=0

        DO I=1,N
            AMAX=0.d0
            DO J=1,N
                IF (DABS(A(I,J)).GT.AMAX) AMAX=DABS(A(I,J))
            END DO ! j loop
            IF(AMAX.LT.TINY) THEN
                CODE = 1
                RETURN
            END IF
            VV(I) = 1.d0 / AMAX
            END DO ! i loop

        DO J=1,N
            DO I=1,J-1
                SUM = A(I,J)
                DO K=1,I-1
                    SUM = SUM - A(I,K)*A(K,J) 
                END DO ! k loop
                A(I,J) = SUM
            END DO ! i loop
            AMAX = 0.d0
            DO I=J,N
                SUM = A(I,J)
                DO K=1,J-1
                    SUM = SUM - A(I,K)*A(K,J) 
                END DO ! k loop
                A(I,J) = SUM
                DUM = VV(I)*DABS(SUM)
                IF(DUM.GE.AMAX) THEN
                    IMAX = I
                    AMAX = DUM
                END IF
            END DO ! i loop  
           
            IF(J.NE.IMAX) THEN
                DO K=1,N
                    DUM = A(IMAX,K)
                    A(IMAX,K) = A(J,K)
                    A(J,K) = DUM
                END DO ! k loop
                D = -D
                VV(IMAX) = VV(J)
            END IF

            INDX(J) = IMAX
            IF(DABS(A(J,J)) < TINY) A(J,J) = TINY

            IF(J.NE.N) THEN
                DUM = 1.d0 / A(J,J)
                DO I=J+1,N
                    A(I,J) = A(I,J)*DUM
                END DO ! i loop
            END IF 
        END DO ! j loop

        RETURN
        END subroutine LUDCMP




program module_particles_test
implicit none
    double precision :: RXIJ, RYIJ, PHII, PHIJ
    double precision :: A1, A2, B1, B2, C1, C2, D1, D2, E1, E2, K, y, dy, dphi, kappa, Interval
    double precision :: a_axis, b_axis, b_axis_4, inv_a_axis_4, pi, j, a_axis_j, j_sq, b_sq
    double precision :: sn_dphi_sq, cs_phii_sq, cs_phii_pr_sq, cs_phij_sq, cs_phij_pr_sq
    double precision :: ju1, ju2, ju3, ju4, ju1_2, ju2_2, ju3_2, ju4_2, ju1_3, ju2_3, ju3_3, ju4_3, ju1_4, ju2_4, ju3_4, ju4_4 
    double precision :: lu1, lu2, lu3, lu4, lu1_2, lu2_2, lu3_2, lu4_2, lu1_3, lu2_3, lu3_3, lu4_3, lu1_4, lu2_4, lu3_4, lu4_4 
    integer :: P, nr_interv
    logical :: ovrlap
    double precision, dimension(8,8) :: mat

    integer :: N=8
    integer :: INDX(8)
    integer :: D, CODE,R
    double precision :: Det


    ovrlap = .false.
   
    ! *** Using Overlap Criterion by Vielliard-Baron *** !
    ! ** from: Phase Transitions of the Classical Hard-Ellipse System - The Journal of Chemical Physics, 1972, Vol 56, Issue 10 ** !                   

    pi=4.*atan(1.d0)
    RXIJ = 0.81469565472048355
    RYIJ = -0.71147581074881250
    PHII = -0.46793242549841968
    PHIJ = 0.86811917629350699
   
    !RXIJ = -8.87083E-002
    !RYIJ = 0.12435
   
    !PHII = -0.53917
    !PHIJ = -8.10069E-002
    a_axis = 1. 
    b_axis = 0.25
    kappa = 4
    nr_interv = 50



    b_axis_4 = b_axis**4.0d0
    inv_a_axis_4 = (1.0d0/a_axis)**4.0d0
    K = b_axis_4 * inv_a_axis_4 - 1.
    j = 0.5*sqrt(RXIJ**2+RYIJ**2)
    
    ju1 = 0.5*( RXIJ*cos(PHII)+RYIJ*sin(PHII))
    ju2 = 0.5*(-RXIJ*sin(PHII)+RYIJ*cos(PHII))
    ju3 = 0.5*( RXIJ*cos(PHIJ)+RYIJ*sin(PHIJ))
    ju4 = 0.5*(-RXIJ*sin(PHIJ)+RYIJ*cos(PHIJ))

    ju1_2 = (0.5*( RXIJ*cos(PHII)+RYIJ*sin(PHII)))**2.0d0
    ju2_2 = (0.5*(-RXIJ*sin(PHII)+RYIJ*cos(PHII)))**2.0d0
    ju3_2 = (0.5*( RXIJ*cos(PHIJ)+RYIJ*sin(PHIJ)))**2.0d0
    ju4_2 = (0.5*(-RXIJ*sin(PHIJ)+RYIJ*cos(PHIJ)))**2.0d0

    ju1_3 = (0.5*( RXIJ*cos(PHII)+RYIJ*sin(PHII)))**3.0d0
    ju2_3 = (0.5*(-RXIJ*sin(PHII)+RYIJ*cos(PHII)))**3.0d0
    ju3_3 = (0.5*( RXIJ*cos(PHIJ)+RYIJ*sin(PHIJ)))**3.0d0
    ju4_3 = (0.5*(-RXIJ*sin(PHIJ)+RYIJ*cos(PHIJ)))**3.0d0

    ju1_4 = (0.5*( RXIJ*cos(PHII)+RYIJ*sin(PHII)))**4.0d0
    ju2_4 = (0.5*(-RXIJ*sin(PHII)+RYIJ*cos(PHII)))**4.0d0
    ju3_4 = (0.5*( RXIJ*cos(PHIJ)+RYIJ*sin(PHIJ)))**4.0d0
    ju4_4 = (0.5*(-RXIJ*sin(PHIJ)+RYIJ*cos(PHIJ)))**4.0d0

    lu1 = 0.5*(-RYIJ*cos(PHII)+RXIJ*sin(PHII))
    lu2 = 0.5*( RYIJ*sin(PHII)+RXIJ*cos(PHII))
    lu3 = 0.5*(-RYIJ*cos(PHIJ)+RXIJ*sin(PHIJ))
    lu4 = 0.5*( RYIJ*sin(PHIJ)+RXIJ*cos(PHIJ))

    lu1_2 = (0.5*(-RYIJ*cos(PHII)+RXIJ*sin(PHII)))**2.0d0
    lu2_2 = (0.5*( RYIJ*sin(PHII)+RXIJ*cos(PHII)))**2.0d0
    lu3_2 = (0.5*(-RYIJ*cos(PHIJ)+RXIJ*sin(PHIJ)))**2.0d0
    lu4_2 = (0.5*( RYIJ*sin(PHIJ)+RXIJ*cos(PHIJ)))**2.0d0

    lu1_3 = (0.5*(-RYIJ*cos(PHII)+RXIJ*sin(PHII)))**3.0d0
    lu2_3 = (0.5*( RYIJ*sin(PHII)+RXIJ*cos(PHII)))**3.0d0
    lu3_3 = (0.5*(-RYIJ*cos(PHIJ)+RXIJ*sin(PHIJ)))**3.0d0
    lu4_3 = (0.5*( RYIJ*sin(PHIJ)+RXIJ*cos(PHIJ)))**3.0d0

    lu1_4 = (0.5*(-RYIJ*cos(PHII)+RXIJ*sin(PHII)))**4.0d0
    lu2_4 = (0.5*( RYIJ*sin(PHII)+RXIJ*cos(PHII)))**4.0d0
    lu3_4 = (0.5*(-RYIJ*cos(PHIJ)+RXIJ*sin(PHIJ)))**4.0d0
    lu4_4 = (0.5*( RYIJ*sin(PHIJ)+RXIJ*cos(PHIJ)))**4.0d0


    A1 = ju1_4 + K*ju1_4 + ju2_4
    A2 = ju3_4 + K*ju3_4 + ju4_4
    B1 = ju1_3*lu1 + K*ju1_3*lu1+ju2_3*lu2
    B2 = ju3_3*lu3 + K*ju3_3*lu3+ju4_3*lu4
    C1 = ju1_2*lu1_2 + K*ju1_2*lu1_2 + ju2_2*lu2_2
    C2 = ju3_2*lu3_2 + K*ju3_2*lu3_2 + ju4_2*lu4_2
    D1 = lu1_4 + K*lu1_4 + lu2_4
    D2 = lu3_4 + K*lu3_4 + lu4_4
    E1 = ju1*lu1_3 + K*ju1*lu1_3 + ju2*lu2_3
    E2 = ju3*lu3_3 + K*ju3*lu3_3 + ju4*lu4_3


    write(*,*) ju1_3, ju2_3, lu1_3, lu2_3
    write(*,*) A1, A2, B1, B2    
    !write(*,*) A1, A2, B1, B2, C1, C2
    

    Interval = (0.5*(1./kappa**2.+kappa**2.)**0.25)/j
    y = -Interval
    dy = 2*Interval * 1./nr_interv

    if (j < b_axis) then
         ovrlap = .true.
    else 
        P=0
        do while ((P .LE. nr_interv) .AND. (ovrlap .EQV. .false.))
           
            mat(1,1) = A1 
            mat(2,2) = A1 
            mat(3,3) = A1 
            mat(4,4) = A1
            mat(2,1) = 4*A1 +  4*y*B1 
            mat(3,2) = mat(2,1)
            mat(4,3) = mat(2,1)
            mat(5,4) = mat(2,1)
            mat(3,1) = 6*A1 + 12*y*B1 +  6*y**2*C1 
            mat(4,2) = mat(3,1)
            mat(5,3) = mat(3,1)
            mat(6,4) = mat(3,1)
            mat(4,1) = 4*A1 + 12*y*B1 + 12*y**2*C1 + 4*y**3*E1 
            mat(5,2) = mat(4,1)
            mat(6,3) = mat(4,1)
            mat(7,4) = mat(4,1)
            mat(5,1) = A1 +  4*y*B1 +  6*y**2*C1 + 4*y**3*E1 + y**4*D1 - b_axis**4
            mat(6,2) = mat(5,1)
            mat(7,3) = mat(5,1)
            mat(8,4) = mat(5,1)
            mat(6,1) = 0.
            mat(7,2) = 0.
            mat(8,3) = 0.
            mat(1,4) = 0.
            mat(7,1) = 0.
            mat(8,2) = 0.
            mat(1,3) = 0.
            mat(2,4) = 0. 
            mat(8,1) = 0.
            mat(1,2) = 0.
            mat(2,3) = 0.
            mat(3,4) = 0. 

            mat(1,5) = A2 
            mat(2,6) = A2
            mat(3,7) = A2 
            mat(4,8) = A2
            mat(2,5) = -4*A2 +  4*y*B2 
            mat(3,6) = mat(2,5)
            mat(4,7) = mat(2,5)
            mat(5,8) = mat(2,5)
            mat(3,5) = 6*A2 - 12*y*B2 +  6*y**2*C2
            mat(4,6) = mat(3,5)
            mat(5,7) = mat(3,5)
            mat(6,8) = mat(3,5)
            mat(4,5) = -4*A2 + 12*y*B2 - 12*y**2*C2 + 4*y**3*E2
            mat(5,6) = mat(4,5)
            mat(6,7) = mat(4,5)
            mat(7,8) = mat(4,5)
            mat(5,5) = A2 -  4*y*B2 +  6*y**2*C2 - 4*y**3*E2 + y**4*D2 - b_axis**4 
            mat(6,6) = mat(5,5)
            mat(7,7) = mat(5,5)
            mat(8,8) = mat(5,5)
            mat(6,5) = 0.
            mat(7,6) = 0.
            mat(8,7) = 0.
            mat(1,8) = 0.
            mat(7,5) = 0.
            mat(8,6) = 0.
            mat(1,7) = 0.
            mat(2,8) = 0.
            mat(8,5) = 0.
            mat(1,6) = 0.
            mat(2,7) = 0.
            mat(3,8) = 0. 

            call ludcmp(mat,N,INDX,D,CODE)
            if (D == 1) then 
                Det=1.
            else if ( D  == -1 ) then 
                Det=-1.
            end if

            do R=1,N
                Det=Det*mat(R,R)
            enddo

            if (Det .LT. 0) then
                ovrlap = .true.
            endif 

            write(*,*) y, Det       
            y = y+dy
            P = P+1   
        end do
    end if
    write (*,*) ovrlap, RXIJ, RYIJ, PHII, PHIJ
            !write (*,*) "A1=", A1, y, Det
    write (*,*) j, Det
    write(*,*) ovrlap

end

!            if ( (psi .gt. 0.0d0) .and. ((f1 .lt. 0.0d0) .or. (f2 .lt. 0.0d0)) ) then
!             if ( (psi .gt. dble_epsilon_2) .and. ((f1 .lt. -dble_epsilon_2) .or. (f2 .lt. -dble_epsilon_2)) ) then
!             if ( (abs(psi) .lt. dble_epsilon .or. psi .lt. 0.0d0) .and. ((abs(f1) .lt. dble_epsilon .or. f1 .lt. 0.0d0) .or. (abs(f2) .lt. dble_epsilon_2 .or. f2 .lt. 0.0d0)) ) then
!                ovrlap = .false.
!            else
!                ovrlap = .true.
!            end if
