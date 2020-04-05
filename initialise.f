        MODULE initialise
        USE     constants
!=========================================================================================
! VARIABLE DICTIONARY
!=========================================================================================
! Simulation parameters
        INTEGER         :: atoms
        REAL(KIND=DBL),DIMENSION(:,:), ALLOCATABLE            :: coord
        CHARACTER(20)       :: session_file='session_in.dat'
        CHARACTER(30)       :: ref_struc_file
        CHARACTER(30)       :: target_struc_file
        CHARACTER(30)       :: output_file
! Global calculation values (coordinates of mu points and mu values)
        REAL(KIND=DBL),DIMENSION(3)  :: centroid                ! centroid
        REAL(KIND=DBL),DIMENSION(3)  :: close_ctd, far_ctd      ! closest and farthest atoms to centroid
        REAL(KIND=DBL),DIMENSION(3)  :: far_far                 ! farthest atom to farthest atom (f2f)
        REAL(KIND=DBL),DIMENSION(16)  :: ref_mu, target_mu      ! array to hold mu values of reference and target structures
        REAL(KIND=DBL)               :: similarity              ! similarity index
        INTEGER                      :: timestep
!=========================================================================================
        CONTAINS

        SUBROUTINE read_session
! read session input file
        IMPLICIT NONE
        CHARACTER       :: dummy
        OPEN(20,file=session_file,status='old')
        READ(20,*) dummy, ref_struc_file
        READ(20,*) dummy, target_struc_file
        READ(20,*) dummy, output_file
        CLOSE(20)
        END SUBROUTINE read_session
        
        SUBROUTINE start_single_session
! start single computation session
        IMPLICIT NONE
        CALL read_coord(ref_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(ref_mu)
        CALL read_coord(target_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(target_mu)
        CALL calc_similarity
        CALL output_results
        END SUBROUTINE start_single_session

        SUBROUTINE start_serial_session
! start serial computation session
        IMPLICIT NONE
        INTEGER         :: ierr !file iostat ierror
        INTEGER         :: iter
        INTEGER         :: counter = 0
        CHARACTER       :: dummy

        CALL read_coord(ref_struc_file)
        CALL compute_struc_moments(ref_mu)
!        timestep = 0 


        OPEN(40,file=target_struc_file,status='old')
        DO
          READ(40,*,iostat=ierr) atoms
          IF(ierr /= 0) EXIT
          counter = counter + 1
!          timestep = timestep + 50
          READ(40,*) dummy, dummy, timestep
          IF (ALLOCATED(coord)) DEALLOCATE(coord)
          ALLOCATE(coord(atoms,3))
          DO iter=1,atoms
            READ(40,*) dummy,coord(iter,1),coord(iter,2),coord(iter,3)
          END DO

          CALL compute_struc_moments(target_mu)
          CALL calc_similarity
          CALL output_results
          
        END DO

        PRINT *, ''
        PRINT *, 'counter', counter
        PRINT *, ''



        END SUBROUTINE start_serial_session

        SUBROUTINE read_coord(struc_file)
! read xyz coordinates file
        IMPLICIT NONE
        CHARACTER(20),INTENT(IN)   :: struc_file
        INTEGER         :: iter
        CHARACTER       :: dummy
        OPEN(21,file=struc_file,status='old')
        READ(21,*) atoms
        READ(21,*) dummy
        IF (ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))
        DO iter=1,atoms
          READ(21,*) dummy, coord(iter,1),coord(iter,2),coord(iter,3)
        END DO
        CLOSE(21)
        END SUBROUTINE read_coord
        
        SUBROUTINE print_coord
! print xyz coordinate files in columns
        IMPLICIT NONE
        INTEGER  :: iter
        PRINT *, ''
        PRINT *, 'Coordinates in xyz'
        PRINT *, ''
        DO iter =1,atoms 
          PRINT *, iter, coord(iter,1), coord(iter,2), coord(iter,3)
        END DO
        PRINT *, ''
        END SUBROUTINE print_coord


        SUBROUTINE calc_centroid
! calculate centroid (1st mu_point)
        IMPLICIT NONE
        INTEGER :: iter
        centroid (:) = 0.0
!        PRINT *, ''
!        PRINT *, 'initial centroid value'
!        PRINT *, centroid(1), centroid(2), centroid(3)
!        PRINT *, ''
        DO iter=1,atoms
          centroid(1) = centroid(1) + coord(iter,1)
          centroid(2) = centroid(2) + coord(iter,2)
          centroid(3) = centroid(3) + coord(iter,3)
        END DO

        centroid(:) = centroid(:) / REAL(atoms)

!        PRINT *, ''
!        PRINT *, 'centroid =', centroid(1), centroid(2), centroid(3)
!        PRINT *, ''
        END SUBROUTINE calc_centroid


        SUBROUTINE locate_mu_points
! locate (1) closest atom to centroid (2nd mu_point)
! locate (2) farthest atom from centroid (3rd mu_point)
! locate (3) farthest atom from the farthest (4th mu_point)
        IMPLICIT NONE
        INTEGER :: iter
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE     :: centroid_diff !difference vectors
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE       :: centroid_dist !distance vectors
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE     :: f2f_diff !difference vectors for far to far
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE       :: f2f_dist !distance vectors for far to far
        REAL(KIND=DBL)  :: close_dist, far_dist     ! closest and farthest distance from centroid
        REAL(KIND=DBL)  :: far_far_dist                 ! farthest to farthest distance
        REAL(KIND=DBL)  :: difference
        REAL(KIND=DBL)  :: interval = 1E-10_DBL

        IF (ALLOCATED(centroid_diff)) DEALLOCATE (centroid_diff)
        IF (ALLOCATED(centroid_dist)) DEALLOCATE (centroid_dist)
        IF (ALLOCATED(f2f_diff)) DEALLOCATE (f2f_diff)
        IF (ALLOCATED(f2f_dist)) DEALLOCATE (f2f_dist)
        ALLOCATE (centroid_diff(atoms,3))
        ALLOCATE (centroid_dist(atoms))
        ALLOCATE (f2f_diff(atoms,3))
        ALLOCATE (f2f_dist(atoms))
!********************************
! Locate closest and farthest atoms from centroid
!********************************
        DO iter=1,3
          centroid_diff(:,iter) = coord(:,iter) - centroid(iter)
        END DO

!        PRINT *, ''
!        PRINT *, 'centroid_diff'
!        PRINT *, ''
!        DO iter =1,atoms
!          PRINT *, centroid_diff(iter,1), centroid_diff(iter,2) &
!     &    , centroid_diff(iter,3)
!        END DO
        DO iter =1,atoms
          centroid_dist(iter) = NORM2(centroid_diff(iter,:))
        END DO
!        PRINT *, ''
!        PRINT *, 'centroid_dist'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, iter, centroid_dist(iter)
!        END DO
        close_dist = MINVAL(centroid_dist)
        far_dist = MAXVAL(centroid_dist)
!        PRINT *, ''
!        PRINT *, 'closest and farthest distance'
!        PRINT *, ''
!        PRINT *, close_dist, far_dist
        DO iter =1,atoms
          difference = close_dist - centroid_dist(iter)
          IF (ABS(difference) < interval) THEN
            close_ctd = coord(iter,:) 
!            PRINT *, 'closest index', iter
            EXIT
          END IF
        END DO
        DO iter =1,atoms
          difference = far_dist - centroid_dist(iter)
          IF (ABS(difference) < interval) THEN
            far_ctd = coord(iter,:) 
!            PRINT *, 'farthest index', iter
!            PRINT *, far_ctd(1), far_ctd(2), far_ctd(3)
            EXIT
          END IF
        END DO

!********************************
! Locate farthest to farthest atom
!********************************
        DO iter=1,3
          f2f_diff(:,iter) = coord(:,iter) - far_ctd(iter) 
        END DO
        DO iter =1,atoms
          f2f_dist(iter) = NORM2(f2f_diff(iter,:))
        END DO
!        PRINT *, ''
!        PRINT *, 'f2f_dist'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, iter, f2f_dist(iter)
!        END DO
        far_far_dist= MAXVAL(f2f_dist)
!        PRINT *, ''
!        PRINT *, 'farthest to farthest distance'
!        PRINT *, ''
!        PRINT *, far_far_dist
        
        DO iter =1,atoms
          difference = far_far_dist - f2f_dist(iter)
          IF (ABS(difference) < interval) THEN
            far_far = coord(iter,:) 
!           PRINT *, 'f2f index', iter
!            PRINT *, far_far(1), far_far(2), far_far(3)
            EXIT
          END IF
        END DO

        DEALLOCATE (centroid_diff)
        DEALLOCATE (centroid_dist)
        DEALLOCATE (f2f_diff)
        DEALLOCATE (f2f_dist)
        END SUBROUTINE locate_mu_points

        SUBROUTINE calc_mu_val(mu_point,mu_val)
! calculate four mu value for a certain mu_point
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(3),INTENT(IN)  :: mu_point        !centroid, closest, farthest or far_far
        REAL(KIND=DBL),DIMENSION(4),INTENT(OUT) :: mu_val
        INTEGER         :: iter
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE     :: diff !difference vectors
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE       :: dist !distance of vectors
        REAL(KIND=DBL)          :: mu_1, mu_2, mu_3, mu_4
        REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE       :: term_1, term_2
        REAL(KIND=DBL)          :: sum_1, sum_2
        REAL(KIND=DBL)          :: std_dev
        IF (ALLOCATED(diff)) DEALLOCATE(diff)
        IF (ALLOCATED(dist)) DEALLOCATE(dist)
        IF (ALLOCATED(term_1)) DEALLOCATE(term_1)
        IF (ALLOCATED(term_2)) DEALLOCATE(term_2)
        ALLOCATE (diff(atoms,3))
        ALLOCATE (dist(atoms))
        ALLOCATE (term_1(atoms))
        ALLOCATE (term_2(atoms))

        DO iter=1,3
          diff(:,iter) = coord(:,iter) - mu_point(iter)
        END DO
        DO iter =1,atoms
          dist(iter) = NORM2(diff(iter,:))
        END DO

!        PRINT *, ''
!        PRINT *, 'mu_point_dist'
!        PRINT *, ''
!        DO iter=1,atoms
!          PRINT *, iter, dist(iter)
!        END DO
 
!********************************
! Calculate mu_1
!********************************
        mu_1 = SUM(dist) / REAL(atoms)
!       PRINT *, 'mu_1', mu_1
!********************************
! Calculate mu_2
!********************************
        term_1(:) = dist(:) - mu_1

        sum_2 = SUM(term_1)
        sum_2 = sum_2**2
        sum_2 = sum_2 / REAL(atoms)


        term_1(:) = term_1(:)**2
        sum_1 = SUM(term_1)

        mu_2 = sum_1 - sum_2
        mu_2 = mu_2 / REAL(atoms-1)

!        PRINT *, 'mu_2', mu_2

!********************************
! Calculate mu_3
!********************************
        std_dev = DSQRT(mu_2)
        
        term_1(:) = dist(:) - mu_1
        term_1(:) = term_1(:) / std_dev
        term_2(:) = term_1(:)           ! for mu_4 calculation
        term_1(:) = term_1(:) ** 3
        sum_1 = SUM(term_1)
        mu_3 = sum_1 / REAL(atoms)

!        PRINT *, 'mu_3', mu_3
!********************************
! Calculate mu_4
!********************************
        term_2(:) = term_2(:)**4
        sum_2 = SUM(term_2)
        mu_4 = sum_2 / REAL(atoms)
        mu_4 = mu_4 - 3.0_DBL

!        PRINT *, 'mu_4', mu_4

        mu_val(1) = mu_1
        mu_val(2) = mu_2
        mu_val(3) = mu_3
        mu_val(4) = mu_4

        END SUBROUTINE calc_mu_val


        SUBROUTINE compute_struc_moments(moment)
! calculate target structure's 16 moment values
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(16),INTENT(OUT)  :: moment       ! array to hold mu values of reference and target structures
        REAL(KIND=DBL),DIMENSION(4)      :: mu_val

        CALL calc_centroid
        CALL locate_mu_points
        CALL calc_mu_val(centroid,mu_val)
        moment(1) = mu_val(1)
        moment(2) = mu_val(2)
        moment(3) = mu_val(3)
        moment(4) = mu_val(4)
        CALL calc_mu_val(close_ctd,mu_val)
        moment(5) = mu_val(1)
        moment(6) = mu_val(2)
        moment(7) = mu_val(3)
        moment(8) = mu_val(4)
        CALL calc_mu_val(far_ctd,mu_val)
        moment(9) = mu_val(1)
        moment(10) = mu_val(2)
        moment(11) = mu_val(3)
        moment(12) = mu_val(4)
        CALL calc_mu_val(far_far,mu_val)
        moment(13) = mu_val(1)
        moment(14) = mu_val(2)
        moment(15) = mu_val(3)
        moment(16) = mu_val(4)

!        PRINT *, ''
!        PRINT *, 'moments'
!        PRINT *, ''
!        DO iter=1,16
!          PRINT *, iter, moment(iter)
!        END DO

        END SUBROUTINE compute_struc_moments

        SUBROUTINE calc_similarity
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(16)  :: denom
        REAL(KIND=DBL)                :: denom_sum

        denom = ABS(ref_mu - target_mu)
        denom_sum = SUM(denom)
        denom_sum = denom_sum / 16.0_DBL
        denom_sum = denom_sum + 1.0_DBL
        denom_sum = 1.0_DBL / denom_sum
        similarity = denom_sum


        PRINT *, ''
        PRINT *, 'similarity index', similarity
        PRINT *, ''
        END SUBROUTINE calc_similarity
        
        SUBROUTINE output_results
        IMPLICIT NONE
        OPEN(30,file=output_file,access='append')
        WRITE(30,*) timestep,similarity
        CLOSE(30)
        END SUBROUTINE output_results

        END MODULE initialise
