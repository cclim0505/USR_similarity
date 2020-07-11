        MODULE simulate
        USE constants       ,ONLY: DBL
        USE initialise      ,ONLY: atoms, coord, read_coord, &
     &    ref_struc_file, target_struc_file
        USE simil           ,ONLY: ref_mu, target_mu, timestep, &
     &    compute_struc_moments, calc_similarity, output_similarity


        PUBLIC   :: start_single_session ! For tests and tryouts
        PUBLIC   :: start_serial_session ! Main subroutine

        CONTAINS

        SUBROUTINE start_single_session
! Start single computation session that doesn't loop over all xyz
! inputs.
! Subroutine for tests and try outs.

        IMPLICIT NONE
        CALL read_coord(ref_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(ref_mu)
        CALL read_coord(target_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(target_mu)
        CALL calc_similarity
        CALL output_similarity
        END SUBROUTINE start_single_session

        SUBROUTINE start_serial_session
! Start serial computation session.
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
          CALL output_similarity
          
        END DO

        PRINT *, ''
        PRINT *, 'counter', counter
        PRINT *, ''
        END SUBROUTINE start_serial_session

        END MODULE simulate
