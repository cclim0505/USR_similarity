        MODULE simulate
        USE constants       ,ONLY: DBL
        USE initialise      ,ONLY: atoms, coord, energy, &
     &    read_coord, &
     &    ref_struc_file, target_struc_file
        USE simil           ,ONLY: ref_mu, target_mu, timestep, &
     &    compute_struc_moments, calc_similarity, output_similarity, &
     &    similarity


        PUBLIC   :: start_serial_session ! Similarity calculation
        PUBLIC   :: start_uniq_session   ! Find unique structures

        PUBLIC   :: start_single_session ! For tests and tryouts

        CONTAINS


        SUBROUTINE start_serial_session
! Start serial computation session.
        IMPLICIT NONE
        INTEGER         :: ierr !file iostat ierror
        INTEGER         :: iter
        INTEGER         :: counter = 0
        CHARACTER       :: dummy
        REAL(KIND=DBL),DIMENSION(16)    :: mu_ref, mu_target

! Calculate mu value of reference structure
        CALL read_coord(ref_struc_file)
        CALL compute_struc_moments(mu_ref)

! Calculate mu value of target structures in a loop
        OPEN(40,file=target_struc_file,status='old')
        DO
          READ(40,*,iostat=ierr) atoms
          IF(ierr /= 0) EXIT
          counter = counter + 1
          READ(40,*) dummy, dummy, dummy, timestep
          IF (ALLOCATED(coord)) DEALLOCATE(coord)
          ALLOCATE(coord(atoms,3))
          DO iter=1,atoms
            READ(40,*) dummy,coord(iter,1),coord(iter,2),coord(iter,3)
          END DO

          CALL compute_struc_moments(mu_target)
          CALL calc_similarity(mu_ref, mu_target)
          CALL output_similarity
          
        END DO

        PRINT *, ''
        PRINT *, 'counter', counter
        PRINT *, ''

        END SUBROUTINE start_serial_session

        SUBROUTINE start_uniq_session
! Simulation that finds all the unique structures in the set.
        USE uniq            ,ONLY: ref_mu_list, init_uniq_var, &
     &    uniq_count, uniq_criteria, save_uniq_struc, &
     &    get_first_uniq
        IMPLICIT NONE
        INTEGER         :: f_struc
        INTEGER         :: ierr !file iostat ierror
        INTEGER         :: iter
        CHARACTER       :: dummy
        LOGICAL         :: found_uniq = .FALSE.
        REAL(KIND=DBL),DIMENSION(16)    :: mu_ref, mu_target

        CALL init_uniq_var
        CALL get_first_uniq

! Open coord file.
        OPEN(NEWUNIT=f_struc,file=target_struc_file,status='old')

        DO
! Read target structures and calculate mu value.
          READ(f_struc,*,iostat=ierr) atoms
          IF(ierr /= 0) EXIT
          READ(f_struc,*) energy, dummy, dummy, timestep
          IF (ALLOCATED(coord)) DEALLOCATE(coord)
          ALLOCATE(coord(atoms,3))
          DO iter=1,atoms
            READ(f_struc,*) dummy, &
     &        coord(iter,1),coord(iter,2),coord(iter,3)
          END DO

          CALL compute_struc_moments(mu_target)

! Calc similarity between target and unique structures.
! Save structure when unique criteria is met.
          found_uniq = .TRUE.
          DO iter=1, uniq_count
            mu_ref = ref_mu_list(iter,:)
            CALL calc_similarity(mu_ref, mu_target)
            IF (similarity >= uniq_criteria) THEN
              found_uniq = .FALSE.
              EXIT
            END IF
          END DO

          IF (found_uniq) THEN
            uniq_count = uniq_count + 1
            ref_mu_list(uniq_count,:) = mu_target
            CALL save_uniq_struc(coord)
          END IF
          
        END DO

        CLOSE(f_struc)

        PRINT*, "Uniq structures are:", uniq_count

        END SUBROUTINE start_uniq_session



        SUBROUTINE start_single_session
! Start single computation session that doesn't loop over all xyz
! inputs.
! Subroutine for tests and try outs.
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(16)    :: mu_ref, mu_target
        CALL read_coord(ref_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(mu_ref)
        CALL read_coord(target_struc_file)
!        CALL print_coord
        CALL compute_struc_moments(mu_target)
        CALL calc_similarity(mu_ref, mu_target)
        CALL output_similarity
        END SUBROUTINE start_single_session

        END MODULE simulate
