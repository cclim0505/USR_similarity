!
! Driver programme for USR-technique
! Update: 5 Mar 2019
!


        PROGRAM main
        USE initialise
        IMPLICIT NONE

        PRINT *, '**********************************'
        PRINT *, 'Driver programme for USR-technique'
        PRINT *, '**********************************'
        PRINT *, ''

        CALL read_session       !read session input file
!        CALL start_single_session
        CALL start_serial_session
!        CALL print_coord

!        CALL reference_struc_moments
!        CALL target_struc_moments
!        CALL similarity_index



!        CALL calc_centroid
!        CALL locate_mu_points
!        CALL calc_mu_val(centroid)
!        CALL calc_mu_val(close_ctd)
!        CALL calc_mu_val(far_ctd)
!        CALL calc_mu_val(far_to_far)
        END PROGRAM main
