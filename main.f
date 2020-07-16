!
! Driver programme for USR-technique
! Update: 16 July 2020
!


        PROGRAM main
        USE initialise      ,ONLY: read_session, job_type
        USE simulate        ,ONLY: start_serial_session, &
     &    start_uniq_session
        IMPLICIT NONE

        PRINT *, '**********************************'
        PRINT *, 'Driver programme for USR-technique'
        PRINT *, '**********************************'
        PRINT *, ''

        CALL read_session       !read session input file

        IF (job_type == 1) THEN
          CALL start_serial_session
        
        ELSE IF (job_type == 2) THEN
          CALL start_uniq_session

        ELSE
          CALL start_serial_session

        END IF


        END PROGRAM main
