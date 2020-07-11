!
! Driver programme for USR-technique
! Update: 5 Mar 2019
!


        PROGRAM main
        USE initialise      ,ONLY: read_session
        USE simulate        ,ONLY: start_serial_session
        IMPLICIT NONE

        PRINT *, '**********************************'
        PRINT *, 'Driver programme for USR-technique'
        PRINT *, '**********************************'
        PRINT *, ''

        CALL read_session       !read session input file
        CALL start_serial_session

        END PROGRAM main
