        MODULE initialise
        USE     constants       ,ONLY: DBL
!=========================================================================================
! VARIABLE DICTIONARY
!=========================================================================================
! Simulation parameters
        INTEGER             :: atoms
        CHARACTER(2)        :: material
        REAL(KIND=DBL),DIMENSION(:,:), ALLOCATABLE            :: coord
        REAL(KIND=DBL)      :: energy
        INTEGER             :: timestep

        CHARACTER(20)       :: session_file='session_in.dat'
        CHARACTER(30)       :: ref_struc_file
        CHARACTER(30)       :: target_struc_file
        CHARACTER(30)       :: output_file

        INTEGER             :: job_type  ! 1: calc simil, 2: find uniq
        INTEGER             :: max_uniq_struc
!=========================================================================================

! Session and initialization subroutines.
        PUBLIC  :: read_session
        PUBLIC  :: read_coord
        PUBLIC  :: print_coord

        CONTAINS

        SUBROUTINE read_session
! Read session input file.
        IMPLICIT NONE
        CHARACTER       :: dummy
        OPEN(20,file=session_file,status='old')
        READ(20,*) dummy, ref_struc_file
        READ(20,*) dummy, target_struc_file
        READ(20,*) dummy, output_file
        READ(20,*) dummy, job_type
        READ(20,*) dummy, max_uniq_struc
        CLOSE(20)
        END SUBROUTINE read_session
        

        SUBROUTINE read_coord(struc_file)
! Read xyz coordinates file.
        IMPLICIT NONE
        CHARACTER(20),INTENT(IN)   :: struc_file
        INTEGER         :: iter
        CHARACTER       :: dummy
        OPEN(21,file=struc_file,status='old')
        READ(21,*) atoms
        READ(21,*) energy, dummy, dummy, timestep
        IF (ALLOCATED(coord)) DEALLOCATE(coord)
        ALLOCATE(coord(atoms,3))
        DO iter=1,atoms
          READ(21,*) material, coord(iter,1),coord(iter,2),coord(iter,3)
        END DO
        CLOSE(21)
        END SUBROUTINE read_coord
        
        SUBROUTINE print_coord
! For debugging: Print xyz coordinate files in columns.
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


        

        END MODULE initialise
