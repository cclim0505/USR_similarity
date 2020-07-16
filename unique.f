        MODULE uniq
        USE constants       ,ONLY: DBL
        USE initialise      ,ONLY: atoms, material, energy, &
     &    read_coord, &
     &    coord, &
     &    timestep, max_uniq_struc, &
     &    target_struc_file
        USE simil           ,ONLY: compute_struc_moments

        CHARACTER(30)                   :: uniq_file='output_uniq.xyz'
        REAL(KIND=DBL),DIMENSION(:,:,:), ALLOCATABLE   :: uniq_coord
        INTEGER                         :: uniq_count
        INTEGER                         :: max_uniq 
        REAL                            :: uniq_criteria = 0.99000E0

        ! n x 16 array of mu values for selected reference structures
        REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE  :: ref_mu_list     

        PUBLIC  :: init_uniq_var
        PUBLIC  :: save_uniq_struc
        PUBLIC  :: get_first_uniq

        CONTAINS

        SUBROUTINE init_uniq_var
! Subroutine to initialise unique module variables.
        max_uniq = max_uniq_struc
        IF (ALLOCATED(ref_mu_list)) DEALLOCATE(ref_mu_list)
        ALLOCATE(ref_mu_list(max_uniq,16))
        IF (ALLOCATED(uniq_coord)) DEALLOCATE(uniq_coord)
        ALLOCATE(uniq_coord(max_uniq,atoms,3))
        uniq_count = 1

        END SUBROUTINE init_uniq_var

        SUBROUTINE save_uniq_struc(in_coord)
! Save structure as uniq structure into output file.
        IMPLICIT NONE
        REAL(KIND=DBL),DIMENSION(:,:), INTENT(IN)   :: in_coord
        INTEGER         :: f_save
        INTEGER         :: iter

        OPEN(NEWUNIT=f_save,FILE=uniq_file,ACCESS='append')
        WRITE(f_save,*) atoms
        WRITE(f_save,*) energy, timestep
        DO iter=1,atoms
          WRITE(f_save,*) material, in_coord(iter,1), &
     &      in_coord(iter,2), in_coord(iter,3)
        END DO
        CLOSE(f_save)
        END SUBROUTINE save_uniq_struc

        SUBROUTINE get_first_uniq
! Read first structure in target set and set as first unique structure.
        IMPLICIT NONE
        CALL read_coord(target_struc_file)
        CALL save_uniq_struc(coord)
        CALL compute_struc_moments(ref_mu_list(1,:))
        END SUBROUTINE get_first_uniq


        END MODULE uniq
