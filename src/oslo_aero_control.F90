module oslo_aero_control

  !-----------------------------------------------------------------------
  ! Provides a control interface to CAM-Oslo packages
  !-----------------------------------------------------------------------

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use spmd_utils,        only: mpicom, mstrid=>masterprocid, masterproc
  use spmd_utils,        only: mpi_logical, mpi_real8, mpi_character, mpi_integer,  mpi_success
  use namelist_utils,    only: find_group_name
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun
  use atm_import_export, only: dms_from_ocn

  implicit none
  private

  public :: oslo_aero_ctl_readnl ! read namelist from file
  public :: oslo_aero_getopts    ! generic query method

  ! Private module data
  character(len=16), parameter :: unset_str = 'UNSET'
  integer,           parameter :: unset_int = huge(1)
  integer, parameter, public   :: dir_string_length=256

  ! Namelist variables:
  real(r8)  :: volc_fraction_coarse = 0.0_r8  !Fraction of volcanic aerosols in coarse mode
  character(len=dir_string_length) :: aerotab_table_dir = unset_str

  ! DMS/Ocean namelist variables
  character(len=20)                :: dms_source       = unset_str
  character(len=32)                :: dms_source_type  = unset_str
  character(len=20)                :: opom_source      = unset_str
  character(len=32)                :: opom_source_type = unset_str
  character(len=dir_string_length) :: ocean_filename   = unset_str
  character(len=dir_string_length) :: ocean_filepath   = unset_str
  integer                          :: dms_cycle_year   = 0 ! =unset_int?
  integer                          :: opom_cycle_year  = 0 ! =unset_int?

!=======================================================================
contains
!=======================================================================

  subroutine oslo_aero_ctl_readnl(nlfile)

    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local variables
    integer :: unitn, ierr
    logical :: fileExists=.false.
    character(len=*), parameter :: subname = 'oslo_ctl_readnl'

    namelist /oslo_ctl_nl/ volc_fraction_coarse, aerotab_table_dir, dms_source, &
                           dms_source_type, opom_source, opom_source_type, &
                           ocean_filename, ocean_filepath, dms_cycle_year, opom_cycle_year
    !-----------------------------------------------------------------------------

    if (masterproc) then
       open (newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'oslo_ctl_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, oslo_ctl_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(volc_fraction_coarse, 1 , mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: volc_fraction_coarse")
    call mpi_bcast(aerotab_table_dir, len(aerotab_table_dir) , mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: aerotab_table_dir")

    ! dms variables
    call mpi_bcast(dms_source, len(dms_source), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: dms_source")
    call mpi_bcast(dms_source_type, len(dms_source_type), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: dms_source_type")
    call mpi_bcast(dms_cycle_year, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: dms_cycle_year")

    ! opom variables
    call mpi_bcast(opom_source, len(opom_source), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: opom_source")
    call mpi_bcast(opom_source_type, len(opom_source_type), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: opom_source_type")
    call mpi_bcast(opom_cycle_year, 1, mpi_integer, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: opom_cycle_year")

    ! ocean variables
    call mpi_bcast(ocean_filename, len(ocean_filename), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: ocean_filename")
    call mpi_bcast(ocean_filepath, len(ocean_filepath), mpi_character, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: ocean_filepath")

    ! Reset dms_source if ocean is sending dms to atm
    if (dms_from_ocn) then
       dms_source = 'ocean_flux'
    end if

    ! Error checking:

    ! Defaults for PBL and microphysics are set in build-namelist.  Check here that
    ! values have been set to guard against problems with hand edited namelists.
    if(volc_fraction_coarse < 0.0_r8 .OR. volc_fraction_coarse > 1.0_r8)then
       write(iulog,*)'cam_oslo: illegal value of volc_fraction_coarse', volc_fraction_coarse
       call endrun('cam_oslo: illegal value of volc_fraction_coarse')
    end if

    if (masterproc) then
       write(iulog,*)"Reading aerosol tables from : " // trim(aerotab_table_dir)
    endif

    ! Error check for OCEAN file
    inquire( file=trim(ocean_filepath)//'/'//trim(ocean_filename), exist=fileExists )
    if(.not. fileExists)then
       call endrun("oslo_aero_control: can not find ocean file "//trim(ocean_filepath)//'/'//trim(ocean_filename))
    else
       if (masterproc) then
          write(iulog,*)"Reading ocean tracers from : " // trim(ocean_filepath)//'/'//trim(ocean_filename)
       end if
    endif

    ! Error check for dms_source from namelist
    if ( dms_source =='ocean_flux' .or. &
         dms_source =='kettle'     .or. &
         dms_source =='lana'       .or. &
         dms_source =='emission_file') then
       if (masterproc) then
          write(iulog,*)"DMS emission source is : "// trim(dms_source)
       end if
    else
       call endrun("oslo_aero_control: no valid dms source from namelist: " //trim(dms_source))
    endif

    ! Error check for opom_source from namelist
    if ( opom_source=='no_file' .or. &
         opom_source=='nilsson' .or. &
         opom_source=='odowd') then
       if (masterproc) then
          write(iulog,*)"Ocean POM emission source is : "// trim(opom_source)
       end if
    else
       call endrun("oslo_aero_control: no valid opom source from namelist: " //trim(opom_source))
    endif

  end subroutine oslo_aero_ctl_readnl

  !==========================================================================
  subroutine oslo_aero_getopts(  &
       volc_fraction_coarse_out, &
       aerotab_table_dir_out,    &
       dms_source_out,           &
       dms_source_type_out,      &
       opom_source_out,          &
       opom_source_type_out,     &
       ocean_filename_out,       &
       ocean_filepath_out,       &
       opom_cycle_year_out,      &
       dms_cycle_year_out  )

    !-----------------------------------------------------------------------
    ! Purpose: Return runtime settings
    !-----------------------------------------------------------------------

    real(r8)         , intent(out), optional :: volc_fraction_coarse_out
    character(len=*) , intent(out), optional :: aerotab_table_dir_out
    character(len=*) , intent(out), optional :: ocean_filename_out
    character(len=*) , intent(out), optional :: ocean_filepath_out
    character(len=*) , intent(out), optional :: dms_source_out
    character(len=*) , intent(out), optional :: dms_source_type_out
    integer          , intent(out), optional :: dms_cycle_year_out
    character(len=*) , intent(out), optional :: opom_source_out
    character(len=*) , intent(out), optional :: opom_source_type_out
    integer          , intent(out), optional :: opom_cycle_year_out

    if ( present(volc_fraction_coarse_out ) ) volc_fraction_coarse_out = volc_fraction_coarse
    if ( present(aerotab_table_dir_out    ) ) aerotab_table_dir_out = aerotab_table_dir
    if ( present(ocean_filename_out       ) ) ocean_filename_out  = ocean_filename
    if ( present(ocean_filepath_out       ) ) ocean_filepath_out  = ocean_filepath
    if ( present(dms_source_out           ) ) dms_source_out      = dms_source
    if ( present(dms_source_type_out      ) ) dms_source_type_out = dms_source_type
    if ( present(dms_cycle_year_out       ) ) dms_cycle_year_out  = dms_cycle_year
    if ( present(opom_source_out          ) ) opom_source_out     = opom_source
    if ( present(opom_source_type_out     ) ) opom_source_type_out= opom_source_type
    if ( present(opom_cycle_year_out      ) ) opom_cycle_year_out = opom_cycle_year

  end subroutine oslo_aero_getopts

end module oslo_aero_control
