module oslo_aero_dust

  ! Calculate emission of all dusts.
  ! Note that the mobilization is calculated in the land model and
  ! the soil erodibility factor is applied here.

  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,           only: pcols, begchunk, endchunk
  use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t
  use spmd_utils,       only: masterproc
  use constituents,     only: cnst_name
  use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
  use mo_constants,     only: pi, d2r
  use cam_logfile,      only: iulog
  use cam_abortutils,   only: endrun
  use cam_pio_utils,    only: cam_pio_openfile
  use ioFileMod,        only: getfil
  use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
  !
  use oslo_aero_share,  only: l_dst_a2, l_dst_a3

  implicit none
  private

  ! public routines
  public :: oslo_aero_dust_readnl
  public :: oslo_aero_dust_init
  public :: oslo_aero_dust_emis

  ! private routines (previously in soil_erod_mod in CAM)
  private :: soil_erod_init

  character(len=6), public :: dust_names(10)

  integer , parameter :: numberOfDustModes = 2  !define in oslo_aero_share?
  real(r8), parameter :: emis_fraction_in_mode(numberOfDustModes) = (/0.13_r8, 0.87_r8 /)
  integer             :: tracerMap(numberOfDustModes) = (/-99, -99/) !index of dust tracers in the modes

  integer , parameter, public :: dust_nbin = numberOfDustModes

  !Related to soil erodibility
  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical, parameter, public :: dust_active = .TRUE.

  real(r8), allocatable ::  soil_erodibility(:,:)  ! soil erodibility factor
  real(r8) :: soil_erod_fact                       ! tuning parameter for dust emissions

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if
#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif
  end subroutine oslo_aero_dust_readnl

  !===============================================================================
  subroutine oslo_aero_dust_init()

    ! local variables
    integer :: i

    call soil_erod_init( dust_emis_fact, soil_erod_file )

    ! Set module variables
    tracerMap(1) = l_dst_a2
    tracerMap(2) = l_dst_a3

    dust_names(:)="      "
    do i=1,numberOfDustModes
       dust_names(i) = cnst_name(tracerMap(i))
    end do

  end subroutine oslo_aero_dust_init

  !===============================================================================
  subroutine oslo_aero_dust_emis(state, cam_in)

    !----------------------------------------------------------------------- 
    ! Purpose: Interface to emission of all dusts.
    ! Notice that the mobilization is calculated in the land model and
    ! the soil erodibility factor is applied here.
    !-----------------------------------------------------------------------

    ! Arguments:
    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

    ! Local variables
    integer           :: lchnk
    integer           :: ncol
    integer           :: i,n
    real(r8)          :: soil_erod_tmp(pcols)
    real(r8)          :: totalEmissionFlux(pcols)
    real(r8), pointer :: cflx(:,:)

    lchnk = state%lchnk
    ncol = state%ncol

    ! Filter away unreasonable values for soil erodibility
    ! (using low values e.g. gives emissions in greenland..)
    where(soil_erodibility(:,lchnk) .lt. 0.1_r8)
       soil_erod_tmp(:)=0.0_r8
    elsewhere
       soil_erod_tmp(:)=soil_erodibility(:,lchnk)
    end where

    totalEmissionFlux(:) = 0.0_r8
    do i=1,ncol
       totalEmissionFlux(i) = totalEmissionFlux(i) + sum(cam_in%dstflx(i,:))
    end do

    ! Note that following CESM use of "dust_emis_fact", the emissions are 
    ! scaled by the INVERSE of the factor!!
    ! There is another random scale factor of 1.15 there. Adapting the exact
    ! same formulation as MAM now and tune later
    ! As of NE-380: Oslo dust emissions are 2/3 of CAM emissions
    ! gives better AOD close to dust sources

    cflx => cam_in%cflx
    do n = 1,numberOfDustModes
       cflx(:ncol, tracerMap(n)) = -1.0_r8*emis_fraction_in_mode(n) &
            *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8  
    end do

  end subroutine oslo_aero_dust_emis

  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )

    ! arguments
    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file

    ! localvaraibles
    real(r8), allocatable :: soil_erodibility_in(:,:)
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
    character(len=cl)     :: infile
    integer               :: did, vid, nlat, nlon
    type(file_desc_t)     :: ncid
    type(interp_type)     :: lon_wgts, lat_wgts
    real(r8)              :: to_lats(pcols), to_lons(pcols)
    integer               :: c, ncols, ierr
    real(r8), parameter   :: zero=0._r8
    real(r8), parameter   :: twopi=2._r8*pi

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
    end if

    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.  
    call getfil(soil_erod_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(soil_erodibility_in(nlon,nlat))

    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
    ierr = pio_get_var( ncid, vid, soil_erodibility_in )

    ! convert to radians and setup regridding
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    soil_erodibility(:,:)=0._r8

    ! regrid
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon, nlat, soil_erodibility(:,c), ncols, lon_wgts, lat_wgts)

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to deallocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to deallocate soil_erodibility_in')
    end if

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

end module oslo_aero_dust
