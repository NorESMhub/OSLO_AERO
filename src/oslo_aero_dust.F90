module oslo_aero_dust

   ! Calculate emission of all dusts.
   ! Note that the mobilization is calculated in the land model and
   ! the soil erodibility factor is applied here.

   use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
   use shr_const_mod,    only: pi => shr_const_pi
   use spmd_utils,       only: mpicom, mstrid=>masterprocid, masterproc
   use spmd_utils,       only: mpi_real8, mpi_character, mpi_success
   use namelist_utils,   only: find_group_name
   use ppgrid,           only: pcols, begchunk, endchunk
   use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
   use constituents,     only: cnst_name, pcnst
   use interpolate_data, only: lininterp_init, lininterp
   use interpolate_data, only: lininterp_finish, interp_type
   use cam_logfile,      only: iulog
   use cam_abortutils,   only: endrun
   use cam_pio_utils,    only: cam_pio_openfile
   use ioFileMod,        only: getfil
   use pio,              only: file_desc_t, pio_inq_dimid, pio_inq_dimlen
   use pio,              only: pio_get_var, pio_inq_varid, PIO_NOWRITE
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
   integer             :: tracerMap(numberOfDustModes) = (/-99, -99/) !index of dust tracers in the modes

   integer , parameter, public :: dust_nbin = numberOfDustModes
   real(r8), parameter :: unset_r8 = huge(1.0_r8)
   !Related to soil erodibility
   real(r8)  :: dust_emis_fact = unset_r8        ! tuning parameter for dust emissions
   character(len=cl)   :: soil_erod_file = 'none' ! full pathname for soil erodibility dataset

   real(r8), allocatable ::  soil_erodibility(:,:) ! soil erodibility factor

   real(r8), parameter, public ::  d2r  = pi/180._r8                 ! radians to degrees

   real(r8) :: emis_fact_in_coarse_mode = unset_r8  ! tuning parameter for distribution of dust emissions between modes
   real(r8), public    :: emis_fraction_in_mode(numberOfDustModes) 

!=============================================================================
contains
!=============================================================================

   subroutine oslo_aero_dust_readnl(nlfile)
      use shr_dust_emis_mod, only: is_dust_emis_zender
      use shr_dust_emis_mod, only: is_zender_soil_erod_from_atm
      use shr_dust_emis_mod, only: shr_dust_emis_readnl

      character(len=*), intent(in) :: nlfile  ! filepath for namelist file

      ! Local variables
      integer                     :: unitn, ierr
      character(len=*), parameter :: subname = 'dust_readnl'

      namelist /dust_nl/ dust_emis_fact, soil_erod_file, emis_fact_in_coarse_mode  
      !---------------------------------------------------------------------------

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

      ! Broadcast namelist variables
      call mpi_bcast(dust_emis_fact, 1, mpi_real8, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) then
         call endrun(subname//" mpi_bcast: dust_emis_fact")
      end if
      call mpi_bcast(soil_erod_file, len(soil_erod_file), mpi_character, &
           mstrid, mpicom, ierr)
      if (ierr /= mpi_success) then
         call endrun(subname//" mpi_bcast: soil_erod_file")
      end if
      call mpi_bcast(emis_fact_in_coarse_mode, 1, mpi_real8, mstrid, &
             mpicom, ierr)
      if (ierr /= mpi_success) then
         call endrun(subname//" mpi_bcast: emis_fact_in_coarse_mode")  
      end if
      
      if (emis_fact_in_coarse_mode == unset_r8) then 
         call endrun(subname//" emis_fact_in_coarse_mode not set")
      endif

      if (dust_emis_fact == unset_r8) then
         call endrun(subname//" dust_emis_fact not set")
      end if

      if (emis_fact_in_coarse_mode < 0.0_r8 .or. emis_fact_in_coarse_mode > 1.0_r8) then
         call endrun(subname//" emis_fact_in_coarse_mode should be between 0 and 1")
      end if
      
      call shr_dust_emis_readnl(mpicom, 'drv_flds_in')

      if ((soil_erod_file /= 'none') .and. &
           (.not. is_zender_soil_erod_from_atm())) then
         call endrun(subname//': should not specify soil_erod_file if Zender soil erosion is not in CAM')
      end if
      
      emis_fraction_in_mode = (/1.0_r8-emis_fact_in_coarse_mode, emis_fact_in_coarse_mode /)
      
      if (masterproc) then
         if (is_dust_emis_zender()) then
            write(iulog,*) subname,': Zender_2003 dust emission method is being used.'
         end if
         if (is_zender_soil_erod_from_atm()) then
            write(iulog, *) subname, ': Zender soil erod file is handled in atm'
            write(iulog, *) subname, ': soil_erod_file = ', trim(soil_erod_file)
            write(iulog, *) subname, ': dust_emis_fact = ', dust_emis_fact
         else
            write(iulog,*) subname,': Leung_2023 dust emission method is being used.'
            write(iulog,*) subname,': dust_emis_fact = ', dust_emis_fact
         end if
            write(iulog,*) 'emis_fact_in_coarse_mode = ', emis_fact_in_coarse_mode
      end if

   end subroutine oslo_aero_dust_readnl

   !=============================================================================
   subroutine oslo_aero_dust_init()
      use shr_dust_emis_mod, only: is_zender_soil_erod_from_atm

      ! local variables
      integer :: imode

      if (is_zender_soil_erod_from_atm()) then
         call soil_erod_init()
      end if

      ! Set module variables
      tracerMap(1) = l_dst_a2
      tracerMap(2) = l_dst_a3

      dust_names(:) = "      "
      do imode = 1, numberOfDustModes
         dust_names(imode) = cnst_name(tracerMap(imode))
      end do

   end subroutine oslo_aero_dust_init

   !=============================================================================
   subroutine oslo_aero_dust_emis(lchnk, ncol, dstflx, cflx)

      !-----------------------------------------------------------------------
      ! Purpose: Interface to emission of all dusts.
      ! Notice that the mobilization is calculated in the land model and
      ! the soil erodibility factor is applied here.
      !-----------------------------------------------------------------------
      use shr_dust_emis_mod, only: is_zender_soil_erod_from_atm

      ! Arguments:
      integer,  intent(in)    :: lchnk
      integer,  intent(in)    :: ncol
      real(r8), intent(in)    :: dstflx(pcols,4)
      real(r8), intent(inout) :: cflx(pcols,pcnst) ! Surface fluxes

      ! Local variables
      integer  :: icol,imode
      real(r8) :: soil_erod_tmp(pcols)
      real(r8) :: totalEmissionFlux(pcols)
      ! Note that following CESM use of "dust_emis_fact", the emissions are
      ! scaled by the INVERSE of the factor!!
      if (is_zender_soil_erod_from_atm()) then
         ! Filter away unreasonable values for soil erodibility
         ! (using low values e.g. gives emissions in greenland..)
         where(soil_erodibility(:,lchnk) < 0.1_r8)
            soil_erod_tmp(:)=0.0_r8
         elsewhere
            soil_erod_tmp(:)=soil_erodibility(:,lchnk)
         end where
  
         totalEmissionFlux(:) = 0.0_r8
         do icol=1,ncol
            totalEmissionFlux(icol) = totalEmissionFlux(icol) + sum(dstflx(icol,:))
         end do
         ! The flux calculations from the Zender_2003 parameterisation also include a 
         ! second scaling factor of 1.15
         do imode = 1,numberOfDustModes
            cflx(:ncol, tracerMap(imode)) = -1.0_r8*emis_fraction_in_mode(imode) &
                 *totalEmissionFlux(:ncol)*soil_erod_tmp(:ncol)/(dust_emis_fact)*1.15_r8
         end do
  
      else ! Leung emissions
  
         totalEmissionFlux(:) = 0.0_r8
         do icol=1,ncol
            totalEmissionFlux(icol) = totalEmissionFlux(icol) + sum(dstflx(icol,:))
         end do
  
         do imode = 1,numberOfDustModes
            cflx(:ncol, tracerMap(imode)) = -1.0_r8*emis_fraction_in_mode(imode) &
                 *totalEmissionFlux(:ncol)/(dust_emis_fact)
         end do
      end if

   end subroutine oslo_aero_dust_emis

   !=============================================================================
   subroutine soil_erod_init()

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

      ! Summary to log file
      if (masterproc) then
         write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
         write(iulog,*) 'soil_erod_mod: dust_emis_fact = ', dust_emis_fact
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
