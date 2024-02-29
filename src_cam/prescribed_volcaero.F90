module prescribed_volcaero

   use shr_kind_mod,     only : r8 => shr_kind_r8, cs => shr_kind_cs
   use cam_abortutils,   only : endrun
   use spmd_utils,       only : mpicom, mstrid=>masterprocid, masterproc
   use spmd_utils,       only : mpi_logical, mpi_real8, mpi_character, mpi_integer, mpi_success
   use tracer_data,      only : trfld, trfile
   use cam_logfile,      only : iulog

   implicit none
   private

   ! Methods
   public :: prescribed_volcaero_readnl
   public :: prescribed_volcaero_register
   public :: prescribed_volcaero_init
   public :: prescribed_volcaero_adv
   public :: init_prescribed_volcaero_restart
   public :: write_prescribed_volcaero_restart
   public :: read_prescribed_volcaero_restart

   ! Module variables
   public :: has_prescribed_volcaero
   public :: solar_bands
   public :: terrestrial_bands

   type(trfld), pointer :: fields(:)
   type(trfile)         :: file

   ! These variables are settable via the namelist (with longer names)
   character(len=cs)  :: filename = ''
   character(len=cs)  :: filelist = ''
   character(len=cs)  :: datapath = ''
   character(len=cs)  :: data_type = 'SERIAL'
   logical            :: rmv_file = .false.
   integer            :: cycle_yr  = 0
   integer            :: fixed_ymd = 0
   integer            :: fixed_tod = 0
   logical            :: has_prescribed_volcaero = .true.
   integer, parameter :: solar_bands=14
   integer, parameter :: terrestrial_bands=16

!===============================================================================
contains
!===============================================================================

   subroutine prescribed_volcaero_readnl(nlfile)

      use namelist_utils,  only : find_group_name

      ! Arguments
      character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

      ! Local variables
      integer            :: unitn, ierr
      character(len=16)  :: prescribed_volcaero_name ! not used
      character(len=256) :: prescribed_volcaero_file
      character(len=256) :: prescribed_volcaero_filelist
      character(len=256) :: prescribed_volcaero_datapath
      character(len=32)  :: prescribed_volcaero_type
      logical            :: prescribed_volcaero_rmfile
      integer            :: prescribed_volcaero_cycle_yr
      integer            :: prescribed_volcaero_fixed_ymd
      integer            :: prescribed_volcaero_fixed_tod
      character(len=*), parameter :: subname = 'prescribed_volcaero_readnl'

      namelist /prescribed_volcaero_nl/   &
           prescribed_volcaero_name,      &
           prescribed_volcaero_file,      &
           prescribed_volcaero_filelist,  &
           prescribed_volcaero_datapath,  &
           prescribed_volcaero_type,      &
           prescribed_volcaero_rmfile,    &
           prescribed_volcaero_cycle_yr,  &
           prescribed_volcaero_fixed_ymd, &
           prescribed_volcaero_fixed_tod
      !---------------------------------------------------

      ! Initialize namelist variables from local module variables.
      prescribed_volcaero_file     = filename
      prescribed_volcaero_filelist = filelist
      prescribed_volcaero_datapath = datapath
      prescribed_volcaero_type     = data_type
      prescribed_volcaero_rmfile   = rmv_file
      prescribed_volcaero_cycle_yr = cycle_yr
      prescribed_volcaero_fixed_ymd= fixed_ymd
      prescribed_volcaero_fixed_tod= fixed_tod

      ! Read namelist
      if (masterproc) then
         open( newunit=unitn, file=trim(nlfile), status='old' )
         call find_group_name(unitn, 'prescribed_volcaero_nl', status=ierr)
         if (ierr == 0) then
            read(unitn, prescribed_volcaero_nl, iostat=ierr)
            if (ierr /= 0) then
               call endrun(subname // ':: ERROR reading namelist')
            end if
         end if
         close(unitn)
      end if

      call mpi_bcast(prescribed_volcaero_file,     len(prescribed_volcaero_file), mpi_character, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_file")
      call mpi_bcast(prescribed_volcaero_filelist, len(prescribed_volcaero_filelist), mpi_character, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_filelist")
      call mpi_bcast(prescribed_volcaero_datapath, len(prescribed_volcaero_datapath), mpi_character, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_datapath")
      call mpi_bcast(prescribed_volcaero_type,     len(prescribed_volcaero_type), mpi_character, mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_type")
      call mpi_bcast(prescribed_volcaero_rmfile,   1, mpi_logical,  mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_rmfile")
      call mpi_bcast(prescribed_volcaero_cycle_yr, 1, mpi_integer,  mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_cycle_yr")
      call mpi_bcast(prescribed_volcaero_fixed_ymd,1, mpi_integer,  mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_fixed_ymd")
      call mpi_bcast(prescribed_volcaero_fixed_tod,1, mpi_integer,  mstrid, mpicom, ierr)
      if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: prescribed_volcaero_fixed_tod")

      ! Update module variables with user settings.
      filename   = prescribed_volcaero_file
      filelist   = prescribed_volcaero_filelist
      datapath   = prescribed_volcaero_datapath
      data_type  = prescribed_volcaero_type
      rmv_file   = prescribed_volcaero_rmfile
      cycle_yr   = prescribed_volcaero_cycle_yr
      fixed_ymd  = prescribed_volcaero_fixed_ymd
      fixed_tod  = prescribed_volcaero_fixed_tod

      ! Turn on prescribed volcanics if user has specified an input dataset.
      if (len_trim(filename) > 0 .and. filename.ne.'NONE') then
         has_prescribed_volcaero = .true.
      end if

   end subroutine prescribed_volcaero_readnl

   !===============================================================================
   subroutine prescribed_volcaero_register()

      use ppgrid,         only: pver,pcols
      use physics_buffer, only: pbuf_add_field, dtype_r8

      ! Local variables
      integer :: idx
      integer :: band
      character(len=3) :: c3
      !---------------------------------------------------

      if (has_prescribed_volcaero) then
         do band=1,solar_bands
            write(c3,'(i3)') band
            call pbuf_add_field('ext_sun'   //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
            call pbuf_add_field('omega_sun' //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
            call pbuf_add_field('g_sun'     //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
         enddo
         do band=1,terrestrial_bands
            write(c3,'(i3)') band
            call pbuf_add_field('ext_earth'   //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
            call pbuf_add_field('omega_earth' //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
            call pbuf_add_field('g_earth'     //trim(adjustl(c3)),'physpkg',dtype_r8,(/pcols,pver/),idx)
         enddo
      endif

   endsubroutine prescribed_volcaero_register

   !===============================================================================
   subroutine prescribed_volcaero_init()

      use tracer_data, only : trcdata_init
      use cam_history, only : addfld

      ! Local variables
      integer           :: band
      character(len=3)  :: c3
      character(len=32) :: specifier(3*(solar_bands+terrestrial_bands))
      !---------------------------------------------------

      if ( has_prescribed_volcaero) then
         if ( masterproc ) then
            write(iulog,*) 'volcanic aerosol is prescribed in :'//trim(filename)
         endif

         do band=1,solar_bands
            write(c3,'(i3)') band
            specifier(band*3-2) = 'ext_sun'   //trim(adjustl(c3))//':'//'ext_sun'//trim(adjustl(c3))
            specifier(band*3-1) = 'omega_sun' //trim(adjustl(c3))//':'//'omega_sun'//trim(adjustl(c3))
            specifier(band*3-0) = 'g_sun'     //trim(adjustl(c3))//':'//'g_sun'//trim(adjustl(c3))
            call addfld('ext_sun'   //trim(adjustl(c3)),(/ 'lev' /), 'I', '1/km', 'Extinction coefficient of solar bands' )
            call addfld('omega_sun' //trim(adjustl(c3)),(/ 'lev' /), 'I', '1'   , 'Single scattering albedo of solar bands' )
            call addfld('g_sun'     //trim(adjustl(c3)),(/ 'lev' /), 'I', '1'   , 'Asymmetry factor of solar bands' )
         enddo
         do band=1,terrestrial_bands
            write(c3,'(i3)') band
            specifier((solar_bands+band)*3-2) = 'ext_earth'   //trim(adjustl(c3))//':'//'ext_earth'//trim(adjustl(c3))
            specifier((solar_bands+band)*3-1) = 'omega_earth' //trim(adjustl(c3))//':'//'omega_earth'//trim(adjustl(c3))
            specifier((solar_bands+band)*3-0) = 'g_earth'     //trim(adjustl(c3))//':'//'g_earth'//trim(adjustl(c3))
            call addfld('ext_earth'   //trim(adjustl(c3)),(/ 'lev' /), 'I', '1/km', 'Extinction coefficient of terrestrial bands' )
            call addfld('omega_earth' //trim(adjustl(c3)),(/ 'lev' /), 'I', '1'   , 'Single scattering albedo of terrestrial bands' )
            call addfld('g_earth'     //trim(adjustl(c3)),(/ 'lev' /), 'I', '1'   , 'Asymmetry factor of terrestrial bands' )
         enddo

         allocate(file%in_pbuf(size(specifier)))
         file%in_pbuf(:) = .true.
         call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
              rmv_file, cycle_yr, fixed_ymd, fixed_tod, data_type)
      endif

   end subroutine prescribed_volcaero_init

   !===============================================================================
   subroutine prescribed_volcaero_adv( state, pbuf2d)

      use tracer_data,    only : advance_trcdata
      use physics_types,  only : physics_state
      use physics_buffer, only : physics_buffer_desc, pbuf_get_field, pbuf_get_chunk
      use ppgrid,         only : begchunk, endchunk, pcols, pver
      use cam_history,    only : outfld
      use tropopause,     only : tropopause_find, TROP_ALG_TWMO, TROP_ALG_CLIMATE

      ! Arguments
      type(physics_state), intent(in)    :: state(begchunk:endchunk)
      type(physics_buffer_desc), pointer :: pbuf2d(:,:)

      ! Local variables
      integer           :: c,ncol,i,k
      integer           :: band
      real(r8), pointer :: data(:,:)
      character(len=3)  :: c3
      integer           :: tropLev(pcols)
      type(physics_buffer_desc), pointer :: pbuf_chnk(:)
      !---------------------------------------------------

      if ( has_prescribed_volcaero) then
         call advance_trcdata( fields, file, state, pbuf2d )

         do c = begchunk,endchunk
            pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
            call tropopause_find(state(c), tropLev, primary=TROP_ALG_TWMO, backup=TROP_ALG_CLIMATE)
            ncol = state(c)%ncol
            do band=1,solar_bands
               write(c3,'(i3)') band
               call pbuf_get_field(pbuf_chnk, fields(band*3-2)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0._r8
                  enddo
               enddo
               call outfld('ext_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
               call pbuf_get_field(pbuf_chnk, fields(band*3-1)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0.999_r8
                  enddo
               enddo
               call outfld('omega_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
               call pbuf_get_field(pbuf_chnk, fields(band*3-0)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0.5_r8
                  enddo
               enddo
               call outfld('g_sun'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
            enddo
            do band=1,terrestrial_bands
               write(c3,'(i3)') band
               call pbuf_get_field(pbuf_chnk, fields((solar_bands+band)*3-2)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0._r8
                  enddo
               enddo
               call outfld('ext_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
               call pbuf_get_field(pbuf_chnk, fields((solar_bands+band)*3-1)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0.999_r8
                  enddo
               enddo
               call outfld('omega_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
               call pbuf_get_field(pbuf_chnk, fields((solar_bands+band)*3-0)%pbuf_ndx, data)
               do i = 1,ncol
                  do k = 1,pver
                     if ( k >= tropLev(i) ) data(i,k) = 0.5_r8
                  enddo
               enddo
               call outfld('g_earth'//trim(adjustl(c3)),data(:,:), pcols, state(c)%lchnk)
            enddo
         enddo
      endif

   end subroutine prescribed_volcaero_adv

   !===============================================================================
   subroutine init_prescribed_volcaero_restart( piofile )
      use tracer_data, only : init_trc_restart
      use pio, only : file_desc_t
      type(file_desc_t),intent(inout) :: pioFile     ! pio File pointer
      call init_trc_restart( 'prescribed_volcaero', piofile, file )
   end subroutine init_prescribed_volcaero_restart

   !===============================================================================
   subroutine write_prescribed_volcaero_restart( piofile )
      use tracer_data, only : write_trc_restart
      use pio, only : file_desc_t
      type(file_desc_t) :: piofile
      call write_trc_restart( piofile, file )
   end subroutine write_prescribed_volcaero_restart

   !===============================================================================
   subroutine read_prescribed_volcaero_restart( pioFile )
      use tracer_data, only : read_trc_restart
      use pio, only : file_desc_t
      type(file_desc_t) :: piofile
      call read_trc_restart( 'prescribed_volcaero', piofile, file )
   end subroutine read_prescribed_volcaero_restart

end module prescribed_volcaero
