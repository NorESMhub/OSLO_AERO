module mo_srf_emissions
  !---------------------------------------------------------------
  ! 	... surface emissions module
  !---------------------------------------------------------------

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use chem_mods,     only : gas_pcnst
  use spmd_utils,    only : masterproc
  use cam_abortutils,only : endrun
  use ioFileMod,     only : getfil
  use cam_logfile,   only : iulog
  use tracer_data,   only : trfld,trfile
#ifdef OSLO_AERO
  use oslo_aero_ocean, only: oslo_aero_dms_inq
#endif

  implicit none

  type :: emission
     integer           :: spc_ndx
     real(r8)          :: mw
     real(r8)          :: scalefactor
     character(len=256):: filename
     character(len=16) :: species
     character(len=8)  :: units
     integer                   :: nsectors
     character(len=32),pointer :: sectors(:)
     type(trfld), pointer      :: fields(:)
     type(trfile)              :: file
  end type emission

  private

  public  :: srf_emissions_inti, set_srf_emissions, set_srf_emissions_time

  logical, public, protected :: has_emis(gas_pcnst) = .false.

  real(r8), parameter :: amufac = 1.65979e-23_r8         ! 1.e4* kg / amu
  type(emission), allocatable :: emissions(:)
  integer                     :: n_emis_files
  integer :: c10h16_ndx, isop_ndx
#ifdef OSLO_AERO
  integer :: dms_ndx
#endif

contains

  subroutine srf_emissions_inti( srf_emis_specifier, emis_type_in, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod )

    !-----------------------------------------------------------------------
    ! 	... initialize the surface emissions
    !-----------------------------------------------------------------------

    use chem_mods,        only : adv_mass
    use mo_chem_utls,     only : get_spc_ndx
    use tracer_data,      only : trcdata_init
    use cam_pio_utils,    only : cam_pio_openfile
    use pio,              only : pio_inquire, pio_nowrite, pio_closefile, pio_inq_varndims
    use pio,              only : pio_inq_varname, pio_inq_vardimid, pio_inq_dimid
    use pio,              only : file_desc_t, pio_get_att, PIO_NOERR, PIO_GLOBAL
    use pio,              only : pio_seterrorhandling, PIO_BCAST_ERROR,PIO_INTERNAL_ERROR
    use chem_surfvals,    only : flbc_list
    use string_utils,     only : GLC
    use m_MergeSorts,     only : IndexSort

    implicit none

    !-----------------------------------------------------------------------
    ! 	... dummy arguments
    !-----------------------------------------------------------------------
    character(len=*), intent(in) :: srf_emis_specifier(:)
    character(len=*), intent(in) :: emis_type_in
    integer,          intent(in) :: emis_cycle_yr
    integer,          intent(in) :: emis_fixed_ymd
    integer,          intent(in) :: emis_fixed_tod

    !-----------------------------------------------------------------------
    ! 	... local variables
    !-----------------------------------------------------------------------
    integer  :: astat
    integer  :: j, l, m, n, i, nn                     ! Indices
    character(len=16)  :: spc_name
    character(len=256) :: filename

    character(len=16)  :: emis_species(size(srf_emis_specifier))
    character(len=256) :: emis_filenam(size(srf_emis_specifier))
    integer  :: emis_indexes(size(srf_emis_specifier))
    integer  :: indx(size(srf_emis_specifier))
    real(r8) :: emis_scalefactor(size(srf_emis_specifier))

    integer :: vid, nvars, isec, num_dims_emis
    integer :: vndims
    logical, allocatable :: is_sector(:)
    type(file_desc_t) :: ncid
    character(len=32)  :: varname
    character(len=256) :: locfn
    integer :: ierr
    character(len=1), parameter :: filelist = ''
    character(len=1), parameter :: datapath = ''
    logical         , parameter :: rmv_file = .false.
    logical :: unstructured
    character(len=32) :: emis_type = ' '
    character(len=80) :: file_interp_type = ' '
    character(len=256) :: tmp_string = ' '
    character(len=32) :: xchr = ' '
    real(r8) :: xdbl
    integer :: time_dimid, ncol_dimid
    integer, allocatable :: dimids(:)

    has_emis(:) = .false.
    nn = 0
    indx(:) = 0

    count_emis: do n=1,size(srf_emis_specifier)
       if ( len_trim(srf_emis_specifier(n) ) == 0 ) then
          exit count_emis
       endif

       i = scan(srf_emis_specifier(n),'->')
       spc_name = trim(adjustl(srf_emis_specifier(n)(:i-1)))

       ! need to parse out scalefactor ...
       tmp_string = adjustl(srf_emis_specifier(n)(i+2:))
       j = scan( tmp_string, '*' )
       if (j>0) then
          xchr = tmp_string(1:j-1) ! get the multipler (left of the '*')
          read( xchr, * ) xdbl   ! convert the string to a real
          tmp_string = adjustl(tmp_string(j+1:)) ! get the filepath name (right of the '*')
       else
          xdbl = 1._r8
       endif
       filename = trim(tmp_string)

       m = get_spc_ndx(spc_name)

       if (m > 0) then
          has_emis(m) = .true.
       else
          write(iulog,*) 'srf_emis_inti: spc_name ',spc_name,' is not included in the simulation'
          call endrun('srf_emis_inti: invalid surface emission specification')
       endif

       if (any( flbc_list == spc_name )) then
          call endrun('srf_emis_inti: ERROR -- cannot specify both fixed LBC ' &
                    //'and emissions for the same species: '//trim(spc_name))
       endif

       nn = nn+1
       emis_species(nn) = spc_name
       emis_filenam(nn) = filename
       emis_indexes(nn) = m
       emis_scalefactor(nn) = xdbl

       indx(n)=n

    enddo count_emis

    n_emis_files = nn

    if (masterproc) write(iulog,*) 'srf_emis_inti: n_emis_files = ',n_emis_files

    allocate( emissions(n_emis_files), stat=astat )
    if( astat/= 0 ) then
       write(iulog,*) 'srf_emis_inti: failed to allocate emissions array; error = ',astat
       call endrun('srf_emis_inti: failed to allocate emissions array')
    end if

    !-----------------------------------------------------------------------
    ! Sort the input files so that the emissions sources are summed in the
    ! same order regardless of the order of the input files in the namelist
    !-----------------------------------------------------------------------
    if (n_emis_files > 0) then
      call IndexSort(n_emis_files, indx, emis_filenam)
    end if

    !-----------------------------------------------------------------------
    ! 	... setup the emission type array
    !-----------------------------------------------------------------------
    do m=1,n_emis_files
       emissions(m)%spc_ndx          = emis_indexes(indx(m))
       emissions(m)%units            = 'Tg/y'
       emissions(m)%species          = emis_species(indx(m))
       emissions(m)%mw               = adv_mass(emis_indexes(indx(m)))                     ! g / mole
       emissions(m)%filename         = emis_filenam(indx(m))
       emissions(m)%scalefactor      = emis_scalefactor(indx(m))
    enddo

    !-----------------------------------------------------------------------
    ! read emis files to determine number of sectors
    !-----------------------------------------------------------------------
    spc_loop: do m = 1, n_emis_files

       emissions(m)%nsectors = 0

       if (masterproc) then
          write(iulog,'(a,i3,a)') 'srf_emissions_inti m: ',m,' init file : '//trim(emissions(m)%filename)
       endif

       call getfil (emissions(m)%filename, locfn, 0)
       call cam_pio_openfile ( ncid, trim(locfn), PIO_NOWRITE)
       ierr = pio_inquire (ncid, nVariables=nvars)

       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_inq_dimid( ncid, 'ncol', ncol_dimid )
       unstructured = ierr==PIO_NOERR
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)

       allocate(is_sector(nvars))
       is_sector(:) = .false.

       if (unstructured) then
          ierr = pio_inq_dimid( ncid, 'time', time_dimid )
       end if

       do vid = 1,nvars

          ierr = pio_inq_varndims (ncid, vid, vndims)

          if (unstructured) then
             num_dims_emis = 2
          else
             num_dims_emis = 3
          endif

          if( vndims < num_dims_emis ) then
             cycle
          elseif( vndims > num_dims_emis ) then
             ierr = pio_inq_varname (ncid, vid, varname)
             write(iulog,*) 'srf_emis_inti: Skipping variable ', trim(varname),', ndims = ',vndims, &
                  ' , species=',trim(emissions(m)%species)
             cycle
          end if

          if (unstructured) then
             allocate( dimids(vndims) )
             ierr = pio_inq_vardimid( ncid, vid, dimids )
             if ( any(dimids(:)==ncol_dimid) .and. any(dimids(:)==time_dimid) ) then
                emissions(m)%nsectors = emissions(m)%nsectors+1
                is_sector(vid)=.true.
             endif
             deallocate(dimids)
          else
             emissions(m)%nsectors = emissions(m)%nsectors+1
             is_sector(vid)=.true.
          end if

       enddo

       allocate( emissions(m)%sectors(emissions(m)%nsectors), stat=astat )
       if( astat/= 0 ) then
         write(iulog,*) 'srf_emis_inti: failed to allocate emissions(m)%sectors array; error = ',astat
         call endrun
       end if

       isec = 1

       do vid = 1,nvars
          if( is_sector(vid) ) then
             ierr = pio_inq_varname(ncid, vid, emissions(m)%sectors(isec))
             isec = isec+1
          endif
       enddo
       deallocate(is_sector)

       ! Global attribute 'input_method' overrides the srf_emis_type namelist setting on
       ! a file-by-file basis.  If the emis file does not contain the 'input_method'
       ! attribute then the srf_emis_type namelist setting is used.
       call pio_seterrorhandling(ncid, PIO_BCAST_ERROR)
       ierr = pio_get_att(ncid, PIO_GLOBAL, 'input_method', file_interp_type)
       call pio_seterrorhandling(ncid, PIO_INTERNAL_ERROR)
       if ( ierr == PIO_NOERR) then
          l = GLC(file_interp_type)
          emis_type(1:l) = file_interp_type(1:l)
          emis_type(l+1:) = ' '
       else
          emis_type = trim(emis_type_in)
       endif

       call pio_closefile (ncid)

       allocate(emissions(m)%file%in_pbuf(size(emissions(m)%sectors)))
       emissions(m)%file%in_pbuf(:) = .false.

       call trcdata_init( emissions(m)%sectors, &
                          emissions(m)%filename, filelist, datapath, &
                          emissions(m)%fields,  &
                          emissions(m)%file, &
                          rmv_file, emis_cycle_yr, emis_fixed_ymd, emis_fixed_tod, trim(emis_type) )

    enddo spc_loop

    c10h16_ndx = get_spc_ndx('C10H16')
    isop_ndx = get_spc_ndx('ISOP')
#ifdef OSLO_AERO
    dms_ndx = get_spc_ndx('DMS')
#endif
  end subroutine srf_emissions_inti

  subroutine set_srf_emissions_time( pbuf2d, state )
    !-----------------------------------------------------------------------
    !       ... check serial case for time span
    !-----------------------------------------------------------------------

    use physics_types,only : physics_state
    use ppgrid,       only : begchunk, endchunk
    use tracer_data,  only : advance_trcdata
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in):: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    !-----------------------------------------------------------------------
    !       ... local variables
    !-----------------------------------------------------------------------
    integer :: m

    do m = 1,n_emis_files
       call advance_trcdata( emissions(m)%fields, emissions(m)%file, state, pbuf2d  )
    end do

  end subroutine set_srf_emissions_time

  ! adds surf flux specified in file to sflx
  subroutine set_srf_emissions( lchnk, ncol, sflx )
    !--------------------------------------------------------
    !	... form the surface fluxes for this latitude slice
    !--------------------------------------------------------

    use mo_constants, only : pi
    use time_manager, only : get_curr_calday
    use string_utils, only : to_lower, GLC
    use phys_grid,    only : get_rlat_all_p, get_rlon_all_p

    implicit none

    !--------------------------------------------------------
    !	... Dummy arguments
    !--------------------------------------------------------
    integer,  intent(in)  :: ncol                  ! columns in chunk
    integer,  intent(in)  :: lchnk                 ! chunk index
    real(r8), intent(out) :: sflx(:,:) ! surface emissions ( kg/m^2/s )

    !--------------------------------------------------------
    !	... local variables
    !--------------------------------------------------------
    integer  ::  i, m, n
    real(r8) ::  factor
    real(r8) ::  dayfrac            ! fration of day in light
    real(r8) ::  iso_off            ! time iso flux turns off
    real(r8) ::  iso_on             ! time iso flux turns on

    logical  :: polar_day,polar_night
    real(r8) :: doy_loc
    real(r8) :: sunon,sunoff
    real(r8) :: loc_angle
    real(r8) :: latitude
    real(r8) :: declination
    real(r8) :: tod
    real(r8) :: calday

    real(r8), parameter :: dayspy = 365._r8
    real(r8), parameter :: twopi = 2.0_r8 * pi
    real(r8), parameter :: pid2  = 0.5_r8 * pi
    real(r8), parameter :: dec_max = 23.45_r8 * pi/180._r8

    real(r8) :: flux(ncol)
    real(r8) :: mfactor
    integer  :: isec

    character(len=12),parameter :: mks_units(4) = (/ "kg/m2/s     ", &
                                                     "kg/m2/sec   ", &
                                                     "kg/m^2/s    ", &
                                                     "kg/m^2/sec  " /)
    character(len=12) :: units

    real(r8), dimension(ncol) :: rlats, rlons

    sflx(:,:) = 0._r8

    !--------------------------------------------------------
    !	... set non-zero emissions
    !--------------------------------------------------------
    emis_loop : do m = 1,n_emis_files

       n = emissions(m)%spc_ndx

       flux(:) = 0._r8
       do isec = 1,emissions(m)%nsectors
          flux(:ncol) = flux(:ncol) + emissions(m)%scalefactor*emissions(m)%fields(isec)%data(:ncol,1,lchnk)
       enddo

       units = to_lower(trim(emissions(m)%fields(1)%units(:GLC(emissions(m)%fields(1)%units))))

       if ( any( mks_units(:) == units ) ) then
          sflx(:ncol,n) = sflx(:ncol,n) + flux(:ncol)
       else
          mfactor = amufac * emissions(m)%mw
          sflx(:ncol,n) = sflx(:ncol,n) + flux(:ncol) * mfactor
       endif

    end do emis_loop

    call get_rlat_all_p( lchnk, ncol, rlats )
    call get_rlon_all_p( lchnk, ncol, rlons )

    calday = get_curr_calday()
    doy_loc     = aint( calday )
    declination = dec_max * cos((doy_loc - 172._r8)*twopi/dayspy)
    tod = (calday - doy_loc) + .5_r8

#ifdef OSLO_AERO
    ! Zero DMS emissions if option is not "from file"
    ! oslo_aero_dms_inq() Returns "true" if "emissions from file"
    if (.not. oslo_aero_dms_inq()) then
       if (dms_ndx > 0) then
          sflx(:ncol,dms_ndx) = 0.0_r8
       end if
    end if
#endif

    do i = 1,ncol
       !
       polar_day   = .false.
       polar_night = .false.
       !
       loc_angle = tod * twopi + rlons(i)
       loc_angle = mod( loc_angle,twopi )
       latitude =  rlats(i)
       !
       !------------------------------------------------------------------
       !        determine if in polar day or night
       !        if not in polar day or night then
       !        calculate terminator longitudes
       !------------------------------------------------------------------
       if( abs(latitude) >= (pid2 - abs(declination)) ) then
          if( sign(1._r8,declination) == sign(1._r8,latitude) ) then
             polar_day = .true.
             sunoff = 2._r8*twopi
             sunon  = -twopi
          else
             polar_night = .true.
          end if
       else
          sunoff = acos( -tan(declination)*tan(latitude) )
          sunon  = twopi - sunoff
       end if

       !--------------------------------------------------------
       !	... adjust alpha-pinene for diurnal variation
       !--------------------------------------------------------
       if( c10h16_ndx > 0 ) then
          if( has_emis(c10h16_ndx) ) then
             if( .not. polar_night .and. .not. polar_day ) then
                dayfrac = sunoff / pi
                sflx(i,c10h16_ndx) = sflx(i,c10h16_ndx) / (.7_r8 + .3_r8*dayfrac)
                if( loc_angle >= sunoff .and. loc_angle <= sunon ) then
                   sflx(i,c10h16_ndx) = sflx(i,c10h16_ndx) * .7_r8
                endif
             end if
          end if
       end if

       !--------------------------------------------------------
       !	... adjust isoprene for diurnal variation
       !--------------------------------------------------------
       if( isop_ndx > 0 ) then
          if( has_emis(isop_ndx) ) then
             if( .not. polar_night ) then
                if( polar_day ) then
                   iso_off = .8_r8 * pi
                   iso_on  = 1.2_r8 * pi
                else
                   iso_off = .8_r8 * sunoff
                   iso_on  = 2._r8 * pi - iso_off
                end if
                if( loc_angle >= iso_off .and. loc_angle <= iso_on ) then
                   sflx(i,isop_ndx) = 0._r8
                else
                   factor = loc_angle - iso_on
                   if( factor <= 0._r8 ) then
                      factor = factor + 2._r8*pi
                   end if
                   factor = factor / (2._r8*iso_off + 1.e-6_r8)
                   sflx(i,isop_ndx) = sflx(i,isop_ndx) * 2._r8 / iso_off * pi * (sin(pi*factor))**2
                end if
             else
                sflx(i,isop_ndx) = 0._r8
             end if
          end if
       end if

    end do

  end subroutine set_srf_emissions

end module mo_srf_emissions
