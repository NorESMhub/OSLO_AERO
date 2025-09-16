module aero_model

  !===============================================================================
  ! Oslo Aerosol Model
  ! Note: SPCAM not supported here
  !===============================================================================

  use shr_kind_mod,             only: r8 => shr_kind_r8
  use spmd_utils,               only: mpicom, mstrid=>masterprocid, masterproc
  use spmd_utils,               only: mpi_logical, mpi_real8, mpi_character, mpi_integer,  mpi_success
  use namelist_utils,           only: find_group_name
  use constituents,             only: pcnst, cnst_name, cnst_get_ind
  use ppgrid,                   only: pcols, pver, pverp, begchunk, endchunk
  use phys_control,             only: phys_getopts, cam_physpkg_is
  use phys_control,             only: history_aerosol_base
  use cam_abortutils,           only: endrun
  use cam_logfile,              only: iulog
  use perf_mod,                 only: t_startf, t_stopf
  use camsrfexch,               only: cam_in_t, cam_out_t
  use aerodep_flx,              only: aerodep_flx_prescribed
  use physics_types,            only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,           only: physics_buffer_desc, pbuf_get_field, pbuf_get_index, pbuf_set_field
  use physconst,                only: gravit, rair, rhoh2o, pi
  use time_manager,             only: get_nstep
  use cam_history,              only: outfld, fieldname_len, addfld, add_default, horiz_only
  use chem_mods,                only: gas_pcnst, adv_mass
  use mo_tracname,              only: solsym
  use mo_setsox,                only: setsox
  use mo_mass_xforms,           only: vmr2mmr, mmr2vmr, mmr2vmri
  use mo_chem_utls,             only: get_rxt_ndx, get_spc_ndx
  use ref_pres,                 only: top_lev => clim_modal_aero_top_lev
  use wv_saturation,            only: qsat_water
  !
  use oslo_aero_share,          only: nmodes_oslo=>nmodes, nbmodes
  use oslo_aero_share,          only: originalSigma, originalNumberMedianRadius
  use oslo_aero_share,          only: init_interp_constants
  use oslo_aero_share,          only: rTabMin, rTabMax, nk, normnk, rBinEdge, rBinMidpoint
  use oslo_aero_share,          only: volumeToNumber, numberToSurface, nBinsTab, rMinAquousChemistry
  use oslo_aero_share,          only: calculateLognormalCDF, calculatedNdLogR, calculateNumberConcentration
  use oslo_aero_share,          only: calculateEquivalentDensityOfFractalMode
  use oslo_aero_share,          only: numberFractionAvailableAqChem, isTracerInMode
  use oslo_aero_share,          only: chemistryIndex, physicsIndex
  use oslo_aero_share,          only: qqcw_get_field, numberOfProcessModeTracers
  use oslo_aero_share,          only: lifeCycleNumberMedianRadius, rhopart, lifeCycleSigma
  use oslo_aero_share,          only: l_so4_a2, l_bc_n, l_bc_ax, l_dms, l_isoprene, l_monoterp
  use oslo_aero_share,          only: MODE_IDX_BC_NUC, MODE_IDX_BC_EXT_AC
  use oslo_aero_share,          only: getNumberofTracersInMode, getCloudTracerIndexDirect, getCloudTracerName
  use oslo_aero_share,          only: getTracerIndex
  use oslo_aero_control,        only: oslo_aero_ctl_readnl, use_aerocom
  use oslo_aero_depos,          only: oslo_aero_depos_init
  use oslo_aero_depos,          only: oslo_aero_depos_dry, oslo_aero_depos_wet, oslo_aero_wetdep_init
  use oslo_aero_coag,           only: initializeCoagulation, coagtend, clcoag
  use oslo_aero_condtend,       only: N_COND_VAP, COND_VAP_ORG_SV, COND_VAP_ORG_LV, COND_VAP_H2SO4
  use oslo_aero_condtend,       only: initializeCondensation, condtend
  use oslo_aero_seasalt,        only: oslo_aero_seasalt_init, oslo_aero_seasalt_emis, seasalt_active
  use oslo_aero_dust,           only: oslo_aero_dust_init, oslo_aero_dust_emis
  use oslo_aero_ocean,          only: oslo_aero_ocean_init, oslo_aero_dms_emis
  use oslo_aero_sox_cldaero,    only: sox_cldaero_init
  use oslo_aero_microp,         only: oslo_aero_microp_readnl
  use oslo_aero_sw_tables,      only: initopt
  use oslo_aero_aerodry_tables, only: initdry
  use oslo_aero_aerocom_tables, only: initaeropt
  use oslo_aero_logn_tables,    only: initlogn

  implicit none
  private

  public :: aero_model_readnl
  public :: aero_model_register
  public :: aero_model_init
  public :: aero_model_gasaerexch     ! create, grow, change, and shrink aerosols.
  public :: aero_model_drydep         ! aerosol dry deposition and sediment
  public :: aero_model_wetdep         ! aerosol wet removal
  public :: aero_model_emissions      ! aerosol emissions
  public :: aero_model_surfarea       ! tropopspheric aerosol wet surface area for chemistry
  public :: aero_model_strat_surfarea ! stratospheric aerosol wet surface area for chemistry

  private :: aero_model_constants

  ! GS and AQ arrays used for source fields summed in oslo_aero_share.F90-summation_fields_writeout()
  ! GS_SOA = GS_SOA_LV + GS_SOA_SV, we use l_soa_lv and l_soa_sv for logic
  ! GS_H2SO4 and AQ_H2SO4, we use l_h2so4 for logic
  ! AQ_SO4_A2_OCW, we use l_so4_a2 for logic inside cloud tracer loop
  ! GS_DMS, we use l_dms for logic
  ! GS_SO2 and AQ_SO2, we us e l_so2 for logic
  real(r8), public, protected, allocatable   :: GS_SOA(:,:)
  real(r8), public, protected, allocatable   :: GS_H2SO4(:,:)
  real(r8), public, protected, allocatable   :: GS_DMS(:,:)
  real(r8), public, protected, allocatable   :: GS_SO2(:,:)
  real(r8), public, protected, allocatable   :: GS_isoprene(:,:)
  real(r8), public, protected, allocatable   :: GS_monoterp(:,:)
  real(r8), public, protected, allocatable   :: AQ_H2SO4(:,:)
  real(r8), public, protected, allocatable   :: AQ_SO4_A2_OCW(:,:)
  real(r8), public, protected, allocatable   :: AQ_SO2(:,:)

  ! Misc private data
  integer :: pblh_idx= 0
  integer :: ndx_h2so4 = -1  ! for surf_area_dens
  integer :: ndx_soa_lv = -1 ! for surf_area_dens
  integer :: ndx_soa_sv = -1 ! for surf_area_dens
  logical :: convproc_do_aer = .false.

  ! Namelist variables
  real(r8) :: sol_facti_cloud_borne   = 1._r8
  real(r8) :: sol_factb_interstitial  = 0.1_r8
  real(r8) :: sol_factic_interstitial = 0.4_r8

!=============================================================================
contains
!=============================================================================

  subroutine aero_model_readnl(nlfilename)
    ! OSLO_AERO begin
    use oslo_aero_dust,   only: oslo_aero_dust_readnl
    ! OSLO_AERO end

    ! read aerosol namelist options

    character(len=*), intent(in) :: nlfilename  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'aero_model_readnl'

    namelist /aerosol_nl/ sol_facti_cloud_borne, sol_factb_interstitial, sol_factic_interstitial
    !-----------------------------------------------------------------------------

    ! Read namelist
    if (masterproc) then
       open(newunit=unitn, file=trim(nlfilename), status='old' )
       call find_group_name(unitn, 'aerosol_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, aerosol_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if
    call mpi_bcast(sol_facti_cloud_borne, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_facti_cloud_borne")
    call mpi_bcast(sol_factb_interstitial, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_factb_interstitial")
    call mpi_bcast(sol_factic_interstitial, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: sol_factic_interstitial")

    call oslo_aero_ctl_readnl(nlfilename)
    call oslo_aero_microp_readnl(nlfilename)
   ! OSLO_AERO begin
   call oslo_aero_dust_readnl(nlfilename)
   ! OSLO_AERO end

  end subroutine aero_model_readnl

  !=============================================================================
  subroutine aero_model_register()
      use oslo_aero_share, only: aero_register,                   &
                                 sulfur_mass_fraction_register,   &
                                 soa_yield_register
      use oslo_aero_depos, only: oslo_aero_depos_register

      call aero_register()

      call oslo_aero_depos_register()

      call sulfur_mass_fraction_register()

      call soa_yield_register()

  end subroutine aero_model_register

  !=============================================================================
  subroutine aero_model_init( pbuf2d )

   use string_utils, only: int2str
   use mo_setsox,    only: sox_inti

   ! args
   type(physics_buffer_desc), pointer :: pbuf2d(:,:)

   ! local vars
   integer           :: icnst,id
   integer           :: l_aero
   character(len=20) :: dummy
   logical           :: history_aerosol ! Output MAM or SECT aerosol tendencies
   character(len=2)  :: unit_basename  ! Units 'kg' or '1'
   integer :: astat
   !------------------------------------
   ! allocate module variables
   allocate( GS_SOA(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_SOA array; error = '//int2str(astat))
   end if
   GS_SOA(:,:) = 0.0_r8
   allocate( GS_H2SO4(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_H2SO4 array; error = '//int2str(astat))
   end if
   GS_H2SO4(:,:) = 0.0_r8
   allocate( GS_DMS(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_DMS array; error = '//int2str(astat))
   end if
   GS_DMS(:,:) = 0.0_r8
   allocate( GS_SO2(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_SO2 array; error = '//int2str(astat))
   end if
   GS_SO2(:,:) = 0.0_r8
   allocate( GS_isoprene(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_isoprene array; error = '//int2str(astat))
   end if
   GS_isoprene(:,:) = 0.0_r8
   allocate( GS_monoterp(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate GS_monoterp array; error = '//int2str(astat))
   end if
   GS_monoterp(:,:) = 0.0_r8
   allocate( AQ_H2SO4(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate AQ_H2SO4 array; error = '//int2str(astat))
   end if
   AQ_H2SO4(:,:) = 0.0_r8
   allocate( AQ_SO4_A2_OCW(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate AQ_SO4_A2_OCW array; error = '//int2str(astat))
   end if
   AQ_SO4_A2_OCW(:,:) = 0.0_r8
   allocate( AQ_SO2(pcols, begchunk:endchunk), stat=astat )
   if( astat/= 0 ) then
      call endrun('aero_model_init: failed to allocate AQ_SO2 array; error = '//int2str(astat))
   end if
   AQ_SO2(:,:) = 0.0_r8
    !------------------------------------

    call phys_getopts(history_aerosol_out=history_aerosol, convproc_do_aer_out=convproc_do_aer)

    ! aqueous chem initialization
    call sox_inti()

    call aero_model_constants
    call init_interp_constants() ! table initialization constants
    call initopt()               ! table initialization
    call initlogn()              ! table initialization
    if (use_aerocom) then
       call initdry()               ! table initialization
       call initaeropt()            ! table initialization
    end if
    call initializeCondensation()
    call oslo_aero_ocean_init()
    call oslo_aero_depos_init(pbuf2d)
    call oslo_aero_dust_init()
    call oslo_aero_seasalt_init()
    call oslo_aero_wetdep_init()

    dummy = 'RAM1'
    call addfld (dummy,horiz_only, 'A','frac','RAM1')
    if ( history_aerosol ) then
       call add_default (dummy, 1, ' ')
    endif

    dummy = 'airFV'
    call addfld (dummy,horiz_only, 'A','frac','FV')
    if ( history_aerosol ) then
       call add_default (dummy, 1, ' ')
    endif

    ! Get height of boundary layer for boundary layer nucleation
    pblh_idx = pbuf_get_index('pblh')

    call cnst_get_ind ( "H2SO4", ndx_h2so4, abort=.true. )
    ndx_h2so4 = chemistryIndex(ndx_h2so4)
    call cnst_get_ind ( "SOA_LV", ndx_soa_lv,abort=.true.)
    ndx_soa_lv = chemistryIndex(ndx_soa_lv)
    call cnst_get_ind ( "SOA_SV", ndx_soa_sv, abort=.true.)
    ndx_soa_sv = chemistryIndex(ndx_soa_sv)

    do icnst = 1,gas_pcnst
       unit_basename = 'kg'  ! Units 'kg' or '1'

       call addfld( 'GS_'//trim(solsym(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
            trim(solsym(icnst))//' gas chemistry/wet removal (for gas species)')

       call addfld( 'AQ_'//trim(solsym(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
            trim(solsym(icnst))//' aqueous chemistry (for gas species)')

       if(physicsIndex(icnst)<=pcnst) then
          if (getCloudTracerIndexDirect(physicsIndex(icnst)) > 0)then
             call addfld( 'AQ_'//getCloudTracerName(physicsIndex(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
                  trim(solsym(icnst))//' aqueous chemistry (for cloud species)')
          end if
       end if

       if ( history_aerosol ) then
          call add_default( 'GS_'//trim(solsym(icnst)), 1, ' ')
          call add_default( 'AQ_'//trim(solsym(icnst)), 1, ' ')
          if(physicsIndex(icnst)<=pcnst) then
             if(getCloudTracerIndexDirect(physicsIndex(icnst))>0)then
                call add_default( 'AQ_'//getCloudTracerName(physicsIndex(icnst)),1,' ')
             end if
          end if
       endif

      call cnst_get_ind(solsym(icnst), l_aero, abort=.false. )
      if ( l_aero == l_dms .or. l_aero == l_isoprene .or. l_aero == l_monoterp ) then
         call addfld( 'sink_'//trim(solsym(icnst)),horiz_only, 'A', unit_basename//'/m2/s ', &
            trim(solsym(icnst))//' sink')
         if ( l_aero == l_dms ) then
            call addfld( 'sink_'//trim(solsym(icnst))//'_S',horiz_only, 'A', 'kg*S/m2/s ', &
               trim(solsym(icnst))//' sink, sulfur mass only')
         endif
         if ( history_aerosol_base ) then
            call add_default( 'sink_'//trim(solsym(icnst)), 1, ' ')
            if ( l_aero == l_dms ) then
               call add_default( 'sink_'//trim(solsym(icnst))//'_S', 1, ' ')
            endif
         endif
      endif
    enddo

    call addfld ('NUCLRATE',(/'lev'/), 'A','#/cm3/s','Nucleation rate')
    call addfld ('FORMRATE',(/'lev'/), 'A','#/cm3/s','Formation rate of 12nm particles')
    call addfld ('COAGNUCL',(/'lev'/), 'A', '/s','Coagulation sink for nucleating particles')
    call addfld ('GRH2SO4',(/'lev'/), 'A', 'nm/hour','Growth rate H2SO4')
    call addfld ('GRSOA',(/'lev'/),'A','nm/hour','Growth rate SOA')
    call addfld ('GR',(/'lev'/), 'A', 'nm/hour','Growth rate, H2SO4+SOA')
    call addfld ('NUCLSOA',(/'lev'/),'A','kg/kg','SOA nucleate')
    call addfld ('ORGNUCL',(/'lev'/),'A','kg/kg','Organic gas available for nucleation')

    if(history_aerosol)then
       call add_default ('NUCLRATE', 1, ' ')
       call add_default ('FORMRATE', 1, ' ')
       call add_default ('COAGNUCL', 1, ' ')
       call add_default ('GRH2SO4', 1, ' ')
       call add_default ('GRSOA', 1, ' ')
       call add_default ('GR', 1, ' ')
       call add_default ('NUCLSOA', 1, ' ')
       call add_default ('ORGNUCL', 1, ' ')
    end if

    call addfld( 'XPH_LWC',    (/ 'lev' /), 'A','kg/kg',   'pH value multiplied by lwc')
    call addfld ('AQSO4_H2O2', horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to H2O2')
    call addfld ('AQSO4_O3',   horiz_only,  'A','kg/m2/s', 'SO4 aqueous phase chemistry due to O3')

    if ( history_aerosol ) then
       call add_default ('XPH_LWC', 1, ' ')
       call add_default ('AQSO4_H2O2', 1, ' ')
       call add_default ('AQSO4_O3', 1, ' ')
    endif

  end subroutine aero_model_init

  !=============================================================================
  subroutine aero_model_drydep  ( state, pbuf, obklen, ustar, cam_in, dt, cam_out, ptend )

    ! args
    type(physics_state),    intent(in)    :: state    ! Physics state variables
    real(r8),               intent(in)    :: obklen(:)
    real(r8),               intent(in)    :: ustar(:) ! sfc fric vel
    type(cam_in_t), target, intent(in)    :: cam_in   ! import state
    real(r8),               intent(in)    :: dt       ! time step
    type(cam_out_t),        intent(inout) :: cam_out  ! export state
    type(physics_ptend),    intent(out)   :: ptend    ! indivdual parameterization tendencies
    type(physics_buffer_desc),    pointer :: pbuf(:)

    ! local vars
    real(r8) :: oslo_dgnumwet(pcols, pver, 0:nmodes_oslo)
    real(r8) :: oslo_wetdens(pcols, pver, 0:nmodes_oslo)
    real(r8) :: oslo_dgnumwet_processmodes(pcols, pver, numberOfProcessModeTracers)
    real(r8) :: oslo_wetdens_processmodes(pcols, pver, numberOfProcessModeTracers)

    oslo_wetdens(:,:,:) = 0._r8
    call calcaersize_sub(state%ncol, state%t, state%q(1,1,1), state%pmid, state%pdel, &
         oslo_dgnumwet, oslo_wetdens, oslo_dgnumwet_processmodes, oslo_wetdens_processmodes)

    call oslo_aero_depos_dry(state%lchnk, state%ncol, state%psetcols, &
         state%t, state%pmid, state%pdel, state%pint, state%q, &
         cam_in%landfrac, cam_in%icefrac, cam_in%ocnfrac, cam_in%fv, cam_in%ram1, cam_in%cflx, &
         pbuf, obklen, ustar, dt, &
         oslo_dgnumwet, oslo_wetdens, oslo_dgnumwet_processmodes, oslo_wetdens_processmodes, &
         cam_out, ptend)

  end subroutine aero_model_drydep

  !=============================================================================
  subroutine aero_model_wetdep( state, dt, dlf, cam_out, ptend, pbuf)

    type(physics_state), intent(in)    :: state       ! Physics state variables
    real(r8),            intent(in)    :: dt          ! time step
    real(r8),            intent(in)    :: dlf(:,:)    ! shallow+deep convective detrainment [kg/kg/s]
    type(cam_out_t),     intent(inout) :: cam_out     ! export state
    type(physics_ptend), intent(out)   :: ptend       ! indivdual parameterization tendencies
    type(physics_buffer_desc), pointer :: pbuf(:)

    call oslo_aero_depos_wet(state%lchnk, state%ncol, state%psetcols, state%pmid, state%pdel, state%q, state%t, &
         dt, dlf, cam_out, ptend, pbuf)

  end subroutine aero_model_wetdep

  !============================================================================
  subroutine aero_model_surfarea(state, mmr, radmean, relhum, pmid, temp, &
       strato_sad, sulfate, rho, ltrop, dlat, het1_ndx, pbuf, ncol, sfc, &
       dm_aer, sad_trop, reff_trop )

    !-------------------------------------------------------------------------
    ! provides wet tropospheric aerosol surface area info for modal aerosols
    ! called from mo_usrrxt
    !-------------------------------------------------------------------------

    ! arguments
    type(physics_state), intent(in) :: state        ! Physics state variables
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: radmean      ! mean radii in cm
    real(r8), intent(in)    :: strato_sad(:,:)
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: ltrop(:)
    real(r8), intent(in)    :: dlat(:)                    ! degrees latitude
    integer,  intent(in)    :: het1_ndx
    real(r8), intent(in)    :: relhum(:,:)
    real(r8), intent(in)    :: rho(:,:) ! total atm density (/cm^3)
    real(r8), intent(in)    :: sulfate(:,:)
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(inout) :: sfc(:,:,:)
    real(r8), intent(inout) :: dm_aer(:,:,:)
    real(r8), intent(inout) :: sad_trop(:,:)
    real(r8), intent(out)   :: reff_trop(:,:)

    ! local vars
    ! HAVE TO GET RID OF THIS MODE 0!! MESSES UP EVERYTHING!!
    real(r8)         :: numberConcentration(pcols,pver,0:nmodes_oslo)
    real(r8), target :: sad_mode(pcols,pver,nmodes_oslo)
    real(r8)         :: rho_air(pcols,pver)
    integer          :: itrac,imode,icol,ilev

    ! Get air density
    do ilev = 1,pver
       do icol=1,ncol
          rho_air(icol,ilev) = pmid(icol,ilev)/(temp(icol,ilev)*287.04_r8)
       end do
    end do

    ! Get number concentrations
    call calculateNumberConcentration(ncol, mmr, rho_air, numberConcentration)

    ! Convert to area using lifecycle-radius
    sad_mode = 0._r8
    sad_trop = 0._r8
    do imode=1,nmodes_oslo
       do ilev=1,pver
          sad_mode(:ncol,ilev,imode) = numberConcentration(:ncol,ilev,imode)*numberToSurface(imode)*1.e-2_r8 !m2/m3 ==> cm2/cm3
          sad_trop(:ncol,ilev) = sad_trop(:ncol,ilev) + sad_mode(:ncol,ilev,imode)
       end do
    end do

    do imode=1,nmodes_oslo
       do ilev=1,pver
          sfc(:ncol,ilev,imode) = sad_mode(:ncol,ilev,imode)     ! aitken_idx:aitken_idx)
          dm_aer(:ncol,ilev,imode) = 2.0_r8*lifeCycleNumberMedianRadius(imode)
       end do
    end do

    ! Need to implement reff_trop here
    reff_trop(:,:) = 1.0e-6_r8

  end subroutine aero_model_surfarea

  !=============================================================================
  subroutine aero_model_strat_surfarea( state, ncol, mmr, pmid, temp, &
       ltrop, pbuf, strato_sad, reff_strat )

    !-------------------------------------------------------------------------
    ! provides WET stratospheric aerosol surface area info for modal aerosols
    ! if modal_strat_sulfate = TRUE -- called from mo_gas_phase_chemdr
    !-------------------------------------------------------------------------

    ! arguments
    type(physics_state), intent(in) :: state        ! Physics state variables
    integer,  intent(in)    :: ncol
    real(r8), intent(in)    :: mmr(:,:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: temp(:,:)
    integer,  intent(in)    :: ltrop(:) ! tropopause level indices
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8), intent(out)   :: strato_sad(:,:)
    real(r8), intent(out)   :: reff_strat(:,:)

    reff_strat = 0.1e-6_r8
    strato_sad = 0._r8

  end subroutine aero_model_strat_surfarea

  !===========================================================================
  subroutine aero_model_gasaerexch( state, loffset, ncol, lchnk, troplev,     &
       delt, reaction_rates, tfld, pmid, pdel, mbar, relhum, zm,  qh2o, cwat, &
       cldfr, cldnum, airdens, invariants, del_h2so4_gasprod, vmr0, vmr, pbuf )

    ! use
    use oslo_aero_share,    only : sulfurMassFraction
    use oslo_aero_share,    only : l_soa_lv, l_soa_sv, l_h2so4, l_so4_a2, l_dms, l_isoprene, l_monoterp, l_so2

    ! arguments
    type(physics_state), intent(in) :: state        ! Physics state variables
    integer,  intent(in)    :: loffset                ! offset applied to modal aero "pointers"
    integer,  intent(in)    :: ncol                   ! number columns in chunk
    integer,  intent(in)    :: lchnk                  ! chunk index
    integer,  intent(in)    :: troplev(pcols)
    real(r8), intent(in)    :: delt                   ! time step size (sec)
    real(r8), intent(in)    :: reaction_rates(:,:,:)  ! reaction rates
    real(r8), intent(in)    :: tfld(:,:)              ! temperature (K)
    real(r8), intent(in)    :: pmid(:,:)              ! pressure at model levels (Pa)
    real(r8), intent(in)    :: pdel(:,:)              ! pressure thickness of levels (Pa)
    real(r8), intent(in)    :: mbar(:,:)              ! mean wet atmospheric mass ( amu )
    real(r8), intent(in)    :: relhum(:,:)            ! relative humidity
    real(r8), intent(in)    :: airdens(:,:)           ! total atms density (molec/cm**3)
    real(r8), intent(in)    :: invariants(:,:,:)
    real(r8), intent(in)    :: zm(:,:)
    real(r8), intent(in)    :: qh2o(:,:)
    real(r8), intent(in)    :: cwat(:,:)              ! cloud liquid water content (kg/kg)
    real(r8), intent(in)    :: cldfr(:,:)
    real(r8), intent(in)    :: cldnum(:,:)            ! droplet number concentration (#/kg)
    real(r8), intent(inout) :: del_h2so4_gasprod(:,:) ! [molec/molec/sec]
    real(r8), intent(in)    :: vmr0(:,:,:)            ! initial mixing ratios (before gas-phase chem changes)
    real(r8), intent(inout) :: vmr(:,:,:)             ! mixing ratios ( vmr )
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! local vars
    integer, parameter :: nmodes_aq_chem = 1
    integer  :: icol,ilev
    integer  :: l_aero
    integer  :: imode,icnst,itrac
    integer  :: nstep
    real(r8) :: wrk(ncol)
    real(r8) :: dvmrcwdt(ncol,pver,gas_pcnst)
    real(r8) :: dvmrdt(ncol,pver,gas_pcnst)
    real(r8) :: vmrcw(ncol,pver,gas_pcnst)   ! cloud-borne aerosol (vmr)
    real(r8) :: del_h2so4_aeruptk(ncol,pver)
    real(r8) :: del_h2so4_aqchem(ncol,pver)
    real(r8) :: mmr_cond_vap_start_of_timestep(pcols,pver,N_COND_VAP)
    real(r8) :: mmr_cond_vap_gasprod(pcols,pver,N_COND_VAP)
    real(r8) :: del_soa_lv_gasprod(ncol,pver)
    real(r8) :: del_soa_sv_gasprod(ncol,pver)
    real(r8) :: dvmrdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: dvmrcwdt_sv1(ncol,pver,gas_pcnst)
    real(r8) :: mmr_tend_ncols(ncol, pver, gas_pcnst)
    real(r8) :: mmr_tend_pcols(pcols, pver, gas_pcnst)
    integer  :: cond_vap_idx
    real(r8) :: aqso4(ncol,nmodes_aq_chem)   ! aqueous phase chemistry
    real(r8) :: aqh2so4(ncol,nmodes_aq_chem) ! aqueous phase chemistry
    real(r8) :: aqso4_h2o2(ncol)             ! SO4 aqueous phase chemistry due to H2O2
    real(r8) :: aqso4_o3(ncol)               ! SO4 aqueous phase chemistry due to O3
    real(r8) :: xphlwc(ncol,pver)            ! pH value multiplied by lwc
    real(r8) :: delt_inverse                 ! 1 / timestep
    real(r8), pointer :: pblh(:)
    character(len=32) :: name

    nstep = get_nstep()

    delt_inverse = 1.0_r8 / delt
    GS_SOA(:ncol,lchnk) = 0._r8
    GS_H2SO4(:ncol,lchnk) = 0._r8
    GS_DMS(:ncol,lchnk) = 0._r8
    GS_SO2(:ncol,lchnk) = 0._r8
    GS_isoprene(:ncol,lchnk) = 0._r8
    GS_monoterp(:ncol,lchnk) = 0._r8
    AQ_H2SO4(:ncol,lchnk) = 0._r8
    AQ_SO4_A2_OCW(:ncol,lchnk) = 0._r8
    AQ_SO2(:ncol,lchnk) = 0._r8

    ! Get height of boundary layer (needed for boundary layer nucleation)
    call pbuf_get_field(pbuf, pblh_idx, pblh)

    ! calculate tendency due to gas phase chemistry and processes
    dvmrdt(:ncol,:,:) = (vmr(:ncol,:,:) - vmr0(:ncol,:,:)) / delt
    do icnst = 1, gas_pcnst
       wrk(:) = 0._r8
       do ilev = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
       end do

       call cnst_get_ind(trim(solsym(icnst)), l_aero, abort=.false.)
       if ( l_aero == l_soa_lv .or. l_aero == l_soa_sv ) then
          GS_SOA(:ncol,lchnk) = GS_SOA(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_h2so4 ) then
          GS_H2SO4(:ncol,lchnk) = GS_H2SO4(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_dms ) then
          GS_DMS(:ncol, lchnk) = GS_DMS(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_so2) then
          GS_SO2(:ncol, lchnk) = GS_SO2(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_isoprene ) then
          GS_isoprene(:ncol, lchnk) = GS_isoprene(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_monoterp ) then
          GS_monoterp(:ncol, lchnk) = GS_monoterp(:ncol,lchnk) + wrk(:ncol)
       endif

       name = 'GS_'//trim(solsym(icnst))
       call outfld( name, wrk(:ncol), ncol, lchnk )

       if ( l_aero == l_dms .or. l_aero == l_isoprene .or. l_aero == l_monoterp) then
         call outfld( 'sink_'//trim(solsym(icnst)), wrk(:ncol), ncol, lchnk )
         if ( l_aero == l_dms ) then
            call outfld( 'sink_'//trim(solsym(icnst))//'_S', ( wrk(:ncol) * sulfurMassFraction(l_aero) ) , ncol, lchnk )
         endif
      endif
    enddo

    ! Get mass mixing ratios at start of time step
    call vmr2mmr( vmr0, mmr_tend_ncols, mbar, ncol )
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_H2SO4)  = mmr_tend_ncols(1:ncol,:,ndx_h2so4)
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_LV) = mmr_tend_ncols(1:ncol,:,ndx_soa_lv)
    mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_SV) = mmr_tend_ncols(1:ncol,:,ndx_soa_sv)
    !
    ! Aerosol processes ...
    call qqcw2vmr( lchnk, vmrcw, mbar, ncol, loffset, pbuf )

    ! save h2so4 change by gas phase chem (for later new particle nucleation)
    if (ndx_h2so4 > 0) then
       del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - vmr0(1:ncol,:,ndx_h2so4)
    endif

    del_soa_lv_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_soa_lv) - vmr0(1:ncol,:,ndx_soa_lv)
    del_soa_sv_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_soa_sv) - vmr0(1:ncol,:,ndx_soa_sv)

    dvmrdt(:ncol,:,:) = vmr(:ncol,:,:)
    dvmrcwdt(:ncol,:,:) = vmrcw(:ncol,:,:)

    !Save intermediate concentrations
    dvmrdt_sv1 = vmr
    dvmrcwdt_sv1 = vmrcw

    ! aqueous chemistry ...
    call setsox( state, ncol, lchnk, loffset, delt, pmid, pdel, tfld, mbar, cwat, &
         cldfr, cldnum, airdens, invariants, vmrcw, vmr, xphlwc, &
         aqso4, aqh2so4, aqso4_h2o2, aqso4_o3)

    call outfld( 'AQSO4_H2O2', aqso4_h2o2(:ncol), ncol, lchnk)
    call outfld( 'AQSO4_O3',   aqso4_o3(:ncol),   ncol, lchnk)
    call outfld( 'XPH_LWC',    xphlwc(:ncol,:),   ncol, lchnk )

    ! vmr tendency from aqchem and soa routines
    dvmrdt_sv1 = (vmr - dvmrdt_sv1)/delt
    dvmrcwdt_sv1 = (vmrcw - dvmrcwdt_sv1)/delt

    if (ndx_h2so4 > 0)then
       del_h2so4_aqchem(:ncol,:) = dvmrdt_sv1(:ncol,:,ndx_h2so4)*delt !"production rate" of H2SO4
    else
       del_h2so4_aqchem(:ncol,:) = 0.0_r8
    end if

    do icnst = 1,gas_pcnst
       wrk(:ncol) = 0._r8
       do ilev = 1,pver
          wrk(:ncol) = wrk(:ncol) + dvmrdt_sv1(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
       end do

       call cnst_get_ind(trim(solsym(icnst)), l_aero, abort=.false.)
       if ( l_aero == l_h2so4) then
          AQ_H2SO4(:ncol,lchnk) = AQ_H2SO4(:ncol,lchnk) + wrk(:ncol)
       else if ( l_aero == l_so2 ) then
          AQ_SO2(:ncol,lchnk) = AQ_SO2(:ncol,lchnk) + wrk(:ncol)
       endif

       name = 'AQ_'//trim(solsym(icnst))
       call outfld( name, wrk(:ncol), ncol, lchnk )

       !In oslo aero also write out the tendencies for the
       !cloud borne aerosols...
       itrac = physicsIndex(icnst)
       if (itrac <= pcnst) then
          if (getCloudTracerIndexDirect(itrac) > 0)then
             name = 'AQ_'//trim(getCloudTracerName(itrac))
             wrk(:ncol)=0.0_r8
             do ilev=1,pver
                wrk(:ncol) = wrk(:ncol) + dvmrcwdt_sv1(:ncol,ilev,icnst)*adv_mass(icnst)/mbar(:ncol,ilev)*pdel(:ncol,ilev)/gravit
             end do

             if ( l_aero == l_so4_a2 ) then
               AQ_SO4_A2_OCW(:ncol,lchnk) = AQ_SO4_A2_OCW(:ncol,lchnk) + wrk(:ncol)
             endif
             call outfld( name, wrk(:ncol), ncol, lchnk )
          end if
       end if
    enddo

    ! condensation
    call vmr2mmr( vmr, mmr_tend_ncols, mbar, ncol )
    do ilev = 1,pver
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_H2SO4) = adv_mass(ndx_h2so4) &
            * (del_h2so4_gasprod(:ncol,ilev)+del_h2so4_aqchem(:ncol,ilev)) / mbar(:ncol,ilev)/delt
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_ORG_LV) = adv_mass(ndx_soa_lv) &
            * del_soa_lv_gasprod(:ncol,ilev) / mbar(:ncol,ilev)/delt
       mmr_cond_vap_gasprod(:ncol,ilev,COND_VAP_ORG_SV) = adv_mass(ndx_soa_sv) &
            * del_soa_sv_gasprod(:ncol,ilev) / mbar(:ncol,ilev)/delt
    end do

    ! This should not happen since there are only production terms for these gases! !
    do cond_vap_idx=1,N_COND_VAP
       where(mmr_cond_vap_gasprod(:ncol,:,cond_vap_idx) < 0.0_r8)
          mmr_cond_vap_gasprod(:ncol,:,cond_vap_idx) = 0.0_r8
       end where
    end do
    mmr_tend_ncols(:ncol,:,ndx_h2so4)  = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_H2SO4)
    mmr_tend_ncols(:ncol,:,ndx_soa_lv) = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_LV)
    mmr_tend_ncols(:ncol,:,ndx_soa_sv) = mmr_cond_vap_start_of_timestep(:ncol,:,COND_VAP_ORG_SV)

    ! Rest of microphysics have pcols dimension
    mmr_tend_pcols(:ncol,:,:) = mmr_tend_ncols(:ncol,:,:)

    ! Condensation
    ! Note use of "zm" here. In CAM5.3-implementation "zi" was used..
    ! zm is passed through the generic interface, and it should not change much
    ! to check if "zm" is below boundary layer height instead of zi
    call condtend( lchnk, mmr_tend_pcols, mmr_cond_vap_gasprod,tfld, pmid, &
         pdel, delt, ncol, pblh, zm, qh2o)  ! cka

    ! Coagulation
    ! OS 280415  Concentratiions in cloud water is in vmr space and as a
    ! temporary variable  (vmrcw) Coagulation between aerosol and cloud
    ! droplets moved to after vmrcw is moved into qqcw (in mmr spac)
    call coagtend( mmr_tend_pcols, pmid, pdel, tfld, delt_inverse, ncol, lchnk)

    ! Convert cloud water to mmr again ==> values in buffer
    call vmr2qqcw( vmrcw, mbar, ncol, loffset, pbuf )

    ! Call cloud coagulation routines (all in mass mixing ratios)
    call clcoag( mmr_tend_pcols, pmid, pdel, tfld, cldnum ,cldfr, delt_inverse, ncol, lchnk, loffset, pbuf)

    ! Make sure mmr==> vmr is done correctly
    mmr_tend_ncols(:ncol,:,:) = mmr_tend_pcols(:ncol,:,:)

    ! Go back to volume mixing ratio for chemistry
    call mmr2vmr( mmr_tend_ncols, vmr, mbar, ncol )

  end subroutine aero_model_gasaerexch

  !=============================================================================
  subroutine aero_model_emissions( state, cam_in )
     use oslo_aero_control, only: dms_from_ocn
     use constituents     , only: cnst_get_ind, sflxnam

    ! Arguments:
    type(physics_state), intent(in)    :: state   ! Physics state variables
    type(cam_in_t),      intent(inout) :: cam_in  ! import state

    ! Local variables
    integer :: icol
    integer :: pndx_fdms  ! DMS surface flux physics index

    call oslo_aero_dust_emis( state%lchnk, state%ncol, cam_in%dstflx, cam_in%cflx)

    if (seasalt_active) then
       call oslo_aero_seasalt_emis(state%ncol, state%lchnk, &
            state%u(:,pver), state%v(:,pver), state%zm(:,pver), &
            cam_in%ocnfrac, cam_in%icefrac, cam_in%sst, cam_in%cflx)
    end if

    ! Pick up correct DMS emissions (replace values from file if requested)
    ! Update cam_in%clfx for pndx_dms when dms is read in or obtained from the ocean
    if (.not. dms_from_ocn) then
       call oslo_aero_dms_emis(state%ncol, state%lchnk, &
            state%u(:,pver), state%v(:,pver), state%zm(:,pver), &
            cam_in%ocnfrac, cam_in%icefrac, cam_in%sst, cam_in%cflx)
    else
       call cnst_get_ind('DMS', pndx_fdms, abort=.true.)
       do icol = 1,state%ncol
          cam_in%cflx(icol,pndx_fdms) = cam_in%fdms(icol)
          ! The addfld call for 'odms' below is in the routine
          ! oslo_aero_ocean_init in module oslo_aero_ocean.F90.
          call outfld('odms', cam_in%fdms(:state%ncol), state%ncol, state%lchnk)
       end do
    end if

  end subroutine aero_model_emissions

  !=============================================================================
  ! private methods
  !=============================================================================

  subroutine qqcw2vmr(lchnk, vmr, mbar, ncol, im, pbuf)

    !-----------------------------------------------------------------
    !	... Xfrom from mass to volume mixing ratio
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !	... Dummy args
    !-----------------------------------------------------------------
    integer, intent(in)     :: lchnk, ncol, im
    real(r8), intent(in)    :: mbar(ncol,pver)
    real(r8), intent(inout) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    !-----------------------------------------------------------------
    !	... Local variables
    !-----------------------------------------------------------------
    integer :: ilev, icnst
    real(r8), pointer :: fldcw(:,:)

    do icnst=1,gas_pcnst
       if( adv_mass(icnst) /= 0._r8 ) then
          fldcw => qqcw_get_field(pbuf, icnst+im)
          if(associated(fldcw)) then
             do ilev = 1,pver
                vmr(:ncol,ilev,icnst) = mbar(:ncol,ilev) * fldcw(:ncol,ilev) / adv_mass(icnst)
             end do
          else
             vmr(:,:,icnst) = 0.0_r8
          end if
       end if
    end do
  end subroutine qqcw2vmr

  !=============================================================================
  subroutine vmr2qqcw( vmr, mbar, ncol, im, pbuf )

    !-----------------------------------------------------------------
    ! Convert from volume to mass mixing ratio
    !-----------------------------------------------------------------

    ! Arguments
    integer , intent(in) :: ncol, im
    real(r8), intent(in) :: mbar(ncol,pver)
    real(r8), intent(in) :: vmr(ncol,pver,gas_pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    ! Local variables
    integer :: ilev ! level index
    integer :: icnst ! tracer index
    real(r8), pointer :: fldcw(:,:)

    ! The non-group species
    do icnst = 1,gas_pcnst
       fldcw => qqcw_get_field(pbuf, icnst+im)
       if( adv_mass(icnst) /= 0._r8 .and. associated(fldcw)) then
          do ilev = 1,pver
             fldcw(:ncol,ilev) = adv_mass(icnst) * vmr(:ncol,ilev,icnst) / mbar(:ncol,ilev)
          end do
       end if
    end do
  end subroutine vmr2qqcw

  !=============================================================================
  subroutine aero_model_constants()

    ! A number of constants used in the emission and size-calculation in CAM-Oslo

    ! local variables
    integer  :: imode,ibin
    real(r8) :: rhob(0:nmodes_oslo) !density of background aerosol in mode
    real(r8) :: rhorbc         !This has to do with fractal dimensions of bc, come back to this!!
    real(r8) :: sumnormnk
    real(r8) :: totalLogDelta
    real(r8) :: logDeltaBin
    real(r8) :: logNextEdge
    !---------------------------------------

    rhob(:)            =-1.0_r8
    volumeToNumber(:)  =-1.0_r8
    numberToSurface(:) =-1.0_r8

    !Prepare modal properties
    do imode =0,nmodes_oslo
       if(getNumberOfTracersInMode(imode) > 0)then

          !Approximate density of mode
          !density of mode is density of first species in mode
          rhob(imode)  = rhopart(getTracerIndex(imode,1,.false.))

          !REPLACE THE EFACT-VARIABLE WITH THIS!!
          volumeToNumber(imode) = 1.0_r8 / &
               ( DEXP ( 4.5_r8 * ( log(originalSigma(imode)) * log(originalSigma(imode)) ) ) &
               *(4.0_r8/3.0_r8)*pi*(originalNumberMedianRadius(imode))**3 )

          numberToSurface(imode) = 4.0_r8*pi*lifeCycleNumberMedianRadius(imode)*lifeCycleNumberMedianRadius(imode)&
               *DEXP(log(lifeCycleSigma(imode))*log(lifeCycleSigma(imode)))
       end if
    end do

    !Find radius in edges and midpoints of bin
    rBinEdge(1) = rTabMin
    totalLogDelta = log(rTabMax/rTabMin)
    logDeltaBin = totalLogDelta / nBinsTab
    do ibin = 2,nBinsTab+1
       logNextEdge = log(rBinEdge(ibin-1)) + logDeltaBin
       rBinEdge(ibin) = DEXP(logNextEdge)
       rBinMidPoint(ibin-1) = sqrt(rBinEdge(ibin)*rBinEdge(ibin-1))
    end do

    !Calculate the fraction of a mode which goes to aquous chemstry
    numberFractionAvailableAqChem(:)=0.0_r8
    do imode = 1,nbmodes
       if(isTracerInMode(imode,l_so4_a2))then
          numberFractionAvailableAqChem(imode) =  1.0_r8 -  &
               calculateLognormalCDF(rMinAquousChemistry,originalNumberMedianRadius(imode), originalSigma(imode))
       end if
    end do

    !Set the density of the fractal mode ==> we get lesser density
    !than the emitted density, so for a given mass emitted, we get
    !more number-concentration!! This is a way of simulating that the
    !aerosols take up more space
    rhorbc = calculateEquivalentDensityOfFractalMode(    &
         rhopart(l_bc_n),                                & !emitted density
         originalNumberMedianRadius(MODE_IDX_BC_NUC),    & !emitted size
         2.5_r8,                                         & !fractal dim
         originalNumberMedianRadius(MODE_IDX_BC_EXT_AC), & !diameter of mode
         originalSigma(MODE_IDX_BC_EXT_AC))                !sigma mode

    rhopart(l_bc_ax) = rhorbc
    !fxm: not the right place for this change of value,
    !but anyway.. this re-calculateion of tracer density
    !influences density of mode used in coagulation
    rhob(MODE_IDX_BC_EXT_AC)=rhorbc

    !Size distribution of the modes!
    !Unclear if this should use the radii assuming growth or not!
    !Mostly used in code where it is sensible to assume some growth has
    !happened, so it is used here
    do imode = 0,nmodes_oslo
       do ibin=1,nBinsTab
          !dN/dlogR (does not sum to one over size range)
          nk(imode,ibin) = calculatedNdLogR(rBinMidPoint(ibin), lifeCycleNumberMedianRadius(imode), lifeCycleSigma(imode))

          !dN (sums to one) over the size range
          normnk(imode,ibin) =logDeltaBin*nk(imode,ibin)
       enddo
    enddo  ! imode

    ! Normalized size distribution must sum to one (accept 2% error)
    do imode=0,nmodes_oslo
       sumNormNk = sum(normnk(imode,:))
       if(abs(sum(normnk(imode,:)) - 1.0_r8) > 2.0e-2_r8)then
          print*, "sum normnk", sum(normnk(imode,:))
          call endrun()
       endif
    enddo

    ! Initialize coagulation including Calculate the coagulation coefficients
    ! Note: Inaccurate density used!
    call initializeCoagulation(rhob, lifeCycleNumberMedianRadius)

  end subroutine aero_model_constants

  !=============================================================================
  subroutine calcaersize_sub( ncol, t, h2ommr, pmid, pdel,wetnumberMedianDiameter,wetrho,   &
       wetNumberMedianDiameter_processmode, wetrho_processmode)

     use oslo_aero_share, only: max_tracers_per_mode, rhtab, rdivr0
     use oslo_aero_share, only: getNumberOfBackgroundTracersInMode
     use oslo_aero_share, only: getDryDensity, tracerInProcessMode
     use oslo_aero_share, only: processModeNumberMedianRadius, processModeSigma

    ! Seland Calculates mean volume size and hygroscopic growth for use in  dry deposition

    integer,  intent(in)  :: ncol               ! number of columns
    real(r8), intent(in)  :: t(pcols,pver)      ! layer temperatures (K)
    real(r8), intent(in)  :: h2ommr(pcols,pver) ! layer specific humidity
    real(r8), intent(in)  :: pmid(pcols,pver)   ! layer pressure (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)  ! layer pressure thickness (Pa)
    real(r8), intent(out) :: wetNumberMedianDiameter(pcols,pver,0:nmodes_oslo)
    real(r8), intent(out) :: wetrho(pcols,pver,0:nmodes_oslo) ! wet aerosol density
    real(r8), intent(out) :: wetNumberMedianDiameter_processmode(pcols,pver,numberOfProcessModeTracers)
    real(r8), intent(out) :: wetrho_processmode(pcols,pver,numberOfProcessModeTracers)

    !     local variables
    integer  :: icol,itrac,ilev,imode,irelh,mm, tracerCounter
    real(r8) :: xrh(pcols,pver)
    real(r8) :: relhum(pcols,pver) ! Relative humidity
    real(r8) :: qs(pcols,pver)     ! saturation specific humidity
    real(r8) :: rmeanvol           ! Mean radius with respect to volume
    integer  :: irh1(pcols,pver),irh2(pcols,pver)
    integer  :: t_irh1,t_irh2
    real(r8) :: t_rh1,t_rh2,t_xrh,rr1,rr2
    real(r8) :: volumeFractionAerosol   !with respect to total (aerosol + water)
    real(r8) :: tmp1
    real(r8) :: wetrad_tmp(max_tracers_per_mode)
    real(r8) :: dry_rhopart_tmp(max_tracers_per_mode)
    real(r8) :: mixed_dry_rho


    !Get the tabulated rh in all grid cells
    do ilev=1,pver
       do icol=1,ncol
          call qsat_water(t(icol,ilev), pmid(icol,ilev), tmp1, qs(icol,ilev))
          xrh(icol,ilev) = h2ommr(icol,ilev)/qs(icol,ilev)
          xrh(icol,ilev) = max(xrh(icol,ilev),0.0_r8)
          xrh(icol,ilev) = min(xrh(icol,ilev),1.0_r8)
          relhum(icol,ilev)=xrh(icol,ilev)
          xrh(icol,ilev)=min(xrh(icol,ilev),rhtab(10))
       end do
    end do

    !Find the relh-index in all grid-points
    do irelh=1,SIZE(rhtab) - 1
       do ilev=1,pver
          do icol=1,ncol
             if (xrh(icol,ilev) >= rhtab(irelh) .and. xrh(icol,ilev) <= rhtab(irelh+1)) then
                irh1(icol,ilev)=irelh                !lower index
                irh2(icol,ilev)=irelh+1              !higher index
             end if
          end do
       end do
    end do

    do ilev=1,pver
       do icol=1,ncol

          !Get the indexes out as floating point single numbers
          t_irh1 = irh1(icol,ilev)
          t_irh2 = irh2(icol,ilev)
          t_rh1  = rhtab(t_irh1)
          t_rh2  = rhtab(t_irh2)
          t_xrh  = xrh(icol,ilev)

          do imode = 0, nmodes_oslo
             !Do some weighting to mass mean property
             !weighting by 1.5 is number median ==> volumetric mean
             !http://dust.ess.uci.edu/facts/psd/psd.pdf
             rmeanvol = lifeCycleNumberMedianRadius(imode)*DEXP(1.5_r8*(log(lifeCycleSigma(imode)))**2)
             wetNumberMedianDiameter(icol,ilev,imode) =  0.1e-6_r8 !Initialize to something..
             mixed_dry_rho = 1.e3_r8

             tracerCounter = 0
             do itrac = 1,getNumberOfBackgroundTracersInMode(imode)

                tracerCounter = tracerCounter + 1

                !which tracer is this?
                mm = getTracerIndex(imode,itrac,.false.)

                !radius of lower rh-bin for this tracer
                rr1=rdivr0(t_irh1,mm)

                !radius of upper rh-bin for this tracer
                rr2=rdivr0(t_irh2,mm)

                !linear interpolate dry ==> wet radius for this tracer
                wetrad_tmp(tracerCounter) = (((t_rh2-t_xrh)*rr1+(t_xrh-t_rh1)*rr2)/ &
                     (t_rh2-t_rh1))*rmeanvol

                !mixed density of dry particle
                dry_rhopart_tmp(tracerCounter) = getDryDensity(imode,itrac)

             end do

             !Find the average growth of this mode
             !(still not taking into account how much we have!!)
             if(TracerCounter > 0)then

                !Convert to diameter and take average (note: This is MASS median diameter)
                wetNumberMedianDiameter(icol,ilev,imode) = 2.0_r8 * SUM(wetrad_tmp(1:tracerCounter))/dble(tracerCounter)

                !Take average density
                mixed_dry_rho = SUM(dry_rhopart_tmp(1:tracerCounter))/dble(tracerCounter)

                !At this point the radius is in "mass mean" space
                volumeFractionAerosol = MIN(1.0_r8, ( 2.0_r8*rmeanVol / wetNumberMedianDiameter(icol,ilev,imode) )**3)

                !wet density
                wetrho(icol,ilev,imode) = mixed_dry_rho * volumeFractionAerosol + (1._r8-volumeFractionAerosol)*rhoh2o

                !convert back to number median diameter (wet)
                wetNumberMedianDiameter(icol,ilev,imode) = &
                     wetNumberMedianDiameter(icol,ilev,imode)*DEXP(-1.5_r8*(log(lifeCycleSigma(imode)))**2)
             endif

          end do     !modes

          !Same thing for the process modes
          do itrac = 1,numberOfProcessModeTracers

             mm = tracerInProcessMode(itrac)   !process mode tracer (physics space)

             !weighting by 1.5 is number median ==> volumetric mean
             !http://dust.ess.uci.edu/facts/psd/psd.pdf
             rmeanvol = processModeNumberMedianRadius(itrac)*DEXP(1.5_r8*(log(processModeSigma(itrac)))**2)

             !radius of lower rh-bin for this tracer
             rr1=rdivr0(t_irh1,mm)

             !radius of upper rh-bin for this tracer
             rr2=rdivr0(t_irh2,mm)

             !Note this is MASS median diameter
             wetNumberMedianDiameter_processmode(icol,ilev,itrac) = (((t_rh2-t_xrh)*rr1+(t_xrh-t_rh1)*rr2)/ &
                  (t_rh2-t_rh1))*rmeanvol*2.0_r8

             volumeFractionAerosol = MIN(1.0, (2.0_r8*rmeanVol/wetnumberMedianDiameter_processmode(icol,ilev,itrac))**3)

             wetrho_processmode(icol,ilev,itrac) = volumeFractionAerosol*rhopart(mm) &
                  + (1.0_r8 - volumeFractionAerosol)*rhoh2o

             !convert back to number median diameter (wet)
             wetNumberMedianDiameter_processMode(icol,ilev,itrac) = &
                  wetNumberMedianDiameter_processMode(icol,ilev,itrac)*DEXP(-1.5_r8*(log(processModeSigma(itrac)))**2)
          end do !process modes
       end do !horizontal points
    end do !layers

  end subroutine calcaersize_sub

end module aero_model
