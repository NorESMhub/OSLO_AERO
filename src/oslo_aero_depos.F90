module oslo_aero_depos

  !------------------------------------------------------------------------------------------------
  ! Compute the contributions from oslo aero modal components of wet and dry
  ! deposition at the surface into the fields passed to the coupler.
  ! Wet deposition routines for both aerosols and gas phase constituents.
  !------------------------------------------------------------------------------------------------

  use shr_kind_mod,            only: r8 => shr_kind_r8
  use ppgrid,                  only: pcols, pver, pverp, begchunk, endchunk
  use constituents,            only: pcnst, cnst_name, cnst_get_ind
  use phys_control,            only: phys_getopts, cam_physpkg_is
  use phys_control,            only: history_aerosol_base,        &
                                     history_aerosol_decomposed,  &
                                     history_aerosol_debug_output
  use cam_abortutils,          only: endrun
  use cam_logfile,             only: iulog
  use camsrfexch,              only: cam_out_t
  use time_manager,            only: is_first_step
  use aerodep_flx,             only: aerodep_flx_prescribed
  use mo_drydep,               only: n_land_type, fraction_landuse
  use physics_types,           only: physics_ptend, physics_ptend_init
  use physics_buffer,          only: physics_buffer_desc, pbuf_get_chunk, pbuf_get_field, pbuf_get_index
  use physics_buffer,          only: pbuf_old_tim_idx
  use physconst,               only: gravit, rair, rhoh2o, boltz, pi, tmelt
  use cam_history,             only: outfld, fieldname_len, addfld, add_default, horiz_only
  use ref_pres,                only: top_lev => clim_modal_aero_top_lev
  !
  use oslo_aero_share,         only: nmodes, max_tracers_per_mode
  use oslo_aero_share,         only: numberOfProcessModeTracers, getNumberOfTracersInMode, getTracerIndex
  use oslo_aero_share,         only: is_process_mode, processModeMap, processModeSigma, lifeCycleSigma
  use oslo_aero_share,         only: belowCloudScavengingCoefficientProcessModes, belowCloudScavengingCoefficient
  use oslo_aero_share,         only: getCloudTracerIndex, GetCloudTracerIndexDirect, getCloudTracerName, qqcw_get_field
  use oslo_aero_share,         only: aerosol_type_name, N_AEROSOL_TYPES, aerosolType, AEROSOL_TYPE_SULFATE
  use oslo_aero_share,         only: l_bc_ax, l_bc_ni, l_bc_ai, l_bc_a, l_bc_ac
  use oslo_aero_share,         only: l_bc_n, l_om_ni, l_om_ai, l_om_ac, l_dst_a2, l_dst_a3
  use oslo_aero_share,         only: l_ss_a2, l_ss_a3, l_so4_a2
  use oslo_aero_dust_sediment, only: oslo_aero_dust_sediment_tend, oslo_aero_dust_sediment_vel

  implicit none
  private          ! Make default type private to the module

  ! Public interfaces
  public :: oslo_aero_depos_register
  public :: oslo_aero_depos_init
  public :: oslo_aero_depos_dry ! dry deposition
  public :: oslo_aero_depos_wet ! wet deposition
  public :: oslo_aero_wetdep_init

  ! Private interfaces
  private :: oslo_aero_depvel_part
  private :: oslo_set_srf_drydep
  private :: oslo_set_srf_wetdep
  private :: calcram
  private :: wetdep_inputs_t
  private :: wetdep_inputs_set
  private :: wetdepa_v2  ! scavenging codes for very soluble aerosols -- CAM5 version
  private :: wetdepg     ! scavenging of gas phase constituents by henry's law
  private :: clddiag     ! calc of cloudy volume and rain mixing ratio

  real(r8), public :: sol_facti_cloud_borne

  real(r8), parameter :: cmftau = 3600._r8
  real(r8), parameter :: molwta = 28.97_r8 ! molecular weight dry air gm/mole

  type wetdep_inputs_t
     real(r8), pointer :: cldt(:,:)  => null()  ! cloud fraction
     real(r8), pointer :: qme(:,:)   => null()
     real(r8), pointer :: prain(:,:) => null()
     real(r8), pointer :: evapr(:,:) => null()
     real(r8) :: cldcu(pcols,pver)     ! convective cloud fraction, currently empty
     real(r8) :: evapc(pcols,pver)     ! Evaporation rate of convective precipitation
     real(r8) :: cmfdqr(pcols,pver)    ! convective production of rain
     real(r8) :: conicw(pcols,pver)    ! convective in-cloud water
     real(r8) :: totcond(pcols, pver)  ! total condensate
     real(r8) :: cldv(pcols,pver)      ! cloudy volume undergoing wet chem and scavenging
     real(r8) :: cldvcu(pcols,pver)    ! Convective precipitation area at the top interface of current layer
     real(r8) :: cldvst(pcols,pver)    ! Stratiform precipitation area at the top interface of current layer
  end type wetdep_inputs_t

  logical :: convproc_do_aer = .FALSE.
  logical :: drydep_lq(pcnst)
  logical :: wetdep_lq(pcnst)

  integer :: fracis_idx      = 0
  integer :: prain_idx       = 0
  integer :: cld_idx         = 0
  integer :: qme_idx         = 0
  integer :: nevapr_idx      = 0
  integer :: icwmrdp_idx     = 0
  integer :: icwmrsh_idx     = 0
  integer :: rprddp_idx      = 0
  integer :: rprdsh_idx      = 0
  integer :: sh_frac_idx     = 0
  integer :: dp_frac_idx     = 0
  integer :: nevapr_shcu_idx = 0
  integer :: nevapr_dpcu_idx = 0
  integer :: ixcldice, ixcldliq

  integer :: idx_wd_a_h2so4 = -1
  integer :: tracer_index(0:nmodes, max_tracers_per_mode)
  integer :: num_tracers_in_mode(0:nmodes)

!===============================================================================
contains
!===============================================================================

   subroutine oslo_aero_depos_register()
      use physics_buffer,  only: pbuf_add_field, dtype_r8
      use ppgrid,          only: pcols

      ! Register a pbuf field for
      call pbuf_add_field('WD_A_H2SO4', 'global', dtype_r8, (/pcols/), idx_wd_a_h2so4)
   end subroutine oslo_aero_depos_register

  subroutine oslo_aero_depos_init( pbuf2d )
    use time_manager,   only: is_first_step
    use physics_buffer, only: pbuf_set_field

    ! Set oslo aeroslo deposition history output

    ! arguments
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! local variables
    integer            :: imode, itrac
    integer            :: lchnk
    integer            :: l_atype
    integer            :: tracerIndex
    integer            :: astat, id
    real(r8), pointer  :: qqcw(:,:)
    logical            :: history_aerosol    ! Output the aerosol tendencies
    character(len=2)   :: unit_basename='kg' ! Units 'kg' or '1'
    character(len=100) :: aName              ! tracer name
    logical            :: is_in_output(pcnst)
    !-----------------------------------------------------------------------

    fracis_idx      = pbuf_get_index('FRACIS')
    prain_idx       = pbuf_get_index('PRAIN')
    nevapr_shcu_idx = pbuf_get_index('NEVAPR_SHCU')

    call phys_getopts( history_aerosol_out = history_aerosol )
    if (is_first_step()) then
       call pbuf_set_field(pbuf2d, idx_wd_a_h2so4, 0.0_r8)
    end if

    is_in_output(:) =.false.
    drydep_lq(:) =.false.
    wetdep_lq(:) =.false.

    do imode = 0,nmodes
       num_tracers_in_mode(imode) = getNumberOfTracersInMode(imode)
       do itrac = 1,num_tracers_in_mode(imode)
          tracer_index(imode,itrac) = getTracerIndex(imode,itrac,.false.)
       end do
    end do

    ! Mode 0 is not subject to wet deposition? (check noresm1 code..)
    do imode=0,nmodes
       do itrac=1,num_tracers_in_mode(imode)

          tracerIndex = tracer_index(imode,itrac)
          drydep_lq(tracerIndex)=.true.
          wetdep_lq(tracerIndex)=.true.

          if(is_in_output(tracerIndex))then
             cycle
          endif

          aName = cnst_name(tracerIndex)

          call addfld (trim(aName)//'SFWET',horiz_only, 'A', unit_basename//'/m2/s',  &
               'Wet deposition flux at surface')
          call addfld (trim(aName)//'SFSIC',horiz_only, 'A', unit_basename//'/m2/s ', &
               'Wet deposition flux (incloud, convective) at surface')
          call addfld (trim(aName)//'SFSIS',horiz_only, 'A', unit_basename//'/m2/s ', &
               'Wet deposition flux (incloud, stratiform) at surface')
          call addfld (trim(aName)//'SFSBC',horiz_only, 'A', unit_basename//'/m2/s ', &
               'Wet deposition flux (belowcloud, convective) at surface')
          call addfld (trim(aName)//'SFSBS',horiz_only, 'A', unit_basename//'/m2/s ', &
               'Wet deposition flux (belowcloud, stratiform) at surface')
          call addfld (trim(aName)//'WET',(/'lev'/), 'A', unit_basename//'/kg/s ',&
               'wet deposition tendency')
          call addfld (trim(aName)//'SIC',(/'lev'/), 'A', unit_basename//'/kg/s ', &
               trim(aName)//' ic wet deposition')
          call addfld (trim(aName)//'SIS',(/'lev'/), 'A', unit_basename//'/kg/s ', &
               trim(aName)//' is wet deposition')
          call addfld (trim(aName)//'SBC',(/'lev'/), 'A', unit_basename//'/kg/s ', &
               trim(aName)//' bc wet deposition')
          call addfld (trim(aName)//'SBS',(/'lev'/), 'A', unit_basename//'/kg/s ', &
               trim(aName)//' bs wet deposition')

          ! Extra wd ouptut
          if ( history_aerosol ) then
             call add_default (trim(aName)//'SFSIC', 1, ' ')
             call add_default (trim(aName)//'SFSIS', 1, ' ')
             call add_default (trim(aName)//'SFSBC', 1, ' ')
             call add_default (trim(aName)//'SFSBS', 1, ' ')
          endif

         if ( history_aerosol .or. history_aerosol_decomposed ) then
            call add_default (trim(aName)//'SFWET', 1, ' ')
         endif

          ! Dry deposition fluxes and velocity
          call addfld (trim(aName)//'DDF',horiz_only, 'A', unit_basename//'/m2/s ', &
               trim(aName)//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(aName)//'TBF',horiz_only, 'A' ,unit_basename//'/m2/s', &
               trim(aName)//' turbulent dry deposition flux')
          call addfld (trim(aName)//'GVF',horiz_only, 'A', unit_basename//'/m2/s ', &
               trim(aName)//' gravitational dry deposition flux')
          call addfld (trim(aName)//'DTQ',(/'lev'/), 'A', unit_basename//'/kg/s ', &
               trim(aName)//' dry deposition')
          call addfld (trim(aName)//'DDV',(/'lev'/), 'A', 'm/s', &
               trim(aName)//' deposition velocity')

          ! extra drydep output
          if ( history_aerosol ) then
             call add_default (trim(aName)//'TBF', 1, ' ')
             call add_default (trim(aName)//'GVF', 1, ' ')
             !call add_default (trim(aName)//'DDV', 1, ' ')
          endif

          if ( history_aerosol .or. history_aerosol_decomposed ) then
            call add_default (trim(aName)//'DDF', 1, ' ')
          endif

          ! some tracers are not in cloud water
          if(getCloudTracerIndexDirect(tracerIndex) < 0)then
             cycle
          endif

          aName = trim(getCloudTracerName(tracerIndex))

          ! Cloud water fields (from mo_chm_diags.F90)
          call addfld (trim(aName)//'SFWET', horiz_only, 'A', unit_basename//'/m2/s', &
               trim(aName)//' wet deposition flux at surface')
          call addfld (trim(aName)//'SFSIC', horiz_only, 'A',unit_basename//'/m2/s ', &
               trim(aName)//' wet deposition flux (incloud, convective) at surface')
          call addfld (trim(aName)//'SFSIS', horiz_only, 'A', unit_basename//'/m2/s ', &
               trim(aName)//' wet deposition flux (incloud, stratiform) at surface')
          call addfld (trim(aName)//'SFSBC', horiz_only, 'A', unit_basename//'/m2/s ' , &
               trim(aName)//' wet deposition flux (belowcloud, convective) at surface')
          call addfld (trim(aName)//'SFSBS', horiz_only, 'A', unit_basename//'/m2/s ' , &
               trim(aName)//' wet deposition flux (belowcloud, stratiform) at surface')

          ! dry deposition
          call addfld (trim(aName)//'DDF',   horiz_only, 'A', unit_basename//'/m2/s ',  &
               trim(aName)//' dry deposition flux at bottom (grav + turb)')
          call addfld (trim(aName)//'TBF',   horiz_only, 'A', unit_basename//'/m2/s ',  &
               trim(aName)//' turbulent dry deposition flux')
          call addfld (trim(aName)//'GVF',   horiz_only, 'A', unit_basename//'/m2/s ',  &
               trim(aName)//' gravitational dry deposition flux')

          is_in_output(tracerIndex) = .true.

       end do !tracers
    enddo    !modes

    do l_atype=1,N_AEROSOL_TYPES

      call addfld( 'dry_'//trim(aerosol_type_name(l_atype)), horiz_only, 'A', unit_basename//'/m2/s ',  &
         trim(aerosol_type_name(l_atype))//' dry deposition flux at bottom (grav + turb)')
      call addfld('wet_'//trim(aerosol_type_name(l_atype)), horiz_only, 'A', unit_basename//'/m2/s', &
         trim(aerosol_type_name(l_atype))//' wet deposition flux at surface')

      if ( l_atype == AEROSOL_TYPE_SULFATE ) then
         call addfld( 'dry_'//trim(aerosol_type_name(l_atype))//'_S', horiz_only, 'A', unit_basename//'*S/m2/s ',  &
            trim(aerosol_type_name(l_atype))//' dry deposition flux at bottom (grav + turb), sulfur mass only')
         call addfld('wet_'//trim(aerosol_type_name(l_atype))//'_S', horiz_only, 'A', unit_basename//'*S/m2/s', &
            trim(aerosol_type_name(l_atype))//' wet deposition flux at surface, sulfur mass only')
         call addfld('wd_a_h2so4_debug', horiz_only, 'A', unit_basename//'*S/m2/s', &
            'wd_a_h2so4_debug - used to debug in-model summation of sulfur mass in aerosols')
      endif

      ! we require history_aerosol_base flag to add the fields to default output
      if ( history_aerosol_base ) then

         call add_default('dry_'//trim(aerosol_type_name(l_atype)), 1, ' ')
         call add_default('wet_'//trim(aerosol_type_name(l_atype)), 1, ' ')

         if ( l_atype == AEROSOL_TYPE_SULFATE ) then
            call add_default('dry_'//trim(aerosol_type_name(l_atype))//'_S', 1, ' ')
            call add_default('wet_'//trim(aerosol_type_name(l_atype))//'_S', 1, ' ')
         endif

      endif
      if ( history_aerosol_debug_output .and. l_atype == AEROSOL_TYPE_SULFATE ) then
         call add_default('wd_a_h2so4_debug', 1, ' ')
      endif
   end do

    !initialize cloud concentrations (initialize cloud bourne constituents in physics buffer)
    if (is_first_step()) then
       do itrac = 1, pcnst
          do lchnk = begchunk, endchunk
             qqcw => qqcw_get_field(pbuf_get_chunk(pbuf2d,lchnk), itrac)
             if (associated(qqcw)) then
                qqcw = 1.e-38_r8
             end if
          end do
       end do
    end if

  end subroutine oslo_aero_depos_init

  !===============================================================================
  subroutine oslo_aero_depos_dry  (lchnk, ncol, psetcols, &
       t, pmid, pdel, pint, q, &
       landfrac, icefrac, ocnfrac, fvin, ram1in, cflx, &
       pbuf, obklen, ustar, dt, &
       dgncur_awet, wetdens, dgncur_awet_processmode, wetdens_processmode, &
       cam_out, ptend)

   ! use
   use oslo_aero_share, only          : sulfurMassFraction

    ! Arguments:
    integer  ,           intent(in)    :: lchnk
    integer  ,           intent(in)    :: ncol
    integer  ,           intent(in)    :: psetcols
    real(r8) ,           intent(in)    :: t(pcols,pver)     ! Model level temperatures (K)
    real(r8) ,           intent(in)    :: q(pcols,pver,pcnst)
    real(r8) ,           intent(in)    :: pmid(pcols,pver)  ! Model level pressures (Pa)
    real(r8) ,           intent(in)    :: pdel(pcols,pver)  ! model layer thickness (Pa)
    real(r8) ,           intent(in)    :: pint(pcols,pverp) ! Model interface pressures (10*Pa)
    real(r8) ,           intent(in)    :: landfrac(pcols)   ! land fraction
    real(r8) ,           intent(in)    :: icefrac(pcols)    ! ice fraction
    real(r8) ,           intent(in)    :: ocnfrac(pcols)    ! ocean fraction
    real(r8) ,           intent(in)    :: fvin(pcols)       !
    real(r8) ,           intent(in)    :: ram1in(pcols)     ! for dry dep velocities from land model for progseasalts
    real(r8) ,           intent(in)    :: cflx(pcols,pcnst) ! Surface fluxes
    type(physics_buffer_desc), pointer :: pbuf(:)
    real(r8),            intent(in)    :: obklen(:)
    real(r8),            intent(in)    :: ustar(:)          ! sfc fric vel
    real(r8),            intent(in)    :: dt                ! time step
    real(r8),            intent(in)    :: dgncur_awet(pcols,pver,0:nmodes)
    real(r8),            intent(in)    :: wetdens(pcols,pver,0:nmodes)
    real(r8),            intent(in)    :: dgncur_awet_processmode(pcols, pver, numberOfProcessModeTracers)
    real(r8),            intent(in)    :: wetdens_processmode(pcols, pver, numberOfProcessModeTracers)
    type(cam_out_t),     intent(inout) :: cam_out           ! export state
    type(physics_ptend), intent(out)   :: ptend             ! indivdual parameterization tendencies

    ! local vars
    real(r8) :: fv(pcols)                ! for dry dep velocities, from land modified over ocean & ice
    real(r8) :: ram1(pcols)              ! for dry dep velocities, from land modified over ocean & ice

    integer :: jvlc                      ! index for last dimension of vlc_xxx arrays
    integer :: lphase                    ! index for interstitial / cloudborne aerosol
    integer :: lspec                     ! index for aerosol number / chem-mass / water-mass
    integer :: imode                     ! aerosol mode index
    integer :: itrac                        ! tracer index
    integer :: icol
    integer :: l_atype                   ! aerosol type index

    real(r8) :: tvs(pcols,pver)
    real(r8) :: rho(pcols,pver)          ! air density in kg/m3
    real(r8) :: sflx(pcols)              ! deposition flux
    real(r8) :: sflx_DDF_arosol_type(pcols, N_AEROSOL_TYPES) ! deposition flux for the aerosol types, no weighted sums
    real(r8) :: sflx_DDF_SULFATE_S(pcols) ! deposition flux for the sulfate aerosol type, sulfur mass only
    real(r8)::  dep_trb(pcols)           ! kg/m2/s
    real(r8)::  dep_grv(pcols)           ! kg/m2/s (total of grav and trb)
    real(r8) :: pvmzaer(pcols,pverp)     ! sedimentation velocity in Pa
    real(r8) :: dqdt_tmp(pcols,pver)     ! temporary array to hold tendency for 1 species

    real(r8) :: rad_drop(pcols,pver)
    real(r8) :: dens_drop(pcols,pver)
    real(r8) :: sg_drop(pcols,pver)
    real(r8) :: rad_aer(pcols,pver)
    real(r8) :: dens_aer(pcols,pver)
    real(r8) :: sg_aer(pcols,pver)

    real(r8) :: vlc_dry(pcols,pver,4)    ! dep velocity
    real(r8) :: vlc_grv(pcols,pver,4)    ! dep velocity
    real(r8)::  vlc_trb(pcols,4)         ! dep velocity
    real(r8) :: aerdepdryis(pcols,pcnst) ! aerosol dry deposition (interstitial)
    real(r8) :: aerdepdrycw(pcols,pcnst) ! aerosol dry deposition (cloud water)
    real(r8), pointer :: fldcw(:,:)

    ! oslo aerosols
    real(r8) :: logSigma
    logical  :: is_done(pcnst,2)
    !-----------------------------------------------------------------------

    aerdepdryis(:,:)=0._r8
    aerdepdrycw(:,:)=0._r8
    sflx_DDF_arosol_type(:,:) = 0._r8
    sflx_DDF_SULFATE_S(:) = 0._r8

    ! calc ram and fv over ocean and sea ice ...
    call calcram( ncol,landfrac, icefrac, ocnfrac, obklen, &
         ustar, ram1in, ram1, t(:,pver), pmid(:,pver), pdel(:,pver), fvin, fv)

    call outfld( 'airFV', fv(:ncol), ncol, lchnk )
    call outfld( 'RAM1', ram1(:ncol), ncol, lchnk )

    ! note that tendencies are not only in sfc layer (because of sedimentation)
    ! and that ptend is updated within each subroutine for different species

    call physics_ptend_init(ptend, psetcols, 'aero_model_drydep', lq=drydep_lq)

    tvs(:ncol,:) = t(:ncol,:) !*(1+q(:ncol,ilev)
    rho(:ncol,:)=  pmid(:ncol,:)/(rair*t(:ncol,:))
    is_done(:,:) = .false.

    ! calc settling/deposition velocities for cloud droplets (and cloud-borne aerosols)
    ! *** mean drop radius should eventually be computed from ndrop and qcldwtr
    rad_drop(:,:) = 5.0e-6_r8
    dens_drop(:,:) = rhoh2o
    sg_drop(:,:) = 1.46_r8

    !jvlc = 3
    jvlc = 4
    call oslo_aero_depvel_part( ncol, t(:,:), pmid(:,:), ram1, fv,  &
         vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
         rad_drop(:,:), dens_drop(:,:), sg_drop(:,:), 3, lchnk)

    !At this point we really need to distribute the lifecycle-tracers over
    !the actual modes (maybe according to surface available of background tracers?)

    !in mam3, jvlc = 1 means number-concentration
    !in oslo_aero, jvlc = 1 means process-modes
    !The following logic is based on that process-mode tracers
    !always follow AFTER the actual tracers!!

    dens_aer(:,:) = 0._r8
    do imode = 0, nmodes   ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial / cloud-borne forms

          if (lphase == 1) then   ! interstial aerosol - calc settling/dep velocities of mode

             logSigma = log(lifeCycleSigma(imode))

             ! rad_aer = volume mean wet radius (imode)
             ! dgncur_awet = geometric mean wet diameter for number distribution (imode)
             if(top_lev > 1) then
                rad_aer(1:ncol,:top_lev-1) = 0._r8
             end if
             rad_aer(1:ncol,top_lev:) = 0.5_r8*dgncur_awet(1:ncol,top_lev:,imode)*exp(1.5_r8*(logSigma**2))

             ! dens_aer(1:ncol,:) = wet density (kg/m3)
             if(top_lev>1)then
                dens_aer(1:ncol,:top_lev-1) = 0._r8
             end if
             dens_aer(1:ncol,top_lev:) = wetdens(1:ncol,top_lev:,imode)

             sg_aer(1:ncol,:) = lifecycleSigma(imode)

             jvlc = 2
             call oslo_aero_depvel_part( ncol, t(:,:), pmid(:,:), ram1, fv,  &
                  vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                  rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, lchnk)
          end if

          do lspec = 1, num_tracers_in_mode(imode)   ! loop over number + constituents

             itrac = tracer_index(imode,lspec)
             if(is_done(itrac,lphase)) then
                cycle
             endif
             is_done(itrac,lphase)=.true.

             ! Calculate sediment velocity

             if (lphase == 1) then
                jvlc = 2 ! mass in clean air tracers
                !Process tracers have their own velocity based on fixed size / density
                !Calculate the velocity to use for this specie..
                if ( is_process_mode(itrac, .false.) ) then
                   jvlc = 1
                   logSigma = log(processModeSigma(processModeMap(itrac)))
                   if(top_lev>1)then
                      rad_aer(1:ncol, top_lev-1) = 0.0_r8
                   end if
                   rad_aer(1:ncol,top_lev:) = &
                        0.5_r8*dgncur_awet_processmode(1:ncol,top_lev:,processModeMap(itrac))*exp(1.5_r8*(logSigma**2))

                   call oslo_aero_depvel_part( ncol, t(:,:), pmid(:,:), ram1, fv,  &
                        vlc_dry(:,:,jvlc), vlc_trb(:,jvlc), vlc_grv(:,:,jvlc),  &
                        rad_aer(:,:), dens_aer(:,:), sg_aer(:,:), 3, lchnk)
                endif
             else
                jvlc = 4              !mass in cloud tracers
             endif

             if (itrac <= 0) cycle

             ! Calculate sediment tendency

             if ((lphase == 1) .and. (lspec <= num_tracers_in_mode(imode))) then
                ptend%lq(itrac) = .TRUE.

                ! use pvprogseasalts instead (means making the top level 0)
                pvmzaer(:ncol,1)=0._r8
                pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

                call outfld( trim(cnst_name(itrac))//'DDV', pvmzaer(:ncol,2:pverp), ncol, lchnk )

                ! use phil's method
                ! convert from meters/sec to pascals/sec, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

                ! calculate the tendencies and sfc fluxes from the above velocities
                call oslo_aero_dust_sediment_tend(ncol, dt, pint(:,:), pmid, pdel, t , &
                     q(:,:,itrac),  pvmzaer,  ptend%q(:,:,itrac), sflx)

                ! apportion dry deposition into turb and gravitational settling for tapes
                dep_trb = 0._r8
                dep_grv = 0._r8
                do icol=1,ncol
                   if (vlc_dry(icol,pver,jvlc) /= 0._r8) then
                      dep_trb(icol)=sflx(icol)*vlc_trb(icol,jvlc)/vlc_dry(icol,pver,jvlc)
                      dep_grv(icol)=sflx(icol)*vlc_grv(icol,pver,jvlc)/vlc_dry(icol,pver,jvlc)
                   endif
                enddo

                call outfld( trim(cnst_name(itrac))//'DDF', sflx(:ncol), ncol, lchnk)
                call outfld( trim(cnst_name(itrac))//'TBF', dep_trb(:ncol), ncol, lchnk )
                call outfld( trim(cnst_name(itrac))//'GVF', dep_grv(:ncol), ncol, lchnk )
                call outfld( trim(cnst_name(itrac))//'DTQ', ptend%q(:ncol,:,itrac), ncol, lchnk)
                aerdepdryis(:ncol,itrac) = sflx(:ncol)

             else  ! lphase == 2

                !Pick up the cloud tracers (oslo)
                fldcw => qqcw_get_field(pbuf, itrac)
                if( .not. associated(fldcw))then
                   cycle
                end if

                ! use pvprogseasalts instead (means making the top level 0)
                pvmzaer(:ncol,1)=0._r8
                pvmzaer(:ncol,2:pverp) = vlc_dry(:ncol,:,jvlc)

                ! Hardwire the method from Phil
                ! convert from meters/sec to pascals/sec
                ! pvprogseasalts(:,1) is assumed zero, use density from layer above in conversion
                pvmzaer(:ncol,2:pverp) = pvmzaer(:ncol,2:pverp) * rho(:ncol,:)*gravit

                ! calculate the tendencies and sfc fluxes from the above velocities
                call oslo_aero_dust_sediment_tend(ncol, dt, pint(:,:), pmid, pdel, t, &
                     fldcw(:,:), pvmzaer, dqdt_tmp(:,:), sflx)

                ! apportion dry deposition into turb and gravitational settling for tapes
                dep_trb = 0._r8
                dep_grv = 0._r8
                do icol=1,ncol
                   if (vlc_dry(icol,pver,jvlc) /= 0._r8) then
                      dep_trb(icol)=sflx(icol)*vlc_trb(icol,jvlc)/vlc_dry(icol,pver,jvlc)
                      dep_grv(icol)=sflx(icol)*vlc_grv(icol,pver,jvlc)/vlc_dry(icol,pver,jvlc)
                   end if
                enddo

                fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                call outfld( trim(getCloudTracerName(itrac))//'DDF', sflx(:ncol), ncol, lchnk)
                call outfld( trim(getCloudTracerName(itrac))//'TBF', dep_trb(:ncol), ncol, lchnk )
                call outfld( trim(getCloudTracerName(itrac))//'GVF', dep_grv(:ncol), ncol, lchnk )
                aerdepdrycw(:ncol,itrac) = sflx(:ncol)

             endif

            ! accumulate the deposition flux for the aerosol type
            ! All will have a version without weighted sum, that is ...DDF
            ! sulfate will have a version with weighted sum, that is ...SDDF
            sflx_DDF_arosol_type(:ncol, aerosolType(itrac)) = sflx_DDF_arosol_type(:ncol, aerosolType(itrac)) + sflx(:ncol)
            if ( aerosolType(itrac) ==  AEROSOL_TYPE_SULFATE ) then
               sflx_DDF_SULFATE_S(:ncol) = sflx_DDF_SULFATE_S(:ncol) + ( sflx(:ncol) * sulfurMassFraction(itrac) )
            endif

          enddo   ! lspec = 0, nspec_amode(imode)+1
       enddo   ! lphase = 1, 2
    enddo   ! imode = 1, ntot_amode

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols

    if (.not. aerodep_flx_prescribed()) then
       call  oslo_set_srf_drydep(ncol, aerdepdryis, aerdepdrycw, &
            cam_out%bcphidry, cam_out%bcphodry, cam_out%ocphidry, cam_out%ocphodry, &
            cam_out%dstdry1, cam_out%dstdry2, cam_out%dstdry3, cam_out%dstdry4)
    endif

   do l_atype=1,N_AEROSOL_TYPES
      ! add the dry deposition rate of the compound aerosols to output
      call outfld('dry_'//trim(aerosol_type_name(l_atype)), sflx_DDF_arosol_type(:ncol,l_atype), ncol, lchnk)
      if ( l_atype == AEROSOL_TYPE_SULFATE ) then
         call outfld('dry_'//trim(aerosol_type_name(l_atype))//'_S', sflx_DDF_SULFATE_S(:ncol), ncol, lchnk)
      endif
   end do

  end subroutine oslo_aero_depos_dry

  !===============================================================================
  subroutine oslo_aero_depos_wet ( lchnk, ncol, psetcols, pmid, pdel, q, t, &
       dt, dlf, cam_out, ptend, pbuf)

   ! use
   use oslo_aero_share, only          : sulfurMassFraction, l_h2so4

    integer ,            intent(in)    :: lchnk            ! chunk identifier
    integer ,            intent(in)    :: ncol             ! number of atmospheri columns
    integer ,            intent(in)    :: psetcols
    real(r8) ,           intent(in)    :: pmid(pcols,pver) ! Model level pressures (Pa)
    real(r8) ,           intent(in)    :: pdel(pcols,pver) ! model layer thickness (Pa)
    real(r8) ,           intent(in)    :: q(pcols,pver,pcnst)
    real(r8) ,           intent(in)    :: t(pcols,pver)    ! Model level temperatures (K)
    real(r8),            intent(in)    :: dt               ! time step
    real(r8),            intent(in)    :: dlf(:,:)         ! shallow+deep convective detrainment [kg/kg/s]
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(cam_out_t),     intent(inout) :: cam_out          ! export state
    type(physics_ptend), intent(out)   :: ptend            ! indivdual parameterization tendencies

    ! Local variables
    integer  :: imode                             ! tracer index
    integer  :: l_atype                           ! aerosol type index
    integer  :: icol,ilev,itrac
    real(r8) :: iscavt(pcols, pver)
    real(r8) :: icscavt(pcols, pver)
    real(r8) :: isscavt(pcols, pver)
    real(r8) :: bcscavt(pcols, pver)
    real(r8) :: bsscavt(pcols, pver)
    real(r8) :: sol_factb, sol_facti
    real(r8) :: sol_factic(pcols,pver)
    real(r8) :: sflx(pcols)                   ! deposition flux
    real(r8) :: sflx_SFWET_arosol_type(pcols, N_AEROSOL_TYPES) ! deposition flux for the aerosol types, aerosol mass
    real(r8) :: sflx_SFWET_SULFATE_S(pcols)   ! deposition flux for sulfate, sulfur mass only
    real(r8) :: scavcoef(pcols,pver)          ! Dana and Hales coefficient (/mm) (0.1)
    integer  :: jnv                           ! index for scavcoefnv 3rd dimension
    integer  :: lphase                        ! index for interstitial / cloudborne aerosol
    integer  :: lspec                         ! index for aerosol number / chem-mass / water-mass
    real(r8) :: dqdt_tmp(pcols,pver)          ! temporary array to hold tendency for 1 species
    real(r8) :: f_act_conv(pcols,pver)        ! prescribed aerosol activation fraction for convective cloud   ! rce 2010/05/01
    real(r8) :: f_act_conv_coarse(pcols,pver) ! similar but for coarse mode                                   ! rce 2010/05/02
    real(r8) :: f_act_conv_coarse_dust, f_act_conv_coarse_nacl ! rce 2010/05/02
    real(r8) :: fracis_cw(pcols,pver)
    real(r8) :: prec(pcols)                    ! precipitation rate
    real(r8) :: q_tmp(pcols,pver)              ! temporary array to hold "most current" mixing ratio for 1 species
    real(r8) :: qqcw_tmp(pcols,pver)           ! temporary array to hold qqcw   ! rce 2010/05/01
    real(r8) :: scavcoefnv(pcols,pver,0:2)     ! Dana and Hales coefficient (/mm) for
    real(r8) :: water_old, water_new           ! temporary old/new aerosol water mix-rat
    logical  :: isprx(pcols,pver)              ! true if precipation
    real(r8) :: aerdepwetis(pcols,pcnst)       ! aerosol wet deposition (interstitial)
    real(r8) :: aerdepwetcw(pcols,pcnst)       ! aerosol wet deposition (cloud water)
    real(r8) :: sflxic(pcols) ! deposition flux
    real(r8) :: sflxbc(pcols) ! deposition flux
    real(r8) :: rcscavt(pcols, pver)
    real(r8) :: rsscavt(pcols, pver)
    real(r8) :: qqcw_in(pcols,pver), qqcw_sav(pcols,pver,pcnst) ! temporary array to hold qqcw for the current mode
    logical  :: is_done(pcnst,2)
    real(r8) :: zeroAerosolConcentration(pcols,pver)
    real(r8), pointer :: fldcw(:,:)
    real(r8), pointer :: fracis(:,:,:)   ! fraction of transported species that are insoluble
    real(r8), pointer :: wd_a_h2so4(:)
    type(wetdep_inputs_t) :: dep_inputs
    !-----------------------------------------------------------------------

    call physics_ptend_init(ptend, psetcols, 'aero_model_wetdep', lq=wetdep_lq)

    is_done(:,:) = .false.

    zeroAerosolConcentration(:,:)=0.0_r8
    sflx_SFWET_arosol_type(:,:) = 0._r8
    sflx_SFWET_SULFATE_S(:) = 0._r8

    ! Wet deposition of mozart aerosol species.
    ptend%name  = ptend%name//'+mz_aero_wetdep'

    call wetdep_inputs_set(ncol, pmid, pdel, t, q, pbuf, dep_inputs)
    call pbuf_get_field(pbuf, fracis_idx, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/) )

    prec(:ncol)=0._r8
    do ilev=1,pver
       where (prec(:ncol) >= 1.e-7_r8)
          isprx(:ncol,ilev) = .true.
       elsewhere
          isprx(:ncol,ilev) = .false.
       endwhere
       prec(:ncol) = prec(:ncol) + &
            (dep_inputs%prain(:ncol,ilev) + dep_inputs%cmfdqr(:ncol,ilev) - dep_inputs%evapr(:ncol,ilev)) * pdel(:ncol,ilev)/gravit
    end do


    ! calculate the mass-weighted sol_factic for coarse mode species
    !    sol_factic_coarse(:,:) = 0.30_r8   ! tuned 1/4
    f_act_conv_coarse(:,:) = 0.60_r8   ! rce 2010/05/02
    f_act_conv_coarse_dust = 0.40_r8   ! rce 2010/05/02
    f_act_conv_coarse_nacl = 0.80_r8   ! rce 2010/05/02
    f_act_conv_coarse(:,:) = 0.5_r8

    scavcoefnv(:,:,0) = 0.0_r8   ! below-cloud scavcoef = 0.0 for cloud-borne species

    do imode = 0, nmodes  ! main loop over aerosol modes

       do lphase = 1, 2   ! loop over interstitial (1) and cloud-borne (2) forms

          ! sol_factb and sol_facti values
          ! sol_factb - currently this is basically a tuning factor
          ! sol_facti & sol_factic - currently has a physical basis, and reflects activation fraction
          !
          ! 2008-mar-07 rce - sol_factb (interstitial) changed from 0.3 to 0.1
          ! - sol_factic (interstitial, dust modes) changed from 1.0 to 0.5
          ! - sol_factic (cloud-borne, pcarb modes) no need to set it to 0.0
          ! because the cloud-borne pcarbon == 0 (no activation)
          !
          ! rce 2010/05/02
          ! prior to this date, sol_factic was used for convective in-cloud wet removal,
          ! and its value reflected a combination of an activation fraction (which varied between modes)
          ! and a tuning factor
          ! from this date forward, two parameters are used for convective in-cloud wet removal
          ! f_act_conv is the activation fraction
          ! note that "non-activation" of aerosol in air entrained into updrafts should
          ! be included here
          ! eventually we might use the activate routine (with w ~= 1 m/s) to calculate
          ! this, but there is still the entrainment issue
          ! sol_factic is strictly a tuning factor
          !
          if (lphase == 1) then   ! interstial aerosol

             scavcoefnv(:,:,1) = 0.1_r8  !Used by MAM for number concentration

             sol_factb  = 0.1_r8   ! all below-cloud scav ON (0.1 "tuning factor")
             ! sol_factb  = 0.03_r8   ! all below-cloud scav ON (0.1 "tuning factor")  ! tuned 1/6

             sol_facti  = 0.0_r8   ! strat  in-cloud scav totally OFF for institial

             sol_factic = 0.4_r8      ! xl 2010/05/20

             !fxm: simplified relative to MAM
             f_act_conv = 0.8 !ag: Introduce tuning per component later
          else   ! cloud-borne aerosol (borne by stratiform cloud drops)
             !default 100 % is scavenged by cloud -borne
             sol_facti_cloud_borne = 1.0_r8

             sol_factb  = 0.0_r8                ! all below-cloud scav OFF (anything cloud-borne is located "in-cloud")
             sol_facti  = sol_facti_cloud_borne ! strat  in-cloud scav cloud-borne tuning factor
             sol_factic = 0.0_r8                ! conv   in-cloud scav OFF (having this on would mean
             !        that conv precip collects strat droplets)
             f_act_conv = 0.0_r8                ! conv   in-cloud scav OFF (having this on would mean
          end if

          if (convproc_do_aer .and. lphase == 1) then
             ! if modal aero convproc is turned on for aerosols, then
             !    turn off the convective in-cloud removal for interstitial aerosols
             !    (but leave the below-cloud on, as convproc only does in-cloud)
             !    and turn off the outfld SFWET, SFSIC, SFSID, SFSEC, and SFSED calls
             ! for (stratiform)-cloudborne aerosols, convective wet removal
             !    (all forms) is zero, so no action is needed
             sol_factic = 0.0_r8
          endif

          do lspec = 1,num_tracers_in_mode(imode)   ! loop over number + chem constituents + water
             itrac = tracer_index(imode,lspec)
             if(is_done(itrac,lphase)) then
                cycle
             endif
             is_done(itrac,lphase)=.true.

             if (lphase == 1) then
                jnv = 2
                !Set correct below cloud scaveing coefficients
                !Hard-coded values per mode in NorESM
                if(is_process_mode(itrac,.FALSE.))then
                   scavcoefnv(:,:,jnv) = belowCloudScavengingCoefficientProcessModes(processModeMap(itrac))
                else
                   scavcoefnv(:,:,jnv) = belowCloudScavengingCoefficient(imode)
                end if
             else
                jnv = 0  !==> below cloud scavenging coefficients are zero (see above)
             endif

             ! Increase scavenging efficiency for large soluble particles.
             if ((lphase==1).and.((itrac==l_ss_a2).or.(itrac==l_ss_a3).or.(itrac==l_so4_a2))) then
                 sol_factic=1.0_r8
                 f_act_conv=1.0_r8
             end if

             if ((lphase == 1) .and. (lspec <= num_tracers_in_mode(imode))) then
                ptend%lq(itrac) = .TRUE.
                dqdt_tmp(:,:) = 0.0_r8
                ! q_tmp reflects changes from modal_aero_calcsize and is the "most current" q
                q_tmp(1:ncol,:) = q(1:ncol,:,itrac) + ptend%q(1:ncol,:,itrac)*dt
                if(convproc_do_aer) then
                   !Feed in the saved cloudborne mixing ratios from phase 2
                   qqcw_in(:,:) = qqcw_sav(:,:,itrac)
                   !Not implemented for oslo aerosols
                else
                   fldcw => qqcw_get_field(pbuf, itrac)
                   if(.not. associated(fldcw))then
                      qqcw_in(:,:) = zeroAerosolConcentration(:,:)
                   else
                      qqcw_in(:,:) = fldcw(:,:)
                   end if
                endif

                call wetdepa_v2( pmid, q(:,:,1), pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, q_tmp, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis(:,:,itrac), sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.false., &
                     qqcw=qqcw_in(:,:),  &
                     f_act_conv=f_act_conv, &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     convproc_do_aer=.false., rcscavt=rcscavt, rsscavt=rsscavt,  &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic )

                ptend%q(1:ncol,:,itrac) = ptend%q(1:ncol,:,itrac) + dqdt_tmp(1:ncol,:)

                call outfld(trim(cnst_name(itrac))//'WET', dqdt_tmp(:ncol,:), ncol, lchnk)
                call outfld(trim(cnst_name(itrac))//'SIC', icscavt(:ncol,:),  ncol, lchnk)
                call outfld(trim(cnst_name(itrac))//'SIS', isscavt(:ncol,:),  ncol, lchnk)
                call outfld(trim(cnst_name(itrac))//'SBC', bcscavt(:ncol,:),  ncol, lchnk)
                call outfld(trim(cnst_name(itrac))//'SBS', bsscavt(:ncol,:),  ncol, lchnk)

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+dqdt_tmp(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                if (.not. convproc_do_aer) call outfld( trim(cnst_name(itrac))//'SFWET', sflx(:ncol), ncol, lchnk)
                aerdepwetis(:ncol,itrac) = sflx(:ncol)

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+icscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                if (.not. convproc_do_aer) call outfld( trim(cnst_name(itrac))//'SFSIC', sflx(:ncol), ncol, lchnk)
                if (convproc_do_aer) sflxic = sflx

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+isscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(itrac))//'SFSIS', sflx(:ncol), ncol, lchnk)

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+bcscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(itrac))//'SFSBC', sflx(:ncol), ncol, lchnk)
                if (convproc_do_aer)sflxbc = sflx

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+bsscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(cnst_name(itrac))//'SFSBS', sflx(:ncol), ncol, lchnk)

             else   ! lphase == 2

                dqdt_tmp(:,:) = 0.0_r8
                qqcw_tmp(:,:) = 0.0_r8    ! rce 2010/05/01

                if (convproc_do_aer) then
                   fldcw => qqcw_get_field(pbuf,itrac)
                   if (.not. associated(fldcw)) then
                      call endrun('attempt to access undefined qqcw_sav for fld_cw')
                   end if
                   qqcw_sav(1:ncol,:,itrac) = fldcw(1:ncol,:)
                   !This option yet not implemented for OSLO_AERO
                else
                   fldcw => qqcw_get_field(pbuf, itrac)
                   if(.not. associated(fldcw))then
                      cycle
                   end if
                endif

                call wetdepa_v2(pmid, q(:,:,1), pdel, &
                     dep_inputs%cldt, dep_inputs%cldcu, dep_inputs%cmfdqr, &
                     dep_inputs%evapc, dep_inputs%conicw, dep_inputs%prain, dep_inputs%qme, &
                     dep_inputs%evapr, dep_inputs%totcond, fldcw, dt, &
                     dqdt_tmp, iscavt, dep_inputs%cldvcu, dep_inputs%cldvst, &
                     dlf, fracis_cw, sol_factb, ncol, &
                     scavcoefnv(:,:,jnv), &
                     is_strat_cloudborne=.true.,  &
                     icscavt=icscavt, isscavt=isscavt, bcscavt=bcscavt, bsscavt=bsscavt, &
                     convproc_do_aer=.false., rcscavt=rcscavt, rsscavt=rsscavt,  &
                     sol_facti_in=sol_facti, sol_factic_in=sol_factic )

                fldcw(1:ncol,:) = fldcw(1:ncol,:) + dqdt_tmp(1:ncol,:) * dt

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+dqdt_tmp(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(getCloudTracerName(itrac))//'SFWET', sflx(:ncol), ncol, lchnk)
                aerdepwetcw(:ncol,itrac) = sflx(:ncol)

                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+icscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(getCloudTracerName(itrac))//'SFSIC', sflx(:ncol), ncol, lchnk)
                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+isscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(getCloudTracerName(itrac))//'SFSIS', sflx(:ncol), ncol, lchnk)
                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+bcscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(getCloudTracerName(itrac))//'SFSBC', sflx(:ncol), ncol, lchnk)
                sflx(:)=0._r8
                do ilev=1,pver
                   do icol=1,ncol
                      sflx(icol)=sflx(icol)+bsscavt(icol,ilev)*pdel(icol,ilev)/gravit
                   enddo
                enddo
                call outfld( trim(getCloudTracerName(itrac))//'SFSBS', sflx(:ncol), ncol, lchnk)

             endif

            ! accumulate the deposition flux for the aerosol type
            sflx(:ncol) = 0._r8
            if ( lphase == 1 ) then
               sflx(:ncol) = aerdepwetis(:ncol,itrac)
            else ! lphase == 2
               sflx(:ncol) = aerdepwetcw(:ncol,itrac)
            endif

            ! accumulate the deposition flux for the aerosol type, aerosol mass
            ! if it is sulfate we accumulate in sflx_SFWET_SULFATE_S in addition
            sflx_SFWET_arosol_type(:ncol, aerosolType(itrac)) = sflx_SFWET_arosol_type(:ncol, aerosolType(itrac)) + sflx(:ncol)
            if ( aerosolType(itrac) == AEROSOL_TYPE_SULFATE ) then
               sflx_SFWET_SULFATE_S(:ncol) = sflx_SFWET_SULFATE_S(:ncol) + ( sflx(:ncol) * sulfurMassFraction(itrac) )
            endif

          enddo   ! lspec = 0, nspec_amode(imode)+1
       enddo   ! lphase = 1, 2
    enddo   ! imode = 1, ntot_amode

   ! add the wet deposition rate of the compound aerosols except sulfur to output
    do l_atype=1,N_AEROSOL_TYPES
      if ( l_atype /= AEROSOL_TYPE_SULFATE ) then
         call outfld('wet_'//trim(aerosol_type_name(l_atype)), -1.0_r8 * (sflx_SFWET_arosol_type(:ncol,l_atype)), ncol, lchnk)
      else if ( l_atype == AEROSOL_TYPE_SULFATE ) then
         ! Add in wd_a_h2so4
         call pbuf_get_field(pbuf, idx_wd_a_h2so4, wd_a_h2so4)

         call outfld('wet_'//trim(aerosol_type_name(l_atype)),                &
         -1.0_r8 * ( sflx_SFWET_arosol_type(:ncol,l_atype) + wd_a_h2so4(:ncol) ),    &
            ncol, lchnk)
         call outfld('wet_'//trim(aerosol_type_name(l_atype))//'_S',                                        &
            -1.0_r8 * ( sflx_SFWET_SULFATE_S(:ncol) + ( wd_a_h2so4(:ncol) * sulfurMassFraction(l_h2so4) )),     &
            ncol, lchnk)
         call outfld('wd_a_h2so4_debug', wd_a_h2so4(:ncol), ncol, lchnk)
      endif
   end do

    ! if the user has specified prescribed aerosol dep fluxes then
    ! do not set cam_out dep fluxes according to the prognostic aerosols
    if (.not. aerodep_flx_prescribed()) then
       call oslo_set_srf_wetdep(ncol, aerdepwetis, aerdepwetcw, &
            cam_out%bcphiwet, cam_out%ocphiwet, &
            cam_out%dstwet1, cam_out%dstwet2, cam_out%dstwet3, cam_out%dstwet4)
    endif

  end subroutine oslo_aero_depos_wet

  !===============================================================================
  subroutine oslo_aero_depvel_part( ncol, t, pmid, ram1, fv, vlc_dry, vlc_trb, vlc_grv,  &
       radius_part, density_part, sig_part, moment, lchnk )

    !    calculates surface deposition velocity of particles
    !    L. Zhang, S. Gong, J. Padro, and L. Barrie
    !    A size-seggregated particle dry deposition scheme for an atmospheric aerosol module
    !    Atmospheric Environment, 35, 549-560, 2001.
    !
    !    Authors: X. Liu

    ! !ARGUMENTS:
    real(r8), intent(in) :: t(pcols,pver)       !atm temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)    !atm pressure (Pa)
    real(r8), intent(in) :: fv(pcols)           !friction velocity (m/s)
    real(r8), intent(in) :: ram1(pcols)         !aerodynamical resistance (s/m)
    real(r8), intent(in) :: radius_part(pcols,pver)    ! mean (volume/number) particle radius (m)
    real(r8), intent(in) :: density_part(pcols,pver)   ! density of particle material (kg/m3)
    real(r8), intent(in) :: sig_part(pcols,pver)       ! geometric standard deviation of particles
    integer,  intent(in) :: moment ! moment of size distribution (0 for number, 2 for surface area, 3 for volume)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: lchnk

    real(r8), intent(out) :: vlc_trb(pcols)       !Turbulent deposn velocity (m/s)
    real(r8), intent(out) :: vlc_grv(pcols,pver)       !grav deposn velocity (m/s)
    real(r8), intent(out) :: vlc_dry(pcols,pver)       !dry deposn velocity (m/s)
    !------------------------------------------------------------------------

    !------------------------------------------------------------------------
    ! Local Variables
    integer  :: imode,icol,ilev,ix                  !indices
    real(r8) :: rho                       !atm density (kg/m**3)
    real(r8) :: vsc_dyn_atm(pcols,pver)   ![kg m-1 s-1] Dynamic viscosity of air
    real(r8) :: vsc_knm_atm(pcols,pver)   ![m2 s-1] Kinematic viscosity of atmosphere
    real(r8) :: shm_nbr                   ![frc] Schmidt number
    real(r8) :: stk_nbr                   ![frc] Stokes number
    real(r8) :: mfp_atm(pcols,pver)       ![m] Mean free path of air
    real(r8) :: dff_aer                   ![m2 s-1] Brownian diffusivity of particle
    real(r8) :: slp_crc(pcols,pver)       ![frc] Slip correction factor
    real(r8) :: rss_trb                   ![s m-1] Resistance to turbulent deposition
    real(r8) :: rss_lmn                   ![s m-1] Quasi-laminar layer resistance
    real(r8) :: brownian                  ! collection efficiency for Browning diffusion
    real(r8) :: impaction                 ! collection efficiency for impaction
    real(r8) :: interception              ! collection efficiency for interception
    real(r8) :: stickfrac                 ! fraction of particles sticking to surface
    real(r8) :: radius_moment(pcols,pver) ! median radius (m) for moment
    real(r8) :: lnsig                     ! ln(sig_part)
    real(r8) :: dispersion                ! accounts for influence of size dist dispersion on bulk settling velocity
    ! assuming radius_part is number mode radius * exp(1.5 ln(sigma))
    integer  :: ltype
    real(r8) :: lnd_frc
    real(r8) :: wrk1, wrk2, wrk3

    ! constants
    real(r8) gamma(11)      ! exponent of schmidt number
    !   data gamma/0.54d+00,  0.56d+00,  0.57d+00,  0.54d+00,  0.54d+00, &
    !              0.56d+00,  0.54d+00,  0.54d+00,  0.54d+00,  0.56d+00, &
    !              0.50d+00/
    data gamma/0.56e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.56e+00_r8,  0.56e+00_r8, &
         0.56e+00_r8,  0.50e+00_r8,  0.54e+00_r8,  0.54e+00_r8,  0.54e+00_r8, &
         0.54e+00_r8/
    save gamma

    real(r8) alpha(11)      ! parameter for impaction
    !   data alpha/50.00d+00,  0.95d+00,  0.80d+00,  1.20d+00,  1.30d+00, &
    !               0.80d+00, 50.00d+00, 50.00d+00,  2.00d+00,  1.50d+00, &
    !             100.00d+00/
    data alpha/1.50e+00_r8,   1.20e+00_r8,  1.20e+00_r8,  0.80e+00_r8,  1.00e+00_r8, &
         0.80e+00_r8, 100.00e+00_r8, 50.00e+00_r8,  2.00e+00_r8,  1.20e+00_r8, &
         50.00e+00_r8/
    save alpha

    real(r8) radius_collector(11) ! radius (m) of surface collectors
    !   data radius_collector/-1.00d+00,  5.10d-03,  3.50d-03,  3.20d-03, 10.00d-03, &
    !                          5.00d-03, -1.00d+00, -1.00d+00, 10.00d-03, 10.00d-03, &
    !                         -1.00d+00/
    data radius_collector/10.00e-03_r8,  3.50e-03_r8,  3.50e-03_r8,  5.10e-03_r8,  2.00e-03_r8, &
         5.00e-03_r8, -1.00e+00_r8, -1.00e+00_r8, 10.00e-03_r8,  3.50e-03_r8, &
         -1.00e+00_r8/
    save radius_collector

    integer            :: iwet(11) ! flag for wet surface = 1, otherwise = -1
    !   data iwet/1,   -1,   -1,   -1,   -1,  &
    !            -1,   -1,   -1,    1,   -1,  &
    !             1/
    data iwet/-1,  -1,   -1,   -1,   -1,  &
         -1,   1,   -1,    1,   -1,  &
         -1/
    save iwet


    !------------------------------------------------------------------------

    if(top_lev>1) then
       vlc_grv(:ncol,:top_lev-1) = 0._r8
       vlc_dry(:ncol,:top_lev-1) = 0._r8
    endif

    do ilev=top_lev,pver
       do icol=1,ncol

          lnsig = log(sig_part(icol,ilev))
          ! use a maximum radius of 50 microns when calculating deposition velocity
          radius_moment(icol,ilev) = min(50.0e-6_r8,radius_part(icol,ilev))*   &
               exp((float(moment)-1.5_r8)*lnsig*lnsig)
          dispersion = exp(2._r8*lnsig*lnsig)

          rho=pmid(icol,ilev)/rair/t(icol,ilev)

          ! Quasi-laminar layer resistance: call rss_lmn_get
          ! Size-independent thermokinetic properties
          vsc_dyn_atm(icol,ilev) = 1.72e-5_r8 * ((t(icol,ilev)/273.0_r8)**1.5_r8) * 393.0_r8 / &
               (t(icol,ilev)+120.0_r8)      ![kg m-1 s-1] RoY94 p. 102
          mfp_atm(icol,ilev) = 2.0_r8 * vsc_dyn_atm(icol,ilev) / &   ![m] SeP97 p. 455
               (pmid(icol,ilev)*sqrt(8.0_r8/(pi*rair*t(icol,ilev))))
          vsc_knm_atm(icol,ilev) = vsc_dyn_atm(icol,ilev) / rho ![m2 s-1] Kinematic viscosity of air

          slp_crc(icol,ilev) = 1.0_r8 + mfp_atm(icol,ilev) * &
               (1.257_r8+0.4_r8*exp(-1.1_r8*radius_moment(icol,ilev)/(mfp_atm(icol,ilev)))) / &
               radius_moment(icol,ilev)   ![frc] Slip correction factor SeP97 p. 464
          vlc_grv(icol,ilev) = (4.0_r8/18.0_r8) * radius_moment(icol,ilev)*radius_moment(icol,ilev)*density_part(icol,ilev)* &
               gravit*slp_crc(icol,ilev) / vsc_dyn_atm(icol,ilev) ![m s-1] Stokes' settling velocity SeP97 p. 466
          vlc_grv(icol,ilev) = vlc_grv(icol,ilev) * dispersion

          vlc_dry(icol,ilev)=vlc_grv(icol,ilev)
       enddo
    enddo
    ilev=pver  ! only look at bottom level for next part
    do icol=1,ncol
       dff_aer = boltz * t(icol,ilev) * slp_crc(icol,ilev) / &    ![m2 s-1]
            (6.0_r8*pi*vsc_dyn_atm(icol,ilev)*radius_moment(icol,ilev)) !SeP97 p.474
       shm_nbr = vsc_knm_atm(icol,ilev) / dff_aer                        ![frc] SeP97 p.972

       wrk2 = 0._r8
       wrk3 = 0._r8
       do ltype = 1,n_land_type
          lnd_frc = fraction_landuse(icol,ltype,lchnk)
          if ( lnd_frc /= 0._r8 ) then
             brownian = shm_nbr**(-gamma(ltype))
             if (radius_collector(ltype) > 0.0_r8) then
                !       vegetated surface
                stk_nbr = vlc_grv(icol,ilev) * fv(icol) / (gravit*radius_collector(ltype))
                interception = 2.0_r8*(radius_moment(icol,ilev)/radius_collector(ltype))**2.0_r8
             else
                !       non-vegetated surface
                stk_nbr = vlc_grv(icol,ilev) * fv(icol) * fv(icol) / (gravit*vsc_knm_atm(icol,ilev))  ![frc] SeP97 p.965
                interception = 0.0_r8
             endif
             impaction = (stk_nbr/(alpha(ltype)+stk_nbr))**2.0_r8

             if (iwet(ltype) > 0) then
                stickfrac = 1.0_r8
             else
                stickfrac = exp(-sqrt(stk_nbr))
                if (stickfrac < 1.0e-10_r8) stickfrac = 1.0e-10_r8
             endif
             rss_lmn = 1.0_r8 / (3.0_r8 * fv(icol) * stickfrac * (brownian+interception+impaction))
             rss_trb = ram1(icol) + rss_lmn + ram1(icol)*rss_lmn*vlc_grv(icol,ilev)

             wrk1 = 1.0_r8 / rss_trb
             wrk2 = wrk2 + lnd_frc*( wrk1 )
             wrk3 = wrk3 + lnd_frc*( wrk1 + vlc_grv(icol,ilev) )
          endif
       enddo  ! n_land_type
       vlc_trb(icol) = wrk2
       vlc_dry(icol,ilev) = wrk3
    enddo !ncol

  end subroutine oslo_aero_depvel_part

  !===============================================================================
  subroutine oslo_set_srf_wetdep(ncol, aerdepwetis, aerdepwetcw, &
       bcphiwet, ocphiwet, dstwet1, dstwet2, dstwet3, dstwet4)

    ! Set surface wet deposition fluxes passed to coupler.

    ! Arguments:
    integer , intent(in)  :: ncol
    real(r8), intent(in)  :: aerdepwetis(:,:)  ! aerosol wet deposition (interstitial)
    real(r8), intent(in)  :: aerdepwetcw(:,:)  ! aerosol wet deposition (cloud water)
    real(r8), intent(out) :: bcphiwet(:)
    real(r8), intent(out) :: ocphiwet(:)
    real(r8), intent(out) :: dstwet1(:)
    real(r8), intent(out) :: dstwet2(:)
    real(r8), intent(out) :: dstwet3(:)
    real(r8), intent(out) :: dstwet4(:)

    ! Local variables:
    integer :: icol
    !----------------------------------------------------------------------------

    bcphiwet(:) = 0._r8
    ocphiwet(:) = 0._r8
    dstwet1(:)  = 0._r8
    dstwet2(:)  = 0._r8
    dstwet3(:)  = 0._r8
    dstwet4(:)  = 0._r8

    ! derive cam_out variables from deposition fluxes
    !  note: wet deposition fluxes are negative into surface,
    !        dry deposition fluxes are positive into surface.
    !        srf models want positive definite fluxes.
    do icol = 1,ncol
       ! black carbon fluxes
       ! note: bc_ax is assumed not to exist in cloud water
       bcphiwet(icol) = -(aerdepwetis(icol,l_bc_ni)+aerdepwetcw(icol,l_bc_ni)+ &
                       aerdepwetis(icol,l_bc_ai)+aerdepwetcw(icol,l_bc_ai)+ &
                       aerdepwetis(icol,l_bc_a )+aerdepwetcw(icol,l_bc_a )+ &
                       aerdepwetis(icol,l_bc_ac)+aerdepwetcw(icol,l_bc_ac)+ &
                       aerdepwetis(icol,l_bc_n )+aerdepwetcw(icol,l_bc_n )+ &
                       aerdepwetis(icol,l_bc_ax))

       ! organic carbon fluxes
       ocphiwet(icol) = -(aerdepwetis(icol,l_om_ni)+aerdepwetcw(icol,l_om_ni)+ &
                       aerdepwetis(icol,l_om_ai)+aerdepwetcw(icol,l_om_ai)+ &
                       aerdepwetis(icol,l_om_ac)+aerdepwetcw(icol,l_om_ac))

       ! dust fluxes
       ! bulk bin1 (fine) dust deposition equals accumulation mode deposition:
       !  A. Simple: Assign all coarse-mode dust to bulk size bin 3:
       dstwet1(icol) = -(aerdepwetis(icol,l_dst_a2)+aerdepwetcw(icol,l_dst_a2))
       dstwet2(icol) = 0._r8
       dstwet3(icol) = -(aerdepwetis(icol,l_dst_a3)+aerdepwetcw(icol,l_dst_a3))
       dstwet4(icol) = 0._r8

    enddo

  end subroutine oslo_set_srf_wetdep

  !===============================================================================
  subroutine oslo_set_srf_drydep(ncol, aerdepdryis, aerdepdrycw, &
       bcphidry, bcphodry, ocphidry, ocphodry, dstdry1, dstdry2, dstdry3, dstdry4)

    ! Set surface dry deposition fluxes passed to coupler.

    ! Arguments:
    integer , intent(in)  :: ncol
    real(r8), intent(in)  :: aerdepdryis(:,:)  ! aerosol dry deposition (interstitial)
    real(r8), intent(in)  :: aerdepdrycw(:,:)  ! aerosol dry deposition (cloud water)
    real(r8), intent(out) :: bcphidry(:)
    real(r8), intent(out) :: bcphodry(:)
    real(r8), intent(out) :: ocphidry(:)
    real(r8), intent(out) :: ocphodry(:)
    real(r8), intent(out) :: dstdry1(:)
    real(r8), intent(out) :: dstdry2(:)
    real(r8), intent(out) :: dstdry3(:)
    real(r8), intent(out) :: dstdry4(:)

    ! Local variables:
    integer :: icol
    !----------------------------------------------------------------------------

    bcphidry(:) = 0._r8
    bcphodry(:) = 0._r8
    ocphidry(:) = 0._r8
    ocphodry(:) = 0._r8
    dstdry1(:)  = 0._r8
    dstdry2(:)  = 0._r8
    dstdry3(:)  = 0._r8
    dstdry4(:)  = 0._r8

    ! wet deposition fluxes are negative into surface,
    ! dry deposition fluxes are positive into surface.
    ! srf models want positive definite fluxes.
    do icol = 1, ncol
       ! black carbon fluxes
       bcphidry(icol) = aerdepdryis(icol,l_bc_ni)+aerdepdrycw(icol,l_bc_ni)+ &
                     aerdepdryis(icol,l_bc_ai)+aerdepdrycw(icol,l_bc_ai)+ &
                     aerdepdryis(icol,l_bc_a )+aerdepdrycw(icol,l_bc_a )+ &
                     aerdepdryis(icol,l_bc_ac)+aerdepdrycw(icol,l_bc_ac)
       bcphodry(icol) = aerdepdryis(icol,l_bc_n )+aerdepdrycw(icol,l_bc_n )+ &
                     aerdepdryis(icol,l_bc_ax)+aerdepdrycw(icol,l_bc_ax)

       ! organic carbon fluxes
       ! djlo : skipped the bc_a contribution (was about om !)
       ! ocphidry(icol) = aerdepdryis(icol,l_om_ni)+aerdepdrycw(icol,l_om_ni)+ &
       !               aerdepdryis(icol,l_om_ai)+aerdepdrycw(icol,l_om_ai)+ &
       !               aerdepdryis(icol,l_om_ac)+aerdepdrycw(icol,l_om_ac)
       ocphidry(icol) = 0._r8
       ocphodry(icol) = 0._r8

       ! dust fluxes
       ! bulk bin1 (fine) : dust deposition equals accumulation mode deposition:
       ! bulk bin 3       : A. Simple: Assign all coarse-mode dust to bulk size bin 3:
       ! bulk bins 2-4    : two options for partitioning deposition into bins 2-4:
       dstdry1(icol) = aerdepdryis(icol,l_dst_a2)+aerdepdrycw(icol,l_dst_a2)
       dstdry2(icol) = 0._r8
       dstdry3(icol) = aerdepdryis(icol,l_dst_a3)+aerdepdrycw(icol,l_dst_a3)
       dstdry4(icol) = 0._r8
    enddo

  end subroutine oslo_set_srf_drydep

  !===============================================================================
  subroutine calcram(ncol, landfrac, icefrac, ocnfrac, obklen, &
       ustar, ram1in, ram1, t, pmid, pdel, fvin, fv)

    ! Calc aerodynamic resistance over oceans and sea ice (comes in from land model)
    ! from Seinfeld and Pandis, p.963. Author: Natalie Mahowald

    integer  , intent(in)  :: ncol
    real(r8) , intent(in)  :: ram1in(pcols)   ! aerodynamical resistance (s/m)
    real(r8) , intent(in)  :: fvin(pcols)     ! sfc frc vel from land
    real(r8) , intent(out) :: ram1(pcols)     ! aerodynamical resistance (s/m)
    real(r8) , intent(out) :: fv(pcols)       ! sfc frc vel from land
    real(r8) , intent(in)  :: obklen(pcols)   ! obklen
    real(r8) , intent(in)  :: ustar(pcols)    ! sfc fric vel
    real(r8) , intent(in)  :: landfrac(pcols) ! land fraction
    real(r8) , intent(in)  :: icefrac(pcols)  ! ice fraction
    real(r8) , intent(in)  :: ocnfrac(pcols)  ! ocean fraction
    real(r8) , intent(in)  :: t(pcols)        ! atm temperature (K)
    real(r8) , intent(in)  :: pmid(pcols)     ! atm pressure (Pa)
    real(r8) , intent(in)  :: pdel(pcols)     ! atm pressure (Pa)

    ! local variables
    real(r8), parameter :: zzocen = 0.0001_r8   ! Ocean aerodynamic roughness length
    real(r8), parameter :: zzsice = 0.0400_r8   ! Sea ice aerodynamic roughness length
    real(r8), parameter :: xkar   = 0.4_r8      ! Von Karman constant
    real(r8) :: z,psi,psi0,nu,nu0,temp,ram
    integer  :: icol

    do icol = 1,ncol
       z=pdel(icol)*rair*t(icol)/pmid(icol)/gravit/2.0_r8   !use half the layer height like Ganzefeld and Lelieveld, 1995
       if(obklen(icol) == 0) then
          psi=0._r8
          psi0=0._r8
       else
          psi=min(max(z/obklen(icol),-1.0_r8),1.0_r8)
          psi0=min(max(zzocen/obklen(icol),-1.0_r8),1.0_r8)
       endif
       temp=z/zzocen
       if(icefrac(icol) > 0.5_r8) then
          if(obklen(icol)>0) then
             psi0=min(max(zzsice/obklen(icol),-1.0_r8),1.0_r8)
          else
             psi0=0.0_r8
          endif
          temp=z/zzsice
       endif
       if(psi> 0._r8) then
          ram=1/xkar/ustar(icol)*(log(temp)+4.7_r8*(psi-psi0))
       else
          nu=(1.00_r8-15.000_r8*psi)**(.25_r8)
          nu0=(1.000_r8-15.000_r8*psi0)**(.25_r8)
          if(ustar(icol).ne.0._r8) then
             ram=1/xkar/ustar(icol)*(log(temp) &
                  +log(((nu0**2+1.00_r8)*(nu0+1.0_r8)**2)/((nu**2+1.0_r8)*(nu+1.00_r8)**2)) &
                  +2.0_r8*(atan(nu)-atan(nu0)))
          else
             ram=0._r8
          endif
       endif
       if(landfrac(icol) < 0.000000001_r8) then
          fv(icol)=ustar(icol)
          ram1(icol)=ram
       else
          fv(icol)=fvin(icol)
          ram1(icol)=ram1in(icol)
       endif
    enddo

    ! fvitt -- fv == 0 causes a floating point exception in dry dep of sea salts and dust
    where ( fv(:ncol) == 0._r8 )
       fv(:ncol) = 1.e-12_r8
    endwhere
  end subroutine calcram

  !==============================================================================
  subroutine oslo_aero_wetdep_init()

    ! Initialize module variables for wet deposition

    cld_idx             = pbuf_get_index('CLD')
    qme_idx             = pbuf_get_index('QME')
    prain_idx           = pbuf_get_index('PRAIN')
    nevapr_idx          = pbuf_get_index('NEVAPR')

    icwmrdp_idx         = pbuf_get_index('ICWMRDP')
    rprddp_idx          = pbuf_get_index('RPRDDP')
    icwmrsh_idx         = pbuf_get_index('ICWMRSH')
    rprdsh_idx          = pbuf_get_index('RPRDSH')
    sh_frac_idx         = pbuf_get_index('SH_FRAC' )
    dp_frac_idx         = pbuf_get_index('DP_FRAC')
    nevapr_shcu_idx     = pbuf_get_index('NEVAPR_SHCU')
    nevapr_dpcu_idx     = pbuf_get_index('NEVAPR_DPCU')

    call cnst_get_ind('CLDICE', ixcldice)
    call cnst_get_ind('CLDLIQ', ixcldliq)

  endsubroutine oslo_aero_wetdep_init

  !==============================================================================
  subroutine wetdep_inputs_set(ncol, pmid, pdel, t, q, pbuf, inputs )

    ! gather up the inputs needed for the wetdepa routines

    ! arguments
    integer,               intent(in)  :: ncol
    real(r8) ,             intent(in)  :: pmid(pcols,pver) ! Model level pressures (Pa)
    real(r8) ,             intent(in)  :: pdel(pcols,pver) ! model layer thickness (Pa)
    real(r8) ,             intent(in)  :: t(pcols,pver)    ! Model level temperatures (K)
    real(r8) ,             intent(in)  :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)          ! physics buffer
    type(wetdep_inputs_t), intent(out) :: inputs           ! collection of wetdepa inputs

    ! local variables
    real(r8), pointer :: icwmrdp(:,:)       ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: rprddp(:,:)        ! rain production, deep convection
    real(r8), pointer :: icwmrsh(:,:)       ! in cloud water mixing ratio, deep convection
    real(r8), pointer :: rprdsh(:,:)        ! rain production, deep convection
    real(r8), pointer :: sh_frac(:,:)       ! Shallow convective cloud fraction
    real(r8), pointer :: dp_frac(:,:)       ! Deep convective cloud fraction
    real(r8), pointer :: evapcsh(:,:)       ! Evaporation rate of shallow convective precipitation >=0.
    real(r8), pointer :: evapcdp(:,:)       ! Evaporation rate of deep    convective precipitation >=0.
    real(r8)          :: rainmr(pcols,pver) ! mixing ratio of rain within cloud volume
    real(r8)          :: cldst(pcols,pver)  ! Stratiform cloud fraction
    integer           :: itim

    itim = pbuf_old_tim_idx()

    call pbuf_get_field(pbuf, cld_idx,         inputs%cldt, start=(/1,1,itim/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, qme_idx,         inputs%qme     )
    call pbuf_get_field(pbuf, prain_idx,       inputs%prain   )
    call pbuf_get_field(pbuf, nevapr_idx,      inputs%evapr   )
    call pbuf_get_field(pbuf, icwmrdp_idx,     icwmrdp )
    call pbuf_get_field(pbuf, icwmrsh_idx,     icwmrsh )
    call pbuf_get_field(pbuf, rprddp_idx,      rprddp  )
    call pbuf_get_field(pbuf, rprdsh_idx,      rprdsh  )
    call pbuf_get_field(pbuf, sh_frac_idx,     sh_frac )
    call pbuf_get_field(pbuf, dp_frac_idx,     dp_frac )
    call pbuf_get_field(pbuf, nevapr_shcu_idx, evapcsh )
    call pbuf_get_field(pbuf, nevapr_dpcu_idx, evapcdp )

    inputs%cldcu(:ncol,:)  = dp_frac(:ncol,:) + sh_frac(:ncol,:)
    cldst(:ncol,:)  = inputs%cldt(:ncol,:) - inputs%cldcu(:ncol,:)  ! Stratiform cloud fraction
    inputs%evapc(:ncol,:)  = evapcsh(:ncol,:) + evapcdp(:ncol,:)
    inputs%cmfdqr(:ncol,:) = rprddp(:ncol,:)  + rprdsh(:ncol,:)

    ! sum deep and shallow convection contributions
    if (cam_physpkg_is('cam5') .or. cam_physpkg_is('cam6')) then
       inputs%conicw(:ncol,:) = (icwmrdp(:ncol,:)*dp_frac(:ncol,:) + icwmrsh(:ncol,:)*sh_frac(:ncol,:))/ &
            max(0.01_r8, sh_frac(:ncol,:) + dp_frac(:ncol,:))
    else
       inputs%conicw(:ncol,:) = icwmrdp(:ncol,:) + icwmrsh(:ncol,:)
    end if

    inputs%totcond(:ncol,:) = q(:ncol,:,ixcldliq) + q(:ncol,:,ixcldice)

    call clddiag( t, pmid, pdel, inputs%cmfdqr, inputs%evapc, &
         inputs%cldt, inputs%cldcu, cldst, inputs%qme, inputs%evapr, &
         inputs%prain, inputs%cldv, inputs%cldvcu, inputs%cldvst, rainmr, ncol )

  end subroutine wetdep_inputs_set

  !==============================================================================
  subroutine clddiag(t, pmid, pdel, cmfdqr, evapc, &
       cldt, cldcu, cldst, cme, evapr, &
       prain, cldv, cldvcu, cldvst, rain, ncol)

    ! ------------------------------------------------------------------------------------
    ! Estimate the cloudy volume which is occupied by rain or cloud water as
    ! the max between the local cloud amount or the
    ! sum above of (cloud*positive precip production)      sum total precip from above
    !              ----------------------------------   x ------------------------
    ! sum above of     (positive precip           )        sum positive precip from above
    ! Author: P. Rasch,  Sungsu Park. Mar.2010
    ! ------------------------------------------------------------------------------------

    ! Input arguments:
    real(r8) , intent(in) :: t(pcols,pver)      ! temperature (K)
    real(r8) , intent(in) :: pmid(pcols,pver)   ! pressure at layer midpoints
    real(r8) , intent(in) :: pdel(pcols,pver)   ! pressure difference across layers
    real(r8) , intent(in) :: cmfdqr(pcols,pver) ! dq/dt due to convective rainout
    real(r8) , intent(in) :: evapc(pcols,pver)  ! Evaporation rate of convective precipitation ( >= 0 )
    real(r8) , intent(in) :: cldt(pcols,pver)   ! total cloud fraction
    real(r8) , intent(in) :: cldcu(pcols,pver)  ! Cumulus cloud fraction
    real(r8) , intent(in) :: cldst(pcols,pver)  ! Stratus cloud fraction
    real(r8) , intent(in) :: cme(pcols,pver)    ! rate of cond-evap within the cloud
    real(r8) , intent(in) :: evapr(pcols,pver)  ! rate of evaporation of falling precipitation (kg/kg/s)
    real(r8) , intent(in) :: prain(pcols,pver)  ! rate of conversion of condensate to precipitation (kg/kg/s)
    integer  , intent(in) :: ncol

    ! Output arguments:
    real(r8), intent(out) :: cldv(pcols,pver)     ! fraction occupied by rain or cloud water
    real(r8), intent(out) :: cldvcu(pcols,pver)   ! Convective precipitation volume
    real(r8), intent(out) :: cldvst(pcols,pver)   ! Stratiform precipitation volume
    real(r8), intent(out) :: rain(pcols,pver)     ! mixing ratio of rain (kg/kg)

    ! Local variables:
    integer  icol, ilev
    real(r8) convfw            ! used in fallspeed calculation; taken from findmcnew
    real(r8) sumppr(pcols)     ! precipitation rate (kg/m2-s)
    real(r8) sumpppr(pcols)    ! sum of positive precips from above
    real(r8) cldv1(pcols)      ! precip weighted cloud fraction from above
    real(r8) lprec             ! local production rate of precip (kg/m2/s)
    real(r8) lprecp            ! local production rate of precip (kg/m2/s) if positive
    real(r8) rho               ! air density
    real(r8) vfall
    real(r8) sumppr_cu(pcols)  ! Convective precipitation rate (kg/m2-s)
    real(r8) sumpppr_cu(pcols) ! Sum of positive convective precips from above
    real(r8) cldv1_cu(pcols)   ! Convective precip weighted convective cloud fraction from above
    real(r8) lprec_cu          ! Local production rate of convective precip (kg/m2/s)
    real(r8) lprecp_cu         ! Local production rate of convective precip (kg/m2/s) if positive
    real(r8) sumppr_st(pcols)  ! Stratiform precipitation rate (kg/m2-s)
    real(r8) sumpppr_st(pcols) ! Sum of positive stratiform precips from above
    real(r8) cldv1_st(pcols)   ! Stratiform precip weighted stratiform cloud fraction from above
    real(r8) lprec_st          ! Local production rate of stratiform precip (kg/m2/s)
    real(r8) lprecp_st         ! Local production rate of stratiform precip (kg/m2/s) if positive
    ! -----------------------------------------------------------------------

    convfw = 1.94_r8*2.13_r8*sqrt(rhoh2o*gravit*2.7e-4_r8)
    do icol=1,ncol
       sumppr(icol) = 0._r8
       cldv1(icol) = 0._r8
       sumpppr(icol) = 1.e-36_r8
       sumppr_cu(icol)  = 0._r8
       cldv1_cu(icol)   = 0._r8
       sumpppr_cu(icol) = 1.e-36_r8
       sumppr_st(icol)  = 0._r8
       cldv1_st(icol)   = 0._r8
       sumpppr_st(icol) = 1.e-36_r8
    end do

    do ilev = 1,pver
       do icol = 1,ncol
          cldv(icol,ilev) = &
               max(min(1._r8, &
               cldv1(icol)/sumpppr(icol) &
               )*sumppr(icol)/sumpppr(icol), &
               cldt(icol,ilev) &
               )
          lprec = pdel(icol,ilev)/gravit  * (prain(icol,ilev)+cmfdqr(icol,ilev)-evapr(icol,ilev))
          lprecp = max(lprec,1.e-30_r8)
          cldv1(icol) = cldv1(icol)  + cldt(icol,ilev)*lprecp
          sumppr(icol) = sumppr(icol) + lprec
          sumpppr(icol) = sumpppr(icol) + lprecp

          ! For convective precipitation volume at the top interface of each layer. Neglect the current layer.
          cldvcu(icol,ilev)   = max(min(1._r8,cldv1_cu(icol)/sumpppr_cu(icol))*(sumppr_cu(icol)/sumpppr_cu(icol)),0._r8)
          lprec_cu      = (pdel(icol,ilev)/gravit)*(cmfdqr(icol,ilev)-evapc(icol,ilev))
          lprecp_cu     = max(lprec_cu,1.e-30_r8)
          cldv1_cu(icol)   = cldv1_cu(icol) + cldcu(icol,ilev)*lprecp_cu
          sumppr_cu(icol)  = sumppr_cu(icol) + lprec_cu
          sumpppr_cu(icol) = sumpppr_cu(icol) + lprecp_cu

          ! For stratiform precipitation volume at the top interface of each layer. Neglect the current layer.
          cldvst(icol,ilev)   = max(min(1._r8,cldv1_st(icol)/sumpppr_st(icol))*(sumppr_st(icol)/sumpppr_st(icol)),0._r8)
          lprec_st      = (pdel(icol,ilev)/gravit)*(prain(icol,ilev)-evapr(icol,ilev))
          lprecp_st     = max(lprec_st,1.e-30_r8)
          cldv1_st(icol)   = cldv1_st(icol) + cldst(icol,ilev)*lprecp_st
          sumppr_st(icol)  = sumppr_st(icol) + lprec_st
          sumpppr_st(icol) = sumpppr_st(icol) + lprecp_st

          rain(icol,ilev) = 0._r8
          if(t(icol,ilev) > tmelt) then
             rho = pmid(icol,ilev)/(rair*t(icol,ilev))
             vfall = convfw/sqrt(rho)
             rain(icol,ilev) = sumppr(icol)/(rho*vfall)
             if (rain(icol,ilev)<1.e-14_r8) rain(icol,ilev) = 0._r8
          endif
       end do
    end do

  end subroutine clddiag

  !==============================================================================
  subroutine wetdepa_v2(                                    &
       p, q, pdel, cldt, cldc,                              &
       cmfdqr, evapc, conicw, precs, conds,                 &
       evaps, cwat, tracer, deltat, scavt,                  &
       iscavt, cldvcu, cldvst, dlf, fracis,                 &
       sol_fact, ncol, scavcoef, is_strat_cloudborne, qqcw, &
       f_act_conv, icscavt, isscavt, bcscavt, bsscavt,      &
       convproc_do_aer, rcscavt, rsscavt,                   &
       sol_facti_in, sol_factic_in )

    !-----------------------------------------------------------------------
    ! scavenging code for very soluble aerosols
    ! This is the CAM5 version of wetdepa.
    !-----------------------------------------------------------------------

    real(r8), intent(in) ::&
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cmfdqr(pcols,pver),   &! rate of production of convective precip
         evapc(pcols,pver),    &! Evaporation rate of convective precipitation
         conicw(pcols,pver),   &! convective cloud water
         cwat(pcols,pver),     &! cloud water amount
         precs(pcols,pver),    &! rate of production of stratiform precip
         conds(pcols,pver),    &! rate of production of condensate
         evaps(pcols,pver),    &! rate of evaporation of precip
         cldvcu(pcols,pver),   &! Convective precipitation area at the top interface of each layer
         cldvst(pcols,pver),   &! Stratiform precipitation area at the top interface of each layer
         dlf(pcols,pver),      &! Detrainment of convective condensate [kg/kg/s]
         deltat,               &! time step
         tracer(pcols,pver)     ! trace species

    ! If subroutine is called with just sol_fact:
    !    sol_fact is used for both in- and below-cloud scavenging
    ! If subroutine is called with optional argument sol_facti_in:
    !    sol_fact  is used for below cloud scavenging
    !    sol_facti is used for in cloud scavenging

    real(r8), intent(in)  :: sol_fact
    integer,  intent(in)  :: ncol
    real(r8), intent(in)  :: scavcoef(pcols,pver) ! Dana and Hales coefficient (/mm) (0.1 if not MODAL_AERO)
    real(r8), intent(out) :: scavt(pcols,pver)    ! scavenging tend
    real(r8), intent(out) :: iscavt(pcols,pver)   ! incloud scavenging tends
    real(r8), intent(out) :: fracis(pcols,pver)   ! fraction of species not scavenged

    ! Setting is_strat_cloudborne=.true. indicates that tracer is stratiform-cloudborne aerosol.
    !   This is only used by MAM code.  The optional args qqcw and f_act_conv are not referenced
    !   in this case.
    ! Setting is_strat_cloudborne=.false. is being used to indicate that the tracers are the
    !   interstitial modal aerosols.  In this case the optional qqcw (the cloud borne mixing ratio
    !   corresponding to the interstitial aerosol) must be provided, as well as the optional f_act_conv.

    logical,  intent(in), optional :: is_strat_cloudborne
    real(r8), intent(in), optional :: qqcw(pcols,pver)
    real(r8), intent(in), optional :: f_act_conv(pcols,pver)

    real(r8), intent(in), optional :: sol_facti_in   ! solubility factor (frac of aerosol scavenged in cloud)
    real(r8), intent(in), optional :: sol_factic_in(pcols,pver)  ! sol_facti_in for convective clouds

    real(r8), intent(out), optional :: icscavt(pcols,pver)     ! incloud, convective
    real(r8), intent(out), optional :: isscavt(pcols,pver)     ! incloud, stratiform
    real(r8), intent(out), optional :: bcscavt(pcols,pver)     ! below cloud, convective
    real(r8), intent(out), optional :: bsscavt(pcols,pver)     ! below cloud, stratiform

    ! Setting convproc_do_aer=.true. removes the resuspension term from bcscavt and
    ! bsscavt and returns those terms as rcscavt and rsscavt respectively.
    logical,  intent(in),  optional :: convproc_do_aer
    real(r8), intent(out), optional :: rcscavt(pcols,pver)     ! resuspension, convective
    real(r8), intent(out), optional :: rsscavt(pcols,pver)     ! resuspension, stratiform

    ! local variables
    integer  :: icol, ilev
    logical  :: out_resuspension
    real(r8) :: omsm                 ! 1 - (a small number)
    real(r8) :: clds(pcols)          ! stratiform cloud fraction
    real(r8) :: fracev(pcols)        ! fraction of precip from above that is evaporating
    real(r8) :: fracev_cu(pcols)     ! Fraction of convective precip from above that is evaporating
    real(r8) :: fracp(pcols)         ! fraction of cloud water converted to precip
    real(r8) :: pdog(pcols)          ! work variable (pdel/gravit)
    real(r8) :: rpdog(pcols)         ! work variable (gravit/pdel)
    real(r8) :: precabc(pcols)       ! conv precip from above (work array)
    real(r8) :: precabs(pcols)       ! strat precip from above (work array)
    real(r8) :: rat(pcols)           ! ratio of amount available to amount removed
    real(r8) :: scavab(pcols)        ! scavenged tracer flux from above (work array)
    real(r8) :: scavabc(pcols)       ! scavenged tracer flux from above (work array)
    real(r8) :: srcc(pcols)          ! tend for convective rain
    real(r8) :: srcs(pcols)          ! tend for stratiform rain
    real(r8) :: srct(pcols)          ! work variable

    real(r8) :: fins(pcols)          ! fraction of rem. rate by strat rain
    real(r8) :: finc(pcols)          ! fraction of rem. rate by conv. rain
    real(r8) :: conv_scav_ic(pcols)  ! convective scavenging incloud
    real(r8) :: conv_scav_bc(pcols)  ! convective scavenging below cloud
    real(r8) :: st_scav_ic(pcols)    ! stratiform scavenging incloud
    real(r8) :: st_scav_bc(pcols)    ! stratiform scavenging below cloud

    real(r8) :: odds(pcols)          ! limit on removal rate (proportional to prec)
    real(r8) :: dblchek(pcols)
    logical  :: found

    real(r8) :: trac_qqcw(pcols)
    real(r8) :: tracer_incu(pcols)
    real(r8) :: tracer_mean(pcols)

    ! For stratiform cloud, cloudborne aerosol is treated explicitly,
    !    and sol_facti is 1.0 for cloudborne, 0.0 for interstitial.
    ! For convective cloud, cloudborne aerosol is not treated explicitly,
    !    and sol_factic is 1.0 for both cloudborne and interstitial.

    real(r8) :: sol_facti              ! in cloud fraction of aerosol scavenged
    real(r8) :: sol_factb              ! below cloud fraction of aerosol scavenged
    real(r8) :: sol_factic(pcols,pver) ! in cloud fraction of aerosol scavenged for convective clouds

    real(r8) :: rdeltat
    ! ------------------------------------------------------------------------

    omsm = 1._r8-2*epsilon(1._r8) ! used to prevent roundoff errors below zero

    ! default (if other sol_facts aren't in call, set all to required sol_fact)
    sol_facti = sol_fact
    sol_factb = sol_fact

    if ( present(sol_facti_in) )  sol_facti = sol_facti_in

    sol_factic  = sol_facti
    if ( present(sol_factic_in ) )  sol_factic  = sol_factic_in

    ! Determine whether resuspension fields are output.
    out_resuspension = .false.
    if (present(convproc_do_aer)) then
       if (convproc_do_aer) then
          if (present(bcscavt) .and. present(bsscavt) .and. &
              present(rcscavt) .and. present(rsscavt) ) then
             out_resuspension = .true.
          else
             call endrun('wetdepa_v2: bcscavt, bsscavt, rcscavt, rsscavt'// &
                  ' must be present when convproc_do_aero true')
          end if
       end if
    end if

    ! this section of code is for highly soluble aerosols,
    ! the assumption is that within the cloud that
    ! all the tracer is in the cloud water
    !
    ! for both convective and stratiform clouds,
    ! the fraction of cloud water converted to precip defines
    ! the amount of tracer which is pulled out.

    precabs(:ncol) = 0.0_r8
    precabc(:ncol) = 0.0_r8
    scavab(:ncol)  = 0.0_r8
    scavabc(:ncol) = 0.0_r8

    do ilev = 1, pver
       do icol = 1, ncol

          clds(icol)  = cldt(icol,ilev) - cldc(icol,ilev)
          pdog(icol)  = pdel(icol,ilev)/gravit
          rpdog(icol) = gravit/pdel(icol,ilev)
          rdeltat  = 1.0_r8/deltat

          ! ****************** Evaporation **************************
          ! calculate the fraction of strat precip from above
          !                 which evaporates within this layer
          fracev(icol) = evaps(icol,ilev)*pdog(icol) &
               /max(1.e-12_r8,precabs(icol))

          ! trap to ensure reasonable ratio bounds
          fracev(icol) = max(0._r8,min(1._r8,fracev(icol)))

          ! Same as above but convective precipitation part
          fracev_cu(icol) = evapc(icol,ilev)*pdog(icol)/max(1.e-12_r8,precabc(icol))
          fracev_cu(icol) = max(0._r8,min(1._r8,fracev_cu(icol)))

          ! ****************** Convection ***************************
          !
          ! set odds proportional to fraction of the grid box that is swept by the
          ! precipitation =precabc/rhoh20*(area of sphere projected on plane
          !                                /volume of sphere)*deltat
          ! assume the radius of a raindrop is 1 e-3 m from Rogers and Yau,
          ! unless the fraction of the area that is cloud is less than odds, in which
          ! case use the cloud fraction (assumes precabs is in kg/m2/s)
          ! is really: precabs*3/4/1000./1e-3*deltat
          ! here I use .1 from Balkanski
          !
          ! use a local rate of convective rain production for incloud scav
          !
          ! Fraction of convective cloud water converted to rain.  This version is used
          ! in 2 of the 3 branches below before fracp is reused in the stratiform calc.
          ! NB: In below formula for fracp conicw is a LWC/IWC that has already
          !     precipitated out, i.e., conicw does not contain precipitation

          fracp(icol) = cmfdqr(icol,ilev)*deltat / &
               max( 1.e-12_r8, cldc(icol,ilev)*conicw(icol,ilev) + (cmfdqr(icol,ilev)+dlf(icol,ilev))*deltat )
          fracp(icol) = max( min( 1._r8, fracp(icol)), 0._r8 )

          if ( present(is_strat_cloudborne) ) then

             if ( is_strat_cloudborne ) then

                ! convective scavenging

                conv_scav_ic(icol) = 0._r8

                conv_scav_bc(icol) = 0._r8

                ! stratiform scavenging

                fracp(icol) = precs(icol,ilev)*deltat / &
                     max( 1.e-12_r8, cwat(icol,ilev) + precs(icol,ilev)*deltat )
                fracp(icol) = max( 0._r8, min(1._r8, fracp(icol)) )
                st_scav_ic(icol) = sol_facti *fracp(icol)*tracer(icol,ilev)*rdeltat

                st_scav_bc(icol) = 0._r8

             else

                ! convective scavenging

                trac_qqcw(icol) = min(qqcw(icol,ilev), &
                     tracer(icol,ilev)*( clds(icol)/max( 0.01_r8, 1._r8-clds(icol) ) ) )

                tracer_incu(icol) = f_act_conv(icol,ilev)*(tracer(icol,ilev) + trac_qqcw(icol))

                conv_scav_ic(icol) = sol_factic(icol,ilev)*cldc(icol,ilev)*fracp(icol)*tracer_incu(icol)*rdeltat

                tracer_mean(icol) = tracer(icol,ilev)*(1._r8 - cldc(icol,ilev)*f_act_conv(icol,ilev)) - &
                     cldc(icol,ilev)*f_act_conv(icol,ilev)*trac_qqcw(icol)
                tracer_mean(icol) = max(0._r8,tracer_mean(icol))

                odds(icol) = precabc(icol)/max(cldvcu(icol,ilev),1.e-5_r8)*scavcoef(icol,ilev)*deltat
                odds(icol) = max(min(1._r8,odds(icol)),0._r8)
                conv_scav_bc(icol) = sol_factb *cldvcu(icol,ilev)*odds(icol)*tracer_mean(icol)*rdeltat


                ! stratiform scavenging

                st_scav_ic(icol) = 0._r8

                odds(icol) = precabs(icol)/max(cldvst(icol,ilev),1.e-5_r8)*scavcoef(icol,ilev)*deltat
                odds(icol) = max(min(1._r8,odds(icol)),0._r8)
                st_scav_bc(icol) = sol_factb *cldvst(icol,ilev)*odds(icol)*tracer_mean(icol)*rdeltat

             end if

          else

             ! convective scavenging

             conv_scav_ic(icol) = sol_factic(icol,ilev)*cldc(icol,ilev)*fracp(icol)*tracer(icol,ilev)*rdeltat

             odds(icol) = precabc(icol)/max(cldvcu(icol,ilev), 1.e-5_r8)*scavcoef(icol,ilev)*deltat
             odds(icol) = max( min(1._r8, odds(icol)), 0._r8)
             conv_scav_bc(icol) = sol_factb*cldvcu(icol,ilev)*odds(icol)*tracer(icol,ilev)*rdeltat

             ! stratiform scavenging

             ! fracp is the fraction of cloud water converted to precip
             ! NB: In below formula for fracp cwat is a LWC/IWC that has already
             !     precipitated out, icol.e., cwat does not contain precipitation
             fracp(icol) = precs(icol,ilev)*deltat / &
                  max( 1.e-12_r8, cwat(icol,ilev) + precs(icol,ilev)*deltat )
             fracp(icol) = max( 0._r8, min( 1._r8, fracp(icol) ) )

             ! assume the corresponding amnt of tracer is removed
             st_scav_ic(icol) = sol_facti*clds(icol)*fracp(icol)*tracer(icol,ilev)*rdeltat

             odds(icol) = precabs(icol)/max(cldvst(icol,ilev),1.e-5_r8)*scavcoef(icol,ilev)*deltat
             odds(icol) = max(min(1._r8,odds(icol)),0._r8)
             st_scav_bc(icol) =sol_factb*(cldvst(icol,ilev)*odds(icol)) *tracer(icol,ilev)*rdeltat

          end if

          ! total convective scavenging
          srcc(icol) = conv_scav_ic(icol) + conv_scav_bc(icol)
          finc(icol) = conv_scav_ic(icol)/(srcc(icol) + 1.e-36_r8)

          ! total stratiform scavenging
          srcs(icol) = st_scav_ic(icol) + st_scav_bc(icol)
          fins(icol) = st_scav_ic(icol)/(srcs(icol) + 1.e-36_r8)

          ! make sure we dont take out more than is there
          ! ratio of amount available to amount removed
          rat(icol) = tracer(icol,ilev)/max(deltat*(srcc(icol)+srcs(icol)),1.e-36_r8)
          if (rat(icol)<1._r8) then
             srcs(icol) = srcs(icol)*rat(icol)
             srcc(icol) = srcc(icol)*rat(icol)
          endif
          srct(icol) = (srcc(icol)+srcs(icol))*omsm


          ! fraction that is not removed within the cloud
          ! (assumed to be interstitial, and subject to convective transport)
          fracp(icol) = deltat*srct(icol)/max(cldvst(icol,ilev)*tracer(icol,ilev),1.e-36_r8)  ! amount removed
          fracp(icol) = max(0._r8,min(1._r8,fracp(icol)))
          fracis(icol,ilev) = 1._r8 - fracp(icol)

          ! tend is all tracer removed by scavenging, plus all re-appearing from evaporation above
          ! Sungsu added cumulus contribution in the below 3 blocks
          scavt(icol,ilev) = -srct(icol) + (fracev(icol)*scavab(icol)+fracev_cu(icol)*scavabc(icol))*rpdog(icol)
          iscavt(icol,ilev) = -(srcc(icol)*finc(icol) + srcs(icol)*fins(icol))*omsm

          if ( present(icscavt) ) icscavt(icol,ilev) = -(srcc(icol)*finc(icol)) * omsm
          if ( present(isscavt) ) isscavt(icol,ilev) = -(srcs(icol)*fins(icol)) * omsm

          if (.not. out_resuspension) then
             if (present(bcscavt)) bcscavt(icol,ilev) = -(srcc(icol) * (1-finc(icol))) * omsm +  &
                  fracev_cu(icol)*scavabc(icol)*rpdog(icol)

             if (present(bsscavt)) bsscavt(icol,ilev) = -(srcs(icol) * (1-fins(icol))) * omsm +  &
                  fracev(icol)*scavab(icol)*rpdog(icol)
          else
             bcscavt(icol,ilev) = -(srcc(icol) * (1-finc(icol))) * omsm
             rcscavt(icol,ilev) =  fracev_cu(icol)*scavabc(icol)*rpdog(icol)

             bsscavt(icol,ilev) = -(srcs(icol) * (1-fins(icol))) * omsm
             rsscavt(icol,ilev) =  fracev(icol)*scavab(icol)*rpdog(icol)
          end if

          dblchek(icol) = tracer(icol,ilev) + deltat*scavt(icol,ilev)

          ! now keep track of scavenged mass and precip
          scavab(icol) = scavab(icol)*(1-fracev(icol)) + srcs(icol)*pdog(icol)
          precabs(icol) = precabs(icol) + (precs(icol,ilev) - evaps(icol,ilev))*pdog(icol)
          scavabc(icol) = scavabc(icol)*(1-fracev_cu(icol)) + srcc(icol)*pdog(icol)
          precabc(icol) = precabc(icol) + (cmfdqr(icol,ilev) - evapc(icol,ilev))*pdog(icol)

       end do ! End of icol = 1, ncol

       found = .false.
       do icol = 1,ncol
          if ( dblchek(icol) < 0._r8 ) then
             found = .true.
             exit
          end if
       end do

       if ( found ) then
          do icol = 1,ncol
             if (dblchek(icol) < 0._r8) then
                write(iulog,*) ' wetdapa: negative value ', icol, ilev, tracer(icol,ilev), &
                     dblchek(icol), scavt(icol,ilev), srct(icol), rat(icol), fracev(icol)
             endif
          end do
       endif

    end do ! End of ilev = 1, pver

  end subroutine wetdepa_v2

  !==============================================================================
  subroutine wetdepg( t, p, q, pdel, &
       cldt, cldc, cmfdqr, evapc, precs, evaps, &
       rain, cwat, tracer, deltat, molwt, &
       solconst, scavt, iscavt, cldv, icwmr1, &
       icwmr2, fracis, ncol )

    !-----------------------------------------------------------------------
    ! scavenging of gas phase constituents by henry's law ( Author: P. Rasch)
    !-----------------------------------------------------------------------

    real(r8), intent(in) ::&
         t(pcols,pver),        &! temperature
         p(pcols,pver),        &! pressure
         q(pcols,pver),        &! moisture
         pdel(pcols,pver),     &! pressure thikness
         cldt(pcols,pver),     &! total cloud fraction
         cldc(pcols,pver),     &! convective cloud fraction
         cmfdqr(pcols,pver),   &! rate of production of convective precip
         rain (pcols,pver),    &! total rainwater mixing ratio
         cwat(pcols,pver),     &! cloud water amount
         precs(pcols,pver),    &! rate of production of stratiform precip
         evaps(pcols,pver),    &! rate of evaporation of precip
         evapc(pcols,pver),    &! Rate of evaporation of convective precipitation
         cldv(pcols,pver),     &! estimate of local volume occupied by clouds
         icwmr1 (pcols,pver),  &! in cloud water mixing ration for zhang scheme
         icwmr2 (pcols,pver),  &! in cloud water mixing ration for hack  scheme
         deltat,               &! time step
         tracer(pcols,pver),   &! trace species
         molwt                  ! molecular weights

    integer, intent(in) :: ncol

    real(r8) &
         solconst(pcols,pver)   ! Henry's law coefficient

    real(r8), intent(out) ::&
         scavt(pcols,pver),    &! scavenging tend
         iscavt(pcols,pver),   &! incloud scavenging tends
         fracis(pcols, pver)    ! fraction of constituent that is insoluble

    ! local variables
    integer icol               ! x index
    integer ilev               ! z index
    real(r8) adjfac         ! factor stolen from cmfmca
    real(r8) aqfrac         ! fraction of tracer in aqueous phase
    real(r8) cwatc          ! local convective total water amount
    real(r8) cwats          ! local stratiform total water amount
    real(r8) cwatl          ! local cloud liq water amount
    real(r8) cwatp          ! local water amount falling from above precip
    real(r8) cwatpl         ! local water amount falling from above precip (liq)
    real(r8) cwatt          ! local sum of strat + conv total water amount
    real(r8) cwatti         ! cwatt/cldv = cloudy grid volume mixing ratio
    real(r8) fracev         ! fraction of precip from above that is evaporating
    real(r8) fracp          ! fraction of cloud water converted to precip
    real(r8) gafrac         ! fraction of tracer in gas phasea
    real(r8) hconst         ! henry's law solubility constant when equation is expressed
                            ! in terms of mixing ratios
    real(r8) mpla           ! moles / liter H2O entering the layer from above
    real(r8) mplb           ! moles / liter H2O leaving the layer below
    real(r8) omsm           ! 1 - (a small number)
    real(r8) part           !  partial pressure of tracer in atmospheres
    real(r8) patm           ! total pressure in atmospheres
    real(r8) pdog           ! work variable (pdel/gravit)
    real(r8) precab(pcols)  ! precip from above (work array)
    real(r8) precbl         ! precip work variable
    real(r8) precxx         ! precip work variable
    real(r8) precxx2        !
    real(r8) precic         ! precip work variable
    real(r8) rat            ! ratio of amount available to amount removed
    real(r8) scavab(pcols)  ! scavenged tracer flux from above (work array)
    real(r8) scavabc(pcols) ! scavenged tracer flux from above (work array)

    real(r8) scavmax        ! an estimate of the max tracer avail for removal
    real(r8) scavbl         ! flux removed at bottom of layer
    real(r8) fins           ! in cloud fraction removed by strat rain
    real(r8) finc           ! in cloud fraction removed by conv rain
    real(r8) rate           ! max removal rate estimate
    real(r8) scavlimt       ! limiting value 1
    real(r8) scavt1         ! limiting value 2
    real(r8) scavin         ! scavenging by incloud processes
    real(r8) scavbc         ! scavenging by below cloud processes
    real(r8) tc
    real(r8) weight         ! ice fraction
    real(r8) wtpl           ! work variable
    real(r8) cldmabs(pcols) ! maximum cloud at or above this level
    real(r8) cldmabc(pcols) ! maximum cloud at or above this level
    !-----------------------------------------------------------

    omsm = 1._r8-2*epsilon(1._r8)   ! used to prevent roundoff errors below zero

    adjfac = deltat/(max(deltat,cmftau)) ! adjustment factor from hack scheme

    ! zero accumulators
    do icol = 1,pcols
       precab(icol) = 1.e-36_r8
       scavab(icol) = 0._r8
       cldmabs(icol) = 0._r8
    end do

    do ilev = 1,pver
       do icol = 1,ncol
          tc = t(icol,ilev) - tmelt
          weight = max(0._r8,min(-tc*0.05_r8,1.0_r8)) ! fraction of condensate that is ice

          cldmabs(icol) = max(cldmabs(icol),cldt(icol,ilev))

          ! partitioning coefs for gas and aqueous phase
          ! take as a cloud water amount, the sum of the stratiform amount
          ! plus the convective rain water amount

          ! convective amnt is just the local precip rate from the hack scheme
          ! since there is no storage of water, this ignores that falling from above
          cwatc = (icwmr1(icol,ilev) + icwmr2(icol,ilev)) * (1._r8-weight)

          ! strat cloud water amount and also ignore the part falling from above
          cwats = cwat(icol,ilev)

          ! cloud water as liq
          cwatl = (1._r8-weight)*cwats

          ! cloud water as ice total suspended condensate as liquid
          cwatt = cwatl + rain(icol,ilev)

          ! incloud version
          cwatti = cwatt/max(cldv(icol,ilev), 0.00001_r8) + cwatc

          ! partitioning terms
          patm = p(icol,ilev)/1.013e5_r8 ! pressure in atmospheres
          hconst = molwta*patm*solconst(icol,ilev)*cwatti/rhoh2o
          aqfrac = hconst/(1._r8+hconst)
          gafrac = 1/(1._r8+hconst)
          fracis(icol,ilev) = gafrac

          ! partial pressure of the tracer in the gridbox in atmospheres
          part = patm*gafrac*tracer(icol,ilev)*molwta/molwt

          ! use henrys law to give moles tracer /liter of water  in this volume
          ! then convert to kg tracer /liter of water (kg tracer / kg water)
          mplb = solconst(icol,ilev)*part*molwt/1000._r8

          pdog = pdel(icol,ilev)/gravit

          ! this part of precip will be carried downward but at a new molarity of mpl
          precic = pdog*(precs(icol,ilev) + cmfdqr(icol,ilev))

          ! we cant take out more than entered, plus that available in the cloud
          scavmax = scavab(icol)+tracer(icol,ilev)*cldv(icol,ilev)/deltat*pdog

          ! flux of tracer by incloud processes
          scavin = precic*(1._r8-weight)*mplb

          ! fraction of precip which entered above that leaves below
          if (cam_physpkg_is('cam5') .or. cam_physpkg_is('cam6')) then
             ! Sungsu added evaporation of convective precipitation below.
             precxx = precab(icol)-pdog*(evaps(icol,ilev)+evapc(icol,ilev))
          else
             precxx = precab(icol)-pdog*evaps(icol,ilev)
          end if
          precxx = max (precxx,0.0_r8)

          ! flux of tracer by below cloud processes
          if (tc > 0) then
             scavbc = precxx*mplb ! if liquid
          else
             precxx2=max(precxx,1.e-36_r8)
             scavbc = scavab(icol)*precxx2/(precab(icol)) ! if ice
          endif

          scavbl = min(scavbc + scavin, scavmax)

          ! first guess assuming that henries law works
          scavt1 = (scavab(icol)-scavbl)/pdog*omsm

          ! pjr this should not be required, but we put it in to make sure we cant remove too much
          ! remember, scavt1 is generally negative (indicating removal)
          scavt1 = max(scavt1,-tracer(icol,ilev)*cldv(icol,ilev)/deltat)

          ! instead just set scavt to scavt1
          scavt(icol,ilev) = scavt1

          ! now update the amount leaving the layer
          scavbl = scavab(icol) - scavt(icol,ilev)*pdog

          ! in cloud amount is that formed locally over the total flux out bottom
          fins = scavin/(scavin + scavbc + 1.e-36_r8)
          iscavt(icol,ilev) = scavt(icol,ilev)*fins

          scavab(icol) = scavbl
          precab(icol) = max(precxx + precic,1.e-36_r8)

       end do
    end do

  end subroutine wetdepg

end module oslo_aero_depos
