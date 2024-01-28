module oslo_aero_nucleate_ice

  !---------------------------------------------------------------------------------
  !  A parameterization of ice nucleation.
  !
  ! Method:
  !  The current method is based on Liu & Penner (2005) & Liu et al. (2007)
  !  It related the ice nucleation with the aerosol number, temperature and the
  !  updraft velocity. It includes homogeneous freezing of sulfate & immersion
  !  freezing on mineral dust (soot disabled) in cirrus clouds, and
  !  Meyers et al. (1992) deposition nucleation in mixed-phase clouds
  !
  !  The effect of preexisting ice crystals on ice nucleation in cirrus clouds is included,
  !  and also consider the sub-grid variability of temperature in cirrus clouds,
  !  following X. Shi et al. ACP (2014).
  !
  !  Ice nucleation in mixed-phase clouds now uses classical nucleation theory (CNT),
  !  follows Y. Wang et al. ACP (2014), Hoose et al. (2010).
  !
  ! Authors:
  !  Xiaohong Liu, 01/2005, modifications by A. Gettelman 2009-2010
  !  Xiangjun Shi & Xiaohong Liu, 01/2014.
  !
  !  With help from C. C. Chen and B. Eaton (2014)
  !---------------------------------------------------------------------------------

  use shr_kind_mod,      only: r8=>shr_kind_r8
  use spmd_utils,        only: masterproc
  use ppgrid,            only: pcols, pver
  use constituents,      only: pcnst, cnst_get_ind
  use physconst,         only: pi, rair, tmelt
  use phys_control,      only: phys_getopts, use_hetfrz_classnuc, cam_physpkg_is
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
  use physics_buffer,    only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
  use cam_history,       only: addfld, add_default, outfld
  use ref_pres,          only: top_lev => trop_cloud_top_lev
  use wv_saturation,     only: qsat_water, svp_water, svp_ice
  use tropopause,        only: tropopause_findChemTrop
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun
  use time_manager,      only: is_first_step
  use nucleate_ice,      only: nucleati_init, nucleati  ! portable module
  !
  use oslo_aero_share,   only: l_dst_a2, l_dst_a3, MODE_IDX_DST_A2, MODE_IDX_DST_A3, rhopart, qqcw_get_field
  use oslo_aero_share,   only: MODE_IDX_DST_A2, MODE_IDX_DST_A3, MODE_IDX_SO4_AC,MODE_IDX_OMBC_INTMIX_COAT_AIT
  use oslo_aero_share,   only: volumeToNumber
  use oslo_aero_share,   only: nmodes

  implicit none
  private

  public :: nucleate_ice_oslo_readnl
  public :: nucleate_ice_oslo_register
  public :: nucleate_ice_oslo_init
  public :: nucleate_ice_oslo_calc

  private :: nucleati

  ! Namelist variables
  logical, public, protected :: use_preexisting_ice = .false.
  logical  :: hist_preexisting_ice = .false.
  logical  :: nucleate_ice_incloud = .false.
  logical  :: nucleate_ice_use_troplev = .false.
  real(r8) :: nucleate_ice_subgrid = -1._r8
  real(r8) :: nucleate_ice_subgrid_strat = -1._r8
  real(r8) :: nucleate_ice_strat = 0.0_r8

  ! Vars set via init method.
  real(r8) :: mincld      ! minimum allowed cloud fraction
  real(r8) :: bulk_scale  ! prescribed aerosol bulk sulfur scale factor

  logical  :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
  logical  :: use_incloud_nuc
  real(r8) :: ci

  ! constituent indices
  integer :: cldliq_idx   = -1
  integer :: cldice_idx   = -1
  integer :: numice_idx   = -1
  integer :: naai_idx     = -1
  integer :: naai_hom_idx = -1
  integer :: ast_idx      = -1
  integer :: qsatfac_idx  = -1

  real(r8), parameter :: Shet   = 1.3_r8     ! het freezing threshold
  real(r8), parameter :: rhoice = 0.5e3_r8   ! kg/m3, Wpice is not sensitive to rhoice
  real(r8), parameter :: minweff= 0.001_r8   ! m/s
  real(r8), parameter :: gamma4=6.0_r8

!===============================================================================
contains
!===============================================================================

  subroutine nucleate_ice_oslo_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use spmd_utils,      only: mpicom, masterprocid, mpi_logical, mpi_real8

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'nucleate_ice_cam_readnl'

    namelist /nucleate_ice_nl/ use_preexisting_ice, hist_preexisting_ice, &
         nucleate_ice_subgrid, nucleate_ice_subgrid_strat, nucleate_ice_strat, &
         nucleate_ice_incloud, nucleate_ice_use_troplev

    !-----------------------------------------------------------------------------

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'nucleate_ice_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, nucleate_ice_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    ! Broadcast namelist variables
    call mpi_bcast(use_preexisting_ice,  1, mpi_logical,masterprocid, mpicom, ierr)
    call mpi_bcast(hist_preexisting_ice, 1, mpi_logical,masterprocid, mpicom, ierr)
    call mpi_bcast(nucleate_ice_subgrid, 1, mpi_real8,  masterprocid, mpicom, ierr)
    call mpi_bcast(nucleate_ice_subgrid_strat, 1, mpi_real8,  masterprocid, mpicom, ierr)
    call mpi_bcast(nucleate_ice_strat,   1, mpi_real8,  masterprocid, mpicom, ierr)
    call mpi_bcast(nucleate_ice_incloud, 1, mpi_logical,masterprocid, mpicom, ierr)
    call mpi_bcast(nucleate_ice_use_troplev, 1, mpi_logical,masterprocid, mpicom, ierr)

    ! Set module variable
    use_incloud_nuc = nucleate_ice_incloud

  end subroutine nucleate_ice_oslo_readnl

  !================================================================================================

  subroutine nucleate_ice_oslo_register()

    call pbuf_add_field('NAAI',     'global', dtype_r8, (/pcols,pver/), naai_idx)
    call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

  end subroutine nucleate_ice_oslo_register

  !================================================================================================

  subroutine nucleate_ice_oslo_init(mincld_in, bulk_scale_in)

    ! arguments
    real(r8)                  , intent(in) :: mincld_in
    real(r8)                  , intent(in) :: bulk_scale_in

    ! local variables
    integer :: ierr
    integer :: m, n
    logical :: history_cesm_forcing
    character(len=*), parameter :: routine = 'nucleate_ice_oslo_init'
    !--------------------------------------------------------------------------------------------

    call phys_getopts(history_cesm_forcing_out = history_cesm_forcing)

    mincld     = mincld_in
    bulk_scale = bulk_scale_in

    ! TODO: IS this necessary?
    ! if (is_first_step()) then
    !    call pbuf_set_field(pbuf, naai_idx, 0.0_r8)
    ! end if

    ! Initialize naai.

    if( masterproc ) then
       write(iulog,*) 'nucleate_ice parameters:'
       write(iulog,*) '  mincld                     = ', mincld_in
       write(iulog,*) '  bulk_scale                 = ', bulk_scale_in
       write(iulog,*) '  use_preexisiting_ice       = ', use_preexisting_ice
       write(iulog,*) '  hist_preexisiting_ice      = ', hist_preexisting_ice
       write(iulog,*) '  nucleate_ice_subgrid       = ', nucleate_ice_subgrid
       write(iulog,*) '  nucleate_ice_subgrid_strat = ', nucleate_ice_subgrid_strat
       write(iulog,*) '  nucleate_ice_strat         = ', nucleate_ice_strat
       write(iulog,*) '  nucleate_ice_incloud       = ', nucleate_ice_incloud
       write(iulog,*) '  nucleate_ice_use_troplev   = ', nucleate_ice_use_troplev
    end if

    call cnst_get_ind('CLDLIQ', cldliq_idx)
    call cnst_get_ind('CLDICE', cldice_idx)
    call cnst_get_ind('NUMICE', numice_idx)
    qsatfac_idx = pbuf_get_index('QSATFAC', ierr)

    if (((nucleate_ice_subgrid .eq. -1._r8) .or. (nucleate_ice_subgrid_strat .eq. -1._r8)) .and. (qsatfac_idx .eq. -1)) then
       call endrun(routine//': ERROR qsatfac is required when subgrid = -1 or subgrid_strat = -1')
    end if

    if (cam_physpkg_is("cam_dev")) then
       ! Updates for PUMAS v1.21+
       call addfld('NIHFTEN',  (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to homogenous freezing')
       call addfld('NIDEPTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to deposition nucleation')
       call addfld('NIIMMTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to immersion freezing')
       call addfld('NIMEYTEN', (/ 'lev' /), 'A', '1/m3/s', 'Activated Ice Number Concentration tendency due to meyers deposition')
    else
       call addfld('NIHF',  (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to homogenous freezing')
       call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to deposition nucleation')
       call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to immersion freezing')
       call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentration due to meyers deposition')
    endif

    call addfld('NIREGM',(/ 'lev' /), 'A', 'C', 'Ice Nucleation Temperature Threshold for Regime')
    call addfld('NISUBGRID',(/ 'lev' /), 'A', '', 'Ice Nucleation subgrid saturation factor')
    call addfld('NITROP_PD',(/ 'lev' /), 'A', '', 'Chemical Tropopause probability')
    if ( history_cesm_forcing ) then
       call add_default('NITROP_PD',8,' ')
    endif

    if (use_preexisting_ice) then
       call addfld('fhom',      (/ 'lev' /), 'A','fraction', 'Fraction of cirrus where homogeneous freezing occur'   )
       call addfld ('WICE',     (/ 'lev' /), 'A','m/s','Vertical velocity Reduction caused by preexisting ice'  )
       call addfld ('WEFF',     (/ 'lev' /), 'A','m/s','Effective Vertical velocity for ice nucleation' )

       if (cam_physpkg_is("cam_dev")) then
          ! Updates for PUMAS v1.21+
          call addfld ('INnso4TEN',   (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency so4 (in) to ice_nucleation')
          call addfld ('INnbcTEN',    (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency bc  (in) to ice_nucleation')
          call addfld ('INndustTEN',  (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency dust (in) ice_nucleation')
          call addfld ('INondustTEN',  (/ 'lev' /), 'A','1/m3/s','Number Concentration tendency dust (out) from ice_nucleation')
          call addfld ('INhetTEN',    (/ 'lev' /), 'A','1/m3/s', &
               'Tendency for contribution for in-cloud ice number density increase by het nucleation in ice cloud')
          call addfld ('INhomTEN',    (/ 'lev' /), 'A','1/m3/s', &
               'Tendency for contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
       else
          call addfld ('INnso4',   (/ 'lev' /), 'A','1/m3','Number Concentration so4 (in) to ice_nucleation')
          call addfld ('INnbc',    (/ 'lev' /), 'A','1/m3','Number Concentration bc  (in) to ice_nucleation')
          call addfld ('INndust',  (/ 'lev' /), 'A','1/m3','Number Concentration dust (in) ice_nucleation')
          call addfld ('INondust',  (/ 'lev' /), 'A','1/m3','Number Concentration dust (out) from ice_nucleation')
          call addfld ('INhet',    (/ 'lev' /), 'A','1/m3', &
               'contribution for in-cloud ice number density increase by het nucleation in ice cloud')
          call addfld ('INhom',    (/ 'lev' /), 'A','1/m3', &
               'contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
       endif

       call addfld ('INFrehom', (/ 'lev' /), 'A','frequency','hom IN frequency ice cloud')
       call addfld ('INFreIN',  (/ 'lev' /), 'A','frequency','frequency of ice nucleation occur')

       if (hist_preexisting_ice) then
          call add_default ('WSUBI   ', 1, ' ')  ! addfld/outfld calls are in microp_aero

          call add_default ('fhom    ', 1, ' ')
          call add_default ('WICE    ', 1, ' ')
          call add_default ('WEFF    ', 1, ' ')
          call add_default ('INnso4  ', 1, ' ')
          call add_default ('INnbc   ', 1, ' ')
          call add_default ('INndust ', 1, ' ')
          call add_default ('INhet   ', 1, ' ')
          call add_default ('INhom   ', 1, ' ')
          call add_default ('INFrehom', 1, ' ')
          call add_default ('INFreIN ', 1, ' ')
       end if
    end if

    lq(l_dst_a2) = .true.
    lq(l_dst_a3) = .true.

    call nucleati_init(use_preexisting_ice, use_hetfrz_classnuc, nucleate_ice_incloud, iulog, pi, &
         mincld)

    ! get indices for fields in the physics buffer
    ast_idx  = pbuf_get_index('AST')

    ci = rhoice*pi/6._r8

  end subroutine nucleate_ice_oslo_init

  !================================================================================================

  subroutine nucleate_ice_oslo_calc( state, wsubi, pbuf, dtime, ptend, numberConcentration)

    ! arguments
    type(physics_state), target, intent(in)  :: state
    real(r8),                    intent(in)  :: wsubi(:,:)
    type(physics_buffer_desc),   pointer     :: pbuf(:)
    real(r8),                    intent(in)  :: dtime
    type(physics_ptend),         intent(out) :: ptend
    real(r8),                    intent(in)  :: numberConcentration(pcols,pver,0:nmodes)

    ! local workspace

    ! naai and naai_hom are the outputs shared with the microphysics
    real(r8), pointer :: naai(:,:)           ! number of activated aerosol for ice nucleation
    real(r8), pointer :: naai_hom(:,:)       ! number of activated aerosol for ice nucleation (homogeneous freezing only)

    integer :: lchnk, ncol
    integer :: itim_old
    integer :: i, k, m

    real(r8), pointer :: t(:,:)              ! input temperature (K)
    real(r8), pointer :: qn(:,:)             ! input water vapor mixing ratio (kg/kg)
    real(r8), pointer :: qc(:,:)             ! cloud water mixing ratio (kg/kg)
    real(r8), pointer :: qi(:,:)             ! cloud ice mixing ratio (kg/kg)
    real(r8), pointer :: ni(:,:)             ! cloud ice number conc (1/kg)
    real(r8), pointer :: pmid(:,:)           ! pressure at layer midpoints (pa)

    real(r8), pointer :: ast(:,:)
    real(r8)          :: icecldf(pcols,pver) ! ice cloud fraction
    real(r8), pointer :: qsatfac(:,:)        ! Subgrid cloud water saturation scaling factor.

    real(r8), pointer :: cld_dst_a2(:,:)     ! mmr cld dst a2
    real(r8), pointer :: cld_dst_a3(:,:)     ! mass m.r. of coarse dust

    real(r8) :: rho(pcols,pver)              ! air density (kg m-3)

    real(r8) :: qs(pcols)                    ! liquid-ice weighted sat mixing rat (kg/kg)
    real(r8) :: es(pcols)                    ! liquid-ice weighted sat vapor press (pa)
    real(r8) :: gammas(pcols)                ! parameter for cond/evap of cloud water
    integer  :: troplev(pcols)               ! tropopause level

    real(r8) :: relhum(pcols,pver)           ! relative humidity
    real(r8) :: icldm(pcols,pver)            ! ice cloud fraction

    real(r8) :: dst_num                      ! total dust aerosol number (#/cm^3)
    real(r8) :: dso4_num                     ! tuning factor for increased so4
    real(r8) :: so4_num                      ! so4 aerosol number (#/cm^3)
    real(r8) :: soot_num                     ! soot (hydrophilic) aerosol number (#/cm^3)
    real(r8) :: wght
    real(r8) :: oso4_num
    real(r8) :: odst_num
    real(r8) :: osoot_num
    real(r8) :: so4_num_st_cr_tot
    real(r8) :: ramp

    real(r8) :: dmc
    real(r8) :: ssmc
    real(r8) :: dust_coarse_fraction                  ! fraction of dust in coarse (a3) mode
    real(r8) :: masslost                              ! [kg/kg] tmp variable for mass lost
    real(r8) :: numberFromSmallDustMode               ! [#/cm3] number of dust activated from small mode

    real(r8) :: subgrid(pcols,pver)
    real(r8) :: trop_pd(pcols,pver)

    ! For pre-existing ice
    real(r8) :: fhom(pcols,pver)     ! how much fraction of cloud can reach Shom
    real(r8) :: wice(pcols,pver)     ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
    real(r8) :: weff(pcols,pver)     ! effective Vertical velocity for ice nucleation (m/s); weff=wsubi-wice
    real(r8) :: INnso4(pcols,pver)   ! #/m3, so4 aerosol number used for ice nucleation
    real(r8) :: INnbc(pcols,pver)    ! #/m3, bc aerosol number used for ice nucleation
    real(r8) :: INndust(pcols,pver)  ! #/m3, dust aerosol number used for ice nucleation
    real(r8) :: INondust(pcols,pver) ! #/m3, dust aerosol number used for ice nucleation
    real(r8) :: INhet(pcols,pver)    ! #/m3, ice number from het freezing
    real(r8) :: INhom(pcols,pver)    ! #/m3, ice number from hom freezing
    real(r8) :: INFrehom(pcols,pver) ! hom freezing occurence frequency.  1 occur, 0 not occur.
    real(r8) :: INFreIN(pcols,pver)  ! ice nucleation occerence frequency.   1 occur, 0 not occur.

    ! history output for ice nucleation
    real(r8) :: nihf(pcols,pver)  !output number conc of ice nuclei due to heterogenous freezing (1/m3)
    real(r8) :: niimm(pcols,pver) !output number conc of ice nuclei due to immersion freezing (hetero nuc) (1/m3)
    real(r8) :: nidep(pcols,pver) !output number conc of ice nuclei due to deoposion nucleation (hetero nuc) (1/m3)
    real(r8) :: nimey(pcols,pver) !output number conc of ice nuclei due to meyers deposition (1/m3)
    real(r8) :: regm(pcols,pver)  !output temperature thershold for nucleation regime

    real(r8) :: so4_num_ac
    real(r8) :: so4_num_cr
    !-------------------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    t     => state%t
    qn    => state%q(:,:,1)
    qc    => state%q(:,:,cldliq_idx)
    qi    => state%q(:,:,cldice_idx)
    ni    => state%q(:,:,numice_idx)
    pmid  => state%pmid

    rho(:ncol,:) = pmid(:ncol,:)/(rair*t(:ncol,:))

    call physics_ptend_init(ptend, state%psetcols, 'nucleatei', lq=lq)

    cld_dst_a2 => qqcw_get_field(pbuf, l_dst_a2)
    cld_dst_a3 => qqcw_get_field(pbuf, l_dst_a2)

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

    icecldf(:ncol,:pver) = ast(:ncol,:pver)

    ! naai and naai_hom are the outputs from this parameterization
    call pbuf_get_field(pbuf, naai_idx, naai)
    call pbuf_get_field(pbuf, naai_hom_idx, naai_hom)
    naai(1:ncol,1:pver)     = 0._r8
    naai_hom(1:ncol,1:pver) = 0._r8

    ! Use the same criteria that is used in chemistry and in CLUBB (for cloud fraction)
    ! to determine whether to use tropospheric or stratospheric settings. Include the
    ! tropopause level so that the cold point tropopause will use the stratospheric values.
    call tropopause_findChemTrop(state, troplev)

    if ((nucleate_ice_subgrid .eq. -1._r8) .or. (nucleate_ice_subgrid_strat .eq. -1._r8)) then
       call pbuf_get_field(pbuf, qsatfac_idx, qsatfac)
    end if

    trop_pd(:,:) = 0._r8

    do k = top_lev, pver
       do i = 1, ncol
          trop_pd(i, troplev(i)) = 1._r8

          if (k <= troplev(i)) then
             if (nucleate_ice_subgrid_strat .eq. -1._r8) then
                subgrid(i, k) = 1._r8 / qsatfac(i, k)
             else
                subgrid(i, k) = nucleate_ice_subgrid_strat
             end if
          else
             if (nucleate_ice_subgrid .eq. -1._r8) then
                subgrid(i, k) = 1._r8 / qsatfac(i, k)
             else
                subgrid(i, k) = nucleate_ice_subgrid
             end if
          end if
       end do
    end do

    ! initialize history output fields for ice nucleation
    nihf(1:ncol,1:pver)  = 0._r8
    niimm(1:ncol,1:pver) = 0._r8
    nidep(1:ncol,1:pver) = 0._r8
    nimey(1:ncol,1:pver) = 0._r8

    regm(1:ncol,1:pver) = 0._r8

    if (use_preexisting_ice) then
       fhom(:,:)     = 0.0_r8
       wice(:,:)     = 0.0_r8
       weff(:,:)     = 0.0_r8
       INnso4(:,:)   = 0.0_r8
       INnbc(:,:)    = 0.0_r8
       INndust(:,:)  = 0.0_r8
       INondust(:,:) = 0.0_r8
       INhet(:,:)    = 0.0_r8
       INhom(:,:)    = 0.0_r8
       INFrehom(:,:) = 0.0_r8
       INFreIN(:,:)  = 0.0_r8
    endif

    do k = top_lev, pver
       ! Get humidity and saturation vapor pressures
       call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol), ncol, gam=gammas(:ncol))

       do i = 1, ncol
          relhum(i,k) = qn(i,k)/qs(i)
          icldm(i,k) = max(icecldf(i,k), mincld) ! get cloud fraction, check for minimum
       end do
    end do

    kloop: do k = top_lev, pver
       iloop: do i = 1, ncol

          so4_num_st_cr_tot = 0._r8
          freezing: if (t(i,k) < tmelt - 5._r8) then

             ! compute aerosol number for so4, soot, and dust with units #/cm^3
             ! soot = accumulation mode
             ! sulfate = aiken mode
             ! dust = coarse mode
             ! since modal has internal mixtures.
             ! Oslo aerosols have two modes.. Need mode-fractions
             so4_num = (numberConcentration(i,k,MODE_IDX_SO4_AC))*1.0e-6_r8
             dst_num = (numberConcentration(i,k,MODE_IDX_DST_A2) + numberConcentration(i,k,MODE_IDX_DST_A3))*1.0e-6_r8
             !soot_num =  numberConcentration(i,k,MODE_IDX_OMBC_INTMIX_COAT_AIT)*1.0e-6_r8
             dust_coarse_fraction = numberConcentration(i,k,MODE_IDX_DST_A3)*1.e-6_r8 / (dst_num+1.e-100_r8)

             ! *** Turn off soot nucleation ***
             soot_num = 0.0_r8

             if (cam_physpkg_is("cam_dev")) then
                call nucleati( &
                     wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
                     qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
                     so4_num, dst_num, soot_num, subgrid(i,k),                 &
                     naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                     wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
                     oso4_num, odst_num, osoot_num, &
                     call_frm_zm_in = .false., add_preexisting_ice_in = .false.)
             else
                call nucleati( &
                     wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
                     qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
                     so4_num, dst_num, soot_num, subgrid(i,k),                 &
                     naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                     wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
                     oso4_num, odst_num, osoot_num)
             end if

             ! Move aerosol used for nucleation from interstial to cloudborne,
             ! otherwise the same coarse mode aerosols will be available again
             ! in the next timestep and will supress homogeneous freezing.
             if (use_preexisting_ice) then

                numberFromSmallDustMode = 0.0_r8

                !Assume the coarse aerosols were activated first
                !so only remove small ones if more than large ones are activated
                if(odst_num .gt. dst_num*dust_coarse_fraction)then

                   !A2-mode
                   numberFromSmallDustMode = odst_num - dst_num*dust_coarse_fraction

                   masslost = (odst_num                  & !all removed
                        - dst_num*dust_coarse_fraction)  & !fraction to coarse mode
                        / volumeToNumber(MODE_IDX_DST_A2) &
                        * rhopart(l_dst_a2) &
                        /rho(i,k)*1e6_r8

                   ptend%q(i,k,l_dst_a2) = -masslost*icldm(i,k)/ dtime
                   cld_dst_a2(i,k) = cld_dst_a2(i,k) + masslost*icldm(i,k)

                end if

                ! Coarse mode (is always lost)
                masslost = (odst_num - numberFromSmallDustMode) &
                     / volumeToNumber(MODE_IDX_DST_A3)    &
                     * rhopart(l_dst_a3)                  &
                     / rho(i,k)*1e6_r8

                ptend%q(i,k,l_dst_a3) = -masslost * icldm(i,k) / dtime
                cld_dst_a3(i,k) = cld_dst_a3(i,k) + masslost*icldm(i,k)

             end if

             !Oslo aerosols do not have explicit treatment of coarse sulfate
             so4_num_cr = 0.0_r8

             ! Liu&Penner does not generate enough nucleation in the polar winter
             ! stratosphere, which affects surface area density, dehydration and
             ! ozone chemistry. Part of this is that there are a larger number of
             ! particles in the accumulation mode than in the Aitken mode. In volcanic
             ! periods, the coarse mode may also be important. As a short
             ! term work around, include the accumulation and coarse mode particles
             ! and assume a larger fraction of the sulfates nucleate in the polar
             ! stratosphere.
             !
             ! Do not include the tropopause level, as stratospheric aerosols
             ! only exist above the tropopause level.
             !
             ! NOTE: This may still not represent the proper particles that
             ! participate in nucleation, because it doesn't include STS and NAT
             ! particles. It may not represent the proper saturation threshold for
             ! nucleation, and wsubi from CLUBB is probably not representative of
             ! wave driven varaibility in the polar stratosphere.
             if (nucleate_ice_use_troplev) then
                if ((k < troplev(i)) .and. (nucleate_ice_strat > 0._r8) .and. (oso4_num > 0._r8)) then
                   so4_num_ac = so4_num*rho(i,k)*1.0e-6_r8 !This is maximum sulfate which can activate
                   dso4_num = max(0._r8, (nucleate_ice_strat * (so4_num_cr + so4_num_ac)) - oso4_num) * 1e6_r8 / rho(i,k)
                   naai(i,k) = naai(i,k) + dso4_num
                   nihf(i,k) = nihf(i,k) + dso4_num
                end if
             else
                ! This maintains backwards compatibility with the previous version.
                if (pmid(i,k) <= 12500._r8 .and. pmid(i,k) > 100._r8 .and. abs(state%lat(i)) >= 60._r8 * pi / 180._r8) then
                   ramp = 1._r8 - min(1._r8, max(0._r8, (pmid(i,k) - 10000._r8) / 2500._r8))

                   if (oso4_num > 0._r8) then
                      dso4_num = (max(oso4_num, ramp * nucleate_ice_strat * so4_num) - oso4_num) * 1e6_r8 / rho(i,k)
                      naai(i,k) = naai(i,k) + dso4_num
                      nihf(i,k) = nihf(i,k) + dso4_num
                   end if
                end if
             end if

             if (cam_physpkg_is("cam_dev")) then
                !Updates for pumas v1.21+

                naai_hom(i,k) = nihf(i,k)/dtime
                naai(i,k)= naai(i,k)/dtime

                ! output activated ice (convert from #/kg -> #/m3/s)
                nihf(i,k)  = nihf(i,k) *rho(i,k)/dtime
                niimm(i,k) = niimm(i,k)*rho(i,k)/dtime
                nidep(i,k) = nidep(i,k)*rho(i,k)/dtime
                nimey(i,k) = nimey(i,k)*rho(i,k)/dtime

                if (use_preexisting_ice) then
                   INnso4(i,k)   = so4_num*1e6_r8            ! (convert from #/cm3 -> #/m3)
                   INnbc(i,k)    = soot_num*1e6_r8
                   INndust(i,k)  = dst_num*1e6_r8
                   INondust(i,k) = odst_num*1e6_r8
                   INFreIN(i,k)  = 1.0_r8                    ! 1,ice nucleation occur
                   INhet(i,k)    = (niimm(i,k) + nidep(i,k)) ! #/m3, nimey not in cirrus
                   INhom(i,k)    = nihf(i,k)                 ! #/m3
                   if (INhom(i,k).gt.1e3_r8)   then          ! > 1/L
                      INFrehom(i,k)=1.0_r8                   ! 1, hom freezing occur
                   endif

                   ! exclude  no ice nucleaton
                   if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then
                      INnso4(i,k) =0.0_r8
                      INnbc(i,k)  =0.0_r8
                      INndust(i,k)=0.0_r8
                      INondust(i,k)=0.0_r8
                      INFreIN(i,k)=0.0_r8
                      INhet(i,k) = 0.0_r8
                      INhom(i,k) = 0.0_r8
                      INFrehom(i,k)=0.0_r8
                      wice(i,k) = 0.0_r8
                      weff(i,k) = 0.0_r8
                      fhom(i,k) = 0.0_r8
                   endif
                end if

             else ! Not cam_dev

                naai_hom(i,k) = nihf(i,k)

                ! output activated ice (convert from #/kg -> #/m3/s)
                nihf(i,k)     = nihf(i,k) *rho(i,k)
                niimm(i,k)    = niimm(i,k)*rho(i,k)
                nidep(i,k)    = nidep(i,k)*rho(i,k)
                nimey(i,k)    = nimey(i,k)*rho(i,k)

                if (use_preexisting_ice) then
                   INnso4(i,k) =so4_num*1e6_r8 ! (convert from #/cm3 -> #/m3/s)
                   INnbc(i,k)  =soot_num*1e6_r8
                   INndust(i,k)=dst_num*1e6_r8
                   INondust(i,k)=odst_num*1e6_r8
                   INFreIN(i,k)=1.0_r8          ! 1,ice nucleation occur
                   INhet(i,k) = (niimm(i,k) + nidep(i,k))   ! #/m3, nimey not in cirrus
                   INhom(i,k) = nihf(i,k)                 ! #/m3
                   if (INhom(i,k).gt.1e3_r8)   then ! > 1/L
                      INFrehom(i,k)=1.0_r8       ! 1, hom freezing occur
                   endif

                   ! exclude  no ice nucleaton
                   if ((INFrehom(i,k) < 0.5_r8) .and. (INhet(i,k) < 1.0_r8))   then
                      INnso4(i,k) =0.0_r8
                      INnbc(i,k)  =0.0_r8
                      INndust(i,k)=0.0_r8
                      INondust(i,k)=0.0_r8
                      INFreIN(i,k)=0.0_r8
                      INhet(i,k) = 0.0_r8
                      INhom(i,k) = 0.0_r8
                      INFrehom(i,k)=0.0_r8
                      wice(i,k) = 0.0_r8
                      weff(i,k) = 0.0_r8
                      fhom(i,k) = 0.0_r8
                   endif
                end if

             end if ! cam_dev

          end if freezing
       end do iloop
    end do kloop

    if (cam_physpkg_is("cam_dev")) then
       ! Updates for PUMAS v1.21+
       call outfld('NIHFTEN',   nihf(:ncol,:), ncol, lchnk)
       call outfld('NIIMMTEN', niimm(:ncol,:), ncol, lchnk)
       call outfld('NIDEPTEN', nidep(:ncol,:), ncol, lchnk)
       call outfld('NIMEYTEN', nimey(:ncol,:), ncol, lchnk)
    else
       call outfld('NIHF',   nihf(:ncol,:), ncol, lchnk)
       call outfld('NIIMM', niimm(:ncol,:), ncol, lchnk)
       call outfld('NIDEP', nidep(:ncol,:), ncol, lchnk)
       call outfld('NIMEY', nimey(:ncol,:), ncol, lchnk)
    end if
    call outfld('NIREGM',    regm(:ncol,:),    ncol, lchnk)
    call outfld('NISUBGRID', subgrid(:ncol,:), ncol, lchnk)
    call outfld('NITROP_PD', trop_pd(:ncol,:), ncol, lchnk)

    if (use_preexisting_ice) then
       call outfld('fhom', fhom(:ncol,:), ncol, lchnk)
       call outfld('WICE', wice(:ncol,:), ncol, lchnk)
       call outfld('WEFF', weff(:ncol,:), ncol, lchnk)
       if (cam_physpkg_is("cam_dev")) then
          ! Updates for PUMAS v1.21+
          call outfld('INnso4TEN',     INnso4(:ncol,:), ncol, lchnk)
          call outfld('INnbcTEN',       INnbc(:ncol,:), ncol, lchnk)
          call outfld('INndustTEN',   INndust(:ncol,:), ncol, lchnk)
          call outfld('INondustTEN', INondust(:ncol,:), ncol, lchnk)
          call outfld('INhetTEN',       INhet(:ncol,:), ncol, lchnk)
          call outfld('INhomTEN',       INhom(:ncol,:), ncol, lchnk)
       else
          call outfld('INnso4  ',    INnso4(:ncol,:), ncol, lchnk)
          call outfld('INnbc   ',     INnbc(:ncol,:), ncol, lchnk)
          call outfld('INndust ',   INndust(:ncol,:), ncol, lchnk)
          call outfld('INondust ', INondust(:ncol,:), ncol, lchnk)
          call outfld('INhet   ',     INhet(:ncol,:), ncol, lchnk)
          call outfld('INhom   ',     INhom(:ncol,:), ncol, lchnk)
       end if
       call outfld('INFrehom', INFrehom(:ncol,:), ncol, lchnk)
       call outfld('INFreIN ',  INFreIN(:ncol,:), ncol, lchnk)
    end if

  end subroutine nucleate_ice_oslo_calc

end module oslo_aero_nucleate_ice
