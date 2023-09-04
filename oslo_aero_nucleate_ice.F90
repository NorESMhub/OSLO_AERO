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
  use phys_control,      only: phys_getopts, use_hetfrz_classnuc
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
  use physics_buffer,    only: pbuf_add_field, dtype_r8, pbuf_old_tim_idx, pbuf_get_index, pbuf_get_field
  use cam_history,       only: addfld, add_default, outfld
  use ref_pres,          only: top_lev => trop_cloud_top_lev
  use wv_saturation,     only: qsat_water, svp_water, svp_ice
  use tropopause,        only: tropopause_findChemTrop
  use cam_logfile,       only: iulog
  use cam_abortutils,    only: endrun
  !
  use oslo_aero_share,   only: l_dst_a2, l_dst_a3, MODE_IDX_DST_A2, MODE_IDX_DST_A3, rhopart, qqcw_get_field
  use oslo_aero_share,   only: MODE_IDX_DST_A2, MODE_IDX_DST_A3, MODE_IDX_SO4_AC,MODE_IDX_OMBC_INTMIX_COAT_AIT 
  use oslo_aero_const,   only: volumeToNumber
  use oslo_aero_params,  only: nmodes

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

  logical :: clim_modal_aero = .true.
  logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
  logical :: use_incloud_nuc
  real(r8) :: ci

  ! constituent indices
  integer :: &
       cldliq_idx = -1, &
       cldice_idx = -1, &
       numice_idx = -1

  integer :: &
       naai_idx,     &
       naai_hom_idx

  integer :: &
       ast_idx   = -1

  integer :: &
       qsatfac_idx

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
    use_incloud_nuc =  nucleate_ice_incloud

  end subroutine nucleate_ice_oslo_readnl

  !================================================================================================

  subroutine nucleate_ice_oslo_register()

    call pbuf_add_field('NAAI',     'physpkg', dtype_r8, (/pcols,pver/), naai_idx)
    call pbuf_add_field('NAAI_HOM', 'physpkg', dtype_r8, (/pcols,pver/), naai_hom_idx)

  end subroutine nucleate_ice_oslo_register

  !================================================================================================

  subroutine nucleate_ice_oslo_init(mincld_in, bulk_scale_in)

    ! arguments
    real(r8), intent(in) :: mincld_in
    real(r8), intent(in) :: bulk_scale_in

    ! local variables
    integer :: ierr
    integer :: m, n
    logical :: history_cesm_forcing
    character(len=*), parameter :: routine = 'nucleate_ice_cam_init'
    !--------------------------------------------------------------------------------------------

    call phys_getopts(history_cesm_forcing_out = history_cesm_forcing)

    mincld     = mincld_in
    bulk_scale = bulk_scale_in

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

    call addfld('NIHF',  (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to homogenous freezing')
    call addfld('NIDEP', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to deposition nucleation')
    call addfld('NIIMM', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to immersion freezing')
    call addfld('NIMEY', (/ 'lev' /), 'A', '1/m3', 'Activated Ice Number Concentation due to meyers deposition')

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
       call addfld ('INnso4',   (/ 'lev' /), 'A','1/m3','Number Concentation so4 (in) to ice_nucleation')
       call addfld ('INnbc',    (/ 'lev' /), 'A','1/m3','Number Concentation bc  (in) to ice_nucleation')
       call addfld ('INndust',  (/ 'lev' /), 'A','1/m3','Number Concentation dust (in) ice_nucleation')
       call addfld ('INondust',  (/ 'lev' /), 'A','1/m3','Number Concentation dust (out) from ice_nucleation')
       call addfld ('INhet',    (/ 'lev' /), 'A','1/m3', &
            'contribution for in-cloud ice number density increase by het nucleation in ice cloud')
       call addfld ('INhom',    (/ 'lev' /), 'A','1/m3', &
            'contribution for in-cloud ice number density increase by hom nucleation in ice cloud')
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

    lq(l_dst_a2) = .TRUE.
    lq(l_dst_a3) = .TRUE.

    ! get indices for fields in the physics buffer
    ast_idx  = pbuf_get_index('AST')

    ci = rhoice*pi/6._r8

  end subroutine nucleate_ice_oslo_init

  !================================================================================================

  subroutine nucleate_ice_oslo_calc( state, wsubi, pbuf, dtime, ptend, numberConcentration)

    ! arguments
    real(r8), intent(in)                       :: numberConcentration(pcols,pver,0:nmodes)
    type(physics_state), target, intent(in)    :: state
    real(r8),                    intent(in)    :: wsubi(:,:)
    type(physics_buffer_desc),   pointer       :: pbuf(:)
    real(r8),                    intent(in)    :: dtime
    type(physics_ptend),         intent(out)   :: ptend

    ! local workspace

    ! naai and naai_hom are the outputs shared with the microphysics
    real(r8), pointer :: naai(:,:)       ! number of activated aerosol for ice nucleation 
    real(r8), pointer :: naai_hom(:,:)   ! number of activated aerosol for ice nucleation (homogeneous freezing only)

    integer :: lchnk, ncol
    integer :: itim_old
    integer :: i, k, m

    real(r8), pointer :: t(:,:)          ! input temperature (K)
    real(r8), pointer :: qn(:,:)         ! input water vapor mixing ratio (kg/kg)
    real(r8), pointer :: qc(:,:)         ! cloud water mixing ratio (kg/kg)
    real(r8), pointer :: qi(:,:)         ! cloud ice mixing ratio (kg/kg)
    real(r8), pointer :: ni(:,:)         ! cloud ice number conc (1/kg)
    real(r8), pointer :: pmid(:,:)       ! pressure at layer midpoints (pa)
    real(r8), pointer :: cld_dst_a2(:,:) ! mmr cld dst a2
    real(r8), pointer :: cld_dst_a3(:,:) ! mass m.r. of coarse dust
    real(r8), pointer :: ast(:,:)
    real(r8), pointer :: qsatfac(:,:)      ! Subgrid cloud water saturation scaling factor.

    real(r8) :: icecldf(pcols,pver)  ! ice cloud fraction
    real(r8) :: rho(pcols,pver)      ! air density (kg m-3)
    real(r8) :: qs(pcols)            ! liquid-ice weighted sat mixing rat (kg/kg)
    real(r8) :: es(pcols)            ! liquid-ice weighted sat vapor press (pa)
    real(r8) :: gammas(pcols)        ! parameter for cond/evap of cloud water
    integer  :: troplev(pcols)       ! tropopause level

    real(r8) :: relhum(pcols,pver)  ! relative humidity
    real(r8) :: icldm(pcols,pver)   ! ice cloud fraction

    real(r8) :: so4_num                               ! so4 aerosol number (#/cm^3)
    real(r8) :: soot_num                              ! soot (hydrophilic) aerosol number (#/cm^3)
    real(r8) :: dst1_num,dst2_num,dst3_num,dst4_num   ! dust aerosol number (#/cm^3)
    real(r8) :: dst_num                               ! total dust aerosol number (#/cm^3)
    real(r8) :: wght
    real(r8) :: dmc
    real(r8) :: ssmc
    real(r8) :: oso4_num
    real(r8) :: odst_num
    real(r8) :: osoot_num
    real(r8) :: dso4_num                              ! tuning factor for increased so4
    real(r8) :: ramp                                  ! ---------- " ----------------
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
    real(r8) :: INFrehom(pcols,pver) !  hom freezing occurence frequency.  1 occur, 0 not occur.
    real(r8) :: INFreIN(pcols,pver)  !  ice nucleation occerence frequency.   1 occur, 0 not occur.

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

    do k = top_lev, pver
       do i = 1, ncol
          rho(i,k) = pmid(i,k)/(rair*t(i,k))
       end do
    end do

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

    if (use_preexisting_ice) then
       fhom(:,:)     = 0.0_r8
       wice(:,:)     = 0.0_r8
       weff(:,:)     = 0.0_r8
       INnso4(:,:)   = 0.0_r8
       INnbc(:,:)    = 0.0_r8
       INndust(:,:)  = 0.0_r8
       INondust(:,:)  = 0.0_r8
       INhet(:,:)    = 0.0_r8
       INhom(:,:)    = 0.0_r8
       INFrehom(:,:) = 0.0_r8
       INFreIN(:,:)  = 0.0_r8
    endif

    do k = top_lev, pver
       ! Get humidity and saturation vapor pressures
       call qsat_water(t(:ncol,k), pmid(:ncol,k), es(:ncol), qs(:ncol), gam=gammas(:ncol))

       do i = 1, ncol
          relhum(i,k) = qn(i,k)/qs(i)
          icldm(i,k) = max(icecldf(i,k), mincld) ! get cloud fraction, check for minimum
       end do
    end do

    do k = top_lev, pver
       do i = 1, ncol

          if (t(i,k) < tmelt - 5._r8) then

             ! compute aerosol number for so4, soot, and dust with units #/cm^3
             so4_num  = 0._r8
             soot_num = 0._r8
             dst1_num = 0._r8
             dst2_num = 0._r8
             dst3_num = 0._r8
             dst4_num = 0._r8
             dst_num  = 0._r8

             if (clim_modal_aero) then
                !For modal aerosols, assume for the upper troposphere:
                ! soot = accumulation mode
                ! sulfate = aiken mode
                ! dust = coarse mode
                ! since modal has internal mixtures.
                soot_num =  numberConcentration(i,k,MODE_IDX_OMBC_INTMIX_COAT_AIT)*1.0e-6_r8

                dst_num = (numberConcentration(i,k,MODE_IDX_DST_A2) &
                     + numberConcentration(i,k,MODE_IDX_DST_A3))*1.0e-6_r8
                !Oslo aerosols have two modes.. Need mode-fractions
                dust_coarse_fraction = numberConcentration(i,k,MODE_IDX_DST_A3)*1.e-6_r8 / (dst_num+1.e-100_r8)


                so4_num = (numberConcentration(i,k,MODE_IDX_SO4_AC))*1.0e-6_r8 

             end if !clim modal aero
             ! *** Turn off soot nucleation ***
             soot_num = 0.0_r8

             call nucleati( &
                  wsubi(i,k), t(i,k), pmid(i,k), relhum(i,k), icldm(i,k),   &
                  qc(i,k), qi(i,k), ni(i,k), rho(i,k),                      &
                  so4_num, dst_num, soot_num, subgrid(i,k),                 &
                  naai(i,k), nihf(i,k), niimm(i,k), nidep(i,k), nimey(i,k), &
                  wice(i,k), weff(i,k), fhom(i,k), regm(i,k),               &
                  oso4_num, odst_num, osoot_num)

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

                   masslost = (odst_num               & !all removed
                        - dst_num*dust_coarse_fraction) & !fraction to coarse mode
                        / volumeToNumber(MODE_IDX_DST_A2) &
                        * rhopart(l_dst_a2) & 
                        /rho(i,k)*1e6_r8 

                   ptend%q(i,k,l_dst_a2) = -masslost*icldm(i,k)/ dtime
                   cld_dst_a2(i,k) = cld_dst_a2(i,k) + masslost*icldm(i,k)

                end if

                ! Coarse mode (is always lost)  
                masslost = (odst_num - numberFromSmallDustMode) &
                     / volumeToNumber(MODE_IDX_DST_A3) &
                     * rhopart(l_dst_a3) & 
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
                if ((k < troplev(i)) .and. (nucleate_ice_strat > 0._r8)) then
                   if (oso4_num > 0._r8) then
                      so4_num_ac = so4_num*rho(i,k)*1.0e-6_r8 !This is maximum sulfate which can activate
                      ! NCAR/MAM4-version
                      ! so4_num_ac = num_accum(i,k)*rho(i,k)*1.0e-6_r8
                      ! NCAR/MAM4-version
                      dso4_num = max(0._r8, (nucleate_ice_strat * (so4_num_cr + so4_num_ac)) - oso4_num) * 1e6_r8 / rho(i,k)
                      naai(i,k) = naai(i,k) + dso4_num
                      nihf(i,k) = nihf(i,k) + dso4_num
                   end if
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

             naai_hom(i,k) = nihf(i,k)

             ! output activated ice (convert from #/kg -> #/m3)
             nihf(i,k)     = nihf(i,k) *rho(i,k)
             niimm(i,k)    = niimm(i,k)*rho(i,k)
             nidep(i,k)    = nidep(i,k)*rho(i,k)
             nimey(i,k)    = nimey(i,k)*rho(i,k)

             if (use_preexisting_ice) then
                INnso4(i,k) =so4_num*1e6_r8  ! (convert from #/cm3 -> #/m3)
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

          end if
       end do
    end do


    call outfld('NIHF',   nihf, pcols, lchnk)
    call outfld('NIIMM', niimm, pcols, lchnk)
    call outfld('NIDEP', nidep, pcols, lchnk)
    call outfld('NIMEY', nimey, pcols, lchnk)
    call outfld('NIREGM', regm, pcols, lchnk)
    call outfld('NISUBGRID', subgrid, pcols, lchnk)
    call outfld('NITROP_PD', trop_pd, pcols, lchnk)

    if (use_preexisting_ice) then
       call outfld( 'fhom' , fhom, pcols, lchnk)
       call outfld( 'WICE' , wice, pcols, lchnk)
       call outfld( 'WEFF' , weff, pcols, lchnk)
       call outfld('INnso4  ',INnso4 , pcols,lchnk)
       call outfld('INnbc   ',INnbc  , pcols,lchnk)
       call outfld('INndust ',INndust, pcols,lchnk)
       call outfld('INondust ',INondust, pcols,lchnk)
       call outfld('INhet   ',INhet  , pcols,lchnk)
       call outfld('INhom   ',INhom  , pcols,lchnk)
       call outfld('INFrehom',INFrehom,pcols,lchnk)
       call outfld('INFreIN ',INFreIN, pcols,lchnk)
    end if

  end subroutine nucleate_ice_oslo_calc

  !===============================================================================

  subroutine nucleati(  &
       wbar, tair, pmid, relhum, cldn,      &
       qc, qi, ni_in, rhoair,               &
       so4_num, dst_num, soot_num, subgrid, &
       nuci, onihf, oniimm, onidep, onimey, &
       wpice, weff, fhom, regm, &
       oso4_num, odst_num, osoot_num)

    ! Input Arguments
    real(r8), intent(in) :: wbar        ! grid cell mean vertical velocity (m/s)
    real(r8), intent(in) :: tair        ! temperature (K)
    real(r8), intent(in) :: pmid        ! pressure at layer midpoints (pa)
    real(r8), intent(in) :: relhum      ! relative humidity with respective to liquid
    real(r8), intent(in) :: cldn        ! new value of cloud fraction    (fraction)
    real(r8), intent(in) :: qc          ! liquid water mixing ratio (kg/kg)
    real(r8), intent(in) :: qi          ! grid-mean preexisting cloud ice mass mixing ratio (kg/kg)
    real(r8), intent(in) :: ni_in       ! grid-mean preexisting cloud ice number conc (#/kg) 
    real(r8), intent(in) :: rhoair      ! air density (kg/m3)
    real(r8), intent(in) :: so4_num     ! so4 aerosol number (#/cm^3)
    real(r8), intent(in) :: dst_num     ! total dust aerosol number (#/cm^3)
    real(r8), intent(in) :: soot_num    ! soot (hydrophilic) aerosol number (#/cm^3)
    real(r8), intent(in) :: subgrid     ! subgrid saturation scaling factor

    ! Output Arguments
    real(r8), intent(out) :: nuci       ! ice number nucleated (#/kg)
    real(r8), intent(out) :: onihf      ! nucleated number from homogeneous freezing of so4
    real(r8), intent(out) :: oniimm     ! nucleated number from immersion freezing
    real(r8), intent(out) :: onidep     ! nucleated number from deposition nucleation
    real(r8), intent(out) :: onimey     ! nucleated number from deposition nucleation  (meyers: mixed phase)
    real(r8), intent(out) :: wpice      ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at Shom
    real(r8), intent(out) :: weff       ! effective Vertical velocity for ice nucleation (m/s); weff=wbar-wpice
    real(r8), intent(out) :: fhom       ! how much fraction of cloud can reach Shom
    real(r8), intent(out) :: regm       ! nucleation regime indiator
    real(r8), intent(out) :: oso4_num   ! so4 aerosol number (#/cm^3)
    real(r8), intent(out) :: odst_num   ! total dust aerosol number (#/cm^3)
    real(r8), intent(out) :: osoot_num  ! soot (hydrophilic) aerosol number (#/cm^3)

    ! Local workspace
    real(r8) :: nihf                      ! nucleated number from homogeneous freezing of so4
    real(r8) :: niimm                     ! nucleated number from immersion freezing
    real(r8) :: nidep                     ! nucleated number from deposition nucleation
    real(r8) :: nimey                     ! nucleated number from deposition nucleation (meyers)
    real(r8) :: n1, ni                    ! nucleated number
    real(r8) :: tc, A, B                  ! work variable
    real(r8) :: esl, esi, deles           ! work variable
    real(r8) :: wbar1, wbar2

    ! used in SUBROUTINE Vpreice
    real(r8) :: Ni_preice        ! cloud ice number conc (1/m3)   
    real(r8) :: lami,Ri_preice   ! mean cloud ice radius (m)
    real(r8) :: Shom             ! initial ice saturation ratio; if <1, use hom threshold Si
    real(r8) :: detaT,RHimean    ! temperature standard deviation, mean cloudy RHi
    real(r8) :: wpicehet   ! diagnosed Vertical velocity Reduction caused by preexisting ice (m/s), at shet
    real(r8) :: weffhet    ! effective Vertical velocity for ice nucleation (m/s)  weff=wbar-wpicehet 
    !-------------------------------------------------------------------------------

    RHimean = relhum*svp_water(tair)/svp_ice(tair)*subgrid

    ! temp variables that depend on use_preexisting_ice
    wbar1 = wbar
    wbar2 = wbar

    ! If not using prexisting ice, the homogeneous freezing happens in the
    ! entire gridbox.
    fhom = 1._r8

    if (use_preexisting_ice) then

       Ni_preice = ni_in*rhoair                    ! (convert from #/kg -> #/m3)
       Ni_preice = Ni_preice / max(mincld,cldn)   ! in-cloud ice number density 

       if (Ni_preice > 10.0_r8 .and. qi > 1.e-10_r8) then    ! > 0.01/L = 10/m3   
          Shom = -1.5_r8   ! if Shom<1 , Shom will be recalculated in SUBROUTINE Vpreice, according to Ren & McKenzie, 2005
          lami = (gamma4*ci*ni_in/qi)**(1._r8/3._r8)
          Ri_preice = 0.5_r8/lami                  ! radius
          Ri_preice = max(Ri_preice, 1e-8_r8)       ! >0.01micron
          call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shom, wpice)
          call Vpreice(pmid, tair, Ri_preice, Ni_preice, Shet, wpicehet)
       else
          wpice    = 0.0_r8
          wpicehet = 0.0_r8
       endif

       weff     = max(wbar-wpice, minweff)
       wpice    = min(wpice, wbar)
       weffhet  = max(wbar-wpicehet,minweff)
       wpicehet = min(wpicehet, wbar)

       wbar1 = weff
       wbar2 = weffhet

       detaT   = wbar/0.23_r8
       if (use_incloud_nuc) then
          call frachom(tair, 1._r8, detaT, fhom)
       else
          call frachom(tair, RHimean, detaT, fhom)
       end if
    end if

    ni = 0._r8
    tc = tair - 273.15_r8

    ! initialize
    niimm = 0._r8
    nidep = 0._r8
    nihf  = 0._r8
    deles = 0._r8
    esi   = 0._r8
    regm  = 0._r8

    oso4_num  = 0._r8
    odst_num  = 0._r8
    osoot_num = 0._r8

    if ((so4_num >= 1.0e-10_r8 .or. (soot_num+dst_num) >= 1.0e-10_r8) .and. cldn > 0._r8) then

       if (RHimean.ge.1.2_r8) then

          if ( ((tc.le.0.0_r8).and.(tc.ge.-37.0_r8).and.(qc.lt.1.e-12_r8)).or.(tc.le.-37.0_r8)) then

             A = -1.4938_r8 * log(soot_num+dst_num) + 12.884_r8
             B = -10.41_r8  * log(soot_num+dst_num) - 67.69_r8
             regm = A * log(wbar1) + B

             ! heterogeneous nucleation only
             if (tc .gt. regm .or. so4_num < 1.0e-10_r8) then

                if(tc.lt.-40._r8 .and. wbar1.gt.1._r8 .and. so4_num >= 1.0e-10_r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                   call hf(tc,wbar1,relhum*subgrid,so4_num,nihf)
                   niimm=0._r8
                   nidep=0._r8

                   ! If some homogeneous nucleation happened, assume all of the that heterogeneous
                   ! and coarse mode sulfate particles nucleated.
                   if (nihf.gt.1e-3_r8) then ! hom occur,  add preexisting ice
                      niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                      odst_num  = dst_num
                      osoot_num = soot_num

                      oso4_num  = nihf
                   endif

                   nihf      = nihf * fhom
                   oso4_num  = oso4_num * fhom

                   n1        = nihf + niimm
                else

                   call hetero(tc,wbar2,soot_num+dst_num,niimm,nidep)

                   nihf = 0._r8
                   n1   = niimm + nidep

                   osoot_num = soot_num * (niimm + nidep) / (soot_num + dst_num) 
                   odst_num  = dst_num  * (niimm + nidep) / (soot_num + dst_num) 
                endif

                ! homogeneous nucleation only
             else if (tc.lt.regm-5._r8 .or. (soot_num+dst_num) < 1.0e-10_r8) then

                call hf(tc,wbar1,relhum*subgrid,so4_num,nihf)
                niimm=0._r8
                nidep=0._r8

                ! If some homogeneous nucleation happened, assume all of the that
                ! heterogeneous and coarse mode sulfate particles nucleated.
                if (nihf.gt.1e-3_r8) then !  hom occur,  add preexisting ice
                   niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                   odst_num  = dst_num
                   osoot_num = soot_num

                   oso4_num  = nihf
                endif

                nihf      = nihf * fhom
                oso4_num  = oso4_num * fhom

                n1        = nihf + niimm

                ! transition between homogeneous and heterogeneous: interpolate in-between
             else

                if (tc.lt.-40._r8 .and. wbar1.gt.1._r8) then ! exclude T<-40 & W>1m/s from hetero. nucleation

                   call hf(tc, wbar1, relhum*subgrid, so4_num, nihf)
                   niimm = 0._r8
                   nidep = 0._r8

                   ! If some homogeneous nucleation happened, assume all of the
                   ! that heterogeneous and coarse mode sulfate particles nucleated.
                   if (nihf.gt.1e-3_r8) then ! hom occur,  add preexisting ice
                      niimm     = dst_num + soot_num       ! assuming dst_num freeze firstly
                      odst_num  = dst_num
                      osoot_num = soot_num

                      oso4_num  = nihf
                   endif

                   nihf      = nihf * fhom
                   oso4_num  = oso4_num * fhom

                   n1        = nihf + niimm

                else

                   call hf(regm-5._r8,wbar1,relhum*subgrid,so4_num,nihf)
                   call hetero(regm,wbar2,soot_num+dst_num,niimm,nidep)

                   ! If some homogeneous nucleation happened, assume all of the
                   ! heterogeneous particles nucleated and add in a fraction of
                   ! the homogeneous freezing.
                   if (nihf.gt.1e-3_r8) then ! hom occur,  add preexisting ice
                      oso4_num  = nihf
                   endif

                   osoot_num = soot_num * (niimm + nidep) / (soot_num + dst_num)
                   odst_num  = dst_num  * (niimm + nidep) / (soot_num + dst_num)

                   nihf      = nihf      * fhom * ((regm - tc) / 5._r8)**2
                   oso4_num  = oso4_num  * fhom * ((regm - tc) / 5._r8)**2

                   n1 = niimm + nidep + nihf

                end if
             end if

             ! Scale the rates for in-cloud number, since this is what
             ! MG is expecting to find.
             ni = n1

             ! If using prexsiting ice, then add it to the total.
             if (use_preexisting_ice) then
                ni = ni + Ni_preice * 1e-6_r8
             end if
          end if
       end if
    end if

    ! deposition/condensation nucleation in mixed clouds (-37<T<0C) (Meyers, 1992)
    if(tc.lt.0._r8 .and. tc.gt.-37._r8 .and. qc.gt.1.e-12_r8) then
       esl = svp_water(tair)     ! over water in mixed clouds
       esi = svp_ice(tair)     ! over ice
       deles = (esl - esi)
       nimey=1.e-3_r8*exp(12.96_r8*deles/esi - 0.639_r8) 
    else
       nimey=0._r8
    endif

    if (use_hetfrz_classnuc) nimey = 0._r8

    nuci=ni + nimey

    if(nuci.gt.9999._r8.or.nuci.lt.0._r8) then
       write(iulog, *) 'Warning: incorrect ice nucleation number (nuci reset =0)'
       write(iulog, *) ni, tair, relhum, wbar, nihf, niimm, nidep,deles,esi,dst_num,so4_num
       nuci=0._r8
    endif

    nuci   = nuci*1.e+6_r8/rhoair    ! change unit from #/cm3 to #/kg
    onimey = nimey*1.e+6_r8/rhoair
    onidep = nidep*1.e+6_r8/rhoair
    oniimm = niimm*1.e+6_r8/rhoair
    onihf  = nihf*1.e+6_r8/rhoair
  end subroutine nucleati

  !===============================================================================

  subroutine hetero(T,ww,Ns,Nis,Nid)

    real(r8), intent(in)  :: T, ww, Ns
    real(r8), intent(out) :: Nis, Nid

    real(r8) A11,A12,A21,A22,B11,B12,B21,B22
    real(r8) B,C

    ! parameters
    A11 = 0.0263_r8
    A12 = -0.0185_r8
    A21 = 2.758_r8
    A22 = 1.3221_r8
    B11 = -0.008_r8
    B12 = -0.0468_r8
    B21 = -0.2667_r8
    B22 = -1.4588_r8

    !     ice from immersion nucleation (cm^-3)

    B = (A11+B11*log(Ns)) * log(ww) + (A12+B12*log(Ns))
    C =  A21+B21*log(Ns)

    Nis = exp(A22) * Ns**B22 * exp(B*T) * ww**C
    Nis = min(Nis,Ns)

    Nid = 0.0_r8    ! don't include deposition nucleation for cirrus clouds when T<-37C

  end subroutine hetero

  !===============================================================================

  subroutine hf(T,ww,RH,Na,Ni)

    real(r8), intent(in)  :: T, ww, RH, Na
    real(r8), intent(out) :: Ni

    real(r8)    A1_fast,A21_fast,A22_fast,B1_fast,B21_fast,B22_fast
    real(r8)    A2_fast,B2_fast
    real(r8)    C1_fast,C2_fast,k1_fast,k2_fast
    real(r8)    A1_slow,A2_slow,B1_slow,B2_slow,B3_slow
    real(r8)    C1_slow,C2_slow,k1_slow,k2_slow
    real(r8)    regm
    real(r8)    A,B,C
    real(r8)    RHw

    ! parameters
    A1_fast  =0.0231_r8
    A21_fast =-1.6387_r8  !(T>-64 deg)
    A22_fast =-6.045_r8   !(T<=-64 deg)
    B1_fast  =-0.008_r8
    B21_fast =-0.042_r8   !(T>-64 deg)
    B22_fast =-0.112_r8   !(T<=-64 deg)
    C1_fast  =0.0739_r8
    C2_fast  =1.2372_r8

    A1_slow  =-0.3949_r8
    A2_slow  =1.282_r8
    B1_slow  =-0.0156_r8
    B2_slow  =0.0111_r8
    B3_slow  =0.0217_r8
    C1_slow  =0.120_r8
    C2_slow  =2.312_r8

    Ni = 0.0_r8

    !RHw parameters
    A = 6.0e-4_r8*log(ww)+6.6e-3_r8
    B = 6.0e-2_r8*log(ww)+1.052_r8
    C = 1.68_r8  *log(ww)+129.35_r8
    RHw=(A*T*T+B*T+C)*0.01_r8

    if((T.le.-37.0_r8) .and. ((RH).ge.RHw)) then

       regm = 6.07_r8*log(ww)-55.0_r8

       if(T.ge.regm) then    ! fast-growth regime

          if(T.gt.-64.0_r8) then
             A2_fast=A21_fast
             B2_fast=B21_fast
          else
             A2_fast=A22_fast
             B2_fast=B22_fast
          endif

          k1_fast = exp(A2_fast + B2_fast*T + C2_fast*log(ww))
          k2_fast = A1_fast+B1_fast*T+C1_fast*log(ww)

          Ni = k1_fast*Na**(k2_fast)
          Ni = min(Ni,Na)

       else       ! slow-growth regime

          k1_slow = exp(A2_slow + (B2_slow+B3_slow*log(ww))*T + C2_slow*log(ww))
          k2_slow = A1_slow+B1_slow*T+C1_slow*log(ww)

          Ni = k1_slow*Na**(k2_slow)
          Ni = min(Ni,Na)

       endif

    end if

  end subroutine hf

  !===============================================================================

  subroutine Vpreice(P_in, T_in, R_in, C_in, S_in, V_out)

    !  based on  Karcher et al. (2006)
    !  VERTICAL VELOCITY CALCULATED FROM DEPOSITIONAL LOSS TERM

    ! arguments
    REAL(r8), INTENT(in)  :: P_in       ! [Pa],INITIAL AIR pressure 
    REAL(r8), INTENT(in)  :: T_in       ! [K] ,INITIAL AIR temperature 
    REAL(r8), INTENT(in)  :: R_in       ! [m],INITIAL MEAN  ICE CRYSTAL NUMBER RADIUS 
    REAL(r8), INTENT(in)  :: C_in       ! [m-3],INITIAL TOTAL ICE CRYSTAL NUMBER DENSITY, [1/cm3]
    REAL(r8), INTENT(in)  :: S_in       ! [-],INITIAL ICE SATURATION RATIO;; if <1, use hom threshold Si 
    REAL(r8), INTENT(out) :: V_out      ! [m/s], VERTICAL VELOCITY REDUCTION (caused by preexisting ice)

    ! parameters
    REAL(r8), PARAMETER :: ALPHAc  = 0.5_r8 ! density of ice (g/cm3), !!!V is not related to ALPHAc 
    REAL(r8), PARAMETER :: FA1c    = 0.601272523_r8        
    REAL(r8), PARAMETER :: FA2c    = 0.000342181855_r8
    REAL(r8), PARAMETER :: FA3c    = 1.49236645E-12_r8        
    REAL(r8), PARAMETER :: WVP1c   = 3.6E+10_r8   
    REAL(r8), PARAMETER :: WVP2c   = 6145.0_r8
    REAL(r8), PARAMETER :: FVTHc   = 11713803.0_r8
    REAL(r8), PARAMETER :: THOUBKc = 7.24637701E+18_r8
    REAL(r8), PARAMETER :: SVOLc   = 3.23E-23_r8    ! SVOL=XMW/RHOICE
    REAL(r8), PARAMETER :: FDc     = 249.239822_r8
    REAL(r8), PARAMETER :: FPIVOLc = 3.89051704E+23_r8         
    REAL(r8) :: T,P,S,R,C
    REAL(r8) :: A1,A2,A3,B1,B2
    REAL(r8) :: T_1,PICE,FLUX,ALP4,CISAT,DLOSS,VICE

    T = T_in          ! K  , K
    P = P_in*1e-2_r8  ! Pa , hpa

    IF (S_in.LT.1.0_r8) THEN
       S = 2.349_r8 - (T/259.0_r8) ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
    ELSE
       S = S_in                    ! INPUT ICE SATURATION RATIO, -,  >1
    ENDIF

    R     = R_in*1e2_r8   ! m  => cm
    C     = C_in*1e-6_r8  ! m-3 => cm-3
    T_1   = 1.0_r8/ T
    PICE  = WVP1c * EXP(-(WVP2c*T_1))
    ALP4  = 0.25_r8 * ALPHAc      
    FLUX  = ALP4 * SQRT(FVTHc*T)
    CISAT = THOUBKc * PICE * T_1   
    A1    = ( FA1c * T_1 - FA2c ) * T_1 
    A2    = 1.0_r8/ CISAT      
    A3    = FA3c * T_1 / P
    B1    = FLUX * SVOLc * CISAT * ( S-1.0_r8 ) 
    B2    = FLUX * FDc * P * T_1**1.94_r8 
    DLOSS = FPIVOLc * C * B1 * R**2 / ( 1.0_r8+ B2 * R )         
    VICE  = ( A2 + A3 * S ) * DLOSS / ( A1 * S )  ! 2006,(19)
    V_out = VICE*1e-2_r8  ! cm/s => m/s

  end subroutine Vpreice

  !===============================================================================

  subroutine frachom(Tmean,RHimean,detaT,fhom)

    ! How much fraction of cirrus might reach Shom  
    ! base on "A cirrus cloud scheme for general circulation models",
    ! B. Karcher and U. Burkhardt 2008

    real(r8), intent(in)  :: Tmean, RHimean, detaT
    real(r8), intent(out) :: fhom

    real(r8), parameter :: seta = 6132.9_r8  ! K
    integer,  parameter :: Nbin=200          ! (Tmean - 3*detaT, Tmean + 3*detaT)

    real(r8) :: PDF_T(Nbin)    ! temperature PDF;  ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
    real(r8) :: Sbin(Nbin)     ! the fluctuations of Si that are driven by the T variations 
    real(r8) :: Sihom, deta
    integer  :: i

    Sihom = 2.349_r8-Tmean/259.0_r8   ! homogeneous freezing threshold, according to Ren & McKenzie, 2005
    fhom  = 0.0_r8

    do i = Nbin, 1, -1

       deta     = (i - 0.5_r8 - Nbin/2)*6.0_r8/Nbin   ! PDF_T=0  outside (Tmean-3*detaT, Tmean+3*detaT)
       Sbin(i)  = RHimean*exp(deta*detaT*seta/Tmean**2.0_r8)
       PDF_T(i) = exp(-deta**2.0_r8/2.0_r8)*6.0_r8/(sqrt(2.0_r8*Pi)*Nbin)


       if (Sbin(i).ge.Sihom) then
          fhom = fhom + PDF_T(i)
       else
          exit
       end if
    end do

    fhom = min(1.0_r8, fhom/0.997_r8)   ! accounting for the finite limits (-3 , 3)
  end subroutine frachom

end module oslo_aero_nucleate_ice
