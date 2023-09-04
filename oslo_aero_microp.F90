module oslo_aero_microp

  !---------------------------------------------------------------------------------
  ! Oslo-aero driver layer for aerosol activation processes.
  ! Author: Andrew Gettelman
  ! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan! May 2010
  ! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
  !    Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
  ! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
  !    microphysics and macrophysics and concentrate aerosol information here
  !---------------------------------------------------------------------------------

  use shr_kind_mod,           only: r8=>shr_kind_r8
  use spmd_utils,             only: masterproc
  use ppgrid,                 only: pcols, pver, pverp
  use ref_pres,               only: top_lev => trop_cloud_top_lev
  use physconst,              only: rair
  use constituents,           only: cnst_get_ind, pcnst
  use physics_types,          only: physics_state, physics_ptend, physics_ptend_init, physics_ptend_sum
  use physics_types,          only: physics_state_copy, physics_update
  use physics_buffer,         only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
  use phys_control,           only: phys_getopts, use_hetfrz_classnuc
  use rad_constituents,       only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, rad_cnst_get_mode_num
  use ndrop_bam,              only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn
  use cam_history,            only: addfld, add_default, outfld
  use cam_logfile,            only: iulog
  !
  use oslo_aero_ndrop,        only: ndrop_init_oslo, dropmixnuc_oslo
  use oslo_aero_conc,         only: oslo_aero_conc_calc
  use oslo_aero_hetfrz,       only: hetfrz_classnuc_oslo_register, hetfrz_classnuc_oslo_init, hetfrz_classnuc_oslo_readnl
  use oslo_aero_hetfrz,       only: hetfrz_classnuc_oslo_calc, hetfrz_classnuc_oslo_save_cbaero
  use oslo_aero_nucleate_ice, only: nucleate_ice_oslo_register, nucleate_ice_oslo_init, nucleate_ice_oslo_readnl
  use oslo_aero_nucleate_ice, only: nucleate_ice_oslo_calc, use_preexisting_ice
  use oslo_aero_params,       only: nmodes_oslo => nmodes
  use oslo_aero_share,        only: MODE_IDX_DST_A2, MODE_IDX_DST_A3, MODE_IDX_SO4_AC, MODE_IDX_OMBC_INTMIX_COAT_AIT
  use oslo_aero_share,        only: lifeCycleNumberMedianRadius, l_dst_a2, l_dst_a3, l_bc_ai
  use oslo_aero_share,        only: getNumberOfTracersInMode, getTracerIndex, getCloudTracerIndex

  implicit none
  private

  public :: oslo_aero_microp_init, oslo_aero_microp_run, oslo_aero_microp_readnl, oslo_aero_microp_register

  ! Private module data

  character(len=16)   :: eddy_scheme

  ! contact freezing due to dust, dust number mean radius (m), 
  ! Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
  real(r8), parameter :: rn_dst1 = 0.258e-6_r8
  real(r8), parameter :: rn_dst2 = 0.717e-6_r8
  real(r8), parameter :: rn_dst3 = 1.576e-6_r8
  real(r8), parameter :: rn_dst4 = 3.026e-6_r8

  ! smallest mixing ratio considered in microphysics
  real(r8), parameter :: qsmall = 1.e-18_r8

  ! minimum allowed cloud fraction
  real(r8), parameter :: mincld = 0.0001_r8

  ! indices in state%q and pbuf structures
  integer :: cldliq_idx   = -1
  integer :: cldice_idx   = -1
  integer :: numliq_idx   = -1
  integer :: numice_idx   = -1
  integer :: kvh_idx      = -1
  integer :: tke_idx      = -1
  integer :: wp2_idx      = -1
  integer :: ast_idx      = -1
  integer :: cldo_idx     = -1
  integer :: dgnumwet_idx = -1

  ! prescribed aerosol bulk sulfur scale factor
  real(r8) :: bulk_scale    

  integer :: npccn_idx, rndst_idx, nacon_idx

!=========================================================================================
contains
!=========================================================================================

  subroutine oslo_aero_microp_readnl(nlfile)

    use namelist_utils, only: find_group_name
    use cam_abortutils, only: endrun
    use mpishorthand

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Namelist variables
    real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'microp_aero_readnl'

    namelist /microp_aero_nl/ microp_aero_bulk_scale
    !-----------------------------------------------------------------------------

    if (masterproc) then
       open(newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'microp_aero_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, microp_aero_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if
#ifdef SPMD
    call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
#endif

    ! set local variables
    bulk_scale = microp_aero_bulk_scale

    call nucleate_ice_oslo_readnl(nlfile)
    call hetfrz_classnuc_oslo_readnl(nlfile)

  end subroutine oslo_aero_microp_readnl

  !=========================================================================================
  subroutine oslo_aero_microp_register
    !----------------------------------------------------------------------- 
    ! Register pbuf fields for aerosols needed by microphysics
    ! Author: Cheryl Craig October 2012
    !-----------------------------------------------------------------------

    use physics_buffer, only: pbuf_add_field, dtype_r8

    call pbuf_add_field('NPCCN', 'physpkg',dtype_r8,(/pcols,pver/)  , npccn_idx)
    call pbuf_add_field('RNDST', 'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
    call pbuf_add_field('NACON', 'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)

    call nucleate_ice_oslo_register()
    call hetfrz_classnuc_oslo_register()

  end subroutine oslo_aero_microp_register

  !=========================================================================================

  subroutine oslo_aero_microp_init

    !----------------------------------------------------------------------- 
    ! Initialize constants for aerosols needed by microphysics
    ! Author: Andrew Gettelman May 2010
    !-----------------------------------------------------------------------

    ! local variables
    integer  :: iaer, ierr
    integer  :: m, n, nmodes, nspec

    character(len=32) :: str32
    character(len=*), parameter :: routine = 'oslo_aero_microp_init'
    logical :: history_amwg
    !-----------------------------------------------------------------------

    ! Query the PBL eddy scheme
    call phys_getopts(eddy_scheme_out=eddy_scheme, history_amwg_out=history_amwg )

    ! Access the physical properties of the aerosols that are affecting the climate
    ! by using routines from the rad_constituents module.

    ! get indices into state and pbuf structures
    call cnst_get_ind('CLDLIQ', cldliq_idx)
    call cnst_get_ind('CLDICE', cldice_idx)
    call cnst_get_ind('NUMLIQ', numliq_idx)
    call cnst_get_ind('NUMICE', numice_idx)

    select case(trim(eddy_scheme))
    case ('diag_TKE')
       tke_idx = pbuf_get_index('tke')   
    case ('CLUBB_SGS')
       wp2_idx = pbuf_get_index('WP2_nadv')
    case default
       kvh_idx = pbuf_get_index('kvh')
    end select
    ast_idx = pbuf_get_index('AST')
    cldo_idx = pbuf_get_index('CLDO')

    call addfld('LCLOUD', (/ 'lev' /), 'A', ' ',   'Liquid cloud fraction used in stratus activation')
    call addfld('WSUB',   (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity'                   )
    call addfld('WSUBI',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice'           )
    if (history_amwg) then
       call add_default ('WSUB     ', 1, ' ')
    end if

    call ndrop_init_oslo()
    call nucleate_ice_oslo_init(mincld, bulk_scale)
    call hetfrz_classnuc_oslo_init(mincld)

  end subroutine oslo_aero_microp_init

  !=========================================================================================
  subroutine oslo_aero_microp_run (state, ptend_all, deltatin, pbuf)

    ! arguments
    type(physics_state),         intent(in)    :: state
    type(physics_ptend),         intent(out)   :: ptend_all
    real(r8),                    intent(in)    :: deltatin     ! time step (s)
    type(physics_buffer_desc),   pointer       :: pbuf(:)

    ! local workspace
    ! all units mks unless otherwise stated
    integer :: i, k, m
    integer :: itim_old
    integer :: nmodes
    type(physics_state) :: state1                             ! Local copy of state variable
    type(physics_ptend) :: ptend_loc
    real(r8), pointer :: ast(:,:)        
    real(r8), pointer :: npccn(:,:)                           ! number of CCN (liquid activated)
    real(r8), pointer :: rndst(:,:,:)                         ! radius of 4 dust bins for contact freezing
    real(r8), pointer :: nacon(:,:,:)                         ! number in 4 dust bins for contact freezing
    real(r8), pointer :: num_coarse(:,:)                      ! number m.r. of coarse mode
    real(r8), pointer :: coarse_dust(:,:)                     ! mass m.r. of coarse dust
    real(r8), pointer :: coarse_nacl(:,:)                     ! mass m.r. of coarse nacl
    real(r8), pointer :: coarse_so4(:,:)                      ! mass m.r. of coarse sulfate
    real(r8), pointer :: kvh(:,:)                             ! vertical eddy diff coef (m2 s-1)
    real(r8), pointer :: tke(:,:)                             ! TKE from the UW PBL scheme (m2 s-2)
    real(r8), pointer :: wp2(:,:)                             ! CLUBB vertical velocity variance
    real(r8), pointer :: cldn(:,:)                            ! cloud fraction
    real(r8), pointer :: cldo(:,:)                            ! old cloud fraction
    real(r8), pointer :: dgnumwet(:,:,:)                      ! aerosol mode diameter
    real(r8), pointer :: aer_mmr(:,:)                         ! aerosol mass mixing ratio
    real(r8) :: rho(pcols,pver)                               ! air density (kg m-3)
    real(r8) :: lcldm(pcols,pver)                             ! liq cloud fraction
    real(r8) :: lcldn(pcols,pver)                             ! fractional coverage of new liquid cloud
    real(r8) :: lcldo(pcols,pver)                             ! fractional coverage of old liquid cloud
    real(r8) :: cldliqf(pcols,pver)                           ! fractional of total cloud that is liquid
    real(r8) :: qcld                                          ! total cloud water
    real(r8) :: nctend_mixnuc(pcols,pver)
    real(r8) :: dum, dum2                                     ! temporary dummy variable
    real(r8) :: dmc, ssmc, so4mc                              ! variables for modal scheme.
    integer  :: dst_idx, num_idx
    real(r8) :: wsub(pcols,pver)                              ! diagnosed sub-grid vertical velocity st. dev. (m/s)
    real(r8) :: wsubi(pcols,pver)                             ! diagnosed sub-grid vertical velocity ice (m/s)
    real(r8) :: nucboas
    real(r8) :: wght
    integer  :: lchnk, ncol
    real(r8) :: factnum(pcols,pver,0:nmodes_oslo)             ! activation fraction for aerosol number
    real(r8) :: qaercwpt(pcols,pver,pcnst)
    logical  :: hasAerosol(pcols, pver, nmodes_oslo)
    real(r8) :: f_acm(pcols,pver, nmodes_oslo)
    real(r8) :: f_bcm(pcols,pver, nmodes_oslo)
    real(r8) :: f_aqm(pcols, pver, nmodes_oslo)
    real(r8) :: f_so4_condm(pcols, pver, nmodes_oslo)         !Needed in "get component fraction"
    real(r8) :: f_soam(pcols, pver, nmodes_oslo)              !Needed in "get component fraction"
    real(r8) :: numberConcentration(pcols,pver,0:nmodes_oslo) ![#/m3] number concentraiton
    real(r8) :: volumeConcentration(pcols,pver,nmodes_oslo)   ![m3/m3] volume concentration
    real(r8) :: hygroscopicity(pcols,pver,nmodes_oslo)        ![mol_{aer}/mol_{water}] hygroscopicity
    real(r8) :: lnsigma(pcols,pver,nmodes_oslo)               ![-] log(base e) sigma
    real(r8) :: CProcessModes(pcols,pver)
    real(r8) :: cam(pcols,pver,nmodes_oslo)
    real(r8) :: f_c(pcols, pver)
    real(r8) :: f_aq(pcols,pver)
    real(r8) :: f_bc(pcols,pver)
    real(r8) :: f_so4_cond(pcols,pver)
    real(r8) :: f_soa(pcols,pver)
    real(r8) :: volumeCore(pcols,pver,nmodes_oslo)
    real(r8) :: volumeCoat(pcols,pver,nmodes_oslo)
    !-------------------------------------------------------------------------------

    call physics_state_copy(state,state1)

    lchnk = state1%lchnk
    ncol  = state1%ncol

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))
    call pbuf_get_field(pbuf, npccn_idx, npccn)
    call pbuf_get_field(pbuf, nacon_idx, nacon)
    call pbuf_get_field(pbuf, rndst_idx, rndst)

    call physics_ptend_init(ptend_all, state%psetcols, 'microp_aero')

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
    call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ! initialize output
    npccn(1:ncol,1:pver)    = 0._r8  
    nacon(1:ncol,1:pver,:)  = 0._r8

    ! set default or fixed dust bins for contact freezing
    rndst(1:ncol,1:pver,1) = rn_dst1
    rndst(1:ncol,1:pver,2) = rn_dst2
    rndst(1:ncol,1:pver,3) = rn_dst3
    rndst(1:ncol,1:pver,4) = rn_dst4

    ! save copy of cloud borne aerosols for use in heterogeneous freezing
    if (use_hetfrz_classnuc) then
       call hetfrz_classnuc_oslo_save_cbaero(state, pbuf)
    end if

    ! initialize time-varying parameters
    do k = top_lev, pver
       do i = 1, ncol
          rho(i,k) = state1%pmid(i,k)/(rair*state1%t(i,k))
       end do
    end do

    factnum(1:ncol,1:pver,0:nmodes_oslo) = 0._r8
    cam(:,:,:) = 0._r8

    ! More refined computation of sub-grid vertical velocity 
    ! Set to be zero at the surface by initialization.

    select case (trim(eddy_scheme))
    case ('diag_TKE')
       call pbuf_get_field(pbuf, tke_idx, tke)
    case ('CLUBB_SGS')
       itim_old = pbuf_old_tim_idx()
       call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))
       allocate(tke(pcols,pverp))
       tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

    case default
       call pbuf_get_field(pbuf, kvh_idx, kvh)
    end select

    ! Set minimum values above top_lev.
    wsub(:ncol,:top_lev-1)  = 0.20_r8
    wsubi(:ncol,:top_lev-1) = 0.001_r8

    do k = top_lev, pver
       do i = 1, ncol

          select case (trim(eddy_scheme))
          case ('diag_TKE', 'CLUBB_SGS')
             wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
             wsub(i,k) = min(wsub(i,k),10._r8)
          case default 
             ! get sub-grid vertical velocity from diff coef.
             ! following morrison et al. 2005, JAS
             ! assume mixing length of 30 m
             dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
             ! use maximum sub-grid vertical vel of 10 m/s
             dum = min(dum, 10._r8)
             ! set wsub to value at current vertical level
             wsub(i,k)  = dum
          end select

          wsubi(i,k) = max(0.001_r8, wsub(i,k))
          if (.not. use_preexisting_ice) then
             wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
          endif
          wsub(i,k)  = max(0.20_r8, wsub(i,k))

       end do
    end do

    call outfld('WSUB',   wsub, pcols, lchnk)
    call outfld('WSUBI', wsubi, pcols, lchnk)

    if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

    ! Get size distributed interstitial aerosol
    call oslo_aero_conc_calc(ncol, state%q, rho, CProcessModes, &
         f_c, f_bc, f_aq, f_so4_cond, f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam, &
         numberConcentration, volumeConcentration, hygroscopicity, lnsigma, hasAerosol, volumeCore, volumeCoat)

    ! -----------------
    ! ICE Nucleation
    ! -----------------
    call nucleate_ice_oslo_calc(state1, wsubi, pbuf, deltatin, ptend_loc, numberConcentration)

    call physics_ptend_sum(ptend_loc, ptend_all, ncol)
    call physics_update(state1, ptend_loc, deltatin)

    ! get liquid cloud fraction, check for minimum
    do k = top_lev, pver
       do i = 1, ncol
          lcldm(i,k) = max(ast(i,k), mincld)
       end do
    end do

    ! -----------------
    ! Droplet Activation
    ! -----------------

    ! partition cloud fraction into liquid water part
    lcldn = 0._r8
    lcldo = 0._r8
    cldliqf = 0._r8
    do k = top_lev, pver
       do i = 1, ncol
          qcld = state1%q(i,k,cldliq_idx) + state1%q(i,k,cldice_idx)
          if (qcld > qsmall) then
             lcldn(i,k)   = cldn(i,k)*state1%q(i,k,cldliq_idx)/qcld
             lcldo(i,k)   = cldo(i,k)*state1%q(i,k,cldliq_idx)/qcld
             cldliqf(i,k) = state1%q(i,k,cldliq_idx)/qcld
          end if
       end do
    end do

    call outfld('LCLOUD', lcldn, pcols, lchnk)

    ! If not using preexsiting ice, then only use cloudbourne aerosol for the
    ! liquid clouds. This is the same behavior as CAM5.
    if (use_preexisting_ice) then
       call dropmixnuc_oslo(                            &
            state1, ptend_loc, deltatin, pbuf, wsub,    &  ! Input
            cldn, cldo, cldliqf,                        &
            hasAerosol,                                 &
            CProcessModes, f_c, f_bc, f_aq, f_so4_cond, &
            f_soa,                                      &
            cam, f_acm, f_bcm, f_aqm, f_so4_condm,      &
            f_soam,                                     &
            numberConcentration, volumeConcentration,   &
            hygroscopicity, lnsigma,                    &
            nctend_mixnuc,                              &  ! Output
            factnum )
    else
       ! Note difference in arguments lcldn, lcldo
       cldliqf = 1._r8
       call dropmixnuc_oslo(                            &
            state1, ptend_loc, deltatin, pbuf, wsub,    &  ! Input
            lcldn, lcldo, cldliqf,                      &
            hasAerosol,                                 &
            CProcessModes, f_c, f_bc, f_aq, f_so4_cond, &
            f_soa,                                      &
            cam, f_acm, f_bcm, f_aqm, f_so4_condm,      &
            f_soam,                                     &
            numberConcentration, volumeConcentration,   &
            hygroscopicity, lnsigma,                    &
            nctend_mixnuc,                              &  ! Output
            factnum )
    end if
    npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

    call physics_ptend_sum(ptend_loc, ptend_all, ncol)
    call physics_update(state1, ptend_loc, deltatin)

    ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
    ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics
    do k = top_lev, pver
       do i = 1, ncol
          if (state1%t(i,k) < 269.15_r8) then
             !fxm: I think model uses bins, not modes.. But to get it 
             !approximately correct, use mode radius in first version
             nacon(i,k,2) = numberConcentration(i,k,MODE_IDX_DST_A2)
             nacon(i,k,3) = numberConcentration(i,k,MODE_IDX_DST_A3) 
             rndst(i,k,2) = lifeCycleNumberMedianRadius(MODE_IDX_DST_A2)
             rndst(i,k,3) = lifeCycleNumberMedianRadius(MODE_IDX_DST_A3)
             nacon(i,k,1) = 0.0_r8 !Set to zero to make sure
             nacon(i,k,4) = 0.0_r8 !Set to zero to make sure
          end if
       end do
    end do

    ! heterogeneous freezing
    if (use_hetfrz_classnuc) then
       call hetfrz_classnuc_oslo_calc(state1, deltatin, factnum, pbuf, &
            numberConcentration, volumeConcentration, &
            f_acm, f_bcm, f_aqm, f_so4_condm, f_soam, &
            hygroscopicity, lnsigma, cam, volumeCore, volumeCoat)
    end if

  end subroutine oslo_aero_microp_run

end module oslo_aero_microp
