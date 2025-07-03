module oslo_aero_hetfrz

  !-----------------------------------------------------------------------
  ! Calculate heterogeneous freezing rates from classical nucleation theory
  !
  ! Author:
  !   Corinna Hoose, UiO, May 2009
  !   Yong Wang and Xiaohong Liu, UWyo, 12/2012,
  !   implement in CAM5 and constrain uncertain parameters using natural dust and
  !   BC(soot) datasets.
  !   Yong Wang and Xiaohong Liu, UWyo, 05/2013, implement the PDF-contact angle
  !   approach: Y. Wang et al., Atmos. Chem. Phys., 2014.
  !   Jack Chen, NCAR, 09/2015, modify calculation of dust activation fraction.
  !---------------------------------------------------------------------------------

  use shr_kind_mod,      only: r8=>shr_kind_r8
  use shr_spfn_mod,      only: erf => shr_spfn_erf
  use spmd_utils,        only: mpicom, mstrid=>masterprocid, masterproc
  use spmd_utils,        only: mpi_logical, mpi_real8, mpi_character, mpi_success
  use ppgrid,            only: pcols, pver, begchunk, endchunk
  use physconst,         only: rair, cpair, rh2o, rhoh2o, mwh2o, tmelt, pi
  use constituents,      only: cnst_get_ind, pcnst
  use physics_types,     only: physics_state
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
  use physics_buffer,    only: pbuf_add_field, dtype_r8
  use phys_control,      only: phys_getopts, use_hetfrz_classnuc
  use cam_history,       only: addfld, add_default, outfld
  use ref_pres,          only: top_lev => trop_cloud_top_lev
  use wv_saturation,     only: svp_water, svp_ice
  use cam_logfile,       only: iulog
  use error_messages,    only: handle_errmsg, alloc_err
  use cam_abortutils,    only: endrun
  !
  use oslo_aero_share,   only: nmodes_oslo => nmodes
  use oslo_aero_share,   only: max_tracers_per_mode, n_tracers_in_mode, n_background_tracers_in_mode
  use oslo_aero_share,   only: CalculateNumberConcentration, calculateNumberMedianRadius
  use oslo_aero_share,   only: MODE_IDX_DST_A2, MODE_IDX_DST_A3, MODE_IDX_OMBC_INTMIX_COAT_AIT
  use oslo_aero_share,   only: getNumberOfTracersInMode, getTracerIndex
  use oslo_aero_share,   only: qqcw_get_field
  use oslo_aero_share,   only: l_dst_a2, l_dst_a3, l_bc_ai, l_bc_ac
  use oslo_aero_share,   only: lifeCycleNumberMedianRadius, lifeCycleSigma

  implicit none
  private

  ! The following are called by microp_aero
  public :: hetfrz_classnuc_oslo_readnl
  public :: hetfrz_classnuc_oslo_register
  public :: hetfrz_classnuc_oslo_init
  public :: hetfrz_classnuc_oslo_calc
  public :: hetfrz_classnuc_oslo_save_cbaero

  private :: get_aer_num
  private :: hetfrz_classnuc_calc
  private :: collkernel
  private :: hetfrz_classnuc_init_pdftheta

  ! Namelist variables
  logical  :: hist_hetfrz_classnuc = .false.
  real(r8) :: hetfrz_bc_scalfac = -huge(1._r8) ! scaling factor for BC  (NOT USED)
  real(r8) :: hetfrz_dust_scalfac = -huge(1._r8) ! scaling factor for dust (NOT USED)

  ! Vars set via init method.
  real(r8) :: mincld      ! minimum allowed cloud fraction

  ! constituent indices
  integer :: cldliq_idx = -1
  integer :: cldice_idx = -1
  integer :: numliq_idx = -1
  integer :: numice_idx = -1

  ! pbuf indices for fields provided by heterogeneous freezing
  integer :: frzimm_idx
  integer :: frzcnt_idx
  integer :: frzdep_idx

  ! pbuf indices for fields needed by heterogeneous freezing
  integer :: ast_idx = -1

  ! Copy of cloud borne aerosols before modification by droplet nucleation
  ! The basis is converted from mass to volume.
  real(r8), allocatable :: aer_cb(:,:,:,:)

  ! PDF theta model
  ! some variables for PDF theta model
  ! immersion freezing
  !
  ! With the original value of pdf_n_theta set to 101 the dust activation
  ! fraction between -15 and 0 C could be overestimated.  This problem was
  ! eliminated by increasing pdf_n_theta to 301.  To reduce the expense of
  ! computing the dust activation fraction the integral is only evaluated
  ! where dim_theta is non-zero.  This was determined to be between
  ! dim_theta index values of 53 through 113.  These loop bounds are
  ! hardcoded in the variables i1 and i2.

  integer, parameter :: pdf_n_theta = 301
  integer, parameter :: i1 = 53
  integer, parameter :: i2 = 113
  real(r8) :: dim_theta(pdf_n_theta) = 0.0_r8
  real(r8) :: pdf_imm_theta(pdf_n_theta) = 0.0_r8
  real(r8) :: pdf_d_theta
  real(r8) :: dim_f_imm_dust_a1(pdf_n_theta) = 0.0_r8
  real(r8) :: dim_f_imm_dust_a3(pdf_n_theta) = 0.0_r8

  !********************************************************
  ! Wang et al., 2014 fitting parameters
  !********************************************************
  ! freezing parameters for immersion freezing
  real(r8),parameter :: theta_imm_bc = 48.0_r8      ! contact angle [deg], converted to rad later !DeMott et al (1990)
  real(r8),parameter :: dga_imm_bc = 14.15E-20_r8   ! activation energy [J]
  real(r8),parameter :: theta_imm_dust = 46.0_r8    ! contact angle [deg], converted to rad later !DeMott et al (2011) SD
  real(r8),parameter :: dga_imm_dust = 14.75E-20_r8 ! activation energy [J]

  ! freezing parameters for deposition nucleation
  real(r8),parameter :: theta_dep_dust = 20.0_r8   ! contact angle [deg], converted to rad later !Koehler et al (2010) SD
  real(r8),parameter :: dga_dep_dust = -8.1E-21_r8 ! activation energy [J]
  real(r8),parameter :: theta_dep_bc = 28._r8      ! contact angle [deg], converted to rad later !Moehler et al (2005), soot
  real(r8),parameter :: dga_dep_bc = -2.E-19_r8    ! activation energy [J]

  real(r8) :: f_imm_bc
  real(r8) :: f_imm_dust_a1, f_imm_dust_a3
  real(r8) :: f_dep_bc
  real(r8) :: f_dep_dust_a1, f_dep_dust_a3
  real(r8) :: f_cnt_bc
  real(r8) :: f_cnt_dust_a1,f_cnt_dust_a3

  real(r8) :: sqroot_f_imm_bc
  real(r8) :: sqroot_dim_f_imm_dust_a1(i1:i2)
  real(r8) :: sqroot_dim_f_imm_dust_a3(i1:i2)
  real(r8) :: sqroot_f_dep_bc
  real(r8) :: sqroot_f_dep_dust_a1
  real(r8) :: sqroot_f_dep_dust_a3
  real(r8) :: sqroot_f_imm_dust_a1
  real(r8) :: sqroot_f_imm_dust_a3

  integer  :: num_bc_idx, num_dst1_idx, num_dst3_idx    ! mode indices
  real(r8) :: fac_volsfc_bc
  real(r8) :: fac_volsfc_dust_a1
  real(r8) :: fac_volsfc_dust_a3

  real(r8) :: con1

  integer  :: tracer_index(0:nmodes_oslo,max_tracers_per_mode)  ! tracer index

  logical  :: pdf_imm_in = .true.
  logical  :: tot_in = .false.

!===============================================================================
contains
!===============================================================================

  subroutine hetfrz_classnuc_oslo_readnl(nlfile)

    use namelist_utils, only: find_group_name

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'hetfrz_classnuc_cam_readnl'

    namelist /hetfrz_classnuc_nl/ hist_hetfrz_classnuc, hetfrz_bc_scalfac, hetfrz_dust_scalfac
    !-----------------------------------------------------------------------------

    if (masterproc) then
       open( newunit=unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'hetfrz_classnuc_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, hetfrz_classnuc_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
    end if

    call mpi_bcast(hist_hetfrz_classnuc, 1, mpi_logical, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hist_hetfrz_classnuc")
    call mpi_bcast(hetfrz_bc_scalfac, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hetfrz_bc_scalfac")
    call mpi_bcast(hetfrz_dust_scalfac, 1, mpi_real8, mstrid, mpicom, ierr)
    if (ierr /= mpi_success) call endrun(subname//" mpi_bcast: hetfrz_dust_scalfac")

    if (masterproc) then
       write(iulog,*) subname,': hist_hetfrz_classnuc = ',hist_hetfrz_classnuc
       write(iulog,*) subname,': hetfrz_bc_scalfac = ',hetfrz_bc_scalfac
       write(iulog,*) subname,': hetfrz_dust_scalfac = ',hetfrz_dust_scalfac
    end if

  end subroutine hetfrz_classnuc_oslo_readnl

  !================================================================================================

  subroutine hetfrz_classnuc_oslo_register()

    if (.not. use_hetfrz_classnuc) return

    ! pbuf fields provided by hetfrz_classnuc
    call pbuf_add_field('FRZIMM', 'physpkg', dtype_r8, (/pcols,pver/), frzimm_idx)
    call pbuf_add_field('FRZCNT', 'physpkg', dtype_r8, (/pcols,pver/), frzcnt_idx)
    call pbuf_add_field('FRZDEP', 'physpkg', dtype_r8, (/pcols,pver/), frzdep_idx)

  end subroutine hetfrz_classnuc_oslo_register

  !================================================================================================

  subroutine hetfrz_classnuc_oslo_init(mincld_in)

    real(r8), intent(in) :: mincld_in

    ! local variables
    integer  :: imode, nspec, itrac
    integer  :: istat
    real(r8) :: sigma_logr_aer
    real(r8) :: sigmag_amode(3)
    real(r8) :: cos_angle
    character(len=32) :: str32
    character(len=*), parameter :: routine = 'hetfrz_classnuc_cam_init'
    !--------------------------------------------------------------------------------------------

    ! This parameterization currently assumes that prognostic modal aerosols are on.  Check...

    if (.not. use_hetfrz_classnuc) return

    mincld = mincld_in

    call cnst_get_ind('CLDLIQ', cldliq_idx)
    call cnst_get_ind('CLDICE', cldice_idx)
    call cnst_get_ind('NUMLIQ', numliq_idx)
    call cnst_get_ind('NUMICE', numice_idx)

    ! pbuf fields used by hetfrz_classnuc
    ast_idx = pbuf_get_index('AST')

    call addfld('bc_num',        (/ 'lev' /), 'A', '#/cm3', 'total bc number')
    call addfld('dst1_num',      (/ 'lev' /), 'A', '#/cm3', 'total dst1 number')
    call addfld('dst3_num',      (/ 'lev' /), 'A', '#/cm3', 'total dst3 number')
    call addfld('bcc_num',       (/ 'lev' /), 'A', '#/cm3', 'coated bc number')
    call addfld('dst1c_num',     (/ 'lev' /), 'A', '#/cm3', 'coated dst1 number')
    call addfld('dst3c_num',     (/ 'lev' /), 'A', '#/cm3', 'coated dst3 number')
    call addfld('bcuc_num',      (/ 'lev' /), 'A', '#/cm3', 'uncoated bc number')
    call addfld('dst1uc_num',    (/ 'lev' /), 'A', '#/cm3', 'uncoated dst1 number')
    call addfld('dst3uc_num',    (/ 'lev' /), 'A', '#/cm3', 'uncoated dst3 number')

    call addfld('bc_a1_num',     (/ 'lev' /), 'A', '#/cm3', 'interstitial bc number')
    call addfld('dst_a1_num',    (/ 'lev' /), 'A', '#/cm3', 'interstitial dst1 number')
    call addfld('dst_a3_num',    (/ 'lev' /), 'A', '#/cm3', 'interstitial dst3 number')
    call addfld('bc_c1_num',     (/ 'lev' /), 'A', '#/cm3', 'cloud borne bc number')
    call addfld('dst_c1_num',    (/ 'lev' /), 'A', '#/cm3', 'cloud borne dst1 number')
    call addfld('dst_c3_num',    (/ 'lev' /), 'A', '#/cm3', 'cloud borne dst3 number')

    call addfld('fn_bc_c1_num',  (/ 'lev' /), 'A', '#/cm3', 'cloud borne bc number derived from fn')
    call addfld('fn_dst_c1_num', (/ 'lev' /), 'A', '#/cm3', 'cloud borne dst1 number derived from fn')
    call addfld('fn_dst_c3_num', (/ 'lev' /), 'A', '#/cm3', 'cloud borne dst3 number derived from fn')

    call addfld('na500',         (/ 'lev' /), 'A', '#/cm3', 'interstitial aerosol number with D>500 nm')
    call addfld('totna500',      (/ 'lev' /), 'A', '#/cm3', 'total aerosol number with D>500 nm')

    call addfld('FREQIMM', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of immersion  freezing')
    call addfld('FREQCNT', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of contact    freezing')
    call addfld('FREQDEP', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of deposition freezing')
    call addfld('FREQMIX', (/ 'lev' /), 'A', 'fraction', 'Fractional occurance of mixed-phase clouds' )

    call addfld('DSTFREZIMM', (/ 'lev' /), 'A', 'm-3s-1', 'dust immersion  freezing rate')
    call addfld('DSTFREZCNT', (/ 'lev' /), 'A', 'm-3s-1', 'dust contact    freezing rate')
    call addfld('DSTFREZDEP', (/ 'lev' /), 'A', 'm-3s-1', 'dust deposition freezing rate')

    call addfld('BCFREZIMM', (/ 'lev' /), 'A', 'm-3s-1', 'bc immersion  freezing rate')
    call addfld('BCFREZCNT', (/ 'lev' /), 'A', 'm-3s-1', 'bc contact    freezing rate')
    call addfld('BCFREZDEP', (/ 'lev' /), 'A', 'm-3s-1', 'bc deposition freezing rate')

    call addfld('NIMIX_IMM', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to het immersion freezing in Mixed Clouds')
    call addfld('NIMIX_CNT', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to het contact freezing in Mixed Clouds')
    call addfld('NIMIX_DEP', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to het deposition freezing in Mixed Clouds')

    call addfld('DSTNIDEP', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to dst dep freezing in Mixed Clouds')
    call addfld('DSTNICNT', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to dst cnt freezing in Mixed Clouds')
    call addfld('DSTNIIMM', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to dst imm freezing in Mixed Clouds')

    call addfld('BCNIDEP', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to bc dep freezing in Mixed Clouds')
    call addfld('BCNICNT', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to bc cnt freezing in Mixed Clouds')
    call addfld('BCNIIMM', (/ 'lev' /), 'A', '#/m3', &
         'Activated Ice Number Concentration due to bc imm freezing in Mixed Clouds')

    call addfld('NUMICE10s', (/ 'lev' /), 'A', '#/m3', &
         'Ice Number Concentration due to het freezing in Mixed Clouds during 10-s period')
    call addfld('NUMIMM10sDST', (/ 'lev' /), 'A', '#/m3', &
         'Ice Number Concentration due to imm freezing by dst in Mixed Clouds during 10-s period')
    call addfld('NUMIMM10sBC', (/ 'lev' /), 'A', '#/m3', &
         'Ice Number Concentration due to imm freezing by bc in Mixed Clouds during 10-s period')

    if (hist_hetfrz_classnuc) then

       call add_default('bc_num', 1, ' ')
       call add_default('dst1_num', 1, ' ')
       call add_default('dst3_num', 1, ' ')
       call add_default('bcc_num', 1, ' ')
       call add_default('dst1c_num', 1, ' ')
       call add_default('dst3c_num', 1, ' ')
       call add_default('bcuc_num', 1, ' ')
       call add_default('dst1uc_num', 1, ' ')
       call add_default('dst3uc_num', 1, ' ')

       call add_default('bc_a1_num', 1, ' ')
       call add_default('dst_a1_num', 1, ' ')
       call add_default('dst_a3_num', 1, ' ')
       call add_default('bc_c1_num', 1, ' ')
       call add_default('dst_c1_num', 1, ' ')
       call add_default('dst_c3_num', 1, ' ')

       call add_default('fn_bc_c1_num', 1, ' ')
       call add_default('fn_dst_c1_num', 1, ' ')
       call add_default('fn_dst_c3_num', 1, ' ')

       call add_default('na500', 1, ' ')
       call add_default('totna500', 1, ' ')

       call add_default('FREQIMM', 1, ' ')
       call add_default('FREQCNT', 1, ' ')
       call add_default('FREQDEP', 1, ' ')
       call add_default('FREQMIX', 1, ' ')

       call add_default('DSTFREZIMM', 1, ' ')
       call add_default('DSTFREZCNT', 1, ' ')
       call add_default('DSTFREZDEP', 1, ' ')

       call add_default('BCFREZIMM', 1, ' ')
       call add_default('BCFREZCNT', 1, ' ')
       call add_default('BCFREZDEP', 1, ' ')

       call add_default('NIMIX_IMM', 1, ' ')
       call add_default('NIMIX_CNT', 1, ' ')
       call add_default('NIMIX_DEP', 1, ' ')

       call add_default('DSTNIDEP', 1, ' ')
       call add_default('DSTNICNT', 1, ' ')
       call add_default('DSTNIIMM', 1, ' ')

       call add_default('BCNIDEP', 1, ' ')
       call add_default('BCNICNT', 1, ' ')
       call add_default('BCNIIMM', 1, ' ')

       call add_default('NUMICE10s', 1, ' ')
       call add_default('NUMIMM10sDST', 1, ' ')
       call add_default('NUMIMM10sBC', 1, ' ')

    end if

    ! The following code sets indices of the mode specific species used
    ! in the module.  Having a list of the species needed allows us to
    ! allocate temporary space for just those species rather than for all the
    ! CAM species (pcnst) which may be considerably more than needed.
    !
    ! The indices set below are for use with the CAM rad_constituents
    ! interfaces.  Using the rad_constituents interfaces isolates the physics
    ! parameterization which requires constituent information from the chemistry
    ! code which provides that information.

    ! Allocate space for copy of cloud borne aerosols before modification by droplet nucleation.
    allocate(aer_cb(pcols,pver,pcnst,begchunk:endchunk), stat=istat)
    call alloc_err(istat, routine, 'aer_cb', pcols*pver*pcnst*(endchunk-begchunk+1))

    ! Initialize all the PDF theta variables:
    ! With the original value of pdf_n_theta set to 101 the dust activation
    ! fraction between -15 and 0 C could be overestimated.  This problem was
    ! eliminated by increasing pdf_n_theta to 301.  To reduce the expense of
    ! computing the dust activation fraction the integral is only evaluated
    ! where dim_theta is non-zero.  This was determined to be between
    ! dim_theta index values of 53 through 113.  These loop bounds are
    ! hardcoded in the variables i1 and i2.

    if (pdf_imm_in) then
       call hetfrz_classnuc_init_pdftheta()
    end if

    ! form factor
    ! only consider flat surfaces due to uncertainty of curved surfaces

    cos_angle = COS(theta_imm_bc*pi/180._r8)
    f_imm_bc = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8
    if (.not. pdf_imm_in) then
       cos_angle = COS(theta_imm_dust*pi/180._r8)
       f_imm_dust_a1 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

       cos_angle = COS(theta_imm_dust*pi/180._r8)
       f_imm_dust_a3 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8
    end if


    ! form factor
    cos_angle = COS(theta_dep_bc*pi/180._r8)
    f_dep_bc = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    cos_angle = COS(theta_dep_dust*pi/180._r8)
    f_dep_dust_a1 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    cos_angle = COS(theta_dep_dust*pi/180._r8)
    f_dep_dust_a3 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    ! form factor
    cos_angle = COS(theta_dep_bc*pi/180._r8)
    f_cnt_bc = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    cos_angle = COS(theta_dep_dust*pi/180._r8)
    f_cnt_dust_a1 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    cos_angle = COS(theta_dep_dust*pi/180._r8)
    f_cnt_dust_a3 = (2._r8+cos_angle)*(1._r8-cos_angle)**2/4._r8

    sqroot_f_imm_bc = SQRT(f_imm_bc)
    sqroot_dim_f_imm_dust_a1(i1:i2) = SQRT(dim_f_imm_dust_a1(i1:i2))
    sqroot_dim_f_imm_dust_a3(i1:i2) = SQRT(dim_f_imm_dust_a3(i1:i2))
    sqroot_f_dep_bc = SQRT(f_dep_bc)
    sqroot_f_dep_dust_a1 = SQRT(f_dep_dust_a1)
    sqroot_f_dep_dust_a3 = SQRT(f_dep_dust_a3)

    num_bc_idx   = MODE_IDX_OMBC_INTMIX_COAT_AIT
    num_dst1_idx = MODE_IDX_DST_A2
    num_dst3_idx = MODE_IDX_DST_A3

    sigmag_amode(1) = lifeCycleSigma(num_bc_idx)
    sigmag_amode(2) = lifeCycleSigma(num_dst1_idx)
    sigmag_amode(3) = lifeCycleSigma(num_dst3_idx)

    fac_volsfc_bc      = exp(2.5*(log(sigmag_amode(1))**2))
    fac_volsfc_dust_a1 = exp(2.5*(log(sigmag_amode(2))**2))
    fac_volsfc_dust_a3 = exp(2.5*(log(sigmag_amode(3))**2))

    con1 = 1._r8/(1.333_r8*pi)**0.333_r8

    do imode = 0,nmodes_oslo
       do itrac = 1,n_tracers_in_mode(imode)
          tracer_index(imode,itrac) = getTracerIndex(imode,itrac,.false.)
       end do
    end do

  end subroutine hetfrz_classnuc_oslo_init

  !================================================================================================

  subroutine hetfrz_classnuc_oslo_calc(          &
       state, deltatin, factnum, pbuf,           &
       numberConcentration, volumeConcentration, &
       f_acm, f_bcm, f_aqm, f_so4_condm, f_soam, &
       hygroscopicity, lnsigma, cam, volumeCore, volumeCoat)

    ! arguments
    type(physics_state), target, intent(in) :: state
    real(r8),                    intent(in) :: deltatin       ! time step (s)
    real(r8),                    intent(in) :: factnum(:,:,:) ! activation fraction for aerosol number
    type(physics_buffer_desc),   pointer    :: pbuf(:)
    real(r8),                    intent(in) :: numberConcentration(pcols,pver,0:nmodes_oslo)
    real(r8),                    intent(in) :: volumeConcentration(pcols,pver,nmodes_oslo)
    real(r8),                    intent(in) :: f_acm(pcols,pver, nmodes_oslo)
    real(r8),                    intent(in) :: f_bcm(pcols,pver, nmodes_oslo)
    real(r8),                    intent(in) :: f_aqm(pcols, pver, nmodes_oslo)
    real(r8),                    intent(in) :: f_so4_condm(pcols, pver, nmodes_oslo)            !Needed in "get component fraction"
    real(r8),                    intent(in) :: f_soam(pcols, pver, nmodes_oslo)
    real(r8),                    intent(in) :: hygroscopicity(pcols,pver,nmodes_oslo)        ![mol_{aer}/mol_{water}] hygroscopicity
    real(r8),                    intent(in) :: lnsigma(pcols,pver,nmodes_oslo)              ![-] log(base e) sigma
    real(r8),                    intent(in) :: cam(pcols,pver,nmodes_oslo)
    real(r8),                    intent(in) :: volumeCore(pcols,pver,nmodes_oslo)
    real(r8),                    intent(in) :: volumeCoat(pcols,pver,nmodes_oslo)

    ! local workspace
    real(r8), pointer :: frzimm(:,:) ! output shared with the microphysics via the pbuf
    real(r8), pointer :: frzcnt(:,:) ! output shared with the microphysics via the pbuf
    real(r8), pointer :: frzdep(:,:) ! output shared with the microphysics via the pbuf
    real(r8), pointer :: ast(:,:)
    integer  :: itim_old
    integer  :: icol, ilev, cos_angle
    real(r8) :: rho(pcols,pver)          ! air density (kg m-3)
    real(r8) :: lcldm(pcols,pver)
    real(r8) :: fn(3)
    real(r8) :: awcam(pcols,pver,3)
    real(r8) :: awfacm(pcols,pver,3)
    real(r8) :: hetraer(pcols,pver,3)
    real(r8) :: dstcoat(pcols,pver,3)
    real(r8) :: total_interstitial_aer_num(pcols,pver,3)
    real(r8) :: total_cloudborne_aer_num(pcols,pver,3)
    real(r8) :: total_aer_num(pcols,pver,3)
    real(r8) :: coated_aer_num(pcols,pver,3)
    real(r8) :: uncoated_aer_num(pcols,pver,3)
    real(r8) :: fn_cloudborne_aer_num(pcols,pver,3)
    real(r8) :: r3lx, supersatice(pcols,pver)
    real(r8) :: eswtr(pcols,pver)
    real(r8) :: esice(pcols,pver)
    real(r8) :: qcic
    real(r8) :: ncic
    real(r8) :: frzbcimm(pcols,pver), frzduimm(pcols,pver)
    real(r8) :: frzbccnt(pcols,pver), frzducnt(pcols,pver)
    real(r8) :: frzbcdep(pcols,pver), frzdudep(pcols,pver)
    real(r8) :: freqimm(pcols,pver), freqcnt(pcols,pver), freqdep(pcols,pver), freqmix(pcols,pver)
    real(r8) :: nnuccc_bc(pcols,pver), nnucct_bc(pcols,pver), nnudep_bc(pcols,pver)
    real(r8) :: nnuccc_dst(pcols,pver), nnucct_dst(pcols,pver), nnudep_dst(pcols,pver)
    real(r8) :: niimm_bc(pcols,pver), nicnt_bc(pcols,pver), nidep_bc(pcols,pver)
    real(r8) :: niimm_dst(pcols,pver), nicnt_dst(pcols,pver), nidep_dst(pcols,pver)
    real(r8) :: numice10s(pcols,pver)
    real(r8) :: numice10s_imm_dst(pcols,pver)
    real(r8) :: numice10s_imm_bc(pcols,pver)
    real(r8) :: CloudnumberConcentration(pcols,pver,0:nmodes_oslo) ! oslo aerosol specific
    real(r8) :: numberMedianRadius(pcols,pver,nmodes_oslo) ! oslo aerosol specific
    real(r8) :: na500(pcols,pver)
    real(r8) :: tot_na500(pcols,pver)
    character(128) :: errstring   ! Error status
    !-------------------------------------------------------------------------------

    associate( &
         lchnk => state%lchnk,             &
         ncol  => state%ncol,              &
         t     => state%t,                 &
         qc    => state%q(:pcols,:pver,cldliq_idx), &
         nc    => state%q(:pcols,:pver,numliq_idx), &
         pmid  => state%pmid               )

    itim_old = pbuf_old_tim_idx()
    call pbuf_get_field(pbuf, ast_idx, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

    rho(:,:) = 0._r8

    do ilev = top_lev, pver
       do icol = 1, ncol
          rho(icol,ilev) = pmid(icol,ilev)/(rair*t(icol,ilev))
       end do
    end do

    do ilev = top_lev, pver
       do icol = 1, ncol
          lcldm(icol,ilev) = max(ast(icol,ilev), mincld)
       end do
    end do

    ! Convert interstitial and cloud borne aerosols from a mass to a volume basis before
    ! being used in get_aer_num
    do icol = 1, pcnst
       aer_cb(:ncol,:,icol,lchnk) = aer_cb(:ncol,:,icol,lchnk) * rho(:ncol,:)
    end do

    !Get estimate of number of aerosols inside clouds
    call calculateNumberConcentration(ncol, aer_cb(:,:,:,lchnk), rho, CloudnumberConcentration)
    call calculateNumberMedianRadius(numberConcentration, volumeConcentration, lnSigma, numberMedianRadius, ncol)

    ! output aerosols as reference information for heterogeneous freezing
    call get_aer_num(ncol, numberConcentration, CloudnumberConcentration, rho, &
         f_acm, f_so4_condm, cam, volumeCore, volumeCoat,                &
         total_aer_num, coated_aer_num, uncoated_aer_num,                &
         total_interstitial_aer_num, total_cloudborne_aer_num,           &
         hetraer, awcam, awfacm, dstcoat,                                &
         na500, tot_na500)

    do ilev = top_lev, pver
       do icol = 1, ncol
          fn_cloudborne_aer_num(icol,ilev,1) = total_aer_num(icol,ilev,1)*factnum(icol,ilev,MODE_IDX_OMBC_INTMIX_COAT_AIT)  ! bc
          fn_cloudborne_aer_num(icol,ilev,2) = total_aer_num(icol,ilev,2)*factnum(icol,ilev,MODE_IDX_DST_A2)
          fn_cloudborne_aer_num(icol,ilev,3) = total_aer_num(icol,ilev,3)*factnum(icol,ilev,MODE_IDX_DST_A3)
       end do
    end do

    call outfld('bc_num',        total_aer_num(:,:,1),    pcols, lchnk)
    call outfld('dst1_num',      total_aer_num(:,:,2),    pcols, lchnk)
    call outfld('dst3_num',      total_aer_num(:,:,3),    pcols, lchnk)

    call outfld('bcc_num',       coated_aer_num(:,:,1),   pcols, lchnk)
    call outfld('dst1c_num',     coated_aer_num(:,:,2),   pcols, lchnk)
    call outfld('dst3c_num',     coated_aer_num(:,:,3),   pcols, lchnk)

    call outfld('bcuc_num',      uncoated_aer_num(:,:,1), pcols, lchnk)
    call outfld('dst1uc_num',    uncoated_aer_num(:,:,2), pcols, lchnk)
    call outfld('dst3uc_num',    uncoated_aer_num(:,:,3), pcols, lchnk)

    call outfld('bc_a1_num',     total_interstitial_aer_num(:,:,1), pcols, lchnk)
    call outfld('dst_a1_num',    total_interstitial_aer_num(:,:,2), pcols, lchnk)
    call outfld('dst_a3_num',    total_interstitial_aer_num(:,:,3), pcols, lchnk)

    call outfld('bc_c1_num',     total_cloudborne_aer_num(:,:,1),   pcols, lchnk)
    call outfld('dst_c1_num',    total_cloudborne_aer_num(:,:,2),   pcols, lchnk)
    call outfld('dst_c3_num',    total_cloudborne_aer_num(:,:,3),   pcols, lchnk)

    call outfld('fn_bc_c1_num',  fn_cloudborne_aer_num(:,:,1),      pcols, lchnk)
    call outfld('fn_dst_c1_num', fn_cloudborne_aer_num(:,:,2),      pcols, lchnk)
    call outfld('fn_dst_c3_num', fn_cloudborne_aer_num(:,:,3),      pcols, lchnk)

    call outfld('na500',         na500,     pcols, lchnk)
    call outfld('totna500',      tot_na500, pcols, lchnk)

    ! frzimm, frzcnt, frzdep are the outputs of this parameterization used by the microphysics
    call pbuf_get_field(pbuf, frzimm_idx, frzimm)
    call pbuf_get_field(pbuf, frzcnt_idx, frzcnt)
    call pbuf_get_field(pbuf, frzdep_idx, frzdep)

    frzimm(:ncol,:) = 0._r8
    frzcnt(:ncol,:) = 0._r8
    frzdep(:ncol,:) = 0._r8

    frzbcimm(:ncol,:) = 0._r8
    frzduimm(:ncol,:) = 0._r8
    frzbccnt(:ncol,:) = 0._r8
    frzducnt(:ncol,:) = 0._r8
    frzbcdep(:ncol,:) = 0._r8
    frzdudep(:ncol,:) = 0._r8

    freqimm(:ncol,:) = 0._r8
    freqcnt(:ncol,:) = 0._r8
    freqdep(:ncol,:) = 0._r8
    freqmix(:ncol,:) = 0._r8

    numice10s(:ncol,:)         = 0._r8
    numice10s_imm_dst(:ncol,:) = 0._r8
    numice10s_imm_bc(:ncol,:)  = 0._r8

    nnuccc_bc(:,:) = 0._r8
    nnucct_bc(:,:) = 0._r8
    nnudep_bc(:,:) = 0._r8

    nnuccc_dst(:,:) = 0._r8
    nnucct_dst(:,:) = 0._r8
    nnudep_dst(:,:) = 0._r8

    niimm_bc(:,:) = 0._r8
    nicnt_bc(:,:) = 0._r8
    nidep_bc(:,:) = 0._r8

    niimm_dst(:,:) = 0._r8
    nicnt_dst(:,:) = 0._r8
    nidep_dst(:,:) = 0._r8

    do ilev = top_lev, pver
       do icol = 1,ncol
          if (t(icol,ilev) > 235.15_r8 .and. t(icol,ilev) < 269.15_r8) then
             eswtr(icol,ilev) = svp_water(t(icol,ilev))
             esice(icol,ilev) = svp_ice(t(icol,ilev))
             supersatice(icol,ilev) = eswtr(icol,ilev)/esice(icol,ilev)
          end if
       end do
    end do

    do ilev = top_lev, pver
       do icol = 1, ncol

          if (t(icol,ilev) > 235.15_r8 .and. t(icol,ilev) < 269.15_r8) then
             qcic = min(qc(icol,ilev)/lcldm(icol,ilev), 5.e-3_r8)
             ncic = max(nc(icol,ilev)/lcldm(icol,ilev), 0._r8)

             r3lx = con1*(rho(icol,ilev)*qcic/(rhoh2o*max(ncic*rho(icol,ilev), 1.0e6_r8)))**0.333_r8 ! in m
             r3lx = max(4.e-6_r8, r3lx)
             fn(1) = factnum(icol,ilev,MODE_IDX_OMBC_INTMIX_COAT_AIT)  ! bc accumulation mode
             fn(2) = factnum(icol,ilev,MODE_IDX_DST_A2)                ! dust_a1 accumulation mode
             fn(3) = factnum(icol,ilev,MODE_IDX_DST_A3)                ! dust_a3 coarse mode

             call hetfrz_classnuc_calc( &
                  deltatin,  t(icol,ilev),  pmid(icol,ilev),  supersatice(icol,ilev), eswtr(icol,ilev), esice(icol,ilev),  &
                  fn,  r3lx,  ncic*rho(icol,ilev)*1.0e-6_r8,  frzbcimm(icol,ilev),  frzduimm(icol,ilev),   &
                  frzbccnt(icol,ilev),  frzducnt(icol,ilev),  frzbcdep(icol,ilev),  frzdudep(icol,ilev),  hetraer(icol,ilev,:), &
                  awcam(icol,ilev,:), awfacm(icol,ilev,:), dstcoat(icol,ilev,:), total_aer_num(icol,ilev,:),  &
                  coated_aer_num(icol,ilev,:), uncoated_aer_num(icol,ilev,:), total_interstitial_aer_num(icol,ilev,:), &
                  total_cloudborne_aer_num(icol,ilev,:), errstring)

             call handle_errmsg(errstring, subname="hetfrz_classnuc_calc")

             frzimm(icol,ilev) = frzbcimm(icol,ilev) + frzduimm(icol,ilev)
             frzcnt(icol,ilev) = frzbccnt(icol,ilev) + frzducnt(icol,ilev)
             frzdep(icol,ilev) = frzbcdep(icol,ilev) + frzdudep(icol,ilev)

             if (frzimm(icol,ilev) > 0._r8) freqimm(icol,ilev) = 1._r8
             if (frzcnt(icol,ilev) > 0._r8) freqcnt(icol,ilev) = 1._r8
             if (frzdep(icol,ilev) > 0._r8) freqdep(icol,ilev) = 1._r8
             if ((frzimm(icol,ilev) + frzcnt(icol,ilev) + frzdep(icol,ilev)) > 0._r8) freqmix(icol,ilev) = 1._r8
          else
             frzimm(icol,ilev) = 0._r8
             frzcnt(icol,ilev) = 0._r8
             frzdep(icol,ilev) = 0._r8
          end if

          nnuccc_bc(icol,ilev) = frzbcimm(icol,ilev)*1.0e6_r8*ast(icol,ilev)
          nnucct_bc(icol,ilev) = frzbccnt(icol,ilev)*1.0e6_r8*ast(icol,ilev)
          nnudep_bc(icol,ilev) = frzbcdep(icol,ilev)*1.0e6_r8*ast(icol,ilev)

          nnuccc_dst(icol,ilev) = frzduimm(icol,ilev)*1.0e6_r8*ast(icol,ilev)
          nnucct_dst(icol,ilev) = frzducnt(icol,ilev)*1.0e6_r8*ast(icol,ilev)
          nnudep_dst(icol,ilev) = frzdudep(icol,ilev)*1.0e6_r8*ast(icol,ilev)

          niimm_bc(icol,ilev) = frzbcimm(icol,ilev)*1.0e6_r8*deltatin
          nicnt_bc(icol,ilev) = frzbccnt(icol,ilev)*1.0e6_r8*deltatin
          nidep_bc(icol,ilev) = frzbcdep(icol,ilev)*1.0e6_r8*deltatin

          niimm_dst(icol,ilev) = frzduimm(icol,ilev)*1.0e6_r8*deltatin
          nicnt_dst(icol,ilev) = frzducnt(icol,ilev)*1.0e6_r8*deltatin
          nidep_dst(icol,ilev) = frzdudep(icol,ilev)*1.0e6_r8*deltatin

          numice10s(icol,ilev) = (frzimm(icol,ilev)+frzcnt(icol,ilev)+frzdep(icol,ilev))*1.0e6_r8*deltatin*(10._r8/deltatin)
          numice10s_imm_dst(icol,ilev) = frzduimm(icol,ilev)*1.0e6_r8*deltatin*(10._r8/deltatin)
          numice10s_imm_bc(icol,ilev) = frzbcimm(icol,ilev)*1.0e6_r8*deltatin*(10._r8/deltatin)
       end do
    end do

    call outfld('FREQIMM', freqimm, pcols, lchnk)
    call outfld('FREQCNT', freqcnt, pcols, lchnk)
    call outfld('FREQDEP', freqdep, pcols, lchnk)
    call outfld('FREQMIX', freqmix, pcols, lchnk)

    call outfld('DSTFREZIMM', nnuccc_dst, pcols, lchnk)
    call outfld('DSTFREZCNT', nnucct_dst, pcols, lchnk)
    call outfld('DSTFREZDEP', nnudep_dst, pcols, lchnk)

    call outfld('BCFREZIMM', nnuccc_bc, pcols, lchnk)
    call outfld('BCFREZCNT', nnucct_bc, pcols, lchnk)
    call outfld('BCFREZDEP', nnudep_bc, pcols, lchnk)

    call outfld('NIMIX_IMM', niimm_bc+niimm_dst, pcols, lchnk)
    call outfld('NIMIX_CNT', nicnt_bc+nicnt_dst, pcols, lchnk)
    call outfld('NIMIX_DEP', nidep_bc+nidep_dst, pcols, lchnk)

    call outfld('DSTNICNT', nicnt_dst, pcols, lchnk)
    call outfld('DSTNIDEP', nidep_dst, pcols, lchnk)
    call outfld('DSTNIIMM', niimm_dst, pcols, lchnk)

    call outfld('BCNICNT', nicnt_bc, pcols, lchnk)
    call outfld('BCNIDEP', nidep_bc, pcols, lchnk)
    call outfld('BCNIIMM', niimm_bc, pcols, lchnk)

    call outfld('NUMICE10s', numice10s, pcols, lchnk)
    call outfld('NUMIMM10sDST', numice10s_imm_dst, pcols, lchnk)
    call outfld('NUMIMM10sBC', numice10s_imm_bc, pcols, lchnk)

    end associate

  end subroutine hetfrz_classnuc_oslo_calc

  !====================================================================================================

  subroutine hetfrz_classnuc_oslo_save_cbaero(state, pbuf)

    ! Save the required cloud borne aerosol constituents.
    type(physics_state),         intent(in)    :: state
    type(physics_buffer_desc),   pointer       :: pbuf(:)

    ! local variables
    integer :: icol, ilev, lchnk, kk, ncol, imode, itrac
    type qqcw_type
       real(r8), pointer :: fldcw(:,:)
    end type qqcw_type
    type(qqcw_type) :: qqcw(pcnst)
    !-------------------------------------------------------------------------------

    ! loop over the cloud borne constituents required by this module and save
    ! a local copy

    lchnk = state%lchnk
    ncol = state%ncol
    aer_cb(1:ncol,1:pver,:,lchnk) = 0.0_r8
    do imode = 1,nmodes_oslo
       do itrac = 1,n_tracers_in_mode(imode)
          kk = tracer_index(imode,itrac) ! This gives the tracer index used in the q-array
          qqcw(kk)%fldcw => qqcw_get_field(pbuf,kk)
          if (associated(qqcw(kk)%fldcw)) then
             do ilev = 1,pver
                do icol = 1,ncol
                   aer_cb(icol,ilev,kk,lchnk) = qqcw(kk)%fldcw(icol,ilev)
                end do
             end do
          end if
       end do
    end do

  end subroutine hetfrz_classnuc_oslo_save_cbaero

  !====================================================================================================

  subroutine get_aer_num(                                    &
       ncol, qaerpt, qaercwpt, rhoair,                       &   ! input
       f_acm, f_condm,                                       &
       cam, volumeCore, volumeCoat,                          &
       total_aer_num, coated_aer_num,  uncoated_aer_num,     &   ! output
       total_interstitial_aer_num, total_cloudborne_aer_num, &
       hetraer, awcam, awfacm, dstcoat,                      &
       na500, tot_na500)

    ! input
    integer , intent(in) :: ncol
    real(r8), intent(in) :: qaerpt(pcols,pver,0:nmodes_oslo)   ! aerosol number and mass mixing ratios(instertitial)
    real(r8), intent(in) :: qaercwpt(pcols,pver,0:nmodes_oslo) ! cloud borne aerosol number and mass mixing ratios
    real(r8), intent(in) :: rhoair(pcols,pver)                 ! air density (kg/m3)
    real(r8), intent(in) :: f_acm(pcols,pver,nmodes_oslo)
    real(r8), intent(in) :: f_condm(pcols,pver,nmodes_oslo)
    real(r8), intent(in) :: cam(pcols,pver,nmodes_oslo)
    real(r8), intent(in) :: volumeCoat(pcols,pver,nmodes_oslo)
    real(r8), intent(in) :: volumeCore(pcols,pver,nmodes_oslo)

    ! output
    real(r8), intent(out) :: total_aer_num(pcols,pver,3)              ! #/cm^3
    real(r8), intent(out) :: coated_aer_num(pcols,pver,3)             ! #/cm^3
    real(r8), intent(out) :: uncoated_aer_num(pcols,pver,3)           ! #/cm^3
    real(r8), intent(out) :: total_interstitial_aer_num(pcols,pver,3) ! #/cm^3
    real(r8), intent(out) :: total_cloudborne_aer_num(pcols,pver,3)   ! #/cm^3
    real(r8), intent(out) :: hetraer(pcols,pver,3)                    ! BC and Dust mass mean radius [m]
    real(r8), intent(out) :: awcam(pcols,pver,3)                      ! modal added mass [mug m-3]
    real(r8), intent(out) :: awfacm(pcols,pver,3)                     ! (OC+BC)/(OC+BC+SO4)
    real(r8), intent(out) :: dstcoat(pcols,pver,3)                    ! coated fraction
    real(r8), intent(out) :: na500(pcols,pver)                        ! #/cm^3 interstitial aerosol number with D>500 nm (#/cm^3)
    real(r8), intent(out) :: tot_na500(pcols,pver)                    ! #/cm^3 total aerosol number with D>500 nm (#/cm^3)

    ! local variables
    real(r8), parameter :: n_so4_monolayers_dust = 1.0_r8 ! number of so4(+nh4) monolayers needed to coat a dust particle
    real(r8), parameter :: dr_so4_monolayers_dust = n_so4_monolayers_dust * 4.76e-10
    real(r8) :: tmp1, tmp2
    real(r8) :: bc_num                                    ! bc number in accumulation mode
    real(r8) :: dst1_num, dst3_num                        ! dust number in accumulation and corase mode
    real(r8) :: dst1_num_imm, dst3_num_imm, bc_num_imm
    real(r8) :: r_bc                         ! model radii of BC modes [m]
    real(r8) :: r_dust_a1, r_dust_a3         ! model radii of dust modes [m]
    integer  :: icol,ilev,imode

    ! calculate mass mean radius
    r_dust_a1 = lifeCycleNumberMedianRadius(num_dst1_idx)
    r_dust_a3 = lifeCycleNumberMedianRadius(num_dst3_idx)
    r_bc      = lifeCycleNumberMedianRadius(num_bc_idx)

    ! set hetraer
    hetraer(:,:,1) = r_bc
    hetraer(:,:,2) = r_dust_a1
    hetraer(:,:,3) = r_dust_a3

    ! init output
    total_aer_num(:,:,:)              = 0._r8
    coated_aer_num(:,:,:)             = 0._r8
    uncoated_aer_num(:,:,:)           = 0._r8
    total_interstitial_aer_num(:,:,:) = 0._r8
    total_cloudborne_aer_num(:,:,:)   = 0._r8
    awcam(:,:,:)                      = 0._r8
    awfacm(:,:,:)                     = 0._r8
    dstcoat(:,:,:)                    = 0._r8
    na500(:,:)                        = 0._r8
    tot_na500(:,:)                    = 0._r8

    do ilev = top_lev, pver
       do icol = 1, ncol

          ! calculate intersitial aerosol
          dst1_num = qaerpt(icol,ilev,num_dst1_idx)*1.0e-6_r8    ! #/cm3
          dst3_num = qaerpt(icol,ilev,num_dst3_idx)*1.0e-6_r8    ! #/cm3
          bc_num   = qaerpt(icol,ilev,num_bc_idx)*1.0e-6_r8      ! #/cm3

          ! calculate cloud borne aerosol
          dst1_num_imm = qaercwpt(icol,ilev,num_dst1_idx)*1.0e-6_r8  ! #/cm3
          dst3_num_imm = qaercwpt(icol,ilev,num_dst3_idx)*1.0e-6_r8  ! #/cm3
          bc_num_imm   = qaercwpt(icol,ilev,num_bc_idx)*1.0e-6_r8    ! #/cm3

          ! calculate coated fraction
          ! volumeCore and volumeCoat from subroutine calculateHygroscopicity in paramix_progncdnc.f90

          tmp1 = volumeCoat(icol,ilev,num_bc_idx)*(r_bc*2._r8)*fac_volsfc_bc
          tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(icol,ilev,num_bc_idx), 0.0_r8)
          dstcoat(icol,ilev,1) = tmp1/tmp2
          if (dstcoat(icol,ilev,1) > 1._r8) then
             dstcoat(icol,ilev,1) = 1._r8
          else if (dstcoat(icol,ilev,1) < 0.001_r8) then
             dstcoat(icol,ilev,1) = 0.001_r8
          end if

          tmp1 = volumeCoat(icol,ilev,num_dst1_idx)*(r_dust_a1*2._r8)*fac_volsfc_dust_a1
          tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(icol,ilev,num_dst1_idx), 0.0_r8)
          dstcoat(icol,ilev,2) = tmp1/tmp2
          if (dstcoat(icol,ilev,2) > 1._r8) then
             dstcoat(icol,ilev,2) = 1._r8
          else if (dstcoat(icol,ilev,2) < 0.001_r8) then
             dstcoat(icol,ilev,2) = 0.001_r8
          end if

          tmp1 = volumeCoat(icol,ilev,num_dst3_idx)*(r_dust_a3*2._r8)*fac_volsfc_dust_a3
          tmp2 = max(6.0_r8*dr_so4_monolayers_dust*volumeCore(icol,ilev,num_dst3_idx), 0.0_r8)
          dstcoat(icol,ilev,3) = tmp1/tmp2
          if (dstcoat(icol,ilev,3) > 1._r8) then
             dstcoat(icol,ilev,3) = 1._r8
          else if (dstcoat(icol,ilev,3) < 0.001_r8) then
             dstcoat(icol,ilev,3) = 0.001_r8
          end if

          ! prepare some variables for water activity

          ! cam ([kg/m3] added mass distributed to modes)
          ! accumulation mode for dust_a1
          if (qaerpt(icol,ilev,num_dst1_idx) > 0._r8) then
             awcam(icol,ilev,2) = cam(icol,ilev,num_dst1_idx)*1.e9_r8    ! kg/m3 -> ug/m3
          else
             awcam(icol,ilev,2) = 0._r8
          end if
          if (awcam(icol,ilev,2) >0._r8) then
             awfacm(icol,ilev,2) = f_acm(icol,ilev,num_dst1_idx)
          else
             awfacm(icol,ilev,2) = 0._r8
          end if

          ! accumulation mode for dust_a3
          if (qaerpt(icol,ilev,num_dst3_idx) > 0._r8) then
             awcam(icol,ilev,3) = cam(icol,ilev,num_dst3_idx)*1.e9_r8    ! kg/m3 -> ug/m3
          else
             awcam(icol,ilev,3) = 0._r8
          end if
          if (awcam(icol,ilev,3) >0._r8) then
             awfacm(icol,ilev,3) = f_acm(icol,ilev,num_dst3_idx)
          else
             awfacm(icol,ilev,3) = 0._r8
          end if

          ! accumulation mode for bc
          if (qaerpt(icol,ilev,num_bc_idx) > 0._r8) then
             awcam(icol,ilev,1) = cam(icol,ilev,num_bc_idx)*1.e9_r8    ! kg/m3 -> ug/m3
          else
             awcam(icol,ilev,1) = 0._r8
          end if
          if (awcam(icol,ilev,1) >0._r8) then
             awfacm(icol,ilev,1) = f_acm(icol,ilev,num_bc_idx)
          else
             awfacm(icol,ilev,1) = 0._r8
          end if

          ! prepare output
          total_interstitial_aer_num(icol,ilev,1) = bc_num
          total_interstitial_aer_num(icol,ilev,2) = dst1_num
          total_interstitial_aer_num(icol,ilev,3) = dst3_num

          total_cloudborne_aer_num(icol,ilev,1) = bc_num_imm
          total_cloudborne_aer_num(icol,ilev,2) = dst1_num_imm
          total_cloudborne_aer_num(icol,ilev,3) = dst3_num_imm

          do imode = 1, 3
             total_aer_num(icol,ilev,imode)    = &
                  total_interstitial_aer_num(icol,ilev,imode) + total_cloudborne_aer_num(icol,ilev,imode)
             coated_aer_num(icol,ilev,imode)   = &
                  total_interstitial_aer_num(icol,ilev,imode) * dstcoat(icol,ilev,imode)
             uncoated_aer_num(icol,ilev,imode) = &
                  total_interstitial_aer_num(icol,ilev,imode) * (1._r8-dstcoat(icol,ilev,imode))
          end do

          ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
          tot_na500(icol,ilev) = total_aer_num(icol,ilev,1)*0.0256_r8  &
                          +total_aer_num(icol,ilev,3)

          ! scaled for D>0.5 um using Clarke et al., 1997; 2004; 2007: rg=0.1um, sig=1.6
          na500(icol,ilev) = total_interstitial_aer_num(icol,ilev,1)*0.0256_r8   &
                      +total_interstitial_aer_num(icol,ilev,3)

       end do
    end do

  end subroutine get_aer_num

  !===================================================================================================

  subroutine hetfrz_classnuc_calc(                                   &
       deltat, t, p, supersatice, eswtr, esice,                      &
       fn, r3lx, icnlx, frzbcimm, frzduimm,                          &
       frzbccnt, frzducnt,                                           &
       frzbcdep, frzdudep, hetraer,                                  &
       awcam, awfacm, dstcoat, total_aer_num,                        &
       coated_aer_num, uncoated_aer_num, total_interstitial_aer_num, &
       total_cloudborne_aer_num, errstring)

    real(r8), intent(in) :: deltat                        ! timestep [s]
    real(r8), intent(in) :: t                             ! temperature [K]
    real(r8), intent(in) :: p                             ! pressure [Pa]
    real(r8), intent(in) :: supersatice                   ! supersaturation ratio wrt ice at 100%rh over water [ ]
    real(r8), intent(in) :: eswtr                         ! saturation vapor pressure water [Pa]
    real(r8), intent(in) :: esice                         ! saturation vapor pressure ice [Pa]
    real(r8), intent(in) :: r3lx                          ! volume mean drop radius [m]
    real(r8), intent(in) :: icnlx                         ! in-cloud droplet concentration [cm-3]
    real(r8), intent(in) :: fn(3)                         ! fraction activated [ ] for cloud borne aerosol number
                                                          ! index values are 1:bc, 2:dust_a1, 3:dust_a3
    real(r8), intent(in) :: hetraer(3)                    ! bc and dust mass mean radius [m]
    real(r8), intent(in) :: awcam(3)                      ! modal added mass [mug m-3]
    real(r8), intent(in) :: awfacm(3)                     ! (OC+BC)/(OC+BC+SO4)
    real(r8), intent(in) :: dstcoat(3)                    ! coated fraction
    real(r8), intent(in) :: total_aer_num(3)              ! total bc and dust number concentration(interstitial+cloudborne) [#/cm^3]
    real(r8), intent(in) :: coated_aer_num(3)             ! coated bc and dust number concentration(interstitial)
    real(r8), intent(in) :: uncoated_aer_num(3)           ! uncoated bc and dust number concentration(interstitial)
    real(r8), intent(in) :: total_interstitial_aer_num(3) ! total bc and dust concentration(interstitial)
    real(r8), intent(in) :: total_cloudborne_aer_num(3)   ! total bc and dust concentration(cloudborne)

    real(r8), intent(out) :: frzbcimm                     ! het. frz by BC immersion nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzduimm                     ! het. frz by dust immersion nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzbccnt                     ! het. frz by BC contact nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzducnt                     ! het. frz by dust contact nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzbcdep                     ! het. frz by BC deposition nucleation [cm-3 s-1]
    real(r8), intent(out) :: frzdudep                     ! het. frz by dust deposition nucleation [cm-3 s-1]
    character(len=*), intent(out) :: errstring

    ! local variables
    real(r8) , parameter :: Mso4 = 96.06_r8
    integer  , parameter :: id_bc   = 1
    integer  , parameter :: id_dst1 = 2
    integer  , parameter :: id_dst3 = 3
    real(r8) , parameter :: n1 = 1.e19_r8        ! number of water molecules in contact with unit area of substrate [m-2]
    real(r8) , parameter :: kboltz = 1.38e-23_r8
    real(r8) , parameter :: hplanck = 6.63e-34_r8
    real(r8) , parameter :: rhplanck = 1._r8/hplanck
    real(r8) , parameter :: amu = 1.66053886e-27_r8
    real(r8) , parameter :: nus = 1.e13_r8       ! frequ. of vibration [s-1] higher freq. (as in P&ILEV, consistent with Anupam's data)
    real(r8) , parameter :: taufrz = 195.435_r8  ! time constant for falloff of freezing rate [s]
    real(r8) , parameter :: rhwincloud = 0.98_r8 ! 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
    real(r8) , parameter :: limfacbc = 0.01_r8   ! max. ice nucleating fraction soot
    real(r8) :: aw(3)                            ! water activity [ ]
    real(r8) :: molal(3)                         ! molality [moles/kg]
    logical  :: do_bc, do_dst1, do_dst3
    real(r8) :: tc
    real(r8) :: vwice
    real(r8) :: rhoice
    real(r8) :: sigma_iw                        ! [J/m2]
    real(r8) :: sigma_iv                        ! [J/m2]
    real(r8) :: rgimm
    real(r8) :: rgdep
    real(r8) :: dg0dep
    real(r8) :: Adep
    real(r8) :: dg0cnt
    real(r8) :: Acnt
    real(r8) :: rgimm_bc
    real(r8) :: rgimm_dust_a1, rgimm_dust_a3
    real(r8) :: dg0imm_bc
    real(r8) :: dg0imm_dust_a1, dg0imm_dust_a3
    real(r8) :: Aimm_bc
    real(r8) :: Aimm_dust_a1, Aimm_dust_a3
    real(r8) :: q,phi
    real(r8) :: r_bc                            ! model radii of BC modes [m]
    real(r8) :: r_dust_a1, r_dust_a3            ! model radii of dust modes [m]
    real(r8) :: Jimm_bc
    real(r8) :: Jimm_dust_a1, Jimm_dust_a3
    real(r8) :: Jdep_bc
    real(r8) :: Jdep_dust_a1, Jdep_dust_a3
    real(r8) :: Jcnt_bc
    real(r8) :: Jcnt_dust_a1,Jcnt_dust_a3
    integer  :: icol

    real(r8) :: Kcoll_bc      ! collision kernel [cm3 s-1]
    real(r8) :: Kcoll_dust_a1 ! collision kernel [cm3 s-1]
    real(r8) :: Kcoll_dust_a3 ! collision kernel [cm3 s-1]
    real(r8) :: dim_Jimm_dust_a1(pdf_n_theta), dim_Jimm_dust_a3(pdf_n_theta)
    real(r8) :: sum_imm_dust_a1, sum_imm_dust_a3
    !------------------------------------------------------------------------------------------------

    tc = t - tmelt
    rhoice = 916.7_r8-0.175_r8*tc-5.e-4_r8*tc**2
    vwice = mwh2o*amu/rhoice
    sigma_iw = (28.5_r8+0.25_r8*tc)*1E-3_r8
    sigma_iv = (76.1_r8-0.155_r8*tc + 28.5_r8+0.25_r8*tc)*1E-3_r8

    ! get mass mean radius
    r_bc = hetraer(1)
    r_dust_a1 = hetraer(2)
    r_dust_a3 = hetraer(3)

    ! calculate collision kernels as a function of environmental parameters and aerosol/droplet sizes
    call collkernel(t, p, eswtr, rhwincloud, r3lx,         &
         r_bc,                                  &  ! BC modes
         r_dust_a1, r_dust_a3,                  &  ! dust modes
         Kcoll_bc,                              &  ! collision kernel [cm3 s-1]
         Kcoll_dust_a1, Kcoll_dust_a3)

    !*****************************************************************************
    !                take water activity into account
    !*****************************************************************************
    ! solute effect
    aw(:) = 1._r8
    molal(:) = 0._r8

    ! The heterogeneous ice freezing temperatures of all IN generally decrease with
    ! increasing total solute mole fraction. Therefore, the large solution concentration
    ! will cause the freezing point depression and the ice freezing temperatures of all
    ! IN will get close to the homogeneous ice freezing temperatures. Since we take into
    ! account water activity for three heterogeneous freezing modes(immersion, deposition,
    ! and contact), we utilize interstitial aerosols(not cloudborne aerosols) to calculate
    ! water activity.
    ! If the index of IN is 0, it means three freezing modes of this aerosol are depressed.

    do icol = 1, 3
       !calculate molality
       if ( total_interstitial_aer_num(icol) > 0._r8 ) then
          molal(icol) = (1.e-6_r8*awcam(icol)*(1._r8-awfacm(icol))/(Mso4*total_interstitial_aer_num(icol)*1.e6_r8))/ &
               (4*pi/3*rhoh2o*(MAX(r3lx,4.e-6_r8))**3)
          aw(icol) = 1._r8/(1._r8 + 2.9244948e-2_r8*molal(icol) + 2.3141243e-3_r8*molal(icol)**2 + 7.8184854e-7_r8*molal(icol)**3)
       end if
    end do

    !*****************************************************************************
    !                immersion freezing begin
    !*****************************************************************************

    frzbcimm = 0._r8
    frzduimm = 0._r8
    frzbccnt = 0._r8
    frzducnt = 0._r8
    frzbcdep = 0._r8
    frzdudep = 0._r8

    ! critical germ size
    rgimm = 2*vwice*sigma_iw/(kboltz*t*LOG(supersatice))

    ! take solute effect into account
    rgimm_bc = rgimm
    rgimm_dust_a1 = rgimm
    rgimm_dust_a3 = rgimm

    do_bc = .false.

    if (aw(id_dst1)*supersatice > 1._r8 ) then
       do_dst1 = .true.
       rgimm_dust_a1 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst1)*supersatice))
    else
       do_dst1 = .false.
    end if

    if (aw(id_dst3)*supersatice > 1._r8 ) then
       do_dst3 = .true.
       rgimm_dust_a3 = 2*vwice*sigma_iw/(kboltz*t*LOG(aw(id_dst3)*supersatice))
    else
       do_dst3 = .false.
    end if

    ! homogeneous energy of germ formation
    dg0imm_bc      = 4*pi/3._r8*sigma_iw*rgimm_bc**2
    dg0imm_dust_a1 = 4*pi/3._r8*sigma_iw*rgimm_dust_a1**2
    dg0imm_dust_a3 = 4*pi/3._r8*sigma_iw*rgimm_dust_a3**2

    ! prefactor
    Aimm_bc = n1*((vwice*rhplanck)/(rgimm_bc**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_bc))
    Aimm_dust_a1 = n1*((vwice*rhplanck)/(rgimm_dust_a1**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_dust_a1))
    Aimm_dust_a3 = n1*((vwice*rhplanck)/(rgimm_dust_a3**3)*SQRT(3._r8/pi*kboltz*T*dg0imm_dust_a3))

    ! nucleation rate per particle
    Jimm_bc = Aimm_bc*r_bc**2/sqroot_f_imm_bc*EXP((-dga_imm_bc-f_imm_bc*dg0imm_bc)/(kboltz*T))
    if (.not. pdf_imm_in) then
       ! 1/sqrt(f)
       ! the expression of Chen et al. (sqrt(f)) may however lead to unphysical
       ! behavior as it implies J->0 when f->0 (i.e. ice nucleation would be
       ! more difficult on easily wettable materials).
       Jimm_dust_a1 = Aimm_dust_a1*r_dust_a1**2/sqroot_f_imm_dust_a1*EXP((-dga_imm_dust-f_imm_dust_a1*dg0imm_dust_a1)/(kboltz*T))
       Jimm_dust_a3 = Aimm_dust_a3*r_dust_a3**2/sqroot_f_imm_dust_a3*EXP((-dga_imm_dust-f_imm_dust_a3*dg0imm_dust_a3)/(kboltz*T))
    end if

    if (pdf_imm_in) then
       dim_Jimm_dust_a1 = 0.0_r8
       dim_Jimm_dust_a3 = 0.0_r8
       do icol = i1,i2
          ! 1/sqrt(f)
          dim_Jimm_dust_a1(icol) = Aimm_dust_a1*r_dust_a1**2/sqroot_dim_f_imm_dust_a1(icol) &
               *EXP((-dga_imm_dust-dim_f_imm_dust_a1(icol)*dg0imm_dust_a1)/(kboltz*T))
          dim_Jimm_dust_a1(icol) = max(dim_Jimm_dust_a1(icol), 0._r8)

          dim_Jimm_dust_a3(icol) = Aimm_dust_a3*r_dust_a3**2/sqroot_dim_f_imm_dust_a1(icol) &
               *EXP((-dga_imm_dust-dim_f_imm_dust_a3(icol)*dg0imm_dust_a3)/(kboltz*T))
          dim_Jimm_dust_a3(icol) = max(dim_Jimm_dust_a3(icol), 0._r8)
       end do
    end if

    ! Limit to 1% of available potential IN (for BC), no limit for dust
    if (pdf_imm_in) then
       sum_imm_dust_a1 = 0._r8
       sum_imm_dust_a3 = 0._r8
       do icol = i1,i2-1
          sum_imm_dust_a1 = sum_imm_dust_a1 &
               + 0.5_r8*((pdf_imm_theta(icol)*exp(-dim_Jimm_dust_a1(icol)*deltat) &
               +          pdf_imm_theta(icol+1)*exp(-dim_Jimm_dust_a1(icol+1)*deltat)))*pdf_d_theta
          sum_imm_dust_a3 = sum_imm_dust_a3 + &
               0.5_r8*((pdf_imm_theta(icol)*exp(-dim_Jimm_dust_a3(icol)*deltat)+ &
                        pdf_imm_theta(icol+1)*exp(-dim_Jimm_dust_a3(icol+1)*deltat)))*pdf_d_theta
       end do
       do icol = i1,i2
          if (sum_imm_dust_a1 > 0.99_r8) then
             sum_imm_dust_a1 = 1.0_r8
          end if
          if (sum_imm_dust_a3 > 0.99_r8) then
             sum_imm_dust_a3 = 1.0_r8
          end if
       end do

    end if

    if (.not.tot_in) then
       if (do_bc) frzbcimm = frzbcimm+MIN(limfacbc*total_cloudborne_aer_num(id_bc)/deltat, &
            total_cloudborne_aer_num(id_bc)/deltat*(1._r8-exp(-Jimm_bc*deltat)))

       if (.not. pdf_imm_in) then
          if (do_dst1) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat, &
               total_cloudborne_aer_num(id_dst1)/deltat*(1._r8-exp(-Jimm_dust_a1*deltat)))
          if (do_dst3) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat, &
               total_cloudborne_aer_num(id_dst3)/deltat*(1._r8-exp(-Jimm_dust_a3*deltat)))
       else
          if (do_dst1) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst1)/deltat,        &
               total_cloudborne_aer_num(id_dst1)/deltat*(1._r8-sum_imm_dust_a1))
          if (do_dst3) frzduimm = frzduimm+MIN(1*total_cloudborne_aer_num(id_dst3)/deltat,        &
               total_cloudborne_aer_num(id_dst3)/deltat*(1._r8-sum_imm_dust_a3))
       end if

    else
       if (do_bc) frzbcimm = frzbcimm+MIN(limfacbc*fn(id_bc)*total_aer_num(id_bc)/deltat, &
            fn(id_bc)*total_aer_num(id_bc)/deltat*(1._r8-exp(-Jimm_bc*deltat)))

       if (.not. pdf_imm_in) then
          if (do_dst1) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat, &
               fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1._r8-exp(-Jimm_dust_a1*deltat)))
          if (do_dst3) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat, &
               fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1._r8-exp(-Jimm_dust_a3*deltat)))
       else
          if (do_dst1) frzduimm = frzduimm+MIN(1*fn(id_dst1)*total_aer_num(id_dst1)/deltat,        &
               fn(id_dst1)*total_aer_num(id_dst1)/deltat*(1._r8-sum_imm_dust_a1))
          if (do_dst3) frzduimm = frzduimm+MIN(1*fn(id_dst3)*total_aer_num(id_dst3)/deltat,        &
               fn(id_dst3)*total_aer_num(id_dst3)/deltat*(1._r8-sum_imm_dust_a3))
       end if
    end if

    if (t > 263.15_r8) then
       frzduimm = 0._r8
       frzbcimm = 0._r8
    end if

    !----------------------------------
    !   Deposition nucleation
    !----------------------------------

    ! critical germ size
    ! assume 98% RH in mixed-phase clouds (Korolev & Isaac, JAS 2006)
    rgdep=2*vwice*sigma_iv/(kboltz*t*LOG(rhwincloud*supersatice))

    ! homogeneous energy of germ formation
    dg0dep = 4*pi/3._r8*sigma_iv*rgdep**2

    ! prefactor
    ! attention: division of small numbers
    Adep = (rhwincloud*eswtr)**2*(vwice/(mwh2o*amu))/(kboltz*T*nus)*SQRT(sigma_iv/(kboltz*T))

    ! nucleation rate per particle
    if (rgdep > 0) then
       Jdep_bc = Adep*r_bc**2/sqroot_f_dep_bc*EXP((-dga_dep_bc-f_dep_bc*dg0dep)/(kboltz*T))
       Jdep_dust_a1 = Adep*r_dust_a1**2/sqroot_f_dep_dust_a1*EXP((-dga_dep_dust-f_dep_dust_a1*dg0dep)/(kboltz*T))
       Jdep_dust_a3 = Adep*r_dust_a3**2/sqroot_f_dep_dust_a3*EXP((-dga_dep_dust-f_dep_dust_a3*dg0dep)/(kboltz*T))
    else
       Jdep_bc = 0._r8
       Jdep_dust_a1 = 0._r8
       Jdep_dust_a3 = 0._r8
    end if

    ! Limit to 1% of available potential IN (for BC), no limit for dust
    if (.not.tot_in) then
       if (do_bc) frzbcdep = frzbcdep+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
            uncoated_aer_num(id_bc)/deltat &
            *(1._r8-exp(-Jdep_bc*deltat)))
       if (do_dst1) frzdudep = frzdudep+MIN(uncoated_aer_num(id_dst1)/deltat, &
            uncoated_aer_num(id_dst1)/deltat &
            *(1._r8-exp(-Jdep_dust_a1*deltat)))
       if (do_dst3) frzdudep = frzdudep+MIN(uncoated_aer_num(id_dst3)/deltat, &
            uncoated_aer_num(id_dst3)/deltat &
            *(1._r8-exp(-Jdep_dust_a3*deltat)))
    else
       if (do_bc) frzbcdep = frzbcdep+MIN(limfacbc*(1._r8-fn(id_bc)) &
            *(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat, &
            (1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat &
            *(1._r8-exp(-Jdep_bc*deltat)))
       if (do_dst1) frzdudep = frzdudep+MIN((1._r8-fn(id_dst1)) &
            *(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
            (1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat &
            *(1._r8-exp(-Jdep_dust_a1*deltat)))
       if (do_dst3) frzdudep = frzdudep+MIN((1._r8-fn(id_dst3)) &
            *(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
            (1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat &
            *(1._r8-exp(-Jdep_dust_a3*deltat)))
    end if

    ! ---------------------------
    ! contact nucleation
    ! ---------------------------

    ! homogeneous energy of germ formation
    dg0cnt = 4*pi/3._r8*sigma_iv*rgimm**2

    ! prefactor
    ! attention: division of small numbers
    Acnt = rhwincloud*eswtr*4*pi/(nus*SQRT(2*pi*mwh2o*amu*kboltz*T))

    ! nucleation rate per particle
    Jcnt_bc      = Acnt*r_bc**2     *EXP((-dga_dep_bc-f_cnt_bc*dg0cnt)       /(kboltz*T))*Kcoll_bc*icnlx
    Jcnt_dust_a1 = Acnt*r_dust_a1**2*EXP((-dga_dep_dust-f_cnt_dust_a1*dg0cnt)/(kboltz*T))*Kcoll_dust_a1*icnlx
    Jcnt_dust_a3 = Acnt*r_dust_a3**2*EXP((-dga_dep_dust-f_cnt_dust_a3*dg0cnt)/(kboltz*T))*Kcoll_dust_a3*icnlx

    ! Limit to 1% of available potential IN (for BC), no limit for dust
    if (.not.tot_in) then
       if (do_bc) frzbccnt = frzbccnt+MIN(limfacbc*uncoated_aer_num(id_bc)/deltat, &
            uncoated_aer_num(id_bc)/deltat &
            *(1._r8-exp(-Jcnt_bc*deltat)))
       if (do_dst1) frzducnt = frzducnt+MIN(uncoated_aer_num(id_dst1)/deltat, &
            uncoated_aer_num(id_dst1)/deltat &
            *(1._r8-exp(-Jcnt_dust_a1*deltat)))
       if (do_dst3) frzducnt = frzducnt+MIN(uncoated_aer_num(id_dst3)/deltat, &
            uncoated_aer_num(id_dst3)/deltat &
            *(1._r8-exp(-Jcnt_dust_a3*deltat)))
    else
       if (do_bc) frzbccnt = frzbccnt+MIN(limfacbc*(1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat, &
            (1._r8-fn(id_bc))*(1._r8-dstcoat(1))*total_aer_num(id_bc)/deltat &
            *(1._r8-exp(-Jcnt_bc*deltat)))
       if (do_dst1) frzducnt = frzducnt+MIN((1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat, &
            (1._r8-fn(id_dst1))*(1._r8-dstcoat(2))*total_aer_num(id_dst1)/deltat &
            *(1._r8-exp(-Jcnt_dust_a1*deltat)))
       if (do_dst3) frzducnt = frzducnt+MIN((1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat, &
            (1._r8-fn(id_dst3))*(1._r8-dstcoat(3))*total_aer_num(id_dst3)/deltat &
            *(1._r8-exp(-Jcnt_dust_a3*deltat)))
    end if

    errstring = ' '
    if (frzducnt <= -1._r8) then
       write(iulog,*) 'hetfrz_classnuc_calc: frzducnt', frzducnt, Jcnt_dust_a1,Jcnt_dust_a3, &
            Kcoll_dust_a1, Kcoll_dust_a3
       errstring = 'ERROR in hetfrz_classnuc_calc::frzducnt'
       return
    end if

  end subroutine  hetfrz_classnuc_calc

  !===================================================================================================

  subroutine collkernel( &
       t, pres, eswtr, rhwincloud, r3lx,       &
       r_bc,                                   &  ! BC modes
       r_dust_a1, r_dust_a3,                   &  ! dust modes
       Kcoll_bc,                               &  ! collision kernel [cm3 s-1]
       Kcoll_dust_a1, Kcoll_dust_a3)

    !-----------------------------------------------------------------------
    ! Purpose: calculate collision kernels as a function of
    ! environmental parameters and aerosol/droplet sizes
    ! Author: Corinna Hoose, UiO, October 2009
    ! Modifications: Yong Wang and Xiaohong Liu, UWyo, 12/2012
    !-----------------------------------------------------------------------

    real(r8), intent(in) :: t                ! temperature [K]
    real(r8), intent(in) :: pres             ! pressure [Pa]
    real(r8), intent(in) :: eswtr            ! saturation vapor pressure of water [Pa]
    real(r8), intent(in) :: r3lx             ! volume mean drop radius [m]
    real(r8), intent(in) :: rhwincloud       ! in-cloud relative humidity over water [ ]
    real(r8), intent(in) :: r_bc             ! model radii of BC modes [m]
    real(r8), intent(in) :: r_dust_a1        ! model radii of dust modes [m]
    real(r8), intent(in) :: r_dust_a3        ! model radii of dust modes [m]
    real(r8), intent(out) :: Kcoll_bc        ! collision kernel [cm3 s-1]
    real(r8), intent(out) :: Kcoll_dust_a1
    real(r8), intent(out) :: Kcoll_dust_a3

    ! local variables
    real(r8) :: a, b, c, a_f, b_f, c_f, f
    real(r8) :: tc          ! temperature [deg C]
    real(r8) :: rho_air     ! air density [kg m-3]
    real(r8) :: viscos_air  ! dynamic viscosity of air [kg m-1 s-1]
    real(r8) :: Ktherm_air  ! thermal conductivity of air [J/(m s K)]
    real(r8) :: lambda      ! mean free path [m]
    real(r8) :: Kn          ! Knudsen number [ ]
    real(r8) :: Re          ! Reynolds number [ ]
    real(r8) :: Pr          ! Prandtl number [ ]
    real(r8) :: Sc          ! Schmidt number [ ]
    real(r8) :: vterm       ! terminal velocity [m s-1]
    real(r8) :: Ktherm      ! thermal conductivity of aerosol [J/(m s K)]
    real(r8) :: Dvap        ! water vapor diffusivity [m2 s-1]
    real(r8) :: Daer        ! aerosol diffusivity [m2 s-1]
    real(r8) :: latvap      ! latent heat of vaporization [J kg-1]
    real(r8) :: kboltz      ! Boltzmann constant [J K-1]
    real(r8) :: G           ! thermodynamic function in Cotton et al. [kg m-1 s-1]
    real(r8) :: r_a         ! aerosol radius [m]
    real(r8) :: f_t         ! factor by Waldmann & Schmidt [ ]
    real(r8) :: Q_heat      ! heat flux [J m-2 s-1]
    real(r8) :: Tdiff_cotton ! temperature difference between droplet and environment [K]
    real(r8) :: K_brownian,K_thermo_cotton,K_diffusio_cotton   ! collision kernels [m3 s-1]
    real(r8) :: K_total     ! total collision kernel [cm3 s-1]
    integer  :: imode
    !------------------------------------------------------------------------------------------------

    Kcoll_bc      = 0._r8
    Kcoll_dust_a1 = 0._r8
    Kcoll_dust_a3 = 0._r8

    tc     = t - tmelt
    kboltz = 1.38065e-23_r8

    ! air viscosity for tc<0, from depvel_part.F90
    viscos_air = (1.718_r8+0.0049_r8*tc-1.2e-5_r8*tc*tc)*1.e-5_r8

    ! air density
    rho_air = pres/(rair*t)

    ! mean free path: Seinfeld & Pandis 8.6
    lambda = 2*viscos_air/(pres*SQRT(8/(pi*rair*t)))

    ! latent heat of vaporization, varies with T
    latvap = 1000*(-0.0000614342_r8*tc**3 + 0.00158927_r8*tc**2 - 2.36418_r8*tc + 2500.79_r8)

    ! droplet terminal velocity after Chen & Liu, QJRMS 2004
    a = 8.8462e2_r8
    b = 9.7593e7_r8
    c = -3.4249e-11_r8
    a_f = 3.1250e-1_r8
    b_f = 1.0552e-3_r8
    c_f = -2.4023_r8
    f = EXP(EXP(a_f + b_f*(LOG(r3lx))**3 + c_f*rho_air**1.5_r8))
    vterm = (a+ (b + c*r3lx)*r3lx)*r3lx*f

    ! Reynolds number
    Re = 2*vterm*r3lx*rho_air/viscos_air

    ! thermal conductivity of air: Seinfeld & Pandis eq. 15.75
    Ktherm_air = 1.e-3_r8*(4.39_r8+0.071_r8*t)  !J/(m s K)

    ! Prandtl number
    Pr = viscos_air*cpair/Ktherm_air

    ! water vapor diffusivity: Pruppacher & Klett 13-3
    Dvap = 0.211e-4_r8*(t/273.15_r8)*(101325._r8/pres)

    ! G-factor = rhoh2o*Xi in Rogers & Yau, p. 104
    G = rhoh2o/((latvap/(rh2o*t) - 1)*latvap*rhoh2o/(Ktherm_air*t) + rhoh2o*rh2o*t/(Dvap*eswtr))

    ! variables depending on aerosol radius
    ! loop over 3 aerosol modes
    do imode = 1, 3
       if (imode == 1) r_a = r_bc
       if (imode == 2) r_a = r_dust_a1
       if (imode == 3) r_a = r_dust_a3
       ! Knudsen number (Seinfeld & Pandis 8.1)
       Kn = lambda/r_a
       ! aerosol diffusivity
       Daer = kboltz*t*(1 + Kn)/(6*pi*r_a*viscos_air)
       ! Schmidt number
       Sc = viscos_air/(Daer*rho_air)

       ! Young (1974) first equ. on page 771
       K_brownian = 4*pi*r3lx*Daer*(1 + 0.3_r8*Re**0.5_r8*Sc**0.33_r8)

       ! thermal conductivities from Seinfeld & Pandis, Table 8.6
       if (imode == 1) Ktherm = 4.2_r8 ! Carbon
       if (imode == 2 .or. imode == 3) Ktherm = 0.72_r8 ! clay

       ! form factor
       f_t = 0.4_r8*(1._r8 + 1.45_r8*Kn + 0.4_r8*Kn*EXP(-1._r8/Kn))      &
            *(Ktherm_air + 2.5_r8*Kn*Ktherm)                      &
            /((1._r8 + 3._r8*Kn)*(2._r8*Ktherm_air + 5._r8*Kn*Ktherm+Ktherm))

       ! calculate T-Tc as in Cotton et al.
       Tdiff_cotton = -G*(rhwincloud - 1._r8)*latvap/Ktherm_air
       Q_heat = Ktherm_air/r3lx*(1._r8 + 0.3_r8*Re**0.5_r8*Pr**0.33_r8)*Tdiff_cotton
       K_thermo_cotton = 4._r8*pi*r3lx*r3lx*f_t*Q_heat/pres
       K_diffusio_cotton = -(1._r8/f_t)*(rh2o*t/latvap)*K_thermo_cotton
       K_total = 1.e6_r8*(K_brownian + K_thermo_cotton + K_diffusio_cotton)  ! convert m3/s -> cm3/s

       ! set K to 0 if negative
       if (K_total < 0._r8) K_total = 0._r8

       if (imode == 1) Kcoll_bc = K_total
       if (imode == 2) Kcoll_dust_a1 = K_total
       if (imode == 3) Kcoll_dust_a3 = K_total
    end do

  end subroutine collkernel

  !===================================================================================================

  subroutine hetfrz_classnuc_init_pdftheta()

    ! Local variables:
    real(r8) :: theta_min, theta_max
    real(r8) :: x1_imm, x2_imm
    real(r8) :: norm_theta_imm
    real(r8) :: imm_dust_mean_theta
    real(r8) :: imm_dust_var_theta
    integer  :: index
    real(r8) :: cos_angle
    real(r8) :: temp
    !----------------------------------------------------------------------------

    theta_min           = pi/180._r8
    theta_max           = 179._r8/180._r8*pi
    imm_dust_mean_theta = 46.0_r8/180.0_r8*pi
    imm_dust_var_theta  = 0.01_r8

    pdf_d_theta = (179._r8-1._r8)/180._r8*pi/(pdf_n_theta-1)

    x1_imm = (LOG(theta_min) - LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
    x2_imm = (LOG(theta_max) - LOG(imm_dust_mean_theta))/(sqrt(2.0_r8)*imm_dust_var_theta)
    norm_theta_imm = (ERF(x2_imm) - ERF(x1_imm))*0.5_r8
    dim_theta = 0.0_r8
    pdf_imm_theta  = 0.0_r8
    do index = i1, i2
       dim_theta(index) = 1._r8/180._r8*pi + (index-1)*pdf_d_theta
       pdf_imm_theta(index) = exp(-((LOG(dim_theta(index)) - LOG(imm_dust_mean_theta))**2._r8) / &
            (2._r8*imm_dust_var_theta**2._r8) ) /  &
            (dim_theta(index)*imm_dust_var_theta*SQRT(2*pi))/norm_theta_imm
    end do

    do index = i1, i2
       cos_angle = cos(dim_theta(index))
       temp = (2+cos_angle)*(1-cos_angle)**2/4._r8
       dim_f_imm_dust_a1(index) = temp
       dim_f_imm_dust_a3(index) = temp
    end do

  end subroutine hetfrz_classnuc_init_pdftheta

end module oslo_aero_hetfrz
