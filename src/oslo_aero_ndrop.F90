module oslo_aero_ndrop

  !---------------------------------------------------------------------------------
  !  Droplet activation by oslo modal aerosols
  !  Compute vertical diffusion and nucleation of cloud droplets
  !---------------------------------------------------------------------------------

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use spmd_utils,        only: masterproc
  use ppgrid,            only: pcols, pver, pverp
  use physconst,         only: pi, rhoh2o, mwh2o, r_universal, rh2o
  use physconst,         only: gravit, latvap, cpair, rair
  use constituents,      only: pcnst, cnst_get_ind, cnst_name, cnst_spec_class_gas, cnst_species_class
  use physics_types,     only: physics_state, physics_ptend, physics_ptend_init
  use physics_buffer,    only: physics_buffer_desc, pbuf_get_index, pbuf_get_field
  use wv_saturation,     only: qsat
  use phys_control,      only: phys_getopts, use_hetfrz_classnuc
  use ref_pres,          only: top_lev => trop_cloud_top_lev
  use shr_spfn_mod,      only: erf => shr_spfn_erf
  use cam_history,       only: addfld, add_default, horiz_only, fieldname_len, outfld
  use cam_abortutils,    only: endrun
  use cam_logfile,       only: iulog
  use perf_mod,          only: t_startf, t_stopf
  !
  use oslo_aero_share,   only: calculateNumberMedianRadius
  use oslo_aero_share,   only: getNumberOfTracersInMode, getNumberOfAerosolTracers
  use oslo_aero_share,   only: getCloudTracerName, getTracerIndex
  use oslo_aero_share,   only: fillAerosolTracerList, fillInverseAerosolTracerList
  use oslo_aero_share,   only: nmodes, nbmodes, max_tracers_per_mode
  use oslo_aero_share,   only: smallNumber, CloudTracerIndex
  use oslo_aero_share,   only: l_so4_a1, l_so4_ac, l_so4_a2, l_bc_ac, l_om_ac, l_soa_a1

  implicit none
  private

  ! public routines
  public :: ndrop_init_oslo
  public :: dropmixnuc_oslo

  ! private routines
  private :: maxsat_oslo
  private :: ccncalc_oslo
  private :: activate_modal_oslo

  ! private variables
  real(r8) :: t0            ! reference temperature
  real(r8) :: aten
  real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
  real(r8) :: alog2, alog3, alogaten
  real(r8) :: third, twothird, sixth, zero
  real(r8) :: sq2, sqpi

  ! supersaturation (%) to determine ccn concentration
  integer,  parameter :: psat=7    ! number of supersaturations to calc ccn concentration
  real(r8), parameter :: supersat(psat)= (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.15_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)
  character(len=8) :: ccn_name(psat)= (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6','CCN7'/)

  ! indices in state and pbuf structures
  integer :: numliq_idx = -1
  integer :: kvh_idx    = -1

  ! description of modal aerosols
  integer               :: ntot_amode     ! number of aerosol modes
  integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
  real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
  real(r8), allocatable :: dgnumlo_amode(:)
  real(r8), allocatable :: dgnumhi_amode(:)
  real(r8), allocatable :: voltonumblo_amode(:)
  real(r8), allocatable :: voltonumbhi_amode(:)

  logical :: history_aerosol      ! Output the MAM aerosol tendencies
  character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
  character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

  ! local indexing for MAM
  integer, allocatable :: mam_idx(:,:)        ! table for local indexing of modal aero number and mmr
  integer :: ncnst_tot                        ! total number of mode number conc + mode species

  integer  :: tracer_index(0:nmodes,max_tracers_per_mode)  ! tracer index
  real(r8) :: sumFraction2(pcnst,pver)

  logical :: tendencyCounted(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
  integer :: n_aerosol_tracers
  integer :: aerosolTracerList(pcnst)         !List where indexes 1...n_aerosol_tracers are the indexes in pcnst
                                              !..something like (/ l_so4_a1, l_bc_a, .../)etc
  integer :: inverseAerosolTracerList(pcnst)  !List where you can back the place in aerosolTracerList if you know the
                                              !tracer index. So in the example above inverseAerosolTracerList(l_so4_a1) = 1

  ! ptr2d_t is used to create arrays of pointers to 2D fields
  type ptr2d_t
     real(r8), pointer :: fld(:,:)
  end type ptr2d_t

  ! modal aerosols
  logical :: lq(pcnst) = .false.      ! set flags true for constituents with non-zero tendencies

!===============================================================================
contains
!===============================================================================

  subroutine ndrop_init_oslo()

    integer            :: ii, itrac, ispec, imode, mm, ilev, isat
    integer            :: nspec_max    ! max number of species in a mode
    character(len=32)  :: tmpname
    character(len=32)  :: tmpname_cw
    character(len=128) :: long_name
    character(len=8)   :: unit
    character(len=10)  :: modeString
    character(len=20)  :: varname
    !-------------------------------------------------------------------------------

    ! get indices into state%q and pbuf structures
    call cnst_get_ind('NUMLIQ', numliq_idx)

    kvh_idx = pbuf_get_index('kvh')

    zero     = 0._r8
    third    = 1._r8/3._r8
    twothird = 2._r8*third
    sixth    = 1._r8/6._r8
    sq2      = sqrt(2._r8)
    sqpi     = sqrt(pi)

    t0       = 273._r8
    surften  = 0.076_r8
    aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
    alogaten = log(aten)
    alog2    = log(2._r8)
    alog3    = log(3._r8)

    ! get info about the modal aerosols, get ntot_amode
    ntot_amode = nmodes
    allocate( &
         nspec_amode(ntot_amode),  &
         sigmag_amode(ntot_amode), &
         dgnumlo_amode(ntot_amode), &
         dgnumhi_amode(ntot_amode), &
         voltonumblo_amode(ntot_amode), &
         voltonumbhi_amode(ntot_amode)  )

    do imode = 1,ntot_amode
       nspec_amode(imode) = getNumberOfTracersInMode(imode)
    enddo

    do imode = 1,ntot_amode
      do ispec = 1,nspec_amode(imode)
        tracer_index(imode,ispec) = getTracerIndex(imode,ispec,.false.)
      end do
    end do

    sumFraction2(:,:) = 0.0_r8
    do ilev=top_lev, pver
      do itrac=1,pcnst
        do imode=1,ntot_amode
          do ispec=1,nspec_amode(imode)
            if (tracer_index(imode,ispec) == itrac) then
              sumFraction2(itrac,ilev) = sumFraction2(itrac,ilev) + 1.0_r8
            endif
          end do ! tracers in mode
        end do ! mode
      end do
    end do

    ! Init the table for local indexing of mam number conc
    ! This table uses species index 0 for the number conc.

    ! Find max number of species in all the modes, and the total
    ! number of mode number concentrations + mode species
    nspec_max = nspec_amode(1)
    ncnst_tot = nspec_amode(1) + 1
    do imode = 2, ntot_amode
       nspec_max = max(nspec_max, nspec_amode(imode))
       ncnst_tot = ncnst_tot + nspec_amode(imode) + 1
    end do

    allocate(mam_idx(ntot_amode,0:nspec_max))
    allocate(fieldname(ncnst_tot))
    allocate(fieldname_cw(ncnst_tot))

    ! Local indexing compresses the mode and number/mass indicies into one index.
    ! This indexing is used by the pointer arrays used to reference state and pbuf
    ! fields.
    ii = 0
    do imode = 1, ntot_amode
       do ispec = 0, nspec_amode(imode)
          ii = ii + 1
          mam_idx(imode,ispec) = ii
       end do
    end do

    ! Add dropmixnuc tendencies for all modal aerosol species

    call phys_getopts(history_aerosol_out=history_aerosol)

    n_aerosol_tracers = getNumberOfAerosolTracers()
    call fillAerosolTracerList(aerosolTracerList)
    call fillInverseAerosolTracerList(aerosolTracerList, inverseAerosolTracerList, n_aerosol_tracers)
    if (masterproc) then
       do ii=1,n_aerosol_tracers
          write(iulog,*) "aerosolTracerList", ii, aerosolTracerList(ii), inverseAerosolTracerList(aerosolTracerList(ii))
       end do
    end if

    lq(:) = .false.  !Initialize

    !Set up tendencies for tracers (output)
    do imode=1,ntot_amode
       do ispec=1,nspec_amode(imode)
          itrac = tracer_index(imode,ispec)

          if(.NOT. lq(itrac))then
             !add dropmixnuc tendencies
             mm=mam_idx(imode,ispec)
             fieldname(mm)=trim(cnst_name(itrac))//"_mixnuc1"
             fieldname_cw(mm)=trim(getCloudTracerName(itrac))//"_mixnuc1"

             long_name = trim(fieldname(mm)) // ' dropmixnuc column tendency'
             call addfld(trim(fieldname(mm)), horiz_only ,'A', "kg/m2/s",long_name)

             long_name = trim(fieldname_cw(mm)) // ' dropmixnuc column tendency'
             call addfld(trim(fieldname_cw(mm)), horiz_only, 'A', "kg/m2/s",long_name)

             if (history_aerosol) then
                call add_default(trim(fieldname(mm)), 1, ' ')
                call add_default(trim(fieldname_cw(mm)),1,' ')
             endif

             !Do tendencies of this tracer
             lq(itrac)=.TRUE.
          endif
       enddo
    enddo
    do imode=1,ntot_amode
       modeString="  "
       write(modeString,"(I2)"),imode
       if(imode < 10) modeString="0"//adjustl(modeString)
       varName = "NMR"//trim(modeString)
       call addfld(varName, (/ 'lev' /),'A', 'm  ', 'number median radius mode '//modeString)
       if(history_aerosol)call add_default(varName, 1, ' ')

       varName = "NCONC"//trim(modeString)
       call addfld(varName, (/ 'lev' /),'A', '#/m3  ', 'number concentration mode '//modeString)
       if(history_aerosol)call add_default(varName, 1, ' ')

       varName = "VCONC"//trim(modeString)
       call addfld(varName, (/ 'lev' /),'A',  'm3/m3  ','volume concentration mode '//modeString)
       if(history_aerosol)call add_default(varName, 1, ' ')

       varName = "SIGMA"//trim(modeString)
       call addfld(varName, (/ 'lev' /),'A', '-','Std. dev. mode '//modeString)

       if(history_aerosol)call add_default(varName, 1, ' ')
       varName = "HYGRO"//trim(modeString)
       call addfld(varName, (/ 'lev' /),'A','-','Hygroscopicity '//modeString)
       if(history_aerosol)call add_default(varName, 1, ' ')
    end do
    call addfld('CCN1',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.02%')
    call addfld('CCN2',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.05%')
    call addfld('CCN3',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.1%')
    call addfld('CCN4',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.15%')
    call addfld('CCN5',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.2%')
    call addfld('CCN6',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.5%')
    call addfld('CCN7',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=1.0%')

    if(history_aerosol)then
       do isat = 1, psat
          call add_default(ccn_name(isat), 1, ' ')
       enddo
    end if

    call addfld('WTKE',     (/ 'lev' /), 'A', 'm/s', 'Standard deviation of updraft velocity')
    call addfld('NDROPMIX', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number mixing')
    call addfld('NDROPSRC', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number source')
    call addfld('NDROPSNK', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number loss by microphysics')
    call addfld('NDROPCOL', horiz_only,  'A', '#/m2', 'Column droplet number')

  end subroutine ndrop_init_oslo

  !===============================================================================

  subroutine dropmixnuc_oslo( state, ptend, dtmicro, pbuf, wsub,                    &
       cldn, cldo, cldliqf, hasAerosol, CProcessModes, f_c, f_bc, f_aq, f_so4_cond, &
       f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam,                        &
       numberConcentration, volumeConcentration,  hygroscopicity, lnsigma, tendnd, fn_in)

    ! vertical diffusion and nucleation of cloud droplets
    ! assume cloud presence controlled by cloud fraction
    ! doesn't distinguish between warm, cold clouds

    ! arguments
    type(physics_state), target, intent(in)  :: state
    type(physics_ptend),         intent(out) :: ptend
    real(r8),                    intent(in)  :: dtmicro                                  ! time step for microphysics (s)
    type(physics_buffer_desc),   pointer     :: pbuf(:)
    real(r8),                    intent(in)  :: wsub(pcols,pver)                         ! subgrid vertical velocity
    real(r8),                    intent(in)  :: cldn(pcols,pver)                         ! cloud fraction
    real(r8),                    intent(in)  :: cldo(pcols,pver)                         ! cloud fraction on previous time step
    real(r8),                    intent(in)  :: cldliqf(pcols,pver)                      ! liquid cloud fraction (liquid / (liquid + ice))
    logical ,                    intent(in)  :: hasAerosol(pcols, pver, nmodes)
    real(r8),                    intent(in)  :: CProcessModes(pcols,pver)
    real(r8),                    intent(in)  :: f_c(pcols,pver)
    real(r8),                    intent(in)  :: f_bc(pcols,pver)
    real(r8),                    intent(in)  :: f_aq(pcols,pver)
    real(r8),                    intent(in)  :: f_so4_cond(pcols,pver)
    real(r8),                    intent(in)  :: f_soa(pcols,pver)
    real(r8),                    intent(in)  :: cam(pcols,pver,nbmodes)
    real(r8),                    intent(in)  :: f_acm(pcols,pver, nbmodes)
    real(r8),                    intent(in)  :: f_bcm(pcols,pver, nbmodes)
    real(r8),                    intent(in)  :: f_aqm(pcols, pver, nbmodes)
    real(r8),                    intent(in)  :: f_so4_condm(pcols, pver, nbmodes)        !Needed in "get component fraction
    real(r8),                    intent(in)  :: f_soam(pcols,pver,nbmodes)
    real(r8),                    intent(in)  :: numberConcentration(pcols,pver,0:nmodes) ![#/m3] number concentraiton
    real(r8),                    intent(in)  :: volumeConcentration(pcols,pver,nmodes)   ![m3/m3] volume concentration
    real(r8),                    intent(in)  :: hygroscopicity(pcols,pver,nmodes)        ![-] hygroscopicity
    real(r8),                    intent(in)  :: lnsigma(pcols,pver,nmodes)               ![-] log(base e) sigma
    real(r8),                    intent(out) :: tendnd(pcols,pver)                       ! change in droplet number concentration (#/kg/s)

    ! Local variables
    integer  :: lchnk                           ! chunk identifier
    integer  :: ncol                            ! number of columns
    real(r8), pointer :: ncldwtr(:,:)           ! droplet number concentration (#/kg)
    real(r8), pointer :: temp(:,:)              ! temperature (K)
    real(r8), pointer :: omega(:,:)             ! vertical velocity (Pa/s)
    real(r8), pointer :: pmid(:,:)              ! mid-level pressure (Pa)
    real(r8), pointer :: pint(:,:)              ! pressure at layer interfaces (Pa)
    real(r8), pointer :: pdel(:,:)              ! pressure thickess of layer (Pa)
    real(r8), pointer :: rpdel(:,:)             ! inverse of pressure thickess of layer (/Pa)
    real(r8), pointer :: zm(:,:)                ! geopotential height of level (m)
    real(r8), pointer :: kvh(:,:)               ! vertical diffusivity (m2/s)
    type(ptr2d_t), allocatable :: raer(:)       ! aerosol mass, number mixing ratios
    type(ptr2d_t), allocatable :: qqcw(:)
    real(r8) :: raertend(pver)                  ! tendency of aerosol mass, number mixing ratios
    real(r8) :: qqcwtend(pver)                  ! tendency of cloudborne aerosol mass, number mixing ratios

    real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8
    real(r8), parameter :: wmixmin = 0.1_r8     ! minimum turbulence vertical velocity (m/s)
    real(r8) :: sq2pi

    integer  :: icol, ilev, ispec, imode, mm, isubmix, mode_cnt, itrac, itrac2, isat
    integer  :: km1, kp1
    integer  :: nnew, nsav, ntemp
    integer  :: nsubmix, nsubmix_bnd
    integer, save :: count_submix(100)
    integer  :: phase                           ! phase of aerosol

    real(r8) :: arg
    real(r8) :: dtinv
    real(r8) :: dtmin, tinv, dtt
    real(r8) :: lcldn(pcols,pver)
    real(r8) :: lcldo(pcols,pver)

    real(r8) :: zs(pver)                        ! inverse of distance between levels (m)
    real(r8) :: qcld(pver)                      ! cloud droplet number mixing ratio (#/kg)
    real(r8) :: qncld(pver)                     ! droplet number nucleated on cloud boundaries
    real(r8) :: srcn(pver)                      ! droplet source rate (/s)
    real(r8) :: cs(pcols,pver)                  ! air density (kg/m3)
    real(r8) :: csbot(pver)                     ! air density at bottom (interface) of layer (kg/m3)
    real(r8) :: csbot_cscen(pver)               ! csbot(icol)/cs(icol,ilev)
    real(r8) :: dz(pcols,pver)                  ! geometric thickness of layers (m)

    real(r8) :: wtke(pcols,pver)                ! turbulent vertical velocity at base of layer k (m/s)
    real(r8) :: wtke_cen(pcols,pver)            ! turbulent vertical velocity at center of layer k (m/s)
    real(r8) :: wbar, wmix, wmin, wmax

    real(r8) :: zn(pver)                        ! g/pdel (m2/g) for layer
    real(r8) :: flxconv                         ! convergence of flux into lowest layer

    real(r8) :: wdiab                           ! diabatic vertical velocity
    real(r8) :: ekd(pver)                       ! diffusivity for droplets (m2/s)
    real(r8) :: ekk(0:pver)                     ! density*diffusivity for droplets (kg/m3 m2/s)
    real(r8) :: ekkp(pver)                      ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
    real(r8) :: ekkm(pver)                      ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
                                                ! above layer ilev  (ilev,ilev+1 interface)
    real(r8) :: dum, dumc
    real(r8) :: tmpa
    real(r8) :: dact
    real(r8) :: fluxntot                        ! (#/cm2/s)
    real(r8) :: dtmix
    real(r8) :: alogarg
    real(r8) :: overlapp(pver), overlapm(pver)  ! cloud overlap below, cloud overlap above
    real(r8) :: nsource(pcols,pver)             ! droplet number source (#/kg/s)
    real(r8) :: ndropmix(pcols,pver)            ! droplet number mixing (#/kg/s) (diagnostic)
    real(r8) :: ndropcol(pcols)                 ! column droplet number (#/m2) (diagnostic)
    real(r8) :: cldo_tmp, cldn_tmp
    real(r8) :: tau_cld_regenerate
    real(r8) :: tau_cld_regenerate_exp
    real(r8) :: zeroaer(pver)
    real(r8) :: taumix_internal_pver_inv        ! 1/(internal mixing time scale for ilev=pver) (1/s)

    real(r8), allocatable :: nact(:,:)          ! fractional aero. number  activation rate (/s)
    real(r8), allocatable :: mact(:,:)          ! fractional aero. mass    activation rate (/s)

    real(r8), allocatable :: raercol(:,:,:)     ! single column of aerosol mass, number mixing ratios
    real(r8), allocatable :: raercol_cw(:,:,:)  ! same as raercol but for cloud-borne phase

                                                !to avoid excessive calls to boundary layer scheme
    real(r8), allocatable :: raercol_tracer(:,:,:)
    real(r8), allocatable :: raercol_cw_tracer(:,:,:)
    real(r8), allocatable :: mact_tracer(:,:)
    real(r8), allocatable :: mfullact_tracer(:,:)

    real(r8)              :: na(pcols), va(pcols), hy(pcols)
    real(r8), allocatable :: naermod(:)         ! (1/m3)
    real(r8), allocatable :: hygro(:)           ! hygroscopicity of aerosol mode
    real(r8), allocatable :: vaerosol(:)        ! interstit+activated aerosol volume conc (cm3/cm3)

    real(r8)              :: source(pver)

    real(r8), allocatable :: fn(:)              ! activation fraction for aerosol number
    real(r8), intent(out) :: fn_in(pcols,pver,0:nmodes)
    real(r8), allocatable :: fm(:)              ! activation fraction for aerosol mass

    real(r8), allocatable :: fluxn(:)           ! number  activation fraction flux (cm/s)
    real(r8), allocatable :: fluxm(:)           ! mass    activation fraction flux (cm/s)
    real(r8)              :: flux_fullact(pver) ! 100%    activation fraction flux (cm/s)
    ! note:  activation fraction fluxes are defined as
    ! fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
    !       / [aero. number conc. in updraft, just below cloudbase (#/cm3)]

    real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
    real(r8), allocatable :: coltend_cw(:,:)    ! column tendency
    real(r8)              :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat

    !for gas species turbulent mixing
    real(r8), pointer     :: rgas(:, :, :)
    real(r8), allocatable :: rgascol(:, :, :)
    real(r8), allocatable :: coltendgas(:)
    real(r8)              :: zerogas(pver)
    character*200         :: fieldnamegas

    real(r8)              :: numberMedianRadius(pcols,pver,nmodes)
    real(r8)              :: sigma(pcols,pver,nmodes)                 ![-] sigma
    real(r8)              :: constituentFraction
    real(r8)              :: volumeCore(pcols,pver,nmodes)
    real(r8)              :: volumeCoat(pcols,pver,nmodes)
    integer               :: tracerIndex
    integer               :: cloud_tracer_index
    integer               :: kcomp
    integer               :: speciesMap(nmodes)
    real(r8), allocatable :: fn_tmp(:), fm_tmp(:)
    real(r8), allocatable :: fluxn_tmp(:), fluxm_tmp(:)
    real(r8)              :: componentFraction
    real(r8)              :: componentFractionOK(nmodes,pcnst,pver)
    real(r8)              :: sumFraction(pcnst,pver)
    logical               :: alert
    real(r8), dimension(pver, pcnst) :: massBalance
    real(r8), dimension(pver, pcnst) :: newMass
    real(r8), dimension(pver,pcnst)  :: newCloud, oldCloud, newAerosol, oldAerosol, deltaCloud
    real(r8), allocatable            :: lnsigman(:)
    character(len=2)                 :: modeString
    character(len=20)                :: varname
    integer                          :: numberOfModes
    !-------------------------------------------------------------------------------

    sq2pi = sqrt(2._r8*pi)

    lchnk = state%lchnk
    ncol  = state%ncol

    ncldwtr  => state%q(:,:,numliq_idx)
    temp     => state%t
    omega    => state%omega
    pmid     => state%pmid
    pint     => state%pint
    pdel     => state%pdel
    rpdel    => state%rpdel
    zm       => state%zm

    call pbuf_get_field(pbuf, kvh_idx, kvh)

    ! Create the liquid weighted cloud fractions that were passsed in
    ! before. This doesn't seem like the best variable, since the cloud could
    ! have liquid condensate, but the part of it that is changing could be the
    ! ice portion; however, this is what was done before.
    lcldo(:ncol,:)  = cldo(:ncol,:)  * cldliqf(:ncol,:)
    lcldn(:ncol,:) = cldn(:ncol,:) * cldliqf(:ncol,:)

    arg = 1.0_r8
    if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
       write(iulog,*) 'erf(1.0) = ',ERF(arg)
       call endrun('dropmixnuc: Error function error')
    endif
    arg = 0.0_r8
    if (erf(arg) /= 0.0_r8) then
       write(iulog,*) 'erf(0.0) = ',erf(arg)
       write(iulog,*) 'dropmixnuc: Error function error'
       call endrun('dropmixnuc: Error function error')
    endif

    dtinv = 1._r8/dtmicro

    allocate(                                         &
         nact(pver,ntot_amode),                       &
         mact(pver,ntot_amode),                       &
         raer(ncnst_tot),                             &
         qqcw(ncnst_tot),                             &
         raercol(pver,ncnst_tot,2),                   &
         raercol_cw(pver,ncnst_tot,2),                &
         coltend(pcols,ncnst_tot),                    &
         coltend_cw(pcols,ncnst_tot),                 &
         naermod(ntot_amode),                         &
         hygro(ntot_amode),                           &
         lnsigman(ntot_amode),                        & !variable std. deviation (CAM-Oslo)
         raercol_tracer(pver,n_aerosol_tracers,2),    &
         raercol_cw_tracer(pver,n_aerosol_tracers,2), &
         mact_tracer(pver,n_aerosol_tracers),         &
         mfullact_tracer(pver,n_aerosol_tracers),     &
         vaerosol(ntot_amode),                        &
         fn(ntot_amode),                              &
         fm(ntot_amode),                              &
         fluxn(ntot_amode),                           &
         fluxm(ntot_amode)               )

    ! Init pointers to mode number and specie mass mixing ratios in
    ! intersitial and cloud borne phases.
    ! Need a list of all aerosol species ==> store in raer (mm)
    ! or qqcw for cloud-borne aerosols (?)
    do imode=1,nmodes  !All aerosol modes

       !NOTE: SEVERAL POINTERS POINT TO SAME FIELD, E.G. CONDENSATE WHICH IS IN SEVERAL MODES
       do ispec = 1, nspec_amode(imode)
          tracerIndex =  tracer_index(imode,ispec)           !Index in q
          cloud_tracer_index = CloudTracerIndex(tracerIndex) !Index in phys-buffer
          mm = mam_idx(imode,ispec)                          !Index in raer/qqcw
          raer(mm)%fld =>  state%q(:,:,tracerIndex)          !NOTE: These are total fields (for example condensate)
          call pbuf_get_field(pbuf, cloud_tracer_index, qqcw(mm)%fld) !NOTE: These are total fields (for example condensate)
       enddo
    enddo
    allocate(                   &
         fn_tmp(ntot_amode),    &
         fm_tmp(ntot_amode),    &
         fluxn_tmp(ntot_amode), &
         fluxm_tmp(ntot_amode))

    ! Initialize turbulent vertical velocity at base of layers (m/s)
    wtke(:,:) = 0._r8

    ! aerosol tendencies
    call physics_ptend_init(ptend, state%psetcols, 'ndrop', lq=lq)

    !Improve this later by using only cloud points ?
    do ilev = top_lev, pver
       do icol=1,ncol
          cs(icol,ilev) = pmid(icol,ilev)/(rair*temp(icol,ilev))        ! air density (kg/m3)
       end do
    end do

    !Output this
    call t_startf('ndrop_calculateNumberMedianRadius')
    call calculateNumberMedianRadius(numberConcentration, volumeConcentration, lnSigma, numberMedianRadius, ncol)
    do imode=1,nmodes
       sigma(:ncol,:,imode) = DEXP(lnSigma(:ncol,:,imode))
       modeString="  "
       write(modeString,"(I2)"),imode
       if(imode < 10) modeString="0"//adjustl(modeString)
       varName = "NMR"//trim(modeString)
       call outfld(varName, numberMedianRadius(:ncol,:,imode), ncol, lchnk)
       varName = "NCONC"//trim(modeString)
       call outfld(varName, numberConcentration(:ncol,:,imode),ncol, lchnk)
       varName = "VCONC"//trim(modeString)
       call outfld(varName, volumeConcentration(:ncol,:,imode), ncol,lchnk)
       varName = "SIGMA"//trim(modeString)
       call outfld(varName, sigma(:ncol,:,imode), ncol,lchnk)
       varName = "HYGRO"//trim(modeString)
       call outfld(varName, hygroscopicity(:ncol,:,imode), ncol,lchnk)
    end do
    call t_stopf('ndrop_calculateNumberMedianRadius')

    alert = .FALSE.
    do ilev=top_lev,pver
       mm = ilev - top_lev + 1
       do imode=1,nmodes
          if(.NOT. alert .and. ANY(numberConcentration(:ncol,ilev,imode) < 0.0_r8 ))then
             alert = .TRUE.
             print*,"STRANGE numberconc", imode, minval(numberConcentration(:ncol,:,:))*1.e-6_r8, "#/cm3", ilev, mm
          endif
       enddo
    enddo
    if (alert)then
       call endrun("strange number concentration")
    endif

    ! tau_cld_regenerate = time scale for regeneration of cloudy air
    ! by (horizontal) exchange with clear air
    tau_cld_regenerate = 3600.0_r8 * 3.0_r8
    tau_cld_regenerate_exp = exp(-dtmicro/tau_cld_regenerate)

    ! overall_main_i_loop
    do icol = 1, ncol

       coltend(icol,:)=0.0_r8
       coltend_cw(icol,:) = 0.0_r8

       do ilev = top_lev, pver-1
          zs(ilev) = 1._r8/(zm(icol,ilev) - zm(icol,ilev+1))
       end do
       zs(pver) = zs(pver-1)

       do ilev = 1, pver
          do imode = 1, ntot_amode
             nact(ilev,imode) = 0._r8
             mact(ilev,imode) = 0._r8
          end do
       end do

       ! load number nucleated into qcld on cloud boundaries
       do ilev = top_lev, pver

          qcld(ilev)  = ncldwtr(icol,ilev)
          qncld(ilev) = 0._r8
          srcn(ilev)  = 0._r8
          cs(icol,ilev)  = pmid(icol,ilev)/(rair*temp(icol,ilev))        ! air density (kg/m3)
          dz(icol,ilev)  = 1._r8/(cs(icol,ilev)*gravit*rpdel(icol,ilev)) ! layer thickness in m

          zn(ilev) = gravit*rpdel(icol,ilev)

          if (ilev < pver) then
             ekd(ilev)   = kvh(icol,ilev+1)
             ekd(ilev)   = max(ekd(ilev), zkmin)
             ekd(ilev)   = min(ekd(ilev), zkmax)
             csbot(ilev) = 2.0_r8*pint(icol,ilev+1)/(rair*(temp(icol,ilev) + temp(icol,ilev+1)))
             csbot_cscen(ilev) = csbot(ilev)/cs(icol,ilev)
          else
             ekd(ilev)   = 0._r8
             csbot(ilev) = cs(icol,ilev)
             csbot_cscen(ilev) = 1.0_r8
          end if

          ! define wtke at layer centers for new-cloud activation
          ! and at layer boundaries for old-cloud activation
          wtke_cen(icol,ilev) = wsub(icol,ilev)
          wtke(icol,ilev)     = wsub(icol,ilev)
          wtke_cen(icol,ilev) = max(wtke_cen(icol,ilev), wmixmin)
          wtke(icol,ilev)     = max(wtke(icol,ilev), wmixmin)
          nsource(icol,ilev) = 0._r8

       end do  ! ilev

       nsav = 1
       nnew = 2

       !get constituent fraction
       call t_startf('ndrop_getConstituentFraction')

       call t_startf('ndrop_getConstituentFraction_calc_aersolFraction')
       componentFractionOK(:,:,:) = 0.0_r8
       do ilev=top_lev, pver
         do imode = 1,ntot_amode
           if (imode <= nbmodes)then
             do ispec = 1, nspec_amode(imode)
               !calculate fraction of component "ispec" in mode "imode" based on concentrations in clear air
               tracerIndex = tracer_index(imode,ispec)

               if (l_so4_a1 == tracerIndex) then !so4 condensation
                  componentFractionOK(imode,tracerIndex,ilev)= (Cam(icol,ilev,imode) &
                       *(1.0_r8 - f_acm(icol,ilev,imode)) & !sulfate fraction
                       *(1.0_r8 - f_aqm(icol,ilev,imode)) & !fraction not from aq phase
                       *(f_so4_condm(icol,ilev,imode)))   & !fraction being condensate
                       /(CProcessModes(icol,ilev)*(1.0_r8-f_c(icol,ilev))*(1.0_r8-f_aq(icol,ilev))*f_so4_cond(icol,ilev)+smallNumber) !total so4 condensate

               else if (l_so4_ac == tracerIndex) then ! so4 coagulation
                  componentFractionOK(imode,tracerIndex,ilev) = (Cam(icol,ilev,imode) &
                       *(1.0_r8 - f_acm(icol,ilev,imode)) &         !sulfate fraction
                       *(1.0_r8 - f_aqm(icol,ilev,imode)) &         !fraction not from aq phase
                       *(1.0_r8 - f_so4_condm(icol,ilev,imode))) &  !fraction not being condensate
                       /(CProcessModes(icol,ilev)*(1.0_r8-f_c(icol,ilev))*(1.0_r8-f_aq(icol,ilev))*(1.0_r8-f_so4_cond(icol,ilev)) + smallNumber)

               else if (l_so4_a2 == tracerIndex) then  !so4 wet phase
                  componentFractionOK(imode,tracerIndex,ilev) = (Cam(icol,ilev,imode) &
                       *(1.0_r8-f_acm(icol,ilev,imode)) & !sulfate fraction
                       *f_aqm(icol,ilev,imode))         & !aq phase fraction of sulfate
                       /(CProcessModes(icol,ilev)*(1.0_r8-f_c(icol,ilev))*(f_aq(icol,ilev))+smallNumber)

               else if (l_bc_ac == tracerIndex)then  !bc coagulated
                  componentFractionOK(imode,tracerIndex,ilev) = (Cam(icol,ilev,imode) &
                       *f_acm(icol,ilev,imode)  & ! carbonaceous fraction
                       *f_bcm(icol,ilev,imode)) & ! bc fraction of carbonaceous
                       /(CProcessModes(icol,ilev)*f_c(icol,ilev)*f_bc(icol,ilev)+smallNumber)

               else if (l_om_ac == tracerIndex ) then  !oc coagulated
                  componentFractionOK(imode,tracerIndex,ilev) = (Cam(icol,ilev,imode) &
                       *f_acm(icol,ilev,imode)            &  ! carbonaceous fraction
                       *(1.0_r8-f_bcm(icol,ilev,imode))   &  ! oc fraction of carbonaceous
                       *(1.0_r8-f_soam(icol,ilev,imode))) &  ! oc fraction which is soa
                       /(CProcessModes(icol,ilev)*f_c(icol,ilev)*(1.0_r8-f_bc(icol,ilev))*(1.0_r8-f_soa(icol,ilev))+smallNumber)

               else if (l_soa_a1 == tracerIndex) then !SOA condensate
                  componentFractionOK(imode,tracerIndex,ilev) = Cam(icol,ilev,imode) &
                       *f_acm(icol,ilev,imode)           &  !carbonaceous fraction
                       *(1.0_r8 -f_bcm(icol,ilev,imode)) &  !om fraction
                       *(f_soam(icol,ilev,imode))        &  !fraction of OM is SOA
                       /(CProcessModes(icol,ilev)*f_c(icol,ilev)*(1.0_r8 -f_bc(icol,ilev))*f_soa(icol,ilev) + smallNumber)
               end if
               if (componentFractionOK(imode,tracerIndex,ilev) >  1.0_r8)then
                  componentFractionOK(imode,tracerIndex,ilev) = 1.0_r8
               endif
             end do
           else
             do ispec = 1, nspec_amode(imode)
               tracerIndex = tracer_index(imode,ispec)
               componentFractionOK(imode,tracerIndex,ilev) = 1.0_r8
             end do
           endif
         end do
       end do
       call t_stopf('ndrop_getConstituentFraction_calc_aersolFraction')

       !Loop over all tracers ==> check that sums to one
       !for all tracers which exist in the oslo-modes

       call t_startf('ndrop_getConstituentFraction_check_trackerNormalization1')
       sumFraction(:,:) = 0.0_r8
       do ilev=top_lev, pver
         do itrac=1,pcnst
           do imode=1,ntot_amode
             sumFraction(itrac,ilev) = sumFraction(itrac,ilev) + componentFractionOK(imode,itrac,ilev)
           end do
         end do
       end do
       call t_stopf('ndrop_getConstituentFraction_check_trackerNormalization1')

       call t_startf('ndrop_getConstituentFraction_check_trackerNormalization2')
       do ilev=top_lev, pver
         do itrac=1,pcnst
           if (sumFraction(itrac,ilev) > 1.e-2_r8) then
             !Just scale what comes out if componentFraction is larger than 1%
             do imode=1,ntot_amode
               componentFractionOK(imode,itrac,ilev) = componentFractionOK(imode,itrac,ilev)/sumFraction(itrac,ilev)
             end do
           else
             ! negative or zero fraction for this species
             ! distribute equal fraction to all receiver modes
             do imode=1,ntot_amode
               componentFractionOK(imode,itrac,ilev)=1.0_r8/max(1.e-30_r8, sumFraction2(itrac,ilev))
             end do !modes
           endif
         end do !tracers
       end do  !levels
       call t_stopf('ndrop_getConstituentFraction_check_trackerNormalization2')

       call t_stopf('ndrop_getConstituentFraction')

       call t_startf('ndrop_getNumberConc')
       do imode = 1, nmodes ! Number of modes
          !Get number concentration of this mode
          mm = mam_idx(imode,0)
          raercol(1:top_lev-1,mm,nsav) = 0.0_r8
          raercol_cw(1:top_lev-1,mm,nsav) = 0.0_r8
          do ilev= top_lev,pver
             raercol(ilev,mm,nsav) = numberConcentration(icol,ilev,imode)/cs(icol,ilev) !#/kg air
             !In oslo model, number concentrations are diagnostics, so
             !approximate number concentration in each mode by total
             !cloud number concentration scaled by how much is available of
             !each mode
             raercol_cw(ilev,mm,nsav) = ncldwtr(icol,ilev)*numberConcentration(icol,ilev,imode) &
                                    /max(1.e-30_r8, sum(numberConcentration(icol,ilev,1:nmodes)))
          enddo

          !These are the mass mixing ratios
          do ispec = 1, nspec_amode(imode)
             mm = mam_idx(imode,ispec)      !index of tracer (all unique)
             raercol(:,mm,nsav) = 0.0_r8
             raercol_cw(:,mm,nsav) = 0.0_r8
             !Several of the fields (raer(mm)%fld point to the same
             !field in q. To avoid double counting, we take into
             !account the component fraction in the mode
             do ilev=top_lev,pver
                if(imode > nbmodes) then
                   componentFraction = 1.0_r8
                else
                   tracerIndex = tracer_index(imode,ispec)
                   componentFraction = componentFractionOK(imode,tracerIndex,ilev)
                endif
                !Assign to the components used here i.e. distribute condensate/coagulate to modes
                raercol_cw(ilev,mm,nsav) = qqcw(mm)%fld(icol,ilev)*componentFraction
                raercol(ilev,mm,nsav)    = raer(mm)%fld(icol,ilev)*componentFraction
             enddo ! ilev (levels)
          end do   ! ispec (species)
       end do      ! imode (modes)
       call t_stopf('ndrop_getNumberConc')

       ! droplet nucleation/aerosol activation
       call t_startf('ndrop_nucleation_activation')

       ! ilev-loop for growing/shrinking cloud calcs .............................
       ! grow_shrink_main_k_loop: &
       do ilev = top_lev, pver

          ! This code was designed for liquid clouds, but the cloudbourne
          ! aerosol can be either from liquid or ice clouds. For the ice clouds,
          ! we do not do regeneration, but as cloud fraction decreases the
          ! aerosols should be returned interstitial. The lack of a liquid cloud
          ! should not mean that all of the aerosol is realease. Therefor a
          ! section has been added for shrinking ice clouds and checks were added
          ! to protect ice cloudbourne aerosols from being released when no
          ! liquid cloud is present.

          ! shrinking ice cloud ......................................................
          cldo_tmp = cldo(icol,ilev) * (1._r8 - cldliqf(icol,ilev))
          cldn_tmp = cldn(icol,ilev) * (1._r8 - cldliqf(icol,ilev))

          if (cldn_tmp < cldo_tmp) then

             ! convert activated aerosol to interstitial in decaying cloud

             dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * (1._r8 - cldliqf(icol,ilev))
             do imode = 1, ntot_amode
                mm = mam_idx(imode,0)
                dact   = raercol_cw(ilev,mm,nsav)*dumc
                raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact   ! cloud-borne aerosol
                raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                do ispec = 1, nspec_amode(imode)
                   mm = mam_idx(imode,ispec)
                   dact    = raercol_cw(ilev,mm,nsav)*dumc
                   raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact  ! cloud-borne aerosol
                   raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                end do
             end do
          end if

          ! shrinking liquid cloud ......................................................
          !    treat the reduction of cloud fraction from when cldn(icol,ilev) < cldo(icol,ilev)
          !    and also dissipate the portion of the cloud that will be regenerated
          cldo_tmp = lcldo(icol,ilev)
          cldn_tmp = lcldn(icol,ilev) * tau_cld_regenerate_exp
          !    alternate formulation
          !    cldn_tmp = cldn(icol,ilev) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

          ! fraction is also provided.
          if (cldn_tmp < cldo_tmp) then
             !  droplet loss in decaying cloud
             nsource(icol,ilev) = nsource(icol,ilev) + qcld(ilev)*(cldn_tmp - cldo_tmp)/cldo_tmp*cldliqf(icol,ilev)*dtinv
             qcld(ilev) = qcld(ilev)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)

             ! convert activated aerosol to interstitial in decaying cloud
             dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * cldliqf(icol,ilev)
             do imode = 1, ntot_amode
                mm = mam_idx(imode,0)
                dact   = raercol_cw(ilev,mm,nsav)*dumc
                raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact   ! cloud-borne aerosol
                raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                do ispec = 1, nspec_amode(imode)
                   mm = mam_idx(imode,ispec)
                   dact    = raercol_cw(ilev,mm,nsav)*dumc
                   raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact  ! cloud-borne aerosol
                   raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                end do
             end do
          end if

          ! growing liquid cloud ......................................................
          !    treat the increase of cloud fraction from when cldn(icol,ilev) > cldo(icol,ilev)
          !    and also regenerate part of the cloud
          cldo_tmp = cldn_tmp
          cldn_tmp = lcldn(icol,ilev)

          if (cldn_tmp-cldo_tmp > 0.01_r8) then

             ! use wtke at layer centers for new-cloud activation
             wbar  = wtke_cen(icol,ilev)
             wmix  = 0._r8
             wmin  = 0._r8
             wmax  = 10._r8
             wdiab = 0._r8

             ! load aerosol properties, assuming external mixtures
             naermod(:) = 0.0_r8
             vaerosol(:) = 0.0_r8
             hygro(:) = 0.0_r8
             lnsigman(:) = alog2

             mode_cnt = 0
             do imode = 1,nmodes
                if(hasAerosol(icol,ilev,imode)) then
                   mode_cnt = mode_cnt + 1
                   naermod(mode_cnt) = numberConcentration(icol,ilev,imode)
                   vaerosol(mode_cnt) = volumeConcentration(icol,ilev,imode)
                   hygro(mode_cnt) = hygroscopicity(icol,ilev,imode)
                   lnsigman(mode_cnt) = lnsigma(icol,ilev,imode)
                   speciesMap(mode_cnt) = imode
                end if
             end do
             numberOfModes = mode_cnt

             ! Call the activation procedure
             if (numberOfModes > 0)then
                if (use_hetfrz_classnuc) then
                   call activate_modal_oslo( wbar, wmix, wdiab, wmin, wmax,       &
                        temp(icol,ilev), cs(icol,ilev), naermod, numberOfModes,          &
                        vaerosol, hygro, fn_in(icol,ilev,1:nmodes), fm, fluxn,     &
                        fluxm, flux_fullact(ilev), lnsigman)
                else
                   call activate_modal_oslo( wbar, wmix, wdiab, wmin, wmax,       &
                        temp(icol,ilev), cs(icol,ilev), naermod, numberOfModes,          &
                        vaerosol, hygro, fn, fm, fluxn,                      &
                        fluxm, flux_fullact(ilev), lnsigman)
                end if
             endif

             dumc = (cldn_tmp - cldo_tmp)

             if (use_hetfrz_classnuc) then
                fn_tmp(:) = fn_in(icol,ilev,1:nmodes)
             else
                fn_tmp(:) = fn(:)
             end if
             fm_tmp(:) = fm(:)
             fluxn_tmp(:) = fluxn(:)
             fluxm_tmp(:) = fluxm(:)
             fn(:) = 0.0_r8
             fn_in(icol,ilev,:) = 0.0_r8
             fm(:) = 0.0_r8
             fluxn(:)=0.0_r8
             fluxm(:)= 0.0_r8
             do imode = 1, numberOfModes   !Number of coexisting modes to be used for activation
                kcomp = speciesMap(imode)       !This is the CAM-oslo mode (modes 1-14 may be activated, mode 0 not)
                if (use_hetfrz_classnuc) then
                   fn_in(icol,ilev,kcomp) = fn_tmp(imode)
                else
                   fn(kcomp) = fn_tmp(imode)
                end if
                fm(kcomp) = fm_tmp(imode)
                fluxn(kcomp) = fluxn_tmp(imode)
                fluxm(kcomp) = fluxm_tmp(imode)
             enddo
             do imode = 1, ntot_amode
                mm = mam_idx(imode,0)
                if (use_hetfrz_classnuc) then
                   dact   = dumc*fn_in(icol,ilev,imode)*numberConcentration(icol,ilev,imode)/cs(icol,ilev) !#/kg_{air}
                else
                   dact   = dumc*fn(imode)*numberConcentration(icol,ilev,imode)/cs(icol,ilev) !#/kg_{air}
                end if
                qcld(ilev) = qcld(ilev) + dact
                nsource(icol,ilev) = nsource(icol,ilev) + dact*dtinv
                raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact  ! cloud-borne aerosol
                raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                dum = dumc*fm(imode)
                do ispec = 1, nspec_amode(imode)
                   mm = mam_idx(imode,ispec)
                   if(imode > nbmodes)then
                      constituentFraction = 1.0_r8
                   else
                      tracerIndex = tracer_index(imode,ispec)
                      constituentFraction = componentFractionOK(imode,tracerIndex,ilev)
                   endif

                   dact    = dum*raer(mm)%fld(icol,ilev)*constituentFraction
                   raercol_cw(ilev,mm,nsav) = raercol_cw(ilev,mm,nsav) + dact  ! cloud-borne aerosol
                   raercol(ilev,mm,nsav)    = raercol(ilev,mm,nsav) - dact
                enddo
             enddo
          endif  ! cldn_tmp-cldo_tmp > 0.01_r8

       enddo  ! grow_shrink_main_k_loop
       ! end of ilev-loop for growing/shrinking cloud calcs ......................
       call t_stopf('ndrop_nucleation_activation')

       ! ......................................................................
       ! start of ilev-loop for calc of old cloud activation tendencies ..........
       !
       ! use current cloud fraction (cldn) exclusively
       ! consider case of cldo(:)=0, cldn(ilev)=1, cldn(ilev+1)=0
       ! previous code (which used cldo below here) would have no cloud-base activation
       ! into layer ilev.  however, activated particles in ilev mix out to ilev+1,
       ! so they are incorrectly depleted with no replacement

       call t_startf('ndrop_oldcloud_activation')
       ! old_cloud_main_k_loop
       do ilev = top_lev, pver
          kp1 = min0(ilev+1, pver)
          taumix_internal_pver_inv = 0.0_r8

          if (lcldn(icol,ilev) > 0.01_r8) then

             wdiab = 0._r8
             wmix  = 0._r8                       ! single updraft
             wbar  = wtke(icol,ilev)                   ! single updraft
             if (ilev == pver) wbar = wtke_cen(icol,ilev) ! single updraft
             wmax  = 10._r8
             wmin  = 0._r8

             if (lcldn(icol,ilev) - lcldn(icol,kp1) > 0.01_r8 .or. ilev == pver) then

                ! cloud base
                ! ekd(ilev) = wtke(icol,ilev)*dz(icol,ilev)/sq2pi
                ! first, should probably have 1/zs(ilev) here rather than dz(icol,ilev) because
                !    the turbulent flux is proportional to ekd(ilev)*zs(ilev),
                !    while the dz(icol,ilev) is used to get flux divergences
                !    and mixing ratio tendency/change
                ! second and more importantly, using a single updraft velocity here
                !    means having monodisperse turbulent updraft and downdrafts.
                !    The sq2pi factor assumes a normal draft spectrum.
                !    The fluxn/fluxm from activate must be consistent with the
                !    fluxes calculated in explmix.
                ekd(ilev) = wbar/zs(ilev)

                alogarg = max(1.e-20_r8, 1._r8/lcldn(icol,ilev) - 1._r8)
                wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)
                phase   = 1   ! interstitial
                naermod(:) = 0.0_r8
                vaerosol(:) = 0.0_r8
                hygro(:) = 0.0_r8
                lnsigman(:) = alog2

                mode_cnt=0
                do imode = 1,nmodes
                   if(hasAerosol(icol,kp1,imode) .eqv. .TRUE.)then
                      mode_cnt = mode_cnt + 1
                      naermod(mode_cnt) = numberConcentration(icol,kp1,imode)
                      vaerosol(mode_cnt) = volumeConcentration(icol,kp1,imode)
                      hygro(mode_cnt) =    hygroscopicity(icol,kp1,imode)
                      lnsigman(mode_cnt) = lnsigma(icol,kp1,imode)
                      speciesMap(mode_cnt) = imode
                   end if
                end do
                numberOfModes = mode_cnt

                if(numberOfModes > 0)then
                   if (use_hetfrz_classnuc) then
                      call activate_modal_oslo(wbar, wmix, wdiab, wmin, wmax, &
                           temp(icol,ilev), cs(icol,ilev), naermod, numberOfModes ,  &
                           vaerosol, hygro, fn_in(icol,ilev,:), fm, fluxn,     &
                           fluxm, flux_fullact(ilev), lnsigman)
                   else
                      call activate_modal_oslo(wbar, wmix, wdiab, wmin, wmax, &
                           temp(icol,ilev), cs(icol,ilev), naermod, numberOfModes ,  &
                           vaerosol, hygro, fn, fm, fluxn,               &
                           fluxm, flux_fullact(ilev), lnsigman)
                   end if
                endif

                !Difference in cloud fraction this layer and above!
                !we are here because there are more clouds above, and some
                !aerosols go into  that layer! ==> calculate additional cloud fraction
                if (ilev < pver) then
                   dumc = lcldn(icol,ilev) - lcldn(icol,kp1)
                else
                   dumc = lcldn(icol,ilev)
                endif

                if (use_hetfrz_classnuc) then
                   fn_tmp(:) = fn_in(icol,ilev,1:nmodes)
                else
                   fn_tmp(:) = fn(:)
                end if
                fm_tmp(:) = fm(:)
                fluxn_tmp(:) = fluxn(:)
                fluxm_tmp(:) = fluxm(:)
                fn(:) = 0.0_r8
                fn_in(icol,ilev,:) = 0.0_r8
                fm(:) = 0.0_r8
                fluxn(:)=0.0_r8
                fluxm(:)= 0.0_r8

                do imode = 1, numberOfModes        !Number of coexisting modes to be used for activation
                   kcomp = speciesMap(imode)       !This is the CAM-oslo mode (modes 1-14 may be activated, mode 0 not)
                   if (use_hetfrz_classnuc) then
                      fn_in(icol,ilev,kcomp) = fn_tmp(imode)
                   else
                      fn(kcomp) = fn_tmp(imode)
                   end if
                   fm(kcomp) = fm_tmp(imode)
                   fluxn(kcomp) = fluxn_tmp(imode)
                   fluxm(kcomp) = fluxm_tmp(imode)
                enddo

                fluxntot = 0.0_r8

                ! flux of activated mass into layer ilev (in kg/m2/s)
                !    = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(ilev)
                ! source of activated mass (in kg/kg/s) = flux divergence
                !    = actmassflux/(cs(icol,ilev)*dz(icol,ilev))
                ! so need factor of csbot_cscen = csbot(ilev)/cs(icol,ilev)
                !                            dum=1./(dz(icol,ilev))
                dum=csbot_cscen(ilev)/(dz(icol,ilev))

                ! code for ilev=pver was changed to use the following conceptual model
                ! in ilev=pver, there can be no cloud-base activation unless one considers
                !    a scenario such as the layer being partially cloudy,
                !    with clear air at bottom and cloudy air at top
                ! assume this scenario, and that the clear/cloudy portions mix with
                !    a timescale taumix_internal = dz(icol,pver)/wtke_cen(icol,pver)
                ! in the absence of other sources/sinks, qact (the activated particle
                !    mixratio) attains a steady state value given by
                !       qact_ss = fcloud*fact*qtot
                !    where fcloud is cloud fraction, fact is activation fraction,
                !    qtot=qact+qint, qint is interstitial particle mixratio
                ! the activation rate (from mixing within the layer) can now be
                !    written as
                !       d(qact)/dt = (qact_ss - qact)/taumix_internal
                !                  = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
                ! note that (fcloud*fact*wtke/dz) is equal to the nact/mact
                ! also, d(qact)/dt can be negative.  in the code below
                !    it is forced to be >= 0
                !
                ! steve --
                !    you will likely want to change this.  i did not really understand
                !       what was previously being done in ilev=pver
                !    in the cam3_5_3 code, wtke(icol,pver) appears to be equal to the
                !       droplet deposition velocity which is quite small
                !    in the cam3_5_37 version, wtke is done differently and is much
                !       larger in ilev=pver, so the activation is stronger there
                !
                if (ilev == pver) then
                   taumix_internal_pver_inv = flux_fullact(ilev)/dz(icol,ilev)
                end if

                do imode = 1, ntot_amode
                   mm = mam_idx(imode,0)
                   fluxn(imode) = fluxn(imode)*dumc
                   fluxm(imode) = fluxm(imode)*dumc
                   nact(ilev,imode) = nact(ilev,imode) + fluxn(imode)*dum
                   mact(ilev,imode) = mact(ilev,imode) + fluxm(imode)*dum
                   if (ilev < pver) then
                      ! note that kp1 is used here
                      fluxntot = fluxntot &
                           + fluxn(imode)*raercol(kp1,mm,nsav)*cs(icol,ilev)
                   else
                      tmpa = raercol(kp1,mm,nsav)*fluxn(imode) &
                           + raercol_cw(kp1,mm,nsav)*(fluxn(imode) &
                           - taumix_internal_pver_inv*dz(icol,ilev))
                      fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(icol,ilev)
                   end if
                end do
                srcn(ilev)      = srcn(ilev) + fluxntot/(cs(icol,ilev)*dz(icol,ilev))
                nsource(icol,ilev) = nsource(icol,ilev) + fluxntot/(cs(icol,ilev)*dz(icol,ilev))
             endif  ! (cldn(icol,ilev) - cldn(icol,kp1) > 0.01 .or. ilev == pver)

          else  ! i.e: cldn(icol,ilev) < 0.01_r8

             ! no liquid cloud
             nsource(icol,ilev) = nsource(icol,ilev) - qcld(ilev)*dtinv
             qcld(ilev) = 0.0_r8

             if (cldn(icol,ilev) < 0.01_r8) then
                ! no ice cloud either
                ! convert activated aerosol to interstitial in decaying cloud
                do imode = 1, ntot_amode
                   mm = mam_idx(imode,0)
                   raercol(ilev,mm,nsav) = raercol(ilev,mm,nsav) + raercol_cw(ilev,mm,nsav)  ! cloud-borne aerosol
                   raercol_cw(ilev,mm,nsav) = 0._r8

                   do ispec = 1, nspec_amode(imode)
                      mm = mam_idx(imode,ispec)
                      raercol(ilev,mm,nsav) = raercol(ilev,mm,nsav) + raercol_cw(ilev,mm,nsav) ! cloud-borne aerosol
                      raercol_cw(ilev,mm,nsav) = 0._r8
                   end do
                end do
             end if
          end if

       end do  ! old_cloud_main_k_loop
       call t_stopf('ndrop_oldcloud_activation')

       ! switch nsav, nnew so that nnew is the updated aerosol
       ntemp = nsav
       nsav  = nnew
       nnew  = ntemp

       call t_startf('ndrop_newdroplets')
       ! load new droplets in layers above, below clouds
       dtmin = dtmicro
       ekk(top_lev-1) = 0.0_r8
       do ilev = top_lev, pver-1
          ! ekd(ilev) is eddy-diffusivity at ilev/ilev+1 interface
          ! want ekk(ilev) = ekd(ilev) * (density at ilev/ilev+1 interface)
          ekk(ilev) = ekd(ilev)*csbot(ilev)
       end do
       ekk(pver) = 0.0_r8

       do ilev = top_lev, pver
          km1 = max0(ilev-1, top_lev)
          ekkp(ilev) = zn(ilev)*ekk(ilev)*zs(ilev)
          ekkm(ilev) = zn(ilev)*ekk(ilev-1)*zs(km1)
          tinv    = ekkp(ilev) + ekkm(ilev)
          ! tinv is the sum of all first-order-loss-rates
          ! for the layer.  for most layers, the activation loss rate
          ! (for interstitial particles) is accounted for by the loss by
          ! turb-transfer to the layer above.
          ! ilev=pver is special, and the loss rate for activation within
          ! the layer must be added to tinv.  if not, the time step
          ! can be too big, and explmix can produce negative values.
          ! the negative values are reset to zero, resulting in an
          ! artificial source.
          if (ilev == pver) then
             tinv = tinv + taumix_internal_pver_inv
          end if
          if (tinv > 1.e-6_r8) then
             dtt   = 1._r8/tinv
             dtmin = min(dtmin, dtt)
          end if
       end do

       dtmix   = 0.9_r8*dtmin
       nsubmix = dtmicro/dtmix + 1
       if (nsubmix > 100) then
          nsubmix_bnd = 100
       else
          nsubmix_bnd = nsubmix
       end if
       count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
       dtmix = dtmicro/nsubmix

       do ilev = top_lev, pver
          kp1 = min(ilev+1, pver)
          km1 = max(ilev-1, top_lev)
          ! maximum overlap assumption
          if (cldn(icol,kp1) > 1.e-10_r8) then
             overlapp(ilev) = min(cldn(icol,ilev)/cldn(icol,kp1), 1._r8)
          else
             overlapp(ilev) = 1._r8
          end if
          if (cldn(icol,km1) > 1.e-10_r8) then
             overlapm(ilev) = min(cldn(icol,ilev)/cldn(icol,km1), 1._r8)
          else
             overlapm(ilev) = 1._r8
          end if
       end do

       !    the activation source(ilev) = mact(ilev,imode)*raercol(kp1,lmass)
       !       should not exceed the rate of transfer of unactivated particles
       !       from kp1 to ilev which = ekkp(ilev)*raercol(kp1,lmass)
       !    however it might if things are not "just right" in subr activate
       !    the following is a safety measure to avoid negatives in explmix
       do imode = 1, ntot_amode
          do ilev = top_lev, pver-1
             nact(ilev,imode) = min( nact(ilev,imode), ekkp(ilev) )
             mact(ilev,imode) = min( mact(ilev,imode), ekkp(ilev) )
          end do
       end do

       !Don't need the mixing per mode in OSLO_AERO ==> only per tracer
       !Note that nsav/nnew is switched above, so operate on nnew here
       !nnew is the updated aerosol
       raercol_tracer(:,:,:) = 0.0_r8
       raercol_cw_tracer(:,:,:) = 0.0_r8
       mact_tracer(:,:) = 0.0_r8
       mfullact_tracer(:,:) = 0.0_r8
       do imode=1,ntot_amode
          do ispec=1,nspec_amode(imode)
             itrac = tracer_index(imode,ispec)  !which tracer are we talking about
             itrac2  = inverseAerosolTracerList(itrac)    !which index is this in the list of aerosol-tracers
             mm = mam_idx(imode,ispec)
             do ilev = top_lev,pver
                raercol_tracer(ilev,itrac2,nnew) = raercol_tracer(ilev,itrac2,nnew) + raercol(ilev,mm,nnew)
                raercol_cw_tracer(ilev,itrac2,nnew) = raercol_cw_tracer(ilev,itrac2,nnew) + raercol_cw(ilev,mm,nnew)
                mact_tracer(ilev,itrac2) = mact_tracer(ilev,itrac2) + mact(ilev,imode)*raercol(ilev,mm,nnew)
                mfullact_tracer(ilev,itrac2) = mfullact_tracer(ilev,itrac2) + raercol(ilev,mm,nnew)
             end do
          end do !ispec
       end do    !imode

       do itrac=1,n_aerosol_tracers
          mact_tracer(:,itrac) = mact_tracer(:,itrac) /(mfullact_tracer(:,itrac) + smallNumber)
       end do
       call t_stopf('ndrop_newdroplets')

       ! old_cloud_nsubmix_loop
       call t_startf('ndrop_oldcloud_nsubmix')
       do isubmix = 1, nsubmix
          qncld(:) = qcld(:)

          ! switch nsav, nnew so that nsav is the updated aerosol
          ntemp   = nsav
          nsav    = nnew
          nnew    = ntemp

          !First mix cloud droplet number concentration
          call t_startf('ndrop_oldcloud_nsubmix_mix_CloudDropNumConc')
          srcn(:) = 0.0_r8
          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             do ilev = top_lev,pver
                if (ilev < pver) then
                   ! update droplet source - activation source in layer ilev involves particles from ilev+1
                   srcn(ilev) = srcn(ilev) + nact(ilev,imode)*(raercol(ilev+1,mm,nsav))
                else
                   ! new formulation for ilev=pver
                   tmpa = raercol(ilev,mm,nsav)*nact(ilev,imode) &
                        + raercol_cw(ilev,mm,nsav)*(nact(ilev,imode)-taumix_internal_pver_inv)
                   srcn(ilev) = srcn(ilev) + max(0.0_r8,tmpa)
                end if
             end do
          end do
          call t_stopf('ndrop_oldcloud_nsubmix_mix_CloudDropNumConc')

          !mixing of cloud droplets
          call t_startf('ndrop_oldcloud_nsubmix_mix_CloudDrop')
          call explmix_oslo(qcld, srcn, ekkp, ekkm, overlapp, overlapm, qncld, zero, zero, pver, dtmix, .false.)

          ! Mix number concentrations consistently!!
          do imode = 1, ntot_amode
             mm = mam_idx(imode,0)
             ! activation source in layer ilev involves particles from ilev+1
             source(top_lev:pver-1) = nact(top_lev:pver-1,imode)*(raercol(top_lev+1:pver,mm,nsav))

             ! new formulation for ilev=pver
             tmpa = raercol(pver,mm,nsav)*nact(pver,imode) &
                  + raercol_cw(pver,mm,nsav)*(nact(pver,imode) - taumix_internal_pver_inv)
             source(pver) = max(0.0_r8, tmpa)
             flxconv = 0._r8

             call explmix_oslo( raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw(:,mm,nsav), zero, zero, pver, dtmix, .false.)

             call explmix_oslo( raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol(:,mm,nsav), zero, flxconv, pver, dtmix, .true., raercol_cw(:,mm,nsav))
          end do

          do itrac=1,n_aerosol_tracers
             source(top_lev:pver-1) = mact_tracer(top_lev:pver-1,itrac) &
                                     *(raercol_tracer(top_lev+1:pver,itrac,nsav))

             tmpa = raercol_tracer(pver,itrac,nsav)*mact_tracer(pver,itrac) &
                  + raercol_cw_tracer(pver,itrac,nsav)*(mact_tracer(pver,itrac) - taumix_internal_pver_inv)

             source(pver) = max(0.0_r8, tmpa)
             flxconv = 0.0_r8

             call explmix_oslo(raercol_cw_tracer(:,itrac,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw_tracer(:,itrac,nsav), zero, zero, pver,  dtmix, .false.)

             call explmix_oslo(raercol_tracer(:,itrac,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol_tracer(:,itrac,nsav), zero, flxconv, pver, dtmix, .true., &
                  raercol_cw_tracer(:,itrac,nsav))
          end do !Number of aerosol tracers
          call t_stopf('ndrop_oldcloud_nsubmix_mix_CloudDrop')

       end do ! old_cloud_nsubmix_loop
       call t_stopf('ndrop_oldcloud_nsubmix')

       call t_startf('ndrop_rest')
       ! Set back to the original framework
       ! Could probably continue in tracer-space from here
       ! but return back to mixture for easier use of std. NCAR code
       tendencyCounted(:)=.FALSE.
       do imode = 1, ntot_amode
          do ispec=1,nspec_amode(imode)
             mm=mam_idx(imode,ispec)
             itrac = tracer_index(imode,ispec)
             itrac2 = inverseAerosolTracerList(itrac)
             !All the tracer-space contains sum of all
             !modes ==> put in first available component
             !and zero in others.
             if(.not.tendencyCounted(itrac))then
                raercol(:,mm,nnew) = raercol_tracer(:,itrac2,nnew)
                raercol_cw(:,mm,nnew) = raercol_cw_tracer(:,itrac2,nnew)
                tendencyCounted(itrac) = .TRUE.
             else
                raercol(:,mm,nnew) = 0.0_r8
                raercol_cw(:,mm,nnew) = 0.0_r8
             end if
          end do
       end do

       ! evaporate particles again if no cloud
       do ilev = top_lev, pver
          if (cldn(icol,ilev) == 0._r8) then
             ! no ice or liquid cloud
             qcld(ilev)=0._r8

             ! convert activated aerosol to interstitial in decaying cloud
             do imode = 1, ntot_amode
                mm = mam_idx(imode,0)
                raercol(ilev,mm,nnew)    = raercol(ilev,mm,nnew) + raercol_cw(ilev,mm,nnew)
                raercol_cw(ilev,mm,nnew) = 0._r8

                do ispec = 1, nspec_amode(imode)
                   mm = mam_idx(imode,ispec)
                   raercol(ilev,mm,nnew)    = raercol(ilev,mm,nnew) + raercol_cw(ilev,mm,nnew)
                   raercol_cw(ilev,mm,nnew) = 0._r8
                end do
             end do
          end if
       end do

       ! droplet number
       ndropcol(icol) = 0._r8

       !Initialize tendnd to zero in all layers since values are set in only top_lev,pver
       !Without this the layers above top_lev would be un-initialized
       tendnd(icol,:) = 0.0_r8

       do ilev = top_lev, pver
          ndropmix(icol,ilev) = (qcld(ilev) - ncldwtr(icol,ilev))*dtinv - nsource(icol,ilev)
          tendnd(icol,ilev) = (max(qcld(ilev), 1.e-6_r8) - ncldwtr(icol,ilev))*dtinv
          ndropcol(icol) = ndropcol(icol) + ncldwtr(icol,ilev)*pdel(icol,ilev)
       end do
       ndropcol(icol) = ndropcol(icol)/gravit

       raertend = 0._r8
       qqcwtend = 0._r8

       coltend_cw(icol,:)=0.0_r8
       coltend(icol,:) = 0.0_r8

       !Need to initialize first because process modes arrive several times
       tendencyCounted(:) = .FALSE.
       do imode=1,ntot_amode
         do ispec = 1,nspec_amode(imode)
           itrac = tracer_index(imode,ispec)
           mm = mam_idx(imode,ispec)

           !column tendencies for output
           if(.NOT. tendencyCounted(itrac))then
             coltend_cw(icol,itrac) = coltend_cw(icol,itrac) &
                  + sum( pdel(icol,top_lev:pver)*(raercol_cw(top_lev:pver,mm,nnew) & !New, splitted,
                  - qqcw(mm)%fld(icol,top_lev:pver) ) )/gravit*dtinv      !Old, total
             tendencyCounted(itrac) = .TRUE.
           else  !Already subtracted total old value, just add new
             coltend_cw(icol,itrac) = coltend_cw(icol,itrac)  &
                  + sum(pdel(icol,top_lev:pver)*raercol_cw(top_lev:pver,mm,nnew))/gravit*dtinv !total already subtracted
           end if

           ptend%q(icol,:,itrac) = 0.0_r8  !Initialize tendencies
           qqcw(mm)%fld(icol,:) = 0.0_r8  !Throw out old concentrations before summing new ones
         end do  ! Tracers
       end do     ! Modes

       !First, sum up all the tracer mass concentrations
       do imode = 1, ntot_amode
         do ispec = 1, nspec_amode(imode)
           mm   = mam_idx(imode,ispec)                !tracer indices for aerosol mass mixing ratios in raer-arrays
           itrac = tracer_index(imode,ispec)           !index in q-array (1-pcnst)

           !This is a bit tricky since in our scheme the tracers can arrive several times
           !the same tracer can exist in several modes, e.g. condensate!!
           !Here we sum this into "qqcw" and "ptend" so that they contain TOTAL of those tracers

           !raercol and raercol_cw do not have totals, they have process-tracers splitted onto modes

           !Tendency at this point is the sum (original value subtracted below)
           ptend%q(icol,top_lev:pver,itrac) = ptend%q(icol,top_lev:pver,itrac) + raercol(top_lev:pver,mm,nnew)
           !for cloud water concentrations, we don't get tendency , only new concentration
           qqcw(mm)%fld(icol,top_lev:pver) = qqcw(mm)%fld(icol,top_lev:pver) + raercol_cw(top_lev:pver,mm,nnew)

         end do
       end do

       !Need this check due to some tracers (e.g. condensate) several times
       tendencyCounted(:) = .FALSE.

       ! Recalculating cloud-borne aerosol number mixing ratios
       do imode=1,ntot_amode
         !Now that all new aerosol masses are summed up, we subtract the original concentrations to obtain the tendencies
         do ispec= 1,nspec_amode(imode)
           mm = mam_idx(imode,ispec)
           itrac = tracer_index(imode,ispec)
           if(.NOT. tendencyCounted(itrac)) then
             ptend%q(icol,top_lev:pver,itrac) = (ptend%q(icol,top_lev:pver,itrac) - raer(mm)%fld(icol,top_lev:pver))*dtinv
             coltend(icol,itrac) = sum(pdel(icol,top_lev:pver)*ptend%q(icol,top_lev:pver,itrac))/gravit !Save column tendency
             tendencyCounted(itrac) = .TRUE.
           endif
         end do !species
       end do !modes
       call t_stopf('ndrop_rest')

    end do  ! overall_main_icol_loop

    ! end of main loop over icol/longitude ....................................

    call outfld('NDROPCOL', ndropcol(:ncol),   ncol, lchnk)
    call outfld('NDROPSRC', nsource(:ncol,:),  ncol, lchnk)
    call outfld('NDROPMIX', ndropmix(:ncol,:), ncol, lchnk)
    call outfld('WTKE    ', wtke(:ncol,:),     ncol, lchnk)

   call ccncalc_oslo(state, pbuf, cs, hasAerosol, numberConcentration, volumeConcentration, &
      hygroscopicity, lnSigma, ccn)
   do isat = 1, psat
      call outfld(ccn_name(isat), ccn(:ncol,:,isat), ncol, lchnk)
   enddo

    tendencyCounted(:)=.FALSE.
    do imode = 1, ntot_amode
       do ispec = 1, nspec_amode(imode)
          mm = mam_idx(imode,ispec)
          itrac = tracer_index(imode,ispec)
          if(.NOT. tendencyCounted(itrac))then
             call outfld(fieldname(mm), coltend(:ncol,itrac), ncol, lchnk)
             call outfld(fieldname_cw(mm), coltend_cw(:ncol,itrac), ncol, lchnk)
             tendencyCounted(itrac)=.TRUE.
          endif
       end do
    end do

    deallocate(nact)
    deallocate(mact)
    deallocate(raer)
    deallocate(qqcw)
    deallocate(raercol)
    deallocate(raercol_cw)
    deallocate(coltend)
    deallocate(coltend_cw)
    deallocate(naermod)
    deallocate(hygro)
    deallocate(lnsigman)  !Variable std. dev (CAM-Oslo)
    deallocate(vaerosol)
    deallocate(fn)
    deallocate(fm)
    deallocate(fluxn)
    deallocate(fluxm)
    deallocate(fluxm_tmp)
    deallocate(fluxn_tmp)
    deallocate(fm_tmp)
    deallocate(fn_tmp)
    deallocate(raercol_tracer)
    deallocate(raercol_cw_tracer)
    deallocate(mact_tracer)
    deallocate(mfullact_tracer)

  end subroutine dropmixnuc_oslo

  !===============================================================================

  subroutine explmix_oslo( q, src, ekkp, ekkm, overlapp, overlapm, &
       qold, surfrate, flxconv, pver, dt, is_unact, qactold )

    ! explicit integration of droplet/aerosol mixing with source due to activation/nucleation

    integer,  intent(in) :: pver           ! number of levels
    real(r8), intent(out):: q(pver)        ! mixing ratio to be updated
    real(r8), intent(in) :: qold(pver)     ! mixing ratio from previous time step
    real(r8), intent(in) :: src(pver)      ! source due to activation/nucleation (/s)
    real(r8), intent(in) :: ekkp(pver)     ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
                                           ! below layer ilev  (ilev,ilev+1 interface)
    real(r8), intent(in) :: ekkm(pver)     ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
                                           ! above layer ilev  (ilev,ilev+1 interface)
    real(r8), intent(in) :: overlapp(pver) ! cloud overlap below
    real(r8), intent(in) :: overlapm(pver) ! cloud overlap above
    real(r8), intent(in) :: surfrate       ! surface exchange rate (/s)
    real(r8), intent(in) :: flxconv        ! convergence of flux from surface
    real(r8), intent(in) :: dt             ! time step (s)
    logical,  intent(in) :: is_unact       ! true if this is an unactivated species
    real(r8), intent(in),optional :: qactold(pver) ! mixing ratio of ACTIVATED species from previous step
                                                   ! *** this should only be present if the current species
                                                   ! is unactivated number/sfc/mass

    integer ilev,kp1,km1

    if ( is_unact ) then
       ! the qactold*(1-overlap) terms are resuspension of activated material
       do ilev=top_lev,pver
          kp1=min(ilev+1,pver)
          km1=max(ilev-1,top_lev)
          q(ilev) = qold(ilev) + dt*(-src(ilev) &
               + ekkp(ilev)*(qold(kp1) - qold(ilev) + qactold(kp1)*(1.0_r8-overlapp(ilev))) &
               + ekkm(ilev)*(qold(km1) - qold(ilev) + qactold(km1)*(1.0_r8-overlapm(ilev))) )
          q(ilev) = max(q(ilev),0._r8)
       end do

       ! diffusion loss at base of lowest layer
       q(pver) = q(pver) - surfrate*qold(pver)*dt + flxconv*dt
       q(pver) = max(q(pver),0._r8)
    else
       do ilev=top_lev,pver
          kp1=min(ilev+1,pver)
          km1=max(ilev-1,top_lev)
          q(ilev) = qold(ilev) + dt*(src(ilev) &
               + ekkp(ilev)*(overlapp(ilev)*qold(kp1)-qold(ilev)) &
               + ekkm(ilev)*(overlapm(ilev)*qold(km1)-qold(ilev)) )
          q(ilev) = max(q(ilev),0._r8) ! force to non-negative if (q(ilev)<-1.e-30) then
       end do

       ! diffusion loss at base of lowest layer
       q(pver) = q(pver) - surfrate*qold(pver)*dt + flxconv*dt ! diffusion loss at base of lowest layer
       q(pver) = max(q(pver),0._r8) ! force to non-negative if(q(pver)<-1.e-30)then
    end if

  end subroutine explmix_oslo

  !===============================================================================

  subroutine activate_modal_oslo(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
       na, nmode, volume, hygro, fn, fm, fluxn, fluxm, flux_fullact, lnsigman )

    ! calculates number, surface, and mass fraction of aerosols activated as CCN
    ! calculates flux of cloud droplets, surface area, and aerosol mass into cloud
    ! assumes an internal mixture within each of up to nmode multiple aerosol modes
    ! a gaussiam spectrum of updrafts can be treated.

    ! mks units

    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.


    ! arguments
    real(r8) , intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
    real(r8) , intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
    real(r8) , intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
    real(r8) , intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
    real(r8) , intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
    real(r8) , intent(in) :: tair          ! air temperature (K)
    real(r8) , intent(in) :: rhoair        ! air density (kg/m3)
    real(r8) , intent(in) :: na(:)         ! aerosol number concentration (/m3)
    integer  , intent(in) :: nmode         ! number of aerosol modes
    real(r8) , intent(in) :: volume(:)     ! aerosol volume concentration (m3/m3)
    real(r8) , intent(in) :: hygro(:)      ! hygroscopicity of aerosol mode
    real(r8) , intent(in) :: lnsigman(:)
    real(r8) , intent(out) :: fn(:)        ! number fraction of aerosols activated
    real(r8) , intent(out) :: fm(:)        ! mass fraction of aerosols activated
    real(r8) , intent(out) :: fluxn(:)     ! flux of activated aerosol number fraction into cloud (cm/s)
    real(r8) , intent(out) :: fluxm(:)     ! flux of activated aerosol mass fraction into cloud (cm/s)
    real(r8) , intent(out) :: flux_fullact ! flux of activated aerosol fraction assuming 100% activation (cm/s)

    ! used for consistency check -- this should match (ekd(ilev)*zs(ilev))
    ! also, fluxm/flux_fullact gives fraction of aerosol mass flux!that is activated

    ! local
    integer, parameter:: nx=200
    integer  :: iquasisect_option, isectional
    real(r8) :: integ,integf
    real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
    real(r8) :: xmin(nmode),xmax(nmode) ! ln(r) at section interfaces
    real(r8) :: volmin(nmode),volmax(nmode) ! volume at interfaces
    real(r8) :: tmass ! total aerosol mass concentration (g/cm3)
    real(r8) :: sign(nmode)    ! geometric standard deviation of size distribution
    real(r8) :: rm ! number mode radius of aerosol at max supersat (cm)
    real(r8) :: pres ! pressure (Pa)
    real(r8) :: path ! mean free path (m)
    real(r8) :: diff ! diffusivity (m2/s)
    real(r8) :: conduct ! thermal conductivity (Joule/m/sec/deg)
    real(r8) :: diff0,conduct0
    real(r8) :: es ! saturation vapor pressure
    real(r8) :: qs ! water vapor saturation mixing ratio
    real(r8) :: dqsdt ! change in qs with temperature
    real(r8) :: dqsdp ! change in qs with pressure
    real(r8) :: g ! thermodynamic function (m2/s)
    real(r8) :: zeta(nmode), eta(nmode)
    real(r8) :: lnsmax ! ln(smax)
    real(r8) :: alpha
    real(r8) :: gamma
    real(r8) :: beta
    real(r8) :: sqrtg
    real(r8) :: amcube(nmode) ! cube of dry mode radius (m)
    real(r8) :: lnsm(nmode) ! ln(smcrit)
    real(r8) :: smc(nmode) ! critical supersaturation for number mode radius
    real(r8) :: sumflx_fullact
    real(r8) :: sumflxn(nmode)
    real(r8) :: sumflxm(nmode)
    real(r8) :: sumfn(nmode)
    real(r8) :: sumfm(nmode)
    real(r8) :: fnold(nmode)   ! number fraction activated
    real(r8) :: fmold(nmode)   ! mass fraction activated
    real(r8) :: exp45logsig_var(nmode)  !variable std. dev (CAM-Oslo)
    real(r8), target :: f1_var(nmode), f2_var(nmode)
    real(r8) :: wold,gold
    real(r8) :: alogam
    real(r8) :: rlo,rhi,xint1,xint2,xint3,xint4
    real(r8) :: wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
    real(r8) :: dfmin,dfmax,fnew,fold,fnmin,fnbar,fsbar,fmbar
    real(r8) :: alw,sqrtalw
    real(r8) :: smax
    real(r8) :: x,arg
    real(r8) :: xmincoeff,xcut,volcut,surfcut
    real(r8) :: z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
    real(r8) :: etafactor1,etafactor2(nmode),etafactor2max
    real(r8) :: grow
    integer  :: imode,n
    character(len=*), parameter :: subname='activate_modal'
    ! numerical integration parameters
    real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8
    real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)
    integer ndist(nx)  ! accumulates frequency distribution of integration bins required
    data ndist/nx*0/
    save ndist

    fn(:)=0._r8
    fm(:)=0._r8
    fluxn(:)=0._r8
    fluxm(:)=0._r8
    flux_fullact=0._r8

    if(nmode == 1.and. na(1) < 1.e-20_r8)return

    if(sigw <= 1.e-5_r8 .and. wbar <= 0._r8)return

    pres = rair*rhoair*tair
    diff0 = 0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
    conduct0 = (5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg

    call qsat(tair, pres, es, qs)

    dqsdt = latvap/(rh2o*tair*tair)*qs
    alpha = gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
    gamma = (1.0_r8+latvap/cpair*dqsdt)/(rhoair*qs)
    etafactor2max = 1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

    grow  =  1._r8/(rhoh2o/(diff0*rhoair*qs) + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))
    sqrtg = sqrt(grow)
    beta  = 2._r8*pi*rhoh2o*grow*gamma

    do imode=1,nmode

       if(volume(imode) > 1.e-39_r8 .and. na(imode) > 1.e-39_r8) then
          ! number mode radius (m)
          exp45logsig_var(imode) = exp(4.5_r8*lnsigman(imode)*lnsigman(imode))
          amcube(imode) = (3._r8*volume(imode)/(4._r8*pi*exp45logsig_var(imode)*na(imode)))  ! only if variable size dist
          f1_var(imode) = 0.5_r8*exp(2.5_r8*lnsigman(imode)*lnsigman(imode))
          f2_var(imode) = 1._r8 + 0.25_r8*lnsigman(imode)

          ! growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
          ! should depend on mean radius of mode to account for gas kinetic effects
          ! see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
          ! for approriate size to use for effective diffusivity.
          etafactor2(imode) = 1._r8/(na(imode)*beta*sqrtg)
          if(hygro(imode) > 1.e-10_r8)then
             smc(imode) = 2._r8*aten*sqrt(aten/(27._r8*hygro(imode)*amcube(imode))) ! only if variable size dist
          else
             smc(imode) = 100._r8
          endif
       else
          smc(imode) = 1._r8
          etafactor2(imode) = etafactor2max ! this should make eta big if na is very small.
       endif
       lnsm(imode) = log(smc(imode)) ! only if variable size dist
    enddo

    if(sigw > 1.e-5_r8)then ! spectrum of updrafts

       wmax = min(wmaxf,wbar+sds*sigw)
       wmin = max(wminf,-wdiab)
       wmin = max(wmin,wbar-sds*sigw)
       w = wmin
       dwmax = eps*sigw
       dw = dwmax
       dfmax = 0.2_r8
       dfmin = 0.1_r8
       if (wmax <= w) return
       do imode=1,nmode
          sumflxn(imode) = 0._r8
          sumfn(imode) = 0._r8
          fnold(imode) = 0._r8
          sumflxm(imode) = 0._r8
          sumfm(imode) = 0._r8
          fmold(imode) = 0._r8
       enddo
       sumflx_fullact = 0._r8

       fold = 0._r8
       wold = 0._r8
       gold = 0._r8

       dwmin = min( dwmax, 0.01_r8 )
       do n = 1, nx

100       wnuc=w+wdiab
          alw=alpha*wnuc
          sqrtalw=sqrt(alw)
          etafactor1=alw*sqrtalw

          do imode=1,nmode
             eta(imode)=etafactor1*etafactor2(imode)
             zeta(imode)=twothird*sqrtalw*aten/sqrtg
          enddo

          call maxsat_oslo(zeta,eta,nmode,smc,smax,f1_var,f2_var)

          lnsmax=log(smax)

          x=twothird*(lnsm(nmode)-lnsmax)/(sq2*lnsigman(nmode))
          fnew=0.5_r8*(1._r8-erf(x))


          dwnew = dw
          if (fnew-fold > dfmax .and. n>1)then
             ! reduce updraft increment for greater accuracy in integration
             if (dw > 1.01_r8*dwmin) then
                dw=0.7_r8*dw
                dw=max(dw,dwmin)
                w=wold+dw
                go to 100
             else
                dwnew = dwmin
             endif
          endif

          if(fnew-fold<dfmin)then
             ! increase updraft increment to accelerate integration
             dwnew=min(1.5_r8*dw,dwmax)
          endif
          fold=fnew

          z=(w-wbar)/(sigw*sq2)
          g=exp(-z*z)
          fnmin=1._r8
          xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

          do imode=1,nmode
              ! modal
             x=twothird*(lnsm(imode)-lnsmax)/(sq2*lnsigman(imode))
             fn(imode)=0.5_r8*(1._r8-erf(x))
             fnmin=min(fn(imode),fnmin)
             ! integration is second order accurate
             ! assumes linear variation of f*g with w
             fnbar=(fn(imode)*g+fnold(imode)*gold)
             arg=x-1.5_r8*sq2*lnsigman(imode)
             fm(imode)=0.5_r8*(1._r8-erf(arg))
             fmbar=(fm(imode)*g+fmold(imode)*gold)
             wb=(w+wold)
             if(w > 0._r8)then
                sumflxn(imode)=sumflxn(imode)+sixth*(wb*fnbar           &
                     +(fn(imode)*g*w+fnold(imode)*gold*wold))*dw
                sumflxm(imode)=sumflxm(imode)+sixth*(wb*fmbar           &
                     +(fm(imode)*g*w+fmold(imode)*gold*wold))*dw
             endif
             sumfn(imode)=sumfn(imode)+0.5_r8*fnbar*dw
             fnold(imode)=fn(imode)
             sumfm(imode)=sumfm(imode)+0.5_r8*fmbar*dw
             fmold(imode)=fm(imode)
          enddo
          ! same form as sumflxm but replace the fm with 1.0
          sumflx_fullact = sumflx_fullact &
               + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
          ! sumg=sumg+0.5_r8*(g+gold)*dw
          gold=g
          wold=w
          dw=dwnew
          if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit
          w=w+dw
          if (n == nx) then
             write(iulog,*)'do loop is too short in activate'
             write(iulog,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
             write(iulog,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
             write(iulog,*)'wnuc=',wnuc
             write(iulog,*)'na=',(na(imode),imode=1,nmode)
             write(iulog,*)'fn=',(fn(imode),imode=1,nmode)
             ! dump all subr parameters to allow testing with standalone code
             ! (build a driver that will read input and call activate)
             write(iulog,*)'wbar,sigw,wdiab,tair,rhoair,nmode='
             write(iulog,*) wbar,sigw,wdiab,tair,rhoair,nmode
             write(iulog,*)'na=',na
             write(iulog,*)'volume=', (volume(imode),imode=1,nmode)
             write(iulog,*)'hydro='
             write(iulog,*) hygro
             call endrun(subname)
          end if

       enddo

       ndist(n)=ndist(n)+1
       if(w<wmaxf)then

          ! contribution from all updrafts stronger than wmax
          ! assuming constant f (close to fmax)
          wnuc=w+wdiab

          z1=(w-wbar)/(sigw*sq2)
          z2=(wmaxf-wbar)/(sigw*sq2)
          g=exp(-z1*z1)
          integ=sigw*0.5_r8*sq2*sqpi*(erf(z2)-erf(z1))
          ! consider only upward flow into cloud base when estimating flux
          wf1=max(w,zero)
          zf1=(wf1-wbar)/(sigw*sq2)
          gf1=exp(-zf1*zf1)
          wf2=max(wmaxf,zero)
          zf2=(wf2-wbar)/(sigw*sq2)
          gf2=exp(-zf2*zf2)
          gf=(gf1-gf2)
          integf=wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

          do imode=1,nmode
             sumflxn(imode)=sumflxn(imode)+integf*fn(imode)
             sumfn(imode)=sumfn(imode)+fn(imode)*integ
             sumflxm(imode)=sumflxm(imode)+integf*fm(imode)
             sumfm(imode)=sumfm(imode)+fm(imode)*integ
          enddo
          ! same form as sumflxm but replace the fm with 1.0
          sumflx_fullact = sumflx_fullact + integf
          ! sumg=sumg+integ
       endif


       do imode=1,nmode
          fn(imode)=sumfn(imode)/(sq2*sqpi*sigw)
          ! fn(imode)=sumfn(imode)/(sumg)
          if(fn(imode) > 1.01_r8)then
             write(iulog,*)'fn=',fn(imode),' > 1 in activate'
             write(iulog,*)'w,imode,na,amcube=',w,imode,na(imode),amcube(imode)
             write(iulog,*)'integ,sumfn,sigw=',integ,sumfn(imode),sigw
             call endrun('activate')
          endif
          fluxn(imode)=sumflxn(imode)/(sq2*sqpi*sigw)
          fm(imode)=sumfm(imode)/(sq2*sqpi*sigw)
          ! fm(imode)=sumfm(imode)/(sumg)
          if(fm(imode) > 1.01_r8)then
             write(iulog,*)'fm=',fm(imode),' > 1 in activate'
          endif
          fluxm(imode)=sumflxm(imode)/(sq2*sqpi*sigw)
       enddo
       ! same form as fluxm
       flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

    else

       ! single updraft
       wnuc=wbar+wdiab

       if(wnuc > 0._r8)then
          w=wbar
          alw=alpha*wnuc
          sqrtalw=sqrt(alw)
          etafactor1=alw*sqrtalw

          do imode = 1,nmode
             eta(imode) = etafactor1*etafactor2(imode)
             zeta(imode) = twothird*sqrtalw*aten/sqrtg
             f1_var(imode) = 0.5_r8*exp(2.5_r8*lnsigman(imode)*lnsigman(imode))
             f2_var(imode) = 1._r8 + 0.25_r8*lnsigman(imode)
          enddo

          call maxsat_oslo(zeta,eta,nmode,smc,smax,f1_var, f2_var)

          lnsmax=log(smax)
          xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3
          do imode = 1,nmode
             x = twothird*(lnsm(imode)-lnsmax)/(sq2*lnsigman(imode))
             fn(imode) = 0.5_r8*(1._r8-erf(x))
             arg = x-1.5_r8*sq2*lnsigman(imode)
             fm(imode) = 0.5_r8*(1._r8-erf(arg))
             if (wbar > 0._r8)then
                fluxn(imode) = fn(imode)*w
                fluxm(imode) = fm(imode)*w
             endif
          enddo
          flux_fullact = w
       endif

    endif

  end subroutine activate_modal_oslo

  !===============================================================================
  subroutine maxsat_oslo(zeta, eta, nmode, smc, smax, f1_in, f2_in)

    ! calculates maximum supersaturation for multiple competing aerosol modes.
    ! Abdul-Razzak and Ghan, A parameterization of aerosol activation.
    ! 2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

    ! arguments
    real(r8), intent(in)  :: zeta(nmode)
    real(r8), intent(in)  :: eta(nmode)
    integer,  intent(in)  :: nmode ! number of modes
    real(r8), intent(in)  :: smc(nmode) ! critical supersaturation for number mode radius
    real(r8), intent(in), target :: f1_in(:)
    real(r8), intent(in), target :: f2_in(:)
    real(r8), intent(out) :: smax ! maximum supersaturation

    ! local variables
    integer  :: imode  ! mode index
    real(r8) :: sum, g1, g2, g1sqrt, g2sqrt
    real(r8), pointer :: f1_used(:), f2_used(:)

    f1_used => f1_in
    f2_used => f2_in

    do imode=1,nmode
       if(zeta(imode) > 1.e5_r8*eta(imode) .or. smc(imode)*smc(imode) > 1.e5_r8*eta(imode))then
          ! weak forcing. essentially none activated
          smax=1.e-20_r8
       else
          ! significant activation of this mode. calc activation all modes.
          exit
       endif
       ! No significant activation in any mode.  Do nothing.
       if (imode == nmode) return
    enddo

    sum = 0.0_r8
    do imode = 1,nmode
       if(eta(imode) > 1.e-20_r8)then
          g1 = zeta(imode)/eta(imode)
          g1sqrt = sqrt(g1)
          g1 = g1sqrt*g1
          g2 = smc(imode)/sqrt(eta(imode)+3._r8*zeta(imode))
          g2sqrt = sqrt(g2)
          g2 = g2sqrt*g2
          sum = sum+(f1_used(imode)*g1+f2_used(imode)*g2)/(smc(imode)*smc(imode))
       else
          sum = 1.e20_r8
       endif
    enddo
    smax = 1._r8/sqrt(sum)

  end subroutine maxsat_oslo

  !===============================================================================

  subroutine ccncalc_oslo(state, pbuf, cs, hasAerosol, numberConcentration, volumeConcentration, &
       hygroscopicity, lnSigma, ccn)

    ! calculates number concentration of aerosols activated as CCN at
    ! supersaturation supersat.
    ! assumes an internal mixture of a multiple externally-mixed aerosol modes cgs units

    ! This was used in the BACCHUS-project where it was agreed that
    ! CCN would not include cloud-borne aerosols. It is possible to
    ! calculate cloud-borne aerosols, but it is complicated, and it was
    ! not needed when this code was made.

    ! arguments
    type(physics_state), target, intent(in)    :: state
    type(physics_buffer_desc),   pointer       :: pbuf(:)
    real(r8) , intent(in)  :: cs(pcols,pver)                           ! air density (kg/m3)
    logical  , intent(in)  :: hasAerosol(pcols, pver, nmodes)
    real(r8) , intent(in)  :: numberConcentration(pcols,pver,0:nmodes) ! interstit+activated aerosol number conc (/m3)
    real(r8) , intent(in)  :: volumeConcentration(pcols,pver,nmodes)   ! interstit+activated aerosol volume conc (m3/m3)
    real(r8) , intent(in)  :: hygroscopicity(pcols,pver,nmodes)
    real(r8) , intent(in)  :: lnSigma(pcols,pver,nmodes)
    real(r8) , intent(out) :: ccn(pcols,pver,psat)                     ! number conc of aerosols activated at supersat (#/m3)

    ! local
    integer  :: lchnk             ! chunk index
    integer  :: ncol              ! number of columns
    real(r8) :: super(psat)       ! supersaturation
    real(r8) :: surften_coef      ! Coefficient in ARGI / ARGII
    real(r8) :: amcube            ! number median radius qubed
    real(r8) :: a                 ! surface tension parameter
    real(r8) :: sm                ! critical supersaturation at mode radius
    real(r8) :: arg               ! factor in eqn 15 ARGII
    real(r8) :: argfactor         ! Coefficient in ARGI/ARGII
    real(r8) :: exp45logsig_var   ! mathematical constants
    integer  :: lsat,imode,icol,ilev        ! mathematical constants
    real(r8) :: smcoefcoef,smcoef ! mathematical constants
    real(r8), pointer   :: tair(:,:)        ! air temperature (K)
    real(r8), parameter :: twothird=2.0_r8/3.0_r8
    real(r8), parameter :: sq2=sqrt(2.0_r8)
    real(r8), parameter :: surften=0.076_r8 !surface tension of water (J/m2)
    !-------------------------------------------------------------------------------

    lchnk = state%lchnk
    ncol  = state%ncol
    tair  => state%t

    super(:) = supersat(:)*0.01_r8

    !This is curvature effect (A) in ARGI eqn 5 in ARG1 (missing division by temperature, see below)
    surften_coef = 2._r8*mwh2o*surften/(r_universal*rhoh2o)

    !This is part of eqn 9 in ARGII where A smcoefcoef is 2/3^(3/2)
    smcoefcoef = 2._r8/sqrt(27._r8)

    ccn(:,:,:) = 0._r8

    do imode=1,nmodes
       do ilev=top_lev,pver
          do icol=1,ncol
             if (hasAerosol(icol,ilev,imode)) then

                !Curvature-parameter "A" in ARGI (eqn 5)
                a = surften_coef/tair(icol,ilev)

                !standard factor for transforming size distr, volume ==> number (google psd.pdf by zender)
                exp45logsig_var = exp(4.5_r8*lnsigma(icol,ilev,imode)*lnsigma(icol,ilev,imode))

                ! Numbe rmedian radius (power of three)
                ! By definition of lognormal distribution only if variable size dist
                amcube =(3._r8*volumeConcentration(icol,ilev,imode) /(4._r8*pi*exp45logsig_var*numberConcentration(icol,ilev,imode)))

                !This is part of eqn 9 in ARGII where A smcoefcoef is 2/3^(3/2)
                smcoef = smcoefcoef * a * sqrt(a)

                !This is finally solving eqn 9 (solve for critical supersat of mode)
                sm = smcoef / sqrt(hygroscopicity(icol,ilev,imode)*amcube) ! critical supersaturation

                !Solve eqn 13 in ARGII
                do lsat = 1,psat

                   !eqn 15 in ARGII
                   argfactor = twothird/(sq2*lnSigma(icol,ilev,imode))

                   !eqn 15 in ARGII
                   arg = argfactor*log(sm/super(lsat))

                   !eqn 13 icol ARGII
                   ccn(icol,ilev,lsat) = ccn(icol,ilev,lsat) + numberConcentration(icol,ilev,imode)*0.5_r8*(1._r8-erf(arg))

                end do
             end if
          end do
       end do
    end do

    ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3

  end subroutine ccncalc_oslo

end module oslo_aero_ndrop
