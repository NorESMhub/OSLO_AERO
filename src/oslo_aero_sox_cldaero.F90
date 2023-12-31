module oslo_aero_sox_cldaero

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver
  use cam_abortutils,  only: endrun
  use mo_chem_utls,    only: get_spc_ndx
  use mo_constants,    only: pi
  use chem_mods,       only: adv_mass
  use physconst,       only: gravit
  use chem_mods,       only: gas_pcnst
  !
  use oslo_aero_share, only: l_so4_a2, chemistryIndex

  implicit none
  private

  public :: sox_cldaero_init
  public :: sox_cldaero_create_obj
  public :: sox_cldaero_update
  public :: sox_cldaero_destroy_obj

  private :: cldaero_uptakerate
  private :: cldaero_allocate
  private :: cldaero_deallocate

  type :: cldaero_conc_t
     real(r8), pointer :: so4c(:,:)
     real(r8), pointer :: nh4c(:,:)
     real(r8), pointer :: no3c(:,:)
     real(r8), pointer :: xlwc(:,:)
     real(r8) :: so4_fact
  end type cldaero_conc_t
  public :: cldaero_conc_t

  integer :: id_msa, id_h2so4, id_so2, id_h2o2, id_nh3
  integer :: id_so4_1a

  real(r8), parameter :: small_value = 1.e-20_r8

!===============================================================================
contains
!===============================================================================

  subroutine sox_cldaero_init()

    ! module variables
    id_msa   = get_spc_ndx( 'MSA' )
    id_h2so4 = get_spc_ndx( 'H2SO4' )
    id_so2   = get_spc_ndx( 'SO2' )
    id_h2o2  = get_spc_ndx( 'H2O2' )
    id_nh3   = get_spc_ndx( 'NH3' )

    if (id_h2so4<1 .or. id_so2<1 .or. id_h2o2<1) then
       call endrun('sox_cldaero_init: oslo aero does not include necessary species' &
            //' -- should not invoke sox_cldaero_init ')
    endif

    id_so4_1a = chemistryIndex(l_so4_a2) 

  end subroutine sox_cldaero_init

  !===============================================================================

  function sox_cldaero_create_obj(cldfrc, qcw, lwc, cfact, ncol, loffset) result( conc_obj )

    ! arguments
    real(r8), intent(in) :: cldfrc(:,:)
    real(r8), intent(in) :: qcw(:,:,:)
    real(r8), intent(in) :: lwc(:,:)
    real(r8), intent(in) :: cfact(:,:)
    integer,  intent(in) :: ncol
    integer,  intent(in) :: loffset
    type(cldaero_conc_t), pointer :: conc_obj

    ! local variables
    integer :: l,n
    integer :: i,k

    conc_obj => cldaero_allocate()

    do k = 1,pver
       do i = 1,ncol
          if (cldfrc(i,k) >0._r8) then
             conc_obj%xlwc(i,k) = lwc(i,k) *cfact(i,k) ! cloud water L(water)/L(air)
             conc_obj%xlwc(i,k) = conc_obj%xlwc(i,k) / cldfrc(i,k) ! liquid water in the cloudy fraction of cell
          else
             conc_obj%xlwc(i,k) = 0._r8
          endif
       enddo
    enddo

    conc_obj%no3c(:,:) = 0._r8

    ! Set concenctration of cloud so4
    conc_obj%so4c(:ncol,:) = qcw(:ncol,:,id_so4_1a) 

    ! current version does not have nh3/nh4 tracers -  so so4 is assumed to be nh4hso4
    ! the partial neutralization of so4 is handled by using a 
    ! -1 charge (instead of -2) in the electro-neutrality equation
    conc_obj%nh4c(:ncol,:) = 0._r8

    ! with 3-mode, assume so4 is nh4hso4, and so half-neutralized
    conc_obj%so4_fact = 1._r8

  end function sox_cldaero_create_obj

  !===============================================================================

  subroutine sox_cldaero_update( &
       ncol, lchnk, loffset, dtime, mbar, pdel, press, tfld, cldnum, cldfrc, cfact, xlwc, &
       delso4_hprxn, xh2so4, xso4, xso4_init, nh3g, hno3g, xnh3, xhno3, xnh4c,  xno3c, xmsa, xso2, xh2o2, qcw, qin, &
       aqso4, aqh2so4, aqso4_h2o2, aqso4_o3, aqso4_h2o2_3d, aqso4_o3_3d)

    !----------------------------------------------------------------------------------
    ! Update the mixing ratios
    !----------------------------------------------------------------------------------

    ! arguments
    integer,  intent(in)    :: ncol
    integer,  intent(in)    :: lchnk         ! chunk id
    integer,  intent(in)    :: loffset
    real(r8), intent(in)    :: dtime         ! time step (sec)
    real(r8), intent(in)    :: mbar(:,:)     ! mean wet atmospheric mass ( amu )
    real(r8), intent(in)    :: pdel(:,:) 
    real(r8), intent(in)    :: press(:,:)
    real(r8), intent(in)    :: tfld(:,:)
    real(r8), intent(in)    :: cldnum(:,:)
    real(r8), intent(in)    :: cldfrc(:,:)
    real(r8), intent(in)    :: cfact(:,:)
    real(r8), intent(in)    :: xlwc(:,:)
    real(r8), intent(in)    :: delso4_hprxn(:,:)
    real(r8), intent(in)    :: xh2so4(:,:)
    real(r8), intent(in)    :: xso4(:,:)
    real(r8), intent(in)    :: xso4_init(:,:)
    real(r8), intent(in)    :: nh3g(:,:)
    real(r8), intent(in)    :: hno3g(:,:)
    real(r8), intent(in)    :: xnh3(:,:)
    real(r8), intent(in)    :: xhno3(:,:)
    real(r8), intent(in)    :: xnh4c(:,:)
    real(r8), intent(in)    :: xmsa(:,:)
    real(r8), intent(in)    :: xso2(:,:)
    real(r8), intent(in)    :: xh2o2(:,:)
    real(r8), intent(in)    :: xno3c(:,:)
    real(r8), intent(inout) :: qcw(:,:,:)                 ! cloud-borne aerosol (vmr)
    real(r8), intent(inout) :: qin(:,:,:)                 ! xported species ( vmr )
    real(r8), intent(out)   :: aqso4(:,:)                 ! aqueous phase chemistry
    real(r8), intent(out)   :: aqh2so4(:,:)               ! aqueous phase chemistry
    real(r8), intent(out)   :: aqso4_h2o2(:)              ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out)   :: aqso4_o3(:)                ! SO4 aqueous phase chemistry due to O3 (kg/m2)
    real(r8), intent(out), optional :: aqso4_h2o2_3d(:,:) ! SO4 aqueous phase chemistry due to H2O2 (kg/m2)
    real(r8), intent(out), optional :: aqso4_o3_3d(:,:)   ! SO4 aqueous phase chemistry due to O3 (kg/m2)


    ! local variables
    real(r8) :: dqdt_aqso4(ncol,pver,gas_pcnst)
    real(r8) :: dqdt_aqh2so4(ncol,pver,gas_pcnst)
    real(r8) :: dqdt_aqhprxn(ncol,pver), dqdt_aqo3rxn(ncol,pver)
    real(r8) :: sflx(1:ncol)
    real(r8) :: delso4_o3rxn
    real(r8) :: dso4dt_aqrxn, dso4dt_hprxn
    real(r8) :: dso4dt_gasuptk, dmsadt_gasuptk
    real(r8) :: dmsadt_gasuptk_tomsa, dmsadt_gasuptk_toso4
    real(r8) :: dqdt_aq, dqdt_wr, dqdt
    real(r8) :: fwetrem, sumf, uptkrate
    real(r8) :: delnh3, delnh4
    integer  :: l, n, m, i,k
    integer  :: ntot_msa_c
    real(r8) :: xl

    ! make sure dqdt is zero initially, for budgets
    dqdt_aqso4(:,:,:) = 0.0_r8
    dqdt_aqh2so4(:,:,:) = 0.0_r8
    dqdt_aqhprxn(:,:) = 0.0_r8
    dqdt_aqo3rxn(:,:) = 0.0_r8

    lev_loop: do k = 1,pver
       col_loop: do i = 1,ncol
          cloud: if (cldfrc(i,k) >= 1.0e-5_r8) then
             xl = xlwc(i,k) ! / cldfrc(i,k)

             IF (XL .ge. 1.e-8_r8) THEN !! WHEN CLOUD IS PRESENTED

                delso4_o3rxn = xso4(i,k) - xso4_init(i,k)

                if (id_nh3>0) then
                   delnh3 = nh3g(i,k) - xnh3(i,k)
                   delnh4 = - delnh3
                endif

                !In the case of OSLO-AEROSOLS, 
                !set no MSA in cloud droplets
                ntot_msa_c = 0

                !   average uptake rate over dtime
                uptkrate = cldaero_uptakerate( xl, cldnum(i,k), cfact(i,k), cldfrc(i,k), tfld(i,k),  press(i,k) )

                ! average uptake rate over dtime
                uptkrate = (1.0_r8 - exp(-min(100._r8,dtime*uptkrate))) / dtime

                !   dso4dt_gasuptk = so4_c tendency from h2so4 gas uptake (mol/mol/s)
                !   dmsadt_gasuptk = msa_c tendency from msa   gas uptake (mol/mol/s)
                dso4dt_gasuptk = xh2so4(i,k) * uptkrate
                if (id_msa > 0) then
                   dmsadt_gasuptk = xmsa(i,k) * uptkrate
                else
                   dmsadt_gasuptk = 0.0_r8
                end if

                !   if no modes have msa aerosol, then "rename" scavenged msa gas to so4
                dmsadt_gasuptk_toso4 = 0.0_r8
                dmsadt_gasuptk_tomsa = dmsadt_gasuptk
                if (ntot_msa_c == 0) then
                   dmsadt_gasuptk_tomsa = 0.0_r8
                   dmsadt_gasuptk_toso4 = dmsadt_gasuptk
                end if

                !-----------------------------------------------------------------------
                !      now compute TMR tendencies
                !      this includes the above aqueous so2 chemistry AND
                !      the uptake of highly soluble aerosol precursor gases (h2so4, msa, ...)
                !      AND the wetremoval of dissolved, unreacted so2 and h2o2

                dso4dt_aqrxn = (delso4_o3rxn + delso4_hprxn(i,k)) / dtime
                dso4dt_hprxn = delso4_hprxn(i,k) / dtime

                !   fwetrem = fraction of in-cloud-water material that is wet removed
                !       fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*clwlrat(i,k)))) )
                fwetrem = 0.0_r8    ! don't have so4 & msa wet removal here

                !Update so4 in cloud water
                l = id_so4_1a       !We only have one aq-phase tracer in CAM_OSLO

                dqdt_aqso4(i,k,l) = dso4dt_aqrxn*cldfrc(i,k)
                dqdt_aqh2so4(i,k,l) = (dso4dt_gasuptk + dmsadt_gasuptk_toso4)*cldfrc(i,k)
                dqdt_aq = dqdt_aqso4(i,k,l) + dqdt_aqh2so4(i,k,l)
                dqdt_wr = -fwetrem*dqdt_aq  !wet removal set to zero above
                dqdt= dqdt_aq + dqdt_wr
                qcw(i,k,l) = qcw(i,k,l) + dqdt*dtime

                !Additional updates for MSA??
                !      For gas species, tendency includes
                !      reactive uptake to cloud water that essentially transforms the gas to
                !      a different species.  Wet removal associated with this is applied
                !      to the "new" species (e.g., so4_c) rather than to the gas.
                !      wet removal of the unreacted gas that is dissolved in cloud water.
                !      Need to multiply both these parts by cldfrc

                !   h2so4 (g) & msa (g)
                qin(i,k,id_h2so4) = qin(i,k,id_h2so4) - dso4dt_gasuptk * dtime * cldfrc(i,k)
                if (id_msa > 0) qin(i,k,id_msa) = qin(i,k,id_msa) - dmsadt_gasuptk * dtime * cldfrc(i,k)


                ! so2 -- the first order loss rate for so2 is frso2_c*clwlrat(i,k)
                ! fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frso2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8   ! don't include so2 wet removal here

                dqdt_wr = -fwetrem*xso2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_aqrxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_so2) = qin(i,k,id_so2) + dqdt * dtime

                !   h2o2 -- the first order loss rate for h2o2 is frh2o2_c*clwlrat(i,k)
                !       fwetrem = max( 0.0_r8, (1.0_r8-exp(-min(100._r8,dtime*frh2o2_c*clwlrat(i,k)))) )
                fwetrem = 0.0_r8   ! don't include h2o2 wet removal here

                dqdt_wr = -fwetrem*xh2o2(i,k)/dtime*cldfrc(i,k)
                dqdt_aq = -dso4dt_hprxn*cldfrc(i,k)
                dqdt = dqdt_aq + dqdt_wr
                qin(i,k,id_h2o2) = qin(i,k,id_h2o2) + dqdt * dtime

                ! NH3
                if (id_nh3>0) then
                   dqdt_aq = delnh3/dtime*cldfrc(i,k)
                   dqdt = dqdt_aq
                   qin(i,k,id_nh3) = qin(i,k,id_nh3) + dqdt * dtime
                endif

                !   for SO4 from H2O2/O3 budgets
                dqdt_aqhprxn(i,k) = dso4dt_hprxn*cldfrc(i,k)
                dqdt_aqo3rxn(i,k) = (dso4dt_aqrxn - dso4dt_hprxn)*cldfrc(i,k)

             ENDIF !! WHEN CLOUD IS PRESENTED

          endif cloud
       enddo col_loop
    enddo lev_loop

    ! Update the mixing ratios
    do k = 1,pver
       qcw(:,k,id_so4_1a) =  MAX( qcw(:,k,id_so4_1a), small_value )
       qin(:,k,id_so2) =  MAX( qin(:,k,id_so2), small_value )
       if ( id_nh3 > 0 ) then
          qin(:,k,id_nh3) =  MAX( qin(:,k,id_nh3), small_value )
       endif
    end do

    ! diagnostics
    l = id_so4_1a !Index of the a2-tracer in cloud water
    n = 1         !Only distribute to one "mode" 
    aqso4(:,n)=0._r8
    do k=1,pver
       do i=1,ncol
          aqso4(i,n)=aqso4(i,n)+dqdt_aqso4(i,k,l)*adv_mass(l)/mbar(i,k) &
               *pdel(i,k)/gravit ! kg/m2/s
       enddo
    enddo

    aqh2so4(:,n)=0._r8
    do k=1,pver
       do i=1,ncol
          aqh2so4(:,n)=aqh2so4(:,n)+dqdt_aqh2so4(i,k,l)*adv_mass(l)/mbar(i,k) &
               *pdel(i,k)/gravit ! kg/m2/s
       enddo
    enddo

    aqso4_h2o2(:) = 0._r8
    do k=1,pver
       do i=1,ncol
          aqso4_h2o2(i)=aqso4_h2o2(i)+dqdt_aqhprxn(i,k)*adv_mass(l)/mbar(i,k) &
               *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo

    if (present(aqso4_h2o2_3d)) then 
       aqso4_h2o2_3d(:,:) = 0._r8
       do k=1,pver
          do i=1,ncol
             aqso4_h2o2_3d(i,k)=dqdt_aqhprxn(i,k)*adv_mass(l)/mbar(i,k) &
                  *pdel(i,k)/gravit ! kg SO4 /m2/s
          enddo
       enddo
    end if

    aqso4_o3(:)=0._r8
    do k=1,pver
       do i=1,ncol
          aqso4_o3(i)=aqso4_o3(i)+dqdt_aqo3rxn(i,k)*adv_mass(l)/mbar(i,k) &
               *pdel(i,k)/gravit ! kg SO4 /m2/s
       enddo
    enddo

    if (present(aqso4_o3_3d)) then
       aqso4_o3_3d(:,:)=0._r8
       do k=1,pver
          do i=1,ncol
             aqso4_o3_3d(i,k)=dqdt_aqo3rxn(i,k)*adv_mass(l)/mbar(i,k) &
                  *pdel(i,k)/gravit ! kg SO4 /m2/s
          enddo
       enddo
    end if

  end subroutine sox_cldaero_update

  !===============================================================================

  subroutine sox_cldaero_destroy_obj( conc_obj )
    type(cldaero_conc_t), pointer :: conc_obj
    call cldaero_deallocate( conc_obj )
  end subroutine sox_cldaero_destroy_obj

  !===============================================================================

  function cldaero_allocate( ) result( cldconc )
    type(cldaero_conc_t), pointer:: cldconc

    allocate( cldconc )
    allocate( cldconc%so4c(pcols,pver) )
    allocate( cldconc%nh4c(pcols,pver) )
    allocate( cldconc%no3c(pcols,pver) )
    allocate( cldconc%xlwc(pcols,pver) )

    cldconc%so4c(:,:) = 0._r8
    cldconc%nh4c(:,:) = 0._r8
    cldconc%no3c(:,:) = 0._r8
    cldconc%xlwc(:,:) = 0._r8
    cldconc%so4_fact  = 2._r8

  end function cldaero_allocate

  !===============================================================================

  subroutine cldaero_deallocate( cldconc )
    type(cldaero_conc_t), pointer :: cldconc

    if ( associated(cldconc%so4c) ) then
       deallocate(cldconc%so4c)
       nullify(cldconc%so4c)
    endif

    if ( associated(cldconc%nh4c) ) then
       deallocate(cldconc%nh4c)
       nullify(cldconc%nh4c)
    endif

    if ( associated(cldconc%no3c) ) then
       deallocate(cldconc%no3c)
       nullify(cldconc%no3c)
    endif

    if ( associated(cldconc%xlwc) ) then
       deallocate(cldconc%xlwc)
       nullify(cldconc%xlwc)
    endif

    deallocate( cldconc )
    nullify( cldconc )

  end subroutine cldaero_deallocate

  !===============================================================================

  function cldaero_uptakerate( xl, cldnum, cfact, cldfrc, tfld,  press ) result( uptkrate )

    ! compute uptake of h2so4 and msa to cloud water
    ! first-order uptake rate is
    ! 4*pi*(drop radius)*(drop number conc)
    ! *(gas diffusivity)*(fuchs sutugin correction)

    ! arguments / output
    real(r8), intent(in) :: xl, cldnum, cfact, cldfrc, tfld,  press
    real(r8) :: uptkrate

    ! local variables
    real(r8) :: rad_cd, radxnum_cd, num_cd
    real(r8) :: gasdiffus, gasspeed, knudsen
    real(r8) :: fuchs_sutugin, volx34pi_cd

    ! num_cd = (drop number conc in 1/cm^3)
    num_cd = 1.0e-3_r8*cldnum*cfact/cldfrc
    num_cd = max( num_cd, 0.0_r8 )

    ! rad_cd = (drop radius in cm), computed from liquid water and drop number,
    ! then bounded by 0.5 and 50.0 micrometers
    ! radxnum_cd = (drop radius)*(drop number conc)
    ! volx34pi_cd = (3/4*pi) * (liquid water volume in cm^3/cm^3)

    volx34pi_cd = xl*0.75_r8/pi

    ! following holds because volx34pi_cd = num_cd*(rad_cd**3)
    radxnum_cd = (volx34pi_cd*num_cd*num_cd)**0.3333333_r8

    ! apply bounds to rad_cd to avoid the occasional unphysical value
    if (radxnum_cd .le. volx34pi_cd*4.0e4_r8) then
       radxnum_cd = volx34pi_cd*4.0e4_r8
       rad_cd = 50.0e-4_r8
    else if (radxnum_cd .ge. volx34pi_cd*4.0e8_r8) then
       radxnum_cd = volx34pi_cd*4.0e8_r8
       rad_cd = 0.5e-4_r8
    else
       rad_cd = radxnum_cd/num_cd
    end if

    ! gasdiffus = h2so4 gas diffusivity from mosaic code (cm^2/s) (pmid must be Pa)
    gasdiffus = 0.557_r8 * (tfld**1.75_r8) / press

    ! gasspeed = h2so4 gas mean molecular speed from mosaic code (cm/s)
    gasspeed = 1.455e4_r8 * sqrt(tfld/98.0_r8)

    ! knudsen number
    knudsen = 3.0_r8*gasdiffus/(gasspeed*rad_cd)

    ! following assumes accomodation coefficient = 0.65
    ! (Adams & Seinfeld, 2002, JGR, and references therein)
    ! fuchs_sutugin = (0.75*accom*(1. + knudsen)) /
    ! (knudsen*(1.0 + knudsen + 0.283*accom) + 0.75*accom)
    fuchs_sutugin = (0.4875_r8*(1._r8 + knudsen)) / (knudsen*(1.184_r8 + knudsen) + 0.4875_r8)

    ! instantaneous uptake rate
    uptkrate = 12.56637_r8*radxnum_cd*gasdiffus*fuchs_sutugin

  end function cldaero_uptakerate

end module oslo_aero_sox_cldaero
