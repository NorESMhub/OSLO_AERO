module oslo_aero_optical_params

  ! Optical parameters for a composite aerosol is calculated by interpolation
  ! from the tables kcomp1.out-kcomp14.out.

  use shr_kind_mod,    only: r8 => shr_kind_r8
  use ppgrid,          only: pcols, pver, pverp
  use constituents,    only: pcnst
  use cam_history,     only: outfld
  use physconst,       only: rair,pi
  use physics_types,   only: physics_state
  use wv_saturation,   only: qsat_water
  !
  use oslo_aero_utils, only: calculateNumberConcentration
  use oslo_aero_conc,  only: calculateBulkProperties, partitionMass
  use oslo_aero_sw_tables
  use oslo_aero_params
  use oslo_aero_const
  use oslo_aero_share

  implicit none
  private

  public :: oslo_aero_optical_params_calc

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_optical_params_calc(lchnk, ncol, pint, pmid,                &
       coszrs, state, t, cld, qm1, Nnatk,                                          &
       per_tau, per_tau_w, per_tau_w_g, per_tau_w_f, per_lw_abs,                   &
       volc_ext_sun, volc_omega_sun, volc_g_sun, volc_ext_earth, volc_omega_earth, &
       aodvis, absvis)

    ! Input arguments
    integer , intent(in) :: lchnk                                 ! chunk identifier
    integer , intent(in) :: ncol                                  ! number of atmospheric columns
    real(r8), intent(in) :: coszrs(pcols)                         ! Cosine solar zenith angle
    real(r8), intent(in) :: pint(pcols,pverp)                     ! Model interface pressures (10*Pa)
    real(r8), intent(in) :: pmid(pcols,pver)                      ! Model level pressures (Pa)
    real(r8), intent(in) :: t(pcols,pver)                         ! Model level temperatures (K)
    real(r8), intent(in) :: cld(pcols,pver)                       ! cloud fraction
    real(r8), intent(in) :: qm1(pcols,pver,pcnst)                 ! Specific humidity and tracers (kg/kg)
    real(r8), intent(in) :: volc_ext_sun(pcols,pver,nbands)       ! volcanic aerosol extinction for solar bands, CMIP6
    real(r8), intent(in) :: volc_omega_sun(pcols,pver,nbands)     ! volcanic aerosol SSA for solar bands, CMIP6
    real(r8), intent(in) :: volc_g_sun(pcols,pver,nbands)         ! volcanic aerosol g for solar bands, CMIP6
    real(r8), intent(in) :: volc_ext_earth(pcols,pver,nlwbands)   ! volcanic aerosol extinction for terrestrial bands, CMIP6
    real(r8), intent(in) :: volc_omega_earth(pcols,pver,nlwbands) ! volcanic aerosol SSA for terrestrial bands, CMIP6
    type(physics_state), intent(in), target :: state

    ! Input-output arguments
    real(r8), intent(inout) :: Nnatk(pcols,pver,0:nmodes)         ! aerosol mode number concentration

    ! Output arguments
    ! AOD and absorptive AOD for visible wavelength closest to 0.55 um (0.442-0.625)
    ! Note that aodvis and absvis output should be devided by dayfoc to give physical (A)AOD values
    real(r8), intent(out) :: per_tau    (pcols,0:pver,nbands)     ! aerosol extinction optical depth
    real(r8), intent(out) :: per_tau_w  (pcols,0:pver,nbands)     ! aerosol single scattering albedo * tau
    real(r8), intent(out) :: per_tau_w_g(pcols,0:pver,nbands)     ! aerosol assymetry parameter * w * tau
    real(r8), intent(out) :: per_tau_w_f(pcols,0:pver,nbands)     ! aerosol forward scattered fraction * w * tau
    real(r8), intent(out) :: per_lw_abs (pcols,pver,nlwbands)     ! aerosol absorption optical depth (LW)
    real(r8), intent(out) :: aodvis(pcols)                        ! AOD vis
    real(r8), intent(out) :: absvis(pcols)                        ! AAOD vis

    ! Local variables
    integer  :: i, k, ib, icol, mplus10
    integer  :: iloop
    logical  :: daylight(pcols)        ! SW calculations also at (polar) night in interpol* if daylight=.true.
    real(r8) :: aodvisvolc(pcols)      ! AOD vis for CMIP6 volcanic aerosol
    real(r8) :: absvisvolc(pcols)      ! AAOD vis for CMIP6 volcanic aerosol
    real(r8) :: bevisvolc(pcols,pver)  ! Extinction in vis wavelength band for CMIP6 volcanic aerosol
    real(r8) :: rhum(pcols,pver)       ! (trimmed) relative humidity for the aerosol calculations
    real(r8) :: deltah_km(pcols,pver)  ! Layer thickness, unit km
    real(r8) :: deltah, airmassl(pcols,pver), airmass(pcols) !akc6
    real(r8) :: Ca(pcols,pver), f_c(pcols,pver), f_bc(pcols,pver), f_aq(pcols,pver)
    real(r8) :: fnbc(pcols,pver), faitbc(pcols,pver), f_so4_cond(pcols,pver)
    real(r8) :: f_soa(pcols,pver),f_soana(pcols,pver)
    real(r8) :: v_soana(pcols,pver)
    real(r8) :: dCtot(pcols,pver), Ctot(pcols,pver)
    real(r8) :: Cam(pcols,pver,nbmodes), fbcm(pcols,pver,nbmodes), fcm(pcols,pver,nbmodes)
    real(r8) :: faqm(pcols,pver,nbmodes), f_condm(pcols,pver,nbmodes)
    real(r8) :: f_soam(pcols, pver,nbmodes), faqm4(pcols,pver)
    real(r8) :: focm(pcols,pver,4)
    real(r8) :: ssa(pcols,pver,0:nmodes,nbands), asym(pcols,pver,0:nmodes,nbands)
    real(r8) :: be(pcols,pver,0:nmodes,nbands), ke(pcols,pver,0:nmodes,nbands)
    real(r8) :: betotvis(pcols,pver), batotvis(pcols,pver)
    real(r8) :: ssatot(pcols,pver,nbands)     ! spectral aerosol single scattering albedo
    real(r8) :: asymtot(pcols,pver,nbands)    ! spectral aerosol asymmetry factor
    real(r8) :: betot(pcols,pver,nbands)      ! spectral aerosol extinction coefficient
    real(r8) :: batotlw(pcols,pver,nlwbands)  ! spectral aerosol absportion extinction in LW
    real(r8) :: kalw(pcols,pver,0:nmodes,nlwbands)
    real(r8) :: balw(pcols,pver,0:nmodes,nlwbands)
    real(r8) :: volc_balw(pcols,0:pver,nlwbands) ! volcanic aerosol absorption coefficient for terrestrial bands, CMIP6
    real(r8) :: rh0(pcols,pver), rhoda(pcols,pver)
    real(r8) :: ssavis(pcols,pver), asymmvis(pcols,pver), extvis(pcols,pver), dayfoc(pcols,pver)
    real(r8) :: n_aer(pcols,pver)
    real(r8) :: es(pcols,pver)      ! saturation vapor pressure
    real(r8) :: qs(pcols,pver)      ! saturation specific humidity
    real(r8) :: rht(pcols,pver)     ! relative humidity (fraction) (rh is already used in opptab)
    real(r8) :: rh_temp(pcols,pver) ! relative humidity (fraction) for input to LUT
    real(r8) :: xrh(pcols,pver)
    integer  :: irh1(pcols,pver)
    real(r8) :: xfombg(pcols,pver)
    integer  :: ifombg1(pcols,pver), ifombg2(pcols,pver)
    real(r8) :: xct(pcols,pver,nmodes)
    integer  :: ict1(pcols,pver,nmodes)
    real(r8) :: xfac(pcols,pver,nbmodes)
    integer  :: ifac1(pcols,pver,nbmodes)
    real(r8) :: xfbc(pcols,pver,nbmodes)
    integer  :: ifbc1(pcols,pver,nbmodes)
    real(r8) :: xfaq(pcols,pver,nbmodes)
    integer  :: ifaq1(pcols,pver,nbmodes)
    real(r8) :: xfbcbg(pcols,pver)
    integer  :: ifbcbg1(pcols,pver)
    real(r8) :: xfbcbgn(pcols,pver)
    integer  :: ifbcbgn1(pcols,pver)
    logical  :: lw_on   ! LW calculations are performed in interpol* if true
    !-------------------------------------------------------------------------

    ! calculate relative humidity for table lookup into rh grid
    call qsat_water(state%t(1:ncol,1:pver), state%pmid(1:ncol,1:pver), es(1:ncol,1:pver), qs(1:ncol,1:pver))

    rht(1:ncol,1:pver) = state%q(1:ncol,1:pver,1) / qs(1:ncol,1:pver)
    rh_temp(1:ncol,1:pver) = min(rht(1:ncol,1:pver),1._r8)

    do k=1,pver
       do icol=1,ncol
          ! Set upper and lower relative humidity for the aerosol calculations
          rhum(icol,k) = min(0.995_r8, max(rh_temp(icol,k), 0.01_r8))
          rhoda(icol,k) = pmid(icol,k)/(rair*t(icol,k))      ! unit kg/m^3
          if (cld(icol,k) .lt. 1.0_r8) then
             rhum(icol,k) = (rhum(icol,k) - cld(icol,k)) / (1.0_r8 - cld(icol,k))  ! clear portion
          end if
          rhum(icol,k) = min(0.995_r8, max(rhum(icol,k), 0.01_r8))
       end do
    end do

    ! Layer thickness with unit km
    do icol=1,ncol
       do k=1,pver
          deltah_km(icol,k)=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
       end do
    end do

    ! interpol-calculations only when daylight or not:
       do icol=1,ncol
          if (coszrs(icol) > 0.0_r8) then
             daylight(icol) = .true.
          else
             daylight(icol) = .false.
          endif
       end do

    ! Set SO4, BC and OC concentrations:

    ! initialize concentration fields
    do i=0,nmodes
       do k=1,pver
          do icol=1,ncol
             Nnatk(icol,k,i)  = 0.0_r8
          end do
       end do
    end do
    do k=1,pver
       do icol=1,ncol
          n_aer(icol,k)     = 0.0_r8
       end do
    end do
    kalw(:,:,:,:)=0._r8
    be(:,:,:,:)=0._r8
    ke(:,:,:,:)=0._r8
    asym(:,:,:,:)=0._r8
    ssa(:,:,:,:)=0._r8
    ! Find process tagged bulk aerosol properies (from the life cycle module):

    call calculateBulkProperties(ncol, qm1, rhoda, Nnatk, Ca, f_c, f_bc, &
         f_aq, f_so4_cond, f_soa, faitbc, fnbc, f_soana)

    ! calculating vulume fractions from mass fractions:
    do k=1,pver
       do icol=1,ncol
          v_soana(icol,k) = f_soana(icol,k)/(f_soana(icol,k) &
               +(1.0_r8-f_soana(icol,k))*rhopart(l_soa_na)/rhopart(l_so4_na))
       end do
    end do

    ! Avoid very small numbers
    do k=1,pver
       do icol=1,ncol
          Ca(icol,k)     = max(eps,Ca(icol,k))
          f_c(icol,k)    = max(eps,f_c(icol,k))
          f_bc(icol,k)   = max(eps,f_bc(icol,k))
          f_aq(icol,k)   = max(eps,f_aq(icol,k))
          fnbc(icol,k)   = max(eps,fnbc(icol,k))
          faitbc(icol,k) = max(eps,faitbc(icol,k))
       end do
    end do

    ! Calculation of the apportionment of internally mixed SO4, BC and OC
    ! mass between the various background modes.

    !==> calls modalapp to partition the mass
    call partitionMass(ncol, nnatk, Ca, f_c, f_bc, f_aq, f_so4_cond, f_soa , &
         cam, fcm, fbcm, faqm, f_condm, f_soam )

    !The following uses non-standard units, #/cm3 and ug/m3
    Nnatk(:ncol,:,:) = Nnatk(:ncol,:,:)*1.e-6_r8
    cam(:ncol,:,:)=cam(:ncol,:,:)*1.e9_r8

    ! Calculate fraction of added mass which is either SOA condensate or OC coagulate,
    ! which in AeroTab are both treated as condensate for kcomp=1-4.
    do i=1,4
       do k=1,pver
          do icol=1,ncol
             focm(icol,k,i) = fcm(icol,k,i)*(1.0_r8-fbcm(icol,k,i))
          enddo
       enddo
    enddo
    do k=1,pver
       do icol=1,ncol
          faqm4(icol,k) = faqm(icol,k,4)
       end do
    enddo

    ! find common input parameters for use in the interpolation routines
    call inputForInterpol (lchnk, ncol, rhum, xrh, irh1,   &
         f_soana, xfombg, ifombg1, faitbc, xfbcbg, ifbcbg1,  &
         fnbc, xfbcbgn, ifbcbgn1, Nnatk, Cam, xct, ict1,     &
         focm, fcm, xfac, ifac1, fbcm, xfbc, ifbc1, faqm, xfaq, ifaq1)

    ! (Wet) Optical properties for each of the aerosol modes:
    lw_on = .true.  ! No LW optics needed for RH=0 (interpol returns 0-values)

    ! BC(ax) mode (dry only):
    call interpol0 (lchnk, ncol, daylight, Nnatk, ssa, asym, be, ke, lw_on, kalw)

    mplus10=0
    ! SO4/SOA(Ait) mode:
    call interpol1 (lchnk, ncol, daylight, xrh, irh1, mplus10, &
         Nnatk, xfombg, ifombg1, xct, ict1,    &
         xfac, ifac1, ssa, asym, be, ke, lw_on, kalw)

    ! BC(Ait) and OC(Ait) modes:
    call interpol2to3 (lchnk, ncol, daylight, xrh, irh1, mplus10, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(Ait) mode:   ------ fcm invalid here (=0). Using faitbc instead
    call interpol4 (lchnk, ncol, daylight, xrh, irh1, mplus10, &
         Nnatk, xfbcbg, ifbcbg1, xct, ict1,    &
         xfac, ifac1, xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)

    ! SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
    call interpol5to10 (lchnk, ncol, daylight, xrh, irh1, &
         Nnatk, xct, ict1, xfac, ifac1, &
         xfbc, ifbc1, xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)

    ! total aerosol number concentrations
    do i=0,nmodes    ! mode 0 to 14
       do k=1,pver
          do icol=1,ncol
             n_aer(icol,k)=n_aer(icol,k)+Nnatk(icol,k,i)
          end do
       enddo
    enddo
    call outfld('N_AER   ',n_aer ,pcols,lchnk)

    mplus10=1
    ! SO4/SOA(Ait) mode:
    !does no longer exist as an externally mixed mode

    ! BC(Ait) and OC(Ait) modes:
    call interpol2to3 (lchnk, ncol, daylight, xrh, irh1, mplus10, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(n) mode:    ------ fcm not valid here (=0). Use fnbc instead
    call interpol4 (lchnk, ncol, daylight, xrh, irh1, mplus10, &
         Nnatk, xfbcbgn, ifbcbgn1, xct, ict1,  &
         xfac, ifac1, xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)

    do k=1,pver
       do icol=1,ncol
          Ctot(icol,k)=0.0_r8
       end do
    enddo

    do i=0,nmodes    ! mode 0 to 14
       do k=1,pver
          do icol=1,ncol
             dCtot(icol,k)=1.e3_r8*be(icol,k,i,4)/(ke(icol,k,i,4)+eps)
             Ctot(icol,k)=Ctot(icol,k)+dCtot(icol,k)*Nnatk(icol,k,i)
          end do
       enddo
    enddo

    ! SW Optical properties of total aerosol:
    do ib=1,nbands
       do k=1,pver
          do icol=1,ncol
             betot(icol,k,ib)=0.0_r8
             ssatot(icol,k,ib)=0.0_r8
             asymtot(icol,k,ib)=0.0_r8
          end do
       enddo
    enddo
    do ib=1,nbands
       do i=0,nmodes
          do k=1,pver
             do icol=1,ncol
                betot(icol,k,ib)=betot(icol,k,ib)+Nnatk(icol,k,i)*be(icol,k,i,ib)
                ssatot(icol,k,ib)=ssatot(icol,k,ib)+Nnatk(icol,k,i) &
                     *be(icol,k,i,ib)*ssa(icol,k,i,ib)
                asymtot(icol,k,ib)=asymtot(icol,k,ib)+Nnatk(icol,k,i) &
                     *be(icol,k,i,ib)*ssa(icol,k,i,ib)*asym(icol,k,i,ib)
             end do
          enddo
       enddo
    enddo

    ! Adding also the volcanic contribution (CMIP6), which is using a CMIP6
    ! band numbering identical to the AeroTab numbering (unlike CAM) both
    ! for SW and LW. I.e., no remapping is required here.
    ! Info from CMIP_CAM6_radiation_v3.nc
    ! wl1_sun = 0.2, 0.263158, 0.344828, 0.441501, 0.625, 0.77821, 1.24224,
    ! 1.2987, 1.62602, 1.94175, 2.15054, 2.5, 3.07692, 3.84615 ;
    ! wl2_sun = 0.263158, 0.344828, 0.441501, 0.625, 0.77821, 1.24224, 1.2987,
    ! 1.62602, 1.94175, 2.15054, 2.5, 3.07692, 3.84615, 12.1951 ;
    ! wl1_earth = 3.07692, 3.84615, 4.20168, 4.44444, 4.80769, 5.55556, 6.75676,
    ! 7.19424, 8.47458, 9.25926, 10.2041, 12.1951, 14.2857, 15.873, 20, 28.5714 ;
    ! wl2_earth = 3.84615, 4.20168, 4.44444, 4.80769, 5.55556, 6.75676, 7.19424,
    ! 8.47458, 9.25926, 10.2041, 12.1951, 14.2857, 15.873, 20, 28.5714, 1000 ;
    do ib=1,nbands
       betot(1:ncol,1:pver,ib) = betot(1:ncol,1:pver,ib) &
            + volc_ext_sun(1:ncol,1:pver,ib)
       ssatot(1:ncol,1:pver,ib) = ssatot(1:ncol,1:pver,ib) &
            + volc_ext_sun(1:ncol,1:pver,ib)*volc_omega_sun(1:ncol,1:pver,ib)
       asymtot(1:ncol,1:pver,ib) = asymtot(1:ncol,1:pver,ib) &
            + volc_ext_sun(1:ncol,1:pver,ib)*volc_omega_sun(1:ncol,1:pver,ib) &
            *volc_g_sun(1:ncol,1:pver,ib)
    enddo
    bevisvolc(1:ncol,1:pver) = volc_ext_sun(1:ncol,1:pver,4)

    ! and then calculate the total bulk optical parameters
    do ib=1,nbands
       do k=1,pver
          do icol=1,ncol
             ssatot(icol,k,ib)=ssatot(icol,k,ib)/(betot(icol,k,ib)+eps)
             asymtot(icol,k,ib)=asymtot(icol,k,ib) &
                  /(betot(icol,k,ib)*ssatot(icol,k,ib)+eps)
          end do
       enddo
    enddo

    !------------------------------------------------------------------------------------------------
    ! Replace CAM5 standard aerosol optics with CAM5-Oslo optics (except top layer: no aerosol)
    ! Remapping from AeroTab to CAM5 SW bands, see p. 167 in the CAM5.0 description:
    ! CAM5 bands         AeroTab bands
    ! 14 3.846 12.195        14
    ! 1 3.077 3.846          13
    ! 2 2.500 3.077          12
    ! 3 2.150 2.500          11
    ! 4 1.942 2.150          10
    ! 5 1.626 1.942           9
    ! 6 1.299 1.626           8
    ! 7 1.242 1.299           7
    ! 8 0.778 1.242           6
    ! 9 0.625 0.778           5
    ! 10 0.442 0.625          4
    ! 11 0.345 0.442          3
    ! 12 0.263 0.345          2
    ! 13 0.200 0.263          1

    do i=1,ncol  ! zero aerosol in the top layer
       do ib=1,14 ! 1-nbands
          per_tau(i,0,ib)= 0._r8
          per_tau_w(i,0,ib)= 0.999_r8
          per_tau_w_g(i,0,ib)= 0.5_r8
          per_tau_w_f(i,0,ib)= 0.25_r8
       end do
       do ib=1,14  ! initialize also for the other layers
          do k=1,pver
             per_tau(i,k,ib)= 0._r8
             per_tau_w(i,k,ib)= 0.999_r8
             per_tau_w_g(i,k,ib)= 0.5_r8
             per_tau_w_f(i,k,ib)= 0.25_r8
          end do
       end do
    end do
    ! Remapping of SW wavelength bands from AeroTab to CAM5
    do i=1,ncol
       do ib=1,13
          do k=1,pver
             per_tau(i,k,ib)=deltah_km(i,k)*betot(i,k,14-ib)
             per_tau_w(i,k,ib)=per_tau(i,k,ib)*max(min(ssatot(i,k,14-ib),0.999999_r8),1.e-6_r8)
             per_tau_w_g(i,k,ib)=per_tau_w(i,k,ib)*asymtot(i,k,14-ib)
             per_tau_w_f(i,k,ib)=per_tau_w_g(i,k,ib)*asymtot(i,k,14-ib)
          end do
       end do
       ib=14
       do k=1,pver
          per_tau(i,k,ib)=deltah_km(i,k)*betot(i,k,ib)
          per_tau_w(i,k,ib)=per_tau(i,k,ib)*max(min(ssatot(i,k,ib),0.999999_r8),1.e-6_r8)
          per_tau_w_g(i,k,ib)=per_tau_w(i,k,ib)*asymtot(i,k,ib)
          per_tau_w_f(i,k,ib)=per_tau_w_g(i,k,ib)*asymtot(i,k,ib)
       end do
    end do  ! ncol
    !------------------------------------------------------------------------------------------------

    ! LW Optical properties of total aerosol:
    do ib=1,nlwbands
       do k=1,pver
          do icol=1,ncol
             batotlw(icol,k,ib)=0.0_r8
          end do
       enddo
    enddo
    do ib=1,nlwbands
       do i=0,nmodes
          do k=1,pver
             do icol=1,ncol
                balw(icol,k,i,ib)=kalw(icol,k,i,ib)*(be(icol,k,i,4)/(ke(icol,k,i,4)+eps))
                batotlw(icol,k,ib)=batotlw(icol,k,ib)+Nnatk(icol,k,i)*balw(icol,k,i,ib)
             end do
          enddo
       enddo
    enddo

    ! Adding also the volcanic contribution (CMIP6), which is also using
    ! AeroTab band numbering, so that a remapping is required here
    do ib=1,nlwbands
       volc_balw(1:ncol,1:pver,ib) = volc_ext_earth(:ncol,1:pver,ib) &
            *(1.0_r8-volc_omega_earth(:ncol,1:pver,ib))
       batotlw(1:ncol,1:pver,ib)=batotlw(1:ncol,1:pver,ib)+volc_balw(1:ncol,1:pver,ib)
    enddo

    ! Remapping of LW wavelength bands from AeroTab to CAM5
    do ib=1,nlwbands
       do i=1,ncol
          do k=1,pver
             per_lw_abs(i,k,ib)=deltah_km(i,k)*batotlw(i,k,17-ib)
             ! if(ib.eq.1.and.k.eq.pver.and.i.eq.1) then
             ! write(*,*) 'per_lw_abs =', per_lw_abs(i,k,ib)
             ! endif
          end do
       end do
    end do

    ! APPROXIMATE aerosol extinction and absorption at 550nm (0.442-0.625 um)
    ! (in the visible wavelength band)
    do k=1,pver
       do icol=1,ncol
          betotvis(icol,k)=betot(icol,k,4)
          batotvis(icol,k)=betotvis(icol,k)*(1.0-ssatot(icol,k,4))
       end do
    enddo

    do k=1,pver
       do icol=1,ncol
          ssavis(icol,k) = 0.0_r8
          asymmvis(icol,k) = 0.0_r8
          extvis(icol,k) = 0.0_r8
          dayfoc(icol,k) = 0.0_r8
       enddo
    end do

    do k=1,pver
       do icol=1,ncol
          ! dayfoc < 1 when looping only over gridcells with daylight
          if(daylight(icol)) then
             dayfoc(icol,k) = 1.0_r8
             ! with the new bands in CAM5, band 4 is now at ca 0.5 um (0.442-0.625)
             ssavis(icol,k) = ssatot(icol,k,4)
             asymmvis(icol,k) = asymtot(icol,k,4)
             extvis(icol,k) = betot(icol,k,4)
          endif
       enddo
    end do

    ! optical parameters in visible light (0.442-0.625um)
    call outfld('SSAVIS  ',ssavis,pcols,lchnk)
    call outfld('ASYMMVIS',asymmvis,pcols,lchnk)
    call outfld('EXTVIS  ',extvis,pcols,lchnk)
    call outfld('DAYFOC  ',dayfoc,pcols,lchnk)

    ! Initialize fields
    do icol=1,ncol
       aodvis(icol)=0.0_r8
       absvis(icol)=0.0_r8
       aodvisvolc(icol)=0.0_r8
       absvisvolc(icol)=0.0_r8
       airmass(icol)=0.0_r8  !akc6
    enddo

    do icol=1,ncol
       if(daylight(icol)) then
          do k=1,pver
             ! Layer thickness, unit km, and layer airmass, unit kg/m2
             deltah=deltah_km(icol,k)
             airmassl(icol,k)=1.e3_r8*deltah*rhoda(icol,k)
             airmass(icol)=airmass(icol)+airmassl(icol,k)  !akc6

             ! Optical depths at ca. 550 nm (0.442-0.625um) all aerosols
             aodvis(icol)=aodvis(icol)+betotvis(icol,k)*deltah
             absvis(icol)=absvis(icol)+batotvis(icol,k)*deltah

             ! Optical depths at ca. 550 nm (0.442-0.625um) CMIP6 volcanic aerosol
             aodvisvolc(icol)=aodvisvolc(icol)+volc_ext_sun(icol,k,4)*deltah
             absvisvolc(icol)=absvisvolc(icol)+volc_ext_sun(icol,k,4) &
                  *(1.0_r8-volc_omega_sun(icol,k,4))*deltah

          end do  ! k
       endif   ! daylight
    end do   ! icol

    ! Extinction and absorption for 0.55 um for the total aerosol, and AODs
    call outfld('AOD_VIS ',aodvis ,pcols,lchnk)
    call outfld('ABSVIS  ',absvis ,pcols,lchnk)
    call outfld('AODVVOLC',aodvisvolc ,pcols,lchnk)
    call outfld('ABSVVOLC',absvisvolc ,pcols,lchnk)
    call outfld('BVISVOLC',bevisvolc ,pcols,lchnk)

  end subroutine oslo_aero_optical_params_calc

end module oslo_aero_optical_params
