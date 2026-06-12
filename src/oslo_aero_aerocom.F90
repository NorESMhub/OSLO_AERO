module oslo_aero_aerocom

  use shr_kind_mod,             only: r8 => shr_kind_r8
  use ppgrid,                   only: pcols, pver, pverp
  use cam_history,              only: outfld
  !
  use oslo_aero_share,          only: cate, cat, fac, faq, fbc, rh, fombg, fbcbg, rh, xrhrf, irhrf1, eps
  use oslo_aero_sw_tables,      only: interpol0, interpol1, interpol2to3, interpol4, interpol5to10
  use oslo_aero_aerodry_tables, only: intdrypar0, intdrypar1, intdrypar2to3, intdrypar4, intdrypar5to10
  use oslo_aero_aerocom_tables, only: intaeropt0, intaeropt1, intaeropt2to3, intaeropt4, intaeropt5to10
  use oslo_aero_linear_interp , only: lininterpol3dim, lininterpol4dim, lininterpol5dim
  use oslo_aero_share,          only: rhopart, l_bc_ni, l_om_ni
  use oslo_aero_share,          only: nmodes, nbmodes, nbands, nlwbands, nbmp1
  use oslo_aero_control,        only: rh_fine_aer_scale_fact_optics

  public  :: aerocom1
  public  :: aerocom2

  private :: opticsAtConstRh

  ! Used by radiation.F90
  real(r8), public, protected :: dod440(pcols)
  real(r8), public, protected :: dod550(pcols)
  real(r8), public, protected :: dod870(pcols)
  real(r8), public, protected :: abs550(pcols)
  real(r8), public, protected :: abs550alt(pcols)

!===============================================================================
contains
!===============================================================================

  subroutine aerocom1(lchnk, ncol, Cam, Nnatk, deltah_km, &
       xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, xfbcbg, ifbcbg1, xfbcbgn, ifbcbgn1, &
       xfombg, ifombg1, Ctotdry)

    ! Called by oslo_aero_optical_params

    ! Arguments
    integer , intent(in)  :: lchnk                        ! chunk identifier
    integer , intent(in)  :: ncol                         ! number of atmospheric columns
    real(r8), intent(in)  :: Cam(pcols,pver,nbmodes)
    real(r8), intent(in)  :: Nnatk(pcols,pver,0:nmodes)   ! modal aerosol number concentration
    real(r8), intent(in)  :: deltah_km(pcols,pver)        ! layer thickness, unit km
    real(r8), intent(in)  :: xct(pcols,pver,nmodes)       ! modal internally mixed SO4+BC+OC conc.
    integer , intent(in)  :: ict1(pcols,pver,nmodes)
    real(r8), intent(in)  :: xfac(pcols,pver,nbmodes)
    integer , intent(in)  :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in)  :: xfbc(pcols,pver,nbmodes)
    integer , intent(in)  :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in)  :: xfaq(pcols,pver,nbmodes)
    integer , intent(in)  :: ifaq1(pcols,pver,nbmodes)
    real(r8), intent(in)  :: xfbcbg(pcols,pver)
    integer , intent(in)  :: ifbcbg1(pcols,pver)
    real(r8), intent(in)  :: xfbcbgn(pcols,pver)
    integer , intent(in)  :: ifbcbgn1(pcols,pver)
    real(r8), intent(in)  :: xfombg(pcols,pver)
    integer , intent(in)  :: ifombg1(pcols,pver)
    real(r8), intent(out) :: Ctotdry(pcols,pver)

    ! Local variables
    integer  :: ilev, icol, imode, kcomp
    real(r8) :: asydry_aer(pcols,pver)     ! dry asymtot in the visible band
    real(r8) :: asymtot(pcols,pver,nbands) ! spectral aerosol asymmetry factor
    real(r8) :: ssatot(pcols,pver,nbands)  ! spectral aerosol single scattering albedo
    real(r8) :: betot(pcols,pver,nbands)   ! spectral aerosol extinction coefficient
    real(r8) :: be(pcols,pver,0:nmodes,nbands)
    real(r8) :: ke(pcols,pver,0:nmodes,nbands)
    real(r8) :: ssa(pcols,pver,0:nmodes,nbands)
    real(r8) :: asym(pcols,pver,0:nmodes,nbands)
    real(r8) :: kalw(pcols,pver,0:nmodes,nlwbands)
    real(r8) :: xrhnull(pcols,pver)
    integer  :: irh1null(pcols,pver)
    real(r8) :: dCtot(pcols,pver)
    real(r8) :: Ctot(pcols,pver)
    real(r8) :: Camrel(pcols,pver,nbmodes)
    real(r8) :: Camtot(pcols,nbmodes)
    real(r8) :: cxsmtot(pcols,nbmodes)
    real(r8) :: cxsmrel(pcols,nbmodes)
    real(r8) :: cxs(pcols,pver)
    real(r8) :: cxstot(pcols,pver)
    real(r8) :: akcxs(pcols)
    logical  :: daylight(pcols) ! SW calculations also at (polar) night in interpol* if daylight=.true.
    logical  :: lw_on           ! LW calculations are performed in interpol* if true
    real(r8) :: xctrel,camdiff,cxsm
    character(len=10) :: modeString
    character(len=20) :: varname
    !-------------------------------------------------------------------------

    ! Initialize overshooting mass summed over all modes
    do ilev=1,pver
       do icol=1,ncol
          cxstot(icol,ilev) = 0.0_r8
       enddo
    enddo

    ! Initializing total and relative exessive (overshooting
    ! w.r.t. look-up table maxima) added mass column:
    do imode=1,nbmodes
       do icol=1,ncol
          Camtot(icol,imode)=0.0_r8
          cxsmtot(icol,imode)=0.0_r8
          cxsmrel(icol,imode)=0.0_r8
       enddo
    enddo

    ! Calculating added internally mixed mass onto each mode 1-10, relative to
    ! maximum mass which can be added w.r.t. the look-up tables (for level ilev),
    ! as well as the relative exessive added mass column:
    do imode=1,4
       do ilev=1,pver
          do icol=1,ncol
             Camrel(icol,ilev,imode) = (Cam(icol,ilev,imode)/(Nnatk(icol,ilev,imode)+eps))/cate(imode,16)
             xctrel=min(max(Camrel(icol,ilev,imode),cate(imode,1)/cate(imode,16)),1.0_r8)
             camdiff=Cam(icol,ilev,imode)-xctrel*cate(imode,16)*(Nnatk(icol,ilev,imode)+eps)
             cxsm=max(0.0_r8,camdiff)
             cxsmtot(icol,imode)=cxsmtot(icol,imode)+cxsm*deltah_km(icol,ilev)
             Camtot(icol,imode)=Camtot(icol,imode)+Cam(icol,ilev,imode)*deltah_km(icol,ilev)
             camdiff=Cam(icol,ilev,imode)-xct(icol,ilev,imode)*(Nnatk(icol,ilev,imode)+eps)
             cxs(icol,ilev)=max(0.0_r8,camdiff)
             cxstot(icol,ilev)= cxstot(icol,ilev)+cxs(icol,ilev)
          enddo
       enddo
    enddo
    do imode=5,nbmodes
       do ilev=1,pver
          do icol=1,ncol
             Camrel(icol,ilev,imode) = (Cam(icol,ilev,imode)/(Nnatk(icol,ilev,imode)+eps))/cat(imode,6)
             xctrel=min(max(Camrel(icol,ilev,imode),cat(imode,1)/cat(imode,6)),1.0_r8)
             camdiff=Cam(icol,ilev,imode)-xctrel*cat(imode,6)*(Nnatk(icol,ilev,imode)+eps)
             cxsm=max(0.0_r8,camdiff)
             cxsmtot(icol,imode)=cxsmtot(icol,imode)+cxsm*deltah_km(icol,ilev)
             Camtot(icol,imode)=Camtot(icol,imode)+Cam(icol,ilev,imode)*deltah_km(icol,ilev)
             camdiff=Cam(icol,ilev,imode)-xct(icol,ilev,imode)*(Nnatk(icol,ilev,imode)+eps)
             cxs(icol,ilev)=max(0.0_r8,camdiff)
             cxstot(icol,ilev)= cxstot(icol,ilev) + cxs(icol,ilev)
          enddo
       enddo
    enddo

    ! Total overshooting mass summed over all modes and all levels
    do icol=1,ncol
       akcxs(icol) = 0.0_r8
       do ilev=1,pver
          akcxs(icol) = akcxs(icol) + cxstot(icol,ilev)*deltah_km(icol,ilev)
       enddo
    enddo
    call outfld('AKCXS   ',akcxs ,pcols,lchnk)

    do imode=1,nbmodes
       do icol=1,ncol
          cxsmrel(icol,imode)=cxsmtot(icol,imode)/(Camtot(icol,imode)+eps)
       enddo
    enddo

    do imode=1,nbmodes
       modeString="  "
       write(modeString,"(I2)"),imode
       if(imode<10) modeString="0"//adjustl(modeString)
       varName = "Camrel"//trim(modeString)
       if(imode.ne.3) call outfld(varName,Camrel(:,:,imode),pcols,lchnk)
    enddo

    do imode=1,nbmodes
       modeString="  "
       write(modeString,"(I2)"),imode
       if(imode<10) modeString="0"//adjustl(modeString)
       varName = "Cxsrel"//trim(modeString)
       if(imode.ne.3) call outfld(varName,cxsmrel(:,imode),pcols,lchnk)
    enddo

    ! Find dry aerosol asymmetry factor and mass for subsequent
    ! calculation of condensed water mass below...
    do ilev=1,pver
       do icol=1,ncol
          asydry_aer(icol,ilev)=0.0_r8
       end do
    enddo

    ! Note: using xrhnull etc as proxy for constant RH input values
    xrhnull(:,:)  = 0.0_r8
    irh1null(:,:) = 1

    ! For aerocom daylight is always .true.
    daylight(:) = .true.

    ! No LW optics needed for RH=0 (interpol returns 0-values)
    lw_on = .false.

    ! BC(ax) mode (dry only):
    call interpol0 (ncol, daylight, Nnatk, ssa, asym, be, ke, lw_on, kalw)

    ! SO4/SOA(Ait) mode:
    call interpol1 (ncol, daylight, xrhnull, irh1null, 1, &
         Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC(Ait) and OC(Ait) modes:
    call interpol2to3 (ncol, daylight, xrhnull, irh1null, 2, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(Ait) mode: fcm not valid here (=0).
    call interpol4 (ncol, daylight, xrhnull, irh1null, 4, &
         Nnatk, xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, &
         xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)

    ! SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
    do kcomp=5,10
       call interpol5to10 (ncol, daylight, xrhnull, irh1null, kcomp, &
            Nnatk, xct, ict1, xfac, ifac1, &
            xfbc, ifbc1, xfaq, ifaq1, &
            ssa, asym, be, ke, lw_on, kalw)
    end do

    ! BC(Ait) and OC(Ait) nucleation modes:
    call interpol2to3 (ncol, daylight, xrhnull, irh1null, 12, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(n) mode:
    call interpol4 (ncol, daylight, xrhnull, irh1null, 14, &
         Nnatk, xfbcbgn, ifbcbgn1, xct, ict1, &
         xfac, ifac1, xfaq, ifaq1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! Compute Ctotdry
    Ctotdry(:,:) = 0.0_r8
    do imode=0,nmodes    ! mode 0 to 14
       do ilev=1,pver
          do icol=1,ncol
             dCtot(icol,ilev)=1.e3_r8*be(icol,ilev,imode,4)/(ke(icol,ilev,imode,4)+eps)
             Ctotdry(icol,ilev) = Ctotdry(icol,ilev)+dCtot(icol,ilev)*Nnatk(icol,ilev,imode)
          end do
       enddo
    enddo

    ! AeroCom Phase III: adding asymmetry factor for dry aerosol, wavelength band 4 only
    ! (and with no CMIP6 volcnic contribution)
    ib=4
    do ilev=1,pver
       do icol=1,ncol
          betot(icol,ilev,ib)=0.0_r8
          ssatot(icol,ilev,ib)=0.0_r8
          asymtot(icol,ilev,ib)=0.0_r8
       end do
    enddo
    do imode=0,nmodes
       do ilev=1,pver
          do icol=1,ncol
             betot(icol,ilev,ib)  =betot(icol,ilev,ib)  &
                  +Nnatk(icol,ilev,imode)*be(icol,ilev,imode,ib)
             ssatot(icol,ilev,ib) =ssatot(icol,ilev,ib) &
                  +Nnatk(icol,ilev,imode)*be(icol,ilev,imode,ib)*ssa(icol,ilev,imode,ib)
             asymtot(icol,ilev,ib)=asymtot(icol,ilev,ib)&
                  +Nnatk(icol,ilev,imode)*be(icol,ilev,imode,ib)*ssa(icol,ilev,imode,ib)*asym(icol,ilev,imode,ib)
          end do
       enddo
    enddo

    do ilev=1,pver
       do icol=1,ncol
          ssatot(icol,ilev,ib) = ssatot(icol,ilev,ib) /(betot(icol,ilev,ib)+eps)
          asymtot(icol,ilev,ib)= asymtot(icol,ilev,ib)/(betot(icol,ilev,ib)*ssatot(icol,ilev,ib)+eps)
          asydry_aer(icol,ilev)= asymtot(icol,ilev,ib)
       end do
    enddo
    call outfld('ASYMMDRY',asydry_aer,pcols,lchnk)

  end subroutine aerocom1

  !===============================================================================

  subroutine aerocom2(lchnk, ncol, Nnatk, pint, deltah_km, faitbc, f_soana, fnbc, rhoda, v_soana, &
       xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, xfbcbg, ifbcbg1, xfbcbgn, ifbcbgn1, &
       xfombg, ifombg1, xrh, irh1)

    ! Arguments
    integer , intent(in) :: lchnk                        ! chunk identifier
    integer , intent(in) :: ncol                         ! number of atmospheric columns
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes)   ! aerosol mode
    real(r8), intent(in) :: pint(pcols,pverp)            ! Model interface pressures (10*Pa)
    real(r8), intent(in) :: deltah_km(pcols,pver)        ! Layer thickness, unit km
    real(r8), intent(in) :: faitbc(pcols,pver)
    real(r8), intent(in) :: f_soana(pcols,pver)
    real(r8), intent(in) :: fnbc(pcols,pver)
    real(r8), intent(in) :: rhoda(pcols,pver)
    real(r8), intent(in) :: v_soana(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)       ! modal internally mixed SO4+BC+OC conc.
    integer , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)
    integer , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)
    integer , intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)
    integer , intent(in) :: ifaq1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbcbg(pcols,pver)
    integer , intent(in) :: ifbcbg1(pcols,pver)
    real(r8), intent(in) :: xfbcbgn(pcols,pver)
    integer , intent(in) :: ifbcbgn1(pcols,pver)
    real(r8), intent(in) :: xfombg(pcols,pver)
    integer , intent(in) :: ifombg1(pcols,pver)
    real(r8), intent(in) :: xrh(pcols,pver)
    integer , intent(in) :: irh1(pcols,pver)

    ! local variables
    integer  :: ilev, icol, imode
    real(r8) :: dload(pcols,0:nmodes)
    real(r8) :: vaercols(pcols)
    real(r8) :: vaercoll(pcols)
    real(r8) :: aaercols(pcols)
    real(r8) :: aaercoll(pcols)
    !
    real(r8) :: bext440n(pcols,pver,0:nbmodes), babs440n(pcols,pver,0:nbmodes)
    real(r8) :: bext500n(pcols,pver,0:nbmodes), babs500n(pcols,pver,0:nbmodes)
    real(r8) :: bext550n(pcols,pver,0:nbmodes), babs550n(pcols,pver,0:nbmodes)
    real(r8) :: bext670n(pcols,pver,0:nbmodes), babs670n(pcols,pver,0:nbmodes)
    real(r8) :: bext870n(pcols,pver,0:nbmodes), babs870n(pcols,pver,0:nbmodes)
    !
    real(r8) :: bebg440n(pcols,pver,0:nbmodes), babg440n(pcols,pver,0:nbmodes)
    real(r8) :: bebg500n(pcols,pver,0:nbmodes), babg500n(pcols,pver,0:nbmodes)
    real(r8) :: bebg550n(pcols,pver,0:nbmodes), babg550n(pcols,pver,0:nbmodes)
    real(r8) :: bebg670n(pcols,pver,0:nbmodes), babg670n(pcols,pver,0:nbmodes)
    real(r8) :: bebg870n(pcols,pver,0:nbmodes), babg870n(pcols,pver,0:nbmodes)
    !
    real(r8) :: bebc440n(pcols,pver,0:nbmodes), babc440n(pcols,pver,0:nbmodes)
    real(r8) :: bebc500n(pcols,pver,0:nbmodes), babc500n(pcols,pver,0:nbmodes)
    real(r8) :: bebc550n(pcols,pver,0:nbmodes), babc550n(pcols,pver,0:nbmodes)
    real(r8) :: bebc670n(pcols,pver,0:nbmodes), babc670n(pcols,pver,0:nbmodes)
    real(r8) :: bebc870n(pcols,pver,0:nbmodes), babc870n(pcols,pver,0:nbmodes)
    !
    real(r8) :: beoc440n(pcols,pver,0:nbmodes), baoc440n(pcols,pver,0:nbmodes)
    real(r8) :: beoc500n(pcols,pver,0:nbmodes), baoc500n(pcols,pver,0:nbmodes)
    real(r8) :: beoc550n(pcols,pver,0:nbmodes), baoc550n(pcols,pver,0:nbmodes)
    real(r8) :: beoc670n(pcols,pver,0:nbmodes), baoc670n(pcols,pver,0:nbmodes)
    real(r8) :: beoc870n(pcols,pver,0:nbmodes), baoc870n(pcols,pver,0:nbmodes)
    !
    real(r8) :: besu440n(pcols,pver,0:nbmodes), basu440n(pcols,pver,0:nbmodes)
    real(r8) :: besu500n(pcols,pver,0:nbmodes), basu500n(pcols,pver,0:nbmodes)
    real(r8) :: besu550n(pcols,pver,0:nbmodes), basu550n(pcols,pver,0:nbmodes)
    real(r8) :: besu670n(pcols,pver,0:nbmodes), basu670n(pcols,pver,0:nbmodes)
    real(r8) :: besu870n(pcols,pver,0:nbmodes), basu870n(pcols,pver,0:nbmodes)
    !
    real(r8) :: bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes)
    real(r8) :: bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes)
    real(r8) :: bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes)
    real(r8) :: bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes)
    real(r8) :: bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes)
    !
    real(r8) :: bebg440(pcols,pver,0:nbmodes), babg440(pcols,pver,0:nbmodes)
    real(r8) :: bebg500(pcols,pver,0:nbmodes), babg500(pcols,pver,0:nbmodes)
    real(r8) :: bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes)
    real(r8) :: bebg670(pcols,pver,0:nbmodes), babg670(pcols,pver,0:nbmodes)
    real(r8) :: bebg870(pcols,pver,0:nbmodes), babg870(pcols,pver,0:nbmodes)
    !
    real(r8) :: bebc440(pcols,pver,0:nbmodes), babc440(pcols,pver,0:nbmodes)
    real(r8) :: bebc500(pcols,pver,0:nbmodes), babc500(pcols,pver,0:nbmodes)
    real(r8) :: bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes)
    real(r8) :: bebc670(pcols,pver,0:nbmodes), babc670(pcols,pver,0:nbmodes)
    real(r8) :: bebc870(pcols,pver,0:nbmodes), babc870(pcols,pver,0:nbmodes)
    !
    real(r8) :: beoc440(pcols,pver,0:nbmodes), baoc440(pcols,pver,0:nbmodes)
    real(r8) :: beoc500(pcols,pver,0:nbmodes), baoc500(pcols,pver,0:nbmodes)
    real(r8) :: beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes)
    real(r8) :: beoc670(pcols,pver,0:nbmodes), baoc670(pcols,pver,0:nbmodes)
    real(r8) :: beoc870(pcols,pver,0:nbmodes), baoc870(pcols,pver,0:nbmodes)
    !
    real(r8) :: besu440(pcols,pver,0:nbmodes), basu440(pcols,pver,0:nbmodes)
    real(r8) :: besu500(pcols,pver,0:nbmodes), basu500(pcols,pver,0:nbmodes)
    real(r8) :: besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes)
    real(r8) :: besu670(pcols,pver,0:nbmodes), basu670(pcols,pver,0:nbmodes)
    real(r8) :: besu870(pcols,pver,0:nbmodes), basu870(pcols,pver,0:nbmodes)
    !
    real(r8) :: be(pcols,pver,0:nmodes,nbands)
    real(r8) :: ke(pcols,pver,0:nmodes,nbands)
    real(r8) :: dod550dry(pcols), abs550dry(pcols)
    !
    real(r8) :: dload3d(pcols,pver,0:nmodes)
    real(r8) :: dload_mi(pcols), dload_ss(pcols)
    real(r8) :: dload_s4(pcols), dload_oc(pcols), dload_bc(pcols)
    real(r8) :: dload_s4_a(pcols), dload_s4_1(pcols), dload_s4_5(pcols)
    real(r8) :: dload_bc_0(pcols), dload_bc_ac(pcols), dload_oc_ac(pcols)
    real(r8) :: dload_bc_2(pcols), dload_bc_4(pcols), dload_bc_12(pcols), dload_bc_14(pcols)
    real(r8) :: dload_oc_4(pcols), dload_oc_14(pcols)
    !
    real(r8) :: cmin(pcols,pver), cseas(pcols,pver)
    !
    real(r8) :: nnat_1(pcols,pver), nnat_2(pcols,pver), nnat_3(pcols,pver)
    real(r8) :: nnat_4(pcols,pver), nnat_5(pcols,pver), nnat_6(pcols,pver)
    real(r8) :: nnat_7(pcols,pver), nnat_8(pcols,pver), nnat_9(pcols,pver)
    real(r8) :: nnat_10(pcols,pver), nnat_12(pcols,pver)
    real(r8) :: nnat_14(pcols,pver), nnat_0(pcols,pver)
    !
    real(r8) :: ck(pcols,pver,0:nmodes), cknorm(pcols,pver,0:nmodes)
    real(r8) :: cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)
    !
    real(r8) :: aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes)
    real(r8) :: vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes)
    real(r8) :: aaeros(pcols,pver,0:nbmodes), aaerol(pcols,pver,0:nbmodes)
    real(r8) :: vaeros(pcols,pver,0:nbmodes), vaerol(pcols,pver,0:nbmodes)
    !
    real(r8) :: cintbg(pcols,pver,0:nbmodes)
    real(r8) :: cintbg05(pcols,pver,0:nbmodes), cintbg125(pcols,pver,0:nbmodes)
    real(r8) :: cintbc(pcols,pver,0:nbmodes)
    real(r8) :: cintbc05(pcols,pver,0:nbmodes), cintbc125(pcols,pver,0:nbmodes)
    real(r8) :: cintoc(pcols,pver,0:nbmodes)
    real(r8) :: cintoc05(pcols,pver,0:nbmodes), cintoc125(pcols,pver,0:nbmodes)
    real(r8) :: cintsc(pcols,pver,0:nbmodes)
    real(r8) :: cintsc05(pcols,pver,0:nbmodes), cintsc125(pcols,pver,0:nbmodes)
    real(r8) :: cintsa(pcols,pver,0:nbmodes)
    real(r8) :: cintsa05(pcols,pver,0:nbmodes), cintsa125(pcols,pver,0:nbmodes)
    !
    real(r8) :: c_mi(pcols,pver), c_mi05(pcols,pver), c_mi125(pcols,pver)
    real(r8) :: c_ss(pcols,pver), c_ss05(pcols,pver), c_ss125(pcols,pver)
    real(r8) :: c_bc(pcols,pver), c_bc05(pcols,pver), c_bc125(pcols,pver)
    real(r8) :: c_oc(pcols,pver), c_oc05(pcols,pver), c_oc125(pcols,pver)
    real(r8) :: c_sa(pcols,pver), c_sa05(pcols,pver), c_sa125(pcols,pver)
    real(r8) :: c_sc(pcols,pver), c_sc05(pcols,pver), c_sc125(pcols,pver)
    real(r8) :: c_s4(pcols,pver), c_s405(pcols,pver), c_s4125(pcols,pver)
    real(r8) :: c_s4_a(pcols,pver), c_s4_1(pcols,pver), c_s4_5(pcols,pver)
    real(r8) :: c_bc_0(pcols,pver), c_bc_ac(pcols,pver), c_oc_ac(pcols,pver)
    real(r8) :: c_bc_2(pcols,pver), c_bc_4(pcols,pver), c_bc_12(pcols,pver), c_bc_14(pcols,pver)
    real(r8) :: c_oc_4(pcols,pver), c_oc_14(pcols,pver)
    real(r8) :: c_tots(pcols), c_tot125s(pcols), c_pm25s(pcols) ! = PM all sizes, PM>2.5um and PM<2.5um (PM2.5)
    real(r8) :: c_tot(pcols,pver), c_tot125(pcols,pver), c_pm25(pcols,pver)
    real(r8) :: c_tot05(pcols,pver), c_pm1(pcols,pver)
    !
    real(r8) :: mmr_pm25(pcols,pver), mmr_pm1(pcols,pver)
    real(r8) :: aaeros_tot(pcols,pver), aaerol_tot(pcols,pver), vaeros_tot(pcols,pver)
    real(r8) :: vaerol_tot(pcols,pver)
    real(r8) :: derlt05(pcols), dergt05(pcols), der(pcols)
    real(r8) :: erlt053d(pcols,pver), ergt053d(pcols,pver), er3d(pcols,pver)
    real(r8) :: bebglt1(pcols,pver,0:nbmodes), bebggt1(pcols,pver,0:nbmodes)
    real(r8) :: bebclt1(pcols,pver,0:nbmodes), bebcgt1(pcols,pver,0:nbmodes)
    real(r8) :: beoclt1(pcols,pver,0:nbmodes), beocgt1(pcols,pver,0:nbmodes)
    real(r8) :: bes4lt1(pcols,pver,0:nbmodes), bes4gt1(pcols,pver,0:nbmodes)
    real(r8) :: backsc550(pcols,pver,0:nbmodes), backsc550x(pcols,pver,nbmp1:nmodes)
    real(r8) :: backsc550tot(pcols,pver), ec550_aer(pcols,pver), abs550_aer(pcols,pver)
    real(r8) :: bs550_aer(pcols,pver)
    real(r8) :: ec550_so4(pcols,pver),ec550_bc(pcols,pver), ec550_pom(pcols,pver)
    real(r8) :: ec550_ss(pcols,pver), ec550_du(pcols,pver)
    real(r8) :: bebglt1n(pcols,pver,0:nbmodes), bebggt1n(pcols,pver,0:nbmodes)
    real(r8) :: bebclt1n(pcols,pver,0:nbmodes), bebcgt1n(pcols,pver,0:nbmodes)
    real(r8) :: beoclt1n(pcols,pver,0:nbmodes), beocgt1n(pcols,pver,0:nbmodes)
    real(r8) :: bes4lt1n(pcols,pver,0:nbmodes), bes4gt1n(pcols,pver,0:nbmodes)
    real(r8) :: backsc550n(pcols,pver,0:nbmodes)
    !
    real(r8) :: bext440tot(pcols,pver), babs440tot(pcols,pver)
    real(r8) :: bext500tot(pcols,pver), babs500tot(pcols,pver)
    real(r8) :: bext550tot(pcols,pver), babs550tot(pcols,pver)
    real(r8) :: bext670tot(pcols,pver), babs670tot(pcols,pver)
    real(r8) :: bext870tot(pcols,pver), babs870tot(pcols,pver)
    !
    real(r8) :: bebg440tot(pcols,pver), babg440tot(pcols,pver)
    real(r8) :: bebg500tot(pcols,pver), babg500tot(pcols,pver)
    real(r8) :: bebg550tot(pcols,pver), babg550tot(pcols,pver)
    real(r8) :: bebg670tot(pcols,pver), babg670tot(pcols,pver)
    real(r8) :: bebg870tot(pcols,pver), babg870tot(pcols,pver)
    !
    real(r8) :: bebc440tot(pcols,pver), babc440tot(pcols,pver)
    real(r8) :: bebc500tot(pcols,pver), babc500tot(pcols,pver)
    real(r8) :: bebc550tot(pcols,pver), babc550tot(pcols,pver)
    real(r8) :: bebc670tot(pcols,pver), babc670tot(pcols,pver)
    real(r8) :: bebc870tot(pcols,pver), babc870tot(pcols,pver)
    !
    real(r8) :: beoc440tot(pcols,pver), baoc440tot(pcols,pver)
    real(r8) :: beoc500tot(pcols,pver), baoc500tot(pcols,pver)
    real(r8) :: beoc550tot(pcols,pver), baoc550tot(pcols,pver)
    real(r8) :: beoc670tot(pcols,pver), baoc670tot(pcols,pver)
    real(r8) :: beoc870tot(pcols,pver), baoc870tot(pcols,pver)
    !
    real(r8) :: besu440tot(pcols,pver), basu440tot(pcols,pver)
    real(r8) :: besu500tot(pcols,pver), basu500tot(pcols,pver)
    real(r8) :: besu550tot(pcols,pver), basu550tot(pcols,pver)
    real(r8) :: besu670tot(pcols,pver), basu670tot(pcols,pver)
    real(r8) :: besu870tot(pcols,pver), basu870tot(pcols,pver)
    !
    real(r8) :: bebglt1t(pcols,pver), bebggt1t(pcols,pver), bebclt1t(pcols,pver)
    real(r8) :: bebcgt1t(pcols,pver), beoclt1t(pcols,pver), beocgt1t(pcols,pver)
    real(r8) :: bes4lt1t(pcols,pver), bes4gt1t(pcols,pver)
    !
    real(r8) :: be440x(pcols,pver,nbmp1:nmodes), ba440x(pcols,pver,nbmp1:nmodes)
    real(r8) :: be500x(pcols,pver,nbmp1:nmodes), ba500x(pcols,pver,nbmp1:nmodes)
    real(r8) :: be550x(pcols,pver,nbmp1:nmodes), ba550x(pcols,pver,nbmp1:nmodes)
    real(r8) :: be670x(pcols,pver,nbmp1:nmodes), ba670x(pcols,pver,nbmp1:nmodes)
    real(r8) :: be870x(pcols,pver,nbmp1:nmodes), ba870x(pcols,pver,nbmp1:nmodes)
    !
    real(r8) :: belt1x(pcols,pver,nbmp1:nmodes), begt1x(pcols,pver,nbmp1:nmodes)
    !
    real(r8) :: bebc440xt(pcols,pver),babc440xt(pcols,pver)
    real(r8) :: bebc500xt(pcols,pver),babc500xt(pcols,pver)
    real(r8) :: bebc550xt(pcols,pver),babc550xt(pcols,pver)
    real(r8) :: bebc670xt(pcols,pver),babc670xt(pcols,pver)
    real(r8) :: bebc870xt(pcols,pver),babc870xt(pcols,pver)
    !
    real(r8) :: beoc440xt(pcols,pver),baoc440xt(pcols,pver)
    real(r8) :: beoc500xt(pcols,pver),baoc500xt(pcols,pver)
    real(r8) :: beoc550xt(pcols,pver),baoc550xt(pcols,pver)
    real(r8) :: beoc670xt(pcols,pver),baoc670xt(pcols,pver)
    real(r8) :: beoc870xt(pcols,pver),baoc870xt(pcols,pver)
    !
    real(r8) :: bbclt1xt(pcols,pver)
    real(r8) :: bbcgt1xt(pcols,pver), boclt1xt(pcols,pver), bocgt1xt(pcols,pver)
    !
    real(r8) :: bint440du(pcols,pver), bint500du(pcols,pver), bint550du(pcols,pver)
    real(r8) :: bint670du(pcols,pver), bint870du(pcols,pver)
    real(r8) :: bint440ss(pcols,pver), bint500ss(pcols,pver), bint550ss(pcols,pver)
    real(r8) :: bint670ss(pcols,pver), bint870ss(pcols,pver)
    real(r8) :: baint550du(pcols,pver), baint550ss(pcols,pver)
    !
    real(r8) :: bedustlt1(pcols,pver), bedustgt1(pcols,pver)
    real(r8) :: besslt1(pcols,pver), bessgt1(pcols,pver)
    !
    real(r8) :: dod4403d(pcols,pver), abs4403d(pcols,pver)
    real(r8) :: dod4403d_ss(pcols,pver)
    real(r8) :: dod4403d_dust(pcols,pver)
    real(r8) :: dod4403d_so4(pcols,pver)
    real(r8) :: dod4403d_bc(pcols,pver)
    real(r8) :: dod4403d_pom(pcols,pver)
    !
    real(r8) :: dod5003d(pcols,pver), abs5003d(pcols,pver)
    real(r8) :: dod5003d_ss(pcols,pver)
    real(r8) :: dod5003d_dust(pcols,pver)
    real(r8) :: dod5003d_so4(pcols,pver)
    real(r8) :: dod5003d_bc(pcols,pver)
    real(r8) :: dod5003d_pom(pcols,pver)
    !
    real(r8) :: dod5503d(pcols,pver), abs5503d(pcols,pver), abs5503dalt(pcols,pver)
    real(r8) :: dod5503d_ss(pcols,pver), abs5503d_ss(pcols,pver)
    real(r8) :: dod5503d_dust(pcols,pver), abs5503d_dust(pcols,pver)
    real(r8) :: dod5503d_so4(pcols,pver), abs5503d_so4(pcols,pver)
    real(r8) :: dod5503d_bc(pcols,pver), abs5503d_bc(pcols,pver)
    real(r8) :: dod5503d_pom(pcols,pver), abs5503d_pom(pcols,pver)
    !
    real(r8) :: abs6703d(pcols,pver)
    real(r8) :: dod6703d(pcols,pver)
    real(r8) :: dod6703d_ss(pcols,pver)
    real(r8) :: dod6703d_dust(pcols,pver)
    real(r8) :: dod6703d_so4(pcols,pver)
    real(r8) :: dod6703d_bc(pcols,pver)
    real(r8) :: dod6703d_pom(pcols,pver)
    !
    real(r8) :: abs8703d(pcols,pver)
    real(r8) :: dod8703d(pcols,pver)
    real(r8) :: dod8703d_ss(pcols,pver)
    real(r8) :: dod8703d_dust(pcols,pver)
    real(r8) :: dod8703d_so4(pcols,pver)
    real(r8) :: dod8703d_bc(pcols,pver)
    real(r8) :: dod8703d_pom(pcols,pver)
    !
    real(r8) :: dod5503dlt1_ss(pcols,pver)
    real(r8) :: dod5503dgt1_ss(pcols,pver)
    real(r8) :: dod5503dlt1_dust(pcols,pver)
    real(r8) :: dod5503dgt1_dust(pcols,pver)
    real(r8) :: dod5503dlt1_so4(pcols,pver)
    real(r8) :: dod5503dgt1_so4(pcols,pver)
    real(r8) :: dod5503dlt1_bc(pcols,pver)
    real(r8) :: dod5503dgt1_bc(pcols,pver)
    real(r8) :: dod5503dlt1_pom(pcols,pver)
    real(r8) :: dod5503dgt1_pom(pcols,pver)
    !
    real(r8) :: dod500(pcols)
    real(r8) :: dod670(pcols)
    !
    real(r8) :: abs440(pcols)
    real(r8) :: abs500(pcols)
    real(r8) :: abs670(pcols)
    real(r8) :: abs870(pcols)
    !
    real(r8) :: dod440_ss(pcols), dod440_dust(pcols), dod440_so4(pcols)
    real(r8) :: dod440_bc(pcols), dod440_pom(pcols)
    !
    real(r8) :: dod500_ss(pcols), dod500_dust(pcols), dod500_so4(pcols)
    real(r8) :: dod500_bc(pcols), dod500_pom(pcols)
    !
    real(r8) :: dod550_ss(pcols), dod550_dust(pcols), dod550_so4(pcols)
    real(r8) :: dod550_bc(pcols), dod550_pom(pcols)
    !
    real(r8) :: dod670_ss(pcols), dod670_dust(pcols), dod670_so4(pcols)
    real(r8) :: dod670_bc(pcols), dod670_pom(pcols)
    !
    real(r8) :: dod870_ss(pcols), dod870_dust(pcols), dod870_so4(pcols)
    real(r8) :: dod870_bc(pcols), dod870_pom(pcols)
    !
    real(r8) :: dod550lt1_ss(pcols)
    real(r8) :: dod550gt1_ss(pcols)
    real(r8) :: dod550lt1_dust(pcols)
    real(r8) :: dod550gt1_dust(pcols)
    real(r8) :: dod550lt1_so4(pcols)
    real(r8) :: dod550gt1_so4(pcols)
    real(r8) :: dod550lt1_bc(pcols)
    real(r8) :: dod550gt1_bc(pcols)
    real(r8) :: dod550lt1_pom(pcols)
    real(r8) :: dod550gt1_pom(pcols)
    !
    real(r8) :: abs550_ss(pcols)
    real(r8) :: abs550_dust(pcols)
    real(r8) :: abs550_so4(pcols)
    real(r8) :: abs550_bc(pcols)
    real(r8) :: abs550_pom(pcols)
    !
    real(r8) :: vnbcarr(pcols,pver), vaitbcarr(pcols,pver)
    real(r8) :: deltah, airmassl(pcols,pver), airmass(pcols)
    real(r8) :: xrhnull(pcols,pver)
    integer  :: irh1null(pcols,pver)
    integer  :: irf
    !-------------------------------------------------------------------------

    ! Initialize fields
    vaercols(:)     = 0.0_r8
    vaercoll(:)     = 0.0_r8
    aaercols(:)     = 0.0_r8
    aaercoll(:)     = 0.0_r8
    dload(:,:)      = 0.0_r8
    bext550n(:,:,:) = 0.0_r8
    babs550n(:,:,:) = 0.0_r8
    bext440n(:,:,:) = 0.0_r8
    babs440n(:,:,:) = 0.0_r8
    bext870n(:,:,:) = 0.0_r8
    babs870n(:,:,:) = 0.0_r8
    babs500n(:,:,:) = 0.0_r8
    babs670n(:,:,:) = 0.0_r8
    vnbcarr(:,:)    = 0.0_r8
    vaitbcarr(:,:)  = 0.0_r8
    cknorm(:,:,:)   = 0.0_r8

    ! AeroCom diagnostics requiring table look-ups with ambient RH.
    irf = 0
    call opticsAtConstRh(                                  &
         lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irf,  &
         xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,          &
         xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,      &
         xfombg, ifombg1, vnbcarr, vaitbcarr, v_soana,     &
         bext440, bext500, bext550, bext670, bext870,      &
         bebg440, bebg500, bebg550, bebg670, bebg870,      &
         bebc440, bebc500, bebc550, bebc670, bebc870,      &
         beoc440, beoc500, beoc550, beoc670, beoc870,      &
         besu440, besu500, besu550, besu670, besu870,      &
         babs440, babs500, babs550, babs670, babs870,      &
         bebglt1, bebggt1, bebclt1, bebcgt1,               &
         beoclt1, beocgt1, bes4lt1, bes4gt1,               &
         backsc550, babg550, babc550, baoc550, basu550,    &
         bext440n, bext500n, bext550n, bext670n, bext870n, &
         bebg440n, bebg500n, bebg550n, bebg670n, bebg870n, &
         bebc440n, bebc500n, bebc550n, bebc670n, bebc870n, &
         beoc440n, beoc500n, beoc550n, beoc670n, beoc870n, &
         besu440n, besu500n, besu550n, besu670n, besu870n, &
         babs440n, babs500n, babs550n, babs670n, babs870n, &
         bebglt1n, bebggt1n, bebclt1n, bebcgt1n,           &
         beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,           &
         backsc550n, babg550n, babc550n, baoc550n, basu550n)

    do ilev=1,pver
       do icol=1,ncol
          bebglt1t(icol,ilev)  =0.0_r8
          bebggt1t(icol,ilev)  =0.0_r8
          bebclt1t(icol,ilev)  =0.0_r8
          bebcgt1t(icol,ilev)  =0.0_r8
          beoclt1t(icol,ilev)  =0.0_r8
          beocgt1t(icol,ilev)  =0.0_r8
          bes4lt1t(icol,ilev)  =0.0_r8
          bes4gt1t(icol,ilev)  =0.0_r8
          bedustlt1(icol,ilev) =0.0_r8
          bedustgt1(icol,ilev) =0.0_r8
          besslt1(icol,ilev)   =0.0_r8
          bessgt1(icol,ilev)   =0.0_r8

          bext440tot(icol,ilev)=0.0_r8
          babs440tot(icol,ilev)=0.0_r8
          bext500tot(icol,ilev)=0.0_r8
          babs500tot(icol,ilev)=0.0_r8
          bext550tot(icol,ilev)=0.0_r8
          babs550tot(icol,ilev)=0.0_r8
          bext670tot(icol,ilev)=0.0_r8
          babs670tot(icol,ilev)=0.0_r8
          bext870tot(icol,ilev)=0.0_r8
          babs870tot(icol,ilev)=0.0_r8

          backsc550tot(icol,ilev)=0.0_r8

          bebg440tot(icol,ilev)=0.0_r8
          bebg500tot(icol,ilev)=0.0_r8
          bebg550tot(icol,ilev)=0.0_r8
          babg550tot(icol,ilev)=0.0_r8
          bebg670tot(icol,ilev)=0.0_r8
          bebg870tot(icol,ilev)=0.0_r8

          bebc440tot(icol,ilev)=0.0_r8
          bebc500tot(icol,ilev)=0.0_r8
          bebc550tot(icol,ilev)=0.0_r8
          babc550tot(icol,ilev)=0.0_r8
          bebc670tot(icol,ilev)=0.0_r8
          bebc870tot(icol,ilev)=0.0_r8

          beoc440tot(icol,ilev)=0.0_r8
          beoc500tot(icol,ilev)=0.0_r8
          beoc550tot(icol,ilev)=0.0_r8
          baoc550tot(icol,ilev)=0.0_r8
          beoc670tot(icol,ilev)=0.0_r8
          beoc870tot(icol,ilev)=0.0_r8

          besu440tot(icol,ilev)=0.0_r8
          besu500tot(icol,ilev)=0.0_r8
          besu550tot(icol,ilev)=0.0_r8
          basu550tot(icol,ilev)=0.0_r8
          besu670tot(icol,ilev)=0.0_r8
          besu870tot(icol,ilev)=0.0_r8
       enddo
    enddo

    do imode=0,nbmodes
       do ilev=1,pver
          do icol=1,ncol
             ! total internal extinction and absorption for 0.44, 0.50, 0.55, 0.68 and 0.87 um
             bext440tot(icol,ilev)   =bext440tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bext440(icol,ilev,imode)
             babs440tot(icol,ilev)   =babs440tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babs440(icol,ilev,imode)
             bext500tot(icol,ilev)   =bext500tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bext500(icol,ilev,imode)
             babs500tot(icol,ilev)   =babs500tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babs500(icol,ilev,imode)
             bext550tot(icol,ilev)   =bext550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bext550(icol,ilev,imode)
             babs550tot(icol,ilev)   =babs550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babs550(icol,ilev,imode)
             bext670tot(icol,ilev)   =bext670tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bext670(icol,ilev,imode)
             babs670tot(icol,ilev)   =babs670tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babs670(icol,ilev,imode)
             bext870tot(icol,ilev)   =bext870tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bext870(icol,ilev,imode)
             babs870tot(icol,ilev)   =babs870tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babs870(icol,ilev,imode)
             backsc550tot(icol,ilev) =backsc550tot(icol,ilev) +Nnatk(icol,ilev,imode)*backsc550(icol,ilev,imode)
             bebg440tot(icol,ilev)   =bebg440tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bebg440(icol,ilev,imode)
             bebg500tot(icol,ilev)   =bebg500tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bebg500(icol,ilev,imode)
             bebg550tot(icol,ilev)   =bebg550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bebg550(icol,ilev,imode)
             babg550tot(icol,ilev)   =babg550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babg550(icol,ilev,imode)
             bebg670tot(icol,ilev)   =bebg670tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bebg670(icol,ilev,imode)
             bebg870tot(icol,ilev)   =bebg870tot(icol,ilev)   +Nnatk(icol,ilev,imode)*bebg870(icol,ilev,imode)
             besu440tot(icol,ilev)   =besu440tot(icol,ilev)   +Nnatk(icol,ilev,imode)*besu440(icol,ilev,imode)
             besu500tot(icol,ilev)   =besu500tot(icol,ilev)   +Nnatk(icol,ilev,imode)*besu500(icol,ilev,imode)
             besu550tot(icol,ilev)   =besu550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*besu550(icol,ilev,imode)
             basu550tot(icol,ilev)   =basu550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*basu550(icol,ilev,imode)
             besu670tot(icol,ilev)   =besu670tot(icol,ilev)   +Nnatk(icol,ilev,imode)*besu670(icol,ilev,imode)
             besu870tot(icol,ilev)   =besu870tot(icol,ilev)   +Nnatk(icol,ilev,imode)*besu870(icol,ilev,imode)

             ! Condensed OC on modes 1-4 and coagulated BC and OC on modes 5-10:
             if(imode>=1) then
                bebc440tot(icol,ilev)=bebc440tot(icol,ilev)+Nnatk(icol,ilev,imode)*bebc440(icol,ilev,imode)
                bebc500tot(icol,ilev)=bebc500tot(icol,ilev)+Nnatk(icol,ilev,imode)*bebc500(icol,ilev,imode)
                bebc550tot(icol,ilev)=bebc550tot(icol,ilev)+Nnatk(icol,ilev,imode)*bebc550(icol,ilev,imode)
                babc550tot(icol,ilev)=babc550tot(icol,ilev)+Nnatk(icol,ilev,imode)*babc550(icol,ilev,imode)
                bebc670tot(icol,ilev)=bebc670tot(icol,ilev)+Nnatk(icol,ilev,imode)*bebc670(icol,ilev,imode)
                bebc870tot(icol,ilev)=bebc870tot(icol,ilev)+Nnatk(icol,ilev,imode)*bebc870(icol,ilev,imode)
                beoc440tot(icol,ilev)=beoc440tot(icol,ilev)+Nnatk(icol,ilev,imode)*beoc440(icol,ilev,imode)
                beoc500tot(icol,ilev)=beoc500tot(icol,ilev)+Nnatk(icol,ilev,imode)*beoc500(icol,ilev,imode)
                beoc550tot(icol,ilev)=beoc550tot(icol,ilev)+Nnatk(icol,ilev,imode)*beoc550(icol,ilev,imode)
                baoc550tot(icol,ilev)=baoc550tot(icol,ilev)+Nnatk(icol,ilev,imode)*baoc550(icol,ilev,imode)
                beoc670tot(icol,ilev)=beoc670tot(icol,ilev)+Nnatk(icol,ilev,imode)*beoc670(icol,ilev,imode)
                beoc870tot(icol,ilev)=beoc870tot(icol,ilev)+Nnatk(icol,ilev,imode)*beoc870(icol,ilev,imode)
             endif  ! imode>=1
             if(imode==6.or.imode==7) then
                bedustlt1(icol,ilev)=bedustlt1(icol,ilev)+Nnatk(icol,ilev,imode)*bebglt1(icol,ilev,imode)
                bedustgt1(icol,ilev)=bedustgt1(icol,ilev)+Nnatk(icol,ilev,imode)*bebggt1(icol,ilev,imode)
             elseif(imode>=8.and.imode<=10) then
                besslt1(icol,ilev)=besslt1(icol,ilev)+Nnatk(icol,ilev,imode)*bebglt1(icol,ilev,imode)
                bessgt1(icol,ilev)=bessgt1(icol,ilev)+Nnatk(icol,ilev,imode)*bebggt1(icol,ilev,imode)
             endif

             ! Condensed/coagulated SO4 on all modes 1-10, and wet-phase SO4 on modes 4-10:
             bes4lt1t(icol,ilev)=bes4lt1t(icol,ilev)+Nnatk(icol,ilev,imode)*bes4lt1(icol,ilev,imode)
             bes4gt1t(icol,ilev)=bes4gt1t(icol,ilev)+Nnatk(icol,ilev,imode)*bes4gt1(icol,ilev,imode)

             ! Condensed OC on mode 1 and coagulated BC and OC on modes 5-10:
             if(imode>=1) then
                bebclt1t(icol,ilev)=bebclt1t(icol,ilev)+Nnatk(icol,ilev,imode)*bebclt1(icol,ilev,imode)
                bebcgt1t(icol,ilev)=bebcgt1t(icol,ilev)+Nnatk(icol,ilev,imode)*bebcgt1(icol,ilev,imode)
                beoclt1t(icol,ilev)=beoclt1t(icol,ilev)+Nnatk(icol,ilev,imode)*beoclt1(icol,ilev,imode)
                beocgt1t(icol,ilev)=beocgt1t(icol,ilev)+Nnatk(icol,ilev,imode)*beocgt1(icol,ilev,imode)
             endif   ! imode>=1
          end do   ! icol
       enddo     ! ilev
    enddo      ! imode

    ! extinction/absorptions (km-1) for each background component
    ! in the internal mixture are
    do ilev=1,pver
       do icol=1,ncol
          bint440du(icol,ilev)=Nnatk(icol,ilev,6)*bebg440(icol,ilev,6)  +Nnatk(icol,ilev,7) *bebg440(icol,ilev,7)
          bint500du(icol,ilev)=Nnatk(icol,ilev,6)*bebg500(icol,ilev,6)  +Nnatk(icol,ilev,7) *bebg500(icol,ilev,7)
          bint550du(icol,ilev)=Nnatk(icol,ilev,6)*bebg550(icol,ilev,6)  +Nnatk(icol,ilev,7) *bebg550(icol,ilev,7)
          bint670du(icol,ilev)=Nnatk(icol,ilev,6)*bebg670(icol,ilev,6)  +Nnatk(icol,ilev,7) *bebg670(icol,ilev,7)
          bint870du(icol,ilev)=Nnatk(icol,ilev,6)*bebg870(icol,ilev,6)  +Nnatk(icol,ilev,7) *bebg870(icol,ilev,7)
          baint550du(icol,ilev)=Nnatk(icol,ilev,6)*babg550(icol,ilev,6) +Nnatk(icol,ilev,7) *babg550(icol,ilev,7)
          !
          bint440ss(icol,ilev)=Nnatk(icol,ilev,8)*bebg440(icol,ilev,8)  +Nnatk(icol,ilev,9) *bebg440(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*bebg440(icol,ilev,10)
          bint500ss(icol,ilev)=Nnatk(icol,ilev,8)*bebg500(icol,ilev,8)  +Nnatk(icol,ilev,9) *bebg500(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*bebg500(icol,ilev,10)
          bint550ss(icol,ilev)=Nnatk(icol,ilev,8)*bebg550(icol,ilev,8)  +Nnatk(icol,ilev,9) *bebg550(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*bebg550(icol,ilev,10)
          bint670ss(icol,ilev)=Nnatk(icol,ilev,8)*bebg670(icol,ilev,8)  +Nnatk(icol,ilev,9) *bebg670(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*bebg670(icol,ilev,10)
          bint870ss(icol,ilev)=Nnatk(icol,ilev,8)*bebg870(icol,ilev,8)  +Nnatk(icol,ilev,9) *bebg870(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*bebg870(icol,ilev,10)
          baint550ss(icol,ilev)=Nnatk(icol,ilev,8)*babg550(icol,ilev,8) +Nnatk(icol,ilev,9) *babg550(icol,ilev,9) &
                                                               +Nnatk(icol,ilev,10)*babg550(icol,ilev,10)
       end do
    enddo

    do imode=11,14
       do ilev=1,pver
          do icol=1,ncol
             be440x(icol,ilev,imode)=bext440n(icol,ilev,imode-10)
             ba440x(icol,ilev,imode)=babs440n(icol,ilev,imode-10)
             be500x(icol,ilev,imode)=bext500n(icol,ilev,imode-10)
             ba500x(icol,ilev,imode)=babs500n(icol,ilev,imode-10)
             be550x(icol,ilev,imode)=bext550n(icol,ilev,imode-10)
             ba550x(icol,ilev,imode)=babs550n(icol,ilev,imode-10)
             be670x(icol,ilev,imode)=bext670n(icol,ilev,imode-10)
             ba670x(icol,ilev,imode)=babs670n(icol,ilev,imode-10)
             be870x(icol,ilev,imode)=bext870n(icol,ilev,imode-10)
             ba870x(icol,ilev,imode)=babs870n(icol,ilev,imode-10)
             belt1x(icol,ilev,imode)=bebglt1n(icol,ilev,imode-10)
             begt1x(icol,ilev,imode)=bebggt1n(icol,ilev,imode-10)
             backsc550x(icol,ilev,imode)=backsc550n(icol,ilev,imode-10)
          end do
       enddo
    enddo

    ! The externally modes' contribution to extinction and absorption:
    do ilev=1,pver
       do icol=1,ncol

          !BC
          vnbcarr(icol,ilev) = fnbc(icol,ilev)/(fnbc(icol,ilev)+(1.0_r8-fnbc(icol,ilev))*rhopart(l_bc_ni)/rhopart(l_om_ni))
          vnbc = vnbcarr(icol,ilev)
          bebc440xt(icol,ilev) =Nnatk(icol,ilev,12)*be440x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*be440x(icol,ilev,14)
          babc440xt(icol,ilev) =Nnatk(icol,ilev,12)*ba440x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*ba440x(icol,ilev,14)
          bebc500xt(icol,ilev) =Nnatk(icol,ilev,12)*be500x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*be500x(icol,ilev,14)
          babc500xt(icol,ilev) =Nnatk(icol,ilev,12)*ba500x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*ba500x(icol,ilev,14)
          bebc550xt(icol,ilev) =Nnatk(icol,ilev,12)*be550x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*be550x(icol,ilev,14)
          babc550xt(icol,ilev) =Nnatk(icol,ilev,12)*ba550x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*ba550x(icol,ilev,14)
          bebc670xt(icol,ilev) =Nnatk(icol,ilev,12)*be670x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*be670x(icol,ilev,14)
          babc670xt(icol,ilev) =Nnatk(icol,ilev,12)*ba670x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*ba670x(icol,ilev,14)
          bebc870xt(icol,ilev) =Nnatk(icol,ilev,12)*be870x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*be870x(icol,ilev,14)
          babc870xt(icol,ilev) =Nnatk(icol,ilev,12)*ba870x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*ba870x(icol,ilev,14)
          bbclt1xt(icol,ilev)  =Nnatk(icol,ilev,12)*belt1x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*belt1x(icol,ilev,14)
          bbcgt1xt(icol,ilev)  =Nnatk(icol,ilev,12)*begt1x(icol,ilev,12)+vnbc*Nnatk(icol,ilev,14)*begt1x(icol,ilev,14)

          !OC
          beoc440xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*be440x(icol,ilev,14)
          baoc440xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*ba440x(icol,ilev,14)
          beoc500xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*be500x(icol,ilev,14)
          baoc500xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*ba500x(icol,ilev,14)
          beoc550xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*be550x(icol,ilev,14)
          baoc550xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*ba550x(icol,ilev,14)
          beoc670xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*be670x(icol,ilev,14)
          baoc670xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*ba670x(icol,ilev,14)
          beoc870xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*be870x(icol,ilev,14)
          baoc870xt(icol,ilev) =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*ba870x(icol,ilev,14)
          boclt1xt(icol,ilev)  =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*belt1x(icol,ilev,14)
          bocgt1xt(icol,ilev)  =(1.0_r8-vnbc)*Nnatk(icol,ilev,14)*begt1x(icol,ilev,14)

          ! Total (for all modes) absorption optical depth and backscattering
          abs550_aer(icol,ilev)=babs550tot(icol,ilev)  &
               +Nnatk(icol,ilev,12)*ba550x(icol,ilev,12)+Nnatk(icol,ilev,14)*ba550x(icol,ilev,14)
          abs550_aer(icol,ilev)=1.e-3_r8*abs550_aer(icol,ilev)

          bs550_aer(icol,ilev)= backsc550tot(icol,ilev)   &
               +Nnatk(icol,ilev,12)*backsc550x(icol,ilev,12)+Nnatk(icol,ilev,14)*backsc550x(icol,ilev,14)
          bs550_aer(icol,ilev)=1.e-3_r8*bs550_aer(icol,ilev)
          !
       end do
    enddo

    ! collect AeroCom-fields for optical depth/absorption of each comp,
    ! 3D and 2D, at 440, 500, 550, 670 and 870 nm, for all d, d<1um and d>1um
    ! initialize 2d-fields
    do icol=1,ncol
       dod440(icol) = 0.0_r8
       abs440(icol) = 0.0_r8
       dod500(icol) = 0.0_r8
       abs500(icol) = 0.0_r8
       dod550(icol) = 0.0_r8
       abs550(icol) = 0.0_r8
       abs550alt(icol) = 0.0_r8
       dod670(icol) = 0.0_r8
       abs670(icol) = 0.0_r8
       dod870(icol) = 0.0_r8
       abs870(icol) = 0.0_r8
       !
       abs550_ss(icol)   = 0.0_r8
       abs550_dust(icol) = 0.0_r8
       abs550_so4(icol)  = 0.0_r8
       abs550_bc(icol)   = 0.0_r8
       abs550_pom(icol)  = 0.0_r8
       !
       dod440_ss(icol)           = 0.0_r8
       dod440_dust(icol)         = 0.0_r8
       dod440_so4(icol)          = 0.0_r8
       dod440_bc(icol)           = 0.0_r8
       dod440_pom(icol)          = 0.0_r8
       dod500_ss(icol)           = 0.0_r8
       dod500_dust(icol)         = 0.0_r8
       dod500_so4(icol)          = 0.0_r8
       dod500_bc(icol)           = 0.0_r8
       dod500_pom(icol)          = 0.0_r8
       dod550_ss(icol)           = 0.0_r8
       dod550_dust(icol)         = 0.0_r8
       dod550_so4(icol)          = 0.0_r8
       dod550_bc(icol)           = 0.0_r8
       dod550_pom(icol)          = 0.0_r8
       dod670_ss(icol)           = 0.0_r8
       dod670_dust(icol)         = 0.0_r8
       dod670_so4(icol)          = 0.0_r8
       dod670_bc(icol)           = 0.0_r8
       dod670_pom(icol)          = 0.0_r8
       dod870_ss(icol)           = 0.0_r8
       dod870_dust(icol)         = 0.0_r8
       dod870_so4(icol)          = 0.0_r8
       dod870_bc(icol)           = 0.0_r8
       dod870_pom(icol)          = 0.0_r8
       dod550lt1_ss(icol)        = 0.0_r8
       dod550gt1_ss(icol)        = 0.0_r8
       dod550lt1_dust(icol)      = 0.0_r8
       dod550gt1_dust(icol)      = 0.0_r8
       dod550lt1_so4(icol)       = 0.0_r8
       dod550gt1_so4(icol)       = 0.0_r8
       dod550lt1_bc(icol)        = 0.0_r8
       dod550gt1_bc(icol)        = 0.0_r8
       dod550lt1_pom(icol)       = 0.0_r8
       dod550gt1_pom(icol)       = 0.0_r8
       do ilev                   =1,pver
          abs4403d(icol,ilev)    = 0.0_r8
          abs5003d(icol,ilev)    = 0.0_r8
          abs5503d(icol,ilev)    = 0.0_r8
          abs6703d(icol,ilev)    = 0.0_r8
          abs8703d(icol,ilev)    = 0.0_r8
          abs5503dalt(icol,ilev) = 0.0_r8
       enddo
    enddo

    do icol=1,ncol
       do ilev=1,pver
          ! Layer thickness, unit km
          deltah=deltah_km(icol,ilev)

          ! 3D optical depths for monthly averages

          !SS
          dod4403d_ss(icol,ilev) = bint440ss(icol,ilev)*deltah
          dod5003d_ss(icol,ilev) = bint500ss(icol,ilev)*deltah
          dod5503d_ss(icol,ilev) = bint550ss(icol,ilev)*deltah
          abs5503d_ss(icol,ilev) = baint550ss(icol,ilev)*deltah
          dod6703d_ss(icol,ilev) = bint670ss(icol,ilev)*deltah
          dod8703d_ss(icol,ilev) = bint870ss(icol,ilev)*deltah

          !DUST
          dod4403d_dust(icol,ilev) = bint440du(icol,ilev)*deltah
          dod5003d_dust(icol,ilev) = bint500du(icol,ilev)*deltah
          dod5503d_dust(icol,ilev) = bint550du(icol,ilev)*deltah
          abs5503d_dust(icol,ilev) = baint550du(icol,ilev)*deltah
          dod6703d_dust(icol,ilev) = bint670du(icol,ilev)*deltah
          dod8703d_dust(icol,ilev) = bint870du(icol,ilev)*deltah

          !SO4 soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          dod4403d_so4(icol,ilev) = (besu440tot(icol,ilev)                          & ! condensate )
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*bebg440(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebg440(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          dod5003d_so4(icol,ilev) = (besu500tot(icol,ilev)                          & ! condensate
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*bebg500(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebg500(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          dod5503d_so4(icol,ilev) = (besu550tot(icol,ilev)                          & ! condensate
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*bebg550(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebg550(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          abs5503d_so4(icol,ilev) = (basu550tot(icol,ilev)                          & ! condensate )
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*babg550(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*babg550(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          dod6703d_so4(icol,ilev) = (besu670tot(icol,ilev)                          & ! condensate
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*bebg670(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebg670(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          dod8703d_so4(icol,ilev) = (besu870tot(icol,ilev)                          & ! condensate
               +(1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*bebg870(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebg870(icol,ilev,5))*deltah                      ! background, SO4(Ait75) mode (5)
          !BC
          vaitbcarr(icol,ilev) = faitbc(icol,ilev)/(faitbc(icol,ilev)          &
               +(1.0_r8-faitbc(icol,ilev))*rhopart(l_bc_ni)/rhopart(l_om_ni))
          vaitbc = vaitbcarr(icol,ilev)
          dod4403d_bc(icol,ilev) = (bebc440tot(icol,ilev)+bebc440xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebg440(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebg440(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebg440(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)
          dod5003d_bc(icol,ilev) = (bebc500tot(icol,ilev)+bebc500xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebg500(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebg500(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebg500(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)
          dod5503d_bc(icol,ilev) = (bebc550tot(icol,ilev)+bebc550xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebg550(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebg550(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebg550(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)
          abs5503d_bc(icol,ilev) = (babc550tot(icol,ilev)+babc550xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*babg550(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*babg550(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*babg550(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)
          dod6703d_bc(icol,ilev) = (bebc670tot(icol,ilev)+bebc670xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebg670(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebg670(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebg670(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)
          dod8703d_bc(icol,ilev) = (bebc870tot(icol,ilev)+bebc870xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebg870(icol,ilev,2)                       & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebg870(icol,ilev,4)                & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebg870(icol,ilev,0))*deltah                 ! background, BC(ax) mode (0)

          !OC soa + v_soana part of mode 11 for the OC volume fraction of that mode v_soana(icol,ilev)
          dod4403d_pom(icol,ilev) = (beoc440tot(icol,ilev)+beoc440xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebg440(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebg440(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)
          dod5003d_pom(icol,ilev) = (beoc500tot(icol,ilev)+beoc500xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebg500(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebg500(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)
          dod5503d_pom(icol,ilev) = (beoc550tot(icol,ilev)+beoc550xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebg550(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebg550(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)
          abs5503d_pom(icol,ilev) = (baoc550tot(icol,ilev)+baoc550xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*babg550(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*babg550(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)
          dod6703d_pom(icol,ilev) = (beoc670tot(icol,ilev)+beoc670xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebg670(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebg670(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)
          dod8703d_pom(icol,ilev) = (beoc870tot(icol,ilev)+beoc870xt(icol,ilev)   & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebg870(icol,ilev,1)*v_soana(icol,ilev)       & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebg870(icol,ilev,4))*deltah    ! background in OC &BC(Ait) mode (4)

          ec550_so4(icol,ilev) = 1.e-3*dod5503d_so4(icol,ilev)/deltah
          ec550_bc(icol,ilev)  = 1.e-3*dod5503d_bc(icol,ilev)/deltah
          ec550_pom(icol,ilev) = 1.e-3*dod5503d_pom(icol,ilev)/deltah
          ec550_ss(icol,ilev)  = 1.e-3*dod5503d_ss(icol,ilev)/deltah
          ec550_du(icol,ilev)  = 1.e-3*dod5503d_dust(icol,ilev)/deltah
          ec550_aer(icol,ilev) = ec550_so4(icol,ilev) + ec550_bc(icol,ilev)&
               + ec550_pom(icol,ilev) + ec550_ss(icol,ilev) + ec550_du(icol,ilev)

          ! Total 3D optical depths/abs. for column integrations
          dod4403d(icol,ilev) = dod4403d_ss(icol,ilev)+dod4403d_dust(icol,ilev) &
               +dod4403d_so4(icol,ilev)+dod4403d_bc(icol,ilev)                  &
               +dod4403d_pom(icol,ilev)
          dod5003d(icol,ilev) = dod5003d_ss(icol,ilev)+dod5003d_dust(icol,ilev) &
               +dod5003d_so4(icol,ilev)+dod5003d_bc(icol,ilev)                  &
               +dod5003d_pom(icol,ilev)
          dod5503d(icol,ilev) = dod5503d_ss(icol,ilev)+dod5503d_dust(icol,ilev) &
               +dod5503d_so4(icol,ilev)+dod5503d_bc(icol,ilev)                  &
               +dod5503d_pom(icol,ilev)
          dod6703d(icol,ilev) = dod6703d_ss(icol,ilev)+dod6703d_dust(icol,ilev) &
               +dod6703d_so4(icol,ilev)+dod6703d_bc(icol,ilev)                  &
               +dod6703d_pom(icol,ilev)
          dod8703d(icol,ilev) = dod8703d_ss(icol,ilev)+dod8703d_dust(icol,ilev) &
               +dod8703d_so4(icol,ilev)+dod8703d_bc(icol,ilev)                  &
               +dod8703d_pom(icol,ilev)
          abs5503d(icol,ilev) = abs5503d_ss(icol,ilev)+abs5503d_dust(icol,ilev) &
               +abs5503d_so4(icol,ilev)+abs5503d_bc(icol,ilev)                  &
               +abs5503d_pom(icol,ilev)

          ! (Note: Local abs550alt is up to 6% larger (annually averaged) in typical b.b.
          ! regions, compared to abs550. This is most likely most correct, but should be checked!)
          do imode=0,10
             abs4403d(icol,ilev) = abs4403d(icol,ilev)+Nnatk(icol,ilev,imode)*babs440(icol,ilev,imode)*deltah
             abs5003d(icol,ilev) = abs5003d(icol,ilev)+Nnatk(icol,ilev,imode)*babs500(icol,ilev,imode)*deltah
             abs6703d(icol,ilev) = abs6703d(icol,ilev)+Nnatk(icol,ilev,imode)*babs670(icol,ilev,imode)*deltah
             abs8703d(icol,ilev) = abs8703d(icol,ilev)+Nnatk(icol,ilev,imode)*babs870(icol,ilev,imode)*deltah
             abs5503dalt(icol,ilev) = abs5503dalt(icol,ilev)+Nnatk(icol,ilev,imode)*babs550(icol,ilev,imode)*deltah
          enddo
          do imode=11,14
             abs4403d(icol,ilev) = abs4403d(icol,ilev)+Nnatk(icol,ilev,imode)*babs440n(icol,ilev,imode-10)*deltah
             abs5003d(icol,ilev) = abs5003d(icol,ilev)+Nnatk(icol,ilev,imode)*babs500n(icol,ilev,imode-10)*deltah
             abs6703d(icol,ilev) = abs6703d(icol,ilev)+Nnatk(icol,ilev,imode)*babs670n(icol,ilev,imode-10)*deltah
             abs8703d(icol,ilev) = abs8703d(icol,ilev)+Nnatk(icol,ilev,imode)*babs870n(icol,ilev,imode-10)*deltah
             abs5503dalt(icol,ilev) = abs5503dalt(icol,ilev)+Nnatk(icol,ilev,imode)*babs550n(icol,ilev,imode-10)*deltah
          enddo

          ! optical depths for d<1um and d>1um (r<0.5um and r>0.5um)

          !SS
          dod5503dlt1_ss(icol,ilev) = besslt1(icol,ilev)*deltah
          dod5503dgt1_ss(icol,ilev) = bessgt1(icol,ilev)*deltah

          !DUST
          dod5503dlt1_dust(icol,ilev) = bedustlt1(icol,ilev)*deltah
          dod5503dgt1_dust(icol,ilev) = bedustgt1(icol,ilev)*deltah

          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          dod5503dlt1_so4(icol,ilev) = (bes4lt1t(icol,ilev)                          & ! condensate
               + Nnatk(icol,ilev,1)*bebglt1(icol,ilev,1)*(1.0_r8-v_soana(icol,ilev)) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebglt1(icol,ilev,5))*deltah                       ! background, SO4(Ait75) mode (5)

          dod5503dgt1_so4(icol,ilev) = (bes4gt1t(icol,ilev)                          & ! condensate + n-mode (11)
               + Nnatk(icol,ilev,1)*bebggt1(icol,ilev,1)*(1.0_r8-v_soana(icol,ilev)) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebggt1(icol,ilev,5))*deltah                       ! background, SO4(Ait75) mode (5)
          !BC
          dod5503dlt1_bc(icol,ilev) =  (bebclt1t(icol,ilev)+bbclt1xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebglt1(icol,ilev,2)                        & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebglt1(icol,ilev,4)                 & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebglt1(icol,ilev,0))*deltah                  ! background, BC(ax) mode (0)

          dod5503dgt1_bc(icol,ilev) =  (bebcgt1t(icol,ilev)+bbcgt1xt(icol,ilev) & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebggt1(icol,ilev,2)                        & ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,ilev,4)*bebggt1(icol,ilev,4)                 & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebggt1(icol,ilev,0))*deltah                  ! background, BC(ax) mode (0)
          !OC
          !soa + v_soana part of mode 11 for the OC volume fraction of that mode
          dod5503dlt1_pom(icol,ilev) = (beoclt1t(icol,ilev)+boclt1xt(icol,ilev)  & ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,ilev,1)*bebglt1(icol,ilev,1)*v_soana(icol,ilev)      & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebglt1(icol,ilev,4))*deltah   ! background in OC &BC(Ait) mode (4)

          dod5503dgt1_pom(icol,ilev) = (beocgt1t(icol,ilev)+bocgt1xt(icol,ilev)  & ! coagulated + n-mode OC&OC (14)
               + Nnatk(icol,ilev,1)*bebggt1(icol,ilev,1)*v_soana(icol,ilev)      & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,ilev,4)*bebggt1(icol,ilev,4))*deltah   ! background in OC &BC(Ait) mode (4)

          ! Column integrated optical depths/abs., total and for each constituent
          dod440(icol)         = dod440(icol)    +dod4403d(icol,ilev)
          abs440(icol)         = abs440(icol)    +abs4403d(icol,ilev)
          dod500(icol)         = dod500(icol)    +dod5003d(icol,ilev)
          abs500(icol)         = abs500(icol)    +abs5003d(icol,ilev)
          dod550(icol)         = dod550(icol)    +dod5503d(icol,ilev)
          abs550(icol)         = abs550(icol)    +abs5503d(icol,ilev)
          abs550alt(icol)      = abs550alt(icol) +abs5503dalt(icol,ilev)
          dod670(icol)         = dod670(icol)    +dod6703d(icol,ilev)
          abs670(icol)         = abs670(icol)    +abs6703d(icol,ilev)
          dod870(icol)         = dod870(icol)    +dod8703d(icol,ilev)
          abs870(icol)         = abs870(icol)    +abs8703d(icol,ilev)

          ! Added abs components
          abs550_ss(icol)      = abs550_ss(icol)      +abs5503d_ss(icol,ilev)
          abs550_dust(icol)    = abs550_dust(icol)    +abs5503d_dust(icol,ilev)
          abs550_so4(icol)     = abs550_so4(icol)     +abs5503d_so4(icol,ilev)
          abs550_bc(icol)      = abs550_bc(icol)      +abs5503d_bc(icol,ilev)
          abs550_pom(icol)     = abs550_pom(icol)     +abs5503d_pom(icol,ilev)
          !
          dod440_ss(icol)      = dod440_ss(icol)      +dod4403d_ss(icol,ilev)
          dod440_dust(icol)    = dod440_dust(icol)    +dod4403d_dust(icol,ilev)
          dod440_so4(icol)     = dod440_so4(icol)     +dod4403d_so4(icol,ilev)
          dod440_bc(icol)      = dod440_bc(icol)      +dod4403d_bc(icol,ilev)
          dod440_pom(icol)     = dod440_pom(icol)     +dod4403d_pom(icol,ilev)
          dod500_ss(icol)      = dod500_ss(icol)      +dod5003d_ss(icol,ilev)
          dod500_dust(icol)    = dod500_dust(icol)    +dod5003d_dust(icol,ilev)
          dod500_so4(icol)     = dod500_so4(icol)     +dod5003d_so4(icol,ilev)
          dod500_bc(icol)      = dod500_bc(icol)      +dod5003d_bc(icol,ilev)
          dod500_pom(icol)     = dod500_pom(icol)     +dod5003d_pom(icol,ilev)
          dod550_ss(icol)      = dod550_ss(icol)      +dod5503d_ss(icol,ilev)
          dod550_dust(icol)    = dod550_dust(icol)    +dod5503d_dust(icol,ilev)
          dod550_so4(icol)     = dod550_so4(icol)     +dod5503d_so4(icol,ilev)
          dod550_bc(icol)      = dod550_bc(icol)      +dod5503d_bc(icol,ilev)
          dod550_pom(icol)     = dod550_pom(icol)     +dod5503d_pom(icol,ilev)
          dod670_ss(icol)      = dod670_ss(icol)      +dod6703d_ss(icol,ilev)
          dod670_dust(icol)    = dod670_dust(icol)    +dod6703d_dust(icol,ilev)
          dod670_so4(icol)     = dod670_so4(icol)     +dod6703d_so4(icol,ilev)
          dod670_bc(icol)      = dod670_bc(icol)      +dod6703d_bc(icol,ilev)
          dod670_pom(icol)     = dod670_pom(icol)     +dod6703d_pom(icol,ilev)
          dod870_ss(icol)      = dod870_ss(icol)      +dod8703d_ss(icol,ilev)
          dod870_dust(icol)    = dod870_dust(icol)    +dod8703d_dust(icol,ilev)
          dod870_so4(icol)     = dod870_so4(icol)     +dod8703d_so4(icol,ilev)
          dod870_bc(icol)      = dod870_bc(icol)      +dod8703d_bc(icol,ilev)
          dod870_pom(icol)     = dod870_pom(icol)     +dod8703d_pom(icol,ilev)
          dod550lt1_ss(icol)   = dod550lt1_ss(icol)   +dod5503dlt1_ss(icol,ilev)
          dod550gt1_ss(icol)   = dod550gt1_ss(icol)   +dod5503dgt1_ss(icol,ilev)
          dod550lt1_dust(icol) = dod550lt1_dust(icol) +dod5503dlt1_dust(icol,ilev)
          dod550gt1_dust(icol) = dod550gt1_dust(icol) +dod5503dgt1_dust(icol,ilev)
          dod550lt1_so4(icol)  = dod550lt1_so4(icol)  +dod5503dlt1_so4(icol,ilev)
          dod550gt1_so4(icol)  = dod550gt1_so4(icol)  +dod5503dgt1_so4(icol,ilev)
          dod550lt1_bc(icol)   = dod550lt1_bc(icol)   +dod5503dlt1_bc(icol,ilev)
          dod550gt1_bc(icol)   = dod550gt1_bc(icol)   +dod5503dgt1_bc(icol,ilev)
          dod550lt1_pom(icol)  = dod550lt1_pom(icol)  +dod5503dlt1_pom(icol,ilev)
          dod550gt1_pom(icol)  = dod550gt1_pom(icol)  +dod5503dgt1_pom(icol,ilev)
       enddo ! ilev
    enddo  ! icol

    ! extinction, absorption (m-1) and backscatter coefficients (m-1 sr-1)
    call outfld('EC550AER',ec550_aer,pcols,lchnk)
    call outfld('ABS550_A',abs550_aer,pcols,lchnk)
    call outfld('BS550AER',bs550_aer,pcols,lchnk)
    !
    ! speciated extinction coefficients (m-1)
    call outfld('EC550SO4',ec550_so4,pcols,lchnk)
    call outfld('EC550BC ',ec550_bc ,pcols,lchnk)
    call outfld('EC550POM',ec550_pom,pcols,lchnk)
    call outfld('EC550SS ',ec550_ss ,pcols,lchnk)
    call outfld('EC550DU ',ec550_du ,pcols,lchnk)
    !
    ! optical depths and absorption as requested by AeroCom
    ! notation: 3=3D, D=DOD, A=ABS, LT=d<1um, GT=d>1um
    call outfld('DOD440  ',dod440         ,pcols,lchnk)
    call outfld('ABS440  ',abs440         ,pcols,lchnk)
    call outfld('DOD500  ',dod500         ,pcols,lchnk)
    call outfld('ABS500  ',abs500         ,pcols,lchnk)
    call outfld('DOD550  ',dod550         ,pcols,lchnk)
    call outfld('ABS550  ',abs550         ,pcols,lchnk)
    call outfld('ABS550AL',abs550alt      ,pcols,lchnk)
    call outfld('DOD670  ',dod670         ,pcols,lchnk)
    call outfld('ABS670  ',abs670         ,pcols,lchnk)
    call outfld('DOD870  ',dod870         ,pcols,lchnk)
    call outfld('ABS870  ',abs870         ,pcols,lchnk)
    call outfld('A550_SS ',abs550_ss      ,pcols,lchnk)
    call outfld('A550_DU ',abs550_dust    ,pcols,lchnk)
    call outfld('A550_SO4',abs550_so4     ,pcols,lchnk)
    call outfld('A550_BC ',abs550_bc      ,pcols,lchnk)
    call outfld('A550_POM',abs550_pom     ,pcols,lchnk)
    !
    call outfld('D440_SS ',dod440_ss      ,pcols,lchnk)
    call outfld('D440_DU ',dod440_dust    ,pcols,lchnk)
    call outfld('D440_SO4',dod440_so4     ,pcols,lchnk)
    call outfld('D440_BC ',dod440_bc      ,pcols,lchnk)
    call outfld('D440_POM',dod440_pom     ,pcols,lchnk)
    call outfld('D500_SS ',dod500_ss      ,pcols,lchnk)
    call outfld('D500_DU ',dod500_dust    ,pcols,lchnk)
    call outfld('D500_SO4',dod500_so4     ,pcols,lchnk)
    call outfld('D500_BC ',dod500_bc      ,pcols,lchnk)
    call outfld('D500_POM',dod500_pom     ,pcols,lchnk)
    call outfld('D550_SS ',dod550_ss      ,pcols,lchnk)
    call outfld('D550_DU ',dod550_dust    ,pcols,lchnk)
    call outfld('D550_SO4',dod550_so4     ,pcols,lchnk)
    call outfld('D550_BC ',dod550_bc      ,pcols,lchnk)
    call outfld('D550_POM',dod550_pom     ,pcols,lchnk)
    call outfld('D670_SS ',dod670_ss      ,pcols,lchnk)
    call outfld('D670_DU ',dod670_dust    ,pcols,lchnk)
    call outfld('D670_SO4',dod670_so4     ,pcols,lchnk)
    call outfld('D670_BC ',dod670_bc      ,pcols,lchnk)
    call outfld('D670_POM',dod670_pom     ,pcols,lchnk)
    call outfld('D870_SS ',dod870_ss      ,pcols,lchnk)
    call outfld('D870_DU ',dod870_dust    ,pcols,lchnk)
    call outfld('D870_SO4',dod870_so4     ,pcols,lchnk)
    call outfld('D870_BC ',dod870_bc      ,pcols,lchnk)
    call outfld('D870_POM',dod870_pom     ,pcols,lchnk)
    call outfld('DLT_SS  ',dod550lt1_ss   ,pcols,lchnk)
    call outfld('DGT_SS  ',dod550gt1_ss   ,pcols,lchnk)
    call outfld('DLT_DUST',dod550lt1_dust ,pcols,lchnk)
    call outfld('DGT_DUST',dod550gt1_dust ,pcols,lchnk)
    call outfld('DLT_SO4 ',dod550lt1_so4  ,pcols,lchnk)
    call outfld('DGT_SO4 ',dod550gt1_so4  ,pcols,lchnk)
    call outfld('DLT_BC  ',dod550lt1_bc   ,pcols,lchnk)
    call outfld('DGT_BC  ',dod550gt1_bc   ,pcols,lchnk)
    call outfld('DLT_POM ',dod550lt1_pom  ,pcols,lchnk)
    call outfld('DGT_POM ',dod550gt1_pom  ,pcols,lchnk)

    ! Dry parameters of each aerosol component
    ! BC(ax) mode
    call intdrypar0(lchnk, ncol, Nnatk,                          &
         cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, &
         cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
         cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
         cknorm,cknlt05,ckngt125)

    ! SO4&SOA(Ait,n) mode
    call intdrypar1(lchnk, ncol, Nnatk, xfombg, ifombg1,         &
         xct, ict1, xfac, ifac1,                                   &
         cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, &
         cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
         cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
         aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)

    ! BC(Ait,n) and OC(Ait,n) modes
    call intdrypar2to3(lchnk, ncol, Nnatk, xct, ict1, xfac, ifac1, &
         cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,   &
         cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,   &
         cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
         aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)

    ! BC&OC(Ait,n) mode   ------ fcm not valid here (=0). Use faitbc or fnbc instead
    call intdrypar4(lchnk, ncol, Nnatk,                          &
         xfbcbg, ifbcbg1, xfbcbgn, ifbcbgn1,                       &
         xct, ict1, xfac, ifac1, xfaq, ifaq1,                      &
         cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, &
         cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
         cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol, &
         aaerosn,aaeroln,vaerosn,vaeroln,cknorm,cknlt05,ckngt125)

    ! SO4(Ait75) (5), mineral (6-7) and Sea-salt (8-10) modes:
    call intdrypar5to10(lchnk, ncol, Nnatk,                      &
         xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1,         &
         cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125, &
         cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125, &
         cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol,&
         cknorm,cknlt05,ckngt125)

    do ilev=1,pver
       do icol=1,ncol
          c_ss(icol,ilev)=0.0_r8
          c_mi(icol,ilev)=0.0_r8
       enddo
    enddo

    do ilev=1,pver
       do icol=1,ncol
          ! mineral and sea-salt background concentrations, internally mixed
          c_mi(icol,ilev)    = Nnatk(icol,ilev,6)*cintbg(icol,ilev,6)    +Nnatk(icol,ilev,7)*cintbg(icol,ilev,7)
          c_mi05(icol,ilev)  = Nnatk(icol,ilev,6)*cintbg05(icol,ilev,6)  +Nnatk(icol,ilev,7)*cintbg05(icol,ilev,7)
          c_mi125(icol,ilev) = Nnatk(icol,ilev,6)*cintbg125(icol,ilev,6) +Nnatk(icol,ilev,7)*cintbg125(icol,ilev,7)
          c_ss(icol,ilev)    = Nnatk(icol,ilev,8)*cintbg(icol,ilev,8)    +Nnatk(icol,ilev,9)*cintbg(icol,ilev,9)    &
                                                                         +Nnatk(icol,ilev,10)*cintbg(icol,ilev,10)
          c_ss05(icol,ilev)  = Nnatk(icol,ilev,8)*cintbg05(icol,ilev,8)  +Nnatk(icol,ilev,9)*cintbg05(icol,ilev,9)  &
                                                                         +Nnatk(icol,ilev,10)*cintbg05(icol,ilev,10)
          c_ss125(icol,ilev) = Nnatk(icol,ilev,8)*cintbg125(icol,ilev,8) +Nnatk(icol,ilev,9)*cintbg125(icol,ilev,9) &
                                                                         +Nnatk(icol,ilev,10)*cintbg125(icol,ilev,10)

          ! internally mixed bc and oc (from coagulation) and so4 concentrations
          ! (sa=so4(aq) and sc=so4(cond+coag), separated because of different density:
          ! necessary for calculation of volume fractions!), and total aerosol surface
          ! areas and volumes.
          c_bc(icol,ilev)       =0.0_r8
          c_bc05(icol,ilev)     =0.0_r8
          c_bc125(icol,ilev)    =0.0_r8
          c_oc(icol,ilev)       =0.0_r8
          c_oc05(icol,ilev)     =0.0_r8
          c_oc125(icol,ilev)    =0.0_r8
          c_s4(icol,ilev)       =0.0_r8
          c_s4_a(icol,ilev)     =0.0_r8
          c_s4_1(icol,ilev)     =0.0_r8
          c_s4_5(icol,ilev)     =0.0_r8
          c_sa(icol,ilev)       =0.0_r8
          c_sa05(icol,ilev)     =0.0_r8
          c_sa125(icol,ilev)    =0.0_r8
          c_sc(icol,ilev)       =0.0_r8
          c_sc05(icol,ilev)     =0.0_r8
          c_sc125(icol,ilev)    =0.0_r8
          aaeros_tot(icol,ilev) =0.0_r8
          aaerol_tot(icol,ilev) =0.0_r8
          vaeros_tot(icol,ilev) =0.0_r8
          vaerol_tot(icol,ilev) =0.0_r8
          c_bc_0(icol,ilev)     =0.0_r8
          c_bc_2(icol,ilev)     =0.0_r8
          c_bc_4(icol,ilev)     =0.0_r8
          c_bc_12(icol,ilev)    =0.0_r8
          c_bc_14(icol,ilev)    =0.0_r8
          c_oc_4(icol,ilev)     =0.0_r8
          c_oc_14(icol,ilev)    =0.0_r8

          c_tot(icol,ilev)    =0.0_r8
          c_tot125(icol,ilev) =0.0_r8
          c_tot05(icol,ilev)  =0.0_r8
          c_pm25(icol,ilev)   =0.0_r8
          c_pm1(icol,ilev)    =0.0_r8
          mmr_pm25(icol,ilev) =0.0_r8
          mmr_pm1(icol,ilev)  =0.0_r8

          do imode=0,nbmodes
             if(imode.ne.3) then
                c_bc(icol,ilev)       = c_bc(icol,ilev)       +Nnatk(icol,ilev,imode)*cintbc(icol,ilev,imode)
                c_bc05(icol,ilev)     = c_bc05(icol,ilev)     +Nnatk(icol,ilev,imode)*cintbc05(icol,ilev,imode)
                c_bc125(icol,ilev)    = c_bc125(icol,ilev)    +Nnatk(icol,ilev,imode)*cintbc125(icol,ilev,imode)
                c_oc(icol,ilev)       = c_oc(icol,ilev)       +Nnatk(icol,ilev,imode)*cintoc(icol,ilev,imode)
                c_oc05(icol,ilev)     = c_oc05(icol,ilev)     +Nnatk(icol,ilev,imode)*cintoc05(icol,ilev,imode)
                c_oc125(icol,ilev)    = c_oc125(icol,ilev)    +Nnatk(icol,ilev,imode)*cintoc125(icol,ilev,imode)
                c_sa(icol,ilev)       = c_sa(icol,ilev)       +Nnatk(icol,ilev,imode)*cintsa(icol,ilev,imode)
                c_sa05(icol,ilev)     = c_sa05(icol,ilev)     +Nnatk(icol,ilev,imode)*cintsa05(icol,ilev,imode)
                c_sa125(icol,ilev)    = c_sa125(icol,ilev)    +Nnatk(icol,ilev,imode)*cintsa125(icol,ilev,imode)
                c_sc(icol,ilev)       = c_sc(icol,ilev)       +Nnatk(icol,ilev,imode)*cintsc(icol,ilev,imode)
                c_sc05(icol,ilev)     = c_sc05(icol,ilev)     +Nnatk(icol,ilev,imode)*cintsc05(icol,ilev,imode)
                c_sc125(icol,ilev)    = c_sc125(icol,ilev)    +Nnatk(icol,ilev,imode)*cintsc125(icol,ilev,imode)
                aaeros_tot(icol,ilev) = aaeros_tot(icol,ilev) +Nnatk(icol,ilev,imode)*aaeros(icol,ilev,imode)
                aaerol_tot(icol,ilev) = aaerol_tot(icol,ilev) +Nnatk(icol,ilev,imode)*aaerol(icol,ilev,imode)
                vaeros_tot(icol,ilev) = vaeros_tot(icol,ilev) +Nnatk(icol,ilev,imode)*vaeros(icol,ilev,imode)
                vaerol_tot(icol,ilev) = vaerol_tot(icol,ilev) +Nnatk(icol,ilev,imode)*vaerol(icol,ilev,imode)
             endif
          enddo
          ! add dry aerosol area and volume of externally mixed modes
          do imode=nbmp1,nmodes
             aaeros_tot(icol,ilev) = aaeros_tot(icol,ilev) +Nnatk(icol,ilev,imode)*aaerosn(icol,ilev,imode)
             aaerol_tot(icol,ilev) = aaerol_tot(icol,ilev) +Nnatk(icol,ilev,imode)*aaeroln(icol,ilev,imode)
             vaeros_tot(icol,ilev) = vaeros_tot(icol,ilev) +Nnatk(icol,ilev,imode)*vaerosn(icol,ilev,imode)
             vaerol_tot(icol,ilev) = vaerol_tot(icol,ilev) +Nnatk(icol,ilev,imode)*vaeroln(icol,ilev,imode)
          end do

          ! Effective radii for particles smaller and greater than 0.5um,
          ! and for all radii, in each layer (er=3*V/A):
          erlt053d(icol,ilev)=3.0_r8*vaeros_tot(icol,ilev)/(aaeros_tot(icol,ilev)+eps)
          ergt053d(icol,ilev)=3.0_r8*vaerol_tot(icol,ilev)/(aaerol_tot(icol,ilev)+eps)
          er3d(icol,ilev)    =3.0_r8*(vaeros_tot(icol,ilev)+vaerol_tot(icol,ilev)) &
               /(aaeros_tot(icol,ilev)+aaerol_tot(icol,ilev)+eps)

          ! column integrated dry aerosol surface areas and volumes
          ! for r<0.5um and r>0.5um (s and l, respectively).
          aaercols(icol)=aaercols(icol)+aaeros_tot(icol,ilev)
          aaercoll(icol)=aaercoll(icol)+aaerol_tot(icol,ilev)
          vaercols(icol)=vaercols(icol)+vaeros_tot(icol,ilev)
          vaercoll(icol)=vaercoll(icol)+vaerol_tot(icol,ilev)

          ! then add background and externally mixed BC, OC and SO4 to mass concentrations
          c_bc_ac(icol,ilev)= c_bc(icol,ilev)
          c_bc_0(icol,ilev) = Nnatk(icol,ilev,0)*cintbg(icol,ilev,0)
          c_bc_2(icol,ilev) = Nnatk(icol,ilev,2)*cintbg(icol,ilev,2)
          c_bc_4(icol,ilev) = Nnatk(icol,ilev,4)*cintbg(icol,ilev,4)*faitbc(icol,ilev)
          c_bc_12(icol,ilev)= Nnatk(icol,ilev,12)*cknorm(icol,ilev,12)
          c_bc_14(icol,ilev)= Nnatk(icol,ilev,14)*cknorm(icol,ilev,14)*fnbc(icol,ilev)
          c_bc(icol,ilev)   = c_bc(icol,ilev)                                      &
               +Nnatk(icol,ilev,2)*cintbg(icol,ilev,2)                   &
               +Nnatk(icol,ilev,4)*cintbg(icol,ilev,4)*faitbc(icol,ilev) &
               +Nnatk(icol,ilev,0)*cintbg(icol,ilev,0)                   &
               +Nnatk(icol,ilev,12)*cknorm(icol,ilev,12)                 &
               +Nnatk(icol,ilev,14)*cknorm(icol,ilev,14)*fnbc(icol,ilev)
          c_bc05(icol,ilev)  = c_bc05(icol,ilev)                                   &
               +Nnatk(icol,ilev,2)*cintbg05(icol,ilev,2)                 &
               +Nnatk(icol,ilev,4)*cintbg05(icol,ilev,4)*faitbc(icol,ilev) &
               +Nnatk(icol,ilev,0)*cintbg05(icol,ilev,0)                 &
               +Nnatk(icol,ilev,12)*cknlt05(icol,ilev,12)                &
               +Nnatk(icol,ilev,14)*cknlt05(icol,ilev,14)*fnbc(icol,ilev)
          c_bc125(icol,ilev) = c_bc125(icol,ilev)                                  &
               +Nnatk(icol,ilev,2)*cintbg125(icol,ilev,2)                &
               +Nnatk(icol,ilev,4)*cintbg125(icol,ilev,4)*faitbc(icol,ilev) &
               +Nnatk(icol,ilev,0)*cintbg125(icol,ilev,0)                &
               +Nnatk(icol,ilev,12)*ckngt125(icol,ilev,12)               &
               +Nnatk(icol,ilev,14)*ckngt125(icol,ilev,14)*fnbc(icol,ilev)
          c_oc_ac(icol,ilev)= c_oc(icol,ilev)
          c_oc_4(icol,ilev)  = Nnatk(icol,ilev,4)*cintbg(icol,ilev,4)*(1.0_r8-faitbc(icol,ilev))
          c_oc_14(icol,ilev) = Nnatk(icol,ilev,14)*cknorm(icol,ilev,14)*(1.0_r8-fnbc(icol,ilev))
          c_oc(icol,ilev)    = c_oc(icol,ilev)                                           &
               +Nnatk(icol,ilev,1)*cintbg(icol,ilev,1)*f_soana(icol,ilev)         &
               +Nnatk(icol,ilev,4)*cintbg(icol,ilev,4)*(1.0_r8-faitbc(icol,ilev))   &
               +Nnatk(icol,ilev,14)*cknorm(icol,ilev,14)*(1.0_r8-fnbc(icol,ilev))
          c_oc05(icol,ilev)  = c_oc05(icol,ilev)                                         &
               +Nnatk(icol,ilev,1)*cintbg05(icol,ilev,1)*f_soana(icol,ilev)       &
               +Nnatk(icol,ilev,4)*cintbg05(icol,ilev,4)*(1.0_r8-faitbc(icol,ilev))  &
               +Nnatk(icol,ilev,14)*cknlt05(icol,ilev,14)*(1.0_r8-fnbc(icol,ilev))
          c_oc125(icol,ilev) = c_oc125(icol,ilev)                                        &
               +Nnatk(icol,ilev,1)*cintbg125(icol,ilev,1)*f_soana(icol,ilev)      &
               +Nnatk(icol,ilev,4)*cintbg125(icol,ilev,4)*(1.0_r8-faitbc(icol,ilev)) &
               +Nnatk(icol,ilev,14)*ckngt125(icol,ilev,14)*(1.0_r8-fnbc(icol,ilev))
          c_s4(icol,ilev)    = c_sa(icol,ilev)+c_sc(icol,ilev)          &
               +Nnatk(icol,ilev,1)*cintbg(icol,ilev,1)*(1.0_r8-f_soana(icol,ilev))   &
               +Nnatk(icol,ilev,5)*cintbg(icol,ilev,5)
          c_s405(icol,ilev)  = c_sa05(icol,ilev)+c_sc05(icol,ilev)      &
               +Nnatk(icol,ilev,1)*cintbg05(icol,ilev,1)*(1.0_r8-f_soana(icol,ilev)) &
               +Nnatk(icol,ilev,5)*cintbg05(icol,ilev,5)
          c_s4125(icol,ilev) = c_sa125(icol,ilev)+c_sc125(icol,ilev)    &
               +Nnatk(icol,ilev,1)*cintbg125(icol,ilev,1)*(1.0_r8-f_soana(icol,ilev)) &
               +Nnatk(icol,ilev,5)*cintbg125(icol,ilev,5)

          c_tot(icol,ilev)    = c_s4(icol,ilev) + c_oc(icol,ilev) + c_bc(icol,ilev) &
               + c_mi(icol,ilev) + c_ss(icol,ilev)
          c_tot125(icol,ilev) = c_s4125(icol,ilev) + c_oc125(icol,ilev) + c_bc125(icol,ilev) &
               + c_mi125(icol,ilev) + c_ss125(icol,ilev)
          c_tot05(icol,ilev) = c_s405(icol,ilev) + c_oc05(icol,ilev) + c_bc05(icol,ilev) &
               + c_mi05(icol,ilev) + c_ss05(icol,ilev)
          c_pm25(icol,ilev)   = c_tot(icol,ilev) - c_tot125(icol,ilev)
          c_pm1(icol,ilev)    = c_tot05(icol,ilev)

          ! mass mixing ratio:
          mmr_pm25(icol,ilev) = 1.e-9*c_pm25(icol,ilev)/rhoda(icol,ilev)
          mmr_pm1(icol,ilev)  = 1.e-9*c_pm1(icol,ilev)/rhoda(icol,ilev)

          c_s4_a(icol,ilev) = c_sa(icol,ilev)+c_sc(icol,ilev)
          c_s4_1(icol,ilev) = Nnatk(icol,ilev,1)*cintbg(icol,ilev,1)*(1.0_r8-f_soana(icol,ilev))
          c_s4_5(icol,ilev) = Nnatk(icol,ilev,5)*cintbg05(icol,ilev,5)

       end do ! icol
    enddo     ! ilev

    ! Total PM and PM2.5 (dry r>1.25um), surface values (ug/m3)
    do icol=1,ncol
       c_tots(icol) = c_tot(icol,pver)
       c_tot125s(icol) = c_tot125(icol,pver)
       c_pm25s(icol) = c_pm25(icol,pver)

    enddo

    ! Effective, column integrated, radii for particles
    ! smaller and greater than 0.5um, and for all radii
    do icol=1,ncol
       derlt05(icol)=3.0_r8*vaercols(icol)/(aaercols(icol)+eps)
       dergt05(icol)=3.0_r8*vaercoll(icol)/(aaercoll(icol)+eps)
       der(icol)=3.0_r8*(vaercols(icol)+vaercoll(icol)) /(aaercols(icol)+aaercoll(icol)+eps)
    enddo

    do icol=1,ncol
       dload_s4(icol)=0.0_r8
       dload_s4_a(icol)=0.0_r8
       dload_s4_1(icol)=0.0_r8
       dload_s4_5(icol)=0.0_r8
       dload_oc(icol)=0.0_r8
       dload_bc(icol)=0.0_r8
       dload_bc_ac(icol)=0.0_r8
       dload_bc_0(icol)=0.0_r8
       dload_bc_2(icol)=0.0_r8
       dload_bc_4(icol)=0.0_r8
       dload_bc_12(icol)=0.0_r8
       dload_bc_14(icol)=0.0_r8
       dload_oc_ac(icol)=0.0_r8
       dload_oc_4(icol)=0.0_r8
       dload_oc_14(icol)=0.0_r8
       do ilev=1,pver
          ! Layer thickness, unit km
          ! deltah=1.e-4_r8*(pint(icol,ilev+1)-pint(icol,ilev))/(rhoda(icol,ilev)*9.8_r8)
          deltah=deltah_km(icol,ilev)

          ! Modal and total mass concentrations for clean and dry aerosol,
          ! i.e. not including coag./cond./Aq. BC,OC,SO4 or condensed water.
          ! Units: ug/m3 for concentrations and mg/m2 (--> kg/m2 later) for mass loading.
          do imode=0,nmodes
             ck(icol,ilev,imode)=cknorm(icol,ilev,imode)*Nnatk(icol,ilev,imode)
             dload3d(icol,ilev,imode)=ck(icol,ilev,imode)*deltah
             dload(icol,imode)=dload(icol,imode)+dload3d(icol,ilev,imode)
          enddo
          nnat_0(icol,ilev) =Nnatk(icol,ilev,0)
          nnat_1(icol,ilev) =Nnatk(icol,ilev,1)
          nnat_2(icol,ilev) =Nnatk(icol,ilev,2)
          nnat_4(icol,ilev) =Nnatk(icol,ilev,4)
          nnat_5(icol,ilev) =Nnatk(icol,ilev,5)
          nnat_6(icol,ilev) =Nnatk(icol,ilev,6)
          nnat_7(icol,ilev) =Nnatk(icol,ilev,7)
          nnat_8(icol,ilev) =Nnatk(icol,ilev,8)
          nnat_9(icol,ilev) =Nnatk(icol,ilev,9)
          nnat_10(icol,ilev)=Nnatk(icol,ilev,10)
          nnat_12(icol,ilev)=Nnatk(icol,ilev,12)
          nnat_14(icol,ilev)=Nnatk(icol,ilev,14)

          ! mineral and sea-salt mass concentrations
          cmin(icol,ilev)=ck(icol,ilev,6)+ck(icol,ilev,7)
          cseas(icol,ilev)=ck(icol,ilev,8)+ck(icol,ilev,9)+ck(icol,ilev,10)

          ! just for checking purposes:
          dload_s4(icol)   =dload_s4(icol)   +c_s4(icol,ilev)*deltah
          dload_s4_a(icol) =dload_s4_a(icol) +c_s4_a(icol,ilev)*deltah
          dload_s4_1(icol) =dload_s4_1(icol) +c_s4_1(icol,ilev)*deltah
          dload_s4_5(icol) =dload_s4_5(icol) +c_s4_5(icol,ilev)*deltah
          dload_oc(icol)   =dload_oc(icol)   +c_oc(icol,ilev)*deltah
          dload_bc(icol)   =dload_bc(icol)   +c_bc(icol,ilev)*deltah
          !
          dload_bc_ac(icol)=dload_bc_ac(icol)+c_bc_ac(icol,ilev)*deltah
          dload_bc_0(icol)=dload_bc_0(icol)+c_bc_0(icol,ilev)*deltah
          dload_bc_2(icol)=dload_bc_2(icol)+c_bc_2(icol,ilev)*deltah
          dload_bc_4(icol)=dload_bc_4(icol)+c_bc_4(icol,ilev)*deltah
          dload_bc_12(icol)=dload_bc_12(icol)+c_bc_12(icol,ilev)*deltah
          dload_bc_14(icol)=dload_bc_14(icol)+c_bc_14(icol,ilev)*deltah
          dload_oc_ac(icol)=dload_oc_ac(icol)+c_oc_ac(icol,ilev)*deltah
          dload_oc_4(icol)=dload_oc_4(icol)+c_oc_4(icol,ilev)*deltah
          dload_oc_14(icol)=dload_oc_14(icol)+c_oc_14(icol,ilev)*deltah
          !
       end do  ! ilev
       dload_mi(icol)=dload(icol,6)+dload(icol,7)
       dload_ss(icol)=dload(icol,8)+dload(icol,9)+dload(icol,10)
    end do   ! icol

    call outfld('PMTOT  ',c_tots   ,pcols,lchnk)
    call outfld('PM25    ',c_pm25s ,pcols,lchnk)
    call outfld('PM2P5   ',c_pm25  ,pcols,lchnk)
    call outfld('MMRPM2P5',mmr_pm25,pcols,lchnk)
    call outfld('MMRPM1  ',mmr_pm1 ,pcols,lchnk)
    call outfld('MMRPM2P5_SRF',mmr_pm25(:pcols,pver),pcols,lchnk)

    call outfld('DLOAD_MI',dload_mi,pcols,lchnk)
    call outfld('DLOAD_SS',dload_ss,pcols,lchnk)
    call outfld('DLOAD_S4',dload_s4,pcols,lchnk)
    call outfld('DLOAD_OC',dload_oc,pcols,lchnk)
    call outfld('DLOAD_BC',dload_bc,pcols,lchnk)

    call outfld('LOADBCAC',dload_bc_ac,pcols,lchnk)
    call outfld('LOADBC0 ',dload_bc_0,pcols,lchnk)
    call outfld('LOADBC2 ',dload_bc_2,pcols,lchnk)
    call outfld('LOADBC4 ',dload_bc_4,pcols,lchnk)
    call outfld('LOADBC12',dload_bc_12,pcols,lchnk)
    call outfld('LOADBC14',dload_bc_14,pcols,lchnk)
    call outfld('LOADOCAC',dload_oc_ac,pcols,lchnk)
    call outfld('LOADOC4 ',dload_oc_4,pcols,lchnk)
    call outfld('LOADOC14',dload_oc_14,pcols,lchnk)

    ! condensed water loading (mg/m2)

    ! number concentrations (1/cm3)
    call outfld('NNAT_0  ',nnat_0 ,pcols,lchnk)
    call outfld('NNAT_1  ',nnat_1 ,pcols,lchnk)
    call outfld('NNAT_2  ',nnat_2 ,pcols,lchnk)
    call outfld('NNAT_4  ',nnat_4 ,pcols,lchnk)
    call outfld('NNAT_5  ',nnat_5 ,pcols,lchnk)
    call outfld('NNAT_6  ',nnat_6 ,pcols,lchnk)
    call outfld('NNAT_7  ',nnat_7 ,pcols,lchnk)
    call outfld('NNAT_8  ',nnat_8 ,pcols,lchnk)
    call outfld('NNAT_9  ',nnat_9 ,pcols,lchnk)
    call outfld('NNAT_10 ',nnat_10,pcols,lchnk)
    call outfld('NNAT_12 ',nnat_12,pcols,lchnk)
    call outfld('NNAT_14 ',nnat_14,pcols,lchnk)

    ! column integrated effective dry radii (um)
    call outfld('DERLT05 ',derlt05,pcols,lchnk)
    call outfld('DERGT05 ',dergt05,pcols,lchnk)
    call outfld('DER     ',der    ,pcols,lchnk)

    ! Extra AeroCom diagnostics requiring table look-ups with RH = constant
    ! Note: using xrhnull etc as proxy for constant RH input values
    irf = 1
    do ilev=1,pver
       do icol=1,ncol
          xrhnull(icol,ilev) = xrhrf(irf)
          irh1null(icol,ilev) = irhrf1(irf)
       end do
    enddo
    call opticsAtConstRh(&
         lchnk, ncol, pint, rhoda, Nnatk, xrhnull, irh1null, irf, &
         xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
         xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
         xfombg, ifombg1, vnbcarr, vaitbcarr, v_soana,      &
         bext440, bext500, bext550, bext670, bext870,       &
         bebg440, bebg500, bebg550, bebg670, bebg870,       &
         bebc440, bebc500, bebc550, bebc670, bebc870,       &
         beoc440, beoc500, beoc550, beoc670, beoc870,       &
         besu440, besu500, besu550, besu670, besu870,       &
         babs440, babs500, babs550, babs670, babs870,       &
         bebglt1, bebggt1, bebclt1, bebcgt1,                &
         beoclt1, beocgt1, bes4lt1, bes4gt1,                &
         backsc550, babg550, babc550, baoc550, basu550,     &
         bext440n, bext500n, bext550n, bext670n, bext870n,  &
         bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,  &
         bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,  &
         beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,  &
         besu440n, besu500n, besu550n, besu670n, besu870n,  &
         babs440n, babs500n, babs550n, babs670n, babs870n,  &
         bebglt1n, bebggt1n, bebclt1n, bebcgt1n,            &
         beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,            &
         backsc550n, babg550n, babc550n, baoc550n, basu550n)

  end subroutine aerocom2

  ! ==========================================================
  subroutine opticsAtConstRh (&
       lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irf,   &
       xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
       xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
       xfombg, ifombg1, vnbc, vaitbc, v_soana,            &
       bext440, bext500, bext550, bext670, bext870,       &
       bebg440, bebg500, bebg550, bebg670, bebg870,       &
       bebc440, bebc500, bebc550, bebc670, bebc870,       &
       beoc440, beoc500, beoc550, beoc670, beoc870,       &
       besu440, besu500, besu550, besu670, besu870,       &
       babs440, babs500, babs550, babs670, babs870,       &
       bebglt1, bebggt1, bebclt1, bebcgt1,                &
       beoclt1, beocgt1, bes4lt1, bes4gt1,                &
       backsc550, babg550, babc550, baoc550, basu550,     &
       bext440n, bext500n, bext550n, bext670n, bext870n,  &
       bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,  &
       bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,  &
       beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,  &
       besu440n, besu500n, besu550n, besu670n, besu870n,  &
       babs440n, babs500n, babs550n, babs670n, babs870n,  &
       bebglt1n, bebggt1n, bebclt1n, bebcgt1n,            &
       beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,            &
       backsc550n, babg550n, babc550n, baoc550n, basu550n)

    ! Extra AeroCom diagnostics requiring table look-ups with constant/fixed RH,
    ! imode.e. for RH = (/"00","40","55","65","75","85" /) (see opttab.F90)

    ! Arguments
    integer,  intent(in) :: lchnk                      ! chunk identifier
    integer,  intent(in) :: ncol                       ! number of atmospheric columns
    real(r8), intent(in) :: pint(pcols,pverp)          ! Model interface pressures (10*Pa)
    real(r8), intent(in) :: rhoda(pcols,pver)          ! Density of dry air (kg/m^3)
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    integer , intent(in) :: irf
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration
    real(r8), intent(in) :: vnbc(pcols,pver)
    real(r8), intent(in) :: vaitbc(pcols,pver)
    real(r8), intent(in) :: v_soana(pcols,pver)
    real(r8), intent(in) :: xfombg(pcols,pver)
    integer,  intent(in) :: ifombg1(pcols,pver)
    real(r8), intent(in) :: xfbcbg(pcols,pver)
    integer,  intent(in) :: ifbcbg1(pcols,pver)
    real(r8), intent(in) :: xfbcbgn(pcols,pver)
    integer,  intent(in) :: ifbcbgn1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! facm for use in the interpolations
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)   ! fbcm for use in the interpolations
    integer,  intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! faqm for use in the interpolations
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)
    real(r8), intent(out) :: bext440(pcols,pver,0:nbmodes), babs440(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext500(pcols,pver,0:nbmodes), babs500(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext550(pcols,pver,0:nbmodes), babs550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext670(pcols,pver,0:nbmodes), babs670(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext870(pcols,pver,0:nbmodes), babs870(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg440(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg500(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg550(pcols,pver,0:nbmodes), babg550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg670(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg870(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc440(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc500(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc550(pcols,pver,0:nbmodes), babc550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc670(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc870(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc440(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc500(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc550(pcols,pver,0:nbmodes), baoc550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc670(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc870(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu440(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu500(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu550(pcols,pver,0:nbmodes), basu550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu670(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu870(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebglt1(pcols,pver,0:nbmodes), bebggt1(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebclt1(pcols,pver,0:nbmodes), bebcgt1(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoclt1(pcols,pver,0:nbmodes), beocgt1(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bes4lt1(pcols,pver,0:nbmodes), bes4gt1(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: backsc550(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext440n(pcols,pver,0:nbmodes), babs440n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext500n(pcols,pver,0:nbmodes), babs500n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext550n(pcols,pver,0:nbmodes), babs550n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext670n(pcols,pver,0:nbmodes), babs670n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bext870n(pcols,pver,0:nbmodes), babs870n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg440n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg500n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg550n(pcols,pver,0:nbmodes), babg550n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg670n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebg870n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc440n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc500n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc550n(pcols,pver,0:nbmodes), babc550n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc670n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebc870n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc440n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc500n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc550n(pcols,pver,0:nbmodes), baoc550n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc670n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoc870n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu440n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu500n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu550n(pcols,pver,0:nbmodes), basu550n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu670n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: besu870n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebglt1n(pcols,pver,0:nbmodes), bebggt1n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bebclt1n(pcols,pver,0:nbmodes), bebcgt1n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: beoclt1n(pcols,pver,0:nbmodes), beocgt1n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: bes4lt1n(pcols,pver,0:nbmodes), bes4gt1n(pcols,pver,0:nbmodes)
    real(r8), intent(out) :: backsc550n(pcols,pver,0:nbmodes)

    ! Local variables
    integer  :: imode, ilev, icol, kcomp, irh, irelh
    integer  :: iloop
    real(r8) :: deltah
    real(r8) :: dod550rh(pcols), abs550rh(pcols)
    real(r8) :: babg440(pcols,pver,0:nbmodes)
    real(r8) :: babg500(pcols,pver,0:nbmodes)
    real(r8) :: babg670(pcols,pver,0:nbmodes)
    real(r8) :: babg870(pcols,pver,0:nbmodes)
    real(r8) :: babc440(pcols,pver,0:nbmodes)
    real(r8) :: babc500(pcols,pver,0:nbmodes)
    real(r8) :: babc670(pcols,pver,0:nbmodes)
    real(r8) :: babc870(pcols,pver,0:nbmodes)
    real(r8) :: baoc440(pcols,pver,0:nbmodes)
    real(r8) :: baoc500(pcols,pver,0:nbmodes)
    real(r8) :: baoc670(pcols,pver,0:nbmodes)
    real(r8) :: baoc870(pcols,pver,0:nbmodes)
    real(r8) :: basu440(pcols,pver,0:nbmodes)
    real(r8) :: basu500(pcols,pver,0:nbmodes)
    real(r8) :: basu670(pcols,pver,0:nbmodes)
    real(r8) :: basu870(pcols,pver,0:nbmodes)
    real(r8) :: ec550rh_aer(pcols,pver)
    real(r8) :: abs550rh_aer(pcols,pver)
    real(r8) :: bebglt1t(pcols,pver)
    real(r8) :: bebclt1t(pcols,pver)
    real(r8) :: beoclt1t(pcols,pver)
    real(r8) :: bes4lt1t(pcols,pver)
    real(r8) :: basu550tot(pcols,pver)
    real(r8) :: babc550tot(pcols,pver)
    real(r8) :: baoc550tot(pcols,pver)
    real(r8) :: babc550xt(pcols,pver)
    real(r8) :: baoc550xt(pcols,pver)
    real(r8) :: ba550x(pcols,pver,nbmp1:nmodes)
    real(r8) :: belt1x(pcols,pver,nbmp1:nmodes)

    ! Additional AeroCom Phase III output:
    real(r8) :: ec440rh_aer(pcols,pver), abs440rh_aer(pcols,pver)
    real(r8) :: ec870rh_aer(pcols,pver), abs870rh_aer(pcols,pver)
    real(r8) :: be550lt1_aer(pcols,pver,0:nbmodes), ec550rhlt1_aer(pcols,pver)
    real(r8) :: abs550rh_bc(pcols,pver), abs550rh_oc(pcols,pver)
    real(r8) :: abs550rh_su(pcols,pver), abs550rh_ss(pcols,pver)
    real(r8) :: abs550rh_du(pcols,pver), ec550rhlt1_bc(pcols,pver)
    real(r8) :: ec550rhlt1_oc(pcols,pver), ec550rhlt1_su(pcols,pver)
    real(r8) :: ec550rhlt1_ss(pcols,pver), ec550rhlt1_du(pcols,pver)
    !
    real(r8) :: babg440n(pcols,pver,0:nbmodes)
    real(r8) :: babg500n(pcols,pver,0:nbmodes)
    real(r8) :: babg670n(pcols,pver,0:nbmodes)
    real(r8) :: babg870n(pcols,pver,0:nbmodes)
    real(r8) :: babc440n(pcols,pver,0:nbmodes)
    real(r8) :: babc500n(pcols,pver,0:nbmodes)
    real(r8) :: babc670n(pcols,pver,0:nbmodes)
    real(r8) :: babc870n(pcols,pver,0:nbmodes)
    real(r8) :: baoc440n(pcols,pver,0:nbmodes)
    real(r8) :: baoc500n(pcols,pver,0:nbmodes)
    real(r8) :: baoc670n(pcols,pver,0:nbmodes)
    real(r8) :: baoc870n(pcols,pver,0:nbmodes)
    real(r8) :: basu440n(pcols,pver,0:nbmodes)
    real(r8) :: basu500n(pcols,pver,0:nbmodes)
    real(r8) :: basu670n(pcols,pver,0:nbmodes)
    real(r8) :: basu870n(pcols,pver,0:nbmodes)
    real(r8) :: bedustlt1(pcols,pver), bedustgt1(pcols,pver)
    real(r8) :: besslt1(pcols,pver), bessgt1(pcols,pver)
    real(r8) :: bbclt1xt(pcols,pver)
    real(r8) :: boclt1xt(pcols,pver)
    real(r8) :: xrh_sc(pcols,pver)
    integer  :: irh1_sc(pcols,pver)
    !-------------------------------------------------------------------------

    belt1x(:,:,:) = 0._r8

    ! scale RH by rh_fine_aer_scale_fact_optics for all modes except dust (6-7) and seasalt (8-10).
    do ilev=1,pver
       do icol=1,ncol
          xrh_sc(icol,ilev) = min(max(xrh(icol,ilev)*rh_fine_aer_scale_fact_optics, rh(1)), rh(10))
       end do
    end do

    irh1_sc(:,:) = 1
    do irelh=1,9
       do ilev=1,pver
          do icol=1,ncol
             if (xrh_sc(icol,ilev) >= rh(irelh) .and. xrh_sc(icol,ilev) <= rh(irelh+1)) then
                irh1_sc(icol,ilev) = irelh
             end if
          end do
       end do
    end do

    ! BC(ax) mode (hydrophobic, so no rhum needed here):
    call intaeropt0(lchnk, ncol, Nnatk,               &
         bext440, bext500, bext550, bext670, bext870,   &
         bebg440, bebg500, bebg550, bebg670, bebg870,   &
         bebc440, bebc500, bebc550, bebc670, bebc870,   &
         beoc440, beoc500, beoc550, beoc670, beoc870,   &
         besu440, besu500, besu550, besu670, besu870,   &
         babs440, babs500, babs550, babs670, babs870,   &
         bebglt1, bebggt1, bebclt1, bebcgt1,            &
         beoclt1, beocgt1, bes4lt1, bes4gt1,            &
         backsc550, babg550, babc550, baoc550, basu550)

    ! SO4(Ait), BC(Ait) and OC(Ait) modes:
    call intaeropt1(lchnk, ncol, xrh_sc, irh1_sc, 1,  &
         Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1,&
         bext440, bext500, bext550, bext670, bext870,   &
         bebg440, bebg500, bebg550, bebg670, bebg870,   &
         bebc440, bebc500, bebc550, bebc670, bebc870,   &
         beoc440, beoc500, beoc550, beoc670, beoc870,   &
         besu440, besu500, besu550, besu670, besu870,   &
         babs440, babs500, babs550, babs670, babs870,   &
         bebglt1, bebggt1, bebclt1, bebcgt1,            &
         beoclt1, beocgt1, bes4lt1, bes4gt1,            &
         backsc550, babg550, babc550, baoc550, basu550)

    call intaeropt2to3(lchnk, ncol, xrh_sc, irh1_sc, 2, &
         Nnatk, xct, ict1, xfac, ifac1,           &
         bext440, bext500, bext550, bext670, bext870,     &
         bebg440, bebg500, bebg550, bebg670, bebg870,     &
         bebc440, bebc500, bebc550, bebc670, bebc870,     &
         beoc440, beoc500, beoc550, beoc670, beoc870,     &
         besu440, besu500, besu550, besu670, besu870,     &
         babs440, babs500, babs550, babs670, babs870,     &
         bebglt1, bebggt1, bebclt1, bebcgt1,              &
         beoclt1, beocgt1, bes4lt1, bes4gt1,              &
         backsc550, babg550, babc550, baoc550, basu550)

    ! BC&OC(Ait) (4), OC&BC(Ait) mode
    call intaeropt4(lchnk, ncol, xrh_sc, irh1_sc, 4, Nnatk,  &
         xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, xfaq, ifaq1, &
         bext440, bext500, bext550, bext670, bext870,          &
         bebg440, bebg500, bebg550, bebg670, bebg870,          &
         bebc440, bebc500, bebc550, bebc670, bebc870,          &
         beoc440, beoc500, beoc550, beoc670, beoc870,          &
         besu440, besu500, besu550, besu670, besu870,          &
         babs440, babs500, babs550, babs670, babs870,          &
         bebglt1, bebggt1, bebclt1, bebcgt1,                   &
         beoclt1, beocgt1, bes4lt1, bes4gt1,                   &
         backsc550, babg550, babc550, baoc550, basu550)

    ! SO4(Ait75) (5) — scaled RH, matching interpol5to10 in oslo_aero_optical_params.F90:
    call intaeropt5to10(lchnk, ncol, xrh_sc, irh1_sc, 5, Nnatk,   &
         xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
         bext440, bext500, bext550, bext670, bext870,      &
         bebg440, bebg500, bebg550, bebg670, bebg870,      &
         bebc440, bebc500, bebc550, bebc670, bebc870,      &
         beoc440, beoc500, beoc550, beoc670, beoc870,      &
         besu440, besu500, besu550, besu670, besu870,      &
         babs440, babs500, babs550, babs670, babs870,      &
         bebglt1, bebggt1, bebclt1, bebcgt1,               &
         beoclt1, beocgt1, bes4lt1, bes4gt1,               &
         backsc550, babg550, babc550, baoc550, basu550)
    ! Mineral (6-7) and Sea-salt (8-10) modes — unscaled RH:
    do kcomp=6,10
       call intaeropt5to10(lchnk, ncol, xrh, irh1, kcomp, Nnatk,   &
            xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
            bext440, bext500, bext550, bext670, bext870,      &
            bebg440, bebg500, bebg550, bebg670, bebg870,      &
            bebc440, bebc500, bebc550, bebc670, bebc870,      &
            beoc440, beoc500, beoc550, beoc670, beoc870,      &
            besu440, besu500, besu550, besu670, besu870,      &
            babs440, babs500, babs550, babs670, babs870,      &
            bebglt1, bebggt1, bebclt1, bebcgt1,               &
            beoclt1, beocgt1, bes4lt1, bes4gt1,               &
            backsc550, babg550, babc550, baoc550, basu550)
    end do

    ! then to the externally mixed SO4(n), BC(n) and OC(n) modes:
    call intaeropt2to3(lchnk, ncol, xrh_sc, irh1_sc, 12,  &
         Nnatk, xct, ict1, xfac, ifac1,                    &
         bext440n, bext500n, bext550n, bext670n, bext870n, &
         bebg440n, bebg500n, bebg550n, bebg670n, bebg870n, &
         bebc440n, bebc500n, bebc550n, bebc670n, bebc870n, &
         beoc440n, beoc500n, beoc550n, beoc670n, beoc870n, &
         besu440n, besu500n, besu550n, besu670n, besu870n, &
         babs440n, babs500n, babs550n, babs670n, babs870n, &
         bebglt1n, bebggt1n, bebclt1n, bebcgt1n,           &
         beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,           &
         backsc550n, babg550n, babc550n, baoc550n, basu550n)

    ! finally the BC&OC(n) mode:
    call intaeropt4(lchnk, ncol, xrh_sc, irh1_sc, 14, Nnatk,    &
         xfbcbgn, ifbcbgn1, xct, ict1, xfac, ifac1, xfaq, ifaq1, &
         bext440n, bext500n, bext550n, bext670n, bext870n,       &
         bebg440n, bebg500n, bebg550n, bebg670n, bebg870n,       &
         bebc440n, bebc500n, bebc550n, bebc670n, bebc870n,       &
         beoc440n, beoc500n, beoc550n, beoc670n, beoc870n,       &
         besu440n, besu500n, besu550n, besu670n, besu870n,       &
         babs440n, babs500n, babs550n, babs670n, babs870n,       &
         bebglt1n, bebggt1n, bebclt1n, bebcgt1n,                 &
         beoclt1n, beocgt1n, bes4lt1n, bes4gt1n,                 &
         backsc550n, babg550n, babc550n, baoc550n, basu550n)

    ! Initialization
    do ilev=1,pver
       do icol=1,ncol
          ec550rh_aer(icol,ilev)    =0.0_r8
          abs550rh_aer(icol,ilev)   =0.0_r8
          ec550rhlt1_aer(icol,ilev) =0.0_r8
          abs550rh_bc(icol,ilev)    =0.0_r8
          abs550rh_oc(icol,ilev)    =0.0_r8
          abs550rh_su(icol,ilev)    =0.0_r8
          abs550rh_ss(icol,ilev)    =0.0_r8
          abs550rh_du(icol,ilev)    =0.0_r8
          ec440rh_aer(icol,ilev)    =0.0_r8
          abs440rh_aer(icol,ilev)   =0.0_r8
          ec870rh_aer(icol,ilev)    =0.0_r8
          abs870rh_aer(icol,ilev)   =0.0_r8
          basu550tot(icol,ilev)     =0.0_r8
          babc550tot(icol,ilev)     =0.0_r8
          baoc550tot(icol,ilev)     =0.0_r8
          bebglt1t(icol,ilev)       =0.0_r8
          bebclt1t(icol,ilev)       =0.0_r8
          beoclt1t(icol,ilev)       =0.0_r8
          bes4lt1t(icol,ilev)       =0.0_r8
          bedustlt1(icol,ilev)      =0.0_r8
          besslt1(icol,ilev)        =0.0_r8
       end do
    end do
    do icol=1,ncol
       dod550rh(icol)=0.0_r8
       abs550rh(icol)=0.0_r8
    end do

    ! Calculation of extinction at given RH and absorption for all r and for r<0.5um
    do ilev=1,pver
       do icol=1,ncol

          do imode = 0,10
             ec550rh_aer(icol,ilev)  = ec550rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext550(icol,ilev,imode)
             abs550rh_aer(icol,ilev) = abs550rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs550(icol,ilev,imode)
             ec440rh_aer(icol,ilev)  = ec440rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext440(icol,ilev,imode)
             abs440rh_aer(icol,ilev) = abs440rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs440(icol,ilev,imode)
             ec870rh_aer(icol,ilev)  = ec870rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext870(icol,ilev,imode)
             abs870rh_aer(icol,ilev) = abs870rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs870(icol,ilev,imode)
             basu550tot(icol,ilev)   = basu550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*basu550(icol,ilev,imode)
             babc550tot(icol,ilev)   = babc550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*babc550(icol,ilev,imode)
             baoc550tot(icol,ilev)   = baoc550tot(icol,ilev)   +Nnatk(icol,ilev,imode)*baoc550(icol,ilev,imode)
             bes4lt1t(icol,ilev)     = bes4lt1t(icol,ilev)     +Nnatk(icol,ilev,imode)*bes4lt1(icol,ilev,imode)
             bebclt1t(icol,ilev)     = bebclt1t(icol,ilev)     +Nnatk(icol,ilev,imode)*bebclt1(icol,ilev,imode)
             beoclt1t(icol,ilev)     = beoclt1t(icol,ilev)     +Nnatk(icol,ilev,imode)*beoclt1(icol,ilev,imode)
          enddo

          do imode = 11,14
             ec550rh_aer(icol,ilev)  = ec550rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext550n(icol,ilev,imode-10)
             abs550rh_aer(icol,ilev) = abs550rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs550n(icol,ilev,imode-10)
             ec440rh_aer(icol,ilev)  = ec440rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext440n(icol,ilev,imode-10)
             abs440rh_aer(icol,ilev) = abs440rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs440n(icol,ilev,imode-10)
             ec870rh_aer(icol,ilev)  = ec870rh_aer(icol,ilev)  +Nnatk(icol,ilev,imode)*bext870n(icol,ilev,imode-10)
             abs870rh_aer(icol,ilev) = abs870rh_aer(icol,ilev) +Nnatk(icol,ilev,imode)*babs870n(icol,ilev,imode-10)
             ba550x(icol,ilev,imode)     = babs550n(icol,ilev,imode-10)
             belt1x(icol,ilev,imode)     = bebglt1n(icol,ilev,imode-10)
          enddo

          do imode=6,7
             bedustlt1(icol,ilev) = bedustlt1(icol,ilev) + Nnatk(icol,ilev,imode)*bebglt1(icol,ilev,imode)
          enddo
          do imode=8,10
             besslt1(icol,ilev) = besslt1(icol,ilev) + Nnatk(icol,ilev,imode)*bebglt1(icol,ilev,imode)
          enddo
          ec550rhlt1_du(icol,ilev) = bedustlt1(icol,ilev)
          ec550rhlt1_ss(icol,ilev) = besslt1(icol,ilev)

          !soa: *(1-v_soan) for the sulfate volume fraction of mode 11
          bbclt1xt(icol,ilev) = Nnatk(icol,ilev,12)*belt1x(icol,ilev,12) &
                           + Nnatk(icol,ilev,14)*belt1x(icol,ilev,14)*vnbc(icol,ilev)

          !soa + v_soan part of mode 11 for the OC volume fraction of that mode
          boclt1xt(icol,ilev) = Nnatk(icol,ilev,13)*belt1x(icol,ilev,13) &
                           + Nnatk(icol,ilev,14)*belt1x(icol,ilev,14)*(1.0_r8-vnbc(icol,ilev))

          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          ec550rhlt1_su(icol,ilev) = bes4lt1t(icol,ilev)                          & ! condensate
               + Nnatk(icol,ilev,1)*bebglt1(icol,ilev,1)*(1.0_r8-v_soana(icol,ilev)) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*bebglt1(icol,ilev,5)                            ! background, SO4(Ait75) mode (5)
          ec550rhlt1_bc(icol,ilev) = bebclt1t(icol,ilev)+bbclt1xt(icol,ilev)         & ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*bebglt1(icol,ilev,2)                          & ! background, BC(Ait) mode (2)
               + Nnatk(icol,ilev,4)*bebglt1(icol,ilev,4)*vaitbc(icol,ilev)           & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*bebglt1(icol,ilev,0)                            ! background, BC(ax) mode (0)

          !soa + v_soan part of mode 11 for the OC volume fraction of that mode
          ec550rhlt1_oc(icol,ilev) = beoclt1t(icol,ilev)+boclt1xt(icol,ilev)         & ! coagulated + n-mode OC (13)
               + Nnatk(icol,ilev,3)*bebglt1(icol,ilev,3)                          & ! background, OC(Ait) mode (3)
               + Nnatk(icol,ilev,4)*bebglt1(icol,ilev,4)*(1.0_r8-vaitbc(icol,ilev))  & ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,1)*bebglt1(icol,ilev,1)*v_soana(icol,ilev)

          ec550rhlt1_aer(icol,ilev) = ec550rhlt1_su(icol,ilev)+ec550rhlt1_bc(icol,ilev) &
               + ec550rhlt1_oc(icol,ilev) + ec550rhlt1_ss(icol,ilev)+ec550rhlt1_du(icol,ilev)
          ec550rhlt1_aer(icol,ilev) = 1.e-3_r8*ec550rhlt1_aer(icol,ilev)

          abs550rh_du(icol,ilev) = Nnatk(icol,ilev,6)*babg550(icol,ilev,6) &
               + Nnatk(icol,ilev,7)*babg550(icol,ilev,7)
          abs550rh_ss(icol,ilev) = Nnatk(icol,ilev,8)*babg550(icol,ilev,8) &
               + Nnatk(icol,ilev,9)*babg550(icol,ilev,9) &
               + Nnatk(icol,ilev,10)*babg550(icol,ilev,10)

          ! soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          abs550rh_su(icol,ilev) = basu550tot(icol,ilev)                   &        ! condensate:w
               + (1.0_r8-v_soana(icol,ilev))*Nnatk(icol,ilev,1)*babg550(icol,ilev,1) & ! background, SO4(Ait) mode (1)
               + Nnatk(icol,ilev,5)*babg550(icol,ilev,5)                            ! background, SO4(Ait75) mode (5)

          ! soa: *(1-v_soan) for the sulfate volume fraction
          babc550xt(icol,ilev) = Nnatk(icol,ilev,12)*ba550x(icol,ilev,12)  &
               + Nnatk(icol,ilev,14)*ba550x(icol,ilev,14)*vnbc(icol,ilev)

          baoc550xt(icol,ilev) = Nnatk(icol,ilev,13)*ba550x(icol,ilev,13) &
               + Nnatk(icol,ilev,14)*ba550x(icol,ilev,14)*(1.0_r8-vnbc(icol,ilev))

          abs550rh_bc(icol,ilev) = babc550tot(icol,ilev)+babc550xt(icol,ilev) &        ! coagulated + n-mode BC (12)
               + Nnatk(icol,ilev,2)*babg550(icol,ilev,2) &                          ! background, BC(Ait) mode (2)
               + vaitbc(icol,ilev)*Nnatk(icol,ilev,4)*babg550(icol,ilev,4) &           ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,ilev,0)*babg550(icol,ilev,0)                            ! background, BC(ax) mode (0)

          abs550rh_oc(icol,ilev) = baoc550tot(icol,ilev)+baoc550xt(icol,ilev) &        ! coagulated + n-mode OC (13)
               + v_soana(icol,ilev)*Nnatk(icol,ilev,1)*babg550(icol,ilev,1) &          ! SOA fraction of mode 1
               + Nnatk(icol,ilev,3)*babg550(icol,ilev,3) &                          ! background, OC(Ait) mode (3)
               + (1.0_r8-vaitbc(icol,ilev))*Nnatk(icol,ilev,4)*babg550(icol,ilev,4)    ! background in OC&BC(Ait) mode (4)

          deltah=1.e-4_r8*(pint(icol,ilev+1)-pint(icol,ilev))/(rhoda(icol,ilev)*9.8_r8)
          dod550rh(icol) = dod550rh(icol)+ec550rh_aer(icol,ilev)*deltah
          abs550rh(icol) = abs550rh(icol)+abs550rh_aer(icol,ilev)*deltah

          ec550rh_aer(icol,ilev)  = 1.e-3_r8*ec550rh_aer(icol,ilev)
          abs550rh_aer(icol,ilev) = 1.e-3_r8*abs550rh_aer(icol,ilev)
          ec440rh_aer(icol,ilev)  = 1.e-3_r8*ec440rh_aer(icol,ilev)
          abs440rh_aer(icol,ilev) = 1.e-3_r8*abs440rh_aer(icol,ilev)
          ec870rh_aer(icol,ilev)  = 1.e-3_r8*ec870rh_aer(icol,ilev)
          abs870rh_aer(icol,ilev) = 1.e-3_r8*abs870rh_aer(icol,ilev)

          abs550rh_bc(icol,ilev)  = 1.e-3_r8*abs550rh_bc(icol,ilev)
          abs550rh_oc(icol,ilev)  = 1.e-3_r8*abs550rh_oc(icol,ilev)
          abs550rh_su(icol,ilev)  = 1.e-3_r8*abs550rh_su(icol,ilev)
          abs550rh_ss(icol,ilev)  = 1.e-3_r8*abs550rh_ss(icol,ilev)
          abs550rh_du(icol,ilev)  = 1.e-3_r8*abs550rh_du(icol,ilev)

       enddo
    enddo

    if (irf == 1) then
       call outfld('ECDRYAER',ec550rh_aer    ,pcols,lchnk)
       call outfld('ABSDRYAE',abs550rh_aer   ,pcols,lchnk)
       call outfld('OD550DRY',dod550rh       ,pcols,lchnk)       ! 2D variable
       call outfld('AB550DRY',abs550rh       ,pcols,lchnk)       ! 2D variable
       call outfld('ECDRY440',ec440rh_aer    ,pcols,lchnk)
       call outfld('ABSDR440',abs440rh_aer   ,pcols,lchnk)
       call outfld('ECDRY870',ec870rh_aer    ,pcols,lchnk)
       call outfld('ABSDR870',abs870rh_aer   ,pcols,lchnk)
       call outfld('ECDRYLT1',ec550rhlt1_aer ,pcols,lchnk)

       ! Since we do not have enough look-up table info to take abs550rhlt1_aer,
       ! instead take out abs550rh for each constituent:
       call outfld('ABSDRYBC',abs550rh_bc    ,pcols,lchnk)
       call outfld('ABSDRYOC',abs550rh_oc    ,pcols,lchnk)
       call outfld('ABSDRYSU',abs550rh_su    ,pcols,lchnk)
       call outfld('ABSDRYSS',abs550rh_ss    ,pcols,lchnk)
       call outfld('ABSDRYDU',abs550rh_du    ,pcols,lchnk)
    end if

  end subroutine opticsAtConstRh

end module oslo_aero_aerocom
