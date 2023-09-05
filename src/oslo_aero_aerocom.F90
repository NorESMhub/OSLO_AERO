module oslo_aero_aerocom

#ifdef AEROCOM
  
  use ppgrid
  use shr_kind_mod,    only: r8 => shr_kind_r8
  use cam_history,     only: outfld
  use physics_types,   only: physics_state
  !
  use oslo_aero_aerocom_opt, only: extinction_coeffs, extinction_coeffsn
  use oslo_aero_aerocom_dry, only: aerodry_prop
  use oslo_aero_sw_tables
  use oslo_aero_share
  use oslo_aero_params
  use oslo_aero_const

  public :: aerocom
  public :: opticsAtConstRh
  public :: intfrh

contains

  subroutine aerocom(daylight, Cam) 

    ! Arguments
    real(r8), intent(in)  :: Cam(pcols,pver,nbmodes) 

    ! Local variables
    integer  i, k, ib, icol, mplus10
    integer iloop
    logical  daylight(pcols)        ! SW calculations also at (polar) night in interpol* if daylight=.true.

    real(r8) Ctotdry(pcols,pver), Cwater(pcols,pver), mmr_aerh2o(pcols,pver), &
         dod550dry(pcols), abs550dry(pcols)

    real(r8) daerh2o(pcols),  dload(pcols,0:nmodes), dload3d(pcols,pver,0:nmodes), &
         dload_mi(pcols), dload_ss(pcols), &
         dload_s4(pcols), dload_oc(pcols), dload_bc(pcols), &
         dload_s4_a(pcols), dload_s4_1(pcols), dload_s4_5(pcols)

    real(r8) dload_bc_0(pcols), dload_bc_ac(pcols), dload_oc_ac(pcols), &
         dload_bc_2(pcols), dload_bc_4(pcols), dload_bc_12(pcols), dload_bc_14(pcols), &
         dload_oc_4(pcols), dload_oc_14(pcols)

    real(r8) cmin(pcols,pver), cseas(pcols,pver)

    real(r8) nnat_1(pcols,pver), nnat_2(pcols,pver), nnat_3(pcols,pver), &
         nnat_4(pcols,pver), nnat_5(pcols,pver), nnat_6(pcols,pver), &
         nnat_7(pcols,pver), nnat_8(pcols,pver), nnat_9(pcols,pver), &
         nnat_10(pcols,pver), nnat_12(pcols,pver), &
         nnat_14(pcols,pver), nnat_0(pcols,pver)

    real(r8) ck(pcols,pver,0:nmodes), cknorm(pcols,pver,0:nmodes), &
         cknlt05(pcols,pver,0:nmodes), ckngt125(pcols,pver,0:nmodes)

    real(r8) aaerosn(pcols,pver,nbmp1:nmodes), aaeroln(pcols,pver,nbmp1:nmodes), &
         vaerosn(pcols,pver,nbmp1:nmodes), vaeroln(pcols,pver,nbmp1:nmodes), &
         aaeros(pcols,pver,0:nbmodes), aaerol(pcols,pver,0:nbmodes), &
         vaeros(pcols,pver,0:nbmodes), vaerol(pcols,pver,0:nbmodes)

    real(r8) cintbg(pcols,pver,0:nbmodes), &
         cintbg05(pcols,pver,0:nbmodes), cintbg125(pcols,pver,0:nbmodes), &
         cintbc(pcols,pver,0:nbmodes), &
         cintbc05(pcols,pver,0:nbmodes), cintbc125(pcols,pver,0:nbmodes), &
         cintoc(pcols,pver,0:nbmodes), &
         cintoc05(pcols,pver,0:nbmodes), cintoc125(pcols,pver,0:nbmodes), &
         cintsc(pcols,pver,0:nbmodes), &
         cintsc05(pcols,pver,0:nbmodes), cintsc125(pcols,pver,0:nbmodes), &
         cintsa(pcols,pver,0:nbmodes), &
         cintsa05(pcols,pver,0:nbmodes), cintsa125(pcols,pver,0:nbmodes)

    real(r8) c_mi(pcols,pver), c_mi05(pcols,pver), c_mi125(pcols,pver), &
         c_ss(pcols,pver), c_ss05(pcols,pver), c_ss125(pcols,pver), &
         c_bc(pcols,pver), c_bc05(pcols,pver), c_bc125(pcols,pver), &
         c_oc(pcols,pver), c_oc05(pcols,pver), c_oc125(pcols,pver), &
         c_sa(pcols,pver), c_sa05(pcols,pver), c_sa125(pcols,pver), &
         c_sc(pcols,pver), c_sc05(pcols,pver), c_sc125(pcols,pver), &
         c_s4(pcols,pver), c_s405(pcols,pver), c_s4125(pcols,pver), &
         c_s4_a(pcols,pver), c_s4_1(pcols,pver), c_s4_5(pcols,pver)

    real(r8) c_bc_0(pcols,pver), c_bc_ac(pcols,pver), c_oc_ac(pcols,pver), &
         c_bc_2(pcols,pver), c_bc_4(pcols,pver), c_bc_12(pcols,pver), c_bc_14(pcols,pver), &
         c_oc_4(pcols,pver), c_oc_14(pcols,pver)

    real(r8) c_tots(pcols), c_tot125s(pcols), c_pm25s(pcols) ! = PM all sizes, PM>2.5um and PM<2.5um (PM2.5)

    real(r8) c_tot(pcols,pver), c_tot125(pcols,pver), c_pm25(pcols,pver), &
         mmr_pm25(pcols,pver), c_tot05(pcols,pver), c_pm1(pcols,pver), mmr_pm1(pcols,pver)

    real(r8) aaeros_tot(pcols,pver), aaerol_tot(pcols,pver), vaeros_tot(pcols,pver), &
         vaerol_tot(pcols,pver), aaercols(pcols), aaercoll(pcols), vaercols(pcols), &
         vaercoll(pcols), derlt05(pcols), dergt05(pcols), der(pcols), &
         erlt053d(pcols,pver), ergt053d(pcols,pver), er3d(pcols,pver)

    real(r8) bebglt1(pcols,pver,0:nbmodes), bebggt1(pcols,pver,0:nbmodes), &
         bebclt1(pcols,pver,0:nbmodes), bebcgt1(pcols,pver,0:nbmodes), &
         beoclt1(pcols,pver,0:nbmodes), beocgt1(pcols,pver,0:nbmodes), &
         bes4lt1(pcols,pver,0:nbmodes), bes4gt1(pcols,pver,0:nbmodes), &
         backsc550(pcols,pver,0:nbmodes), backsc550x(pcols,pver,nbmp1:nmodes), &
         backsc550tot(pcols,pver), ec550_aer(pcols,pver), abs550_aer(pcols,pver), &
         bs550_aer(pcols,pver)

    ! Additional AeroCom Phase III output:
    real(r8) asydry_aer(pcols,pver)    ! dry asymtot in the visible band
    !
    real(r8) ec550_so4(pcols,pver),ec550_bc(pcols,pver), ec550_pom(pcols,pver), &
         ec550_ss(pcols,pver), ec550_du(pcols,pver)

    real(r8) bebglt1n(pcols,pver,0:nbmodes), bebggt1n(pcols,pver,0:nbmodes), &
         bebclt1n(pcols,pver,0:nbmodes), bebcgt1n(pcols,pver,0:nbmodes), &
         beoclt1n(pcols,pver,0:nbmodes), beocgt1n(pcols,pver,0:nbmodes), &
         bes4lt1n(pcols,pver,0:nbmodes), bes4gt1n(pcols,pver,0:nbmodes), &
         backsc550n(pcols,pver,0:nbmodes)

    real(r8) bext440tot(pcols,pver), babs440tot(pcols,pver), &
         bext500tot(pcols,pver), babs500tot(pcols,pver), &
         bext550tot(pcols,pver), babs550tot(pcols,pver), &
         bext670tot(pcols,pver), babs670tot(pcols,pver), &
         bext870tot(pcols,pver), babs870tot(pcols,pver), &
         bebg440tot(pcols,pver), &
         bebg500tot(pcols,pver), &
         bebg550tot(pcols,pver), babg550tot(pcols,pver), &
         bebg670tot(pcols,pver), &
         bebg870tot(pcols,pver), &
         bebc440tot(pcols,pver), &
         bebc500tot(pcols,pver), &
         bebc550tot(pcols,pver), babc550tot(pcols,pver), &
         bebc670tot(pcols,pver), &
         bebc870tot(pcols,pver), &
         beoc440tot(pcols,pver), &
         beoc500tot(pcols,pver), &
         beoc550tot(pcols,pver), baoc550tot(pcols,pver), &
         beoc670tot(pcols,pver), &
         beoc870tot(pcols,pver), &
         besu440tot(pcols,pver), &
         besu500tot(pcols,pver), &
         besu550tot(pcols,pver), basu550tot(pcols,pver), &
         besu670tot(pcols,pver), &
         besu870tot(pcols,pver)

    real(r8) bebglt1t(pcols,pver), bebggt1t(pcols,pver), bebclt1t(pcols,pver), &
         bebcgt1t(pcols,pver), beoclt1t(pcols,pver), beocgt1t(pcols,pver), &
         bes4lt1t(pcols,pver), bes4gt1t(pcols,pver)

    real(r8) be440x(pcols,pver,nbmp1:nmodes), ba440x(pcols,pver,nbmp1:nmodes), &
         be500x(pcols,pver,nbmp1:nmodes), ba500x(pcols,pver,nbmp1:nmodes), &
         be550x(pcols,pver,nbmp1:nmodes), ba550x(pcols,pver,nbmp1:nmodes), &
         be670x(pcols,pver,nbmp1:nmodes), ba670x(pcols,pver,nbmp1:nmodes), &
         be870x(pcols,pver,nbmp1:nmodes), ba870x(pcols,pver,nbmp1:nmodes), &
         belt1x(pcols,pver,nbmp1:nmodes), begt1x(pcols,pver,nbmp1:nmodes)

    real(r8) bebc440xt(pcols,pver),babc440xt(pcols,pver), &
         bebc500xt(pcols,pver),babc500xt(pcols,pver), &
         bebc550xt(pcols,pver),babc550xt(pcols,pver), &
         bebc670xt(pcols,pver),babc670xt(pcols,pver), &
         bebc870xt(pcols,pver),babc870xt(pcols,pver), &
         beoc440xt(pcols,pver),baoc440xt(pcols,pver), &
         beoc500xt(pcols,pver),baoc500xt(pcols,pver), &
         beoc550xt(pcols,pver),baoc550xt(pcols,pver), &
         beoc670xt(pcols,pver),baoc670xt(pcols,pver), &
         beoc870xt(pcols,pver),baoc870xt(pcols,pver)

    real(r8) bbclt1xt(pcols,pver), &
         bbcgt1xt(pcols,pver), boclt1xt(pcols,pver), bocgt1xt(pcols,pver)

    real(r8) bint440du(pcols,pver), bint500du(pcols,pver), bint550du(pcols,pver), &
         bint670du(pcols,pver), bint870du(pcols,pver), &
         bint440ss(pcols,pver), bint500ss(pcols,pver), bint550ss(pcols,pver), &
         bint670ss(pcols,pver), bint870ss(pcols,pver), &
         baint550du(pcols,pver), baint550ss(pcols,pver)

    real(r8) bedustlt1(pcols,pver), bedustgt1(pcols,pver), &
         besslt1(pcols,pver), bessgt1(pcols,pver)

    real(r8) dod4403d(pcols,pver), abs4403d(pcols,pver), &
         dod4403d_ss(pcols,pver),   & ! abs4403d_ss(pcols,pver), &
         dod4403d_dust(pcols,pver), & ! abs4403d_dust(pcols,pver), &
         dod4403d_so4(pcols,pver),  & ! abs4403d_so4(pcols,pver), &
         dod4403d_bc(pcols,pver),   & ! abs4403d_bc(pcols,pver), &
         dod4403d_pom(pcols,pver),  & ! abs4403d_pom(pcols,pver), &
         dod5003d(pcols,pver), abs5003d(pcols,pver), &
         dod5003d_ss(pcols,pver),   & ! abs5003d_ss(pcols,pver), &
         dod5003d_dust(pcols,pver), & ! abs5003d_dust(pcols,pver), &
         dod5003d_so4(pcols,pver),  & ! abs5003d_so4(pcols,pver), &
         dod5003d_bc(pcols,pver),   & ! abs5003d_bc(pcols,pver), &
         dod5003d_pom(pcols,pver),  & ! abs5003d_pom(pcols,pver), &
         dod5503d(pcols,pver), abs5503d(pcols,pver), abs5503dalt(pcols,pver), &
         dod5503d_ss(pcols,pver), abs5503d_ss(pcols,pver), &
         dod5503d_dust(pcols,pver), abs5503d_dust(pcols,pver), &
         dod5503d_so4(pcols,pver), abs5503d_so4(pcols,pver), &
         dod5503d_bc(pcols,pver), abs5503d_bc(pcols,pver), &
         dod5503d_pom(pcols,pver), abs5503d_pom(pcols,pver), &
         dod6703d(pcols,pver), abs6703d(pcols,pver), &
         dod6703d_ss(pcols,pver),   & ! abs6703d_ss(pcols,pver), &
         dod6703d_dust(pcols,pver), & ! abs6703d_dust(pcols,pver), &
         dod6703d_so4(pcols,pver),  & ! abs6703d_so4(pcols,pver), &
         dod6703d_bc(pcols,pver),   & ! abs6703d_bc(pcols,pver), &
         dod6703d_pom(pcols,pver),  & ! abs6703d_pom(pcols,pver), &
         dod8703d(pcols,pver), abs8703d(pcols,pver), &
         dod8703d_ss(pcols,pver),   & ! abs8703d_ss(pcols,pver), &
         dod8703d_dust(pcols,pver), & ! abs8703d_dust(pcols,pver), &
         dod8703d_so4(pcols,pver),  & ! abs8703d_so4(pcols,pver), &
         dod8703d_bc(pcols,pver),   & ! abs8703d_bc(pcols,pver), &
         dod8703d_pom(pcols,pver) ! abs8703d_pom(pcols,pver)

    real(r8) dod5503dlt1_ss(pcols,pver), dod5503dgt1_ss(pcols,pver), &
         dod5503dlt1_dust(pcols,pver), dod5503dgt1_dust(pcols,pver), &
         dod5503dlt1_so4(pcols,pver), dod5503dgt1_so4(pcols,pver), &
         dod5503dlt1_bc(pcols,pver), dod5503dgt1_bc(pcols,pver), &
         dod5503dlt1_pom(pcols,pver), dod5503dgt1_pom(pcols,pver)

    real(r8) abs440(pcols), dod500(pcols), abs500(pcols),  &
         dod670(pcols),&
         abs670(pcols), abs870(pcols),                 &
         dod440_ss(pcols), dod440_dust(pcols), dod440_so4(pcols),     &
         dod440_bc(pcols), dod440_pom(pcols),                         &
         dod500_ss(pcols), dod500_dust(pcols), dod500_so4(pcols),     &
         dod500_bc(pcols), dod500_pom(pcols),                         &
         dod550_ss(pcols), dod550_dust(pcols), dod550_so4(pcols),     &
         dod550_bc(pcols), dod550_pom(pcols),                         &
         dod670_ss(pcols), dod670_dust(pcols), dod670_so4(pcols),     &
         dod670_bc(pcols), dod670_pom(pcols),                         &
         dod870_ss(pcols), dod870_dust(pcols), dod870_so4(pcols),     &
         dod870_bc(pcols), dod870_pom(pcols),                         &
         dod550lt1_ss(pcols), dod550gt1_ss(pcols), dod550lt1_dust(pcols), &
         dod550gt1_dust(pcols), dod550lt1_so4(pcols), &
         dod550gt1_so4(pcols), dod550lt1_bc(pcols), dod550gt1_bc(pcols), &
         dod550lt1_pom(pcols), dod550gt1_pom(pcols)

    real(r8) abs550_ss(pcols), abs550_dust(pcols), &
         abs550_so4(pcols), abs550_bc(pcols), abs550_pom(pcols)

    real(r8) batotsw13(pcols,pver), batotlw01(pcols,pver)
    character(len=10) :: modeString
    character(len=20) :: varname
    integer irf,irfmax
    real(r8) Camrel(pcols,pver,nbmodes)
    real(r8) Camtot(pcols,nbmodes)
    real(r8) cxsmtot(pcols,nbmodes)
    real(r8) cxsmrel(pcols,nbmodes)
    real(r8) xctrel,camdiff,cxsm
    real(r8) cxs(pcols,pver), cxstot(pcols,pver), akcxs(pcols)
    !-------------------------------------------------------------------------

    !     interpol-calculations only when daylight or not:
    do icol=1,ncol
       if (coszrs(icol) > 0.0_r8) then
          daylight(icol) = .true. 
       else
          daylight(icol) = .false.
       endif
    end do

    ! Initialize overshooting mass summed over all modes
    do k=1,pver
       do icol=1,ncol
          cxstot(icol,k)=0.0_r8
       enddo
    enddo
    do icol=1,ncol
       akcxs(icol)=0.0_r8
    enddo

    ! Initializing total and relative exessive (overshooting w.r.t.
    ! look-up table maxima) added mass column:
    do i=1,nbmodes
       do icol=1,ncol
          Camtot(icol,i)=0.0_r8
          cxsmtot(icol,i)=0.0_r8
          cxsmrel(icol,i)=0.0_r8
       enddo
    enddo

    ! Calculating added internally mixed mass onto each mode 1-10, relative to
    ! maximum mass which can be added w.r.t. the look-up tables (for level k),
    ! as well as the relative exessive added mass column:
    do i=1,4
       do k=1,pver
          do icol=1,ncol
             Camrel(icol,k,i) = (Cam(icol,k,i)/(Nnatk(icol,k,i)+eps))/cate(i,16)
             xctrel           = min(max(Camrel(icol,k,i),cate(i,1)/cate(i,16)),1.0_r8)
             camdiff          = Cam(icol,k,i)-xctrel*cate(i,16)*(Nnatk(icol,k,i)+eps)
             cxsm             = max(0.0_r8,camdiff)
             cxsmtot(icol,i)  = cxsmtot(icol,i)+cxsm*deltah_km(icol,k)
             Camtot(icol,i)   = Camtot(icol,i)+Cam(icol,k,i)*deltah_km(icol,k)
             camdiff          = Cam(icol,k,i)-xct(icol,k,i)*(Nnatk(icol,k,i)+eps)
             cxs(icol,k)      = max(0.0_r8,camdiff)
             cxstot(icol,k)   = cxstot(icol,k)+cxs(icol,k)
          enddo
       enddo
    enddo
    do i=5,nbmodes
       do k=1,pver
          do icol=1,ncol
             Camrel(icol,k,i) = (Cam(icol,k,i)/(Nnatk(icol,k,i)+eps))/cat(i,6)
             xctrel           = min(max(Camrel(icol,k,i),cat(i,1)/cat(i,6)),1.0_r8)
             camdiff          = Cam(icol,k,i)-xctrel*cat(i,6)*(Nnatk(icol,k,i)+eps)
             cxsm             = max(0.0_r8,camdiff)
             cxsmtot(icol,i)  = cxsmtot(icol,i)+cxsm*deltah_km(icol,k)
             Camtot(icol,i)   = Camtot(icol,i)+Cam(icol,k,i)*deltah_km(icol,k)
             camdiff          = Cam(icol,k,i)-xct(icol,k,i)*(Nnatk(icol,k,i)+eps)
             cxs(icol,k)      = max(0.0_r8,camdiff)
             cxstot(icol,k)   = cxstot(icol,k)+cxs(icol,k)
          enddo
       enddo
    enddo

    ! Total overshooting mass summed over all modes and all levels
    do icol=1,ncol
       do k=1,pver
          akcxs(icol) =akcxs(icol)+cxstot(icol,k)*deltah_km(icol,k)
       enddo
    enddo
    call outfld('AKCXS   ',akcxs ,pcols,lchnk)

    do i=1,nbmodes
       do icol=1,ncol
          cxsmrel(icol,i)=cxsmtot(icol,i)/(Camtot(icol,i)+eps)
       enddo
    enddo

    do i=1,nbmodes
       modeString="  "
       write(modeString,"(I2)"),i
       if(i.lt.10) modeString="0"//adjustl(modeString)
       varName = "Camrel"//trim(modeString)
       if(i.ne.3) call outfld(varName,Camrel(:,:,i),pcols,lchnk)
    enddo

    do i=1,nbmodes
       modeString="  "
       write(modeString,"(I2)"),i
       if(i.lt.10) modeString="0"//adjustl(modeString)
       varName = "Cxsrel"//trim(modeString)
       if(i.ne.3) call outfld(varName,cxsmrel(:,i),pcols,lchnk)
    enddo

    ! AeroCom: Find dry aerosol asymmetry factor and mass for subsequent
    ! calculation of condensed water mass below...
    do k=1,pver
       do icol=1,ncol
          Ctotdry(icol,k)=0.0_r8
          rh0(icol,k)=0.0_r8
          asydry_aer(icol,k)=0.0_r8
       end do
    enddo

    ! and define the respective RH input variables for dry aerosols
    do k=1,pver
       do icol=1,ncol
          xrhnull(icol,k)=rh(1)
          irh1null(icol,k)=1
       end do
    enddo

    !--------
    lw_on = .false.  ! No LW optics needed for RH=0 (interpol returns 0-values)
    !--------

    ! BC(ax) mode (dry only):
    call interpol0 (lchnk, ncol, daylight, Nnatk, ssa, asym, be, ke, lw_on, kalw)

    mplus10 = 0
    ! SO4/SOA(Ait) mode:
    call interpol1 (lchnk, ncol, daylight, xrhnull, irh1null, mplus10, &
         Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC(Ait) and OC(Ait) modes:
    call interpol2to3 (lchnk, ncol, daylight, xrhnull, irh1null, mplus10, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(Ait) mode:   ------ fcm not valid here (=0). Use faitbc instead
    call interpol4 (lchnk, ncol, daylight, xrhnull, irh1null, mplus10, &
         Nnatk, xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, &
         xfaq, ifaq1, ssa, asym, be, ke, lw_on, kalw)

    ! SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
    call interpol5to10 (lchnk, ncol, daylight, xrhnull, irh1null, &
         Nnatk, xct, ict1, xfac, ifac1, &
         xfbc, ifbc1, xfaq, ifaq1, &
         ssa, asym, be, ke, lw_on, kalw)

    mplus10 = 1
    ! BC(Ait) and OC(Ait) modes:
    call interpol2to3 (lchnk, ncol, daylight, xrhnull, irh1null, mplus10, &
         Nnatk, xct, ict1, xfac, ifac1, &
         ssa, asym, be, ke, lw_on, kalw)

    ! BC&OC(n) mode:   ------ fcm not valid here (=0). Use fnbc instead
    call interpol4 (lchnk, ncol, daylight, xrhnull, irh1null, mplus10, &
         Nnatk, xfbcbgn, ifbcbgn1, xct, ict1, &
         xfac, ifac1, xfaq, ifaq1, &
         ssa, asym, be, ke, lw_on, kalw)

    do i=0,nmodes    ! mode 0 to 14
       do k=1,pver
          do icol=1,ncol
             dCtot(icol,k)=1.e3_r8*be(icol,k,i,4)/(ke(icol,k,i,4)+eps)
             Ctotdry(icol,k)=Ctotdry(icol,k)+dCtot(icol,k)*Nnatk(icol,k,i)
          end do
       enddo
    enddo

    ! AeroCom Phase III: adding asymmetry factor for dry aerosol, wavelength band 4 only
    ! (and with no CMIP6 volcnic contribution)
    ib=4
    do k=1,pver
       do icol=1,ncol
          betot(icol,k,ib)=0.0_r8
          ssatot(icol,k,ib)=0.0_r8
          asymtot(icol,k,ib)=0.0_r8
       end do
    enddo
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
    do k=1,pver
       do icol=1,ncol
          ssatot(icol,k,ib)=ssatot(icol,k,ib)/(betot(icol,k,ib)+eps)
          asymtot(icol,k,ib)=asymtot(icol,k,ib) &
               /(betot(icol,k,ib)*ssatot(icol,k,ib)+eps)
          asydry_aer(icol,k)=asymtot(icol,k,ib)
       end do
    enddo
    call outfld('ASYMMDRY',asydry_aer,pcols,lchnk)

    !..................!

    !     Mass concentration (ug/m3) and mmr (kg/kg) of aerosol condensed water
    do k=1,pver
       do icol=1,ncol
          Cwater(icol,k)=Ctot(icol,k)-Ctotdry(icol,k)
          mmr_aerh2o(icol,k)=1.e-9_r8*Cwater(icol,k)/rhoda(icol,k)
       end do
    enddo

    !..................!

    do i=1,ncol 
       do k=1,pver 
          batotsw13(i,k)=betot(i,k,13)*(1.0_r8-ssatot(i,k,13))
          batotlw01(i,k)=batotlw(i,k,1)
       end do
    end do
    !    These two fields should be close to equal, both representing absorption 
    !    in the 3.077-3.846 um wavelenght band (i.e., a check of LUT for LW vs. SW).
    call outfld('BATSW13 ',batotsw13,pcols,lchnk)
    call outfld('BATLW01 ',batotlw01,pcols,lchnk)

    !..................!

    call outfld('BETOTVIS',betotvis,pcols,lchnk)
    call outfld('BATOTVIS',batotvis,pcols,lchnk)

    !       Initialize fields
    do icol=1,ncol
       daerh2o(icol)=0.0_r8
       vaercols(icol)=0.0_r8
       vaercoll(icol)=0.0_r8
       aaercols(icol)=0.0_r8
       aaercoll(icol)=0.0_r8
       do i=0,nmodes
          dload(icol,i)=0.0_r8
       enddo
    enddo
    vnbcarr(:,:)    = 0.0_r8
    vaitbcarr(:,:)  = 0.0_r8
    cknorm(:,:,:)   = 0.0_r8

    !     AeroCom diagnostics requiring table look-ups with ambient RH. 
    do irf=0,0
       call opticsAtConstRh(lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irf, &
            xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
            xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
            xfombg, ifombg1, vnbcarr, vaitbcarr, v_soana)
    end do ! irf

    do k=1,pver
       do icol=1,ncol

          bebglt1t(icol,k)=0.0_r8
          bebggt1t(icol,k)=0.0_r8
          bebclt1t(icol,k)=0.0_r8
          bebcgt1t(icol,k)=0.0_r8
          beoclt1t(icol,k)=0.0_r8
          beocgt1t(icol,k)=0.0_r8
          bes4lt1t(icol,k)=0.0_r8
          bes4gt1t(icol,k)=0.0_r8 
          bedustlt1(icol,k)=0.0_r8
          bedustgt1(icol,k)=0.0_r8
          besslt1(icol,k)=0.0_r8
          bessgt1(icol,k)=0.0_r8

          bext440tot(icol,k)=0.0_r8 
          babs440tot(icol,k)=0.0_r8 
          bext500tot(icol,k)=0.0_r8 
          babs500tot(icol,k)=0.0_r8 
          bext550tot(icol,k)=0.0_r8 
          babs550tot(icol,k)=0.0_r8 
          bext670tot(icol,k)=0.0_r8 
          babs670tot(icol,k)=0.0_r8 
          bext870tot(icol,k)=0.0_r8 
          babs870tot(icol,k)=0.0_r8 

          backsc550tot(icol,k)=0.0_r8 

          bebg440tot(icol,k)=0.0_r8 
          bebg500tot(icol,k)=0.0_r8 
          bebg550tot(icol,k)=0.0_r8 
          babg550tot(icol,k)=0.0_r8 
          bebg670tot(icol,k)=0.0_r8 
          bebg870tot(icol,k)=0.0_r8 

          bebc440tot(icol,k)=0.0_r8 
          bebc500tot(icol,k)=0.0_r8 
          bebc550tot(icol,k)=0.0_r8 
          babc550tot(icol,k)=0.0_r8 
          bebc670tot(icol,k)=0.0_r8 
          bebc870tot(icol,k)=0.0_r8 

          beoc440tot(icol,k)=0.0_r8 
          beoc500tot(icol,k)=0.0_r8 
          beoc550tot(icol,k)=0.0_r8 
          baoc550tot(icol,k)=0.0_r8 
          beoc670tot(icol,k)=0.0_r8 
          beoc870tot(icol,k)=0.0_r8 

          besu440tot(icol,k)=0.0_r8 
          besu500tot(icol,k)=0.0_r8 
          besu550tot(icol,k)=0.0_r8 
          basu550tot(icol,k)=0.0_r8 
          besu670tot(icol,k)=0.0_r8 
          besu870tot(icol,k)=0.0_r8 

       enddo
    enddo

    do i=0,nbmodes
       do k=1,pver
          do icol=1,ncol
             !      total internal extinction and absorption for 0.44, 0.50, 0.55, 0.68 and 0.87 um
             bext440tot(icol,k)=bext440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bext440(icol,k,i)
             babs440tot(icol,k)=babs440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs440(icol,k,i)
             bext500tot(icol,k)=bext500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bext500(icol,k,i)
             babs500tot(icol,k)=babs500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs500(icol,k,i)
             bext550tot(icol,k)=bext550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bext550(icol,k,i)
             babs550tot(icol,k)=babs550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs550(icol,k,i)
             bext670tot(icol,k)=bext670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bext670(icol,k,i)
             babs670tot(icol,k)=babs670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs670(icol,k,i)
             bext870tot(icol,k)=bext870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bext870(icol,k,i)
             babs870tot(icol,k)=babs870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs870(icol,k,i)
             backsc550tot(icol,k)=backsc550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%backsc550(icol,k,i)

             !      extinction and absorption for 0.44, 0.50, 0.55 (no abs), 0.68 and 0.87 um
             !      for the whole background aerosol (icluding SO4,BC, and OC for modes 0-5)
             bebg440tot(icol,k)=bebg440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebg440(icol,k,i)
             bebg500tot(icol,k)=bebg500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebg500(icol,k,i)
             bebg550tot(icol,k)=bebg550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebg550(icol,k,i)
             babg550tot(icol,k)=babg550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babg550(icol,k,i)
             bebg670tot(icol,k)=bebg670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebg670(icol,k,i)
             bebg870tot(icol,k)=bebg870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebg870(icol,k,i)
             besu440tot(icol,k)=besu440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%besu440(icol,k,i)
             besu500tot(icol,k)=besu500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%besu500(icol,k,i)
             besu550tot(icol,k)=besu550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%besu550(icol,k,i)
             basu550tot(icol,k)=basu550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%basu550(icol,k,i)
             besu670tot(icol,k)=besu670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%besu670(icol,k,i)
             besu870tot(icol,k)=besu870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%besu870(icol,k,i)
             !
             !      Condensed OC on modes 1-4 and coagulated BC and OC on modes 5-10:
             if(i>=1) then
                bebc440tot(icol,k)=bebc440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebc440(icol,k,i)
                bebc500tot(icol,k)=bebc500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebc500(icol,k,i)
                bebc550tot(icol,k)=bebc550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebc550(icol,k,i)
                babc550tot(icol,k)=babc550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babc550(icol,k,i)
                bebc670tot(icol,k)=bebc670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebc670(icol,k,i)
                bebc870tot(icol,k)=bebc870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%bebc870(icol,k,i)
                beoc440tot(icol,k)=beoc440tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%beoc440(icol,k,i)
                beoc500tot(icol,k)=beoc500tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%beoc500(icol,k,i)
                beoc550tot(icol,k)=beoc550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%beoc550(icol,k,i)
                baoc550tot(icol,k)=baoc550tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%baoc550(icol,k,i)
                beoc670tot(icol,k)=beoc670tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%beoc670(icol,k,i)
                beoc870tot(icol,k)=beoc870tot(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%beoc870(icol,k,i)
             endif  ! i>=1
             if(i==6.or.i==7) then
                bedustlt1(icol,k)=bedustlt1(icol,k) +Nnatk(icol,k,i)*bebglt1(icol,k,i)
                bedustgt1(icol,k)=bedustgt1(icol,k) +Nnatk(icol,k,i)*bebggt1(icol,k,i)
             elseif(i>=8.and.i<=10) then
                besslt1(icol,k)=besslt1(icol,k) +Nnatk(icol,k,i)*bebglt1(icol,k,i)
                bessgt1(icol,k)=bessgt1(icol,k) +Nnatk(icol,k,i)*bebggt1(icol,k,i)
             endif
             !      Condensed/coagulated SO4 on all modes 1-10, and wet-phase SO4 on modes 4-10:
             bes4lt1t(icol,k)=bes4lt1t(icol,k) +Nnatk(icol,k,i)*bes4lt1(icol,k,i)
             bes4gt1t(icol,k)=bes4gt1t(icol,k) +Nnatk(icol,k,i)*bes4gt1(icol,k,i)
             !      Condensed OC on mode 1 and coagulated BC and OC on modes 5-10:
             if(i>=1) then
                bebclt1t(icol,k)=bebclt1t(icol,k) +Nnatk(icol,k,i)*bebclt1(icol,k,i)
                bebcgt1t(icol,k)=bebcgt1t(icol,k) +Nnatk(icol,k,i)*bebcgt1(icol,k,i)
                beoclt1t(icol,k)=beoclt1t(icol,k) +Nnatk(icol,k,i)*beoclt1(icol,k,i)
                beocgt1t(icol,k)=beocgt1t(icol,k) +Nnatk(icol,k,i)*beocgt1(icol,k,i)
             endif   ! i>=1
          end do   ! icol
       enddo     ! k
    enddo      ! i

    !      extinction/absorptions (km-1) for each background component 
    !      in the internal mixture are
    do k=1,pver
       do icol=1,ncol
          bint440du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%bebg440(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%bebg440(icol,k,7)
          bint500du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%bebg500(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%bebg500(icol,k,7)
          bint550du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%bebg550(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%bebg550(icol,k,7)
          bint670du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%bebg670(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%bebg670(icol,k,7)
          bint870du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%bebg870(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%bebg870(icol,k,7)
          bint440ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%bebg440(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%bebg440(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%bebg440(icol,k,10)
          bint500ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%bebg500(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%bebg500(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%bebg500(icol,k,10)
          bint550ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%bebg550(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%bebg550(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%bebg550(icol,k,10)
          bint670ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%bebg670(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%bebg670(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%bebg670(icol,k,10)
          bint870ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%bebg870(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%bebg870(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%bebg870(icol,k,10)
          baint550du(icol,k)=Nnatk(icol,k,6)*extinction_coeffs%babg550(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%babg550(icol,k,7)
          baint550ss(icol,k)=Nnatk(icol,k,8)*extinction_coeffs%babg550(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%babg550(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%babg550(icol,k,10)
       end do
    enddo

    ! Need to make the following substitutions
    ! bebglt1 bebglt1n => extinction_coeffs%bebg550lt1
    ! bebggt1 bebggt1n => extinction_coeffs%bebg550gt1
    ! bebclt1 bebclt1n => extinction_coeffs%bebc550lt1
    ! bebcgt1 bebcgt1n => extinction_coeffs%bebc550gt1
    ! beoclt1 beoclt1n => extinction_coeffs%beoc550lt1
    ! beocgt1 beocgt1n => extinction_coeffs%beoc550gt1
    ! bes4lt1 bes4lt1n => extinction_coeffs%besu550lt1
    ! bes4gt1 bes4gt1n => extinction_coeffs%besu550gt1

    do i=11,14
       do k=1,pver
          do icol=1,ncol
             be440x(icol,k,i) = extinction_coeffsn%bext440(icol,k,i-10)
             ba440x(icol,k,i) = extinction_coeffsn%babs440(icol,k,i-10)
             be500x(icol,k,i) = extinction_coeffsn%bext500(icol,k,i-10)
             ba500x(icol,k,i) = extinction_coeffsn%babs500(icol,k,i-10)
             be550x(icol,k,i) = extinction_coeffsn%bext550(icol,k,i-10)
             ba550x(icol,k,i) = extinction_coeffsn%babs550(icol,k,i-10)
             be670x(icol,k,i) = extinction_coeffsn%bext670(icol,k,i-10)
             ba670x(icol,k,i) = extinction_coeffsn%babs670(icol,k,i-10)
             be870x(icol,k,i) = extinction_coeffsn%bext870(icol,k,i-10)
             ba870x(icol,k,i) = extinction_coeffsn%babs870(icol,k,i-10)
             belt1x(icol,k,i) = extinction_coeffsn%bebg550lt1(icol,k,i-10)
             begt1x(icol,k,i) = extinction_coeffsn%bebg550gt1(icol,k,i-10)
             backsc550x(icol,k,i) = extinction_coeffsn%backsc550(icol,k,i-10)
          end do
       enddo
    enddo

    !     The externally modes' contribution to extinction and absorption:
    do k=1,pver
       do icol=1,ncol

          !BC
          vnbcarr(icol,k) = fnbc(icol,k)/(fnbc(icol,k) &
               +(1.0_r8-fnbc(icol,k))*rhopart(l_bc_ni)/rhopart(l_om_ni))
          vnbc = vnbcarr(icol,k)            
          bebc440xt(icol,k) =Nnatk(icol,k,12)*be440x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*be440x(icol,k,14)
          babc440xt(icol,k) =Nnatk(icol,k,12)*ba440x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*ba440x(icol,k,14)
          bebc500xt(icol,k) =Nnatk(icol,k,12)*be500x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*be500x(icol,k,14)
          babc500xt(icol,k) =Nnatk(icol,k,12)*ba500x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*ba500x(icol,k,14)
          bebc550xt(icol,k) =Nnatk(icol,k,12)*be550x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*be550x(icol,k,14)
          babc550xt(icol,k) =Nnatk(icol,k,12)*ba550x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*ba550x(icol,k,14)
          bebc670xt(icol,k) =Nnatk(icol,k,12)*be670x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*be670x(icol,k,14)
          babc670xt(icol,k) =Nnatk(icol,k,12)*ba670x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*ba670x(icol,k,14)
          bebc870xt(icol,k) =Nnatk(icol,k,12)*be870x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*be870x(icol,k,14)
          babc870xt(icol,k) =Nnatk(icol,k,12)*ba870x(icol,k,12)  &
               +vnbc*Nnatk(icol,k,14)*ba870x(icol,k,14)
          bbclt1xt(icol,k)=Nnatk(icol,k,12)*belt1x(icol,k,12) &
               +vnbc*Nnatk(icol,k,14)*belt1x(icol,k,14)
          bbcgt1xt(icol,k)=Nnatk(icol,k,12)*begt1x(icol,k,12) &
               +vnbc*Nnatk(icol,k,14)*begt1x(icol,k,14)
          !OC
          beoc440xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be440x(icol,k,14) 
          baoc440xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba440x(icol,k,14) 
          beoc500xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be500x(icol,k,14) 
          baoc500xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba500x(icol,k,14) 
          beoc550xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be550x(icol,k,14) 
          baoc550xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba550x(icol,k,14) 
          beoc670xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be670x(icol,k,14) 
          baoc670xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba670x(icol,k,14) 
          beoc870xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*be870x(icol,k,14) 
          baoc870xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*ba870x(icol,k,14) 
          boclt1xt(icol,k) =  &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*belt1x(icol,k,14) 
          bocgt1xt(icol,k) = &
               +(1.0_r8-vnbc)*Nnatk(icol,k,14)*begt1x(icol,k,14) 
          !     Total (for all modes) absorption optical depth and backscattering
          abs550_aer(icol,k)=babs550tot(icol,k)  &
               +Nnatk(icol,k,12)*ba550x(icol,k,12) &
               +Nnatk(icol,k,14)*ba550x(icol,k,14)
          abs550_aer(icol,k)=1.e-3_r8*abs550_aer(icol,k)
          bs550_aer(icol,k)= backsc550tot(icol,k)   &
               +Nnatk(icol,k,12)*backsc550x(icol,k,12) &
               +Nnatk(icol,k,14)*backsc550x(icol,k,14)
          bs550_aer(icol,k)=1.e-3_r8*bs550_aer(icol,k)
          !
       end do
    enddo

    !       collect AeroCom-fields for optical depth/absorption of each comp,
    !       3D and 2D, at 440, 500, 550, 670 and 870 nm, for all d, d<1um and d>1um  
    !        initialize 2d-fields
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
       abs550_ss(icol) = 0.0_r8
       abs550_dust(icol) = 0.0_r8
       abs550_so4(icol) = 0.0_r8
       abs550_bc(icol) = 0.0_r8
       abs550_pom(icol) = 0.0_r8
       !
       dod440_ss(icol) = 0.0_r8
       dod440_dust(icol) = 0.0_r8
       dod440_so4(icol) = 0.0_r8
       dod440_bc(icol) = 0.0_r8
       dod440_pom(icol) = 0.0_r8
       dod500_ss(icol) = 0.0_r8
       dod500_dust(icol) = 0.0_r8
       dod500_so4(icol) = 0.0_r8
       dod500_bc(icol) = 0.0_r8
       dod500_pom(icol) = 0.0_r8
       dod550_ss(icol) = 0.0_r8
       dod550_dust(icol) = 0.0_r8
       dod550_so4(icol) = 0.0_r8
       dod550_bc(icol) = 0.0_r8
       dod550_pom(icol) = 0.0_r8
       dod670_ss(icol) = 0.0_r8
       dod670_dust(icol) = 0.0_r8
       dod670_so4(icol) = 0.0_r8
       dod670_bc(icol) = 0.0_r8
       dod670_pom(icol) = 0.0_r8
       dod870_ss(icol) = 0.0_r8
       dod870_dust(icol) = 0.0_r8
       dod870_so4(icol) = 0.0_r8
       dod870_bc(icol) = 0.0_r8
       dod870_pom(icol) = 0.0_r8
       dod550lt1_ss(icol) = 0.0_r8
       dod550gt1_ss(icol) = 0.0_r8
       dod550lt1_dust(icol) = 0.0_r8
       dod550gt1_dust(icol) = 0.0_r8
       dod550lt1_so4(icol) = 0.0_r8
       dod550gt1_so4(icol) = 0.0_r8
       dod550lt1_bc(icol) = 0.0_r8
       dod550gt1_bc(icol) = 0.0_r8
       dod550lt1_pom(icol) = 0.0_r8
       dod550gt1_pom(icol) = 0.0_r8
       do k=1,pver
          abs4403d(icol,k) = 0.0_r8
          abs5003d(icol,k) = 0.0_r8
          abs5503d(icol,k) = 0.0_r8
          abs6703d(icol,k) = 0.0_r8
          abs8703d(icol,k) = 0.0_r8
          abs5503dalt(icol,k) = 0.0_r8
       enddo
    enddo

    do icol=1,ncol
       do k=1,pver
          !          Layer thickness, unit km
          deltah=deltah_km(icol,k)
          !          if(k==pver) write(*,*) 'icol, deltah(pmxsub)=', icol, deltah
          !          3D optical depths for monthly averages
          !SS
          dod4403d_ss(icol,k) = bint440ss(icol,k)*deltah
          dod5003d_ss(icol,k) = bint500ss(icol,k)*deltah
          dod5503d_ss(icol,k) = bint550ss(icol,k)*deltah
          abs5503d_ss(icol,k) = baint550ss(icol,k)*deltah
          dod6703d_ss(icol,k) = bint670ss(icol,k)*deltah
          dod8703d_ss(icol,k) = bint870ss(icol,k)*deltah
          !DUST
          dod4403d_dust(icol,k) = bint440du(icol,k)*deltah
          dod5003d_dust(icol,k) = bint500du(icol,k)*deltah
          dod5503d_dust(icol,k) = bint550du(icol,k)*deltah
          abs5503d_dust(icol,k) = baint550du(icol,k)*deltah
          dod6703d_dust(icol,k) = bint670du(icol,k)*deltah
          dod8703d_dust(icol,k) = bint870du(icol,k)*deltah
          !SO4
          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          dod4403d_so4(icol,k) = (besu440tot(icol,k)                 &       ! condensate )
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%bebg440(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg440(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          dod5003d_so4(icol,k) = (besu500tot(icol,k)                 &       ! condensate 
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%bebg500(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg500(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          dod5503d_so4(icol,k) = (besu550tot(icol,k)                   &     ! condensate 
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%bebg550(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg550(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          abs5503d_so4(icol,k) = (basu550tot(icol,k)                 &       ! condensate )
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%babg550(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%babg550(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          dod6703d_so4(icol,k) = (besu670tot(icol,k)                   &     ! condensate
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%bebg670(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg670(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          dod8703d_so4(icol,k) = (besu870tot(icol,k)                   &     ! condensate 
               +(1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%bebg870(icol,k,1) &       ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg870(icol,k,5))*deltah ! background, SO4(Ait75) mode (5)
          !BC
          vaitbcarr(icol,k) = faitbc(icol,k)/(faitbc(icol,k) &
               +(1.0_r8-faitbc(icol,k))*rhopart(l_bc_ni)/rhopart(l_om_ni))
          vaitbc = vaitbcarr(icol,k)
          dod4403d_bc(icol,k) = (bebc440tot(icol,k)+bebc440xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg440(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%bebg440(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg440(icol,k,0))*deltah ! background, BC(ax) mode (0)
          dod5003d_bc(icol,k) = (bebc500tot(icol,k)+bebc500xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg500(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%bebg500(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg500(icol,k,0))*deltah ! background, BC(ax) mode (0)
          dod5503d_bc(icol,k) = (bebc550tot(icol,k)+bebc550xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg550(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%bebg550(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg550(icol,k,0))*deltah ! background, BC(ax) mode (0)
          abs5503d_bc(icol,k) = (babc550tot(icol,k)+babc550xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%babg550(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%babg550(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%babg550(icol,k,0))*deltah ! background, BC(ax) mode (0)
          dod6703d_bc(icol,k) = (bebc670tot(icol,k)+bebc670xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg670(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%bebg670(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg670(icol,k,0))*deltah ! background, BC(ax) mode (0)
          dod8703d_bc(icol,k) = (bebc870tot(icol,k)+bebc870xt(icol,k)  &     ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg870(icol,k,2) &       ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*extinction_coeffs%bebg870(icol,k,4) &       ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg870(icol,k,0))*deltah ! background, BC(ax) mode (0)
          !OC        
          !soa + v_soana part of mode 11 for the OC volume fraction of that mode
          ! v_soana(icol,k)
          dod4403d_pom(icol,k) = (beoc440tot(icol,k)+beoc440xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg440(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%bebg440(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
          dod5003d_pom(icol,k) = (beoc500tot(icol,k)+beoc500xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg500(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%bebg500(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
          dod5503d_pom(icol,k) = (beoc550tot(icol,k)+beoc550xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg550(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%bebg550(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
          abs5503d_pom(icol,k) = (baoc550tot(icol,k)+baoc550xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%babg550(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%babg550(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
          dod6703d_pom(icol,k) = (beoc670tot(icol,k)+beoc670xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg670(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%bebg670(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)
          dod8703d_pom(icol,k) = (beoc870tot(icol,k)+beoc870xt(icol,k) &     ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg870(icol,k,1)*v_soana(icol,k) & ! SOA fraction of mode 1
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*extinction_coeffs%bebg870(icol,k,4))*deltah ! background in OC&BC(Ait) mode (4)

          ec550_so4(icol,k) = 1.e-3*dod5503d_so4(icol,k)/deltah
          ec550_bc(icol,k)  = 1.e-3*dod5503d_bc(icol,k)/deltah
          ec550_pom(icol,k) = 1.e-3*dod5503d_pom(icol,k)/deltah
          ec550_ss(icol,k)  = 1.e-3*dod5503d_ss(icol,k)/deltah
          ec550_du(icol,k)  = 1.e-3*dod5503d_dust(icol,k)/deltah
          ec550_aer(icol,k) = ec550_so4(icol,k)+ec550_bc(icol,k)+ec550_pom(icol,k) &
               + ec550_ss(icol,k)+ec550_du(icol,k)

          !          Total 3D optical depths/abs. for column integrations
          dod4403d(icol,k) = dod4403d_ss(icol,k)+dod4403d_dust(icol,k) &
               +dod4403d_so4(icol,k)+dod4403d_bc(icol,k)    &
               +dod4403d_pom(icol,k)                     
          dod5003d(icol,k) = dod5003d_ss(icol,k)+dod5003d_dust(icol,k) &
               +dod5003d_so4(icol,k)+dod5003d_bc(icol,k)    &
               +dod5003d_pom(icol,k)                     
          dod5503d(icol,k) = dod5503d_ss(icol,k)+dod5503d_dust(icol,k) &
               +dod5503d_so4(icol,k)+dod5503d_bc(icol,k)    &
               +dod5503d_pom(icol,k)                     
          dod6703d(icol,k) = dod6703d_ss(icol,k)+dod6703d_dust(icol,k) &
               +dod6703d_so4(icol,k)+dod6703d_bc(icol,k)    &
               +dod6703d_pom(icol,k)                     
          dod8703d(icol,k) = dod8703d_ss(icol,k)+dod8703d_dust(icol,k) &
               +dod8703d_so4(icol,k)+dod8703d_bc(icol,k)    &
               +dod8703d_pom(icol,k)                     
          abs5503d(icol,k) = abs5503d_ss(icol,k)+abs5503d_dust(icol,k) &
               +abs5503d_so4(icol,k)+abs5503d_bc(icol,k)    &
               +abs5503d_pom(icol,k)                     
          !   (Note: Local abs550alt is up to 6% larger (annually averaged) in typical b.b.
          !   regions, compared to abs550. This is most likely most correct, but should be checked!)
          do i=0,10
             abs4403d(icol,k) = abs4403d(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs440(icol,k,i)*deltah
             abs5003d(icol,k) = abs5003d(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs500(icol,k,i)*deltah
             abs6703d(icol,k) = abs6703d(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs670(icol,k,i)*deltah
             abs8703d(icol,k) = abs8703d(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs870(icol,k,i)*deltah
             abs5503dalt(icol,k) = abs5503dalt(icol,k)+Nnatk(icol,k,i)*extinction_coeffs%babs550(icol,k,i)*deltah
          enddo
          do i=11,14
             abs4403d(icol,k) = abs4403d(icol,k)+Nnatk(icol,k,i)*extinction_coeffsn%babs440(icol,k,i-10)*deltah
             abs5003d(icol,k) = abs5003d(icol,k)+Nnatk(icol,k,i)*extinction_coeffsn%babs500(icol,k,i-10)*deltah
             abs6703d(icol,k) = abs6703d(icol,k)+Nnatk(icol,k,i)*extinction_coeffsn%babs670(icol,k,i-10)*deltah
             abs8703d(icol,k) = abs8703d(icol,k)+Nnatk(icol,k,i)*extinction_coeffsn%babs870(icol,k,i-10)*deltah
             abs5503dalt(icol,k) = abs5503dalt(icol,k)+Nnatk(icol,k,i)*extinction_coeffsn%babs550(icol,k,i-10)*deltah
          enddo

          !            optical depths for d<1um and d>1um (r<0.5um and r>0.5um)
          !SS
          dod5503dlt1_ss(icol,k) = besslt1(icol,k)*deltah
          dod5503dgt1_ss(icol,k) = bessgt1(icol,k)*deltah
          !DUST
          dod5503dlt1_dust(icol,k) = bedustlt1(icol,k)*deltah
          dod5503dgt1_dust(icol,k) = bedustgt1(icol,k)*deltah

          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          dod5503dlt1_so4(icol,k) = (bes4lt1t(icol,k)                  &   ! condensate
               + Nnatk(icol,k,1)*bebglt1(icol,k,1)*(1.0_r8-v_soana(icol,k)) &   ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*bebglt1(icol,k,5))*deltah             ! background, SO4(Ait75) mode (5)
          dod5503dgt1_so4(icol,k) = (bes4gt1t(icol,k)                  &   ! condensate + n-mode (11)
               + Nnatk(icol,k,1)*bebggt1(icol,k,1)*(1.0_r8-v_soana(icol,k)) &   ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*bebggt1(icol,k,5))*deltah             ! background, SO4(Ait75) mode (5)
          !BC
          dod5503dlt1_bc(icol,k) =  (bebclt1t(icol,k)+bbclt1xt(icol,k) &   ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*bebglt1(icol,k,2)                 &   ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*bebglt1(icol,k,4)                 &   ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*bebglt1(icol,k,0))*deltah             ! background, BC(ax) mode (0)
          dod5503dgt1_bc(icol,k) =  (bebcgt1t(icol,k)+bbcgt1xt(icol,k) &   ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*bebggt1(icol,k,2)                 &   ! background, BC(Ait) mode (2)
               + vaitbc*Nnatk(icol,k,4)*bebggt1(icol,k,4)                 &   ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*bebggt1(icol,k,0))*deltah             ! background, BC(ax) mode (0)
          !OC
          !soa + v_soana part of mode 11 for the OC volume fraction of that mode
          dod5503dlt1_pom(icol,k) = (beoclt1t(icol,k)+boclt1xt(icol,k) &   ! coagulated + n-mode OC&BC (14)
               + Nnatk(icol,k,1)*bebglt1(icol,k,1)*v_soana(icol,k) &   ! SOA fraction of mode 1
                                !-3                    + Nnatk(icol,k,3)*bebglt1(icol,k,3)                 &   ! background, OC(Ait) mode (3)
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebglt1(icol,k,4))*deltah          ! background in OC&BC(Ait) mode (4)
          dod5503dgt1_pom(icol,k) = (beocgt1t(icol,k)+bocgt1xt(icol,k) &   ! coagulated + n-mode OC&OC (14)
               + Nnatk(icol,k,1)*bebggt1(icol,k,1)*v_soana(icol,k) &   ! SOA fraction of mode 1
                                !-3                    + Nnatk(icol,k,3)*bebggt1(icol,k,3)                 &   ! background, OC(Ait) mode (3)
               + (1.0_r8-vaitbc)*Nnatk(icol,k,4)*bebggt1(icol,k,4))*deltah          ! background in OC&BC(Ait) mode (4)

          !          Column integrated optical depths/abs., total and for each constituent
          dod440(icol)         = dod440(icol)+dod4403d(icol,k)
          abs440(icol)         = abs440(icol)+abs4403d(icol,k)
          dod500(icol)         = dod500(icol)+dod5003d(icol,k)
          abs500(icol)         = abs500(icol)+abs5003d(icol,k)
          dod550(icol)         = dod550(icol)+dod5503d(icol,k)
          abs550(icol)         = abs550(icol)+abs5503d(icol,k)
          abs550alt(icol)      = abs550alt(icol)+abs5503dalt(icol,k)
          dod670(icol)         = dod670(icol)+dod6703d(icol,k)
          abs670(icol)         = abs670(icol)+abs6703d(icol,k)
          dod870(icol)         = dod870(icol)+dod8703d(icol,k)
          abs870(icol)         = abs870(icol)+abs8703d(icol,k)
          ! Added abs components
          abs550_ss(icol)      = abs550_ss(icol)+abs5503d_ss(icol,k)
          abs550_dust(icol)    = abs550_dust(icol)+abs5503d_dust(icol,k)
          abs550_so4(icol)     = abs550_so4(icol)+abs5503d_so4(icol,k)
          abs550_bc(icol)      = abs550_bc(icol)+abs5503d_bc(icol,k)
          abs550_pom(icol)     = abs550_pom(icol)+abs5503d_pom(icol,k)
          !
          dod440_ss(icol)      = dod440_ss(icol)+dod4403d_ss(icol,k)
          dod440_dust(icol)    = dod440_dust(icol)+dod4403d_dust(icol,k)
          dod440_so4(icol)     = dod440_so4(icol)+dod4403d_so4(icol,k)
          dod440_bc(icol)      = dod440_bc(icol)+dod4403d_bc(icol,k)
          dod440_pom(icol)     = dod440_pom(icol)+dod4403d_pom(icol,k)
          dod500_ss(icol)      = dod500_ss(icol)+dod5003d_ss(icol,k)
          dod500_dust(icol)    = dod500_dust(icol)+dod5003d_dust(icol,k)
          dod500_so4(icol)     = dod500_so4(icol)+dod5003d_so4(icol,k)
          dod500_bc(icol)      = dod500_bc(icol)+dod5003d_bc(icol,k)
          dod500_pom(icol)     = dod500_pom(icol)+dod5003d_pom(icol,k)
          dod550_ss(icol)      = dod550_ss(icol)+dod5503d_ss(icol,k)
          dod550_dust(icol)    = dod550_dust(icol)+dod5503d_dust(icol,k)
          dod550_so4(icol)     = dod550_so4(icol)+dod5503d_so4(icol,k)
          dod550_bc(icol)      = dod550_bc(icol)+dod5503d_bc(icol,k)
          dod550_pom(icol)     = dod550_pom(icol)+dod5503d_pom(icol,k)
          dod670_ss(icol)      = dod670_ss(icol)+dod6703d_ss(icol,k)
          dod670_dust(icol)    = dod670_dust(icol)+dod6703d_dust(icol,k)
          dod670_so4(icol)     = dod670_so4(icol)+dod6703d_so4(icol,k)
          dod670_bc(icol)      = dod670_bc(icol)+dod6703d_bc(icol,k)
          dod670_pom(icol)     = dod670_pom(icol)+dod6703d_pom(icol,k)
          dod870_ss(icol)      = dod870_ss(icol)+dod8703d_ss(icol,k)
          dod870_dust(icol)    = dod870_dust(icol)+dod8703d_dust(icol,k)
          dod870_so4(icol)     = dod870_so4(icol)+dod8703d_so4(icol,k)
          dod870_bc(icol)      = dod870_bc(icol)+dod8703d_bc(icol,k)
          dod870_pom(icol)     = dod870_pom(icol)+dod8703d_pom(icol,k)
          dod550lt1_ss(icol)   = dod550lt1_ss(icol)+dod5503dlt1_ss(icol,k)
          dod550gt1_ss(icol)   = dod550gt1_ss(icol)+dod5503dgt1_ss(icol,k)
          dod550lt1_dust(icol) = dod550lt1_dust(icol)+dod5503dlt1_dust(icol,k)
          dod550gt1_dust(icol) = dod550gt1_dust(icol)+dod5503dgt1_dust(icol,k)
          dod550lt1_so4(icol)  = dod550lt1_so4(icol)+dod5503dlt1_so4(icol,k) 
          dod550gt1_so4(icol)  = dod550gt1_so4(icol)+dod5503dgt1_so4(icol,k) 
          dod550lt1_bc(icol)   = dod550lt1_bc(icol)+dod5503dlt1_bc(icol,k)  
          dod550gt1_bc(icol)   = dod550gt1_bc(icol)+dod5503dgt1_bc(icol,k)  
          dod550lt1_pom(icol)  = dod550lt1_pom(icol)+dod5503dlt1_pom(icol,k) 
          dod550gt1_pom(icol)  = dod550gt1_pom(icol)+dod5503dgt1_pom(icol,k)
       enddo ! k

    enddo  ! icol

    !       extinction, absorption (m-1) and backscatter coefficients (m-1 sr-1)
    call outfld('EC550AER',ec550_aer,pcols,lchnk)
    call outfld('ABS550_A',abs550_aer,pcols,lchnk)
    call outfld('BS550AER',bs550_aer,pcols,lchnk)
    !
    !       speciated extinction coefficients (m-1)
    call outfld('EC550SO4',ec550_so4,pcols,lchnk)
    call outfld('EC550BC ',ec550_bc ,pcols,lchnk)
    call outfld('EC550POM',ec550_pom,pcols,lchnk)
    call outfld('EC550SS ',ec550_ss ,pcols,lchnk)
    call outfld('EC550DU ',ec550_du ,pcols,lchnk)
    !
    !       optical depths and absorption as requested by AeroCom
    !       notation: 3=3D, D=DOD, A=ABS, LT=d<1um, GT=d>1um
    call outfld('DOD440  ',dod440  ,pcols,lchnk)
    call outfld('ABS440  ',abs440  ,pcols,lchnk)
    call outfld('DOD500  ',dod500  ,pcols,lchnk)
    call outfld('ABS500  ',abs500  ,pcols,lchnk)
    call outfld('DOD550  ',dod550  ,pcols,lchnk)
    call outfld('ABS550  ',abs550  ,pcols,lchnk)
    call outfld('ABS550AL',abs550alt,pcols,lchnk)
    call outfld('DOD670  ',dod670  ,pcols,lchnk)
    call outfld('ABS670  ',abs670  ,pcols,lchnk)
    call outfld('DOD870  ',dod870  ,pcols,lchnk)
    call outfld('ABS870  ',abs870  ,pcols,lchnk)
    call outfld('A550_SS ',abs550_ss  ,pcols,lchnk)
    call outfld('A550_DU ',abs550_dust,pcols,lchnk)
    call outfld('A550_SO4',abs550_so4 ,pcols,lchnk)
    call outfld('A550_BC ',abs550_bc  ,pcols,lchnk)
    call outfld('A550_POM',abs550_pom ,pcols,lchnk)
    !
    call outfld('D440_SS ',dod440_ss  ,pcols,lchnk)
    call outfld('D440_DU ',dod440_dust,pcols,lchnk)
    call outfld('D440_SO4',dod440_so4 ,pcols,lchnk)
    call outfld('D440_BC ',dod440_bc  ,pcols,lchnk)
    call outfld('D440_POM',dod440_pom ,pcols,lchnk)
    call outfld('D500_SS ',dod500_ss  ,pcols,lchnk)
    call outfld('D500_DU ',dod500_dust,pcols,lchnk)
    call outfld('D500_SO4',dod500_so4 ,pcols,lchnk)
    call outfld('D500_BC ',dod500_bc  ,pcols,lchnk)
    call outfld('D500_POM',dod500_pom ,pcols,lchnk)
    call outfld('D550_SS ',dod550_ss  ,pcols,lchnk)
    call outfld('D550_DU ',dod550_dust,pcols,lchnk)
    call outfld('D550_SO4',dod550_so4 ,pcols,lchnk)
    call outfld('D550_BC ',dod550_bc  ,pcols,lchnk)
    call outfld('D550_POM',dod550_pom ,pcols,lchnk)
    call outfld('D670_SS ',dod670_ss  ,pcols,lchnk)
    call outfld('D670_DU ',dod670_dust,pcols,lchnk)
    call outfld('D670_SO4',dod670_so4 ,pcols,lchnk)
    call outfld('D670_BC ',dod670_bc  ,pcols,lchnk)
    call outfld('D670_POM',dod670_pom ,pcols,lchnk)
    call outfld('D870_SS ',dod870_ss  ,pcols,lchnk)
    call outfld('D870_DU ',dod870_dust,pcols,lchnk)
    call outfld('D870_SO4',dod870_so4 ,pcols,lchnk)
    call outfld('D870_BC ',dod870_bc  ,pcols,lchnk)
    call outfld('D870_POM',dod870_pom ,pcols,lchnk)
    call outfld('DLT_SS  ',dod550lt1_ss,pcols,lchnk)
    call outfld('DGT_SS  ',dod550gt1_ss,pcols,lchnk)
    call outfld('DLT_DUST',dod550lt1_dust,pcols,lchnk)
    call outfld('DGT_DUST',dod550gt1_dust,pcols,lchnk)
    call outfld('DLT_SO4 ',dod550lt1_so4,pcols,lchnk)
    call outfld('DGT_SO4 ',dod550gt1_so4,pcols,lchnk)
    call outfld('DLT_BC  ',dod550lt1_bc,pcols,lchnk)
    call outfld('DGT_BC  ',dod550gt1_bc,pcols,lchnk)
    call outfld('DLT_POM ',dod550lt1_pom,pcols,lchnk)
    call outfld('DGT_POM ',dod550gt1_pom,pcols,lchnk)
    !       Dry parameters of each aerosol component 
    !       BC(ax) mode
    call aerodry_prop%intdrypar0(lchnk, ncol, Nnatk)

    !       SO4&SOA(Ait,n) mode
    call aerodry_prop%intdrypar1(lchnk, ncol, Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1)

    !       BC(Ait,n) and OC(Ait,n) modes
    call aerodry_prop%intdrypar2to3(lchnk, ncol, Nnatk, xct, ict1, xfac, ifac1)

    !       BC&OC(Ait,n) mode   ------ fcm not valid here (=0). Use faitbc or fnbc instead
    call aerodry_prop%intdrypar4(lchnk, ncol, Nnatk, xfbcbg, ifbcbg1, xfbcbgn, ifbcbgn1, &
         xct, ict1, xfac, ifac1, xfaq, ifaq1)

    !       SO4(Ait75) (5), mineral (6-7) and Sea-salt (8-10) modes:
    call aerodry_prop%intdrypar5to10(lchnk, ncol, Nnatk, xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1)

    do k=1,pver
       do icol=1,ncol
          c_ss(icol,k)=0.0_r8
          c_mi(icol,k)=0.0_r8
       enddo
    enddo

    do k=1,pver
       do icol=1,ncol
          ! mineral and sea-salt background concentrations, internally mixed
          c_mi(icol,k)    = Nnatk(icol,k,6) * aerodry_prop%cintbg(icol,k,6)   &
               +Nnatk(icol,k,7) * aerodry_prop%cintbg(icol,k,7)
          c_mi05(icol,k)  = Nnatk(icol,k,6) * aerodry_prop%cintbg05(icol,k,6) &
               +Nnatk(icol,k,7) * aerodry_prop%cintbg05(icol,k,7)
          c_mi125(icol,k) = Nnatk(icol,k,6) * aerodry_prop%cintbg125(icol,k,6)& 
               +Nnatk(icol,k,7) * aerodry_prop%cintbg125(icol,k,7) 
          c_ss(icol,k)    = Nnatk(icol,k,8) * aerodry_prop%cintbg(icol,k,8)   &
               +Nnatk(icol,k,9) * aerodry_prop%cintbg(icol,k,9)    &
               +Nnatk(icol,k,10) * aerodry_prop%cintbg(icol,k,10)
          c_ss05(icol,k)  = Nnatk(icol,k,8) * aerodry_prop%cintbg05(icol,k,8) &
               +Nnatk(icol,k,9) * aerodry_prop%cintbg05(icol,k,9)  &
               +Nnatk(icol,k,10) * aerodry_prop%cintbg05(icol,k,10)
          c_ss125(icol,k) = Nnatk(icol,k,8) * aerodry_prop%cintbg125(icol,k,8)&
               +Nnatk(icol,k,9) * aerodry_prop%cintbg125(icol,k,9) &
               +Nnatk(icol,k,10) * aerodry_prop%cintbg125(icol,k,10)

          !           internally mixed bc and oc (from coagulation) and so4 concentrations 
          !           (sa=so4(aq) and sc=so4(cond+coag), separated because of different density: 
          !           necessary for calculation of volume fractions!), and total aerosol surface 
          !           areas and volumes. 
          c_bc(icol,k)=0.0_r8
          c_bc05(icol,k)=0.0_r8
          c_bc125(icol,k)=0.0_r8
          c_oc(icol,k)=0.0_r8
          c_oc05(icol,k)=0.0_r8
          c_oc125(icol,k)=0.0_r8
          c_s4(icol,k)=0.0_r8
          c_s4_a(icol,k)=0.0_r8
          c_s4_1(icol,k)=0.0_r8
          c_s4_5(icol,k)=0.0_r8
          c_sa(icol,k)=0.0_r8
          c_sa05(icol,k)=0.0_r8
          c_sa125(icol,k)=0.0_r8
          c_sc(icol,k)=0.0_r8
          c_sc05(icol,k)=0.0_r8
          c_sc125(icol,k)=0.0_r8
          aaeros_tot(icol,k)=0.0_r8
          aaerol_tot(icol,k)=0.0_r8
          vaeros_tot(icol,k)=0.0_r8
          vaerol_tot(icol,k)=0.0_r8
          c_bc_0(icol,k)=0.0_r8
          c_bc_2(icol,k)=0.0_r8
          c_bc_4(icol,k)=0.0_r8
          c_bc_12(icol,k)=0.0_r8
          c_bc_14(icol,k)=0.0_r8
          c_oc_4(icol,k)=0.0_r8
          c_oc_14(icol,k)=0.0_r8
          c_tot(icol,k)=0.0_r8
          c_tot125(icol,k)=0.0_r8
          c_tot05(icol,k)=0.0_r8
          c_pm25(icol,k)=0.0_r8
          c_pm1(icol,k)=0.0_r8
          mmr_pm25(icol,k)=0.0_r8
          mmr_pm1(icol,k)=0.0_r8

          do i=0,nbmodes
             if(i.ne.3) then
                c_bc(icol,k)       = c_bc(icol,k)       + Nnatk(icol,k,i) * aerodry_prop%cintbc(icol,k,i)
                c_bc05(icol,k)     = c_bc05(icol,k)     + Nnatk(icol,k,i) * aerodry_prop%cintbc05(icol,k,i)
                c_bc125(icol,k)    = c_bc125(icol,k)    + Nnatk(icol,k,i) * aerodry_prop%cintbc125(icol,k,i)
                c_oc(icol,k)       = c_oc(icol,k)       + Nnatk(icol,k,i) * aerodry_prop%cintoc(icol,k,i)
                c_oc05(icol,k)     = c_oc05(icol,k)     + Nnatk(icol,k,i) * aerodry_prop%cintoc05(icol,k,i)
                c_oc125(icol,k)    = c_oc125(icol,k)    + Nnatk(icol,k,i) * aerodry_prop%cintoc125(icol,k,i)
                c_sa(icol,k)       = c_sa(icol,k)       + Nnatk(icol,k,i) * aerodry_prop%cintsa(icol,k,i)
                c_sa05(icol,k)     = c_sa05(icol,k)     + Nnatk(icol,k,i) * aerodry_prop%cintsa05(icol,k,i)
                c_sa125(icol,k)    = c_sa125(icol,k)    + Nnatk(icol,k,i) * aerodry_prop%cintsa125(icol,k,i)
                c_sc(icol,k)       = c_sc(icol,k)       + Nnatk(icol,k,i) * aerodry_prop%cintsc(icol,k,i)
                c_sc05(icol,k)     = c_sc05(icol,k)     + Nnatk(icol,k,i) * aerodry_prop%cintsc05(icol,k,i)
                c_sc125(icol,k)    = c_sc125(icol,k)    + Nnatk(icol,k,i) * aerodry_prop%cintsc125(icol,k,i)
                aaeros_tot(icol,k) = aaeros_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%aaeros(icol,k,i)
                aaerol_tot(icol,k) = aaerol_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%aaerol(icol,k,i)
                vaeros_tot(icol,k) = vaeros_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%vaeros(icol,k,i)
                vaerol_tot(icol,k) = vaerol_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%vaerol(icol,k,i)
             endif
          enddo
          !           add dry aerosol area and volume of externally mixed modes
          do i=nbmp1,nmodes
             aaeros_tot(icol,k) = aaeros_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%aaerosn(icol,k,i)
             aaerol_tot(icol,k) = aaerol_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%aaeroln(icol,k,i)
             vaeros_tot(icol,k) = vaeros_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%vaerosn(icol,k,i)
             vaerol_tot(icol,k) = vaerol_tot(icol,k) + Nnatk(icol,k,i) * aerodry_prop%vaeroln(icol,k,i)
          end do

          !c_er3d           
          !           Effective radii for particles smaller and greater than 0.5um, 
          !           and for all radii, in each layer (er=3*V/A):
          erlt053d(icol,k)=3.0_r8*vaeros_tot(icol,k) /(aaeros_tot(icol,k)+eps)
          ergt053d(icol,k)=3.0_r8*vaerol_tot(icol,k) /(aaerol_tot(icol,k)+eps)
          er3d(icol,k)=3.0_r8*(vaeros_tot(icol,k)+vaerol_tot(icol,k)) /(aaeros_tot(icol,k)+aaerol_tot(icol,k)+eps)

          !c_er3d
          !           column integrated dry aerosol surface areas and volumes
          !           for r<0.5um and r>0.5um (s and l, respectively).
          aaercols(icol)=aaercols(icol)+aaeros_tot(icol,k)
          aaercoll(icol)=aaercoll(icol)+aaerol_tot(icol,k)
          vaercols(icol)=vaercols(icol)+vaeros_tot(icol,k)
          vaercoll(icol)=vaercoll(icol)+vaerol_tot(icol,k)

          !           then add background and externally mixed BC, OC and SO4 to mass concentrations
          c_bc_ac(icol,k)= c_bc(icol,k)
          c_bc_0(icol,k) = Nnatk(icol,k,0) * aerodry_prop%cintbg(icol,k,0)
          c_bc_2(icol,k) = Nnatk(icol,k,2) * aerodry_prop%cintbg(icol,k,2)
          c_bc_4(icol,k) = Nnatk(icol,k,4) * aerodry_prop%cintbg(icol,k,4)*faitbc(icol,k)
          c_bc_12(icol,k)= Nnatk(icol,k,12) * aerodry_prop%cknorm(icol,k,12)
          c_bc_14(icol,k)= Nnatk(icol,k,14) * aerodry_prop%cknorm(icol,k,14)*fnbc(icol,k)
          c_bc(icol,k)   = c_bc(icol,k)                                      &
               +Nnatk(icol,k,2) * aerodry_prop%cintbg(icol,k,2)                   &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg(icol,k,4) * faitbc(icol,k) &
               +Nnatk(icol,k,0) * aerodry_prop%cintbg(icol,k,0)                   &
               +Nnatk(icol,k,12) * aerodry_prop%cknorm(icol,k,12)                 &
               +Nnatk(icol,k,14) * aerodry_prop%cknorm(icol,k,14)*fnbc(icol,k)
          c_bc05(icol,k)  = c_bc05(icol,k)                                   &
               +Nnatk(icol,k,2) * aerodry_prop%cintbg05(icol,k,2)                 &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg05(icol,k,4)*faitbc(icol,k) &
               +Nnatk(icol,k,0) * aerodry_prop%cintbg05(icol,k,0)                 &
               +Nnatk(icol,k,12) * aerodry_prop%cknlt05(icol,k,12)                &
               +Nnatk(icol,k,14) * aerodry_prop%cknlt05(icol,k,14)*fnbc(icol,k)
          c_bc125(icol,k) = c_bc125(icol,k)                                  &
               +Nnatk(icol,k,2) * aerodry_prop%cintbg125(icol,k,2)                &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg125(icol,k,4)*faitbc(icol,k) &
               +Nnatk(icol,k,0) * aerodry_prop%cintbg125(icol,k,0)                &
               +Nnatk(icol,k,12) * aerodry_prop%ckngt125(icol,k,12)               &
               +Nnatk(icol,k,14) * aerodry_prop%ckngt125(icol,k,14)*fnbc(icol,k)
          c_oc_ac(icol,k)= c_oc(icol,k)
          c_oc_4(icol,k)  = Nnatk(icol,k,4) * aerodry_prop%cintbg(icol,k,4)*(1.0_r8-faitbc(icol,k))
          c_oc_14(icol,k) = Nnatk(icol,k,14) * aerodry_prop%cknorm(icol,k,14)*(1.0_r8-fnbc(icol,k))
          c_oc(icol,k)    = c_oc(icol,k)                                           &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg(icol,k,1)*f_soana(icol,k)         &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg(icol,k,4)*(1.0_r8-faitbc(icol,k))   &
               +Nnatk(icol,k,14) * aerodry_prop%cknorm(icol,k,14)*(1.0_r8-fnbc(icol,k))
          c_oc05(icol,k)  = c_oc05(icol,k)                                         &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg05(icol,k,1)*f_soana(icol,k)       &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg05(icol,k,4)*(1.0_r8-faitbc(icol,k))  &
               +Nnatk(icol,k,14) * aerodry_prop%cknlt05(icol,k,14)*(1.0_r8-fnbc(icol,k))
          c_oc125(icol,k) = c_oc125(icol,k)                                        &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg125(icol,k,1)*f_soana(icol,k)      &
               +Nnatk(icol,k,4) * aerodry_prop%cintbg125(icol,k,4)*(1.0_r8-faitbc(icol,k)) &
               +Nnatk(icol,k,14) * aerodry_prop%ckngt125(icol,k,14)*(1.0_r8-fnbc(icol,k))
          c_s4(icol,k)    = c_sa(icol,k)+c_sc(icol,k)          &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg(icol,k,1)*(1.0_r8-f_soana(icol,k))   &
               +Nnatk(icol,k,5) * aerodry_prop%cintbg(icol,k,5)     
          c_s405(icol,k)  = c_sa05(icol,k)+c_sc05(icol,k)      &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg05(icol,k,1)*(1.0_r8-f_soana(icol,k)) &
               +Nnatk(icol,k,5) * aerodry_prop%cintbg05(icol,k,5)    
          c_s4125(icol,k) = c_sa125(icol,k)+c_sc125(icol,k)    &
               +Nnatk(icol,k,1) * aerodry_prop%cintbg125(icol,k,1)*(1.0_r8-f_soana(icol,k)) &
               +Nnatk(icol,k,5) * aerodry_prop%cintbg125(icol,k,5)  

          c_tot(icol,k)  = c_s4(icol,k) + c_oc(icol,k) + c_bc(icol,k) + c_mi(icol,k) + c_ss(icol,k)
          c_tot125(icol,k) = c_s4125(icol,k) + c_oc125(icol,k) + c_bc125(icol,k) + c_mi125(icol,k) + c_ss125(icol,k)
          c_tot05(icol,k) = c_s405(icol,k) + c_oc05(icol,k) + c_bc05(icol,k) + c_mi05(icol,k) + c_ss05(icol,k)
          c_pm25(icol,k) = c_tot(icol,k) - c_tot125(icol,k)
          c_pm1(icol,k) = c_tot05(icol,k)

          !  mass mixing ratio:
          mmr_pm25(icol,k) = 1.e-9*c_pm25(icol,k)/rhoda(icol,k)   
          mmr_pm1(icol,k)  = 1.e-9*c_pm1(icol,k)/rhoda(icol,k)   

          !            converting from S to SO4 concentrations is no longer necessary, since 
          !             sc=H2SO4 and sa=(NH4)2SO4 now, not SO4 as in CAM4-Oslo     
          !             c_s4(icol,k)=c_s4(icol,k)/3._r8
          !             c_s405(icol,k)=c_s405(icol,k)/3._r8
          !             c_s4125(icol,k)=c_s4125(icol,k)/3._r8

          c_s4_a(icol,k) = c_sa(icol,k)+c_sc(icol,k) 
          c_s4_1(icol,k) = Nnatk(icol,k,1) * aerodry_prop%cintbg(icol,k,1)*(1.0_r8-f_soana(icol,k))
          c_s4_5(icol,k) = Nnatk(icol,k,5) * aerodry_prop%cintbg05(icol,k,5) 

       end do ! icol
    enddo     ! k

    !       Total PM and PM2.5 (dry r>1.25um), surface values (ug/m3)
    do icol=1,ncol
       c_tots(icol) = c_tot(icol,pver)
       c_tot125s(icol) = c_tot125(icol,pver)
       c_pm25s(icol) = c_pm25(icol,pver)
    enddo

    !       Effective, column integrated, radii for particles
    !       smaller and greater than 0.5um, and for all radii
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
       do k=1,pver
          !         Layer thickness, unit km
          !-          deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
          deltah=deltah_km(icol,k)
          !         Modal and total mass concentrations for clean and dry aerosol, 
          !         i.e. not including coag./cond./Aq. BC,OC,SO4 or condensed water. 
          !         Units: ug/m3 for concentrations and mg/m2 (--> kg/m2 later) for mass loading.  
          do i=0,nmodes
             ck(icol,k,i)=cknorm(icol,k,i)*Nnatk(icol,k,i)
             dload3d(icol,k,i)=ck(icol,k,i)*deltah
             dload(icol,i)=dload(icol,i)+dload3d(icol,k,i)
          enddo
          nnat_0(icol,k) =Nnatk(icol,k,0)
          nnat_1(icol,k) =Nnatk(icol,k,1)
          nnat_2(icol,k) =Nnatk(icol,k,2)
          nnat_4(icol,k) =Nnatk(icol,k,4)
          nnat_5(icol,k) =Nnatk(icol,k,5)
          nnat_6(icol,k) =Nnatk(icol,k,6)
          nnat_7(icol,k) =Nnatk(icol,k,7)
          nnat_8(icol,k) =Nnatk(icol,k,8)
          nnat_9(icol,k) =Nnatk(icol,k,9)
          nnat_10(icol,k)=Nnatk(icol,k,10)
          nnat_12(icol,k)=Nnatk(icol,k,12)
          nnat_14(icol,k)=Nnatk(icol,k,14)
          !         mineral and sea-salt mass concentrations
          cmin(icol,k)=ck(icol,k,6)+ck(icol,k,7)               
          cseas(icol,k)=ck(icol,k,8)+ck(icol,k,9)+ck(icol,k,10)
          !          Aerocom: Condensed water loading (mg_m2)
          daerh2o(icol)=daerh2o(icol)+Cwater(icol,k)*deltah
          !         just for checking purposes:
          dload_s4(icol)=dload_s4(icol)+c_s4(icol,k)*deltah
          dload_s4_a(icol)=dload_s4_a(icol)+c_s4_a(icol,k)*deltah
          dload_s4_1(icol)=dload_s4_1(icol)+c_s4_1(icol,k)*deltah
          dload_s4_5(icol)=dload_s4_5(icol)+c_s4_5(icol,k)*deltah
          dload_oc(icol)=dload_oc(icol)+c_oc(icol,k)*deltah
          dload_bc(icol)=dload_bc(icol)+c_bc(icol,k)*deltah
          !
          dload_bc_ac(icol)=dload_bc_ac(icol)+c_bc_ac(icol,k)*deltah
          dload_bc_0(icol)=dload_bc_0(icol)+c_bc_0(icol,k)*deltah
          dload_bc_2(icol)=dload_bc_2(icol)+c_bc_2(icol,k)*deltah
          dload_bc_4(icol)=dload_bc_4(icol)+c_bc_4(icol,k)*deltah
          dload_bc_12(icol)=dload_bc_12(icol)+c_bc_12(icol,k)*deltah
          dload_bc_14(icol)=dload_bc_14(icol)+c_bc_14(icol,k)*deltah
          dload_oc_ac(icol)=dload_oc_ac(icol)+c_oc_ac(icol,k)*deltah
          dload_oc_4(icol)=dload_oc_4(icol)+c_oc_4(icol,k)*deltah
          dload_oc_14(icol)=dload_oc_14(icol)+c_oc_14(icol,k)*deltah
          !
       end do  ! k
       dload_mi(icol)=dload(icol,6)+dload(icol,7)
       dload_ss(icol)=dload(icol,8)+dload(icol,9)+dload(icol,10)
    end do   ! icol

    call outfld('PMTOT  ',c_tots   ,pcols,lchnk)
    call outfld('PM25    ',c_pm25s ,pcols,lchnk)
    call outfld('PM2P5   ',c_pm25  ,pcols,lchnk)
    call outfld('MMRPM2P5',mmr_pm25,pcols,lchnk)
    call outfld('MMRPM1  ',mmr_pm1 ,pcols,lchnk)
    call outfld('MMRPM2P5_SRF',mmr_pm25(:pcols,pver),pcols,lchnk) 
    !       total (all r) dry concentrations (ug/m3) and loadings (mg/m2)
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
    !       condensed water mmr (kg/kg)
    call outfld('MMR_AH2O',mmr_aerh2o,pcols,lchnk)
    !       condensed water loading (mg/m2)
    call outfld('DAERH2O ',daerh2o ,pcols,lchnk)
    !       number concentrations (1/cm3)
    call outfld('NNAT_0  ',nnat_0 ,pcols,lchnk)
    call outfld('NNAT_1  ',nnat_1 ,pcols,lchnk)
    call outfld('NNAT_2  ',nnat_2 ,pcols,lchnk)
    !=0        call outfld('NNAT_3  ',nnat_3 ,pcols,lchnk)
    call outfld('NNAT_4  ',nnat_4 ,pcols,lchnk)
    call outfld('NNAT_5  ',nnat_5 ,pcols,lchnk)
    call outfld('NNAT_6  ',nnat_6 ,pcols,lchnk)
    call outfld('NNAT_7  ',nnat_7 ,pcols,lchnk)
    call outfld('NNAT_8  ',nnat_8 ,pcols,lchnk)
    call outfld('NNAT_9  ',nnat_9 ,pcols,lchnk)
    call outfld('NNAT_10 ',nnat_10,pcols,lchnk)
    !=0        call outfld('NNAT_11 ',nnat_11,pcols,lchnk)
    call outfld('NNAT_12 ',nnat_12,pcols,lchnk)
    !=0        call outfld('NNAT_13 ',nnat_13,pcols,lchnk)
    call outfld('NNAT_14 ',nnat_14,pcols,lchnk)
    !akc6        call outfld('AIRMASSL',airmassl,pcols,lchnk)
    call outfld('AIRMASSL',airmassl,pcols,lchnk)
    call outfld('AIRMASS ',airmass,pcols,lchnk)  !akc6

    !c_er3d 
    !       effective dry radii (um) in each layer
    !        call outfld('ERLT053D',erlt053d,pcols,lchnk)
    !        call outfld('ERGT053D',ergt053d,pcols,lchnk)
    !        call outfld('ER3D    ',er3d    ,pcols,lchnk)
    !c_er3d           
    !       column integrated effective dry radii (um)
    call outfld('DERLT05 ',derlt05,pcols,lchnk)
    call outfld('DERGT05 ',dergt05,pcols,lchnk)
    call outfld('DER     ',der    ,pcols,lchnk)
    !
    !     Extra AeroCom diagnostics requiring table look-ups with RH = constant 

#ifdef AEROCOM_INSITU
    irfmax=6
#else
    irfmax=1
#endif  ! AEROCOM_INSITU

    !     Note: using xrhnull etc as proxy for constant RH input values (see oslo_aero_sw_tables.F90)
    do irf=1,irfmax
       do k=1,pver
          do icol=1,ncol
             xrhnull(icol,k)=xrhrf(irf)
             irh1null(icol,k)=irhrf1(irf)
          end do
       enddo
       call opticsAtConstRh(lchnk, ncol, pint, rhoda, Nnatk, xrhnull, irh1null, irf, &
            xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
            xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
            xfombg, ifombg1, vnbcarr, vaitbcarr, v_soana)
    end do ! irf

  end subroutine aerocom

  subroutine opticsAtConstRh (lchnk, ncol, pint, rhoda, Nnatk, xrh, irh1, irf, &
       xct, ict1, xfaq, ifaq1, xfbcbg, ifbcbg1,           &
       xfbcbgn, ifbcbgn1, xfac, ifac1, xfbc, ifbc1,       &
       xfombg, ifombg1, vnbc, vaitbc, v_soana)

    ! Extra AeroCom diagnostics requiring table look-ups with constant/fixed RH,
    ! i.e. for RH = (/"00","40","55","65","75","85" /) (see oslo_aero_sw_tables.F90)

    ! Input arguments
    !
    integer,  intent(in) :: lchnk                      ! chunk identifier
    integer,  intent(in) :: ncol                       ! number of atmospheric columns
    real(r8), intent(in) :: pint(pcols,pverp)          ! Model interface pressures (10*Pa)
    real(r8), intent(in) :: rhoda(pcols,pver)          ! Density of dry air (kg/m^3)
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration  
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    integer,  intent(in) :: irf
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)        
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! faqm for use in the interpolations 
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbcbg(pcols,pver)
    integer,  intent(in) :: ifbcbg1(pcols,pver)
    real(r8), intent(in) :: xfbcbgn(pcols,pver)
    integer,  intent(in) :: ifbcbgn1(pcols,pver)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! facm for use in the interpolations 
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)   ! fbcm for use in the interpolations 
    integer,  intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfombg(pcols,pver)
    integer,  intent(in) :: ifombg1(pcols,pver)
    real(r8), intent(in) :: vnbc(pcols,pver)
    real(r8), intent(in) :: vaitbc(pcols,pver)
    real(r8), intent(in) :: v_soana(pcols,pver)
    !
    ! Local variables
    !
    integer  :: i, k, icol, mplus10, irh
    integer  :: iloop
    real(r8) :: deltah
    real(r8) :: dod550rh(pcols), abs550rh(pcols)
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

    ! Additionl AeroCom Phase III output:   
    real(r8) :: ec440rh_aer(pcols,pver)
    real(r8) :: abs440rh_aer(pcols,pver)
    real(r8) :: ec870rh_aer(pcols,pver)
    real(r8) :: abs870rh_aer(pcols,pver)
    real(r8) :: be550lt1_aer(pcols,pver,0:nbmodes)
    real(r8) :: ec550rhlt1_aer(pcols,pver)
    real(r8) :: abs550rh_bc(pcols,pver)
    real(r8) :: abs550rh_oc(pcols,pver)
    real(r8) :: abs550rh_su(pcols,pver)
    real(r8) :: abs550rh_ss(pcols,pver)
    real(r8) :: abs550rh_du(pcols,pver)
    real(r8) :: ec550rhlt1_bc(pcols,pver)
    real(r8) :: ec550rhlt1_oc(pcols,pver)
    real(r8) :: ec550rhlt1_su(pcols,pver)
    real(r8) :: ec550rhlt1_ss(pcols,pver)
    real(r8) :: ec550rhlt1_du(pcols,pver) 
    real(r8) :: bedustlt1(pcols,pver) 
    real(r8) :: bedustgt1(pcols,pver)
    real(r8) :: besslt1(pcols,pver)
    real(r8) :: bessgt1(pcols,pver)
    real(r8) :: bbclt1xt(pcols,pver)
    real(r8) :: boclt1xt(pcols,pver)
    real(r8) :: bocgt1xt(pcols,pver)

    character(len=10) :: modeString
    character(len=20) :: varname
    !--------------------------------------------------

    belt1x(:,:,:) = 0._r8

    do iloop=1,1

       ! BC(ax) mode (hydrophobic, so no rhum needed here):
       call extinction_coeffs%intaeropt0(lchnk, ncol, Nnatk)

       ! SO4(Ait), BC(Ait) and OC(Ait) modes:
       mplus10=0
       call extinction_coeffs%intaeropt1(lchnk, ncol, xrh, irh1, mplus10,  &
            Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1)

       mplus10=0
       call extinction_coeffs%intaeropt2to3(lchnk, ncol, xrh, irh1, mplus10, &
            Nnatk, xct, ict1, xfac, ifac1)

       ! BC&OC(Ait) (4), OC&BC(Ait) mode
       mplus10=0
       call extinction_coeffs%intaeropt4(lchnk, ncol, xrh, irh1, mplus10, Nnatk,  &
            xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, xfaq, ifaq1)

       ! SO4(Ait75) (5), Mineral (6-7) and Sea-salt (8-10) modes:
       call extinction_coeffs%intaeropt5to10(lchnk, ncol, xrh, irh1, Nnatk,   &
            xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1)

       ! then to the externally mixed SO4(n), BC(n) and OC(n) modes:
       mplus10=1
       call extinction_coeffsn%intaeropt2to3(lchnk, ncol, xrh, irh1, mplus10,  &
            Nnatk, xct, ict1, xfac, ifac1)

       ! and finally the BC&OC(n) mode:
       mplus10=1
       call extinction_coeffsn%intaeropt4(lchnk, ncol, xrh, irh1, mplus10, Nnatk,    &
            xfbcbgn, ifbcbgn1, xct, ict1, xfac, ifac1, xfaq, ifaq1)

    end do ! iloop


    ! Initialization
    do k=1,pver  
       do icol=1,ncol
          ec550rh_aer(icol,k)     = 0.0_r8 
          abs550rh_aer(icol,k)    = 0.0_r8 
          ec550rhlt1_aer(icol,k)  = 0.0_r8 
          abs550rh_bc(icol,k)     = 0.0_r8
          abs550rh_oc(icol,k)     = 0.0_r8
          abs550rh_su(icol,k)     = 0.0_r8
          abs550rh_ss(icol,k)     = 0.0_r8
          abs550rh_du(icol,k)     = 0.0_r8
          ec440rh_aer(icol,k)     = 0.0_r8 
          abs440rh_aer(icol,k)    = 0.0_r8 
          ec870rh_aer(icol,k)     = 0.0_r8 
          abs870rh_aer(icol,k)    = 0.0_r8 
          basu550tot(icol,k)      = 0.0_r8 
          babc550tot(icol,k)      = 0.0_r8 
          baoc550tot(icol,k)      = 0.0_r8 
          bebglt1t(icol,k)        = 0.0_r8
          bebclt1t(icol,k)        = 0.0_r8
          beoclt1t(icol,k)        = 0.0_r8
          bes4lt1t(icol,k)        = 0.0_r8
          bedustlt1(icol,k)       = 0.0_r8
          besslt1(icol,k)         = 0.0_r8
       end do
    end do
    do icol=1,ncol
       dod550rh(icol) = 0.0_r8 
       abs550rh(icol) = 0.0_r8 
    end do

    ! Calculation of extinction at given RH and absorption for all r and for r<0.5um
    do k=1,pver
       do icol=1,ncol

          do i=0,10
             ec550rh_aer(icol,k)  = ec550rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext550(icol,k,i)
             abs550rh_aer(icol,k) = abs550rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs550(icol,k,i)
             ec440rh_aer(icol,k)  = ec440rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext440(icol,k,i)
             abs440rh_aer(icol,k) = abs440rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs440(icol,k,i)
             ec870rh_aer(icol,k)  = ec870rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffs%bext870(icol,k,i)
             abs870rh_aer(icol,k) = abs870rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%babs870(icol,k,i)
             basu550tot(icol,k)   = basu550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%basu550(icol,k,i)
             babc550tot(icol,k)   = babc550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%babc550(icol,k,i)
             baoc550tot(icol,k)   = baoc550tot(icol,k)   + Nnatk(icol,k,i)*extinction_coeffs%baoc550(icol,k,i)
             bes4lt1t(icol,k)     = bes4lt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%besu550lt1(icol,k,i)
             bebclt1t(icol,k)     = bebclt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%bebc550lt1(icol,k,i)
             beoclt1t(icol,k)     = beoclt1t(icol,k)     + Nnatk(icol,k,i)*extinction_coeffs%beoc550lt1(icol,k,i)
          enddo
          do i=11,14
             ec550rh_aer(icol,k)  = ec550rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext550(icol,k,i-10)
             abs550rh_aer(icol,k) = abs550rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs550(icol,k,i-10)
             ec440rh_aer(icol,k)  = ec440rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext440(icol,k,i-10)
             abs440rh_aer(icol,k) = abs440rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs440(icol,k,i-10)
             ec870rh_aer(icol,k)  = ec870rh_aer(icol,k)  + Nnatk(icol,k,i)*extinction_coeffsn%bext870(icol,k,i-10)
             abs870rh_aer(icol,k) = abs870rh_aer(icol,k) + Nnatk(icol,k,i)*extinction_coeffsn%babs870(icol,k,i-10)
             ba550x(icol,k,i)     = extinction_coeffsn%babs550(icol,k,i-10)
             belt1x(icol,k,i)     = extinction_coeffs%bebg550lt1(icol,k,i-10) !???
          enddo
          do i=6,7
             bedustlt1(icol,k) = bedustlt1(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%bebg550lt1(icol,k,i)
          enddo
          do i=8,10
             besslt1(icol,k) = besslt1(icol,k) + Nnatk(icol,k,i)*extinction_coeffs%bebg550lt1(icol,k,i)
          enddo
          ec550rhlt1_du(icol,k) = bedustlt1(icol,k)
          ec550rhlt1_ss(icol,k) = besslt1(icol,k)

          !soa: *(1-v_soan) for the sulfate volume fraction of mode 11
          bbclt1xt(icol,k) = Nnatk(icol,k,12)*belt1x(icol,k,12) &
               + Nnatk(icol,k,14)*belt1x(icol,k,14)*vnbc(icol,k)
          !soa + v_soan part of mode 11 for the OC volume fraction of that mode
          boclt1xt(icol,k) = Nnatk(icol,k,13)*belt1x(icol,k,13) &
               + Nnatk(icol,k,14)*belt1x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          ec550rhlt1_su(icol,k) = bes4lt1t(icol,k)                         &  ! condensate
               + Nnatk(icol,k,1)*extinction_coeffs%bebg550lt1(icol,k,1)*(1.0_r8-v_soana(icol,k))&  ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%bebg550lt1(icol,k,5)                            ! background, SO4(Ait75) mode (5)
          ec550rhlt1_bc(icol,k) = bebclt1t(icol,k)+bbclt1xt(icol,k)        &  ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%bebg550lt1(icol,k,2)                        &  ! background, BC(Ait) mode (2)
               + Nnatk(icol,k,4)*extinction_coeffs%bebg550lt1(icol,k,4)*vaitbc(icol,k)         &  ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%bebg550lt1(icol,k,0)                           ! background, BC(ax) mode (0)

          !soa + v_soan part of mode 11 for the OC volume fraction of that mode
          ec550rhlt1_oc(icol,k) = beoclt1t(icol,k)+boclt1xt(icol,k)        &  ! coagulated + n-mode OC (13)
               + Nnatk(icol,k,3)*extinction_coeffs%bebg550lt1(icol,k,3)                        &  ! background, OC(Ait) mode (3)
               + Nnatk(icol,k,4)*extinction_coeffs%bebg550lt1(icol,k,4)*(1.0_r8-vaitbc(icol,k))&  ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,1)*extinction_coeffs%bebg550lt1(icol,k,1)*v_soana(icol,k)

          ec550rhlt1_aer(icol,k) = ec550rhlt1_su(icol,k)+ec550rhlt1_bc(icol,k) &
               + ec550rhlt1_oc(icol,k) + ec550rhlt1_ss(icol,k)+ec550rhlt1_du(icol,k)
          ec550rhlt1_aer(icol,k) = 1.e-3_r8*ec550rhlt1_aer(icol,k)

          abs550rh_du(icol,k) = Nnatk(icol,k,6)*extinction_coeffs%babg550(icol,k,6) &
               + Nnatk(icol,k,7)*extinction_coeffs%babg550(icol,k,7)
          abs550rh_ss(icol,k) = Nnatk(icol,k,8)*extinction_coeffs%babg550(icol,k,8) &
               + Nnatk(icol,k,9)*extinction_coeffs%babg550(icol,k,9) &
               + Nnatk(icol,k,10)*extinction_coeffs%babg550(icol,k,10)

          !soa: *(1-v_soana) for the sulfate volume fraction of mode 1
          abs550rh_su(icol,k) = basu550tot(icol,k)                   &  ! condensate:w
               + (1.0_r8-v_soana(icol,k))*Nnatk(icol,k,1)*extinction_coeffs%babg550(icol,k,1) &  ! background, SO4(Ait) mode (1)
               + Nnatk(icol,k,5)*extinction_coeffs%babg550(icol,k,5)    ! background, SO4(Ait75) mode (5)

          !soa: *(1-v_soan) for the sulfate volume fraction
          babc550xt(icol,k) = Nnatk(icol,k,12)*ba550x(icol,k,12)  &
               + Nnatk(icol,k,14)*ba550x(icol,k,14)*vnbc(icol,k)

          baoc550xt(icol,k) = Nnatk(icol,k,13)*ba550x(icol,k,13) &
               + Nnatk(icol,k,14)*ba550x(icol,k,14)*(1.0_r8-vnbc(icol,k)) 

          abs550rh_bc(icol,k) = babc550tot(icol,k)+babc550xt(icol,k) &                       ! coagulated + n-mode BC (12)
               + Nnatk(icol,k,2)*extinction_coeffs%babg550(icol,k,2) &                       ! background, BC(Ait) mode (2)
               + vaitbc(icol,k)*Nnatk(icol,k,4)*extinction_coeffs%babg550(icol,k,4) &        ! background in OC&BC(Ait) mode (4)
               + Nnatk(icol,k,0)*extinction_coeffs%babg550(icol,k,0)                         ! background, BC(ax) mode (0)

          abs550rh_oc(icol,k) = baoc550tot(icol,k)+baoc550xt(icol,k) &                       ! coagulated + n-mode OC (13)
               + v_soana(icol,k)*Nnatk(icol,k,1)*extinction_coeffs%babg550(icol,k,1) &       ! SOA fraction of mode 1
               + Nnatk(icol,k,3)*extinction_coeffs%babg550(icol,k,3) &                       ! background, OC(Ait) mode (3)
               + (1.0_r8-vaitbc(icol,k))*Nnatk(icol,k,4)*extinction_coeffs%babg550(icol,k,4) ! background in OC&BC(Ait) mode (4)

          deltah=1.e-4_r8*(pint(icol,k+1)-pint(icol,k))/(rhoda(icol,k)*9.8_r8)
          dod550rh(icol) = dod550rh(icol)+ec550rh_aer(icol,k)*deltah
          abs550rh(icol) = abs550rh(icol)+abs550rh_aer(icol,k)*deltah

          ec550rh_aer(icol,k)  = 1.e-3_r8*ec550rh_aer(icol,k)
          abs550rh_aer(icol,k) = 1.e-3_r8*abs550rh_aer(icol,k)
          ec440rh_aer(icol,k)  = 1.e-3_r8*ec440rh_aer(icol,k)
          abs440rh_aer(icol,k) = 1.e-3_r8*abs440rh_aer(icol,k)
          ec870rh_aer(icol,k)  = 1.e-3_r8*ec870rh_aer(icol,k)
          abs870rh_aer(icol,k) = 1.e-3_r8*abs870rh_aer(icol,k)

          abs550rh_bc(icol,k)  = 1.e-3_r8*abs550rh_bc(icol,k)
          abs550rh_oc(icol,k)  = 1.e-3_r8*abs550rh_oc(icol,k)
          abs550rh_su(icol,k)  = 1.e-3_r8*abs550rh_su(icol,k)
          abs550rh_ss(icol,k)  = 1.e-3_r8*abs550rh_ss(icol,k)
          abs550rh_du(icol,k)  = 1.e-3_r8*abs550rh_du(icol,k)

       enddo
    enddo

    if(irf.eq.1) then

       call outfld('ECDRYAER',ec550rh_aer,pcols,lchnk)
       call outfld('ABSDRYAE',abs550rh_aer,pcols,lchnk)
       call outfld('OD550DRY',dod550rh,pcols,lchnk)       ! 2D variable
       call outfld('AB550DRY',abs550rh,pcols,lchnk)       ! 2D variable
       call outfld('ECDRY440',ec440rh_aer,pcols,lchnk)
       call outfld('ABSDR440',abs440rh_aer,pcols,lchnk)
       call outfld('ECDRY870',ec870rh_aer,pcols,lchnk)
       call outfld('ABSDR870',abs870rh_aer,pcols,lchnk)
       call outfld('ECDRYLT1',ec550rhlt1_aer,pcols,lchnk)
       !         Since we do not have enough look-up table info to take abs550rhlt1_aer,
       !         instead take out abs550rh for each constituent:
       call outfld('ABSDRYBC',abs550rh_bc,pcols,lchnk)
       call outfld('ABSDRYOC',abs550rh_oc,pcols,lchnk)
       call outfld('ABSDRYSU',abs550rh_su,pcols,lchnk)
       call outfld('ABSDRYSS',abs550rh_ss,pcols,lchnk)
       call outfld('ABSDRYDU',abs550rh_du,pcols,lchnk)

    elseif(irf.ge.2) then   ! only happens for AEROCOM_INSITU

       irh=RF(irf)

       modeString="  "
       write(modeString,"(I2)"),irh
       if(RF(irf).eq.0) modeString="00"
       varName = "EC55RH"//trim(modeString)
       call outfld(varName,ec550rh_aer(:,:),pcols,lchnk)
       varName = "AB55RH"//trim(modeString)
       call outfld(varName,abs550rh_aer(:,:),pcols,lchnk)

    end if ! irf

  end subroutine opticsAtConstRh

  subroutine intfrh (lchnk, ncol, v3so4, v3insol, v3oc, v3ss, relh, frh)

    ! Written by Alf Kirkevaag in November 2011, based on interpol1to3 in optinterpol.F90
    ! called by NorESM/physpkg

    ! Input arguments
    integer, intent(in)  :: lchnk                      ! chunk identifier
    integer, intent(in)  :: ncol                       ! number of atmospheric columns
    real(r8), intent(in) :: v3so4(pcols,pver,nmodes)   ! Modal mass fraction of Sulfate
    real(r8), intent(in) :: v3insol(pcols,pver,nmodes) ! Modal mass fraction of BC and dust
    real(r8), intent(in) :: v3oc(pcols,pver,nmodes)    ! Modal mass fraction of OC (POM)
    real(r8), intent(in) :: v3ss(pcols,pver,nmodes)    ! Modal mass fraction of sea-salt
    real(r8), intent(in) :: relh(pcols,pver)           ! Ambient relatve humidity (fraction)
    !
    ! Output arguments
    real(r8), intent(out) :: frh(pcols,pver,nmodes)    ! Modal humidity growth factor 
    !
    ! Local variables
    integer  :: i, ierr, irelh, kcomp, k, icol
    integer  :: irh1(pcols,pver), irh2(pcols,pver)
    real(r8) :: a, b, e, fso4, finsol, foc, fss
    real(r8) :: xrh(pcols,pver)
    integer  :: t_irh1, t_irh2
    real(r8) :: t_xrh, t_rh1, t_rh2
    parameter (e=2.718281828)

    ! Relative humidity intries from oslo_aero_sw_tables
    ! rh = (/ 0.0_r8,  0.37_r8, 0.47_r8, 0.65_r8, 0.75_r8, &
    !         0.8_r8,  0.85_r8, 0.9_r8,  0.95_r8, 0.995_r8 /)
    ! Humidity growth factors which are consistent with the aerosol optics look-up tables:
    real(r8), dimension(10) :: fh_SO4   = &
         (/ 1.00_r8, 1.34_r8, 1.40_r8, 1.53_r8, 1.64_r8, &
         1.71_r8, 1.81_r8, 1.98_r8, 2.39_r8, 5.04_r8 /)
    real(r8), dimension(10) :: fh_insol = &
         (/ 1.00_r8, 1.01_r8, 1.01_r8, 1.02_r8, 1.02_r8, &
         1.02_r8, 1.02_r8, 1.02_r8, 1.02_r8, 1.02_r8 /)
    real(r8), dimension(10) :: fh_OC    = &
         (/ 1.00_r8, 1.02_r8, 1.05_r8, 1.14_r8, 1.19_r8, &
         1.22_r8, 1.27_r8, 1.36_r8, 1.59_r8, 3.18_r8 /)
    real(r8), dimension(10) :: fh_SS    = &
         (/ 1.00_r8, 1.01_r8, 1.02_r8, 1.56_r8, 1.87_r8, &
         1.97_r8, 2.12_r8, 2.35_r8, 2.88_r8, 6.08_r8 /)
    ! -----------------------------------------

    ! write(*,*) 'Before xrh-loop'
    do k=1,pver
       do icol=1,ncol
          !test          xrh(icol,k)  = 0.8
          xrh(icol,k)  = min(max(relh(icol,k),rh(1)),rh(10))
       end do
    end do

    ! write(*,*) 'Before rh-loop'
    do irelh=1,9
       do k=1,pver
          do icol=1,ncol
             if(xrh(icol,k) >= rh(irelh).and. &
                  xrh(icol,k)<=rh(irelh+1)) then
                irh1(icol,k)=irelh
                irh2(icol,k)=irelh+1
             endif
          end do
       end do
    end do

    ! Loop over all relevant modes (kcomp=1,2,4-11,13,14)
    ! (mode 3 is no longer included, and 12 is insoluble)

    do kcomp=1,14

       do icol=1,ncol
          do k=1,pver
             frh(icol,k,kcomp)=0.0_r8
          end do
       end do

       if(kcomp.ne.3.and.kcomp.ne.12) then

          do k=1,pver
             do icol=1,ncol

                ! Collect all the vector elements into temporary storage
                ! to avoid cache conflicts and excessive cross-referencing

                t_irh1 = irh1(icol,k)
                t_irh2 = irh2(icol,k)

                ! write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2

                t_rh1  = rh(t_irh1)
                t_rh2  = rh(t_irh2)

                t_xrh  = xrh(icol,k)

                if(t_xrh <= 0.37) then  ! linear averaging w.r.t. small RH:
                   fso4  = ((t_rh2-t_xrh)*fh_SO4(t_irh1)+(t_xrh-t_rh1)*fh_SO4(t_irh2)) /(t_rh2-t_rh1)    
                   finsol= ((t_rh2-t_xrh)*fh_insol(t_irh1)+(t_xrh-t_rh1)*fh_insol(t_irh2)) /(t_rh2-t_rh1)    
                   foc   = ((t_rh2-t_xrh)*fh_OC(t_irh1)+(t_xrh-t_rh1)*fh_OC(t_irh2)) /(t_rh2-t_rh1)    
                   fss   = ((t_rh2-t_xrh)*fh_SS(t_irh1)+(t_xrh-t_rh1)*fh_SS(t_irh2)) /(t_rh2-t_rh1)    
                else                  ! exponential averaging w.r.t. large RH:   
                   a = (log(fh_SO4(t_irh2))-log(fh_SO4(t_irh1)))/(t_rh2-t_rh1)
                   b = (t_rh2*log(fh_SO4(t_irh1))-t_rh1*log(fh_SO4(t_irh2)))/(t_rh2-t_rh1)
                   fso4 = e**(a*t_xrh+b)
                   a = (log(fh_insol(t_irh2))-log(fh_insol(t_irh1)))/(t_rh2-t_rh1)
                   b = (t_rh2*log(fh_insol(t_irh1))-t_rh1*log(fh_insol(t_irh2)))/(t_rh2-t_rh1)
                   finsol = e**(a*t_xrh+b)
                   a = (log(fh_OC(t_irh2))-log(fh_OC(t_irh1)))/(t_rh2-t_rh1)
                   b = (t_rh2*log(fh_OC(t_irh1))-t_rh1*log(fh_OC(t_irh2)))/(t_rh2-t_rh1)
                   foc = e**(a*t_xrh+b)
                   a = (log(fh_SS(t_irh2))-log(fh_SS(t_irh1)))/(t_rh2-t_rh1)
                   b = (t_rh2*log(fh_SS(t_irh1))-t_rh1*log(fh_SS(t_irh2)))/(t_rh2-t_rh1)
                   fss = e**(a*t_xrh+b)
                endif

                ! linear interpolation w.r.t. mass fractions of each internally mixed component
                ! (this assumption is only used here, while the full Koehler equation are solved
                ! for the look-up tables for log-normal size distributions and aerosol optics):

                frh(icol,k,kcomp) = v3so4(icol,k,kcomp)*fso4 + v3insol(icol,k,kcomp)*finsol &
                     + v3oc(icol,k,kcomp) *foc + v3ss(icol,k,kcomp)*fss

                !  write(*,*) 'frh =', frh(icol,k,kcomp)
             end do ! icol
          end do   ! k
       endif    ! kcomp.ne.3.and.kcomp.ne.12
    end do   ! kcomp

  end subroutine intfrh

#endif

end module oslo_aero_aerocom
