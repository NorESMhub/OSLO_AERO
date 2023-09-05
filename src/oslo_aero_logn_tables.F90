module oslo_aero_logn_tables

  use shr_kind_mod,            only: r8 => shr_kind_r8
  use ppgrid,                  only: pcols
  use cam_logfile,             only: iulog
  use spmd_utils,              only: masterproc
  !
  use oslo_aero_control,       only: oslo_aero_getopts,dir_string_length
  use oslo_aero_sw_tables,     only: cate, fac, faq, fbc, cat
  use oslo_aero_linear_interp, only: lininterpol3dim, lininterpol4dim
  use oslo_aero_params,        only: nmodes, nbmodes
  use oslo_aero_share

  implicit none
  private

  public :: initlogn
  public :: intlog1to3_sub
  public :: intlog4_sub
  public :: intlog5to10_sub

  real(r8) :: rrr1to3 (3,16,6)         ! Modal radius array, mode 1 - 3
  real(r8) :: sss1to3 (3,16,6)         ! Standard deviation array, Mode 1 -3
  real(r8) :: rrr4 (16,6,6)            ! Modal radius array, mode 4
  real(r8) :: sss4 (16,6,6)            ! Modal radius array, mode 4
  real(r8) :: rrr (5:10,6,6,6,6)       ! Modal radius array, mode 5 - 10
  real(r8) :: sss (5:10,6,6,6,6)       ! Standard deviation array, mode 5 - 10

  real(r8) :: calog1to3(3,96)          ! Array for reading catot from file
  real(r8) :: rk1to3 (3,96)            ! Array for reading modal radius from file
  real(r8) :: stdv1to3 (3,96)          ! Array for reading std. dev. from file
  real(r8) :: fraclog1to3 (3,96)       ! Same as frac4, but for initlogn.F90

  real(r8) :: calog4(576)              ! Same as catot4, but for initlogn.F90
  real(r8) :: fraclog4(576)            ! Same as frac4, but for initlogn.F90
  real(r8) :: fraqlog4(576)            ! Same as fraq4, but for initlogn.F90
  real(r8) :: rk4 (576)                ! Array for reading modal radius from file
  real(r8) :: stdv4 (576)              ! Array for reading std. dev. from file

  real(r8) :: calog (5:10,1296)        ! Same as catot, but for initlogn.F90
  real(r8) :: fraclog5to10 (5:10,1296) ! Same as frac5to10, but for initlogn.F90
  real(r8) :: fabclog5to10 (5:10,1296) ! Same as fabc5to10, but for initlogn.F90
  real(r8) :: fraqlog5to10 (5:10,1296) ! Same as fraq5to10, but for initlogn.F90
  real(r8) :: rk5to10 (5:10,1296)      ! Array for reading modal radius from file
  real(r8) :: stdv5to10 (5:10,1296)    ! Array for reading std. dev. from file

!=======================================================
contains
!=======================================================

  subroutine initlogn()

    ! Reads the tabulated parameters for "best lognormal fits" of the
    ! aerosol size distribution wrt CCN activation as calculated by Alf Kirkevaag.

    integer  :: kcomp, ictot, ifac, ifbc, ifaq
    integer  :: ic, ifil, lin
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps4 = 1.e-4_r8
    character(len=dir_string_length) :: aerotab_table_dir

    ! Where are the tables stored??
    call oslo_aero_getopts(aerotab_table_dir_out=aerotab_table_dir)

    open(20,file=trim(aerotab_table_dir)//'/logntilp1.out' ,form='formatted',status='old')  ! SO4&SOA(n/Ait)
    open(21,file=trim(aerotab_table_dir)//'/logntilp2.out' ,form='formatted',status='old')  ! BC(n/Ait)
    open(22,file=trim(aerotab_table_dir)//'/logntilp3.out' ,form='formatted',status='old')  ! OC(n/Ait)
    open(23,file=trim(aerotab_table_dir)//'/logntilp4.out' ,form='formatted',status='old')  ! BC&OC(n/Ait)
    open(24,file=trim(aerotab_table_dir)//'/logntilp5.out' ,form='formatted',status='old')  ! SO4(Ait75)
    open(25,file=trim(aerotab_table_dir)//'/logntilp6.out' ,form='formatted',status='old')  ! MINACC
    open(26,file=trim(aerotab_table_dir)//'/logntilp7.out' ,form='formatted',status='old')  ! MINCOA
    open(27,file=trim(aerotab_table_dir)//'/logntilp8.out' ,form='formatted',status='old')  ! SEASF
    open(28,file=trim(aerotab_table_dir)//'/logntilp9.out' ,form='formatted',status='old')  ! SEASACC
    open(29,file=trim(aerotab_table_dir)//'/logntilp10.out',form='formatted',status='old')  ! SEASCOA
    if (masterproc) then
       write(iulog,*)'nlog open ok'
    end if

    ! Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 20,29
       call checkTableHeader (ifil)
    enddo

    ! ************************************************************************
    ! Mode  1      (SO4&SOA + condesate from H2SO4 and SOA)
    ! Modes 2 to 3 (BC/OC + condesate from H2SO4 and SOA)
    !
    ! These two are treated the same way since there is no dependence on
    ! fombg (SOA fraction in the background) for mode 1.
    ! ************************************************************************

    do ifil = 1,2
       do lin = 1,96   ! 16*6 entries
          read(19+ifil,993) kcomp, calog1to3(ifil,lin), fraclog1to3 (ifil, lin), &
               rk1to3(ifil,lin), stdv1to3(ifil,lin)

          do ic=1,16
             if(abs((calog1to3(ifil,lin)-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
                ictot=ic
                exit
             endif
          end do
          do ic=1,6
             if(abs(fraclog1to3(ifil,lin)-fac(ic))<eps4) then
                ifac=ic
                exit
             endif
          end do
          sss1to3(kcomp,ictot,ifac) = stdv1to3(ifil,lin)
          rrr1to3(kcomp,ictot,ifac) = rk1to3(ifil,lin)
       end do   ! lin
    end do    ! ifil

    ! Prescribed dummy values for kcomp=3
    kcomp=3
    do ictot=1,16
       do ifac=1,6
          sss1to3(kcomp,ictot,ifac)=1.0_r8
          rrr1to3(kcomp,ictot,ifac)=1.0_r8
       enddo
    enddo

    do kcomp=1,2
       do ictot=1,16
          do ifac=1,6
             if(sss1to3(kcomp,ictot,ifac)<=0.0_r8) then
                write(*,*) 'sss1to3 =',  ictot, ifac, sss1to3(kcomp,ictot,ifac)
                write(*,*) 'Error in initialization of sss1to3'
                stop
             endif
             if(rrr1to3(kcomp,ictot,ifac)<=0.0_r8) then
                write(*,*) 'rrr1to3 =', ictot, ifac, rrr1to3(kcomp,ictot,ifac)
                write(*,*) 'Error in initialization of rrr1to3'
                stop
             endif
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'nlog mode 1-3 ok'
    end if

    ! ************************************************************************
    ! Mode 4 (BC&OC + condesate from H2SO4 + wetphase (NH4)2SO4)
    ! ************************************************************************

    ifil = 4
    do lin = 1,576   ! 16 entries
       read(19+ifil,994) kcomp, calog4(lin) &
            ,fraclog4(lin), fraqlog4(lin), rk4(lin), stdv4(lin)

       do ic=1,16
          if(abs((calog4(lin)-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(fraclog4(lin)-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(fraqlog4(lin)-faq(ic))<eps4) then
             ifaq=ic
             exit
          endif
       end do

       rrr4(ictot,ifac,ifaq) = rk4(lin)
       sss4(ictot,ifac,ifaq) = stdv4(lin)
    end do   ! lin

    do ifac=1,6
       do ifaq=1,6
          do ictot=1,16
             if(rrr4(ictot,ifac,ifaq)<=0.0_r8) then
                write(*,*) 'rrr4 =',ictot,ifac,ifaq,rrr4(ictot,ifac,ifaq)
                write(*,*) 'Error in initialization of rrr4'
                stop
             endif
             if(sss4(ictot,ifac,ifaq)<=0.0_r8) then
                write(*,*) 'sss4 =',ictot,ifac,ifaq,sss4(ictot,ifac,ifaq)
                write(*,*) 'Error in initialization of sss4'
                stop
             endif
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'nlog mode 4 ok'
    end if

    ! ************************************************************************
    ! Modes 5 to 10 (SO4(ait75) and mineral and seasalt-modes + cond./coag./aq.)
    ! ************************************************************************

    do ifil = 5,10
       do lin = 1,1296   ! 6**4 entries
          read(19+ifil,995) kcomp, calog(ifil,lin), &
               fraclog5to10(ifil,lin), fabclog5to10(ifil,lin), fraqlog5to10(ifil,lin), &
               rk5to10(ifil,lin), stdv5to10(ifil,lin)

          do ic=1,6
             if(abs((calog(ifil,lin)-cat(kcomp,ic))/cat(kcomp,ic))<eps2) then
                ictot=ic
                exit
             endif
          end do
          do ic=1,6
             if(abs(fraclog5to10(ifil,lin)-fac(ic))<eps4) then
                ifac=ic
                exit
             endif
          end do
          ! PRINT*,'ifac not found',fraclog5to10(ifil,lin) TODO:
          do ic=1,6
             if(abs((fabclog5to10(ifil,lin)-fbc(ic))/fbc(ic))<eps2) then
                ifbc=ic
                exit
             endif
          end do
          ! PRINT*,'ifbc not found',fabclog5to10(ifil,lin) TODO:
          do ic=1,6
             if(abs(fraqlog5to10(ifil,lin)-faq(ic))<eps4) then
                ifaq=ic
                exit
             endif
          end do
          ! PRINT*,'ifaq not found',fraqlog5to10(ifil,lin) TODO:

          rrr(kcomp,ictot,ifac,ifbc,ifaq) = rk5to10(ifil,lin)
          sss(kcomp,ictot,ifac,ifbc,ifaq) = stdv5to10(ifil,lin)

       end do   ! lin
    end do    ! ifil

    do kcomp=5,10
       do ifac=1,6
          do ifbc=1,6
             do ictot=1,6
                if(rrr(kcomp,ictot,ifac,ifbc,ifaq)<=0.0_r8) then
                   write(*,*) 'rrr =',kcomp,ictot,ifac,ifbc,ifaq,rrr(kcomp,ictot,ifac,ifbc,ifaq)
                   write(*,*) 'Error in initialization of rrr'
                   stop
                endif
                if(sss(kcomp,ictot,ifac,ifbc,ifaq)<=0.0_r8) then
                   write(*,*) 'sss =',ictot,ifac,ifbc,ifaq,sss(kcomp,ictot,ifac,ifbc,ifaq)
                   write(*,*) 'Error in initialization of sss'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'nlog mode 5-10 ok'

    do ifil=20,29
       close (ifil)
    end do

993 format(I3,4(x,e12.5))
994 format(I3,5(x,e12.5))
995 format(I3,6(x,e12.5))

  end subroutine initlogn

  !=======================================================
  subroutine intlog1to3_sub (ncol, kcomp, xctin, Nnat, xfacin, cxs, xstdv, xrk)

    ! Created by Trude Storelvmo, fall 2007. This subroutine gives as output
    ! the "new" modal radius and standard deviation for a given aerosol mode, kcomp
    ! 1-3. These parameters are calculated for a best lognormal fit approximation of
    ! the aerosol size distribution. This because the aerosol activation routine
    ! (developed by Abdul-Razzak & Ghan, 2000) requiers the size distribution to be
    ! described by lognormal modes.
    ! Changed by Alf Kirkevåg to take into account condensation of SOA, September 2015,

    integer  , intent(in)  :: ncol
    integer  , intent(in)  :: kcomp
    real(r8) , intent(in)  :: Nnat(pcols)   ! Modal number concentration
    real(r8) , intent(in)  :: xctin(pcols)  ! total internally mixed conc. (ug/m3)
    real(r8) , intent(in)  :: xfacin(pcols) ! SOA/(SOA+H2SO4) for condensated mass
    real(r8) , intent(out) :: xstdv(pcols)  ! log10 of standard deviation for lognormal fit
    real(r8) , intent(out) :: xrk(pcols)    ! Modal radius for lognormal fit
    real(r8) , intent(out) :: cxs(pcols)    ! excess (modal) internally mixed conc.

    ! local variables
    real(r8) :: camdiff
    real(r8) :: xct(pcols)
    real(r8) :: xfac(ncol)
    integer  :: lon, long
    integer  :: i, ictot, ict1, ict2
    real(r8) :: r1, r2, s1, s2
    integer  :: ifac, ifac1, ifac2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xct, t_cat1, t_cat2
    real(r8) :: r11, r12, r21, r22, s11, s12, s21, s22
    real(r8) :: d2mx(2), dxm1(2), invd(2)
    real(r8) :: ess
    real(r8), parameter :: eps= 1.0e-10_r8

    ! Initialize excess mass cxs, wrt. maximum allowed internal mixing
    do lon=1,ncol
       cxs(lon) = 0.0_r8
       xct(lon) = 0.0_r8
       xfac(lon) = 0.0_r8

       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8

       xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))
       xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
       camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

       cxs(lon)  = max(0.0_r8,camdiff)

       ictot=1
       ess = xct(lon)
       do while (ictot.lt.15 .and. (ess.lt.cate(kcomp,ictot) .or. ess.gt.cate(kcomp,ictot+1)))
          ictot=ictot+1
       enddo
       ict1=ictot
       ict2=ictot+1

       ifac=1
       ess = xfac(lon)
       do while (ifac.lt.5 .and. (ess.lt.fac(ifac) .or. ess.gt.fac(ifac+1)))
          ifac=ifac+1
       enddo
       ifac1=ifac
       ifac2=ifac+1

       ! Collect all the vector elements into temporary storage
       ! to avoid cache conflicts and excessive cross-referencing

       t_cat1 = cate(kcomp,ict1)
       t_cat2 = cate(kcomp,ict2)
       t_fac1 = fac(ifac1)
       t_fac2 = fac(ifac2)

       t_xct  = xct(lon)
       t_xfac = xfac(lon)

       ! partial lengths along each dimension (1-2) for interpolation

       d2mx(1) = (t_cat2-t_xct)
       dxm1(1) = (t_xct-t_cat1)
       invd(1) = 1.0_r8/(t_cat2-t_cat1)
       d2mx(2) = (t_fac2-t_xfac)
       dxm1(2) = (t_xfac-t_fac1)
       invd(2) = 1.0_r8/(t_fac2-t_fac1)

       ! interpolated (in 2 dimensions) modal median radius:

       r11=rrr1to3(kcomp,ict1,ifac1)
       r12=rrr1to3(kcomp,ict1,ifac2)
       r21=rrr1to3(kcomp,ict2,ifac1)
       r22=rrr1to3(kcomp,ict2,ifac2)

       r1 =d2mx(2)*r11+dxm1(2)*r12
       r2 =d2mx(2)*r21+dxm1(2)*r22

       xrk(lon) = (d2mx(1)*r1+dxm1(1)*r2)*invd(2)*invd(1)*1.e-6_r8  !Look-up table radii in um

       ! interpolated (in 2 dimensions) modal standard deviation:

       s11=sss1to3(kcomp,ict1,ifac1)
       s12=sss1to3(kcomp,ict1,ifac2)
       s21=sss1to3(kcomp,ict2,ifac1)
       s22=sss1to3(kcomp,ict2,ifac2)

       s1 =d2mx(2)*s11+dxm1(2)*s12
       s2 =d2mx(2)*s21+dxm1(2)*s22

       xstdv(lon) = (d2mx(1)*s1+dxm1(1)*s2)*invd(2)*invd(1)
    end do   ! lon

  end subroutine intlog1to3_sub

  !=======================================================
  subroutine intlog4_sub (ncol, kcomp, xctin, Nnat, xfacin, xfaqin, cxs, xstdv, xrk)

    ! Created by Trude Storelvmo, fall 2007. This subroutine gives as output
    ! the "new" modal radius and standard deviation for aerosol mode kcomp=4.
    ! These parameters are calculated for a best lognormal fit approximation of
    ! the aerosol size distribution. This because the aerosol activation routine
    ! (developed by Abdul-Razzak & Ghan, 2000) requires the size distribution
    ! to be described by lognormal modes.
    ! Changed by Alf Kirkevåg to take into account condensation of SOA, September
    ! 2015, and also rewritten to a more generalized for for interpolations using
    ! common subroutines interpol*dim.

    integer  , intent(in)  :: ncol
    integer  , intent(in)  :: kcomp
    real(r8) , intent(in)  :: Nnat(pcols)   ! Modal number concentration
    real(r8) , intent(in)  :: xctin(pcols)  ! total internally mixed conc. (ug/m3)
    real(r8) , intent(in)  :: xfacin(pcols) ! SOA/(SOA+H2SO4) for condensated mass
    real(r8) , intent(in)  :: xfaqin(pcols) ! = Cso4a2/(Cso4a1+Cso4a2)
    real(r8) , intent(out) :: xstdv(pcols)  ! log10 of standard deviation for lognormal fit
    real(r8) , intent(out) :: xrk(pcols)    ! Modal radius for lognormal fit
    real(r8) , intent(out) :: cxs(pcols)    ! excess (modal) internally mixed conc.

    ! local variables
    real(r8) :: camdiff
    real(r8) :: xct(pcols)
    real(r8) :: xfac(pcols)
    real(r8) :: xfaq(pcols)
    integer  :: lon, long
    integer  :: i, ictot, ifac, ifaq
    integer  :: ict1, ict2, ifac1, ifac2, ifaq1, ifaq2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xct, t_cat1, t_cat2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: r1, r2, s1, s2, tmp, e
    real(r8) :: d2mx(3), dxm1(3), invd(3)
    real(r8) :: sizepar3d(2,2,2)
    real(r8), parameter :: eps=1.0e-60_r8 ! introduced by (or inspired by) Egil Stoeren:

    ! Initialize excess mass cxs, wrt. maximum allowed internal mixing
    do lon=1,ncol
       cxs(lon)  = 0.0_r8
       xct(lon)  = 0.0_r8
       xfac(lon) = 0.0_r8
       xfaq(lon) = 0.0_r8

       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8

       xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cate(kcomp,1)),cate(kcomp,16))
       xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
       xfaq(lon) = min(max(xfaqin(lon),faq(1)),faq(6))

       camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

       cxs(lon)  = max(0.0_r8,camdiff)

       ictot=1
       tmp = xct(lon)
       do while (ictot.lt.15 .and. (tmp.lt.cate(kcomp,ictot) .or. tmp.gt.cate(kcomp,ictot+1)))
          ictot=ictot+1
       enddo
       ict1=ictot
       ict2=ictot+1

       ifac=1
       tmp = xfac(lon)
       do while (ifac.lt.5 .and. (tmp.lt.fac(ifac) .or. tmp.gt.fac(ifac+1)))
          ifac=ifac+1
       enddo
       ifac1=ifac
       ifac2=ifac+1

       ifaq=1
       tmp = xfaq(lon)
       do while (ifaq.lt.5 .and. (tmp.lt.faq(ifaq) .or. tmp.gt.faq(ifaq+1)))
          ifaq=ifaq+1
       enddo
       ifaq1=ifaq
       ifaq2=ifaq+1

       ! Collect all the vector elements into temporary storage
       ! to avoid cache conflicts and excessive cross-referencing
       t_cat1 = cate(kcomp,ict1)
       t_cat2 = cate(kcomp,ict2)
       t_fac1 = fac(ifac1)
       t_fac2 = fac(ifac2)
       t_faq1 = faq(ifaq1)
       t_faq2 = faq(ifaq2)

       t_xct  = xct(lon)
       t_xfac = xfac(lon)
       t_xfaq = xfaq(lon)

       ! partial lengths along each dimension (1-4) for interpolation
       d2mx(1) = (t_cat2-t_xct)
       dxm1(1) = (t_xct-t_cat1)
       invd(1) = 1.0_r8/(t_cat2-t_cat1)
       d2mx(2) = (t_fac2-t_xfac)
       dxm1(2) = (t_xfac-t_fac1)
       invd(2) = 1.0_r8/(t_fac2-t_fac1)
       d2mx(3) = (t_faq2-t_xfaq)
       dxm1(3) = (t_xfaq-t_faq1)
       invd(3) = 1.0_r8/(t_faq2-t_faq1)

       ! Table points as basis for multidimentional linear interpolation,
       ! modal median radius:

       sizepar3d(1,1,1)=rrr4(ict1,ifac1,ifaq1)
       sizepar3d(1,1,2)=rrr4(ict1,ifac1,ifaq2)
       sizepar3d(1,2,1)=rrr4(ict1,ifac2,ifaq1)
       sizepar3d(1,2,2)=rrr4(ict1,ifac2,ifaq2)
       sizepar3d(2,1,1)=rrr4(ict2,ifac1,ifaq1)
       sizepar3d(2,1,2)=rrr4(ict2,ifac1,ifaq2)
       sizepar3d(2,2,1)=rrr4(ict2,ifac2,ifaq1)
       sizepar3d(2,2,2)=rrr4(ict2,ifac2,ifaq2)

       ! interpolation in the faq and fac dimension
       call lininterpol3dim (d2mx, dxm1, invd, sizepar3d, r1, r2)

       ! finally, interpolation in the cate dimension
       xrk(lon)=(d2mx(1)*r1+dxm1(1)*r2)*invd(1)*1.e-6_r8  ! look up table radii in um


       ! Table points as basis for multidimentional linear interpolation,
       ! modal standard deviation:
       sizepar3d(1,1,1)=sss4(ict1,ifac1,ifaq1)
       sizepar3d(1,1,2)=sss4(ict1,ifac1,ifaq2)
       sizepar3d(1,2,1)=sss4(ict1,ifac2,ifaq1)
       sizepar3d(1,2,2)=sss4(ict1,ifac2,ifaq2)
       sizepar3d(2,1,1)=sss4(ict2,ifac1,ifaq1)
       sizepar3d(2,1,2)=sss4(ict2,ifac1,ifaq2)
       sizepar3d(2,2,1)=sss4(ict2,ifac2,ifaq1)
       sizepar3d(2,2,2)=sss4(ict2,ifac2,ifaq2)

       ! interpolation in the faq and fac dimension
       call lininterpol3dim (d2mx, dxm1, invd, sizepar3d, s1, s2)

       ! finally, interpolation in the cate dimension
       xstdv(lon)=(d2mx(1)*s1+dxm1(1)*s2)*invd(1)

    end do   ! lon

  end subroutine intlog4_sub

  !=======================================================
  subroutine intlog5to10_sub (ncol, kcomp, xctin, Nnat, &
       xfacin, xfbcin, xfaqin, cxs, xstdv, xrk)

    !This subroutine gives as output the "new" modal radius and standard deviation
    !for a given aerosol mode, kcomp 1-5. These parameters are calculated for a
    !best lognormal fit approximation of the aerosol size distribution.

    integer  , intent(in)  :: ncol
    integer  , intent(in)  :: kcomp
    real(r8) , intent(in)  :: Nnat(pcols)   ! Modal number concentration
    real(r8) , intent(in)  :: xctin(pcols)  ! total internally mixed conc. (ug/m3)
    real(r8) , intent(in)  :: xfacin(pcols) ! = (Cbc+Coc)/(Cbc+Coc+Cso4)
    real(r8) , intent(in)  :: xfbcin(pcols) ! = Cbc/(Cbc+Coc)
    real(r8) , intent(in)  :: xfaqin(pcols) ! = Cso4a2/(Cso4a1+Cso4a2)
    real(r8) , intent(out) :: xstdv(pcols)  ! log10 of standard deviation of lognormal fit
    real(r8) , intent(out) :: xrk(pcols)    ! Modal radius of lognormal fit
    real(r8) , intent(out) :: cxs(pcols)    ! excess (modal) internally mixed conc.

    ! local variables
    real(r8) :: xctsave, camdiff
    real(r8) :: xct(pcols), xfac(pcols),  xfbc(pcols), xfaq(pcols)
    integer  :: lon, long
    integer  :: i, ictot, ifac, ifbc, ifaq
    integer  :: ict1, ict2, ifac1, ifac2
    integer  :: ifbc1, ifbc2, ifaq1, ifaq2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xct, t_cat1, t_cat2
    real(r8) :: t_faq1, t_faq2, t_xfaq, t_fbc1, t_fbc2, t_xfbc
    real(r8) :: r1, r2, s1, s2, tmp, e
    real(r8) :: d2mx(4), dxm1(4), invd(4)
    real(r8) :: sizepar4d(2,2,2,2)
    real(r8), parameter :: eps=1.0e-10_r8

    ! Initialize excess mass cxs, wrt. maximum allowed internal mixing
    do lon=1,ncol
       cxs(lon) = 0.0_r8
       xct(lon)  = 0.0_r8
       xfac(lon) = 0.0_r8
       xfbc(lon) = 0.0_r8
       xfaq(lon) = 0.0_r8

       xstdv(lon) = 0._r8
       xrk(lon) = 0._r8

       xct(lon)  = min(max(xctin(lon)/(Nnat(lon)+eps),cat(kcomp,1)),cat(kcomp,6))
       xfac(lon) = min(max(xfacin(lon),fac(1)),fac(6))
       xfbc(lon) = min(max(xfbcin(lon),fbc(1)),fbc(6))
       xfaq(lon) = min(max(xfaqin(lon),faq(1)),faq(6))

       camdiff   = xctin(lon)-xct(lon)*(Nnat(lon)+eps)

       cxs(lon)  = max(0.0_r8,camdiff)

       ictot=1
       tmp = xct(lon)
       do while (ictot.lt.5.and.(tmp.lt.cat(kcomp,ictot).or. tmp.gt.cat(kcomp,ictot+1)))
          ictot=ictot+1
       enddo
       ict1=ictot
       ict2=ictot+1

       ifac=1
       tmp = xfac(lon)
       do while (ifac.lt.5.and.(tmp.lt.fac(ifac).or. tmp.gt.fac(ifac+1)))
          ifac=ifac+1
       enddo
       ifac1=ifac
       ifac2=ifac+1

       ifbc=1
       tmp = xfbc(lon)
       do while (ifbc.lt.5.and.(tmp.lt.fbc(ifbc).or. tmp.gt.fbc(ifbc+1)))
          ifbc=ifbc+1
       enddo
       ifbc1=ifbc
       ifbc2=ifbc+1

       ifaq=1
       tmp = xfaq(lon)
       do while (ifaq.lt.5.and.(tmp.lt.faq(ifaq) .or.tmp.gt.faq(ifaq+1)))
          ifaq=ifaq+1
       enddo
       ifaq1=ifaq
       ifaq2=ifaq+1

       ! Collect all the vector elements into temporary storage
       ! to avoid cache conflicts and excessive cross-referencing
       t_cat1 = cat(kcomp,ict1)
       t_cat2 = cat(kcomp,ict2)
       t_fac1 = fac(ifac1)
       t_fac2 = fac(ifac2)
       t_fbc1 = fbc(ifbc1)
       t_fbc2 = fbc(ifbc2)
       t_faq1 = faq(ifaq1)
       t_faq2 = faq(ifaq2)

       t_xct  = xct(lon)
       t_xfac = xfac(lon)
       t_xfbc = xfbc(lon)
       t_xfaq = xfaq(lon)

       ! partial lengths along each dimension (1-4) for interpolation
       d2mx(1) = (t_cat2-t_xct)
       dxm1(1) = (t_xct-t_cat1)
       invd(1) = 1.0_r8/(t_cat2-t_cat1)
       d2mx(2) = (t_fac2-t_xfac)
       dxm1(2) = (t_xfac-t_fac1)
       invd(2) = 1.0_r8/(t_fac2-t_fac1)
       d2mx(3) = (t_fbc2-t_xfbc)
       dxm1(3) = (t_xfbc-t_fbc1)
       invd(3) = 1.0_r8/(t_fbc2-t_fbc1)
       d2mx(4) = (t_faq2-t_xfaq)
       dxm1(4) = (t_xfaq-t_faq1)
       invd(4) = 1.0_r8/(t_faq2-t_faq1)

       ! Table points as basis for multidimentional linear interpolation,
       ! modal median radius:

       sizepar4d(1,1,1,1) = rrr(kcomp,ict1,ifac1,ifbc1,ifaq1)
       sizepar4d(1,1,1,2) = rrr(kcomp,ict1,ifac1,ifbc1,ifaq2)
       sizepar4d(1,1,2,1) = rrr(kcomp,ict1,ifac1,ifbc2,ifaq1)
       sizepar4d(1,1,2,2) = rrr(kcomp,ict1,ifac1,ifbc2,ifaq2)
       sizepar4d(1,2,1,1) = rrr(kcomp,ict1,ifac2,ifbc1,ifaq1)
       sizepar4d(1,2,1,2) = rrr(kcomp,ict1,ifac2,ifbc1,ifaq2)
       sizepar4d(1,2,2,1) = rrr(kcomp,ict1,ifac2,ifbc2,ifaq1)
       sizepar4d(1,2,2,2) = rrr(kcomp,ict1,ifac2,ifbc2,ifaq2)
       sizepar4d(2,1,1,1) = rrr(kcomp,ict2,ifac1,ifbc1,ifaq1)
       sizepar4d(2,1,1,2) = rrr(kcomp,ict2,ifac1,ifbc1,ifaq2)
       sizepar4d(2,1,2,1) = rrr(kcomp,ict2,ifac1,ifbc2,ifaq1)
       sizepar4d(2,1,2,2) = rrr(kcomp,ict2,ifac1,ifbc2,ifaq2)
       sizepar4d(2,2,1,1) = rrr(kcomp,ict2,ifac2,ifbc1,ifaq1)
       sizepar4d(2,2,1,2) = rrr(kcomp,ict2,ifac2,ifbc1,ifaq2)
       sizepar4d(2,2,2,1) = rrr(kcomp,ict2,ifac2,ifbc2,ifaq1)
       sizepar4d(2,2,2,2) = rrr(kcomp,ict2,ifac2,ifbc2,ifaq2)

       ! interpolation in the faq, fbc, fac and cat dimensions
       call lininterpol4dim (d2mx, dxm1, invd, sizepar4d, r1, r2)

       ! finally, interpolation in the cat dimension
       xrk(lon)=(d2mx(1)*r1+dxm1(1)*r2)*invd(1)*1.e-6_r8  ! look-up table radii in um

       ! Table points as basis for multidimentional linear interpolation,
       ! modal standard deviation:

       sizepar4d(1,1,1,1)=sss(kcomp,ict1,ifac1,ifbc1,ifaq1)
       sizepar4d(1,1,1,2)=sss(kcomp,ict1,ifac1,ifbc1,ifaq2)
       sizepar4d(1,1,2,1)=sss(kcomp,ict1,ifac1,ifbc2,ifaq1)
       sizepar4d(1,1,2,2)=sss(kcomp,ict1,ifac1,ifbc2,ifaq2)
       sizepar4d(1,2,1,1)=sss(kcomp,ict1,ifac2,ifbc1,ifaq1)
       sizepar4d(1,2,1,2)=sss(kcomp,ict1,ifac2,ifbc1,ifaq2)
       sizepar4d(1,2,2,1)=sss(kcomp,ict1,ifac2,ifbc2,ifaq1)
       sizepar4d(1,2,2,2)=sss(kcomp,ict1,ifac2,ifbc2,ifaq2)
       sizepar4d(2,1,1,1)=sss(kcomp,ict2,ifac1,ifbc1,ifaq1)
       sizepar4d(2,1,1,2)=sss(kcomp,ict2,ifac1,ifbc1,ifaq2)
       sizepar4d(2,1,2,1)=sss(kcomp,ict2,ifac1,ifbc2,ifaq1)
       sizepar4d(2,1,2,2)=sss(kcomp,ict2,ifac1,ifbc2,ifaq2)
       sizepar4d(2,2,1,1)=sss(kcomp,ict2,ifac2,ifbc1,ifaq1)
       sizepar4d(2,2,1,2)=sss(kcomp,ict2,ifac2,ifbc1,ifaq2)
       sizepar4d(2,2,2,1)=sss(kcomp,ict2,ifac2,ifbc2,ifaq1)
       sizepar4d(2,2,2,2)=sss(kcomp,ict2,ifac2,ifbc2,ifaq2)

       ! interpolation in the faq, fbc, fac and cat dimensions
       call lininterpol4dim (d2mx, dxm1, invd, sizepar4d, s1, s2)

       ! finally, interpolation in the cat dimension
       xstdv(lon)=(d2mx(1)*s1+dxm1(1)*s2)*invd(1)

    end do   ! lon
  end subroutine intlog5to10_sub

  !********************************************************************************************
  subroutine checkTableHeader (ifil)

    ! Read the header-text in a look-up table (in file with iu=ifil).

    ! arguments
    integer, intent(in) :: ifil

    ! local variables
    character*80 :: headertext
    character*12 :: text0, text1

    text0='X-CHECK LUT'
    text1='none       '
    do while (text1(2:12) .ne. text0(2:12))
       read(ifil,'(A)') headertext
       text1 = headertext(2:12)
    enddo
  end subroutine checkTableHeader

end module oslo_aero_logn_tables
