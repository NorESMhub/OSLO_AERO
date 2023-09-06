module oslo_aero_aerocom_opt

#ifdef AEROCOM

  use shr_kind_mod            , only : r8 => shr_kind_r8
  use ppgrid                  , only : pcols, pver
  use cam_logfile             , only : iulog
  !
  use oslo_aero_params       , only : nmodes, nbmodes
  use oslo_aero_sw_tables     , only : cate, cat, fac, faq, fbc, rh, fombg, fbcbg
  use oslo_aero_control       , only : oslo_aero_getopts, dir_string_length
  use oslo_aero_linear_interp , only : lininterpol3dim, lininterpol4dim, lininterpol5dim

  implicit none
  private

  ! Set by init_aeropt Mode0
  real(r8) :: bex440, bax440
  real(r8) :: bex500, bax500, bax550
  real(r8) :: bex670, bax670
  real(r8) :: bex870, bax870
  real(r8) :: bex550lt1, bex550gt1, backscx550

  ! Set by init_aeropt Mode1
  real(r8), public :: bep1(38,10,6,16,6)

  ! Set by init_aeropt Mode2to3
  real(r8), public :: bep2to3 (38,10,16,6,2:3)

  ! Set by init_aeropt Mode4
  real(r8), public :: bep4(38,10,6,16,6,6)

  ! Set by init_aeropt Mode5to10
  real(r8), public :: bep5to10(38,10,6,6,6,6,5:10)

  ! Modal total and absorption extiction coefficients (for AeroCom)
  ! for 440nm, 500nm, 550nm, 670nm and 870nm, and for d<1um (lt1) and d>1um (gt1).
  ! March 2009: + backscatter coefficient, backsc550 (km-1 sr-1).
  type, public :: extinction_coeffs_type
     real(r8) :: bext440(pcols,pver,0:nbmodes)
     real(r8) :: babs440(pcols,pver,0:nbmodes)
     real(r8) :: bext500(pcols,pver,0:nbmodes)
     real(r8) :: babs500(pcols,pver,0:nbmodes)
     real(r8) :: bext550(pcols,pver,0:nbmodes)
     real(r8) :: babs550(pcols,pver,0:nbmodes)
     real(r8) :: bext670(pcols,pver,0:nbmodes)
     real(r8) :: babs670(pcols,pver,0:nbmodes)
     real(r8) :: bext870(pcols,pver,0:nbmodes)
     real(r8) :: babs870(pcols,pver,0:nbmodes)
     real(r8) :: bebg440(pcols,pver,0:nbmodes)
     real(r8) :: bebg500(pcols,pver,0:nbmodes)
     real(r8) :: bebg550(pcols,pver,0:nbmodes)
     real(r8) :: babg550(pcols,pver,0:nbmodes)
     real(r8) :: bebg670(pcols,pver,0:nbmodes)
     real(r8) :: bebg870(pcols,pver,0:nbmodes)
     real(r8) :: bebc440(pcols,pver,0:nbmodes)
     real(r8) :: bebc500(pcols,pver,0:nbmodes)
     real(r8) :: bebc550(pcols,pver,0:nbmodes)
     real(r8) :: babc550(pcols,pver,0:nbmodes)
     real(r8) :: bebc670(pcols,pver,0:nbmodes)
     real(r8) :: bebc870(pcols,pver,0:nbmodes)
     real(r8) :: beoc440(pcols,pver,0:nbmodes)
     real(r8) :: beoc500(pcols,pver,0:nbmodes)
     real(r8) :: beoc550(pcols,pver,0:nbmodes)
     real(r8) :: baoc550(pcols,pver,0:nbmodes)
     real(r8) :: beoc670(pcols,pver,0:nbmodes)
     real(r8) :: beoc870(pcols,pver,0:nbmodes)
     real(r8) :: besu440(pcols,pver,0:nbmodes)
     real(r8) :: besu500(pcols,pver,0:nbmodes)
     real(r8) :: besu550(pcols,pver,0:nbmodes)
     real(r8) :: basu550(pcols,pver,0:nbmodes)
     real(r8) :: besu670(pcols,pver,0:nbmodes)
     real(r8) :: besu870(pcols,pver,0:nbmodes)
     real(r8) :: bebg550lt1(pcols,pver,0:nbmodes)
     real(r8) :: bebg550gt1(pcols,pver,0:nbmodes)
     real(r8) :: bebc550lt1(pcols,pver,0:nbmodes)
     real(r8) :: bebc550gt1(pcols,pver,0:nbmodes)
     real(r8) :: beoc550lt1(pcols,pver,0:nbmodes)
     real(r8) :: beoc550gt1(pcols,pver,0:nbmodes)
     real(r8) :: besu550lt1(pcols,pver,0:nbmodes)
     real(r8) :: besu550gt1(pcols,pver,0:nbmodes)
     real(r8) :: backsc550(pcols,pver,0:nbmodes)

   contains
     procedure :: intaeropt0
     procedure :: intaeropt1
     procedure :: intaeropt2to3
     procedure :: intaeropt4
     procedure :: intaeropt5to10
     procedure :: zero
     procedure :: update
  end type extinction_coeffs_type

  type(extinction_coeffs_type), public :: extinction_coeffs
  type(extinction_coeffs_type), public :: extinction_coeffsn

  public :: initaeropt

! ==========================================================
contains
! ==========================================================

  subroutine initaeropt()

    !Purpose: To read in the AeroCom look-up tables for aerosol optical properties.
    !     The grid for discrete input-values in the look-up tables is defined in opptab.

    !     Tabulating the 'aerocomk'-files to save computing time.
    !     Updated for new kcomp1.out including condensed SOA - Alf Kirkevåg, May 2013
    !     Extended for new SOA treatment - Alf Kirkevaag, September 2015.
    !     Modified for optimized added masses and mass fractions for
    !     concentrations from condensation, coagulation or cloud-processing
    !     - Alf Kirkevaag, May 2016.
    !     Modified for optimized added masses and mass fractions for concentrations from
    !     condensation, coagulation or cloud-processing - Alf Kirkevaag, May 2016.

    integer  :: ic, ifil, lin, iv
    integer  :: kcomp, irelh, ictot, ifac, ifbc, ifaq
    integer  :: ifombg, ifbcbg
    real(r8) :: catot, relh, frombg, frbcbg, frac, fabc, fraq
    real(r8) :: bext440, babs440, bext500, babs500, babs550
    real(r8) :: bext670, babs670, bext870, babs870
    real(r8) :: bebg440, babg440, bebg500, babg500, babg550
    real(r8) :: bebg670, babg670, bebg870, babg870
    real(r8) :: bebc440, babc440, bebc500, babc500, babc550
    real(r8) :: bebc670, babc670, bebc870, babc870
    real(r8) :: beoc440, baoc440, beoc500, baoc500, baoc550
    real(r8) :: beoc670, baoc670, beoc870, baoc870
    real(r8) :: besu440, basu440, besu500, basu500, basu550
    real(r8) :: besu670, basu670, besu870
    real(r8) :: bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1
    real(r8) :: beoc550lt1, beoc550gt1, besu550lt1, besu550gt1
    real(r8) :: backscat550
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps4 = 1.e-4_r8
    real(r8) :: eps6 = 1.e-6_r8
    real(r8) :: eps7 = 1.e-7_r8
    character(len=dir_string_length) :: aerotab_table_dir
    !-----------------------------------------------------------

    call oslo_aero_getopts(aerotab_table_dir_out = aerotab_table_dir)

    open(20,file=trim(aerotab_table_dir)//'/aerocomk0.out' , form='formatted',status='old')
    open(21,file=trim(aerotab_table_dir)//'/aerocomk1.out' , form='formatted',status='old')
    open(11,file=trim(aerotab_table_dir)//'/aerocomk2.out' , form='formatted',status='old')
    open(12,file=trim(aerotab_table_dir)//'/aerocomk3.out' , form='formatted',status='old')
    open(13,file=trim(aerotab_table_dir)//'/aerocomk4.out' , form='formatted',status='old')
    open(14,file=trim(aerotab_table_dir)//'/aerocomk5.out' , form='formatted',status='old')
    open(15,file=trim(aerotab_table_dir)//'/aerocomk6.out' , form='formatted',status='old')
    open(16,file=trim(aerotab_table_dir)//'/aerocomk7.out' , form='formatted',status='old')
    open(17,file=trim(aerotab_table_dir)//'/aerocomk8.out' , form='formatted',status='old')
    open(18,file=trim(aerotab_table_dir)//'/aerocomk9.out' , form='formatted',status='old')
    open(19,file=trim(aerotab_table_dir)//'/aerocomk10.out', form='formatted',status='old')

    ! Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 11,21
       call checkTableHeader (ifil)
    enddo
    !
    !-------------------------------------------
    ! Mode 0, BC
    !-------------------------------------------
    !
    read(20,'(I2,f6.3,12e11.4)')                                 &
         kcomp, relh,                                            &
         bex440, bax440, bex500, bax500, bax550, bex670, bax670, &
         bex870, bax870, bex550lt1, bex550gt1, backscx550

    if(bex440<=0.0_r8) then
       write(*,*) 'bex440 =', bex440
       write(*,*) 'Error in initialization of bex1'
       stop
    endif
    write(iulog,*)'aerocom mode 0 ok'
    !
    !-------------------------------------------
    ! Mode 1 (H2SO4 and SOA + condensate from H2SO4 and SOA)
    !-------------------------------------------
    !
    do lin = 1,5760     ! 10x6x16x6
       read(21,'(I2,f6.3,3e10.3,38e10.3)')                  &
            kcomp, relh, frombg, catot, frac,               &
            bext440, bext500, bext670, bext870,             &
            bebg440, bebg500, bebg670, bebg870,             &
            bebc440, bebc500, bebc670, bebc870,             &
            beoc440, beoc500, beoc670, beoc870,             &
            besu440, besu500, besu670, besu870,             &
            babs440, babs500, babs550, babs670, babs870,    &
            bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
            beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
            backscat550, babg550, babc550, baoc550, basu550

       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frombg-fombg(ic))<eps4) then
             ifombg=ic
             exit
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do

       bep1(1,irelh,ifombg,ictot,ifac) = bext440 ! unit km^-1
       bep1(2,irelh,ifombg,ictot,ifac) = bext500
       bep1(3,irelh,ifombg,ictot,ifac) = bext670
       bep1(4,irelh,ifombg,ictot,ifac) = bext870
       bep1(5,irelh,ifombg,ictot,ifac) = bebg440
       bep1(6,irelh,ifombg,ictot,ifac) = bebg500
       bep1(7,irelh,ifombg,ictot,ifac) = bebg670
       bep1(8,irelh,ifombg,ictot,ifac) = bebg870
       bep1(9,irelh,ifombg,ictot,ifac) = bebc440  ! = 0
       bep1(10,irelh,ifombg,ictot,ifac) = bebc500 ! = 0
       bep1(11,irelh,ifombg,ictot,ifac) = bebc670 ! = 0
       bep1(12,irelh,ifombg,ictot,ifac) = bebc870 ! = 0
       bep1(13,irelh,ifombg,ictot,ifac) = beoc440
       bep1(14,irelh,ifombg,ictot,ifac) = beoc500
       bep1(15,irelh,ifombg,ictot,ifac) = beoc670
       bep1(16,irelh,ifombg,ictot,ifac) = beoc870
       bep1(17,irelh,ifombg,ictot,ifac) = besu440
       bep1(18,irelh,ifombg,ictot,ifac) = besu500
       bep1(19,irelh,ifombg,ictot,ifac) = besu670
       bep1(20,irelh,ifombg,ictot,ifac) = besu870
       bep1(21,irelh,ifombg,ictot,ifac) = babs440
       bep1(22,irelh,ifombg,ictot,ifac) = babs500
       bep1(23,irelh,ifombg,ictot,ifac) = babs550
       bep1(24,irelh,ifombg,ictot,ifac) = babs670
       bep1(25,irelh,ifombg,ictot,ifac) = babs870
       bep1(26,irelh,ifombg,ictot,ifac) = bebg550lt1
       bep1(27,irelh,ifombg,ictot,ifac) = bebg550gt1
       bep1(28,irelh,ifombg,ictot,ifac) = bebc550lt1 ! = 0
       bep1(29,irelh,ifombg,ictot,ifac) = bebc550gt1 ! = 0
       bep1(30,irelh,ifombg,ictot,ifac) = beoc550lt1
       bep1(31,irelh,ifombg,ictot,ifac) = beoc550gt1
       bep1(32,irelh,ifombg,ictot,ifac) = besu550lt1
       bep1(33,irelh,ifombg,ictot,ifac) = besu550gt1
       bep1(34,irelh,ifombg,ictot,ifac) = backscat550
       bep1(35,irelh,ifombg,ictot,ifac) = babg550
       bep1(36,irelh,ifombg,ictot,ifac) = babc550 ! = 0
       bep1(37,irelh,ifombg,ictot,ifac) = baoc550
       bep1(38,irelh,ifombg,ictot,ifac) = basu550
    end do  ! lin

    do irelh=1,10
       do ifombg=1,6
          do ictot=1,16
             do ifac=1,6
                if(bep1(1,irelh,ifombg,ictot,ifac)<=0.0_r8) then
                   write(*,*) 'bep1 =', irelh,ifombg, ictot, ifac, &
                        bep1(1,irelh,ifombg,ictot,ifac)
                   write(*,*) 'Error in initialization of bep1'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 1 ok'
    !
    !-------------------------------------------
    ! Mode 2  (BC/OC + condesate from H2SO4 and SOA)
    ! Note that mode 3 is no longer active
    !-------------------------------------------
    !
    do lin = 1,960     ! 10x16x6
       read(11,'(I2,f6.3,2e10.3,38e10.3)')                  &
            kcomp, relh, catot, frac,                       &
            bext440, bext500, bext670, bext870,             &
            bebg440, bebg500, bebg670, bebg870,             &
            bebc440, bebc500, bebc670, bebc870,             &
            beoc440, beoc500, beoc670, beoc870,             &
            besu440, besu500, besu670, besu870,             &
            babs440, babs500, babs550, babs670, babs870,    &
            bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
            beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
            backscat550, babg550, babc550, baoc550, basu550

       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             exit
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do

       bep2to3(1,irelh,ictot,ifac,kcomp) = bext440 ! unit km^-1
       bep2to3(2,irelh,ictot,ifac,kcomp) = bext500
       bep2to3(3,irelh,ictot,ifac,kcomp) = bext670
       bep2to3(4,irelh,ictot,ifac,kcomp) = bext870
       bep2to3(5,irelh,ictot,ifac,kcomp) = bebg440
       bep2to3(6,irelh,ictot,ifac,kcomp) = bebg500
       bep2to3(7,irelh,ictot,ifac,kcomp) = bebg670
       bep2to3(8,irelh,ictot,ifac,kcomp) = bebg870
       bep2to3(9,irelh,ictot,ifac,kcomp) = bebc440  ! = 0
       bep2to3(10,irelh,ictot,ifac,kcomp) = bebc500 ! = 0
       bep2to3(11,irelh,ictot,ifac,kcomp) = bebc670 ! = 0
       bep2to3(12,irelh,ictot,ifac,kcomp) = bebc870 ! = 0
       bep2to3(13,irelh,ictot,ifac,kcomp) = beoc440
       bep2to3(14,irelh,ictot,ifac,kcomp) = beoc500
       bep2to3(15,irelh,ictot,ifac,kcomp) = beoc670
       bep2to3(16,irelh,ictot,ifac,kcomp) = beoc870
       bep2to3(17,irelh,ictot,ifac,kcomp) = besu440
       bep2to3(18,irelh,ictot,ifac,kcomp) = besu500
       bep2to3(19,irelh,ictot,ifac,kcomp) = besu670
       bep2to3(20,irelh,ictot,ifac,kcomp) = besu870
       bep2to3(21,irelh,ictot,ifac,kcomp) = babs440
       bep2to3(22,irelh,ictot,ifac,kcomp) = babs500
       bep2to3(23,irelh,ictot,ifac,kcomp) = babs550
       bep2to3(24,irelh,ictot,ifac,kcomp) = babs670
       bep2to3(25,irelh,ictot,ifac,kcomp) = babs870
       bep2to3(26,irelh,ictot,ifac,kcomp) = bebg550lt1
       bep2to3(27,irelh,ictot,ifac,kcomp) = bebg550gt1
       bep2to3(28,irelh,ictot,ifac,kcomp) = bebc550lt1 ! = 0
       bep2to3(29,irelh,ictot,ifac,kcomp) = bebc550gt1 ! = 0
       bep2to3(30,irelh,ictot,ifac,kcomp) = beoc550lt1
       bep2to3(31,irelh,ictot,ifac,kcomp) = beoc550gt1
       bep2to3(32,irelh,ictot,ifac,kcomp) = besu550lt1
       bep2to3(33,irelh,ictot,ifac,kcomp) = besu550gt1
       bep2to3(34,irelh,ictot,ifac,kcomp) = backscat550
       bep2to3(35,irelh,ictot,ifac,kcomp) = babg550
       bep2to3(36,irelh,ictot,ifac,kcomp) = babc550 ! = 0
       bep2to3(37,irelh,ictot,ifac,kcomp) = baoc550
       bep2to3(38,irelh,ictot,ifac,kcomp) = basu550
    end do

    ! Prescribed dummy values for unused kcomp=3
    kcomp=3
    do irelh=1,10
       do ictot=1,16
          do ifac=1,6
             do iv=1,38
                bep2to3(iv,irelh,ictot,ifac,kcomp)=1.0_r8
             enddo
          enddo
       enddo
    enddo

    do kcomp=2,3
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                if(bep2to3(1,irelh,ictot,ifac,kcomp)<=0.0_r8) then
                   write(*,*) 'bep2to3 =', irelh, ictot, ifac, bep2to3(1,irelh,ictot,ifac,kcomp)
                   write(*,*) 'Error in initialization of bep2to3'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 2-3 ok'
    !
    !-------------------------------------------
    ! Mode 4 (BC&OC + condesate from H2SO4 and SOA + wetphase (NH4)2SO4)
    !-------------------------------------------
    !
    do lin = 1,34560     ! 10x16x6x6x6
       read(13,'(I2,f6.3,3e10.3,f5.2,38e10.3)')            &
            kcomp, relh, frbcbg, catot, frac, fraq,         &
            bext440, bext500, bext670, bext870,             &
            bebg440, bebg500, bebg670, bebg870,             &
            bebc440, bebc500, bebc670, bebc870,             &
            beoc440, beoc500, beoc670, beoc870,             &
            besu440, besu500, besu670, besu870,             &
            babs440, babs500, babs550, babs670, babs870,    &
            bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
            beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
            backscat550, babg550, babc550, baoc550, basu550

       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs((frbcbg-fbcbg(ic))/fbcbg(ic))<eps2) then
             ifbcbg=ic
             exit
          endif
       end do
       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(fraq-faq(ic))<eps4) then
             ifaq=ic
             exit
          endif
       end do

       bep4(1,irelh,ifbcbg,ictot,ifac,ifaq) = bext440 ! unit km^-1
       bep4(2,irelh,ifbcbg,ictot,ifac,ifaq) = bext500
       bep4(3,irelh,ifbcbg,ictot,ifac,ifaq) = bext670
       bep4(4,irelh,ifbcbg,ictot,ifac,ifaq) = bext870
       bep4(5,irelh,ifbcbg,ictot,ifac,ifaq) = bebg440
       bep4(6,irelh,ifbcbg,ictot,ifac,ifaq) = bebg500
       bep4(7,irelh,ifbcbg,ictot,ifac,ifaq) = bebg670
       bep4(8,irelh,ifbcbg,ictot,ifac,ifaq) = bebg870
       bep4(9,irelh,ifbcbg,ictot,ifac,ifaq) = bebc440
       bep4(10,irelh,ifbcbg,ictot,ifac,ifaq) = bebc500
       bep4(11,irelh,ifbcbg,ictot,ifac,ifaq) = bebc670
       bep4(12,irelh,ifbcbg,ictot,ifac,ifaq) = bebc870
       bep4(13,irelh,ifbcbg,ictot,ifac,ifaq) = beoc440 ! = 0
       bep4(14,irelh,ifbcbg,ictot,ifac,ifaq) = beoc500 ! = 0
       bep4(15,irelh,ifbcbg,ictot,ifac,ifaq) = beoc670 ! = 0
       bep4(16,irelh,ifbcbg,ictot,ifac,ifaq) = beoc870 ! = 0
       bep4(17,irelh,ifbcbg,ictot,ifac,ifaq) = besu440
       bep4(18,irelh,ifbcbg,ictot,ifac,ifaq) = besu500
       bep4(19,irelh,ifbcbg,ictot,ifac,ifaq) = besu670
       bep4(20,irelh,ifbcbg,ictot,ifac,ifaq) = besu870
       bep4(21,irelh,ifbcbg,ictot,ifac,ifaq) = babs440
       bep4(22,irelh,ifbcbg,ictot,ifac,ifaq) = babs500
       bep4(23,irelh,ifbcbg,ictot,ifac,ifaq) = babs550
       bep4(24,irelh,ifbcbg,ictot,ifac,ifaq) = babs670
       bep4(25,irelh,ifbcbg,ictot,ifac,ifaq) = babs870
       bep4(26,irelh,ifbcbg,ictot,ifac,ifaq) = bebg550lt1
       bep4(27,irelh,ifbcbg,ictot,ifac,ifaq) = bebg550gt1
       bep4(28,irelh,ifbcbg,ictot,ifac,ifaq) = bebc550lt1
       bep4(29,irelh,ifbcbg,ictot,ifac,ifaq) = bebc550gt1
       bep4(30,irelh,ifbcbg,ictot,ifac,ifaq) = beoc550lt1 ! = 0
       bep4(31,irelh,ifbcbg,ictot,ifac,ifaq) = beoc550gt1 ! = 0
       bep4(32,irelh,ifbcbg,ictot,ifac,ifaq) = besu550lt1
       bep4(33,irelh,ifbcbg,ictot,ifac,ifaq) = besu550gt1
       bep4(34,irelh,ifbcbg,ictot,ifac,ifaq) = backscat550
       bep4(35,irelh,ifbcbg,ictot,ifac,ifaq) = babg550
       bep4(36,irelh,ifbcbg,ictot,ifac,ifaq) = babc550 ! = 0
       bep4(37,irelh,ifbcbg,ictot,ifac,ifaq) = baoc550 ! = 0
       bep4(38,irelh,ifbcbg,ictot,ifac,ifaq) = basu550
    end do

    do irelh=1,10
       do ifbcbg=1,6
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(bep4(1,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                      write(*,*) 'bep4 =', irelh, ifbcbg, ictot, ifac, ifaq, &
                           bep4(1,irelh,ifbcbg,ictot,ifac,ifaq)
                      write(*,*) 'Error in initialization of bep4'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 4 ok'
    !
    !-------------------------------------------
    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !-------------------------------------------
    !
    do ifil = 5,10
       do lin = 1,12960     ! 10x6x6x6x6
          read(9+ifil,'(I2,f6.3,3e10.3,f5.2,38e10.3)')          &
               kcomp, relh, catot, frac, fabc, fraq,           &
               bext440, bext500, bext670, bext870,             &
               bebg440, bebg500, bebg670, bebg870,             &
               bebc440, bebc500, bebc670, bebc870,             &
               beoc440, beoc500, beoc670, beoc870,             &
               besu440, besu500, besu670, besu870,             &
               babs440, babs500, babs550, babs670, babs870,    &
               bebg550lt1, bebg550gt1, bebc550lt1, bebc550gt1, &
               beoc550lt1, beoc550gt1, besu550lt1, besu550gt1, &
               backscat550, babg550, babc550, baoc550, basu550

       	  do ic=1,10
             if(abs(relh-rh(ic))<eps4) then
                irelh=ic
                exit
             endif
          end do
          do ic=1,6
             if(abs((catot-cat(kcomp,ic))/cat(kcomp,ic))<eps2) then
                ictot=ic
                exit
             end if
          end do
          do ic=1,6
             if(abs(frac-fac(ic))<eps4) then
                ifac=ic
                exit
             endif
          end do
          do ic=1,6
             if(abs((fabc-fbc(ic))/fbc(ic))<eps2) then
                ifbc=ic
                exit
             endif
          end do
          do ic=1,6
             if(abs(fraq-faq(ic))<eps4) then
                ifaq=ic
                exit
             endif
          end do

          bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bext440 ! unit km^-1
          bep5to10(2,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bext500
          bep5to10(3,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bext670
          bep5to10(4,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bext870
          bep5to10(5,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg440
          bep5to10(6,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg500
          bep5to10(7,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg670
          bep5to10(8,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg870
          bep5to10(9,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc440
          bep5to10(10,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc500
          bep5to10(11,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc670
          bep5to10(12,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc870
          bep5to10(13,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc440
          bep5to10(14,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc500
          bep5to10(15,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc670
          bep5to10(16,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc870
          bep5to10(17,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu440
          bep5to10(18,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu500
          bep5to10(19,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu670
          bep5to10(20,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu870
          bep5to10(21,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babs440
          bep5to10(22,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babs500
          bep5to10(23,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babs550
          bep5to10(24,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babs670
          bep5to10(25,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babs870
          bep5to10(26,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg550lt1
          bep5to10(27,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebg550gt1
          bep5to10(28,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc550lt1
          bep5to10(29,irelh,ictot,ifac,ifbc,ifaq,kcomp) = bebc550gt1
          bep5to10(30,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc550lt1
          bep5to10(31,irelh,ictot,ifac,ifbc,ifaq,kcomp) = beoc550gt1
          bep5to10(32,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu550lt1
          bep5to10(33,irelh,ictot,ifac,ifbc,ifaq,kcomp) = besu550gt1
          bep5to10(34,irelh,ictot,ifac,ifbc,ifaq,kcomp) = backscat550
          bep5to10(35,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babg550
          bep5to10(36,irelh,ictot,ifac,ifbc,ifaq,kcomp) = babc550
          bep5to10(37,irelh,ictot,ifac,ifbc,ifaq,kcomp) = baoc550
          bep5to10(38,irelh,ictot,ifac,ifbc,ifaq,kcomp) = basu550
       end do
    end do

    do kcomp=5,10
       do irelh=1,10
          do ictot=1,6
             do ifac=1,6
                do ifaq=1,6
                   if(bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                      write(*,*) 'bep5to10 =', kcomp, irelh, ictot, ifac, ifbc, ifaq, &
                           bep5to10(1,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                      write(*,*) 'Error in initialization of bep5to10'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerocom mode 5-10 ok'

    ! Close files
    do ifil=10,21
       close (ifil)
    end do

  end subroutine initaeropt

  ! ==========================================================
  subroutine intaeropt0 (this, lchnk, ncol, Nnatk)

    ! Arguments
    class(extinction_coeffs_type) :: this
    integer                , intent(in)    :: lchnk                     ! chunk identifier
    integer                , intent(in)    :: ncol                      ! number of atmospheric columns
    real(r8)               , intent(in)    :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration

    ! Local variables
    integer i, iv, ierr, k, kcomp, icol

    kcomp=0
    call this%zero(kcomp, ncol)

    ! Mode 0 BC(ax)
    do k = 1,pver
       do icol = 1,ncol
          if (Nnatk(icol,k,kcomp).gt.0) then
             this%bext440(icol,k,kcomp)=bex440
             this%babs440(icol,k,kcomp)=bax440
             this%bext500(icol,k,kcomp)=bex500
             this%babs500(icol,k,kcomp)=bax500
             this%bext550(icol,k,kcomp)=bex550lt1+bex550gt1
             this%babs550(icol,k,kcomp)=bax550
             this%bext670(icol,k,kcomp)=bex670
             this%babs670(icol,k,kcomp)=bax670
             this%bext870(icol,k,kcomp)=bex870
             this%babs870(icol,k,kcomp)=bax870
             this%bebg440(icol,k,kcomp)=bex440
             this%bebg500(icol,k,kcomp)=bex500
             this%bebg550(icol,k,kcomp)=bex550lt1+bex550gt1
             this%babg550(icol,k,kcomp)=bax550
             this%bebg670(icol,k,kcomp)=bex670
             this%bebg870(icol,k,kcomp)=bex870
             this%bebg550lt1(icol,k,kcomp)=bex550lt1
             this%bebg550gt1(icol,k,kcomp)=bex550gt1
             this%backsc550(icol,k,kcomp)=backscx550
          endif
       end do ! icol
    end do ! k

  end subroutine intaeropt0

  ! ==========================================================
  subroutine intaeropt1 (this, lchnk, ncol, xrh, irh1, mplus10, &
       Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1)

    ! arguments
    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    integer  , intent(in) :: mplus10                    ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xfombg(pcols,pver)         ! SOA/(SOA+H2SO4) for the background mode
    integer  , intent(in) :: ifombg1(pcols,pver)
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, ierr, irelh, ifombg, ictot, ifac, kcomp, k, icol, kc10
    integer  :: t_irh1, t_irh2, t_ifo1, t_ifo2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_rh1, t_rh2, t_fombg1, t_fombg2, t_xfombg
    real(r8) :: t_xct, t_cat1, t_cat2
    real(r8) :: d2mx(4), dxm1(4), invd(4)
    real(r8) :: opt4d(2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)
    !----------------------------------------------

    ! SO4/SOA(Ait) mode:
    kcomp = 1
    call this%zero(kcomp, ncol)

    if(mplus10 == 0) then
       kc10 = kcomp
    else
       write(*,*) "mplus10=1 is no loger an option for kcomp=1."
       stop
    endif

    do k=1,pver
       do icol=1,ncol

          if(Nnatk(icol,k,kc10).gt.0) then

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ifo1 = ifombg1(icol,k)
             t_ifo2 = t_ifo1+1
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_fombg1 = fombg(t_ifo1)
             t_fombg2 = fombg(t_ifo2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)
             t_xfombg = xfombg(icol,k)

             ! partial lengths along each dimension (1-4) for interpolation
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_fombg2-t_xfombg)
             dxm1(2) = (t_xfombg-t_fombg1)
             invd(2) = 1.0_r8/(t_fombg2-t_fombg1)
             d2mx(3) = (t_cat2-t_xct)
             dxm1(3) = (t_xct-t_cat1)
             invd(3) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(4) = (t_fac2-t_xfac)
             dxm1(4) = (t_xfac-t_fac1)
             invd(4) = 1.0_r8/(t_fac2-t_fac1)

             do iv=1,38  ! variable number
                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1) = bep1(iv,t_irh1,t_ifo1,t_ict1,t_ifc1)
                opt4d(1,1,1,2) = bep1(iv,t_irh1,t_ifo1,t_ict1,t_ifc2)
                opt4d(1,1,2,1) = bep1(iv,t_irh1,t_ifo1,t_ict2,t_ifc1)
                opt4d(1,1,2,2) = bep1(iv,t_irh1,t_ifo1,t_ict2,t_ifc2)
                opt4d(1,2,1,1) = bep1(iv,t_irh1,t_ifo2,t_ict1,t_ifc1)
                opt4d(1,2,1,2) = bep1(iv,t_irh1,t_ifo2,t_ict1,t_ifc2)
                opt4d(1,2,2,1) = bep1(iv,t_irh1,t_ifo2,t_ict2,t_ifc1)
                opt4d(1,2,2,2) = bep1(iv,t_irh1,t_ifo2,t_ict2,t_ifc2)
                opt4d(2,1,1,1) = bep1(iv,t_irh2,t_ifo1,t_ict1,t_ifc1)
                opt4d(2,1,1,2) = bep1(iv,t_irh2,t_ifo1,t_ict1,t_ifc2)
                opt4d(2,1,2,1) = bep1(iv,t_irh2,t_ifo1,t_ict2,t_ifc1)
                opt4d(2,1,2,2) = bep1(iv,t_irh2,t_ifo1,t_ict2,t_ifc2)
                opt4d(2,2,1,1) = bep1(iv,t_irh2,t_ifo2,t_ict1,t_ifc1)
                opt4d(2,2,1,2) = bep1(iv,t_irh2,t_ifo2,t_ict1,t_ifc2)
                opt4d(2,2,2,1) = bep1(iv,t_irh2,t_ifo2,t_ict2,t_ifc1)
                opt4d(2,2,2,2) = bep1(iv,t_irh2,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac, cat and fombg dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, opt1, opt2)

                ! finally, interpolation in the rh dimension
                opt(iv)=((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) / (t_rh2-t_rh1)
             end do ! iv=1,38

             ! update extinction coefficient
             call this%update(icol, k, kcomp, opt)

          end if
       end do ! end of icol loop
    end do  ! end of k loop

  end subroutine intaeropt1

  ! ==========================================================
  subroutine intaeropt2to3 (this, lchnk, ncol, xrh, irh1, mplus10, &
       Nnatk, xct, ict1, xfac, ifac1)

    !   Extended by Alf Kirkevaag to include SOA in September 2015

    ! Arguments
    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: lchnk                       ! chunk identifier
    integer  , intent(in) :: ncol                        ! number of atmospheric columns
    integer  , intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol, kc10
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(3), dxm1(3), invd(3)
    real(r8) :: opt3d(2,2,2)
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! SO4(Ait), BC(Ait) and OC(Ait) modes:

    do kcomp=2,3
       call this%zero(kcomp, ncol)
    end do

    kcomp = 2 ! kcomp=3 is no longer used
    do k=1,pver
       do icol=1,ncol

          if(Nnatk(icol,k,kc10).gt.0) then

             !      Collect all the vector elements into temporary storage
             !      to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)

             !     partial lengths along each dimension (1-4) for interpolation
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)

             do iv=1,38  ! variable number

                !  end points as basis for multidimentional linear interpolation
                opt3d(1,1,1)=bep2to3(iv,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=bep2to3(iv,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=bep2to3(iv,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=bep2to3(iv,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=bep2to3(iv,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=bep2to3(iv,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=bep2to3(iv,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=bep2to3(iv,t_irh2,t_ict2,t_ifc2,kcomp)

                !     interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, opt1, opt2)

                !     finally, interpolation in the rh dimension
                opt(iv)=((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)

             end do ! iv=1,38

             ! determine extinction coefficient
             call this%update(icol, k, kcomp, opt)

          end if ! Nnatk > 0
       end do ! icol
    end do ! k

  end subroutine intaeropt2to3

  ! ==========================================================
  subroutine intaeropt4 (this, lchnk, ncol, xrh, irh1, mplus10, Nnatk,   &
       xfbcbg, ifbcbg1, xct, ict1, xfac, ifac1, xfaq, ifaq1)

    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    integer  , intent(in) :: mplus10                    ! mode number (0) or number + 10 (1)
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xfbcbg(pcols,pver)
    integer  , intent(in) :: ifbcbg1(pcols,pver)
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol, kc10
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2,  t_ifa1, t_ifa2
    real(r8) :: t_fbcbg1, t_fbcbg2
    integer  :: t_ifb1, t_ifb2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_xct, t_rh1, t_rh2
    real(r8) :: t_cat1, t_cat2
    real(r8) :: t_xfbcbg
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! BC&OC(Ait) mode:
    kcomp = 4
    call this%zero(kcomp, ncol)

    if(mplus10==0) then
       kc10=kcomp
    else
       kc10=kcomp+10
    endif

    do k=1,pver
       do icol=1,ncol
          if(Nnatk(icol,k,kc10).gt.0) then
             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ifb1 = ifbcbg1(icol,k)
             t_ifb2 = t_ifb1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifa1 = ifaq1(icol,k,kcomp)
             t_ifa2 = t_ifa1+1

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_fbcbg1 = fbcbg(t_ifb1)
             t_fbcbg2 = fbcbg(t_ifb2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_faq1 = faq(t_ifa1)
             t_faq2 = faq(t_ifa2)

             t_xrh  = xrh(icol,k)
             t_xfbcbg = xfbcbg(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)
             t_xfaq = xfaq(icol,k,kcomp)

             ! partial lengths along each dimension (1-5) for interpolation
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_fbcbg2-t_xfbcbg)
             dxm1(2) = (t_xfbcbg-t_fbcbg1)
             invd(2) = 1.0_r8/(t_fbcbg2-t_fbcbg1)
             d2mx(3) = (t_cat2-t_xct)
             dxm1(3) = (t_xct-t_cat1)
             invd(3) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(4) = (t_fac2-t_xfac)
             dxm1(4) = (t_xfac-t_fac1)
             invd(4) = 1.0_r8/(t_fac2-t_fac1)
             d2mx(5) = (t_faq2-t_xfaq)
             dxm1(5) = (t_xfaq-t_faq1)
             invd(5) = 1.0_r8/(t_faq2-t_faq1)


             do iv=1,38  ! variable number

                opt5d(1,1,1,1,1)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=bep4(iv,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=bep4(iv,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=bep4(iv,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=bep4(iv,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=bep4(iv,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=bep4(iv,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=bep4(iv,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=bep4(iv,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, opt1, opt2)

                ! finally, interpolation in the rh dimension
                opt(iv) = ((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)

             end do ! iv=1,38

             ! determine extinction coefficient
             call this%update(icol, k, kcomp, opt)

          end if ! Nnatk > 0
       end do ! icol
    end do ! k

  end subroutine intaeropt4

  ! ==========================================================
  subroutine intaeropt5to10 (this, lchnk, ncol, xrh, irh1, Nnatk,    &
       xct, ict1, xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1)

    ! Arguments
    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    real(r8) , intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer  , intent(in) :: irh1(pcols,pver)
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! modal (OC+BC)/(SO4+BC+OC)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfbc(pcols,pver,nbmodes)   ! modal BC/(OC+BC)
    integer  , intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: i, iv, kcomp, k, icol
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
    integer  :: t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fbc1, t_fbc2, t_xfbc
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xrh, t_xct, t_rh1, t_rh2
    real(r8) :: t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: opt1, opt2, opt(38)
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

    do kcomp=5,10
       ! zero extinction coefficients for this kcomp
       call this%zero(kcomp, ncol)

       do k=1,pver
          do icol=1,ncol
             if(Nnatk(icol,k,kcomp).gt.0) then
                ! Collect all the vector elements into temporary storage
                ! to avoid cache conflicts and excessive cross-referencing

                t_irh1 = irh1(icol,k)
                t_irh2 = t_irh1+1
                t_ict1 = ict1(icol,k,kcomp)
                t_ict2 = t_ict1+1
                t_ifc1 = ifac1(icol,k,kcomp)
                t_ifc2 = t_ifc1+1

                t_ifb1 = ifbc1(icol,k,kcomp)
                t_ifb2 = t_ifb1+1
                t_ifa1 = ifaq1(icol,k,kcomp)
                t_ifa2 = t_ifa1+1

                t_rh1  = rh(t_irh1)
                t_rh2  = rh(t_irh2)
                t_cat1 = cat(kcomp,t_ict1)
                t_cat2 = cat(kcomp,t_ict2)
                t_fac1 = fac(t_ifc1)
                t_fac2 = fac(t_ifc2)
                t_fbc1 = fbc(t_ifb1)
                t_fbc2 = fbc(t_ifb2)
                t_faq1 = faq(t_ifa1)
                t_faq2 = faq(t_ifa2)

                t_xrh  = xrh(icol,k)
                t_xct  = xct(icol,k,kcomp)
                t_xfac = xfac(icol,k,kcomp)
                t_xfbc = xfbc(icol,k,kcomp)
                t_xfaq = xfaq(icol,k,kcomp)

                ! partial lengths along each dimension (1-5) for interpolation
                d2mx(1) = (t_rh2-t_xrh)
                dxm1(1) = (t_xrh-t_rh1)
                invd(1) = 1.0_r8/(t_rh2-t_rh1)
                d2mx(2) = (t_cat2-t_xct)
                dxm1(2) = (t_xct-t_cat1)
                invd(2) = 1.0_r8/(t_cat2-t_cat1)
                d2mx(3) = (t_fac2-t_xfac)
                dxm1(3) = (t_xfac-t_fac1)
                invd(3) = 1.0_r8/(t_fac2-t_fac1)
                d2mx(4) = (t_fbc2-t_xfbc)
                dxm1(4) = (t_xfbc-t_fbc1)
                invd(4) = 1.0_r8/(t_fbc2-t_fbc1)
                d2mx(5) = (t_faq2-t_xfaq)
                dxm1(5) = (t_xfaq-t_faq1)
                invd(5) = 1.0_r8/(t_faq2-t_faq1)


                do iv=1,38  ! variable number
                   opt5d(1,1,1,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=bep5to10(iv,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=bep5to10(iv,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=bep5to10(iv,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=bep5to10(iv,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, opt1, opt2)

                   ! finally, interpolation in the rh dimension
                   opt(iv) = ((t_rh2-t_xrh)*opt1+(t_xrh-t_rh1)*opt2) /(t_rh2-t_rh1)

                end do ! iv=1,38

                ! determine extinction coefficient
                call this%update(icol, k, kcomp, opt)

             end if ! Nnatk > 0
          end do ! icol
       end do ! k
    end do ! kcomp

  end subroutine intaeropt5to10

  ! ==========================================================
  subroutine zero(this, kcomp, ncol)

    class(extinction_coeffs_type) :: this
    integer , intent(in)    :: kcomp
    integer , intent(in)    :: ncol

    integer :: k
    integer :: icol

    ! initialize all output fields to zero
    do k=1,pver
       do icol=1,ncol
          this%bext440(icol,k,kcomp) = 0.0_r8
          this%babs440(icol,k,kcomp) = 0.0_r8
          this%bext500(icol,k,kcomp) = 0.0_r8
          this%babs500(icol,k,kcomp) = 0.0_r8
          this%bext550(icol,k,kcomp) = 0.0_r8
          this%babs550(icol,k,kcomp) = 0.0_r8
          this%bext670(icol,k,kcomp) = 0.0_r8
          this%babs670(icol,k,kcomp) = 0.0_r8
          this%bext870(icol,k,kcomp) = 0.0_r8
          this%babs870(icol,k,kcomp) = 0.0_r8
          this%bebg440(icol,k,kcomp) = 0.0_r8
          this%bebg500(icol,k,kcomp) = 0.0_r8
          this%bebg550(icol,k,kcomp) = 0.0_r8
          this%babg550(icol,k,kcomp) = 0.0_r8
          this%bebg670(icol,k,kcomp) = 0.0_r8
          this%bebg870(icol,k,kcomp) = 0.0_r8
          this%bebc440(icol,k,kcomp) = 0.0_r8
          this%bebc500(icol,k,kcomp) = 0.0_r8
          this%bebc550(icol,k,kcomp) = 0.0_r8
          this%babc550(icol,k,kcomp) = 0.0_r8
          this%bebc670(icol,k,kcomp) = 0.0_r8
          this%bebc870(icol,k,kcomp) = 0.0_r8
          this%beoc440(icol,k,kcomp) = 0.0_r8
          this%beoc500(icol,k,kcomp) = 0.0_r8
          this%beoc550(icol,k,kcomp) = 0.0_r8
          this%baoc550(icol,k,kcomp) = 0.0_r8
          this%beoc670(icol,k,kcomp) = 0.0_r8
          this%beoc870(icol,k,kcomp) = 0.0_r8
          this%besu440(icol,k,kcomp) = 0.0_r8
          this%besu500(icol,k,kcomp) = 0.0_r8
          this%besu550(icol,k,kcomp) = 0.0_r8
          this%basu550(icol,k,kcomp) = 0.0_r8
          this%besu670(icol,k,kcomp) = 0.0_r8
          this%besu870(icol,k,kcomp) = 0.0_r8
          this%bebg550lt1(icol,k,kcomp) = 0.0_r8
          this%bebg550gt1(icol,k,kcomp) = 0.0_r8
          this%bebc550lt1(icol,k,kcomp) = 0.0_r8
          this%bebc550gt1(icol,k,kcomp) = 0.0_r8
          this%beoc550lt1(icol,k,kcomp) = 0.0_r8
          this%beoc550gt1(icol,k,kcomp) = 0.0_r8
          this%besu550lt1(icol,k,kcomp) = 0.0_r8
          this%besu550gt1(icol,k,kcomp) = 0.0_r8
          this%backsc550(icol,k,kcomp) = 0.0_r8
       end do
    end do

  end subroutine zero

  ! ==========================================================
  subroutine update(this, icol, k, kcomp, opt)

    class(extinction_coeffs_type) :: this
    integer  , intent(in) :: icol
    integer  , intent(in) :: k
    integer  , intent(in) :: kcomp
    real(r8) , intent(in) :: opt(:)

    this%bext440(icol,k,kcomp)    = opt(1)
    this%bext500(icol,k,kcomp)    = opt(2)
    this%bext670(icol,k,kcomp)    = opt(3)
    this%bext870(icol,k,kcomp)    = opt(4)
    this%bebg440(icol,k,kcomp)    = opt(5)
    this%bebg500(icol,k,kcomp)    = opt(6)
    this%bebg670(icol,k,kcomp)    = opt(7)
    this%bebg870(icol,k,kcomp)    = opt(8)
    this%bebc440(icol,k,kcomp)    = opt(9)
    this%bebc500(icol,k,kcomp)    = opt(10)
    this%bebc670(icol,k,kcomp)    = opt(11)
    this%bebc870(icol,k,kcomp)    = opt(12)
    this%beoc440(icol,k,kcomp)    = opt(13)
    this%beoc500(icol,k,kcomp)    = opt(14)
    this%beoc670(icol,k,kcomp)    = opt(15)
    this%beoc870(icol,k,kcomp)    = opt(16)
    this%besu440(icol,k,kcomp)    = opt(17)
    this%besu500(icol,k,kcomp)    = opt(18)
    this%besu670(icol,k,kcomp)    = opt(19)
    this%besu870(icol,k,kcomp)    = opt(20)
    this%babs440(icol,k,kcomp)    = opt(21)
    this%babs500(icol,k,kcomp)    = opt(22)
    this%babs550(icol,k,kcomp)    = opt(23)
    this%babs670(icol,k,kcomp)    = opt(24)
    this%babs870(icol,k,kcomp)    = opt(25)
    this%bebg550lt1(icol,k,kcomp) = opt(26)
    this%bebg550gt1(icol,k,kcomp) = opt(27)
    this%bebc550lt1(icol,k,kcomp) = opt(28)
    this%bebc550gt1(icol,k,kcomp) = opt(29)
    this%beoc550lt1(icol,k,kcomp) = opt(30)
    this%beoc550gt1(icol,k,kcomp) = opt(31)
    this%besu550lt1(icol,k,kcomp) = opt(32)
    this%besu550gt1(icol,k,kcomp) = opt(33)
    this%backsc550(icol,k,kcomp)  = opt(34)
    this%babg550(icol,k,kcomp)    = opt(35)
    this%babc550(icol,k,kcomp)    = opt(36)
    this%baoc550(icol,k,kcomp)    = opt(37)
    this%basu550(icol,k,kcomp)    = opt(38)
    this%bebg550(icol,k,kcomp)    = opt(26)+opt(27)
    this%bebc550(icol,k,kcomp)    = opt(28)+opt(29)
    this%beoc550(icol,k,kcomp)    = opt(30)+opt(31)
    this%besu550(icol,k,kcomp)    = opt(32)+opt(33)
    this%bext550(icol,k,kcomp)    = this%bebg550(icol,k,kcomp) + this%bebc550(icol,k,kcomp) &
                                   +this%beoc550(icol,k,kcomp) + this%besu550(icol,k,kcomp)
  end subroutine update

  subroutine checkTableHeader (ifil)
    ! Read the header-text in a look-up table (in file with iu=ifil). 

    integer, intent(in) :: ifil
    character*80 :: headertext
    character*12 :: text0, text1 

    text0='X-CHECK LUT'
    text1='none       '
    do while (text1(2:12) .ne. text0(2:12))
       read(ifil,'(A)') headertext
       text1 = headertext(2:12) 
    enddo
  end subroutine checkTableHeader

#endif

end module oslo_aero_aerocom_opt
