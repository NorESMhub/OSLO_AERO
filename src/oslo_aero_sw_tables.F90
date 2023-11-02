module oslo_aero_sw_tables

  ! Purpose: To read in SW look-up tables for calculation of aerosol optical properties,
  ! and to interpolate between look-up table entries for SW optical aerosol properties.  

  ! Internal mixtures of process-tagged mass
  ! cate : total added mass (µg/m3 per particle per cm3) from condensation
  !        and wet phase chemistry/cloud processing, for kcomp = 1-2.
  !        cate should be scaled up/down whenever the modal parameters (modal
  !        radius and width) are increased/decreased a lot.
  ! cat  : total added mass (µg/m3 per particle per cm-3) from coagulation, condensation
  !        and wet phase chemistry/cloud processing, for kcomp = 5-10.
  !        cat should be scaled up/down whenever the modal parameters (modal
  !        radius and width) are increased/decreased a lot.
  ! fac  : mass fraction of cat or cate from coagulating carbonaceous aerosols (BC+OM).
  !        The remaining mass cate*(1-fac) or cat*(1-fac) is SO4.
  ! fbc  : mass fraction of BC from coagulating carbonaceous aerosols, BC/(BC+OM).
  ! faq  : mass fraction of sulfate which is produced in wet-phase, SO4aq/SO4.
  !        The remaining SO4 mass, SO4*(1-faq), is from condensation. 

  use shr_kind_mod            , only: r8 => shr_kind_r8
  use shr_sys_mod             , only: shr_sys_abort
  use ppgrid                  , only: pcols, pver
  use cam_logfile             , only: iulog
  use spmd_utils              , only: masterproc
  use oslo_aero_control       , only: oslo_aero_getopts, dir_string_length
  use oslo_aero_linear_interp , only: lininterpol3dim, lininterpol4dim, lininterpol5dim
  use oslo_aero_params        , only: nmodes, nbmodes, nbands, nlwbands
  use oslo_aero_const         , only: cate, cat, fac, faq, fbc, rh, fombg, fbcbg, rh, e

  implicit none
  private

  ! Interfaces
  public :: initopt
  public :: interpol0
  public :: interpol1
  public :: interpol2to3
  public :: interpol4
  public :: interpol5to10

  private :: initopt_lw

  ! ----------------------------
  ! Module variables set by table lookup
  ! ----------------------------

  real(r8), public :: om0(nbands)
  real(r8), public :: g0(nbands)
  real(r8), public :: be0(nbands)
  real(r8), public :: ke0(nbands)

  real(r8), public :: om1(nbands,10,6,16,6)
  real(r8), public :: g1 (nbands,10,6,16,6)
  real(r8), public :: be1(nbands,10,6,16,6)
  real(r8), public :: ke1(nbands,10,6,16,6)

  real(r8), public :: om2to3(nbands,10,16,6,2:3)
  real(r8), public :: g2to3 (nbands,10,16,6,2:3)
  real(r8), public :: be2to3(nbands,10,16,6,2:3)
  real(r8), public :: ke2to3(nbands,10,16,6,2:3)

  real(r8), public :: om4(nbands,10,6,16,6,6)
  real(r8), public :: g4 (nbands,10,6,16,6,6)
  real(r8), public :: be4(nbands,10,6,16,6,6)
  real(r8), public :: ke4(nbands,10,6,16,6,6)

  real(r8), public :: om5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: g5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: be5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: ke5to10(nbands,10,6,6,6,6,5:10)

  ! Long wave table info
  real(r8), public :: ka0(nlwbands)
  real(r8), public :: ka1(nlwbands,10,6,16,6)
  real(r8), public :: ka2to3(nlwbands,10,16,6,2:3)
  real(r8), public :: ka4(nlwbands,10,6,16,6,6)
  real(r8), public :: ka5to10(nlwbands,10,6,6,6,6,5:10)

!=============================================================================
contains
!=============================================================================

  subroutine initopt()

    ! Local variables
    integer  :: kcomp, iwl, irelh, ictot, ifac, ifbc, ifaq, i, irf
    integer  :: ifombg, ifbcbg
    integer  :: ik, ic, ifil, lin, linmax
    real(r8) :: catot, relh, frac, fabc, fraq, frombg, frbcbg
    real(r8) :: ssa, ass, ext, spext
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps3 = 1.e-3_r8
    real(r8) :: eps4 = 1.e-4_r8
    real(r8) :: eps6 = 1.e-6_r8
    character(len=dir_string_length) :: aerotab_table_dir
    !-----------------------------------------------------------

    call oslo_aero_getopts(aerotab_table_dir_out= aerotab_table_dir)

    ! Opening the 'kcomp'-files:
    open(50,file=trim(aerotab_table_dir)//'/kcomp0.out' ,form='formatted',status='old')
    open(40,file=trim(aerotab_table_dir)//'/kcomp1.out' ,form='formatted',status='old')
    open(41,file=trim(aerotab_table_dir)//'/kcomp2.out' ,form='formatted',status='old')
    open(42,file=trim(aerotab_table_dir)//'/kcomp3.out' ,form='formatted',status='old')
    open(43,file=trim(aerotab_table_dir)//'/kcomp4.out' ,form='formatted',status='old')
    open(44,file=trim(aerotab_table_dir)//'/kcomp5.out' ,form='formatted',status='old')
    open(45,file=trim(aerotab_table_dir)//'/kcomp6.out' ,form='formatted',status='old')
    open(46,file=trim(aerotab_table_dir)//'/kcomp7.out' ,form='formatted',status='old')
    open(47,file=trim(aerotab_table_dir)//'/kcomp8.out' ,form='formatted',status='old')
    open(48,file=trim(aerotab_table_dir)//'/kcomp9.out' ,form='formatted',status='old')
    open(49,file=trim(aerotab_table_dir)//'/kcomp10.out',form='formatted',status='old')

    ! Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 40,50
       call checkTableHeader (ifil)
    enddo

    !-------------------------------------------
    ! Mode 0, BC(ax)
    !-------------------------------------------

    ifil = 11
    linmax=nbands
    do lin = 1,linmax
       read(39+ifil,'(2I3,f8.3,4(x,e12.5))') kcomp, iwl, relh, ssa, ass, ext, spext
       om0(iwl)=ssa
       g0 (iwl)=ass
       be0(iwl)=ext    ! unit km^-1
       ke0(iwl)=spext  ! unit m^2/g
       ! write(iulog,*) 'kcomp, om =', kcomp, om0(iwl)
       ! write(iulog,*) 'kcomp, g  =', kcomp, g0(iwl)
       ! write(iulog,*) 'kcomp, be =', kcomp, be0(iwl)
       ! write(iulog,*) 'kcomp, ke =', kcomp, ke0(iwl)
    end do

    do iwl=1,nbands
       if(be0(iwl)<=0.0_r8) then
          write(iulog,*) 'be0 =', iwl, be0(iwl)
          call shr_sys_abort('initopt: error in initialization of be0')
       endif
    enddo

    if (masterproc) then
       write(iulog,*)'kcompN tables: mode 0 ok'
    end if

    !-------------------------------------------
    ! Mode 1 (H2SO4 and SOA + condesate from H2SO4 and SOA)
    !-------------------------------------------

    linmax = nbands*10*6*16*6   ! 14*10*6*16*6
    do lin = 1,linmax

       read(40,'(2I3,f8.3,3(x,e10.3),4(x,e12.5))') kcomp, iwl, relh, frombg, catot, frac, ssa, ass, ext, spext

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
          if(abs(frombg-fombg(ic))<eps4) then
             ifombg=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             exit
          endif
       end do

       om1(iwl,irelh,ifombg,ictot,ifac)=ssa
       g1 (iwl,irelh,ifombg,ictot,ifac)=ass
       be1(iwl,irelh,ifombg,ictot,ifac)=ext    ! unit km^-1
       ke1(iwl,irelh,ifombg,ictot,ifac)=spext  ! unit m^2/g

       ! write(iulog,*) 'kcomp, om =', kcomp, om1(iwl,irelh,ifombg,ictot,ifac)
       ! write(iulog,*) 'kcomp, g  =', kcomp, g1(iwl,irelh,ifombg,ictot,ifac)
       ! write(iulog,*) 'kcomp, be =', kcomp, be1(iwl,irelh,ifombg,ictot,ifac)
       ! write(iulog,*) 'kcomp, ke =', kcomp, ke1(iwl,irelh,ifombg,ictot,ifac)

    end do  ! lin

    kcomp=1
    do iwl=1,nbands
       do irelh=1,10
          do ifombg=1,6
             do ictot=1,16
                do ifac=1,6
                   if(be1(iwl,irelh,ifombg,ictot,ifac)<=0.0_r8) then
                      write(iulog,*) 'be1 =', iwl, irelh, ifombg, ictot, be1(iwl,irelh,ifombg,ictot,ifac)
                      call shr_sys_abort('initopt: error in initialization of be1')
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'kcompN tables: mode 1 ok'
    end if

    !-------------------------------------------
    ! Modes 2 to 3 (BC/OC + condensate from H2SO4 and SOA)
    !-------------------------------------------

    linmax=nbands*10*16*6
    do lin = 1,linmax

       read(41,'(2I3,f8.3,2(x,e10.3),4(x,e12.5))') kcomp, iwl, relh, catot, frac, ssa, ass, ext, spext

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

       om2to3(iwl,irelh,ictot,ifac,kcomp)=ssa
       g2to3 (iwl,irelh,ictot,ifac,kcomp)=ass
       be2to3(iwl,irelh,ictot,ifac,kcomp)=ext    ! unit km^-1
       ke2to3(iwl,irelh,ictot,ifac,kcomp)=spext  ! unit m^2/g

    end do  ! lin

    ! Prescribed dummy values for kcomp=3
    kcomp=3
    do iwl=1,nbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                om2to3(iwl,irelh,ictot,ifac,kcomp)=0.999_r8
                g2to3 (iwl,irelh,ictot,ifac,kcomp)=0.5_r8
                be2to3(iwl,irelh,ictot,ifac,kcomp)=0.0001_r8    ! unit km^-1
                ke2to3(iwl,irelh,ictot,ifac,kcomp)=1.0_r8       ! unit m^2/g
             enddo
          enddo
       enddo
    enddo

    do kcomp=2,3
       do iwl=1,nbands
          do irelh=1,10
             do ictot=1,16
                do ifac=1,6
                   if(be2to3(iwl,irelh,ictot,ifac,kcomp)<=0.0_r8) then
                      write(iulog,*) 'be2to3 =', iwl, irelh, ictot, ifac, be2to3(iwl,irelh,ictot,ifac,kcomp)
                      call shr_sys_abort('initopt: error in initialization of be2to3')
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    if (masterproc) then
       write(iulog,*)'kcompN tables: modes 2-3 ok'
    end if

    !-------------------------------------------
    ! Mode 4 (BC&OC + condensate from H2SO4 and SOA + wet phase (NH4)2SO4)
    !-------------------------------------------

    linmax = nbands*10*6*16*6*6
    do lin = 1,linmax
       read(43,'(2I3,f8.3,3(x,e10.3),f7.2,4(x,e12.5))') kcomp, iwl, relh, frbcbg, catot, frac, fraq, &
            ssa, ass, ext, spext

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
       do ic=1,6
          if(abs(fraq-faq(ic))<eps4) then
             ifaq=ic
             exit
          endif
       end do
       do ic=1,6
          if(abs((frbcbg-fbcbg(ic))/fbcbg(ic))<eps2) then
             ifbcbg=ic
             exit
          endif
       end do

       om4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ssa
       g4 (iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ass
       be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=ext    ! unit km^-1
       ke4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=spext  ! unit m^2/g

       ! write(iulog,*) 'kcomp, om =', kcomp, om4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
       ! write(iulog,*) 'kcomp, g  =', kcomp, g4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
       ! write(iulog,*) 'kcomp, be =', kcomp, be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
       ! write(iulog,*) 'kcomp, ke =', kcomp, ke4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)

    end do

    do iwl=1,nbands
       do irelh=1,10
          do ifbcbg=1,6
             do ictot=1,16
                do ifac=1,6
                   do ifaq=1,6
                      if(be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                         write(iulog,*) 'be4 =', iwl, irelh, ifbcbg, ictot, ifac, ifaq, &
                              be4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
                         call shr_sys_abort('initopt: error in initialization of be4')
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    if (masterproc) then
       write(iulog,*)'kcompN tables: mode 4 ok'
    end if

    !-------------------------------------------
    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !-------------------------------------------

    linmax = nbands*10*6*6*6*6     ! 14*10*6*6*6*6
    do ifil = 5,10
       do lin = 1,linmax

          read(39+ifil,'(2I3,f8.3,3(x,e10.3),f7.2,4(x,e12.5))') &
               kcomp, iwl, relh, catot, frac, fabc, fraq, ssa, ass, ext, spext

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
             endif
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

          om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ssa
          g5to10 (iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ass
          be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=ext    ! unit km^-1
          ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spext  ! unit m^2/g

          ! write(iulog,*) 'kcomp, om =', kcomp, om5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
          ! write(iulog,*) 'kcomp, g  =', kcomp, g5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
          ! write(iulog,*) 'kcomp, be =', kcomp, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
          ! write(iulog,*) 'kcomp, ke =', kcomp, ke5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)

       end do  ! ifil
    end do    ! lin

    do kcomp=5,10
       do iwl=1,nbands
          do irelh=1,10
             do ictot=1,6
                do ifac=1,6
                   do ifaq=1,6
                      if(be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                         write(iulog,*) 'be5to10 =', &
                              iwl, irelh, ictot, ifac, ifbc, ifaq, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                         write(iulog,*) 'kcomp, abs((fabc-fbc)/fbc) =', kcomp, abs((fabc-fbc(ic))/fbc(ic))
                         call shr_sys_abort('initopt: error in initialization of be5to10')
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    if (masterproc) then
       write(iulog,*)'kcompN tables: modes 5-10 ok'
    end if

    do ifil=40,50
       close (ifil)
    end do

    call initopt_lw()  ! initialize lw table info used in this module

  end subroutine initopt

  !=============================================================================

  subroutine initopt_lw

    !---------------------------------------------------------------
    ! Read in lwkcompN.out tables and set module variables
    ! om[0,1,2to3,4,5-10], g[0,1,2to3,4,5-10], be[0,1,2to3,4,5-10] and ke[0,1,2to3,4,5-10]
    !
    ! Modified for new aerosol schemes by Alf Kirkevaag in 2006. 
    ! Modified for new wavelength bands and look-up tables by 
    ! Alf Kirkevaag in 2014 and for SOA in 2015.
    !---------------------------------------------------------------

    ! Local variables
    integer  :: kcomp, iwl, irelh, ictot, ifac, ifbc, ifaq
    integer  :: ifombg, ifbcbg
    integer  :: ic, ifil, lin, linmax
    real(r8) :: catot, relh, frac, fabc, fraq, frombg, frbcbg
    real(r8) :: spabs
    real(r8) :: rh2(10)
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps3 = 1.e-3_r8
    real(r8) :: eps4 = 1.e-4_r8
    real(r8) :: eps6 = 1.e-6_r8
    real(r8) :: eps7 = 1.e-7_r8
    character(len=dir_string_length) :: aerotab_table_dir
    !------------------------------------------------------------------------

    call oslo_aero_getopts(aerotab_table_dir_out = aerotab_table_dir)

    open(40,file=trim(aerotab_table_dir)//'/lwkcomp1.out'  ,form="formatted",status="old")
    open(41,file=trim(aerotab_table_dir)//'/lwkcomp2.out'  ,form="formatted",status="old")
    open(42,file=trim(aerotab_table_dir)//'/lwkcomp3.out'  ,form="formatted",status="old")
    open(43,file=trim(aerotab_table_dir)//'/lwkcomp4.out'  ,form="formatted",status="old")
    open(44,file=trim(aerotab_table_dir)//'/lwkcomp5.out'  ,form="formatted",status="old")
    open(45,file=trim(aerotab_table_dir)//'/lwkcomp6.out'  ,form="formatted",status="old")
    open(46,file=trim(aerotab_table_dir)//'/lwkcomp7.out'  ,form="formatted",status="old")
    open(47,file=trim(aerotab_table_dir)//'/lwkcomp8.out'  ,form="formatted",status="old")
    open(48,file=trim(aerotab_table_dir)//'/lwkcomp9.out'  ,form="formatted",status="old")
    open(49,file=trim(aerotab_table_dir)//'/lwkcomp10.out' ,form="formatted",status="old")
    open(50,file=trim(aerotab_table_dir)//'/lwkcomp0.out'  ,form="formatted",status="old")

    !     Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 40,50
       call checkTableHeader (ifil)
    enddo

    !-------------------------------------------
    ! Mode 0, BC(ax)
    !-------------------------------------------

    ifil = 11
    linmax=nlwbands
    do lin = 1,linmax
       read(39+ifil,996) kcomp, iwl, relh, spabs
       ka0(iwl)=spabs  ! unit m^2/g
    end do

    do iwl=1,nlwbands
       if(ka0(iwl)<=0.0_r8) then
          write(iulog,*) 'ka0 =', iwl, ka0(iwl)
          call shr_sys_abort('initopt_lw: error in initialization of ka0')
       endif
    enddo
    if (masterproc) then
       write(iulog,*)'lwkcompN tables: lw mode 0 ok'
    end if

    !-------------------------------------------
    ! Mode 1 (H2SO4 + condesate from H2SO4 and SOA)
    !-------------------------------------------

    ifil = 1
    linmax=nlwbands*10*6*16*6
    do lin = 1,linmax

       read(39+ifil,997) kcomp, iwl, relh, frombg, catot, frac, spabs

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

       ka1(iwl,irelh,ifombg,ictot,ifac)=spabs  ! unit m^2/g

       ! write(*,*) 'kcomp, ka =', kcomp, ka1(iwl,irelh,ifombg,ictot,ifac)
       ! if(ifil==1) write(iulog,*) 'iwl,irelh,ifombg,ictot,ifac,ka =', &
       ! iwl,irelh,ictot,ifac,ka1(iwl,irelh,ifombg,ictot,ifac)

    end do  ! lin

    do iwl=1,nlwbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                if(ka1(iwl,irelh,ifombg,ictot,ifac)<=0.0_r8) then
                   write(iulog,*) 'ka1 =', iwl, irelh, ifombg, ictot, ifac, ka1(iwl,irelh,ifombg,ictot,ifac)
                   call shr_sys_abort('initopt_lw: error in initialization of ka1')
                endif
             enddo
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'lwkcompN tables: lw mode 1 ok'
    end if

    !-------------------------------------------
    ! Modes 2 to 3 (BC or OC + condensate from H2SO4 and SOA)
    !-------------------------------------------

    linmax = nlwbands*10*16*6
    do ifil = 2,2
       do lin = 1,linmax
          read(39+ifil,994) kcomp, iwl, relh, catot, frac, spabs

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

          ka2to3(iwl,irelh,ictot,ifac,kcomp)=spabs  ! unit m^2/g

       end do  ! lin
    end do    ! ifil

    ! Prescribed dummy values for kcomp=3
    kcomp=3
    do iwl=1,nlwbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                ka2to3(iwl,irelh,ictot,ifac,kcomp)=1.0_r8
             enddo
          enddo
       enddo
    enddo

    do kcomp=2,2
       do iwl=1,nlwbands
          do irelh=1,10
             do ictot=1,16
                do ifac=1,6
                   if(ka2to3(iwl,irelh,ictot,ifac,kcomp)<=0.0_r8) then
                      write(iulog,*) 'ka2to3 =', &
                           iwl, irelh, ictot, ifac, ka2to3(iwl,irelh,ictot,ifac,kcomp)
                      call shr_sys_abort('initopt_lw: error in initialization of ka2to3')
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'lwkcompN tables: lw mode 2-3 ok'
    end if

    !-------------------------------------------
    ! Mode 4 (BC&OC + condesate from H2SO4 and SOA + wetphase (NH4)2SO4)
    !-------------------------------------------

    ifil = 4
    linmax = nlwbands*10*6*16*6*6
    do lin = 1,linmax

       read(39+ifil,995) kcomp, iwl, relh, frbcbg, catot, frac, fraq, spabs

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

       ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=spabs  ! unit m^2/g
    end do

    do iwl=1,nlwbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                      write(iulog,*) 'ka4 =', &
                           iwl, irelh, ifbcbg, ictot, ifac, ifaq, ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
                      call shr_sys_abort('initopt_lw: error in initialization of ka4')
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'lwkcompN tables: lw mode 4 ok'
    end if

    !-------------------------------------------
    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !-------------------------------------------

    linmax = nlwbands*10*6*6*6*6
    do ifil = 5,10
       do lin = 1,linmax
          read(39+ifil,993) kcomp, iwl, relh, catot, frac, fabc, fraq, spabs

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
             endif
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

          ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spabs  ! unit m^2/g
       end do
    end do

    do kcomp=5,10
       do iwl=1,nlwbands
          do irelh=1,10
             do ictot=1,6
                do ifac=1,6
                   do ifaq=1,6
                      if(ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                         write(iulog,*) 'ka5to10 =', &
                              iwl, irelh, ictot, ifac, ifbc, ifaq, ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                         call shr_sys_abort('initopt_lw: error in initialization of ka5to10')
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
    if (masterproc) then
       write(iulog,*)'lwkcompN tables: lw mode 5-10 ok'
    end if

993 format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 5-10
994 format(2I3,f8.3,2(x,e10.3),x,e12.5)           ! 2-3
995 format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 4
996 format(2I3,f8.3,x,e12.5)                      ! 0
997 format(2I3,f8.3,3(x,e10.3),x,e12.5)           ! 1

    do ifil=40,50
       close (ifil)
    end do

  end subroutine initopt_lw

  !=============================================================================

  subroutine interpol0 (ncol, daylight, Nnatk, omega, gass, bex, ske, lw_on, kabs)
    !
    ! Arguments
    integer  , intent(in)  :: ncol                               ! number of atmospheric columns
    logical  , intent(in)  :: daylight(pcols)                    ! calculations also at (polar) night if daylight=.true.
    logical  , intent(in)  :: lw_on                              ! LW calculations are performed if true
    real(r8) , intent(in)  :: Nnatk(pcols,pver,0:nmodes)         ! modal aerosol number concentration
    real(r8) , intent(out) :: omega(pcols,pver,0:nmodes,nbands)  ! spectral modal single scattering albedo
    real(r8) , intent(out) :: gass(pcols,pver,0:nmodes,nbands)   ! spectral modal asymmetry factor
    real(r8) , intent(out) :: bex(pcols,pver,0:nmodes,nbands)    ! spectral modal extinction coefficient
    real(r8) , intent(out) :: ske(pcols,pver,0:nmodes,nbands)    ! spectral modal specific extinction coefficient
    real(r8) , intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands) ! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol
    !---------------------------------------

    kcomp=0
    do i=1,nbands
       do icol=1,ncol
          do k=1,pver
             omega(icol,k,kcomp,i)=0.0_r8
             gass(icol,k,kcomp,i)=0.0_r8
             bex(icol,k,kcomp,i)=0.0_r8
             ske(icol,k,kcomp,i)=0.0_r8
          end do
       end do
    end do
    do i=1,nlwbands
       do icol=1,ncol
          do k=1,pver
             kabs(icol,k,kcomp,i)=0.0_r8
          end do
       end do
    end do

    ! SW optical parameters

    do k=1,pver
       do icol=1,ncol
          ! if(Nnatk(icol,k,kcomp)>0.0_r8) then
          if(daylight(icol)) then
             do i=1,nbands   ! i = wavelength index
                omega(icol,k,kcomp,i)=om0(i)
                gass(icol,k,kcomp,i)=g0(i)
                bex(icol,k,kcomp,i)=be0(i)
                ske(icol,k,kcomp,i)=ke0(i)
             end do          ! i
          else  ! daylight
             ! Need be and ke in   nband=4 for lw calculation
             bex(icol,k,kcomp,4)=be0(4)
             ske(icol,k,kcomp,4)=ke0(4)
          end if ! daylight
       end do ! icol
    end do ! k

    ! LW optical parameters

    if(lw_on) then
       do k=1,pver
          do icol=1,ncol
             do i=1,nlwbands   ! i = wavelength index
                kabs(icol,k,kcomp,i)=ka0(i)
             end do            ! i
          end do ! icol
       end do ! k

    endif ! lw_on

  end subroutine interpol0

  !=============================================================================

  subroutine interpol1 (ncol, daylight, xrh, irh1, mplus10, Nnatk, xfombg, ifombg1, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    !
    ! Arguments
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xfombg(pcols,pver)         ! SOA/(SOA+H2SO4) for the background mode
    integer,  intent(in) :: ifombg1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)

    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absoption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol, kc10
    real(r8) a, b
    integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifo1, t_ifo2
    real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2
    real(r8) t_cat1, t_cat2, t_fombg1, t_fombg2, t_xfombg
    real(r8) d2mx(4), dxm1(4), invd(4)
    real(r8) opt4d(2,2,2,2)
    real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    do kcomp=1,1

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifo1 = ifombg1(icol,k)
             t_ifo2 = t_ifo1+1

             t_rh1  = rh(t_irh1)
             !x      t_rh2  = t_rh1+1
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_fombg1 = fombg(t_ifo1)
             t_fombg2 = fombg(t_ifo2)

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kcomp)
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


             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   ! single scattering albedo:

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=om1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=om1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=om1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=om1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=om1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=om1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=om1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=om1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) /(t_rh2-t_rh1)
                   !alt       omega(icol,k,kc10,i)=(d2mx(1)*ome1+dxm1(1)*ome2)*invd(1)

                   ! asymmetry factor

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=g1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=g1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=g1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=g1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=g1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=g1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=g1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=g1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, ge1, ge2)

                   ! finally, interpolation in the rh dimension (dim. 1)
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) /(t_rh2-t_rh1)
                   !alt      gass(icol,k,kc10,i)=(d2mx(1)*ge1+dxm1(1)*ge2)*invd(1)

                   ! aerosol extinction

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             else  ! daylight

                ! aerosol extinction used for size information in LW
                i=4

                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                opt4d(1,1,1,2)=be1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                opt4d(1,1,2,1)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                opt4d(1,1,2,2)=be1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                opt4d(1,2,1,1)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                opt4d(1,2,1,2)=be1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                opt4d(1,2,2,1)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                opt4d(1,2,2,2)=be1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                opt4d(2,1,1,1)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                opt4d(2,1,1,2)=be1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                opt4d(2,1,2,1)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                opt4d(2,1,2,2)=be1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                opt4d(2,2,1,1)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                opt4d(2,2,1,2)=be1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                opt4d(2,2,2,1)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                opt4d(2,2,2,2)=be1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac, cat and fombg dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                ! aerosol specific extinction

                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                opt4d(1,1,1,2)=ke1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                opt4d(1,1,2,1)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                opt4d(1,1,2,2)=ke1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                opt4d(1,2,1,1)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                opt4d(1,2,1,2)=ke1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                opt4d(1,2,2,1)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                opt4d(1,2,2,2)=ke1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                opt4d(2,1,1,1)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                opt4d(2,1,1,2)=ke1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                opt4d(2,1,2,1)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                opt4d(2,1,2,2)=ke1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                opt4d(2,2,1,1)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                opt4d(2,2,1,2)=ke1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                opt4d(2,2,2,1)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                opt4d(2,2,2,2)=ke1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac, cat and fombg dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                   !alt        ske(icol,k,kc10,i)=(d2mx(1)*ske1+dxm1(1)*ske2)*invd(1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                   !alt        a=(log(ske2)-log(ske1))*invd(1)
                   !alt        b=(t_rh2*log(ske1)-t_rh1*log(ske2))*invd(1)
                   !alt        ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i

             if (lw_on) then

                ! LW optical parameters
                do i=1,nlwbands            ! i = wavelength index

                   ! aerosol specific absorption in LW

                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc1)
                   opt4d(1,1,1,2)=ka1(i,t_irh1,t_ifo1,t_ict1,t_ifc2)
                   opt4d(1,1,2,1)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc1)
                   opt4d(1,1,2,2)=ka1(i,t_irh1,t_ifo1,t_ict2,t_ifc2)
                   opt4d(1,2,1,1)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc1)
                   opt4d(1,2,1,2)=ka1(i,t_irh1,t_ifo2,t_ict1,t_ifc2)
                   opt4d(1,2,2,1)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc1)
                   opt4d(1,2,2,2)=ka1(i,t_irh1,t_ifo2,t_ict2,t_ifc2)
                   opt4d(2,1,1,1)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc1)
                   opt4d(2,1,1,2)=ka1(i,t_irh2,t_ifo1,t_ict1,t_ifc2)
                   opt4d(2,1,2,1)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc1)
                   opt4d(2,1,2,2)=ka1(i,t_irh2,t_ifo1,t_ict2,t_ifc2)
                   opt4d(2,2,1,1)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc1)
                   opt4d(2,2,1,2)=ka1(i,t_irh2,t_ifo2,t_ict1,t_ifc2)
                   opt4d(2,2,2,1)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc1)
                   opt4d(2,2,2,2)=ka1(i,t_irh2,t_ifo2,t_ict2,t_ifc2)

                   ! interpolation in the fac, cat and fombg dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30)
                   kabs2=max(kabs2,1.e-30)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k
       ! write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
       ! write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
       ! write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
       ! write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)
    end do  ! kcomp

  end subroutine interpol1

  !=============================================================================

  subroutine interpol2to3 (ncol, daylight, xrh, irh1, mplus10, Nnatk, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    ! Arguments
    integer,  intent(in)  :: ncol                       ! number of atmospheric columns
    integer,  intent(in)  :: mplus10                    ! mode number (0) or number + 10 (1)
    logical,  intent(in)  :: daylight(pcols)            ! only daylight calculations if .true.
    logical,  intent(in)  :: lw_on                      ! LW calculations are performed if true
    real(r8), intent(in)  :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in)  :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in)  :: irh1(pcols,pver)
    real(r8), intent(in)  :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in)  :: ict1(pcols,pver,nmodes)
    real(r8), intent(in)  :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in)  :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer  :: i, kcomp, k, icol, kc10
    real(r8) :: a, b
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2,t_cat1, t_cat2
    real(r8) :: d2mx(3), dxm1(3), invd(3)
    real(r8) :: opt3d(2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: kabs1, kabs2
    !---------------------------------------

    do kcomp=2,2

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing
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

             ! partial lengths along each dimension (1-4) for interpolation
             d2mx(1) = (t_rh2-t_xrh)
             dxm1(1) = (t_xrh-t_rh1)
             invd(1) = 1.0_r8/(t_rh2-t_rh1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)

             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands            ! i = wavelength index

                   ! single scattering albedo:

                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=om2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=om2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=om2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=om2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=om2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=om2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=om2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=om2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2)/(t_rh2-t_rh1)

                   ! asymmetry factor
                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=g2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=g2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=g2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=g2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=g2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=g2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=g2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=g2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2)/(t_rh2-t_rh1)

                   ! aerosol extinction
                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=be2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=be2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=be2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=be2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=be2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=be2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=be2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=be2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, bex1, bex2)

                   bex1=max(bex1,1.e-30)
                   bex2=max(bex2,1.e-30)

                   ! finally, interpolation in the rh dimension
                   if(t_xrh <= 0.37) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2)/(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight

                ! aerosol extinction used for LW size information
                ! end points as basis for multidimentional linear interpolation
                i=4
                opt3d(1,1,1)=be2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=be2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=be2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=be2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=be2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=be2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=be2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=be2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                ! interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, bex1, bex2)

                bex1=max(bex1,1.e-30)
                bex2=max(bex2,1.e-30)

                ! finally, interpolation in the rh dimension
                if(t_xrh <= 0.37) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2)/(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                ! aerosol specific extinction
                ! end points as basis for multidimentional linear interpolation
                opt3d(1,1,1)=ke2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                opt3d(1,1,2)=ke2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                opt3d(1,2,1)=ke2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                opt3d(1,2,2)=ke2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                opt3d(2,1,1)=ke2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                opt3d(2,1,2)=ke2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                opt3d(2,2,1)=ke2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                opt3d(2,2,2)=ke2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                ! interpolation in the (fac and) cat dimension
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, ske1, ske2)

                ske1=max(ske1,1.e-30)
                ske2=max(ske2,1.e-30)

                ! finally, interpolation in the rh dimension
                if(t_xrh <= 0.37) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2)/(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i

             ! LW optical parameters
             if (lw_on) then

                do i=1,nlwbands            ! i = wavelength index

                   ! aerosol specific absorption in LW
                   ! end points as basis for multidimentional linear interpolation
                   opt3d(1,1,1)=ka2to3(i,t_irh1,t_ict1,t_ifc1,kcomp)
                   opt3d(1,1,2)=ka2to3(i,t_irh1,t_ict1,t_ifc2,kcomp)
                   opt3d(1,2,1)=ka2to3(i,t_irh1,t_ict2,t_ifc1,kcomp)
                   opt3d(1,2,2)=ka2to3(i,t_irh1,t_ict2,t_ifc2,kcomp)
                   opt3d(2,1,1)=ka2to3(i,t_irh2,t_ict1,t_ifc1,kcomp)
                   opt3d(2,1,2)=ka2to3(i,t_irh2,t_ict1,t_ifc2,kcomp)
                   opt3d(2,2,1)=ka2to3(i,t_irh2,t_ict2,t_ifc1,kcomp)
                   opt3d(2,2,2)=ka2to3(i,t_irh2,t_ict2,t_ifc2,kcomp)

                   ! interpolation in the (fac and) cat dimension
                   call lininterpol3dim (d2mx, dxm1, invd, opt3d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2)/(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             endif ! lw_on
          end do ! icol
       end do ! k
       ! write(*,*) 'kcomp, omega(1,26,kcomp,4)=', kcomp, omega(1,26,kcomp,4)
       ! write(*,*) 'kcomp, gass(1,26,kcomp,4)=', kcomp, gass(1,26,kcomp,4)
       ! write(*,*) 'kcomp, bex(1,26,kcomp,4)=', kcomp, bex(1,26,kcomp,4)
       ! write(*,*) 'kcomp, ske(1,26,kcomp,4)=', kcomp, ske(1,26,kcomp,4)
    end do  ! kcomp

  end subroutine interpol2to3

  !=============================================================================

  subroutine interpol4 (ncol, daylight, xrh, irh1, mplus10, Nnatk, xfbcbg, ifbcbg1, &
       xct, ict1, xfac, ifac1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Arguments
    integer, intent(in)   :: ncol                       ! number of atmospheric columns
    integer, intent(in)   :: mplus10                    ! mode number (0) or number + 10 (1)
    logical, intent(in)   :: daylight(pcols)            ! only daylight calculations if .true.
    logical, intent(in)   :: lw_on                      ! LW calculations are performed if true
    real(r8), intent(in)  :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in)  :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in)  :: irh1(pcols,pver)
    real(r8), intent(in)  :: xfbcbg(pcols,pver)         ! mass fraction BC/(BC+OC) for the background mode
    integer,  intent(in)  :: ifbcbg1(pcols,pver)
    real(r8), intent(in)  :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in)  :: ict1(pcols,pver,nmodes)
    real(r8), intent(in)  :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in)  :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in)  :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer,  intent(in)  :: ifaq1(pcols,pver,nbmodes)
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer  :: i, kcomp, k, kc10, icol
    real(r8) :: a, b
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2, t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq, t_fbcbg1, t_fbcbg2, t_xfbcbg, t_fac1
    real(r8) :: t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: kabs1, kabs2
    !---------------------------------------

    do kcomp=4,4

       if(mplus10==0) then
          kc10=kcomp
       else
          kc10=kcomp+10
       endif

       ! write(*,*) 'Before init-loop', kc10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kc10,i)=0.0_r8
                gass(icol,k,kc10,i)=0.0_r8
                bex(icol,k,kc10,i)=0.0_r8
                ske(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kc10,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing

             t_irh1 = irh1(icol,k)
             t_irh2 = t_irh1+1
             t_ict1 = ict1(icol,k,kc10)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifb1 = ifbcbg1(icol,k)
             t_ifb2 = t_ifb1+1
             t_ifa1 = ifaq1(icol,k,kcomp)
             t_ifa2 = t_ifa1+1

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_fbcbg1 = fbcbg(t_ifb1)
             t_fbcbg2 = fbcbg(t_ifb2)
             t_faq1 = faq(t_ifa1)
             t_faq2 = faq(t_ifa2)

             t_xrh  = xrh(icol,k)
             t_xct  = xct(icol,k,kc10)
             t_xfac = xfac(icol,k,kcomp)
             t_xfbcbg = xfbcbg(icol,k)
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

             ! SW optical parameters
             if(daylight(icol)) then

                do i=1,nbands ! wavelength index

                   ! single scattering albedo:
                   opt5d(1,1,1,1,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=om4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=om4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=om4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=om4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=om4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=om4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=om4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=om4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) /(t_rh2-t_rh1)

                   ! asymmetry factor
                   opt5d(1,1,1,1,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=g4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=g4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=g4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=g4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=g4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=g4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=g4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=g4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2)/(t_rh2-t_rh1)

                   ! aerosol extinction
                   opt5d(1,1,1,1,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             else  ! daylight

                ! aerosol extinction called for use in size estimate for use in LW
                i=4
                opt5d(1,1,1,1,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=be4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=be4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=be4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=be4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=be4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=be4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=be4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=be4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2)/(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                ! aerosol specific extinction
                opt5d(1,1,1,1,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,1,1,1,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,1,1,2,1)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,1,1,2,2)=ke4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,1,2,1,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,1,2,1,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,1,2,2,1)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,1,2,2,2)=ke4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(1,2,1,1,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(1,2,1,1,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(1,2,1,2,1)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(1,2,1,2,2)=ke4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(1,2,2,1,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(1,2,2,1,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(1,2,2,2,1)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(1,2,2,2,2)=ke4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,1,1,1,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,1,1,1,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,1,1,2,1)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,1,1,2,2)=ke4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,1,2,1,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,1,2,1,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,1,2,2,1)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,1,2,2,2)=ke4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt5d(2,2,1,1,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt5d(2,2,1,1,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt5d(2,2,1,2,1)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt5d(2,2,1,2,2)=ke4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt5d(2,2,2,1,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt5d(2,2,2,1,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt5d(2,2,2,2,1)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt5d(2,2,2,2,2)=ke4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac, cat and fbcbg dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2)/(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif
             end do ! i

             ! LW optical parameters
             if (lw_on) then
                do i=1,nlwbands            ! i = wavelength index

                   ! aerosol specific absorption
                   opt5d(1,1,1,1,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,1,1,1,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,1,1,2,1)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,1,1,2,2)=ka4(i,t_irh1,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,1,2,1,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,1,2,1,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,1,2,2,1)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,1,2,2,2)=ka4(i,t_irh1,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(1,2,1,1,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(1,2,1,1,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(1,2,1,2,1)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(1,2,1,2,2)=ka4(i,t_irh1,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(1,2,2,1,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(1,2,2,1,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(1,2,2,2,1)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(1,2,2,2,2)=ka4(i,t_irh1,t_ifb2,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,1,1,1,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,1,1,1,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,1,1,2,1)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,1,1,2,2)=ka4(i,t_irh2,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,1,2,1,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,1,2,1,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,1,2,2,1)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,1,2,2,2)=ka4(i,t_irh2,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                   opt5d(2,2,1,1,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                   opt5d(2,2,1,1,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                   opt5d(2,2,1,2,1)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                   opt5d(2,2,1,2,2)=ka4(i,t_irh2,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                   opt5d(2,2,2,1,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                   opt5d(2,2,2,1,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                   opt5d(2,2,2,2,1)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                   opt5d(2,2,2,2,2)=ka4(i,t_irh2,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                   ! interpolation in the faq, fac, cat and fbcbg dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kc10,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2)/(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             endif ! lw_on
          end do ! icol
       end do ! k
       ! write(*,*) 'kcomp, omega(1,26,kc10,4)=', kcomp, omega(1,26,kc10,4)
       ! write(*,*) 'kcomp, gass(1,26,kc10,4)=', kcomp, gass(1,26,kc10,4)
       ! write(*,*) 'kcomp, bex(1,26,kc10,4)=', kcomp, bex(1,26,kc10,4)
       ! write(*,*) 'kcomp, ske(1,26,kc10,4)=', kcomp, ske(1,26,kc10,4)
    end do  ! kcomp

  end subroutine interpol4

  !=============================================================================

  subroutine interpol5to10 (ncol, daylight, xrh, irh1, Nnatk, xct, ict1, &
       xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! modal (OC+BC)/(SO4+BC+OC)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfbc(pcols,pver,nbmodes)   ! modal BC/(OC+BC)
    integer,  intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Output arguments
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient

    ! Local variables
    integer  :: i, kcomp, k, icol
    real(r8) :: a, b
    integer  :: t_irh1, t_irh2, t_ict1, t_ict2, t_ifa1, t_ifa2
    integer  :: t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq, t_fbc1, t_fbc2, t_xfbc, t_fac1
    real(r8) :: t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2, t_cat1, t_cat2
    real(r8) :: d2mx(5), dxm1(5), invd(5)
    real(r8) :: opt5d(2,2,2,2,2)
    real(r8) :: ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) :: kabs1, kabs2
    !---------------------------------------

    do kcomp=5,10
       do i=1,nbands
          do icol=1,ncol
             do k=1,pver
                omega(icol,k,kcomp,i)=0.0_r8
                gass(icol,k,kcomp,i)=0.0_r8
                bex(icol,k,kcomp,i)=0.0_r8
                ske(icol,k,kcomp,i)=0.0_r8
             end do
          end do
       end do
       do i=1,nlwbands
          do icol=1,ncol
             do k=1,pver
                kabs(icol,k,kcomp,i)=0.0_r8
             end do
          end do
       end do

       do k=1,pver
          do icol=1,ncol

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

             ! SW optical parameters
             if(daylight(icol)) then
                do i=1,nbands            ! i = wavelength index

                   ! single scattering albedo:
                   opt5d(1,1,1,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=om5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=om5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=om5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=om5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ome1, ome2)

                   ! finally, interpolation in the rh dimension
                   omega(icol,k,kcomp,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2)/(t_rh2-t_rh1)

                   ! asymmetry factor
                   opt5d(1,1,1,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=g5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=g5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=g5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=g5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, ge1, ge2)

                   ! finally, interpolation in the rh dimension
                   gass(icol,k,kcomp,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2)/(t_rh2-t_rh1)

                   ! aerosol extinction
                   opt5d(1,1,1,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                   bex1=max(bex1,1.e-30_r8)
                   bex2=max(bex2,1.e-30_r8)

                   ! finally, interpolation in the rh dimension
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1 + (t_xrh-t_rh1)*bex2)/(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             else  ! daylight

                ! aerosol extinction  used for aerosol size estimate needed for LW calculations
                i=4
                opt5d(1,1,1,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,1,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,1,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,1,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,1,2,1,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,2,1,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,2,2,1)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,2,2,2)=be5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,1,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,1,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,1,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,1,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,2,1,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,2,1,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,2,2,1)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,2,2,2)=be5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,1,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,1,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,1,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,1,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,2,1,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,2,1,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,2,2,1)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,2,2,2)=be5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,1,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,1,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,1,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,1,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,2,1,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,2,1,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,2,2,1)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,2,2,2)=be5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                ! interpolation in the faq, fbc, fac and cat dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, bex1, bex2)

                bex1=max(bex1,1.e-30_r8)
                bex2=max(bex2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37_r8) then
                   bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                ! aerosol specific extinction
                opt5d(1,1,1,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,1,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,1,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,1,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,1,2,1,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,1,2,1,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,1,2,2,1)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,1,2,2,2)=ke5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,1,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,1,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,1,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,1,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(1,2,2,1,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(1,2,2,1,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(1,2,2,2,1)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(1,2,2,2,2)=ke5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,1,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,1,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,1,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,1,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,1,2,1,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,1,2,1,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,1,2,2,1)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,1,2,2,2)=ke5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,1,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,1,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,1,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,1,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                opt5d(2,2,2,1,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                opt5d(2,2,2,1,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                opt5d(2,2,2,2,1)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                opt5d(2,2,2,2,2)=ke5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                ! interpolation in the faq, fbc, fac and cat dimensions
                call lininterpol5dim (d2mx, dxm1, invd, opt5d, ske1, ske2)

                ske1=max(ske1,1.e-30_r8)
                ske2=max(ske2,1.e-30_r8)

                ! finally, interpolation in the rh dimension
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37_r8) then
                   ske(icol,k,kcomp,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kcomp,i)=e**(a*t_xrh+b)
                endif
             end do ! i

             ! LW optical parameters
             if (lw_on) then
                do i=1,nlwbands            ! i = wavelength index
                   ! aerosol specific absorption
                   opt5d(1,1,1,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,1,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,1,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,1,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,1,2,1,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,1,2,1,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,1,2,2,1)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,1,2,2,2)=ka5to10(i,t_irh1,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,1,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,1,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,1,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,1,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(1,2,2,1,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(1,2,2,1,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(1,2,2,2,1)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(1,2,2,2,2)=ka5to10(i,t_irh1,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,1,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,1,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,1,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,1,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,1,2,1,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,1,2,1,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,1,2,2,1)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,1,2,2,2)=ka5to10(i,t_irh2,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,1,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,1,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,1,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,1,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt5d(2,2,2,1,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt5d(2,2,2,1,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt5d(2,2,2,2,1)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt5d(2,2,2,2,2)=ka5to10(i,t_irh2,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, fac and cat dimensions
                   call lininterpol5dim (d2mx, dxm1, invd, opt5d, kabs1, kabs2)

                   kabs1=max(kabs1,1.e-30_r8)
                   kabs2=max(kabs2,1.e-30_r8)

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
                      kabs(icol,k,kcomp,i)=((t_rh2-t_xrh)*kabs1+(t_xrh-t_rh1)*kabs2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(kabs2)-log(kabs1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(kabs1)-t_rh1*log(kabs2))/(t_rh2-t_rh1)
                      kabs(icol,k,kcomp,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

             endif ! lw_on

          end do ! icol
       end do ! k
    end do  ! kcomp

  end subroutine interpol5to10

  !=============================================================================

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

end module oslo_aero_sw_tables
