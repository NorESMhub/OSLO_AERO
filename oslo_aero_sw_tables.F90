module oslo_aero_sw_tables

  ! Purpose: To read in SW look-up tables for calculation of aerosol optical properties,
  ! and to define the grid for discrete input-values in these look-up tables.

  ! Purpose: To interpolate between look-up table entries for SW optical aerosol properties.
  ! Optimized for speed by Arild Burud and Egil Storen (NoSerC), June-July 2002
  ! Updated for new kcomp1.out including condensed SOA - Alf Kirkevaag, May 2013.
  ! Extended for new SOA treatment for  kcomp1-4.out and treating SOA as coagulated OC
  ! for kcomp5-10 - Alf Kirkevaag, August 2015, and also rewritten to a more generalized
  ! for for interpolations using common subroutines interpol*dim.

  ! Modified for new wavelength bands and look-up tables - Alf Kirkevaag Dec. 2013.
  ! Updated for reading input files with extra header info - Alf Kirkevaag, May 2015.
  ! Extended for new SOA treatment - Alf Kirkevaag, August 2015.
  ! Added output (ASCII) Jabuary 2016: #ifdef COLTST4INTCONS -> extinction
  ! koefficients (wrt. all added mass including condensed water vapour) are
  ! written out for checking against the look-up tables (using xmgrace), e.g.
  ! as function of RH (to be changed to whatever parameter the user is interested in)
  ! Modified for optimized added masses and mass fractions for concentrations from
  ! condensation, coagulation or cloud-processing - Alf Kirkevaag, May 2016.
  ! Modified cate values for kcomp=2 (as  in AeroTab) - Alf Kirkevaag October 2016.

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
  use ppgrid                  , only: pcols, pver
  use cam_logfile             , only: iulog
  !
  use oslo_aero_control       , only: oslo_aero_getopts, dir_string_length
  use oslo_aero_linear_interp , only: lininterpol3dim, lininterpol4dim, lininterpol5dim
  use oslo_aero_params        , only: nmodes, nbmodes

  implicit none
  private

  ! Interfaces
  public :: initopt
  public :: initopt_lw
  public :: inputForInterpol
  public :: interpol0
  public :: interpol1
  public :: interpol2to3
  public :: interpol4
  public :: interpol5to10

  integer, public, parameter :: nbands=14    ! number of aerosol spectral bands in SW
  integer, public, parameter :: nbmp1=11     ! number of first non-background mode

  real(r8), public, dimension(10)     :: rh
  real(r8), public, dimension(6)      :: fombg, fbcbg, fac, fbc, faq
  real(r8), public, dimension(4,16)   :: cate
  real(r8), public, dimension(5:10,6) :: cat

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

  real(r8), public :: om0(nbands)
  real(r8), public :: g0(nbands)
  real(r8), public :: be0(nbands)
  real(r8), public :: ke0(nbands)

  real(r8), public :: om5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: g5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: be5to10(nbands,10,6,6,6,6,5:10)
  real(r8), public :: ke5to10(nbands,10,6,6,6,6,5:10)

  ! relative humidity (RH, as integer for output variable names) for use in AeroCom code
  integer, public, dimension(6) :: RF = (/0, 40, 55, 65, 75, 85 /)

  ! AeroCom specific RH input variables for use in opticsAtConstRh.F90
  integer , public :: irhrf1(6)
  real(r8), public :: xrhrf(6)

  real(r8), public :: e, eps
  parameter (e=2.718281828_r8, eps=1.0e-30_r8)

  ! Array bounds in the tabulated optical parameters
  integer, public, parameter :: nlwbands=16    ! number of aerosol spectral bands in LW

  real(r8), public :: ka0(nlwbands)
  real(r8), public :: ka1(nlwbands,10,6,16,6)
  real(r8), public :: ka2to3(nlwbands,10,16,6,2:3)
  real(r8), public :: ka4(nlwbands,10,6,16,6,6)
  real(r8), public :: ka5to10(nlwbands,10,6,6,6,6,5:10)

contains

  subroutine initopt()

    !---------------------------------------------------------------
    ! Modified for new aerosol schemes by Alf Kirkevaag in January
    ! 2006. Modified for new wavelength bands and look-up tables
    ! by Alf Kirkevaag in December 2013, and for SOA in August 2015.
    !---------------------------------------------------------------

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

    ! Defining array bounds for tabulated optical parameters (and r and sigma)
    ! relative humidity (only 0 value used for r and sigma tables):
    rh = (/ 0.0_r8, 0.37_r8, 0.47_r8, 0.65_r8, 0.75_r8, 0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8, 0.995_r8 /)

    ! AeroCom specific RH input variables for use in opticsAtConstRh.F90
    do irf=1,6
       xrhrf(irf)  = real(RF(irf))*0.01_r8
    enddo
    do irelh=1,9
       do irf=1,6
          if(xrhrf(irf)>=rh(irelh).and.xrhrf(irf)<=rh(irelh+1)) then
             irhrf1(irf)=irelh
          endif
       end do
    end do

    ! mass fractions internal mixtures in background (fombg and fbcbg) and mass added to the
    ! background modes (fac, faq, faq)
    fombg = (/ 0.0_r8, 0.2_r8,  0.4_r8, 0.6_r8, 0.8_r8, 1.0_r8  /)
    fac =   (/ 0.0_r8, 0.2_r8,  0.4_r8, 0.6_r8, 0.8_r8, 1.0_r8  /)
    faq =   (/ 0.0_r8, 0.2_r8,  0.4_r8, 0.6_r8, 0.8_r8, 1.0_r8  /)

    ! with more weight on low fractions (thus a logaritmic f axis) for BC,
    ! which is less ambundant than sulfate and OC, and the first value
    ! corresponding to a clean background mode:
    fbcbg(1)=1.e-10_r8
    fbc(1)=1.e-10_r8
    do i=2,6
       fbcbg(i)=10**((i-1)/4.0_r8-1.25_r8)
       fbc(i)=fbcbg(i)
    end do
    ! and most weight on small concentrations for added mass onto the background:
    do kcomp=1,4
       cate(kcomp,1)=1.e-10_r8
       do i=2,16
          if(kcomp.eq.1.or.kcomp.eq.2) then
             cate(kcomp,i)=10.0_r8**((i-1)/3.0_r8-6.222_r8)
          elseif(kcomp.eq.3) then
             cate(kcomp,i)=1.0e-10_r8  ! not used
          else
             cate(kcomp,i)=10.0_r8**((i-1)/3.0_r8-4.301_r8)
          endif
       end do
    end do
    do kcomp=5,10
       cat(kcomp,1) =1.e-10_r8
       do i=2,6
          if(kcomp.eq.5) then
             cat(kcomp,i)=10.0_r8**((i-1)-3.824_r8)
          elseif(kcomp.eq.6) then
             cat(kcomp,i)=10.0_r8**((i-1)-3.523_r8)
          elseif(kcomp.eq.7) then
             cat(kcomp,i)=10.0_r8**((i-1)-3.699_r8)
          elseif(kcomp.eq.8) then
             cat(kcomp,i)=10.0_r8**((i-1)-4.921_r8)
          elseif(kcomp.eq.9) then
             cat(kcomp,i)=10.0_r8**((i-1)-3.301_r8)
          else
             cat(kcomp,i)=10.0_r8**((i-1)-3.699_r8)
          endif
       end do
    end do

    call oslo_aero_getopts(aerotab_table_dir_out= aerotab_table_dir)

    ! Opening the 'kcomp'-files:

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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

    ! Then reading in the look-up table entries for each file (kcomp*.out)

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    ! Mode 0, BC(ax)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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
          write(iulog,*) 'Error in initialization of be0'
          stop
       endif
    enddo

    write(iulog,*)'mode 0 ok'


    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    ! Mode 1 (H2SO4 and SOA + condesate from H2SO4 and SOA)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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
                      write(iulog,*) 'Error in initialization of be1'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'mode 1 ok'

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    ! Modes 2 to 3 (BC/OC + condensate from H2SO4 and SOA)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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
                      write(iulog,*) 'Error in initialization of be2to3'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'modes 2-3 ok'

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    ! Mode 4 (BC&OC + condensate from H2SO4 and SOA + wet phase (NH4)2SO4)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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
                         write(iulog,*) 'Error in initialization of be4'
                         stop
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'mode 4 ok'

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

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
                         write(iulog,*) 'be5to10 =', iwl, irelh, ictot, ifac, ifbc, ifaq, be5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                         write(iulog,*) 'Error in initialization of be5to10'
                         write(iulog,*) 'kcomp, abs((fabc-fbc)/fbc) =', kcomp, abs((fabc-fbc(ic))/fbc(ic))
                         stop
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'modes 5-10 ok'

    do ifil=40,50
       close (ifil)
    end do

  end subroutine initopt

  !********************************************************************************************
  subroutine initopt_lw

    !---------------------------------------------------------------
    !   Modified for new aerosol schemes by Alf Kirkevaag in January
    !   2006. Modified for new wavelength
    !   bands and look-up tables by Alf Kirkevaag in January 2014,
    !   and for SOA in August 2015.
    !---------------------------------------------------------------

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

    !     Then reading in the look-up table entries for each file (lwkcomp*.out)

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    !       Mode 0, BC(ax)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

    ifil = 11
    linmax=nlwbands
    do lin = 1,linmax
       read(39+ifil,996) kcomp, iwl, relh, spabs
       ka0(iwl)=spabs  ! unit m^2/g
    end do

    do iwl=1,nlwbands
       if(ka0(iwl)<=0.0_r8) then
          write(iulog,*) 'ka0 =', iwl, ka0(iwl)
          write(iulog,*) 'Error in initialization of ka0'
          stop
       endif
    enddo
    write(iulog,*)'lw mode 0 ok'

    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    !       Mode 1 (H2SO4 + condesate from H2SO4 and SOA)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

    ifil = 1
    linmax=nlwbands*10*6*16*6
    do lin = 1,linmax

       read(39+ifil,997) kcomp, iwl, relh, frombg, catot, frac, spabs

       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             goto 121
          endif
       end do
121    continue

       do ic=1,6
          if(abs(frombg-fombg(ic))<eps4) then
             ifombg=ic
             goto 122
          endif
       end do
122    continue

       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             goto 131
          endif
       end do
131    continue

       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             goto 141
          endif
       end do
141    continue

       ka1(iwl,irelh,ifombg,ictot,ifac)=spabs  ! unit m^2/g

       !      write(*,*) 'kcomp, ka =', kcomp, ka1(iwl,irelh,ifombg,ictot,ifac)
       !      if(ifil==1) write(iulog,*) 'iwl,irelh,ifombg,ictot,ifac,ka =', &
       !                  iwl,irelh,ictot,ifac,ka1(iwl,irelh,ifombg,ictot,ifac)

    end do  ! lin

    do iwl=1,nlwbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                if(ka1(iwl,irelh,ifombg,ictot,ifac)<=0.0_r8) then
                   write(iulog,*) 'ka1 =', iwl, irelh, ifombg, ictot, ifac, ka1(iwl,irelh,ifombg,ictot,ifac)
                   write(iulog,*) 'Error in initialization of ka1'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'lw new mode 1 ok'


    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    !       Modes 2 to 3 (BC or OC + condensate from H2SO4 and SOA)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

    linmax = nlwbands*10*16*6
    !      do ifil = 2,3
    do ifil = 2,2
       do lin = 1,linmax

          read(39+ifil,994) kcomp, iwl, relh, catot, frac, spabs

       	  do ic=1,10
             if(abs(relh-rh(ic))<eps4) then
                irelh=ic
                goto 61
             endif
	  end do
61        continue

 	  do ic=1,16
             if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
                ictot=ic
                goto 71
             endif
	  end do
71        continue

 	  do ic=1,6
             if(abs(frac-fac(ic))<eps4) then
                ifac=ic
                goto 72
             endif
	  end do
72        continue

          ka2to3(iwl,irelh,ictot,ifac,kcomp)=spabs  ! unit m^2/g

       end do  ! lin
    end do    ! ifil

    !   Prescribed dummy values for kcomp=3
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
                      write(iulog,*) 'ka2to3 =', iwl, irelh, ictot, ifac, ka2to3(iwl,irelh,ictot,ifac,kcomp)
                      write(iulog,*) 'Error in initialization of ka2to3'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'lw mode 2-3 ok'


    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    !       Mode 4 (BC&OC + condesate from H2SO4 and SOA + wetphase (NH4)2SO4)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

    ifil = 4
    linmax = nlwbands*10*6*16*6*6
    do lin = 1,linmax

       read(39+ifil,995) kcomp, iwl, relh, frbcbg, catot, frac, fraq, spabs

       do ic=1,10
          if(abs(relh-rh(ic))<eps4) then
             irelh=ic
             goto 81
          endif
       end do
81     continue

       do ic=1,6
          !	   if(abs(frbcbg-fbcbg(ic))<eps4) then
          !	   if(abs(frbcbg-fbcbg(ic))<eps3) then
          if(abs((frbcbg-fbcbg(ic))/fbcbg(ic))<eps2) then
             ifbcbg=ic
             goto 92
          endif
       end do
92     continue

       do ic=1,16
          if(abs((catot-cate(kcomp,ic))/cate(kcomp,ic))<eps2) then
             ictot=ic
             goto 91
          endif
       end do
91     continue

       do ic=1,6
          if(abs(frac-fac(ic))<eps4) then
             ifac=ic
             goto 101
          endif
       end do
101    continue

       do ic=1,6
          if(abs(fraq-faq(ic))<eps4) then
             ifaq=ic
             goto 111
          endif
       end do
111    continue

       ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)=spabs  ! unit m^2/g

       !      write(*,*) 'kcomp, ka =', kcomp, ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
    end do

    do iwl=1,nlwbands
       do irelh=1,10
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                      write(iulog,*) 'ka4 =', iwl, irelh, ifbcbg, ictot, ifac, ifaq, ka4(iwl,irelh,ifbcbg,ictot,ifac,ifaq)
                      write(iulog,*) 'Error in initialization of ka4'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'lw mode 4 ok'


    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
    !       Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.)
    !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc

    linmax = nlwbands*10*6*6*6*6
    do ifil = 5,10
       do lin = 1,linmax

          read(39+ifil,993) kcomp, iwl, relh, catot, frac, fabc, fraq, spabs

       	  do ic=1,10
             if(abs(relh-rh(ic))<eps4) then
                irelh=ic
                goto 11
             endif
	  end do
11        continue

 	  do ic=1,6
       !	   if(abs(catot-cat(kcomp,ic))<eps6) then
             if(abs((catot-cat(kcomp,ic))/cat(kcomp,ic))<eps2) then
                ictot=ic
                goto 21
             endif
	  end do
21        continue

 	  do ic=1,6
             if(abs(frac-fac(ic))<eps4) then
                ifac=ic
                goto 31
             endif
	  end do
31        continue

 	  do ic=1,6
       !	   if(abs(fabc-fbc(ic))<eps4) then
             if(abs((fabc-fbc(ic))/fbc(ic))<eps2) then
                ifbc=ic
                goto 41
             endif
	  end do
41        continue

	  do ic=1,6
             if(abs(fraq-faq(ic))<eps4) then
                ifaq=ic
                goto 51
             endif
	  end do
51        continue

          ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)=spabs  ! unit m^2/g

          !      write(*,*) 'kcomp, ka =', kcomp, ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
       end do
    end do


    do kcomp=5,10
       do iwl=1,nlwbands
          do irelh=1,10
             do ictot=1,6
                do ifac=1,6
                   do ifaq=1,6
                      if(ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                         write(iulog,*) 'ka5to10 =', iwl, irelh, ictot, ifac, ifbc, ifaq, ka5to10(iwl,irelh,ictot,ifac,ifbc,ifaq,kcomp)
                         write(iulog,*) 'Error in initialization of ka5to10'
                         stop
                      endif
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo

    write(iulog,*)'lw mode 5-10 ok'

993 format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 5-10
994 format(2I3,f8.3,2(x,e10.3),x,e12.5)           ! 2-3
995 format(2I3,f8.3,3(x,e10.3),f7.2,x,e12.5)      ! 4
996 format(2I3,f8.3,x,e12.5)                      ! 0
997 format(2I3,f8.3,3(x,e10.3),x,e12.5)           ! 1

    do ifil=40,50
       close (ifil)
    end do

  end subroutine initopt_lw

  !********************************************************************************************
  subroutine inputForInterpol (lchnk, ncol, rhum, xrh, irh1, &
       f_soana, xfombg, ifombg1, faitbc, xfbcbg, ifbcbg1,    &
       fnbc, xfbcbgn, ifbcbgn1, Nnatk, Cam, xct, ict1,       &
       focm, fcm, xfac, ifac1, fbcm, xfbc, ifbc1, faqm, xfaq, ifaq1)

    !
    ! Input arguments
    integer, intent(in)  :: lchnk                      ! chunk identifier
    integer, intent(in)  :: ncol                       ! number of atmospheric columns
    real(r8), intent(in) :: rhum(pcols,pver)           ! level relative humidity (fraction)
    real(r8), intent(in) :: f_soana(pcols,pver)        ! SOA/(SOA+H2SO4) mass fraction for the background in mode 1
    real(r8), intent(in) :: faitbc(pcols,pver)         ! BC/(BC + OC) mass fraction for the background in mode 4
    real(r8), intent(in) :: fnbc(pcols,pver)           ! BC/(BC + OC) mass fraction for the background in mode 14
    real(r8), intent(in) :: focm(pcols,pver,4)         ! fraction of added mass which is either SOA condensate or OC coagulate
    real(r8), intent(in) :: Cam(pcols,pver,nbmodes)    ! added internally mixed SO4+BC+OC concentration for a normalized mode
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! aerosol mode number concentration
    real(r8), intent(in) :: fcm(pcols,pver,nbmodes)    ! fraction of added mass which is either BC or OC/SOA (carbonaceous)
    real(r8), intent(in) :: fbcm(pcols,pver,nbmodes)   ! fraction of added mass as BC/(BC+OC)
    real(r8), intent(in) :: faqm(pcols,pver,nbmodes)   ! fraction of added sulfate which is from aqueous phase (ammonium sulfate)
    !
    ! Output arguments
    real(r8), intent(out) :: xrh(pcols,pver)           ! rhum for use in the interpolations
    integer,  intent(out) :: irh1(pcols,pver)
    real(r8), intent(out) :: xfombg(pcols,pver)        ! f_soana for use in the interpolations (mode 1)
    integer,  intent(out) :: ifombg1(pcols,pver)
    real(r8), intent(out) :: xfbcbg(pcols,pver)        ! faitbc for use in the interpolations (mode 4)
    integer,  intent(out) :: ifbcbg1(pcols,pver)
    real(r8), intent(out) :: xfbcbgn(pcols,pver)       ! fnbc for use in the interpolations (mode 14)
    integer,  intent(out) :: ifbcbgn1(pcols,pver)
    real(r8), intent(out) :: xct(pcols,pver,nmodes)    ! Cam/Nnatk for use in the interpolations
    integer,  intent(out) :: ict1(pcols,pver,nmodes)
    real(r8), intent(out) :: xfac(pcols,pver,nbmodes)  ! focm (1-4) or fcm (5-10) for use in the interpolations
    integer,  intent(out) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(out) :: xfbc(pcols,pver,nbmodes)  ! fbcm for use in the interpolations
    integer,  intent(out) :: ifbc1(pcols,pver,nbmodes)
    real(r8), intent(out) :: xfaq(pcols,pver,nbmodes)  ! faqm for use in the interpolations
    integer,  intent(out) :: ifaq1(pcols,pver,nbmodes)
    !
    ! Local variables
    integer k, icol, i, irelh
    real(r8) :: eps10 = 1.e-10_r8
    !------------------------------------------------------------------------
    !
    do k=1,pver
       do icol=1,ncol
          xrh(icol,k)  = min(max(rhum(icol,k),rh(1)),rh(10))
       end do
    end do

    do irelh=1,9
       do k=1,pver
          do icol=1,ncol
             if(xrh(icol,k) >= rh(irelh) .and. xrh(icol,k)<=rh(irelh+1)) then
                irh1(icol,k)=irelh
             endif
          end do
       end do
    end do

    do k=1,pver
       do icol=1,ncol
          ! find common xfombg, ifombg1 and ifombg2 for use in the interpolation routines
          xfombg(icol,k) =min(max(f_soana(icol,k),fombg(1)),fombg(6))
          ifombg1(icol,k)=int(5.0_r8*xfombg(icol,k)-eps10)+1
       end do
    enddo

    do k=1,pver
       do icol=1,ncol
          ! find common xfbcbg, ifbcbg1 and ifbcbg2 for use in the interpolation routines
          xfbcbg(icol,k) =min(max(faitbc(icol,k),fbcbg(1)),fbcbg(6))
          ifbcbg1(icol,k)=min(max(int(4*log10(xfbcbg(icol,k))+6),1),5)

          ! find common xfbcbgn, ifbcbgn1 and ifbcbgn2 for use in the interpolation routines
          xfbcbgn(icol,k) =min(max(fnbc(icol,k),fbcbg(1)),fbcbg(6))
          ifbcbgn1(icol,k)=min(max(int(4*log10(xfbcbgn(icol,k))+6),1),5)
       end do
    enddo

    do i=1,4
       do k=1,pver
          do icol=1,ncol
             ! find common xfac, ifac1 and ifac2 for use in the interpolation routines
             xfac(icol,k,i) =min(max(focm(icol,k,i),fac(1)),fac(6))
             ifac1(icol,k,i)=int(5.0_r8*xfac(icol,k,i)-eps10)+1
          end do
       enddo
    enddo
    do i=5,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfac, ifac1 and ifac2 for use in the interpolation routines
             xfac(icol,k,i) =min(max(fcm(icol,k,i),fac(1)),fac(6))
             ifac1(icol,k,i)=int(5.0_r8*xfac(icol,k,i)-eps10)+1
          end do
       enddo
    enddo

    do i=1,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfbc, ifbc1 and ifbc2 for use in the interpolation routines
             xfbc(icol,k,i) =min(max(fbcm(icol,k,i),fbc(1)),fbc(6))
             ifbc1(icol,k,i)=min(max(int(4*log10(xfbc(icol,k,i))+6),1),5)
          end do
       enddo
    enddo

    do i=1,nbmodes
       do k=1,pver
          do icol=1,ncol
             ! find common xfaq, ifaq1 and ifaq2 for use in the interpolation routines
             xfaq(icol,k,i) =min(max(faqm(icol,k,i),faq(1)),faq(6))
             ifaq1(icol,k,i)=int(5.0_r8*xfaq(icol,k,i)-eps10)+1
          end do
       enddo
    enddo

    ! find common xct, ict1 and ict2 for use in the interpolation routines
    do i=1,4
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=min(max(Cam(icol,k,i)/(Nnatk(icol,k,i)+eps),cate(i,1)),cate(i,16))
             if(i.le.2) then
                ict1(icol,k,i)=min(max(int(3*log10(xct(icol,k,i))+19.666_r8),1),15)
             elseif(i.eq.3) then ! mode not used
                xct(icol,k,i)=cate(i,1)
                ict1(icol,k,i)=1
             else
                ict1(icol,k,i)=min(max(int(3*log10(xct(icol,k,i))+13.903_r8),1),15)
             endif
          end do
       end do
    end do

    do i=5,10
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=min(max(Cam(icol,k,i)/(Nnatk(icol,k,i)+eps),cat(i,1)),cat(i,6))
             if(i.eq.5) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.824_r8),1),5)
             elseif(i.eq.6) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.523_r8),1),5)
             elseif(i.eq.7) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.699_r8),1),5)
             elseif(i.eq.8) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+5.921_r8),1),5)
             elseif(i.eq.9) then
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.301_r8),1),5)
             else
                ict1(icol,k,i)=min(max(int(log10(xct(icol,k,i))+4.699_r8),1),5)
             endif
          end do
       end do
    end do

    do i=11,nmodes ! for the externally mixed modes 11-14 (now only 12 and 14)
       do k=1,pver
          do icol=1,ncol
             xct(icol,k,i)=cate(i-10,1)
             ict1(icol,k,i)=1
          end do
       end do
    end do

    return

  end subroutine inputForInterpol

  !********************************************************************************************
  subroutine interpol0 (lchnk, ncol, daylight, Nnatk, omega, gass, bex, ske, lw_on, kabs)
    !
    ! Arguments
    integer  , intent(in)  :: lchnk                              ! chunk identifier
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

  !********************************************************************************************
  subroutine interpol1 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, xfombg, ifombg1, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    !
    ! Arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                      !alt        bex(icol,k,kc10,i)=(d2mx(1)*bex1+dxm1(1)*bex2)*invd(1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                      !alt        a=(log(bex2)-log(bex1))*invd(1)
                      !alt        b=(t_rh2*log(bex1)-t_rh1*log(bex2))*invd(1)
                      !alt        bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i

                ! if(bex(icol,k,kc10,8)<1.e-20_r8) then
                ! write(*,995) 'bex(8)=', kc10, t_xrh, t_xct, t_xfac, t_xfombg, bex(icol,k,kc10,8)
                ! endif
             else  ! daylight


                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

    return
  end subroutine interpol1


  !********************************************************************************************
  subroutine interpol2to3 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, &
       xct, ict1, xfac, ifac1, omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)

    ! Output arguments
    real(r8), intent(out) :: omega(pcols,pver,0:nmodes,nbands) ! spectral modal single scattering albedo
    real(r8), intent(out) :: gass(pcols,pver,0:nmodes,nbands)  ! spectral modal asymmetry factor
    real(r8), intent(out) :: bex(pcols,pver,0:nmodes,nbands)   ! spectral modal extinction coefficient
    real(r8), intent(out) :: ske(pcols,pver,0:nmodes,nbands)   ! spectral modal specific extinction coefficient
    real(r8), intent(out) :: kabs(pcols,pver,0:nmodes,nlwbands)! LW spectral modal specific absorption coefficient
    !
    ! Local variables
    integer i, kcomp, k, icol, kc10
    real(r8) a, b
    integer t_irh1, t_irh2, t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) t_fac1, t_fac2, t_xfac, t_xrh, t_xct, t_rh1, t_rh2,t_cat1, t_cat2
    real(r8) d2mx(3), dxm1(3), invd(3)
    real(r8) opt3d(2,2,2)
    real(r8) ome1, ome2, ge1, ge2, bex1, bex2, ske1, ske2
    real(r8) kabs1, kabs2
    !---------------------------------------

    ! write(*,*) 'Before kcomp-loop'
    ! do kcomp=2,3
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

             ! write(*,*) 't_irh1,t_irh2=',t_irh1,t_irh2
             ! write(*,*) 't_ict1,t_ict2=',t_ict1,t_ict2
             ! write(*,*) 't_ifc1,t_ifc2=',t_ifc1,t_ifc2
             ! write(*,*) 't_ifa1,t_ifa2=',t_ifa1,t_ifa2

             t_rh1  = rh(t_irh1)
             t_rh2  = rh(t_irh2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)

             ! write(*,*) 't_rh1,t_rh2,t_cat1,t_cat2=',t_rh1,t_rh2,t_cat1,t_cat2
             ! write(*,*) 't_fac1,t_fac2=',t_fac1,t_fac2

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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) gass(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37) then
                      bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight



                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
                ! aerosol extinction used for LW size information

                i=4
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
                ! write(*,*) 'Before bex'
                if(t_xrh <= 0.37) then
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight



             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                ! write(*,*) 'Before ske'
                if(t_xrh <= 0.37) then
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i



             if (lw_on) then

                ! LW optical parameters
                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
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

    return
  end subroutine interpol2to3

  !********************************************************************************************

  subroutine interpol4 (lchnk, ncol, daylight, xrh, irh1, mplus10, Nnatk, xfbcbg, ifbcbg1, &
       xct, ict1, xfac, ifac1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
    integer, intent(in) :: ncol                        ! number of atmospheric columns
    integer, intent(in) :: mplus10                     ! mode number (0) or number + 10 (1)
    logical, intent(in) :: daylight(pcols)             ! only daylight calculations if .true.
    logical, intent(in) :: lw_on                       ! LW calculations are performed if true
    real(r8), intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in) :: xrh(pcols,pver)            ! level relative humidity (fraction)
    integer,  intent(in) :: irh1(pcols,pver)
    real(r8), intent(in) :: xfbcbg(pcols,pver)         ! mass fraction BC/(BC+OC) for the background mode
    integer,  intent(in) :: ifbcbg1(pcols,pver)
    real(r8), intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in) :: ict1(pcols,pver,nmodes)
    real(r8), intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8), intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer,  intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! Output arguments
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

    ! write(*,*) 'Before kcomp-loop'
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

                do i=1,nbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kc10,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kc10,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kc10,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   bex(icol,k,kc10,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                   bex(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             endif  ! daylight

             do i=4,4            ! i = wavelength index

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ske(icol,k,kc10,i)=((t_rh2-t_xrh)*ske1+(t_xrh-t_rh1)*ske2) &
                        /(t_rh2-t_rh1)
                else
                   a=(log(ske2)-log(ske1))/(t_rh2-t_rh1)
                   b=(t_rh2*log(ske1)-t_rh1*log(ske2))/(t_rh2-t_rh1)
                   ske(icol,k,kc10,i)=e**(a*t_xrh+b)
                endif

             end do ! i



             if (lw_on) then

                ! LW optical parameters

                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                   ! write(*,*) 'Before kabs'
                   if(t_xrh <= 0.37_r8) then
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

       ! write(*,*) 'kcomp, omega(1,26,kc10,4)=', kcomp, omega(1,26,kc10,4)
       ! write(*,*) 'kcomp, gass(1,26,kc10,4)=', kcomp, gass(1,26,kc10,4)
       ! write(*,*) 'kcomp, bex(1,26,kc10,4)=', kcomp, bex(1,26,kc10,4)
       ! write(*,*) 'kcomp, ske(1,26,kc10,4)=', kcomp, ske(1,26,kc10,4)

    end do  ! kcomp

  end subroutine interpol4

  !********************************************************************************************
  subroutine interpol5to10 (lchnk, ncol, daylight, xrh, irh1, Nnatk, xct, ict1, &
       xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1, &
       omega, gass, bex, ske, lw_on, kabs)

    ! Input arguments
    integer, intent(in) :: lchnk                       ! chunk identifier
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

    ! write(*,*) 'Before kcomp-loop'
    do kcomp=5,10

       ! write(*,*) 'Before init-loop', kcomp
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

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before omega'
                   omega(icol,k,kcomp,i)=((t_rh2-t_xrh)*ome1+(t_xrh-t_rh1)*ome2) &
                        /(t_rh2-t_rh1)
                   ! write(*,*) omega(icol,k,kcomp,i)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before gass'
                   gass(icol,k,kcomp,i)=((t_rh2-t_xrh)*ge1+(t_xrh-t_rh1)*ge2) &
                        /(t_rh2-t_rh1)

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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
                   ! write(*,*) 'Before bex'
                   if(t_xrh <= 0.37_r8) then
                      bex(icol,k,kcomp,i)=((t_rh2-t_xrh)*bex1+(t_xrh-t_rh1)*bex2) &
                           /(t_rh2-t_rh1)
                   else
                      a=(log(bex2)-log(bex1))/(t_rh2-t_rh1)
                      b=(t_rh2*log(bex1)-t_rh1*log(bex2))/(t_rh2-t_rh1)
                      bex(icol,k,kcomp,i)=e**(a*t_xrh+b)
                   endif

                end do ! i
             else  ! daylight


                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

                !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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



             if (lw_on) then

                ! LW optical parameters

                do i=1,nlwbands            ! i = wavelength index

                   !ccccccccc1ccccccccc2ccccccccc3ccccccccc4ccccccccc5ccccccccc6ccccccccc7cc
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

  !********************************************************************************************
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
