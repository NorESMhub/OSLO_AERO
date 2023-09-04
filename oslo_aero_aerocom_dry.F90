module oslo_aero_aerocom_dry

#ifdef AEROCOM

  use shr_kind_mod            , only: r8 => shr_kind_r8
  use ppgrid                  , only: pcols, pver
  use cam_logfile             , only: iulog
  !
  use oslo_aero_params       , only: nmodes, nbmodes
  use oslo_aero_sw_tables     , only: cate, cat, fac, faq, fbc, fombg, fbcbg, nbmp1
  use oslo_aero_linear_interp , only: lininterpol3dim, lininterpol4dim, lininterpol5dim  
  use oslo_aero_control       , only: oslo_aero_getopts, dir_string_length

  implicit none
  private

  ! Set by init_dryp Mode0
  real(r8) :: a0cintbg, a0cintbg05, a0cintbg125
  real(r8) :: a0aaeros, a0aaerol, a0vaeros, a0vaerol

  ! Used by init_dryp Mode1
  real(r8) :: a1var(19,6,16,6)

  ! Used by init_dryp Mode2to3
  real(r8) :: a2to3var(19,16,6,2:3)

  ! Used by init_dryp Mode4
  real(r8) :: a4var(19,6,16,6,6)

  ! Used by init_dryp Mode5
  real(r8) :: a5to10var(19,6,6,6,6,5:10)

  type, public :: aerodry_prop_type
     ! modal mass concentrations (cint), area (aaero) and volume (vaero)
     ! (for AeroCom determination of particle effective radii) of each constituent.
     ! cint*05 and cint*125 are for r<0.5um and r>1.25um, respectively.
     ! aaeros and vaeros are integrated over r<0.5um, and aaerol and vaerol over r>0.5um.

     real(r8) :: cintbg(pcols,pver,0:nbmodes)
     real(r8) :: cintbg05(pcols,pver,0:nbmodes)
     real(r8) :: cintbg125(pcols,pver,0:nbmodes)
     real(r8) :: cintbc(pcols,pver,0:nbmodes)
     real(r8) :: cintbc05(pcols,pver,0:nbmodes)
     real(r8) :: cintbc125(pcols,pver,0:nbmodes)
     real(r8) :: cintoc(pcols,pver,0:nbmodes)
     real(r8) :: cintoc05(pcols,pver,0:nbmodes)
     real(r8) :: cintoc125(pcols,pver,0:nbmodes)
     real(r8) :: cintsc(pcols,pver,0:nbmodes)
     real(r8) :: cintsc05(pcols,pver,0:nbmodes)
     real(r8) :: cintsc125(pcols,pver,0:nbmodes)
     real(r8) :: cintsa(pcols,pver,0:nbmodes)
     real(r8) :: cintsa05(pcols,pver,0:nbmodes)
     real(r8) :: cintsa125(pcols,pver,0:nbmodes)
     real(r8) :: aaeros(pcols,pver,0:nbmodes)
     real(r8) :: aaerol(pcols,pver,0:nbmodes)
     real(r8) :: vaeros(pcols,pver,0:nbmodes)
     real(r8) :: vaerol(pcols,pver,0:nbmodes)

     real(r8) :: aaerosn(pcols,pver,nbmp1:nmodes)
     real(r8) :: aaeroln(pcols,pver,nbmp1:nmodes)
     real(r8) :: vaerosn(pcols,pver,nbmp1:nmodes)
     real(r8) :: vaeroln(pcols,pver,nbmp1:nmodes)
     real(r8) :: cknorm(pcols,pver,0:nmodes)
     real(r8) :: cknlt05(pcols,pver,0:nmodes)
     real(r8) :: ckngt125(pcols,pver,0:nmodes)

   contains
     procedure :: intdrypar0
     procedure :: intdrypar1
     procedure :: intdrypar2to3
     procedure :: intdrypar4
     procedure :: intdrypar5to10
     procedure :: zero
     procedure :: update

  end type aerodry_prop_type

  type(aerodry_prop_type), public :: aerodry_prop

  public :: initdryp

! ==========================================================
contains
! ==========================================================

  subroutine initdryp()

    !Purpose: To read in the AeroCom look-up tables for calculation of dry
    !     aerosol size and mass distribution properties. The grid for discrete
    !     input-values in the look-up tables is defined in opptab.

    !     Tabulating the 'aerodryk'-files to save computing time. Routine
    !     originally made by  Alf Kirkevaag, and modified for new aerosol
    !     schemes in January 2006.
    !     Updated for new kcomp1.out including condensed SOA - Alf Kirkevåg,
    !     May 2013, and extended for new SOA treatment October 2015.
    !     Modified for optimized added masses and mass fractions for
    !     concentrations from condensation, coagulation or cloud-processing
    !     - Alf Kirkevaag, May 2016.
    !     Modified for optimized added masses and mass fractions for concentrations from
    !     condensation, coagulation or cloud-processing - Alf Kirkevaag, May 2016.

    ! local variables
    integer  :: iv, kcomp, ifombg, ifbcbg, ictot, ifac, ifbc, ifaq
    integer  :: ic, ifil, lin
    real(r8) :: frombg, frbcbg, catot, frac, fabc, fraq
    real(r8) :: cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125
    real(r8) :: cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125
    real(r8) :: cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol
    real(r8) :: eps2 = 1.e-2_r8
    real(r8) :: eps4 = 1.e-4_r8
    real(r8) :: eps6 = 1.e-6_r8
    real(r8) :: eps7 = 1.e-7_r8
    character(len=dir_string_length) :: aerotab_table_dir

    call oslo_aero_getopts(aerotab_table_dir_out = aerotab_table_dir)

    open(20,file=trim(aerotab_table_dir)//'/aerodryk0.out'  ,form='formatted',status='old')
    open(21,file=trim(aerotab_table_dir)//'/aerodryk1.out'  ,form='formatted',status='old')
    open(11,file=trim(aerotab_table_dir)//'/aerodryk2.out'  ,form='formatted',status='old')
    open(12,file=trim(aerotab_table_dir)//'/aerodryk3.out'  ,form='formatted',status='old')
    open(13,file=trim(aerotab_table_dir)//'/aerodryk4.out'  ,form='formatted',status='old')
    open(14,file=trim(aerotab_table_dir)//'/aerodryk5.out'  ,form='formatted',status='old')
    open(15,file=trim(aerotab_table_dir)//'/aerodryk6.out'  ,form='formatted',status='old')
    open(16,file=trim(aerotab_table_dir)//'/aerodryk7.out'  ,form='formatted',status='old')
    open(17,file=trim(aerotab_table_dir)//'/aerodryk8.out'  ,form='formatted',status='old')
    open(18,file=trim(aerotab_table_dir)//'/aerodryk9.out'  ,form='formatted',status='old')
    open(19,file=trim(aerotab_table_dir)//'/aerodryk10.out' ,form='formatted',status='old')

    ! Skipping the header-text in all input files (Later: use it to check AeroTab - CAM5-Oslo consistency!)
    do ifil = 11,21
       call checkTableHeader (ifil)
    enddo
    !
    !-------------------------------------------
    ! Mode 0, BC(ax)
    !-------------------------------------------
    !
    read(20,996) kcomp, cintbg, cintbg05, cintbg125, aaeros, aaerol, vaeros, vaerol

    ! no ictot-, ifac-, ifbc- or ifaq-dependency for this mode,
    ! since BC(ax) is purely externally mixed

    a0cintbg=cintbg
    a0cintbg05=cintbg05
    a0cintbg125=cintbg125

    a0aaeros=aaeros
    a0aaerol=aaerol
    a0vaeros=vaeros
    a0vaerol=vaerol
    write(iulog,*)'mode 0 ok'

    !-------------------------------------------
    ! Mode 1 (H2SO4 and SOA + condensate from H2SO4 and SOA)
    !-------------------------------------------

    do lin = 1,576     ! 6x16x6
       read(21,997) kcomp, frombg, catot, frac,               &
            cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
            cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &
            cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

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

       ! no ifombg-dependency for this mode, since all catot
       ! comes from condensate or from wet-phase sulfate

       a1var(1,ifombg,ictot,ifac)  =cintbg
       a1var(2,ifombg,ictot,ifac)  =cintbg05
       a1var(3,ifombg,ictot,ifac)  =cintbg125
       a1var(4,ifombg,ictot,ifac)  =cintbc
       a1var(5,ifombg,ictot,ifac)  =cintbc05
       a1var(6,ifombg,ictot,ifac)  =cintbc125
       a1var(7,ifombg,ictot,ifac)  =cintoc
       a1var(8,ifombg,ictot,ifac)  =cintoc05
       a1var(9,ifombg,ictot,ifac)  =cintoc125
       a1var(10,ifombg,ictot,ifac) =cintsc
       a1var(11,ifombg,ictot,ifac) =cintsc05
       a1var(12,ifombg,ictot,ifac) =cintsc125
       a1var(13,ifombg,ictot,ifac) =cintsa
       a1var(14,ifombg,ictot,ifac) =cintsa05
       a1var(15,ifombg,ictot,ifac) =cintsa125
       a1var(16,ifombg,ictot,ifac) =aaeros
       a1var(17,ifombg,ictot,ifac) =aaerol
       a1var(18,ifombg,ictot,ifac) =vaeros
       a1var(19,ifombg,ictot,ifac) =vaerol

       if(cintsa<cintsa05) then
          write(*,*) 'cintsatot =', ictot, ifac, ifaq, cintsa, cintsa05, cintsa125
       end if
    end do  ! lin

    do iv=1,19
       do ifombg=1,6
          do ictot=1,16
             do ifac=1,6
                if(a1var(iv,ifombg,ictot,ifac)<=0.0_r8) then
                   write(*,*) 'a1var =', iv, ifombg, ictot, ifac, a1var(iv,ifombg,ictot,ifac)
                   write(*,*) 'Error in initialization of a1var'
                   stop
                endif
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'new aerodry mode 1 ok'

    !-------------------------------------------
    ! Modes 2 to 3 (BC/OC + condesate from H2SO4 and SOA)
    ! Note that mode 3 is no longer active
    !-------------------------------------------

    do lin = 1,96     ! 16x6
       read(11,994) &
            kcomp, catot, frac,                                        &
            cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
            cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &
            cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

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

       ! no ifbc- or ifaq-dependency for these modes,
       ! since all catot comes from condensed SOA and H2SO4

       a2to3var(1,ictot,ifac,kcomp) = cintbg
       a2to3var(2,ictot,ifac,kcomp) = cintbg05
       a2to3var(3,ictot,ifac,kcomp) = cintbg125
       a2to3var(4,ictot,ifac,kcomp) = cintbc
       a2to3var(5,ictot,ifac,kcomp) = cintbc05
       a2to3var(6,ictot,ifac,kcomp) = cintbc125
       a2to3var(7,ictot,ifac,kcomp) = cintoc
       a2to3var(8,ictot,ifac,kcomp) = cintoc05
       a2to3var(9,ictot,ifac,kcomp) = cintoc125
       a2to3var(10,ictot,ifac,kcomp) = cintsc
       a2to3var(11,ictot,ifac,kcomp) = cintsc05
       a2to3var(12,ictot,ifac,kcomp) = cintsc125
       a2to3var(13,ictot,ifac,kcomp) = cintsa
       a2to3var(14,ictot,ifac,kcomp) = cintsa05
       a2to3var(15,ictot,ifac,kcomp) = cintsa125
       a2to3var(16,ictot,ifac,kcomp) = aaeros
       a2to3var(17,ictot,ifac,kcomp) = aaerol
       a2to3var(18,ictot,ifac,kcomp) = vaeros
       a2to3var(19,ictot,ifac,kcomp) = vaerol
    end do  ! lin

    ! Prescribed dummy values for unused kcomp=3
    kcomp=3
    do ictot=1,16
       do ifac=1,6
          do iv=1,19
             a2to3var(iv,ictot,ifac,kcomp)=1.0_r8
          enddo
       enddo
    enddo

    do iv=1,19
       do kcomp=2,3
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(a2to3var(iv,ictot,ifac,kcomp)<=0.0_r8) then
                      write(*,*) 'a2to3var =', iv, kcomp, ictot, a2to3var(iv,ictot,ifac,kcomp)
                      write(*,*) 'Error in initialization of a2to3var'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerodry mode 2-3 ok'
    !
    !-------------------------------------------
    !       Mode 4 (BC&OC + condensate from H2SO4 and SOA + wetphase (NH4)2SO4)
    !-------------------------------------------
    !
    do lin = 1,3456     ! 16x6x6x6

       read(13,995) &
            kcomp, frbcbg, catot, frac, fraq,                          &
            cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
            cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &
            cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

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

       ! no ifbc-dependency for this mode, since all catot
       ! comes from condensed or wet-phase sulfate

       a4var(1,ifbcbg,ictot,ifac,ifaq) = cintbg
       a4var(2,ifbcbg,ictot,ifac,ifaq) = cintbg05
       a4var(3,ifbcbg,ictot,ifac,ifaq) = cintbg125
       a4var(4,ifbcbg,ictot,ifac,ifaq) = cintbc
       a4var(5,ifbcbg,ictot,ifac,ifaq) = cintbc05
       a4var(6,ifbcbg,ictot,ifac,ifaq) = cintbc125
       a4var(7,ifbcbg,ictot,ifac,ifaq) = cintoc
       a4var(8,ifbcbg,ictot,ifac,ifaq) = cintoc05
       a4var(9,ifbcbg,ictot,ifac,ifaq) = cintoc125
       a4var(10,ifbcbg,ictot,ifac,ifaq) = cintsc
       a4var(11,ifbcbg,ictot,ifac,ifaq) = cintsc05
       a4var(12,ifbcbg,ictot,ifac,ifaq) = cintsc125
       a4var(13,ifbcbg,ictot,ifac,ifaq) = cintsa
       a4var(14,ifbcbg,ictot,ifac,ifaq) = cintsa05
       a4var(15,ifbcbg,ictot,ifac,ifaq) = cintsa125
       a4var(16,ifbcbg,ictot,ifac,ifaq) = aaeros
       a4var(17,ifbcbg,ictot,ifac,ifaq) = aaerol
       a4var(18,ifbcbg,ictot,ifac,ifaq) = vaeros
       a4var(19,ifbcbg,ictot,ifac,ifaq) = vaerol

       if(cintsa<cintsa05) then
          write(*,*) 'cintsatot =', ictot, ifac, ifaq, cintsa, cintsa05, cintsa125
       end if
    end do  ! lin

    do iv=1,19
       do ifbcbg=1,6
          do ictot=1,16
             do ifac=1,6
                do ifaq=1,6
                   if(a4var(iv,ifbcbg,ictot,ifac,ifaq)<=0.0_r8) then
                      write(*,*) 'a4var =', iv, ifbcbg, ictot, ifac, ifaq, &
                           a4var(iv,ifbcbg,ictot,ifac,ifaq)
                      write(*,*) 'Error in initialization of a4var'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerodry mode 4 ok'
    !
    !-------------------------------------------
    !       Modes 5 to 10 (mineral and seasalt-modes + cond./coag./aq.)
    !-------------------------------------------
    !
    do ifil = 5,10
       do lin = 1,1296     ! 6x6x6x6

          read(9+ifil,993) kcomp, catot, frac, fabc, fraq,            &
               cintbg, cintbg05, cintbg125, cintbc, cintbc05, cintbc125,  &
               cintoc, cintoc05, cintoc125, cintsc, cintsc05, cintsc125,  &
               cintsa, cintsa05, cintsa125, aaeros, aaerol, vaeros, vaerol

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

          a5to10var(1,ictot,ifac,ifbc,ifaq,kcomp)  = cintbg
          a5to10var(2,ictot,ifac,ifbc,ifaq,kcomp)  = cintbg05
          a5to10var(3,ictot,ifac,ifbc,ifaq,kcomp)  = cintbg125
          a5to10var(4,ictot,ifac,ifbc,ifaq,kcomp)  = cintbc
          a5to10var(5,ictot,ifac,ifbc,ifaq,kcomp)  = cintbc05
          a5to10var(6,ictot,ifac,ifbc,ifaq,kcomp)  = cintbc125
          a5to10var(7,ictot,ifac,ifbc,ifaq,kcomp)  = cintoc
          a5to10var(8,ictot,ifac,ifbc,ifaq,kcomp)  = cintoc05
          a5to10var(9,ictot,ifac,ifbc,ifaq,kcomp)  = cintoc125
          a5to10var(10,ictot,ifac,ifbc,ifaq,kcomp) = cintsc
          a5to10var(11,ictot,ifac,ifbc,ifaq,kcomp) = cintsc05
          a5to10var(12,ictot,ifac,ifbc,ifaq,kcomp) = cintsc125
          a5to10var(13,ictot,ifac,ifbc,ifaq,kcomp) = cintsa
          a5to10var(14,ictot,ifac,ifbc,ifaq,kcomp) = cintsa05
          a5to10var(15,ictot,ifac,ifbc,ifaq,kcomp) = cintsa125
          a5to10var(16,ictot,ifac,ifbc,ifaq,kcomp) = aaeros
          a5to10var(17,ictot,ifac,ifbc,ifaq,kcomp) = aaerol
          a5to10var(18,ictot,ifac,ifbc,ifaq,kcomp) = vaeros
          a5to10var(19,ictot,ifac,ifbc,ifaq,kcomp) = vaerol
       end do  ! lin
    end do    ! ifil

    do iv=1,19
       do kcomp=5,10
          do ictot=1,6
             do ifac=1,6
                do ifaq=1,6
                   if(a5to10var(iv,ictot,ifac,ifbc,ifaq,kcomp)<=0.0_r8) then
                      write(*,*) 'a5to10var =', iv, kcomp, ictot, ifac, ifbc, ifaq, &
                           a5to10var(iv,ictot,ifac,ifbc,ifaq,kcomp)
                      write(*,*) 'Error in initialization of a5to10var'
                      stop
                   endif
                enddo
             enddo
          enddo
       enddo
    enddo
    write(iulog,*)'aerodry mode 5-10 ok'

993 format(I2,23e10.3)
994 format(I2,21e10.3)
995 format(I2,23e10.3)
996 format(I2,7e11.4)
997 format(I2,22e10.3)

    do ifil=10,21
       close (ifil)
    end do

  end subroutine initdryp

  ! ==========================================================
  subroutine intdrypar0 (this, lchnk, ncol, Nnatk)

    ! arguments
    class (aerodry_prop_type):: this
    integer  , intent(in)    :: lchnk                      ! chunk identifier
    integer  , intent(in)    :: ncol                       ! number of atmospheric columns
    real(r8) , intent(in)    :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration 

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: i, ierr, kcomp, k, icol
    parameter (eps=1.0e-60_r8)
    !-----------------------------------------------------------

    ! Mode 0, BC(ax):
    kcomp=0
    call aerodry_prop%zero(kcomp,ncol)
    do k=1,pver
       do icol=1,ncol
          if(Nnatk(icol,k,kcomp) > 0.0_r8) then
             this%cintbg(icol,k,kcomp)    = a0cintbg
             this%cintbg05(icol,k,kcomp)  = a0cintbg05
             this%cintbg125(icol,k,kcomp) = a0cintbg125
             this%cintbc(icol,k,kcomp)    = eps
             this%cintbc05(icol,k,kcomp)  = eps
             this%cintbc125(icol,k,kcomp) = eps
             this%cintoc(icol,k,kcomp)    = eps
             this%cintoc05(icol,k,kcomp)  = eps
             this%cintoc125(icol,k,kcomp) = eps
             this%cintsc(icol,k,kcomp)    = eps
             this%cintsc05(icol,k,kcomp)  = eps
             this%cintsc125(icol,k,kcomp) = eps
             this%cintsa(icol,k,kcomp)    = eps
             this%cintsa05(icol,k,kcomp)  = eps
             this%cintsa125(icol,k,kcomp) = eps
             this%aaeros(icol,k,kcomp)    = a0aaeros
             this%aaerol(icol,k,kcomp)    = a0aaerol
             this%vaeros(icol,k,kcomp)    = a0vaeros
             this%vaerol(icol,k,kcomp)    = a0vaerol
          endif
          this%cknorm(icol,k,kcomp)  = a0cintbg
          this%cknlt05(icol,k,kcomp) = a0cintbg05
          this%ckngt125(icol,k,kcomp)= a0cintbg125
       end do ! icol
    end do ! k

  end subroutine intdrypar0

  ! ==========================================================
  subroutine intdrypar1 (this, lchnk, ncol, Nnatk, xfombg, ifombg1, xct, ict1, xfac, ifac1)

    ! Output arguments: Modal mass concentrations (cint), area (aaero) and volume (vaero)
    ! (for AeroCom determination of particle effective radii) of each constituent. cint*05
    ! and cint*125 are  for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
    ! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.

    ! Arguments
    class(aerodry_prop_type) :: this
    integer, intent(in)   :: lchnk                       ! chunk identifier
    integer, intent(in)   :: ncol                        ! number of atmospheric columns
    real(r8), intent(in)  :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8), intent(in)  :: xfombg(pcols,pver)         ! SOA/(SOA+H2SO4) for the background mode (1)
    integer,  intent(in)  :: ifombg1(pcols,pver)
    real(r8), intent(in)  :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer,  intent(in)  :: ict1(pcols,pver,nmodes)
    real(r8), intent(in)  :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer,  intent(in)  :: ifac1(pcols,pver,nbmodes)

    ! Local variables
    real(r8) :: a, b, e, eps
    integer  :: iv, kcomp, k, icol
    integer  :: t_ifo1, t_ifo2
    integer  :: t_ict1, t_ict2, t_ifc1, t_ifc2
    real(r8) :: t_xct,  t_cat1, t_cat2
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_fombg1, t_fombg2, t_xfombg, t_xfombgn
    real(r8) :: d2mx(3), dxm1(3), invd(3)
    real(r8) :: opt3d(2,2,2)
    real(r8) :: opt1, opt2, opt
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    !---------------------
    ! Mode 1, SO4(Ait):
    !---------------------
    kcomp=1
    call this%zero(kcomp,ncol)

    do k=1,pver
       do icol=1,ncol
          if(Nnatk(icol,k,kcomp)>0.0_r8) then

             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing
             t_ifo1 = ifombg1(icol,k)
             t_ifo2 = t_ifo1+1
             t_fombg1 = fombg(t_ifo1)
             t_fombg2 = fombg(t_ifo2)
             t_xfombg = xfombg(icol,k)
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_xct  = xct(icol,k,kcomp)
             t_xfac = xfac(icol,k,kcomp)

             ! partial lengths along each dimension (1-3) for interpolation
             d2mx(1) = (t_fombg2-t_xfombg)
             dxm1(1) = (t_xfombg-t_fombg1)
             invd(1) = 1.0_r8/(t_fombg2-t_fombg1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)

             do iv=1,19  ! variable number

                ! end points as basis for multidimentional linear interpolation
                opt3d(1,1,1)=a1var(iv,t_ifo1,t_ict1,t_ifc1)
                opt3d(1,1,2)=a1var(iv,t_ifo1,t_ict1,t_ifc2)
                opt3d(1,2,1)=a1var(iv,t_ifo1,t_ict2,t_ifc1)
                opt3d(1,2,2)=a1var(iv,t_ifo1,t_ict2,t_ifc2)
                opt3d(2,1,1)=a1var(iv,t_ifo2,t_ict1,t_ifc1)
                opt3d(2,1,2)=a1var(iv,t_ifo2,t_ict1,t_ifc2)
                opt3d(2,2,1)=a1var(iv,t_ifo2,t_ict2,t_ifc1)
                opt3d(2,2,2)=a1var(iv,t_ifo2,t_ict2,t_ifc2)

                ! interpolation in the fac and cat dimensions
                call lininterpol3dim (d2mx, dxm1, invd, opt3d, opt1, opt2)

                ! finally, interpolation in the fombg dimension
                opt = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

                ! update the properties
                call this%update(kcomp, k, icol, iv, opt)

             end do ! iv=1,19
          endif
       end do ! icol
    end do ! k


    ! Dry parameters for externally mixed mode 11, SO4(n):
    kcomp=11
    do k=1,pver
       do icol=1,ncol
          ! Neither total background concentrations (OM + sulfate)
          ! nor areas & volumes depend on fombg:
          this%cknorm(icol,k,kcomp)   = a1var(1,1,1,1)
          this%cknlt05(icol,k,kcomp)  = a1var(2,1,1,1)
          this%ckngt125(icol,k,kcomp) = a1var(3,1,1,1)
          this%aaerosn(icol,k,kcomp)  = a1var(16,1,1,1)
          this%aaeroln(icol,k,kcomp)  = a1var(17,1,1,1)
          this%vaerosn(icol,k,kcomp)  = a1var(18,1,1,1)
          this%vaeroln(icol,k,kcomp)  = a1var(19,1,1,1)
       end do ! icol
    end do ! k

  end subroutine intdrypar1

  ! ==========================================================
  subroutine intdrypar2to3 (this, lchnk, ncol, Nnatk, xct, ict1, xfac, ifac1)

    ! Modal mass concentrations (cint), area (aaero)
    ! and volume (vaero) (for AeroCom determination of particle
    ! effective radii) of each constituent. cint*05 and cint*125 are
    ! for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
    ! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.

    ! arguments
    class(aerodry_prop_type) :: this
    integer  , intent(in) :: lchnk                       ! chunk identifier
    integer  , intent(in) :: ncol                        ! number of atmospheric columns
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: iv, kcomp, k, icol
    integer  :: t_ict1, t_ict2
    real(r8) :: t_xct,  t_cat1, t_cat2
    real(r8) :: t_fac1, t_fac2, t_xfac
    integer  :: t_ifc1, t_ifc2
    real(r8) :: d2mx(2), dxm1(2), invd(2)
    real(r8) :: opt2d(2,2)
    real(r8) :: opt1, opt2, opt
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! Modes 1-3,  SO4(Ait), BC(Ait) and OC(Ait):

    do kcomp=2,3
       call this%zero(kcomp, ncol)
    end do ! kcomp

    kcomp = 1
    do k=1,pver
       do icol=1,ncol
          if(Nnatk(icol,k,kcomp)>0.0_r8) then
             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_xct  = xct(icol,k,kcomp)
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_xfac = xfac(icol,k,kcomp)

             ! partial lengths along each dimension (1-2) for interpolation
             d2mx(1) = (t_cat2-t_xct)
             dxm1(1) = (t_xct-t_cat1)
             invd(1) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(2) = (t_fac2-t_xfac)
             dxm1(2) = (t_xfac-t_fac1)
             invd(2) = 1.0_r8/(t_fac2-t_fac1)

             do iv=1,19  ! variable number

                ! end points as basis for multidimentional linear interpolation
                opt2d(1,1) = a2to3var(iv,t_ict1,t_ifc1,kcomp)
                opt2d(1,2) = a2to3var(iv,t_ict1,t_ifc2,kcomp)
                opt2d(2,1) = a2to3var(iv,t_ict2,t_ifc1,kcomp)
                opt2d(2,2) = a2to3var(iv,t_ict2,t_ifc2,kcomp)

                ! interpolation in the fac dimension
                opt1 = (d2mx(2)*opt2d(1,1)+dxm1(2)*opt2d(1,2))*invd(2)
                opt2 = (d2mx(2)*opt2d(2,1)+dxm1(2)*opt2d(2,2))*invd(2)

                ! finally, interpolation in the cat dimension
                opt = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

                call this%update(kcomp, k, icol, iv, opt)
             end do
          end if
       end do
    end do

    ! Dry parameters for externally mixed modes modes 12-13,
    ! BC(n) and OC(n):

    do kcomp=12,13    ! using dummy initialization for kcomp=3
       do k=1,pver
          do icol=1,ncol
             this%cknorm(icol,k,kcomp)  = a2to3var(1,1,1,kcomp-10)
             this%cknlt05(icol,k,kcomp) = a2to3var(2,1,1,kcomp-10)
             this%ckngt125(icol,k,kcomp)= a2to3var(3,1,1,kcomp-10)
             this%aaerosn(icol,k,kcomp) = a2to3var(16,1,1,kcomp-10)
             this%aaeroln(icol,k,kcomp) = a2to3var(17,1,1,kcomp-10)
             this%vaerosn(icol,k,kcomp) = a2to3var(18,1,1,kcomp-10)
             this%vaeroln(icol,k,kcomp) = a2to3var(19,1,1,kcomp-10)
          end do ! icol
       end do ! k
    end do  ! kcomp

  end subroutine intdrypar2to3

  ! ==========================================================
  subroutine intdrypar4 (this, lchnk, ncol, Nnatk, xfbcbg, ifbcbg1, &
       xfbcbgn, ifbcbgn1, xct, ict1, xfac, ifac1, xfaq, ifaq1)

    ! Output arguments: Modal mass concentrations (cint), area (aaero)
    ! and volume (vaero) (for AeroCom determination of particle
    ! effective radii) of each constituent. cint*05 and cint*125 are
    ! for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
    ! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.

    ! arguments
    class(aerodry_prop_type) :: this
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! condensed SOA/(SOA+H2SO4) (1-4) or added carbonaceous fraction (5-10)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfbcbg(pcols,pver)         ! mass fraction BC/(BC+OC) for the background mode (4)
    integer  , intent(in) :: ifbcbg1(pcols,pver)
    real(r8) , intent(in) :: xfbcbgn(pcols,pver)        ! mass fraction BC/(BC+OC) for the background mode (14)
    integer  , intent(in) :: ifbcbgn1(pcols,pver)

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: iv, kcomp, k, icol
    integer  :: t_ifb1, t_ifb2
    integer  :: t_ict1, t_ict2, t_ifc1, t_ifc2, t_ifa1, t_ifa2
    real(r8) :: t_fbcbg1, t_fbcbg2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xct,  t_cat1, t_cat2
    real(r8) :: t_xfbcbg
    real(r8) :: d2mx(4), dxm1(4), invd(4)
    real(r8) :: opt4d(2,2,2,2)
    real(r8) :: opt1, opt2, opt
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! Mode 4, BC&OC(Ait):
    kcomp=4
    call this%zero(kcomp, ncol)

    do k=1,pver
       do icol=1,ncol
          if(Nnatk(icol,k,kcomp)>0.0_r8) then
             ! Collect all the vector elements into temporary storage
             ! to avoid cache conflicts and excessive cross-referencing
             t_ifb1 = ifbcbg1(icol,k)
             t_ifb2 = t_ifb1+1
             t_ict1 = ict1(icol,k,kcomp)
             t_ict2 = t_ict1+1
             t_ifc1 = ifac1(icol,k,kcomp)
             t_ifc2 = t_ifc1+1
             t_ifa1 = ifaq1(icol,k,kcomp)
             t_ifa2 = t_ifa1+1
             t_fbcbg1 = fbcbg(t_ifb1)
             t_fbcbg2 = fbcbg(t_ifb2)
             t_cat1 = cate(kcomp,t_ict1)
             t_cat2 = cate(kcomp,t_ict2)
             t_fac1 = fac(t_ifc1)
             t_fac2 = fac(t_ifc2)
             t_faq1 = faq(t_ifa1)
             t_faq2 = faq(t_ifa2)
             t_xfbcbg = xfbcbg(icol,k)
             t_xct  = xct(icol,k,kcomp)
             t_xfac = xfac(icol,k,kcomp)
             t_xfaq = xfaq(icol,k,kcomp)

             ! partial lengths along each dimension (1-5) for interpolation
             d2mx(1) = (t_fbcbg2-t_xfbcbg)
             dxm1(1) = (t_xfbcbg-t_fbcbg1)
             invd(1) = 1.0_r8/(t_fbcbg2-t_fbcbg1)
             d2mx(2) = (t_cat2-t_xct)
             dxm1(2) = (t_xct-t_cat1)
             invd(2) = 1.0_r8/(t_cat2-t_cat1)
             d2mx(3) = (t_fac2-t_xfac)
             dxm1(3) = (t_xfac-t_fac1)
             invd(3) = 1.0_r8/(t_fac2-t_fac1)
             d2mx(4) = (t_faq2-t_xfaq)
             dxm1(4) = (t_xfaq-t_faq1)
             invd(4) = 1.0_r8/(t_faq2-t_faq1)

             do iv=1,19  ! variable number

                ! end points as basis for multidimentional linear interpolation
                opt4d(1,1,1,1)=a4var(iv,t_ifb1,t_ict1,t_ifc1,t_ifa1)
                opt4d(1,1,1,2)=a4var(iv,t_ifb1,t_ict1,t_ifc1,t_ifa2)
                opt4d(1,1,2,1)=a4var(iv,t_ifb1,t_ict1,t_ifc2,t_ifa1)
                opt4d(1,1,2,2)=a4var(iv,t_ifb1,t_ict1,t_ifc2,t_ifa2)
                opt4d(1,2,1,1)=a4var(iv,t_ifb1,t_ict2,t_ifc1,t_ifa1)
                opt4d(1,2,1,2)=a4var(iv,t_ifb1,t_ict2,t_ifc1,t_ifa2)
                opt4d(1,2,2,1)=a4var(iv,t_ifb1,t_ict2,t_ifc2,t_ifa1)
                opt4d(1,2,2,2)=a4var(iv,t_ifb1,t_ict2,t_ifc2,t_ifa2)
                opt4d(2,1,1,1)=a4var(iv,t_ifb2,t_ict1,t_ifc1,t_ifa1)
                opt4d(2,1,1,2)=a4var(iv,t_ifb2,t_ict1,t_ifc1,t_ifa2)
                opt4d(2,1,2,1)=a4var(iv,t_ifb2,t_ict1,t_ifc2,t_ifa1)
                opt4d(2,1,2,2)=a4var(iv,t_ifb2,t_ict1,t_ifc2,t_ifa2)
                opt4d(2,2,1,1)=a4var(iv,t_ifb2,t_ict2,t_ifc1,t_ifa1)
                opt4d(2,2,1,2)=a4var(iv,t_ifb2,t_ict2,t_ifc1,t_ifa2)
                opt4d(2,2,2,1)=a4var(iv,t_ifb2,t_ict2,t_ifc2,t_ifa1)
                opt4d(2,2,2,2)=a4var(iv,t_ifb2,t_ict2,t_ifc2,t_ifa2)

                ! interpolation in the faq, fac and cat dimensions
                call lininterpol4dim (d2mx, dxm1, invd, opt4d, opt1, opt2)

                ! finally, interpolation in the fbcbg dimension
                opt = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

                call this%update(kcomp, k, icol, iv, opt)
             end do
          endif
       end do ! icol
    end do ! k

    kcomp=14
    do k=1,pver
       do icol=1,ncol

          t_ifb1 = ifbcbgn1(icol,k)
          t_ifb2 = t_ifb1+1
          t_fbcbg1 = fbcbg(t_ifb1)
          t_fbcbg2 = fbcbg(t_ifb2)
          t_xfbcbg = xfbcbgn(icol,k)

          d2mx(1) = (t_fbcbg2-t_xfbcbg)
          dxm1(1) = (t_xfbcbg-t_fbcbg1)
          invd(1) = 1.0_r8/(t_fbcbg2-t_fbcbg1)

          !  Only interpolation in the fbcbg dimension for mode 14
          opt1 = a4var(1,1,1,1,1)
          opt2 = a4var(1,2,1,1,1)
          this%cknorm(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          opt1 = a4var(2,1,1,1,1)
          opt2 = a4var(2,2,1,1,1)
          this%cknlt05(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          opt1 = a4var(3,1,1,1,1)
          opt2 = a4var(3,2,1,1,1)
          this%ckngt125(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          ! (The remaining variables are actually independent of fbcbg,
          ! but we follow the same procedure anyway:)

          opt1 = a4var(16,1,1,1,1)
          opt2 = a4var(16,2,1,1,1)
          this%aaerosn(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          opt1 = a4var(17,1,1,1,1)
          opt2 = a4var(17,2,1,1,1)
          this%aaeroln(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          opt1 = a4var(18,1,1,1,1)
          opt2 = a4var(18,2,1,1,1)
          this%vaerosn(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

          opt1 = a4var(19,1,1,1,1)
          opt2 = a4var(19,2,1,1,1)
          this%vaeroln(icol,k,kcomp) = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

       end do ! icol
    end do ! k

  end subroutine intdrypar4

  ! ==========================================================
  subroutine intdrypar5to10 (this, lchnk, ncol, Nnatk, xct, ict1, &
       xfac, ifac1, xfbc, ifbc1, xfaq, ifaq1)

    ! Output arguments: Modal mass concentrations (cint), area (aaero)
    ! and volume (vaero) (for AeroCom determination of particle
    ! effective radii) of each constituent. cint*05 and cint*125 are
    ! for r<0.5um and r>1.25um, respectively. aaeros and vaeros are
    ! integrated over r<0.5um, and aaerol and vaerol over r>0.5um.

    ! arguments
    class(aerodry_prop_type) :: this
    integer  , intent(in) :: lchnk                      ! chunk identifier
    integer  , intent(in) :: ncol                       ! number of atmospheric columns
    real(r8) , intent(in) :: Nnatk(pcols,pver,0:nmodes) ! modal aerosol number concentration
    real(r8) , intent(in) :: xct(pcols,pver,nmodes)     ! modal internally mixed SO4+BC+OC conc.
    integer  , intent(in) :: ict1(pcols,pver,nmodes)
    real(r8) , intent(in) :: xfac(pcols,pver,nbmodes)   ! modal (OC+BC)/(SO4+BC+OC)
    integer  , intent(in) :: ifac1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfbc(pcols,pver,nbmodes)   ! modal BC/(OC+BC)
    integer  , intent(in) :: ifbc1(pcols,pver,nbmodes)
    real(r8) , intent(in) :: xfaq(pcols,pver,nbmodes)   ! modal SO4(aq)/SO4
    integer  , intent(in) :: ifaq1(pcols,pver,nbmodes)

    ! local variables
    real(r8) :: a, b, e, eps
    integer  :: iv, kcomp, k, icol
    integer  :: t_ict1, t_ict2, t_ifa1, t_ifa2
    integer  :: t_ifb1, t_ifb2, t_ifc1, t_ifc2
    real(r8) :: t_faq1, t_faq2, t_xfaq
    real(r8) :: t_fbc1, t_fbc2, t_xfbc
    real(r8) :: t_fac1, t_fac2, t_xfac
    real(r8) :: t_xct,  t_cat1, t_cat2
    real(r8) :: d2mx(4), dxm1(4), invd(4)
    real(r8) :: opt4d(2,2,2,2)
    real(r8) :: opt1, opt2, opt
    parameter (e=2.718281828_r8, eps=1.0e-60_r8)

    ! Modes 5 to 10 (SO4(Ait75) and mineral and seasalt-modes + cond./coag./aq.):

    do kcomp=5,10
       call this%zero(kcomp,ncol)

       do k=1,pver
          do icol=1,ncol
             if(Nnatk(icol,k,kcomp)>0.0_r8) then
                ! Collect all the vector elements into temporary storage
                ! to avoid cache conflicts and excessive cross-referencing
                t_ict1 = ict1(icol,k,kcomp)
                t_ict2 = t_ict1+1
                t_ifc1 = ifac1(icol,k,kcomp)
                t_ifc2 = t_ifc1+1
                t_ifb1 = ifbc1(icol,k,kcomp)
                t_ifb2 = t_ifb1+1
                t_ifa1 = ifaq1(icol,k,kcomp)
                t_ifa2 = t_ifa1+1
                t_cat1 = cat(kcomp,t_ict1)
                t_cat2 = cat(kcomp,t_ict2)
                t_fac1 = fac(t_ifc1)
                t_fac2 = fac(t_ifc2)
                t_fbc1 = fbc(t_ifb1)
                t_fbc2 = fbc(t_ifb2)
                t_faq1 = faq(t_ifa1)
                t_faq2 = faq(t_ifa2)
                t_xct  = xct(icol,k,kcomp)
                t_xfac = xfac(icol,k,kcomp)
                t_xfbc = xfbc(icol,k,kcomp)
                t_xfaq = xfaq(icol,k,kcomp)

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
                !soa

                do iv=1,19  ! variable number
                   ! end points as basis for multidimentional linear interpolation
                   opt4d(1,1,1,1)=a5to10var(iv,t_ict1,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt4d(1,1,1,2)=a5to10var(iv,t_ict1,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt4d(1,1,2,1)=a5to10var(iv,t_ict1,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt4d(1,1,2,2)=a5to10var(iv,t_ict1,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt4d(1,2,1,1)=a5to10var(iv,t_ict1,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt4d(1,2,1,2)=a5to10var(iv,t_ict1,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt4d(1,2,2,1)=a5to10var(iv,t_ict1,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt4d(1,2,2,2)=a5to10var(iv,t_ict1,t_ifc2,t_ifb2,t_ifa2,kcomp)
                   opt4d(2,1,1,1)=a5to10var(iv,t_ict2,t_ifc1,t_ifb1,t_ifa1,kcomp)
                   opt4d(2,1,1,2)=a5to10var(iv,t_ict2,t_ifc1,t_ifb1,t_ifa2,kcomp)
                   opt4d(2,1,2,1)=a5to10var(iv,t_ict2,t_ifc1,t_ifb2,t_ifa1,kcomp)
                   opt4d(2,1,2,2)=a5to10var(iv,t_ict2,t_ifc1,t_ifb2,t_ifa2,kcomp)
                   opt4d(2,2,1,1)=a5to10var(iv,t_ict2,t_ifc2,t_ifb1,t_ifa1,kcomp)
                   opt4d(2,2,1,2)=a5to10var(iv,t_ict2,t_ifc2,t_ifb1,t_ifa2,kcomp)
                   opt4d(2,2,2,1)=a5to10var(iv,t_ict2,t_ifc2,t_ifb2,t_ifa1,kcomp)
                   opt4d(2,2,2,2)=a5to10var(iv,t_ict2,t_ifc2,t_ifb2,t_ifa2,kcomp)

                   ! interpolation in the faq, fbc, and fac and dimensions
                   call lininterpol4dim (d2mx, dxm1, invd, opt4d, opt1, opt2)

                   ! finally, interpolation in the cat dimension
                   opt = (d2mx(1)*opt1+dxm1(1)*opt2)*invd(1)

                   call this%update(kcomp, k, icol, iv, opt)
                end do
             endif

             this%cknorm(icol,k,kcomp)  = a5to10var(1,1,1,1,1,kcomp)
             this%cknlt05(icol,k,kcomp) = a5to10var(2,1,1,1,1,kcomp)
             this%ckngt125(icol,k,kcomp)= a5to10var(3,1,1,1,1,kcomp)

          end do ! icol
       end do ! k
    end do  ! kcomp

  end subroutine intdrypar5to10

  ! ==========================================================
  subroutine zero(this, kcomp, ncol)

    class(aerodry_prop_type) :: this
    integer , intent(in)  :: kcomp
    integer , intent(in)  :: ncol

    integer :: k
    integer :: icol

    ! initialize all output fields to zero
    do k=1,pver
       do icol=1,ncol
          this%cintbg(icol,k,kcomp)    = 0.0_r8
          this%cintbg05(icol,k,kcomp)  = 0.0_r8
          this%cintbg125(icol,k,kcomp) = 0.0_r8
          this%cintbc(icol,k,kcomp)    = 0.0_r8
          this%cintbc05(icol,k,kcomp)  = 0.0_r8
          this%cintbc125(icol,k,kcomp) = 0.0_r8
          this%cintoc(icol,k,kcomp)    = 0.0_r8
          this%cintoc05(icol,k,kcomp)  = 0.0_r8
          this%cintoc125(icol,k,kcomp) = 0.0_r8
          this%cintsc(icol,k,kcomp)    = 0.0_r8
          this%cintsc05(icol,k,kcomp)  = 0.0_r8
          this%cintsc125(icol,k,kcomp) = 0.0_r8
          this%cintsa(icol,k,kcomp)    = 0.0_r8
          this%cintsa05(icol,k,kcomp)  = 0.0_r8
          this%cintsa125(icol,k,kcomp) = 0.0_r8
          this%aaeros(icol,k,kcomp)    = 0.0_r8
          this%aaerol(icol,k,kcomp)    = 0.0_r8
          this%vaeros(icol,k,kcomp)    = 0.0_r8
          this%vaerol(icol,k,kcomp)    = 0.0_r8
       end do
    end do
  end subroutine zero

  ! ==========================================================
  subroutine update(this, kcomp, k, icol, iv, opt)

    class(aerodry_prop_type) :: this
    integer , intent(in)  :: kcomp
    integer , intent(in)  :: k
    integer , intent(in)  :: icol
    integer , intent(in)  :: iv
    real(r8), intent(in)  :: opt

    if(iv==1) then
       this%cintbg(icol,k,kcomp)=opt
    elseif(iv==2) then
       this%cintbg05(icol,k,kcomp)=opt
    elseif(iv==3) then
       this%cintbg125(icol,k,kcomp)=opt
    elseif(iv==4) then
       this%cintbc(icol,k,kcomp)=opt
    elseif(iv==5) then
       this%cintbc05(icol,k,kcomp)=opt
    elseif(iv==6) then
       this%cintbc125(icol,k,kcomp)=opt
    elseif(iv==7) then
       this%cintoc(icol,k,kcomp)=opt
    elseif(iv==8) then
       this%cintoc05(icol,k,kcomp)=opt
    elseif(iv==9) then
       this%cintoc125(icol,k,kcomp)=opt
    elseif(iv==10) then
       this%cintsc(icol,k,kcomp)=opt
    elseif(iv==11) then
       this%cintsc05(icol,k,kcomp)=opt
    elseif(iv==12) then
       this%cintsc125(icol,k,kcomp)=opt
    elseif(iv==13) then
       this%cintsa(icol,k,kcomp)=opt
    elseif(iv==14) then
       this%cintsa05(icol,k,kcomp)=opt
    elseif(iv==15) then
       this%cintsa125(icol,k,kcomp)=opt
    elseif(iv==16) then
       this%aaeros(icol,k,kcomp)=opt
    elseif(iv==17) then
       this%aaerol(icol,k,kcomp)=opt
    elseif(iv==18) then
       this%vaeros(icol,k,kcomp)=opt
    elseif(iv==19) then
       this%vaerol(icol,k,kcomp)=opt
    endif

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

end module oslo_aero_aerocom_dry
