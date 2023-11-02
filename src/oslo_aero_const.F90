module oslo_aero_const

  !-----------------------------------------------------------------------------
  ! Module containing oslo_aero constants
  !-----------------------------------------------------------------------------

  use shr_kind_mod,      only: r8 => shr_kind_r8
  use oslo_aero_params,  only: nmodes
  use physconst,         only: pi
  !
  implicit none
  public

  public :: init_interp_constants

  real(r8), parameter :: smallNumber = 1.e-100_r8
  real(r8), parameter :: rTabMin = 1.e-9_r8               ![m] smallest lookup table size
  real(r8), parameter :: rTabMax = 20.e-6_r8              ![m] largest lookup table size
  integer,  parameter :: nBinsTab = 44                    ![nbr] number of tabulated bins
  real(r8), parameter :: rMinAquousChemistry = 0.05e-6_r8 ! Smallest particle which can receive aquous chemistry mass
  real(r8), parameter :: sq2pi = 1._r8/sqrt(2.0_r8*pi)

  real(r8) :: nk(0:nmodes,nbinsTab)                       !dN/dlogr for modes
  real(r8) :: normnk(0:nmodes,nbinsTab)                   !dN for modes (sums to one over size range)
  real(r8) :: rBinEdge(nBinsTab+1)
  real(r8) :: rBinMidpoint(nBinsTab)
  real(r8) :: volumeToNumber(0:nmodes)                    !m3 ==> #
  real(r8) :: numberToSurface(0:nmodes)                   !# ==> m2

  ! Constants used in interpolation
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

  real(r8) :: rh(10)
  real(r8) :: fombg(6), fbcbg(6), fac(6), fbc(6), faq(6)
  real(r8) :: cate(4,16)
  real(r8) :: cat(5:10,6)

  ! relative humidity (RH, as integer for output variable names) for use in AeroCom code
  integer  :: RF(6)

  ! AeroCom specific RH input variables for use in opticsAtConstRh.F90
  integer  :: irhrf1(6)
  real(r8) :: xrhrf(6)

  real(r8), parameter :: e=2.718281828_r8
  real(r8), parameter :: eps=1.0e-30_r8  

!=============================================================================
contains
!=============================================================================

  subroutine init_interp_constants()

    !---------------------------------------------------------------
    ! set module variables
    !---------------------------------------------------------------

    ! Local variables
    integer :: irf, irelh, kcomp, i
    !-----------------------------------------------------------

    ! Defining array bounds for tabulated optical parameters (and r and sigma)
    ! relative humidity (only 0 value used for r and sigma tables):
    rh = (/ 0.0_r8, 0.37_r8, 0.47_r8, 0.65_r8, 0.75_r8, 0.8_r8, 0.85_r8, 0.9_r8, 0.95_r8, 0.995_r8 /)

    ! relative humidity (RH, as integer for output variable names) for use in AeroCom code
    RF = (/0, 40, 55, 65, 75, 85 /)

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
    ! and most weight on small concentrations for added mass onto the background:

    fbcbg(1)=1.e-10_r8
    fbc(1)=1.e-10_r8
    do i=2,6
       fbcbg(i)=10**((i-1)/4.0_r8-1.25_r8)
       fbc(i)=fbcbg(i)
    end do

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

  end subroutine init_interp_constants

end module oslo_aero_const




