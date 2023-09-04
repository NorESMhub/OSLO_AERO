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

end module oslo_aero_const




