module oslo_aero_seasalt

  !-----------------------------------------------------------------------
  ! compute emission of sea salt
  !-----------------------------------------------------------------------

  use shr_kind_mod,    only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,          only: pcols, pver
  use constituents,    only: cnst_name
  use camsrfexch,      only: cam_in_t
  use physics_types,   only: physics_state
  !
  use oslo_aero_const, only: volumeToNumber
  use oslo_aero_ocean, only: oslo_aero_opom_inq, oslo_aero_opom_emis
  use oslo_aero_share, only: rhopart, l_om_ni, l_ss_a1, l_ss_a2, l_ss_a3
  use oslo_aero_share, only: MODE_IDX_SS_A1, MODE_IDX_SS_A2, MODE_IDX_SS_A3

  implicit none
  private

  integer , parameter :: numberOfSaltModes = 3

  character(len=6)  , public :: seasalt_names(10)
  integer, parameter, public :: seasalt_nbin = numberOfSaltModes  ! needed by mo_photo.F90
  logical, parameter, public :: seasalt_active = .true.

  integer :: modeMap(numberOfSaltModes)    ! [idx] which modes are we modifying
  integer :: tracerMap(numberOfSaltModes)  ! [idx] which tracers are we modifying

  public :: oslo_aero_seasalt_init
  public :: oslo_aero_seasalt_emis

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_seasalt_init()

    integer :: i

    modeMap(1) = MODE_IDX_SS_A1
    modeMap(2) = MODE_IDX_SS_A2
    modeMap(3) = MODE_IDX_SS_A3

    tracerMap(1) = l_ss_a1
    tracerMap(2) = l_ss_a2
    tracerMap(3) = l_ss_a3

    seasalt_names(:) = "      "
    do i = 1,numberOfSaltModes
       seasalt_names(i) = cnst_name(tracerMap(i))
    end do

  end subroutine oslo_aero_seasalt_init

  !===============================================================================
  subroutine oslo_aero_seasalt_emis(state, cam_in)

    ! Arguments:
    type(physics_state),    intent(in)    :: state   ! Physics state variables
    type(cam_in_t), target, intent(inout) :: cam_in  ! import state

    ! Local variables
    integer             :: n                                   ![] counter for modes
    integer             :: ncol                                ![nbr] number of columns in use
    integer             :: lchnk                               !chunk index
    real(r8)            :: whiteCapAreaFraction(pcols)         ![fraction]
    real(r8)            :: open_ocean(pcols)                   ![fraction]
    real(r8)            :: numberFlux(pcols,numberofSaltModes) ![#/m2/sec]
    real(r8)            :: u10m(pcols)                         ![m/s]
    real(r8), pointer   :: sst(:)                              ![frc] sea surface temperature
    real(r8), pointer   :: ocnfrc(:)                           ![frc] ocean fraction
    real(r8), pointer   :: icefrc(:)                           ![frc] ice fraction
    real(r8)            :: spracklenOMOceanSource(pcols)       ![kg/m2/s] spracklen ocean source
    real(r8)            :: onOMOceanSource(pcols)              ![kg/m2/s] OM source from Nilsson/O'Dowd
    real(r8)            :: OMOceanSource(pcols)                ![kg/m2/s] new OM ocean source
    real(r8), parameter :: z0= 0.0001_r8                       ![m] roughness length over ocean

    !New numbers are based on Salter et al. (2105):
    !www.atmos-chem-phys-discuss.net/15/13783/2015/doi:10.5194/acpd-15-13783-2015
    !Values from Table 1 in Salter et al. (2015):
    real(r8), parameter :: coeffA(numberOfSaltModes) = (/-5.2168e5_r8,  0.0_r8,      0.0_r8      /)
    real(r8), parameter :: coeffB(numberOfSaltModes) = (/ 3.31725e7_r8, 7.374e5_r8,  1.4210e4_r8 /)
    real(r8), parameter :: coeffC(numberOfSaltModes) = (/-6.95275e8_r8,-2.4803e7_r8, 1.4662e7_r8 /)
    real(r8), parameter :: coeffD(numberOfSaltModes) = (/ 1.0684e10_r8, 7.7373e8_r8, 1.7075e8_r8 /)

    !After discussions with Alf K, it is better to scale with only smallest SS-mode since POM is small
    !and assume same production mechanism. Nudged 1 degree simulations give 2.52 Tg/yr of SS_A1, so
    !to obtain 7.7, we need to scale them by 7.7 / 2.52 ==> 3.03
    !updated value for Salter et al. sea-salt treatment, which gives global annual SS_A1 emissions of
    !2.663 instead of 0.153 ng m-2 s-1 (i.e. ca 17 times more than the old sea-salt treatment):
    real(r8), parameter :: seasaltToSpracklenOM2 = 3.03_r8*0.153_r8/2.663_r8

    !number of columns in use
    ncol = state%ncol
    lchnk = state%lchnk

    !pointers to land model variables
    ocnfrc => cam_in%ocnfrac
    icefrc => cam_in%icefrac
    sst    => cam_in%sst

    !start with midpoint wind speed
    u10m(:ncol)=sqrt(state%u(:ncol,pver)**2+state%v(:ncol,pver)**2)

    ! move the winds to 10m high from the midpoint of the gridbox:
    u10m(:ncol)=u10m(:ncol)*log(10._r8/z0)/log(state%zm(:ncol,pver)/z0)

    ! New whitecap area fraction / air entrainment flux from eqn. 6 in Salter et al. (2015)
    ! JCA & MS Using Hanson & Phillips 99 air entrainment vs. wind speed
    ! (Note the uncertainty in the factor 2, written as 2 pluss/minus 1 in Eq. 6 -> possible tuning factor)
    whitecapAreaFraction(:ncol) = (2.0_r8*10.0_r8**(-8.0_r8))*(u10m(:ncol)**3.74_r8)
    whitecapAreaFraction(:ncol) = ocnfrc(:ncol) * (1._r8-icefrc(:ncol)) * whitecapAreaFraction(:ncol)

    ! Determine open ocean fraction on gridcell
    open_ocean(:ncol) = ocnfrc(:ncol) * (1._r8-icefrc(:ncol))

    ! Eqn. 9 in Salter et al. (2015)
    do n=1,numberOfSaltModes
       numberFlux(:ncol,n) = whitecapAreaFraction(:ncol)*                                    &
            ( coeffA(n)*(sst(:ncol)-273.15_r8)*(sst(:ncol)-273.15_r8)*(sst(:ncol)-273.15_r8) &
            + coeffB(n)*(sst(:ncol)-273.15_r8)*(sst(:ncol)-273.15_r8)                        &
            + coeffC(n)*(sst(:ncol)-273.15_r8)                                               &
            + coeffD(n) )
    end do

    do n=1,numberOfSaltModes
       cam_in%cflx(:ncol, tracerMap(n)) = numberFlux(:ncol,n) & !#/m2/sec
            / volumeToNumber(modeMap(n))                      & !==> m3/m2/sec
            * rhopart(tracerMap(n))                             !==> kg/m2/sec
    end do
    spracklenOMOceanSource(:ncol) = cam_in%cflx(:ncol, tracerMap(1))*seasaltToSpracklenOM2

    if (oslo_aero_opom_inq())then
       call oslo_aero_opom_emis(cam_in%cflx(:ncol, tracerMap(1)), &
            cam_in%cflx(:ncol,tracerMap(2)), cam_in%cflx(:ncol,tracerMap(3)), &
            open_ocean, ncol, lchnk, onOMOceanSource )
       OMOceanSource(:ncol) = onOMOceanSource(:ncol)
    else
       OMOceanSource(:ncol) = spracklenOMOceanSource(:ncol)
    endif

    !Add OM ocean source to cam_in
    cam_in%cflx(:ncol,l_om_ni) = cam_in%cflx(:ncol,l_om_ni) + OMOceanSource(:ncol)

  end subroutine oslo_aero_seasalt_emis

end module oslo_aero_seasalt
