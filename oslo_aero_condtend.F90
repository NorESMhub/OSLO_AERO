module oslo_aero_condtend

  ! Calculate the sulphate nucleation rate, and condensation rate of
  ! aerosols used for parameterising the transfer of externally mixed
  ! aitken mode particles into an internal mixture.
  ! Note the parameterisation for conversion of externally mixed particles
  ! used the h2so4 lifetime onto the particles, and not a given
  ! increase in particle radius. Will be improved in future versions of the model
  
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use ppgrid,             only: pcols, pver, pverp
  use phys_control,       only: phys_getopts
  use chem_mods,          only: gas_pcnst
  use mo_tracname,        only: solsym
  use cam_history,        only: addfld, add_default, fieldname_len, horiz_only, outfld
  use physconst,          only: rair, gravit, pi, avogad
  use chem_mods,          only: adv_mass !molecular weights from mozart
  use wv_saturation,      only: qsat_water
  use m_spc_id,           only: id_H2SO4, id_soa_lv
  !
  use oslo_aero_coag,     only: normalizedCoagulationSink, receiverMode,numberOfCoagulationReceivers
  use oslo_aero_coag,     only: numberOfAddCoagReceivers,addReceiverMode,normCoagSinkAdd
  use constituents,       only: pcnst  ! h2so4 and soa nucleation (cka)
  use oslo_aero_share   ! only: MODE_IDX_SO4SOA_AIT, rhopart, l_so4_a1, l_soa_lv, l_so4_na, l_soa_na
  use oslo_aero_params  ! only: originalNumberMedianRadius
  use oslo_aero_const   ! only: volumeToNumber

  implicit none
  private

  ! public routines
  public :: registerCondensation
  public :: initializeCondensation
  public :: condtend

  ! private routines
  private :: aeronucl
  private :: appformrate

  integer, parameter, public :: N_COND_VAP = 3
  integer, parameter, public :: COND_VAP_H2SO4 = 1
  integer, parameter, public :: COND_VAP_ORG_LV = 2
  integer, parameter, public :: COND_VAP_ORG_SV = 3

  real(r8), public  :: normalizedCondensationSink(0:nmodes,N_COND_VAP) ! [m3/#/s] condensation sink per particle in mode i

  integer , private :: lifeCycleReceiver(gas_pcnst)                    ! [-] array of transformation of life cycle tracers
  real(r8), private :: stickingCoefficient(0:nmodes,N_COND_VAP)        ! [-] stickingCoefficient for H2SO4 on a mode
  integer , private :: cond_vap_map(N_COND_VAP)

  ! Assumed number of monolayers
  real(r8), parameter, private :: n_so4_monolayers_age = 3.0_r8

  ! thickness of the so4 monolayers (m)
  ! for so4(+nh4), use bi-sulfate mw and 1.77 g/cm3 as in MAM
  real(r8), parameter, public :: dr_so4_monolayers_age = n_so4_monolayers_age * 4.76e-10_r8

!===============================================================================
contains
!===============================================================================

  subroutine registerCondensation()

    ! local variables
    integer :: iDonor
    integer :: l_donor
    integer :: tracerIndex
    integer :: mode_index_donor

    !These are the lifecycle-species which receive mass when
    !the externally mixed modes receive condensate,
    !e.g. the receiver of l_so4_n mass is the tracer l_so4_na
    lifeCycleReceiver(:) = -99
    lifeCycleReceiver(chemistryIndex(l_bc_n))   = chemistryIndex(l_bc_a)    !create bc int mix from bc in mode 12
    lifeCycleReceiver(chemistryIndex(l_bc_ni))  = chemistryIndex(l_bc_ai)   !create bc int mix from bc in mode 14
    lifeCycleReceiver(chemistryIndex(l_om_ni))  = chemistryIndex(l_om_ai)

    !!create om int mix from om in mode 14
    lifeCycleReceiver(chemistryIndex(l_bc_ax))  = chemistryIndex(l_bc_ai)
    !!create bc int mix from bc in mode 0. Note Mass is conserved but not number

    !Sticking coeffcients for H2SO4 condensation
    !See table 1 in Kirkevag et al (2013)
    !http://www.geosci-model-dev.net/6/207/2013/gmd-6-207-2013.html
    !Note: In NorESM1, sticking coefficients of the externally mixed modes were
    !used for the internally mixed modes in modallapp. In condtend the internally
    !mixed modes had sticking coefficient = 1.0
    !This might be correct, but is too confusing, so here just
    !assign based on background aerosol and table 1 in Kirkevag et al
    stickingCoefficient(:,:) = 1.0_r8
    stickingCoefficient(MODE_IDX_BC_EXT_AC,:) = 0.3_r8
    stickingCoefficient(MODE_IDX_BC_AIT,:) = 0.3_r8
    stickingCoefficient(MODE_IDX_OMBC_INTMIX_COAT_AIT,:) = 0.5_r8
    stickingCoefficient(MODE_IDX_DST_A2,:) = 0.3_r8
    stickingCoefficient(MODE_IDX_DST_A3,:) = 0.3_r8
    stickingCoefficient(MODE_IDX_BC_NUC,:) = 0.3_r8
    stickingCoefficient(MODE_IDX_OMBC_INTMIX_AIT,:) = 0.5_r8

  end subroutine registerCondensation

  !===============================================================================

  subroutine initializeCondensation()

    !condensation coefficients:
    !Theory: Poling et al, "The properties of gases and liquids"
    !5th edition, eqn 11-4-4

    ! local variables
    real(r8), parameter :: aunit = 1.6606e-27_r8  ![kg] Atomic mass unit
    real(r8), parameter :: boltz = 1.3806e-23_r8   ![J/K/molec]
    real(r8), parameter :: t0 = 273.15_r8         ![K] standard temperature
    real(r8), parameter :: p0 = 101325.0_r8       ! [Pa] Standard pressure
    real(r8), parameter :: radair = 1.73e-10_r8   ![m] Typical air molecule collision radius
    real(r8), parameter :: Mair = 28.97_r8        ![amu/molec] Molecular weight for dry air
    !Diffusion volumes for simple molecules [Poling et al], table 11-1
    real(r8), dimension(N_COND_VAP), parameter :: vad = (/51.96_r8, 208.18_r8, 208.18_r8/) ![cm3/mol]
    real(r8), parameter :: vadAir       = 19.7_r8                                          ![cm3/mol]
    real(r8), parameter :: aThird = 1.0_r8/3.0_r8
    real(r8), parameter :: cm2Tom2 = 1.e-4_r8       !convert from cm2 ==> m2

    real(r8), dimension(0:100,0:nmodes,N_COND_VAP) :: DiffusionCoefficient   ! [m2/s] Diffusion coefficient
    character(len=fieldname_len+3) :: fieldname_donor
    character(len=fieldname_len+3) :: fieldname_receiver
    character(128)                 :: long_name
    character(8)                   :: unit

    integer                        :: nsiz !counter for aerotab sizes
    integer                        :: iChem             !counter for chemical species
    integer                        :: mode_index_donor  !index for mode
    integer                        :: iMode             !Counter for mode
    integer                        :: tracerIndex       !counter for chem. spec

    logical                        :: history_aerosol
    logical                        :: isAlreadyOnList(gas_pcnst)
    integer                        :: cond_vap_idx

    real(r8), dimension(N_COND_VAP) :: mfv  ![m] mean free path
    real(r8), dimension(N_COND_VAP) :: diff ![m2/s] diffusion coefficient for cond. vap
    real(r8) :: molecularWeight !amu/molec molecular weight
    real(r8) :: Mdual ![molec/amu] 1/M_1 + 1/M_2
    real(r8) :: rho   ![kg/m3] density of component in question
    real(r8) :: radmol ![m] radius molecule
    real(r8), dimension(N_COND_VAP) :: th     !thermal velocity

    !Couple the condenseable vapours to chemical species for properties and indexes
    cond_vap_map(COND_VAP_H2SO4) = chemistryIndex(l_h2so4)
    cond_vap_map(COND_VAP_ORG_LV) = chemistryIndex(l_soa_lv)
    cond_vap_map(COND_VAP_ORG_SV) = chemistryIndex(l_soa_sv)

    do cond_vap_idx = 1, N_COND_VAP

       rho = rhopart(physicsIndex(cond_vap_map(cond_vap_idx))) !pick up densities from oslo_aero_share

       molecularWeight=adv_mass(cond_vap_map(cond_vap_idx))    !pick up molecular weights from mozart

       ! https://en.wikipedia.org/wiki/Thermal_velocity
       th(cond_vap_idx) = sqrt(8.0_r8*boltz*t0/(pi*molecularweight*aunit))   ! thermal velocity for H2SO4 in air (m/s)

       ! Radius of molecul (straight forward assuming spherical)
       radmol=(3.0_r8*molecularWeight*aunit/(4.0_r8*pi*rho))**aThird    ! molecule radius

       Mdual=2.0_r8/(1.0_r8/Mair+1.0_r8/molecularWeight) !factor of [1/m_1 + 1_m2]

       ! calculating microphysical parameters from equations in Ch. 8 of Seinfeld & Pandis (1998):
       ! mean free path for molec in air (m)
       mfv(cond_vap_idx)=1.0_r8/(pi*sqrt(1.0_r8+MolecularWeight/Mair)*(radair+radmol)**2*p0/(boltz*t0)) 

       ! Solve eqn 11-4.4 in Poling et al
       ! (A bit hard to follow units here, but result in the book is in cm2/s)..
       ! so scale by "cm2Tom2" to get m2/sec
       diff(cond_vap_idx) = cm2Tom2*0.00143_r8*t0**1.75_r8/((p0/1.0e5_r8)*sqrt(Mdual)   &
            *(((Vad(cond_vap_idx))**aThird+(Vadair)**aThird)**2))
    end do

    do cond_vap_idx = 1, N_COND_VAP
       do imode = 0, nmodes         !all modes receive condensation
          do nsiz = 1, nBinsTab     !aerotab sizes
             !Correct for non-continuum effects, formula is from
             !Chuang and Penner, Tellus, 1995, sticking coeffient from
             !Vignati et al, JGR, 2004
             !fxm: make "diff ==> diff (cond_vap_idx)
             DiffusionCoefficient(nsiz,imode,cond_vap_idx) = diff(cond_vap_idx)  &    !original diffusion coefficient
                  /(                                    &
                  rBinMidPoint(nsiz)/(rBinMidPoint(nsiz)+mfv(cond_vap_idx))  &  !non-continuum correction factor
                  +4.0_r8*diff(cond_vap_idx)/(stickingCoefficient(imode,cond_vap_idx)*th(cond_vap_idx)*rBinMidPoint(nsiz)) &
                  )
          enddo
       end do !receiver modes
    end do

    normalizedCondensationSink(:,:) = 0.0_r8
    !Find sink per particle in mode "imode"
    !Eqn 13 in Kulmala et al, Tellus 53B, 2001, pp 479
    !http://onlinelibrary.wiley.com/doi/10.1034/j.1600-0889.2001.530411.x/abstract
    do cond_vap_idx =1, N_COND_VAP
       do imode = 0, nmodes
          do nsiz = 1, nBinsTab
             normalizedCondensationSink(imode,cond_vap_idx) =     &
                  normalizedCondensationSink(imode,cond_vap_idx)  &
                  + 4.0_r8*pi                                     &
                  * DiffusionCoefficient(nsiz,imode,cond_vap_idx) &    ![m2/s] diffusion coefficient
                  * rBinMidPoint(nsiz)                            &    ![m] look up table radius
                  * normnk(imode,nsiz)                                 ![frc]
          end do
       end do
    end do

    !Initialize output
    call phys_getopts(history_aerosol_out = history_aerosol)

    isAlreadyOnList(:) = .FALSE.
    do iChem = 1,gas_pcnst
       !Does this tracer have a receiver? If yes: It participate in condensation tendencies
       if(lifeCycleReceiver(iChem) .gt. 0)then
          unit = "kg/m2/s"
          fieldname_donor = trim(solsym(iChem))//"condTend"
          fieldname_receiver = trim(solsym(lifeCycleReceiver(iChem)))//"condTend"
          if(.not. isAlreadyOnList(lifeCycleReceiver(iChem)))then
             call addfld( fieldname_receiver, horiz_only, "A", unit, "condensation tendency" )
             isAlreadyOnList(lifeCycleReceiver(iChem))=.TRUE.
          end if
          call addfld( fieldname_donor, horiz_only, "A", unit, "condensation tendency" )
          if(history_aerosol)then
             call add_default( fieldname_receiver, 1, ' ' )
             call add_default( fieldname_donor   , 1, ' ')
          end if
       end if
    end do

    !Need to add so4_a1, soa_na, so4_na, soa_a1 also (which are not parts of the donor-receiver stuff)
    fieldname_receiver = trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
    call addfld( fieldname_receiver, horiz_only, 'A', unit, "condensation tendency")
    if(history_aerosol)then
       call add_default( fieldname_receiver, 1, ' ' )
    end if

    fieldname_receiver = trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
    call addfld( fieldname_receiver, horiz_only, "A", unit, "condensation tendency" )
    if(history_aerosol)then
       call add_default( fieldname_receiver, 1, ' ' )
    end if

    fieldname_receiver = trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
    call addfld( fieldname_receiver, horiz_only, 'A', unit , "condensation tendency" )
    if(history_aerosol)then
       call add_default( fieldname_receiver, 1, ' ' )
    end if

    fieldname_receiver = trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
    call addfld( fieldname_receiver, horiz_only, 'A', unit, "condensation tendency" )
    if(history_aerosol)then
       call add_default( fieldname_receiver, 1, ' ' )
    end if

  end subroutine initializeCondensation

  !===============================================================================

  subroutine condtend(lchnk, q, cond_vap_gasprod, temperature, &
       pmid, pdel, dt, ncol, pblh, zm, qh20)

    ! Calculate the sulphate nucleation rate, and condensation rate of
    ! aerosols used for parameterising the transfer of externally mixed
    ! aitken mode particles into an internal mixture.
    ! Note the parameterisation for conversion of externally mixed particles
    ! used the h2so4 lifetime onto the particles, and not a given
    ! increase in particle radius. Will be improved in future versions of the model
    ! Added input for h2so4 and soa nucleation: soa_lv_gasprod, soa_sv_gasprod, pblh,zm,qh20 (cka)

    ! arguments
    integer,  intent(in) :: lchnk                      ! chunk identifier
    real(r8), intent(inout) :: q(pcols,pver,gas_pcnst) ! TMR [kg/kg] including moisture
    real(r8), intent(in) :: cond_vap_gasprod(pcols,pver,N_COND_VAP) ! TMR [kg/kg/sec]] production rate of H2SO4 (gas prod - aq phase uptake)
    real(r8), intent(in) :: temperature(pcols,pver)    ! Temperature (K)
    real(r8), intent(in) :: pmid(pcols,pver)           ! [Pa] pressure at mid point
    real(r8), intent(in) :: pdel(pcols,pver)           ! [Pa] difference in grid cell
    real(r8), intent(in) :: dt                         ! Time step
    integer,  intent(in) :: ncol                       ! number of columns
    ! Needed for soa nucleation treatment
    real(r8), intent(in) :: pblh(pcols)               ! pbl height (m)
    real(r8), intent(in) :: zm(pcols,pverp)           ! midlayer geopotential height above the surface (m) (pver+1)
    real(r8), intent(in) :: qh20(pcols,pver)          ! specific humidity (kg/kg)

    ! local
    character(len=fieldname_len+3) :: fieldname
    integer  :: i,k,nsiz
    integer  :: mode_index_donor            ![idx] index of mode donating mass
    integer  :: mode_index_receiver         ![idx] index of mode receiving mass
    integer  :: tracerIndex
    integer  :: l_donor
    integer  :: l_receiver
    integer  :: iDonor                                                                       ![idx] counter for externally mixed modes
    real(r8) :: condensationSink(0:nmodes, N_COND_VAP)                                       ![1/s] loss rate per mode (mixture)
    real(r8) :: condensationSinkFraction(pcols,pver,numberOfExternallyMixedModes,N_COND_VAP) ![frc]
    real(r8) :: sumCondensationSink(pcols,pver, N_COND_VAP)                                  ![1/s] sum of condensation sink
    real(r8) :: totalLoss(pcols,pver,gas_pcnst)                                              ![kg/kg] tracer lost
    real(r8) :: numberConcentration(0:nmodes)                                                ![#/m3] number concentration
    real(r8) :: numberConcentrationExtMix(pcols,pver,numberOfExternallyMixedModes)
    real(r8) :: coltend(pcols, gas_pcnst)
    real(r8) :: tracer_coltend(pcols)
    real(r8) :: intermediateConcentration(pcols,pver,N_COND_VAP)
    real(r8) :: rhoAir(pcols,pver)                           ![kg/m3] density of air

    ! Volume of added  material from condensate;  surface area of core particle;
    real(r8)            :: volume_shell, area_core,vol_monolayer
    real (r8)           :: frac_transfer               ! Fraction of hydrophobic material converted to an internally mixed mode
    logical             :: history_aerosol
    character(128)      :: long_name                   ! [-] needed for diagnostics

    ! needed for h2so4 and soa nucleation treatment
    integer             :: modeIndexReceiverCoag       ! Index of modes receiving coagulate
    integer             :: iCoagReceiver               ! counter for species receiving coagulate
    real(r8)            :: coagulationSink(pcols,pver) ! [1/s] coaglation loss for SO4_n and soa_n
    real(r8), parameter :: lvocfrac=0.5                ! Fraction of organic oxidation products with low enough

    !volatility to enter nucleation mode particles (1-24 nm)
    real(r8)            :: soa_lv_forNucleation(pcols,pver) ![kg/kg] soa gas available for nucleation
    real(r8)            :: gasLost(pcols,pver,N_COND_VAP)   ![kg/kg] budget terms on H2SO4 (gas)
    real(r8)            :: fracNucl(pcols,pver,N_COND_VAP)  ![frc] fraction of gas nucleated
    real(r8)            :: firstOrderLossRateNucl(pcols,pver,N_COND_VAP) ![1/s] first order loss rate due to nucleation
    real(r8)            :: nuclso4(pcols,pver) ![kg/kg/s] Nucleated so4 mass tendency from RM's parameterization
    real(r8)            :: nuclsoa(pcols,pver) ![kg/kg/s] Nucleated soa mass tendency from RM's parameterization
    integer             :: cond_vap_idx

    !Initialize h2so4 and soa nucl variables
    coagulationSink(:,:)=0.0_r8
    condensationSinkFraction(:,:,:,:) = 0.0_r8  !Sink to the coming "receiver" of any vapour
    numberConcentrationExtMix(:,:,:) = 0.0_r8

    do k=1,pver
       do i=1,ncol

          condensationSink(:,:) = 0.0_r8  !Sink to the coming "receiver" of any vapour

          !NB: The following is duplicated code, coordinate with oslo_aero_coag!
          !Initialize number concentration for this receiver

          !Air density
          rhoAir(i,k) = pmid(i,k)/rair/temperature(i,k)

          numberConcentration(:) = 0.0_r8

          !Go though all modes receiving condensation
          do mode_index_receiver = 0, nmodes

             !Go through all core species in that mode
             do tracerIndex = 1, getNumberOfBackgroundTracersInMode(mode_index_receiver)

                !Find the lifecycle-specie receiving the condensation
                l_receiver = getTracerIndex(mode_index_receiver, tracerIndex, .true.)

                !Add up the number concentration of the receiving mode [#/m3]
                numberConcentration(mode_index_receiver) = numberConcentration(mode_index_receiver) &  !previous value
                     + q(i,k,l_receiver)                   &  !kg/kg
                     / rhopart(physicsIndex(l_receiver))   &  !m3/kg ==> m3_{aer}/kg_{air}
                     * volumeToNumber(mode_index_receiver) &  !#/m3 ==> #/kg_{air}
                     * rhoAir(i,k)                                 !kg/m3 ==> #/m3_{air}
             end do !Lifecycle "core" species in this mode
          enddo


          !All modes are condensation receivers
          do cond_vap_idx=1,N_COND_VAP
             do mode_index_receiver = 0, nmodes

                !This is the loss rate a gas molecule will see due to aerosol surface area
                condensationSink(mode_index_receiver,cond_vap_idx)   = &             !==> [1/s]
                     normalizedCondensationSink(mode_index_receiver,cond_vap_idx)  & ![m3/#/s]
                   * numberConcentration(mode_index_receiver)                        ![#/m3]

             end do !Loop over receivers
          end do

          !Find concentration after condensation of all condenseable vapours
          do cond_vap_idx=1,N_COND_VAP

             !sum of cond. sink for this vapour [1/s]
             sumCondensationSink(i,k,cond_vap_idx) = sum(condensationSink(:,cond_vap_idx))

             !Solve the intermediate (end of timestep) concentration using
             !euler backward solution C_{old} + P *dt - L*C_{new}*dt = C_{new} ==>
             !Cnew -Cold = prod - loss ==>
             intermediateConcentration(i,k,cond_vap_idx) = &
                  ( q(i,k,cond_vap_map(cond_vap_idx)) + cond_vap_gasprod(i,k,cond_vap_idx)*dt ) &
                  / (1.0_r8 + sumCondensationSink(i,k,cond_vap_idx)*dt)
          end do

          !Save the fraction of condensation sink for the externally mixed modes
          !(Needed below to find volume shell)
          do cond_vap_idx=1,N_COND_VAP

             do iDonor = 1,numberOfExternallyMixedModes
                !Find the mode in question
                mode_index_donor = externallyMixedMode(iDonor)

                !Remember fraction of cond sink for this mode
                condensationSinkFraction(i,k,iDonor,cond_vap_idx) = &
                     condensationSink(mode_index_donor,cond_vap_idx) / sumCondensationSink(i,k,cond_vap_idx)

                !Remember number concentration in this mode
                numberConcentrationExtMix(i,k,iDonor) = numberConcentration(mode_index_donor)
             end do
          end do

          ! Assume only a fraction of ORG_LV left can contribute to nucleation
          ! fraction of soa_lv left that is assumend to have low enough volatility to nucleate.
          soa_lv_forNucleation(i,k) = lvocfrac*intermediateConcentration(i,k,COND_VAP_ORG_LV) 
          
          !Sum coagulation sink for nucleated so4 and soa particles over all receivers of coagulate. Needed for RM's nucleation code
          !OBS - looks like RM's coagulation sink is multiplied by 10^-12??
          modeIndexReceiverCoag = 0
          do iCoagReceiver = 1, numberOfCoagulationReceivers

             modeIndexReceiverCoag = receiverMode(iCoagReceiver)

             coagulationSink(i,k) =   &                                                  ![1/s]
                  coagulationSink(i,k) + &                                               ![1/] previous value
                  normalizedCoagulationSink(modeIndexReceiverCoag,MODE_IDX_SO4SOA_AIT) & ![m3/#/s]
                  * numberConcentration(modeIndexReceiverCoag)                           !numberConcentration (#/m3)
          end do    !coagulation sink

          !Sum coagulation sink for nucleated so4 and soa particles over all additional
          !receivers od coagulate (not directly affecting the life-cycle).
          do iCoagReceiver = 1, numberOfAddCoagReceivers

             modeIndexReceiverCoag = addReceiverMode(iCoagReceiver)

             coagulationSink(i,k) =   &                        ![1/s]
                  coagulationSink(i,k) + &                     ![1/] previous value
                  normCoagSinkAdd(iCoagReceiver) &             ![m3/#/s]
                  * numberConcentration(modeIndexReceiverCoag) !numberConcentration (#/m3)
          end do    !coagulation sink

       end do !index i
    end do !index k

    !Calculate nucleated masses of so4 and soa (nuclso4, nuclsoa)
    !following RM's parameterization (cka)
    call aeronucl(lchnk,ncol,temperature, pmid, qh20, &
         intermediateConcentration(:,:,COND_VAP_H2SO4), soa_lv_forNucleation, &
         coagulationSink, nuclso4, nuclsoa, zm, pblh)


    firstOrderLossRateNucl(:,:,:)=0.0_r8
    do k=1,pver
       do i=1,ncol

          !First order loss rate (1/s) for nucleation
          firstOrderLossRateNucl(i,k,COND_VAP_H2SO4) = nuclSo4(i,k)/intermediateConcentration(i,k,COND_VAP_H2SO4)

          !First order loss rate (1/s) for nucleation
          firstOrderLossRateNucl(i,k,COND_VAP_ORG_LV) = nuclSOA(i,k)/intermediateConcentration(i,k,COND_VAP_ORG_LV)

          do cond_vap_idx = 1,N_COND_VAP
             !Solve implicitly (again)
             !C_new - C_old =  PROD_{gas} - CS*C_new*dt - LR_{nucl}*C_new =>
             intermediateConcentration(i,k,cond_vap_idx) = &
                  ( q(i,k,cond_vap_map(cond_vap_idx)) + cond_vap_gasprod(i,k,cond_vap_idx)*dt ) &
                  / (1.0_r8 + sumCondensationSink(i,k,cond_vap_idx)*dt + firstOrderLossRateNucl(i,k,cond_vap_idx)*dt)

             !fraction nucleated
             fracNucl(i,k,cond_vap_idx) = firstOrderLossRateNucl(i,k,cond_vap_idx) &
                  /(firstOrderLossRateNucl(i,k,cond_vap_idx) + sumCondensationSink(i,k,cond_vap_idx))
             !From budget, we get: lost = prod -cnew + cold
             gasLost(i,k,cond_vap_idx) = cond_vap_gasprod(i,k,cond_vap_idx)*dt   & !Produced
                  + q(i,k,cond_vap_map(cond_vap_idx))            & !cold
                  - intermediateConcentration(i,k,cond_vap_idx)    !cnew

          end do !cond_vap_idx

          !Add nuceated mass to so4_na mode
          q(i,k,chemistryIndex(l_so4_na)) =  q(i,k,chemistryIndex(l_so4_na))       &
               + gasLost(i,k,COND_VAP_H2SO4)*fracNucl(i,k,COND_VAP_H2SO4)

          !H2SO4 condensate
          q(i,k,chemistryIndex(l_so4_a1)) = q(i,k,chemistryIndex(l_so4_a1))         &
               + gasLost(i,k,COND_VAP_H2SO4)*(1.0_r8-fracNucl(i,k,COND_VAP_H2SO4))

          !Add nucleated mass to soa_na mode
          q(i,k,chemistryIndex(l_soa_na)) =  q(i,k,chemistryIndex(l_soa_na))       &
               + gasLost(i,k,COND_VAP_ORG_LV)*fracNucl(i,k,COND_VAP_ORG_LV)

          !Organic condensate (from both soa_lv and soa_sv) goes to the soaCondensateReceiver tracer (cka)
          q(i,k,chemistryIndex(l_soa_a1)) = q(i,k,chemistryIndex(l_soa_a1))         &
               + gasLost(i,k,COND_VAP_ORG_SV)                             &           ! "semi volatile" can not nucleate
               + gasLost(i,k,COND_VAP_ORG_LV)*(1.0_r8-fracNucl(i,k,COND_VAP_ORG_LV))  ! part of low volatile which does not nucleate

          !condenseable vapours
          q(i,k,chemistryIndex(l_h2so4))  = intermediateConcentration(i,k,COND_VAP_H2SO4)
          q(i,k,chemistryIndex(l_soa_lv)) = intermediateConcentration(i,k,COND_VAP_ORG_LV)
          q(i,k,chemistryIndex(l_soa_sv)) = intermediateConcentration(i,k,COND_VAP_ORG_SV)


          !Condensation transfers mass from externally mixed to internally mixed modes
          do iDonor = 1,numberOfExternallyMixedModes

             !Find the mode in question
             mode_index_donor    = externallyMixedMode(iDonor)

             if(getNumberOfTracersInMode(mode_index_donor) .eq. 0)then
                cycle
             end if

             volume_shell = 0.0_r8
             do cond_vap_idx = 1, N_COND_VAP

                !Add up volume shell for this
                !condenseable vapour
                volume_shell = volume_shell                                               &
                     + condensationSinkFraction(i,k,iDonor,cond_vap_idx)                 & ![frc]
                     * gasLost(i,k,cond_vap_idx)*(1.0_r8-fracNucl(i,k,cond_vap_idx))     & ![kg/kg]
                     * invRhoPart(physicsIndex(cond_vap_map(cond_vap_idx)))              & !*[m3/kg] ==> [m3/kg_{air}
                     * rhoAir(i,k)                                                         !*[kg/m3] ==> m3/m3

             end do

             area_core=numberConcentrationExtMix(i,k,iDonor)*numberToSurface(mode_index_donor)   !#/m3 * m2/# ==> m2/m3
             vol_monolayer=area_core*dr_so4_monolayers_age

             ! Small fraction retained to avoid numerical irregularities
             frac_transfer=min((volume_shell/vol_monolayer),0.999_r8)

             !How many tracers exist in donor mode?
             !The "donor" is the externally mixed mode which will soon
             !become internally mixed. The externally mixed is donating mass
             !and the internally mixed is receiving...
             do tracerIndex = 1, getNumberOfTracersInMode(mode_index_donor)

                !Indexes here are in "chemistry space"
                l_donor    = getTracerIndex(mode_index_donor, tracerIndex,.true.)
                l_receiver = lifeCycleReceiver(l_donor)

                if( l_receiver .le. 0)then
                   stop !something wrong
                endif

                !Transfer from donor to receiver takes into account
                !fraction transferred
                totalLoss(i,k,l_donor) = frac_transfer*q(i,k,l_donor)
                q(i,k,l_donor) = q(i,k,l_donor) - totalLoss(i,k,l_donor)
                q(i,k,l_receiver) = q(i,k,l_receiver) + totalLoss(i,k,l_donor)
             end do !tracers in mode
          end do    !loop over receivers
       end do !physical index k
    end do    !physical index i

    !Output for diagnostics
    call phys_getopts(history_aerosol_out = history_aerosol)

    if(history_aerosol)then
       coltend(:ncol,:) = 0.0_r8
       do i=1,gas_pcnst
          !Check if species contributes to condensation
          if(lifeCycleReceiver(i) .gt. 0)then
             !Loss from the donor specie
             tracer_coltend(:ncol) = sum(totalLoss(:ncol, :,i)*pdel(:ncol,:),2)/gravit/dt
             coltend(:ncol,i) = coltend(:ncol,i) - tracer_coltend(:ncol) !negative (loss for donor)
             coltend(:ncol,lifeCycleReceiver(i)) = coltend(:ncol,lifeCycleReceiver(i)) + tracer_coltend(:ncol)
          endif
       end do

       ! Remove so4_n ---> directly into so4_na
       coltend(:ncol,chemistryIndex(l_so4_na)) = coltend(:ncol,chemistryIndex(l_so4_na)) + &
            sum(                                         &
            gasLost(:ncol,:,COND_VAP_H2SO4)           &
            *fracNucl(:ncol,:,COND_VAP_H2SO4)*pdel(:ncol,:) , 2 &
            )/gravit/dt

       !Take into account H2SO4 (gas) condensed in budget
       coltend(:ncol,chemistryIndex(l_so4_a1)) = coltend(:ncol,chemistryIndex(l_so4_a1)) + &
            sum(                                         &
            gasLost(:ncol,:,COND_VAP_H2SO4)           &
            *(1.0_r8 - fracNucl(:ncol,:,COND_VAP_H2SO4))*pdel(:ncol,:) , 2 &
            )/gravit/dt

       !Take into account soa_lv (gas) nucleated in budget
       coltend(:ncol,chemistryIndex(l_soa_na)) = coltend(:ncol,chemistryIndex(l_soa_na)) + &
            sum(                                         &
            gasLost(:ncol,:,COND_VAP_ORG_LV)              &
            *fracNucl(:ncol,:,COND_VAP_ORG_LV)*pdel(:ncol,:) , 2 &
            )/gravit/dt

       !Take into account soa gas condensed in the budget (both LV and SV)
       coltend(:ncol,chemistryIndex(l_soa_a1)) = coltend(:ncol,chemistryIndex(l_soa_a1)) + &
            sum(                                         &
            gasLost(:ncol,:,COND_VAP_ORG_LV)           &
            *(1.0_r8 - fracNucl(:ncol,:,COND_VAP_ORG_LV))*pdel(:ncol,:) , 2 &
            )/gravit/dt                        &
            +                                  &
            sum(                                         &
            gasLost(:ncol,:,COND_VAP_ORG_SV)*pdel(:ncol,:) , 2 &
            )/gravit/dt

       do i=1,gas_pcnst
          if(lifeCycleReceiver(i) .gt. 0 )then
             long_name= trim(solsym(i))//"condTend"
             call outfld(long_name, coltend(:ncol,i), pcols, lchnk)
             long_name= trim(solsym(lifeCycleReceiver(i)))//"condTend"
             call outfld(long_name, coltend(:ncol,lifeCycleReceiver(i)),pcols,lchnk)
          end if
       end do
       long_name=trim(solsym(chemistryIndex(l_so4_a1)))//"condTend"
       call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_a1)),pcols,lchnk)
       long_name=trim(solsym(chemistryIndex(l_soa_a1)))//"condTend"
       call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_a1)),pcols,lchnk)
       long_name=trim(solsym(chemistryIndex(l_so4_na)))//"condTend"
       call outfld(long_name, coltend(:ncol,chemistryIndex(l_so4_na)),pcols,lchnk)
       long_name=trim(solsym(chemistryIndex(l_soa_na)))//"condTend"
       call outfld(long_name, coltend(:ncol,chemistryIndex(l_soa_na)),pcols,lchnk)

    endif

  end subroutine condtend

  !===============================================================================

  subroutine aeronucl(lchnk, ncol, t, pmid, h2ommr, h2so4pc, oxidorg, coagnuc, nuclso4, nuclorg, zm, pblht)

    ! Subroutine to calculate nucleation (formation) rates of new particles
    ! At the moment, the final nucleation rate consists of
    !  (1) Binary sulphuric acid-water nucleation in whole atmosphere (Vehkamaki et al., 2002, JGR)
    !      JGR, vol 107, No D22, http://onlinelibrary.wiley.com/doi/10.1029/2002JD002184/abstract
    !  (2) Boundary-layer nucleation
    !      Paasonen et al (2010), ACP, vol 10, pp 11223: http://www.atmos-chem-phys.net/10/11223/2010/acp-10-11223-2010.html
    !  (3) First version published ACP (Risto Makkonen)
    !      ACP, vol 14, no 10, pp 5127 http://www.atmos-chem-phys.net/14/5127/2014/acp-14-5127-2014.html
    ! Modified Spring 2015, cka

    !-- Arguments
    integer,  intent(in)  :: lchnk                    ! chunk identifier
    integer,  intent(in)  :: ncol                     ! number of atmospheric column
    real(r8), intent(in)  :: pmid(pcols,pver)         ! layer pressure (Pa)
    real(r8), intent(in)  :: h2ommr(pcols,pver)       ! layer specific humidity
    real(r8), intent(in)  :: t(pcols,pver)            ! Temperature (K)
    real(r8), intent(in)  :: h2so4pc(pcols,pver)      ! Sulphuric acid concentration (kg kg-1)
    real(r8), intent(in)  :: oxidorg(pcols,pver)      ! Organic vapour concentration (kg kg-1)
    real(r8), intent(in)  :: coagnuc(pcols,pver)      ! Coagulation sink for nucleating particles [1/s]
    real(r8), intent(out) :: nuclorg(pcols,pver)      ! Nucleated mass (ORG)
    real(r8), intent(out) :: nuclso4(pcols,pver)      ! Nucleated mass (H2SO4)
    real(r8), intent(in)  :: zm(pcols,pver)           ! Height at layer midpoints (m)
    real(r8), intent(in)  :: pblht(pcols)             ! Planetary boundary layer height (m)

    !-- Local variables

    real(r8), parameter   :: pi=3.141592654_r8
    !cka+
    real(r8), parameter   :: gasconst_R=8.314472_r8    ! universal gas constant [J mol-1 K-1]
    real(r8), parameter   :: h2so4_dens=1841._r8       ! h2so4 density [kg m-3]
    real(r8), parameter   :: org_dens=2000._r8         ! density of organics [kg m-3], based on RM assumptions
    !cka -

    integer               :: i,k
    real(r8)              :: qs(pcols,pver)            ! Saturation specific humidity
    real(r8)              :: relhum(pcols,pver)        ! Relative humidity
    real(r8)              :: h2so4(pcols,pver)         ! Sulphuric acid concentration [#/cm3]
    real(r8)              :: nuclvolume(pcols,pver)    ! [m3/m3/s] Nucleated mass (SO4+ORG)
    real(r8)              :: rhoair(pcols,pver)        ! density of air [kg/m3] !cka
    real(r8)              :: pblht_lim(pcols)          ! Planetary boundary layer height (m) (500m<pblht_lim<7000m) (cka)

    real(r8)              :: nuclrate_bin(pcols,pver) ! Binary nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_bin(pcols,pver) ! Binary formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_bin(pcols,pver) ! Binary nucleation critical cluster size (m)
    real(r8)              :: nuclrate_pbl(pcols,pver) ! Boundary layer nucleation rate (# cm-3 s-1)
    real(r8)              :: formrate_pbl(pcols,pver) ! Boundary layer formation rate (12 nm) (# cm-3 s-1)
    real(r8)              :: nuclsize_pbl(pcols,pver) ! Boundary layer nucleation formation size (m)

    real(r8)              :: orgforgrowth(pcols,pver) ! Organic vapour mass available for growth

    real(r8)              :: d_form                   ! Particle size at calculated formation rate [m]
    real(r8)              :: gr(pcols,pver), grh2so4(pcols,pver), grorg(pcols,pver) !growth rates
    real(r8)              :: vmolh2so4, vmolorg       ! [m/s] molecular speed of condenseable gases
    real(r8)              :: frach2so4
    real(r8)              :: dummy

    integer               :: atm_nucleation           ! Nucleation parameterization for the whole atmosphere
    integer               :: pbl_nucleation           ! Nucleation parameterization for the boundary layer
    real(r8)              :: molmass_h2so4            ! molecular mass of h2so4 [g/mol]
    real(r8)              :: molmass_soa              ! molecular mass of soa [g/mol]

    ! Variables for binary nucleation parameterization
    real(r8)              :: zrhoa, zrh, zt, zt2, zt3, zlogrh, zlogrh2, zlogrh3, zlogrhoa, zlogrhoa2, zlogrhoa3, x, zxmole, zix
    real(r8)              :: zjnuc, zntot, zrc, zrxc

    !cka: OBS    call phys_getopts(pbl_nucleation_out=pbl_nucleation, atm_nucleation_out=atm_nucleation)
    !cka: testing by setting these flags:
    pbl_nucleation = 2
    atm_nucleation = 1

    nuclso4(:,:)=0._r8
    nuclorg(:,:)=0._r8
    nuclrate_bin(:,:)=0._r8
    nuclrate_pbl(:,:)=0._r8
    formrate_bin(:,:)=0._r8
    formrate_pbl(:,:)=0._r8
    !-- The highest level in planetary boundary layer
    do i=1,ncol
       pblht_lim(i)=MIN(MAX(pblht(i),500._r8),7000._r8)
    end do

    !-- Get molecular mass of h2so4 and soa_lv (cka)
    molmass_h2so4=adv_mass(id_H2SO4)
    molmass_soa=adv_mass(id_SOA_LV)

    !-- Formation diameters (m). Nucleated particles are inserted to SO4(n), same size used for soa  (cka)
    d_form=2._r8*originalNumberMedianRadius(MODE_IDX_SO4SOA_AIT)

    !-- Conversion of H2SO4 from kg/kg to #/cm3
    !-- and calculation of relative humidity (needed by binary nucleation parameterization)
    do k=1,pver
       do i=1,ncol
          rhoair(i,k)=pmid(i,k)/(t(i,k)*rair)
          !avogad*1.e-3_r8 to get molec/mol instead of molec/kmol
          h2so4(i,k)=(1.e-6_r8*h2so4pc(i,k)*avogad*1.e-3_r8*rhoair(i,k))/(molmass_h2so4*1.E-3_r8)
          orgforgrowth(i,k)=(1.e-6_r8*oxidorg(i,k)*avogad*1.e-3_r8*rhoair(i,k))/(molmass_soa*1.E-3_r8)
          orgforgrowth(i,k)=MAX(MIN(orgforgrowth(i,k),1.E10_r8),0._r8)

          call qsat_water(t(i,k),pmid(i,k),dummy,qs(i,k))

          relhum(i,k) = h2ommr(i,k)/qs(i,k)
          relhum(i,k) = max(relhum(i,k),0.0_r8)
          relhum(i,k) = min(relhum(i,k),1.0_r8)
       end do !ncol
    end do     !layers

    !-- Binary sulphuric acid-water nucleation rate
    if(atm_nucleation .EQ. 1) then
       do k=1,pver
          do i=1,ncol

             ! Calculate nucleation only for valid thermodynamic conditions:
             zrhoa = max(h2so4(i,k),1.E+4_r8)
             zrhoa = min(zrhoa,1.E11_r8)

             zrh   = max(relhum(i,k),1.E-4_r8)
             zrh   = min(zrh,1.0_r8)

             zt    = max(t(i,k),190.15_r8)
             zt    = min(zt,300.15_r8)

             zt2 = zt*zt
             zt3 = zt2*zt

             ! Equation (11) - molefraction of H2SO4 in the critical cluster

             zlogrh  = LOG(zrh)
             zlogrh2 = zlogrh*zlogrh
             zlogrh3 = zlogrh2*zlogrh

             zlogrhoa  = LOG(zrhoa)
             zlogrhoa2 = zlogrhoa*zlogrhoa
             zlogrhoa3 = zlogrhoa2*zlogrhoa

             x=0.7409967177282139_r8 - 0.002663785665140117_r8*zt   &
                  + 0.002010478847383187_r8*zlogrh    &
                  - 0.0001832894131464668_r8*zt*zlogrh    &
                  + 0.001574072538464286_r8*zlogrh2        &
                  - 0.00001790589121766952_r8*zt*zlogrh2    &
                  + 0.0001844027436573778_r8*zlogrh3     &
                  -  1.503452308794887e-6_r8*zt*zlogrh3    &
                  - 0.003499978417957668_r8*zlogrhoa   &
                  + 0.0000504021689382576_r8*zt*zlogrhoa

             zxmole=x

             zix = 1.0_r8/x

             ! Equation (12) - nucleation rate in 1/cm3s

             zjnuc=0.1430901615568665_r8 + 2.219563673425199_r8*zt -   &
                  0.02739106114964264_r8*zt2 +     &
                  0.00007228107239317088_r8*zt3 + 5.91822263375044_r8*zix +     &
                  0.1174886643003278_r8*zlogrh + 0.4625315047693772_r8*zt*zlogrh -     &
                  0.01180591129059253_r8*zt2*zlogrh +     &
                  0.0000404196487152575_r8*zt3*zlogrh +    &
                  (15.79628615047088_r8*zlogrh)*zix -     &
                  0.215553951893509_r8*zlogrh2 -    &
                  0.0810269192332194_r8*zt*zlogrh2 +     &
                  0.001435808434184642_r8*zt2*zlogrh2 -    &
                  4.775796947178588e-6_r8*zt3*zlogrh2 -     &
                  (2.912974063702185_r8*zlogrh2)*zix -   &
                  3.588557942822751_r8*zlogrh3 +     &
                  0.04950795302831703_r8*zt*zlogrh3 -     &
                  0.0002138195118737068_r8*zt2*zlogrh3 +    &
                  3.108005107949533e-7_r8*zt3*zlogrh3 -     &
                  (0.02933332747098296_r8*zlogrh3)*zix +     &
                  1.145983818561277_r8*zlogrhoa -    &
                  0.6007956227856778_r8*zt*zlogrhoa +    &
                  0.00864244733283759_r8*zt2*zlogrhoa -    &
                  0.00002289467254710888_r8*zt3*zlogrhoa -    &
                  (8.44984513869014_r8*zlogrhoa)*zix +    &
                  2.158548369286559_r8*zlogrh*zlogrhoa +   &
                  0.0808121412840917_r8*zt*zlogrh*zlogrhoa -    &
                  0.0004073815255395214_r8*zt2*zlogrh*zlogrhoa -   &
                  4.019572560156515e-7_r8*zt3*zlogrh*zlogrhoa +    &
                  (0.7213255852557236_r8*zlogrh*zlogrhoa)*zix +    &
                  1.62409850488771_r8*zlogrh2*zlogrhoa -    &
                  0.01601062035325362_r8*zt*zlogrh2*zlogrhoa +   &
                  0.00003771238979714162_r8*zt2*zlogrh2*zlogrhoa +    &
                  3.217942606371182e-8_r8*zt3*zlogrh2*zlogrhoa -    &
                  (0.01132550810022116_r8*zlogrh2*zlogrhoa)*zix +    &
                  9.71681713056504_r8*zlogrhoa2 -    &
                  0.1150478558347306_r8*zt*zlogrhoa2 +    &
                  0.0001570982486038294_r8*zt2*zlogrhoa2 +    &
                  4.009144680125015e-7_r8*zt3*zlogrhoa2 +    &
                  (0.7118597859976135_r8*zlogrhoa2)*zix -    &
                  1.056105824379897_r8*zlogrh*zlogrhoa2 +    &
                  0.00903377584628419_r8*zt*zlogrh*zlogrhoa2 -    &
                  0.00001984167387090606_r8*zt2*zlogrh*zlogrhoa2 +    &
                  2.460478196482179e-8_r8*zt3*zlogrh*zlogrhoa2 -    &
                  (0.05790872906645181_r8*zlogrh*zlogrhoa2)*zix -    &
                  0.1487119673397459_r8*zlogrhoa3 +    &
                  0.002835082097822667_r8*zt*zlogrhoa3 -    &
                  9.24618825471694e-6_r8*zt2*zlogrhoa3 +    &
                  5.004267665960894e-9_r8*zt3*zlogrhoa3 -    &
                  (0.01270805101481648_r8*zlogrhoa3)*zix

             zjnuc=EXP(zjnuc)      !   add. Eq. (12) [1/(cm^3s)]

             ! Equation (13) - total number of molecules in the critical cluster

             zntot=-0.002954125078716302_r8 - 0.0976834264241286_r8*zt +   &
                  0.001024847927067835_r8*zt2 - 2.186459697726116e-6_r8*zt3 -    &
                  0.1017165718716887_r8*zix - 0.002050640345231486_r8*zlogrh -   &
                  0.007585041382707174_r8*zt*zlogrh +    &
                  0.0001926539658089536_r8*zt2*zlogrh -   &
                  6.70429719683894e-7_r8*zt3*zlogrh -    &
                  (0.2557744774673163_r8*zlogrh)*zix +   &
                  0.003223076552477191_r8*zlogrh2 +   &
                  0.000852636632240633_r8*zt*zlogrh2 -    &
                  0.00001547571354871789_r8*zt2*zlogrh2 +   &
                  5.666608424980593e-8_r8*zt3*zlogrh2 +    &
                  (0.03384437400744206_r8*zlogrh2)*zix +   &
                  0.04743226764572505_r8*zlogrh3 -    &
                  0.0006251042204583412_r8*zt*zlogrh3 +   &
                  2.650663328519478e-6_r8*zt2*zlogrh3 -    &
                  3.674710848763778e-9_r8*zt3*zlogrh3 -   &
                  (0.0002672510825259393_r8*zlogrh3)*zix -    &
                  0.01252108546759328_r8*zlogrhoa +   &
                  0.005806550506277202_r8*zt*zlogrhoa -    &
                  0.0001016735312443444_r8*zt2*zlogrhoa +   &
                  2.881946187214505e-7_r8*zt3*zlogrhoa +    &
                  (0.0942243379396279_r8*zlogrhoa)*zix -   &
                  0.0385459592773097_r8*zlogrh*zlogrhoa -   &
                  0.0006723156277391984_r8*zt*zlogrh*zlogrhoa +   &
                  2.602884877659698e-6_r8*zt2*zlogrh*zlogrhoa +    &
                  1.194163699688297e-8_r8*zt3*zlogrh*zlogrhoa -   &
                  (0.00851515345806281_r8*zlogrh*zlogrhoa)*zix -    &
                  0.01837488495738111_r8*zlogrh2*zlogrhoa +   &
                  0.0001720723574407498_r8*zt*zlogrh2*zlogrhoa -   &
                  3.717657974086814e-7_r8*zt2*zlogrh2*zlogrhoa -    &
                  5.148746022615196e-10_r8*zt3*zlogrh2*zlogrhoa +    &
                  (0.0002686602132926594_r8*zlogrh2*zlogrhoa)*zix -   &
                  0.06199739728812199_r8*zlogrhoa2 +    &
                  0.000906958053583576_r8*zt*zlogrhoa2 -   &
                  9.11727926129757e-7_r8*zt2*zlogrhoa2 -    &
                  5.367963396508457e-9_r8*zt3*zlogrhoa2 -   &
                  (0.007742343393937707_r8*zlogrhoa2)*zix +    &
                  0.0121827103101659_r8*zlogrh*zlogrhoa2 -   &
                  0.0001066499571188091_r8*zt*zlogrh*zlogrhoa2 +    &
                  2.534598655067518e-7_r8*zt2*zlogrh*zlogrhoa2 -    &
                  3.635186504599571e-10_r8*zt3*zlogrh*zlogrhoa2 +    &
                  (0.0006100650851863252_r8*zlogrh*zlogrhoa2)*zix +   &
                  0.0003201836700403512_r8*zlogrhoa3 -    &
                  0.0000174761713262546_r8*zt*zlogrhoa3 +   &
                  6.065037668052182e-8_r8*zt2*zlogrhoa3 -    &
                  1.421771723004557e-11_r8*zt3*zlogrhoa3 +   &
                  (0.0001357509859501723_r8*zlogrhoa3)*zix

             zntot=EXP(zntot)  !  add. Eq. (13)

             ! Equation (14) - radius of the critical cluster in nm

             zrc=EXP(-1.6524245_r8+0.42316402_r8*x+0.33466487_r8*LOG(zntot))    ! [nm]

             !----1.2) Limiter

             IF(zjnuc<1.e-7_r8 .OR. zntot<4.0_r8) zjnuc=0.0_r8

             ! limitation to 1E+10 [1/cm3s]

             nuclrate_bin(i,k)=MAX(MIN(zjnuc,1.E10_r8),0._r8)
             nuclsize_bin(i,k)=MAX(MIN(zrc,1.E2_r8),0.01_r8)

          end do
       end do
    else   !No atmospheric nucleation
       nuclrate_bin(:,:)=0._r8
       nuclsize_bin(:,:)=1._r8
    end if

    !-- Boundary layer nucleation
    do k=1,pver
       do i=1,ncol

          !-- Nucleation rate #/cm3/s
          if(pblht_lim(i)>zm(i,k) .AND. pbl_nucleation>0) then

             if(pbl_nucleation .EQ. 1) then

                !-- Paasonen et al. (2010), eqn 10, Table 4
                nuclrate_pbl(i,k)=(1.7E-6_r8)*h2so4(i,k)

             else if(pbl_nucleation .EQ. 2) then

                !-- Paasonen et al. (2010)
                !values from Table 3 in Paasonen et al (2010), modified version of eqn 14
                nuclrate_pbl(i,k)=(6.1E-7_r8)*h2so4(i,k)+(0.39E-7_r8)*orgforgrowth(i,k)

             end if

             nuclrate_pbl(i,k)=MAX(MIN(nuclrate_pbl(i,k),1.E10_r8),0._r8)

          else !Not using PBL-nucleation
             nuclrate_pbl(i,k)=0._r8
          end if
          !Size [nm] of particles in PBL
          nuclsize_pbl(i,k)=2._r8

       end do !horizontal points
    end do     !levels

    !-- Calculate total nucleated mass
    do k=1,pver
       do i=1,ncol

          !   Molecular speed and growth rate: H2SO4. Eq. 21 in Kerminen and Kulmala 2002
          vmolh2so4=SQRT(8._r8*gasconst_R*t(i,k)/(pi*molmass_h2so4*1.E-3_r8))
          grh2so4(i,k)=(3.E-9_r8/h2so4_dens)*(vmolh2so4*molmass_h2so4*h2so4(i,k))
          grh2so4(i,k)=MAX(MIN(grh2so4(i,k),10000._r8),1.E-10_r8)

          !   Molecular speed and growth rate: ORG. Eq. 21 in Kerminen and Kulmala 2002
          vmolorg=SQRT(8._r8*gasconst_R*t(i,k)/(pi*molmass_soa*1.E-3_r8))
          grorg(i,k)=(3.E-9_r8/org_dens)*(vmolorg*molmass_soa*orgforgrowth(i,k))
          grorg(i,k)=MAX(MIN(grorg(i,k),10000._r8),1.E-10_r8)

          ! Combined growth rate (cka)
          gr(i,k)=grh2so4(i,k)+grorg(i,k)

          !-- Lehtinen 2007 parameterization for apparent formation rate
          !   diameters in nm, growth rate in nm h-1, coagulation in s-1

          call appformrate(nuclsize_bin(i,k), d_form*1.E9_r8, nuclrate_bin(i,k), formrate_bin(i,k), coagnuc(i,k), gr(i,k))
          call appformrate(nuclsize_pbl(i,k), d_form*1.E9_r8, nuclrate_pbl(i,k), formrate_pbl(i,k), coagnuc(i,k), gr(i,k))

          formrate_bin(i,k)=MAX(MIN(formrate_bin(i,k),1.E3_r8),0._r8)
          formrate_pbl(i,k)=MAX(MIN(formrate_pbl(i,k),1.E3_r8),0._r8)

          !   Number of mol nucleated per g air per second.
          nuclvolume(i,k) = (formrate_bin(i,k) + formrate_pbl(i,k)) & ![particles/cm3]
               *1.0e6_r8                                 & !==> [particles / m3 /]
               /volumeToNumber(MODE_IDX_SO4SOA_AIT)   & !==> [m3_{aer} / m3_{air} / sec]
               / rhoair(i,k)                            !==> m3_{aer} / kg_{air} /sec

          !Estimate how much is organic based on growth-rate
          if(gr(i,k)>1.E-10_r8) then
             frach2so4=grh2so4(i,k)/gr(i,k)
          else
             frach2so4=1._r8
          end if

          ! Nucleated so4 and soa mass mixing ratio per second [kg kg-1 s-1]
          ! used density of particle phase, not of condensing gas
          nuclso4(i,k)=rhopart(l_so4_na)*nuclvolume(i,k)*frach2so4
          nuclorg(i,k)=rhopart(l_soa_na)*nuclvolume(i,k)*(1.0_r8-frach2so4)

       end do
    end do

    !-- Diagnostic output
    call outfld('NUCLRATE', nuclrate_bin+nuclrate_pbl, pcols   ,lchnk)
    call outfld('FORMRATE', formrate_bin+formrate_pbl, pcols   ,lchnk)
    call outfld('COAGNUCL', coagnuc, pcols   ,lchnk)
    call outfld('GRH2SO4', grh2so4, pcols   ,lchnk)
    call outfld('GRSOA', grorg, pcols   ,lchnk)
    call outfld('GR', gr, pcols   ,lchnk)

    return
  end subroutine aeronucl

  !===============================================================================

  subroutine appformrate(d1, dx, j1, jx, CoagS_dx, gr)
    !-- appformrate calculates the formation rate jx of dx sized particles from the nucleation rate j1 (d1 sized particles)
    !-- Formation rate is parameterized according to Lehtinen et al. (2007), JAS 38:988-994
    !-- Parameterization takes into account the loss of particles due to coagulation
    !-- Growth by self-coagulation is not accounted for
    !-- Typically, 1% of 1 nm nuclei make it to 12 nm
    !-- Written by Risto Makkonen
    ! First estimate: 99% of particles are lost during growth from 1 nm to 12 nm

    !-- Arguments

    real(r8), intent(in)  :: d1                  ! Size of nucleation-sized particles (nm)
    real(r8), intent(in)  :: dx                  ! Size of calculated apparent formation rate (nm)
    real(r8), intent(in)  :: j1                  ! Nucleation rate of d1 sized particles (# cm-3 s-1)
    real(r8), intent(out) :: jx                  ! Formation rate of dx sized particles (# cm-3 s-1)
    real(r8), intent(in)  :: CoagS_dx            ! Coagulation term for nucleating particles (s-1)
    real(r8), intent(in)  :: gr                  ! Particle growth rate (nm h-1)

    !-- Local variables

    real(r8)              :: m
    real(r8)              :: gamma
    real(r8)              :: CoagS_d1            ! Coagulation term for nucleating particles, calculated from CoagS_dx

    ! In Hyytiala, typically 80% of the nuclei are scavenged onto larger background particles while they grow from 1 to 3 nm

    !-- (Eq. 6) Exponent m, depends on background distribution
    ! m=log(CoagS_dx/CoagS_d1)/log(dx/d1)
    ! Or, if we dont want to calculate CoagS_d1, lets assume a typical value for m (-1.5 -- -1.9) and calculate CoagS_d1 from Eq.5
    m=-1.6_r8
    CoagS_d1=CoagS_dx*(d1/dx)**m
    CoagS_d1=MAX(MIN(CoagS_d1,1.E2_r8),1.E-10_r8)

    gamma=(1._r8/(m+1._r8))*((dx/d1)**(m+1._r8)-1._r8)
    gamma=MAX(MIN(gamma,1.E2_r8),1.E-10_r8)

    !-- (Eq. 7) CoagS_d1 is multiplied with 3600 to get units h-1
    jx=j1*exp(-gamma*d1*CoagS_d1*3600._r8/gr)

  end subroutine appformrate

end module oslo_aero_condtend
