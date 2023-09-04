module oslo_aero_conc

  ! Calculate concentrations of aerosol modes based on lifecycle species

  use shr_kind_mod ,         only: r8 => shr_kind_r8
  use ppgrid       ,         only: pcols, pver
  use physconst    ,         only: density_water =>rhoh2o, molecularWeightWater=>mwh2o, pi
  use constituents ,         only: pcnst, cnst_name
  !
  use oslo_aero_logn_tables, only: intlog1to3_sub, intlog4_sub, intlog5to10_sub, initlogn
  use oslo_aero_utils,       only: calculateNumberConcentration
  use oslo_aero_coag,        only: normalizedCoagulationSink
  use oslo_aero_condtend,    only: normalizedCondensationSink, COND_VAP_H2SO4, COND_VAP_ORG_SV
  use oslo_aero_const,       only: smallNumber, volumeToNumber,smallNumber
  use oslo_aero_params
  use oslo_aero_share

  implicit none
  private

  public :: oslo_aero_conc_calc
  public :: calculateBulkProperties
  public :: partitionMass

  private :: getAerosolMask
  private :: calculateHygroscopicity
  private :: addModeHygroscopicity
  private :: doLognormalInterpolation
  private :: modalapp2d

  ! Size of molecule-layer which defines when particles are coated
  real(r8), parameter :: coatingLimit = 2.e-9_r8  ![m]

  ! The fraction of soluble material required in a components before it!will add to any coating
  real(r8), parameter :: solubleMassFractionCoatingLimit=0.50_r8
  real(r8), parameter :: aThird       = 1.0_r8/3.0_r8
  real(r8), parameter :: ln10         = log(10.0_r8)

  logical :: init_logn_tables = .false.

contains

  !********************************************************************************************
  subroutine oslo_aero_conc_calc(ncol, mmr, rho_air, CProcessModes, &
       f_c, f_bc, f_aq, f_so4_cond, f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam, &
       numberConcentration, volumeConcentration, hygroscopicity, lnsigma, hasAerosol, volumeCore, volumeCoat)

    !----------------------------------------------
    ! Calculate concentrations of aerosol modes based on lifecycle species
    !----------------------------------------------

    ! arguments
    integer,  intent(in)  :: ncol                                     ! Number of columns used in chunk
    real(r8), intent(in)  :: mmr(pcols,pver,pcnst)                    ! [kg/kg] mass mixing ratio of tracers
    real(r8), intent(in)  :: rho_air(pcols,pver)                      ! [kg/m3] air density
    logical,  intent(out) :: hasAerosol(pcols, pver, nmodes)          ! [t/f] do we have this type of aerosol here?
    real(r8), intent(out) :: f_acm(pcols,pver, nbmodes)               ! [frc] carbon fraction in mode
    real(r8), intent(out) :: f_bcm(pcols,pver, nbmodes)               ! [frc] fraction of c being bc
    real(r8), intent(out) :: f_aqm(pcols, pver, nbmodes)              ! [frc] fraction of sulfate being aquous
    real(r8), intent(out) :: f_so4_condm(pcols, pver, nbmodes)        ! [frc] fraction of non-aquous SO4 being condensate
    real(r8), intent(out) :: f_soam(pcols, pver, nbmodes)             ! Needed in "get component fraction"
    real(r8), intent(out) :: numberConcentration(pcols,pver,0:nmodes) ! [#/m3] number concentraiton
    real(r8), intent(out) :: volumeConcentration(pcols,pver,nmodes)   ! [m3/m3] volume concentration
    real(r8), intent(out) :: hygroscopicity(pcols,pver,nmodes)        ! [mol_{aer}/mol_{water}] hygroscopicity
    real(r8), intent(out) :: lnsigma(pcols,pver,nmodes)               ! [-] log(base e) sigma
    real(r8), intent(out) :: CProcessModes(pcols,pver)
    real(r8), intent(out) :: cam(pcols,pver,nbmodes)
    real(r8), intent(out) :: f_c(pcols, pver)
    real(r8), intent(out) :: f_aq(pcols,pver)
    real(r8), intent(out) :: f_bc(pcols,pver)
    real(r8), intent(out) :: f_so4_cond(pcols,pver)
    real(r8), intent(out) :: f_soa(pcols,pver)
    real(r8), intent(out) :: volumeCore(pcols,pver,nmodes)
    real(r8), intent(out) :: volumeCoat(pcols,pver,nmodes)

    ! local variables
    real(r8) :: f_aitbc(pcols,pver) ! [-] bc fraction in the coated bc-oc mode
    real(r8) :: f_nbc(pcols,pver)   ! [-] mass fraction of bc in uncoated bc/oc mode
    real(r8) :: f_soana(pcols,pver) ! [-]

    !Get mass, number concentration and the total add-ons (previous convaer)
    call calculateBulkProperties(ncol, mmr, rho_air, numberConcentration, CProcessModes, &
         f_c, f_bc, f_aq, f_so4_cond, f_soa, f_aitbc, f_nbc, f_soana)

    ! Find the points where we have aerosol (number concentration)
    call getAerosolMask(ncol, numberConcentration, hasAerosol)

    ! Find out how much is added per size-mode (modalapp)
    call partitionMass( ncol, numberConcentration, CProcessModes, &
         f_c, f_bc, f_aq, f_so4_cond, f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam)

    ! Calculate they hygroscopicity
    call calculateHygroscopicity(  ncol, mmr, numberConcentration, rho_air, Cam, &
         f_acm, f_bcm, f_aqm, hasAerosol, hygroscopicity, &
         volumeConcentration, volumeCore, volumeCoat)

    ! Do the interpolation to new modes
    call doLognormalInterpolation(ncol, numberConcentration, hasAerosol, cam, &
         volumeConcentration, f_c, f_acm, f_bcm, f_aqm, f_aitbc, lnSigma)

  end subroutine oslo_aero_conc_calc

  !******************************************************************
  subroutine calculateBulkProperties( ncol, qm, rho_air, numberConcentration, CProcessModes, &
       f_c, f_bc, f_aq, f_so4_cond, f_soa, f_aitbc, f_nbc, f_soana)

    !----------------------------------------------
    ! Create bulk properties (dependent on tracers, not size modes)
    !----------------------------------------------

    ! arguments
    integer , intent(in)  :: ncol                 ! [nbr] number of columns used
    real(r8), intent(in)  :: qm(pcols,pver,pcnst) ! [kg/kg] mmr for transported tracers
    real(r8), intent(in)  :: rho_air(pcols,pver)  ! [kg/m3] air density
    real(r8), intent(out) :: numberConcentration(pcols,pver,0:nmodes) ! [#/m3] aerosol number concentration
    real(r8), intent(out) :: f_c(pcols,pver)        ![-] mass fraction of process mode being c
    real(r8), intent(out) :: f_bc(pcols,pver)       ![-] mass fraction of c being bc
    real(r8), intent(out) :: f_aq(pcols,pver)       ![-] mass fraction of s being aq phase
    real(r8), intent(out) :: f_so4_cond(pcols,pver) ![-] mass fraction of non-aq s being condensate
    real(r8), intent(out) :: f_soa(pcols,pver)      ![-] mass fraction of OM being SOA
    real(r8), intent(out) :: f_aitbc(pcols,pver)    ![-] mass fraction of bc in bc/oc mixed, coated mode
    real(r8), intent(out) :: f_nbc(pcols,pver)      ![-] mass fraction of bc in bc/oc mixed, un-coated mode
    real(r8), intent(out) :: f_soana(pcols,pver)    ![-] mass fraction of soa in background in int mix ait mode (1)

    !Local variables
    real(r8) :: totalProcessModes(pcols,pver)    ! [kg/kg] Int. mixed (cond./coag./aq.) SO4+BC+OC concentration
    real(r8) :: CProcessModes(pcols,pver)        ! [kg/m3] Int. mixed (cond./coag./aq.) SO4+BC+OC concentration
    integer  :: k  !counter for layers

    ! Total number concentration per mode
    call calculateNumberConcentration(ncol, qm, rho_air, numberConcentration)

    do k=1,pver

       !Total coagulated bc and oc and SO4 (condensate, wet phase and coagulated) (kg/kg)
       !internally mixed with background modes
       totalProcessModes(:ncol,k)  = qm(:ncol,k,l_bc_ac) + qm(:ncol,k,l_om_ac) &
            +  qm(:ncol,k,l_so4_a1) + qm(:ncol,k,l_so4_a2) + qm(:ncol,k,l_so4_ac) + qm(:ncol,k,l_soa_a1)

       CProcessModes(:ncol,k) = rho_air(:ncol,k)*totalProcessModes(:ncol,k)  !==> kg/m3

       !fraction of process-mode being carbonaceous
       f_c(:ncol,k)   = min((qm(:ncol,k,l_bc_ac)+qm(:ncol,k,l_om_ac)+qm(:ncol,k,l_soa_a1) )&
            /(totalProcessModes(:ncol,k)+smallNumber), 1.0_r8)

       !fraction of "c" being bc (total is oc and bc)
       f_bc(:ncol,k)  = min(qm(:ncol,k,l_bc_ac)/(qm(:ncol,k,l_bc_ac)+qm(:ncol,k,l_om_ac)+qm(:ncol,k,l_soa_a1)+smallNumber), 1.0_r8)

       !fraction of non-aqeous phase sulphate being condensate
       f_so4_cond(:ncol,k) = min(qm(:ncol,k,l_so4_a1)/(qm(:ncol,k,l_so4_a1)+qm(:ncol,k,l_so4_ac)+smallNumber), 1.0_r8)

       !fraction of sulphate being aquous phase (total is condensate + aqeous phase + coagulate)
       f_aq(:ncol,k)  = min(qm(:ncol,k,l_so4_a2) &
            /(qm(:ncol,k,l_so4_a1)+qm(:ncol,k,l_so4_a2)+qm(:ncol,k,l_so4_ac)+smallNumber),1.0_r8)

       !fraction of bc in the sulfate-coated bc/oc mode (total background is bc and oc)
       f_aitbc(:ncol,k) = min(qm(:ncol,k,l_bc_ai) / (qm(:ncol,k,l_bc_ai) + qm(:ncol,k,l_om_ai) + smallNumber), 1.0_r8)

       !fraction of bc in the un-coated bc/oc (total is bc and oc)
       f_nbc(:ncol,k) = min(qm(:ncol,k,l_bc_ni) / (qm(:ncol,k,l_bc_ni) + qm(:ncol,k,l_om_ni) + smallNumber),1.0_r8)

       !fraction of OM process-mode which is SOA
       f_soa(:ncol,k) = min(qm(:ncol,k,l_soa_a1) / (qm(:ncol,k,l_om_ac) + qm(:ncol,k,l_soa_a1) + smallNumber), 1.0_r8)

       !fraction of "background" int-mix (mode 1) which is SOA
       f_soana(:ncol,k) = min(qm(:ncol,k,l_soa_na) / (qm(:ncol,k,l_soa_na) + qm(:ncol,k,l_so4_na) + smallNumber), 1.0_r8 )

    end do !k

    return
  end subroutine calculateBulkProperties

  !********************************************************************************
  subroutine partitionMass(  ncol, Nnatk, CProcessModes, &
       f_c, f_bc, f_aq, f_so4_cond, f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam)

    integer , intent(in)  :: ncol                        ! [nbr] number of columns used
    real(r8), intent(in)  :: Nnatk(pcols,pver,0:nmodes)  ! [#/m3] number concentration
    real(r8), intent(in)  :: CProcessModes(pcols,pver)   ! [kg/m3] total added mass
    real(r8), intent(in)  :: f_c(pcols,pver)             ! [frc] fraction of added mass being c
    real(r8), intent(in)  :: f_bc(pcols,pver)            ! [frc] fraction of c being bc
    real(r8), intent(in)  :: f_aq(pcols,pver)            ! [frc] fraction of SO4 being aq
    real(r8), intent(in)  :: f_so4_cond(pcols,pver)      ! [frc] fraction of SO4 coag+cond being cond
    real(r8), intent(in)  :: f_soa(pcols,pver)           ! [frc] fraction of OM being SOA
    real(r8), intent(out) :: cam(pcols, pver, nbmodes)   ! [kg/m3] added mass distributed to modes
    real(r8), intent(out) :: f_acm(pcols,pver,nbmodes)   ! [frc] as f_c per mode
    real(r8), intent(out) :: f_bcm(pcols,pver,nbmodes)   ! [frc] as f_bc per mode
    real(r8), intent(out) :: f_aqm(pcols,pver,nbmodes)   ! [frc] as f_aq per mode
    real(r8), intent(out) :: f_so4_condm(pcols,pver,nbmodes) ! [frc] fraction of non aq sulfate being coagulate
    real(r8), intent(out) :: f_soam(pcols,pver,nbmodes)  ! [frc] fraction of OC being SOA

    call modalapp2d(ncol, Nnatk(1,1,1), CProcessModes, &
         f_c, f_bc, f_aq, f_so4_cond, f_soa, cam, f_acm, f_bcm, f_aqm, f_so4_condm, f_soam)

  end subroutine partitionMass

  !*************************************************************
  subroutine getAerosolMask(ncol,numberConcentration, hasAerosol)

    ! Find out where we have aerosols

    integer, intent(in)   :: ncol         !number of columns used
    real(r8), intent(in)  :: numberConcentration(pcols, pver, 0:nmodes)
    logical, intent(out)  :: hasAerosol(pcols, pver, nmodes)
    integer               :: k !counter for levels
    integer               :: m !counter for modes

    do m=1,nmodes
       do k=1,pver
          where(numberConcentration(:ncol,k,m) .gt. smallNumber)
             hasAerosol(:ncol,k,m)= .true.
          elsewhere
             hasAerosol(:ncol,k,m) = .false.
          end where
       end do  !levels
    end do     !modes

  end subroutine getAerosolMask

  !*************************************************************
  subroutine calculateHygroscopicity(ncol, mmr, numberConcentration, rho_air, Cam, &
       f_acm, f_bcm, f_aqm, hasAerosol, hygroscopicity, volumeConcentration, volumeCore,volumeCoat)

    ! A parameterization of aerosol activation 2. Multiple aerosol types, JGR, vol 105, noD5, pp 6837
    ! http://onlinelibrary.wiley.com/doi/10.1029/1999JD901161/abstract
    ! Abdul-Razzak and S. Ghan:

    ! arguments
    integer  , intent(in)  :: ncol
    real(r8) , intent(in)  :: mmr(pcols,pver,pcnst)                    ! [kg/kg] mass mixing ratios
    real(r8) , intent(in)  :: numberConcentration(pcols,pver,0:nmodes) ! [#/m3] number concentrations
    real(r8) , intent(in)  :: rho_air(pcols,pver)                      ! [kg/m3] air density
    real(r8) , intent(in)  :: Cam(pcols, pver, nbmodes)                ! [kg/m3] total added mass during microphysics
    real(r8) , intent(in)  :: f_acm(pcols,pver,nbmodes)                ! [-] fraction of added mass which is carbon
    real(r8) , intent(in)  :: f_aqm(pcols,pver,nbmodes)                ! [-] fraction of sulfate which is aq. phase
    real(r8) , intent(in)  :: f_bcm(pcols,pver,nbmodes)                ! [-] fraction of C which is bc
    logical  , intent(in)  :: hasAerosol(pcols,pver,nmodes)            ! [t/f] do we have aerosols
    real(r8) , intent(out) :: hygroscopicity(pcols,pver,nmodes)
    real(r8) , intent(out) :: volumeConcentration(pcols,pver,nmodes)
    real(r8) , intent(out) :: volumeCore(pcols,pver,nmodes)            ![m3]
    real(r8) , intent(out) :: volumeCoat(pcols,pver,nmodes)            ![m3]

    ! local variables
    integer  :: kcomp !counter for modes
    integer  :: l     !counter for components
    integer  :: k     !counter for levels
    integer  :: tracerIndex
    integer  :: i
    real(r8) :: hygroscopicityAvg(pcols,pver)
    real(r8) :: hygroscopicityCoat(pcols,pver)
    real(r8) :: massConcentrationTracerInMode(pcols,pver)
    real(r8) :: averageRadiusCore(pcols,pver)  ![m]
    real(r8) :: averageRadiusTotal(pcols,pver) ![m]

    ! initialize
    hygroscopicity(:,:,:) = 0.0_r8
    volumeConcentration(:,:,:)=0.0_r8

    do kcomp=1,nmodes

       !Don't do anything if no tracers in mode
       if(getNumberOfBackgroundTracersInMode(kcomp) .lt. 1)then
          volumeCore(:,:,kcomp)=smallNumber
          volumeCoat(:,:,kcomp)=smallNumber
          volumeConcentration(:,:,kcomp)=smallNumber
          hygroscopicity(:,:,kcomp) = smallNumber
          cycle
       end if

       hygroscopicityAvg(:,:) = 0.0_r8
       hygroscopicityCoat(:,:) = 0.0_r8
       volumeCore(:,:,kcomp) = 0.0_r8
       volumeCoat(:,:,kcomp) = 0.0_r8

       !Loop over tracers in mode
       do l=1,getNumberOfBackgroundTracersInMode(kcomp)

          tracerIndex = getTracerIndex(kcomp,l,.false.) !get index in physcis space

          do k=1,pver
             massConcentrationTracerInMode(:ncol,k) = mmr(:ncol,k,tracerIndex)*rho_air(:ncol,k)
          end do

          ! hasAerosol is true if any concentration in this point
          call addModeHygroscopicity(   ncol, hasAerosol(:,:,kcomp), &
               massConcentrationTracerInMode, volumeCore(:,:,kcomp), volumeCoat(:,:,kcomp), &
               hygroscopicityAvg, hygroscopicityCoat, tracerIndex)

       end do !background tracers in mode (l)

       !The background modes can have tracer mass added to them
       if (kcomp .le. nbmodes)then

          ! added aquous sulfate
          if(isTracerInMode(kcomp,l_so4_a2))then

             do k=1,pver
                massConcentrationTracerInMode(:ncol,k) = Cam(:ncol,k,kcomp)*(1.0_r8 - f_acm(:ncol,k,kcomp))*f_aqm(:ncol,k,kcomp)
             end do

             ! hasAerosol is true if any concentration in this point
             call addModeHygroscopicity( ncol, hasAerosol(:,:,kcomp), &
                  massConcentrationTracerInMode, volumeCore(:,:,kcomp), volumeCoat(:,:,kcomp), &
                  hygroscopicityAvg, hygroscopicityCoat, l_so4_a2)

          endif

          ! added condensate/coagulate
          ! All modes which have coagulate have also condensate, so it is
          ! ok to check for condensate and add the combined mass..
          if (isTracerInMode(kcomp,l_so4_a1))then
             do k=1,pver
                massConcentrationTracerInMode(:ncol,k) = Cam(:ncol,k,kcomp)*(1.0_r8 - f_acm(:ncol,k,kcomp))*(1.0_r8 - f_aqm(:ncol,k,kcomp))
             end do

             call addModeHygroscopicity(ncol, hasAerosol(:,:,kcomp),  &
                  massConcentrationTracerInMode, volumeCore(:,:,kcomp), volumeCoat(:,:,kcomp), &
                  hygroscopicityAvg, hygroscopicityCoat, l_so4_a1)

          endif

          ! Added bc
          if (isTracerInMode(kcomp,l_bc_ac))then
             do k=1,pver
                massConcentrationTracerInMode(:ncol,k) = Cam(:ncol,k,kcomp)*f_acm(:ncol,k,kcomp)*f_bcm(:ncol,k,kcomp)
             end do

             call addModeHygroscopicity( ncol, hasAerosol(:,:,kcomp), &
                  massConcentrationTracerInMode, volumeCore(:,:,kcomp), volumeCoat(:,:,kcomp), &
                  hygroscopicityAvg, hygroscopicityCoat, l_bc_ac )
          endif

          ! Added oc (both POM and SOA), then both have the same
          ! properties, so add combined mass here.
          ! All modes which have condensate also has coagulate, so OK to check
          ! for condensate and distribute the sum..
          if (isTracerInMode(kcomp,l_soa_a1))then
             do k=1,pver
                massConcentrationTracerInMode(:ncol,k) = Cam(:ncol,k,kcomp)*f_acm(:ncol,k,kcomp)*(1.0_r8 -f_bcm(:ncol,k,kcomp))
             end do

             call addModeHygroscopicity( ncol                            &
                  , hasAerosol(:,:,kcomp)         &  !true if any concentration in this point
                  , massConcentrationTracerInMode &
                  , volumeCore(:,:,kcomp)         &
                  , volumeCoat(:,:,kcomp)         &
                  , hygroscopicityAvg             &
                  , hygroscopicityCoat            &
                  , l_om_ac                       &
                  )
          endif
       end if

       !Note: NCAR definitions of molecular weights are kg/kmol. This is used
       !inside "addModeHygroscopicity" and here as in molecularWeightWater. SI units are kg/mol, but
       !the error cancels out since eqn 4 has Mw_water/Mw_tracer

       do k=1,pver

          !Finally, when the sums are calculated, Apply finally eqn 4 here!!

          where (hasAerosol(:ncol,k,kcomp))
             where(VolumeCoat(:ncol,k,kcomp) .gt. 1.e-30_r8)
                !If there is enough soluble material, a coating will be formed: In that case, the
                !volume of the aerosol in question is only the volume of the coating!
                hygroscopicityCoat(:ncol,k) = molecularWeightWater*hygroscopicityCoat(:ncol,k) &
                     /( density_water * volumeCoat(:ncol,k,kcomp)) !Note use of volume Coating here
             elsewhere
                hygroscopicityCoat(:ncol,k) = 1.e-30_r8
             endwhere
             !mode total volume:
             volumeConcentration(:ncol,k,kcomp) = volumeCore(:ncol,k,kcomp) + volumeCoat(:ncol,k,kcomp)

             !hygroscopicity of mixture (Note use of total volume to get average hygroscopicity)
             hygroscopicityAvg(:ncol,k) = molecularWeightWater*hygroscopicityAvg(:ncol,k) &
                  /(density_water * volumeConcentration(:ncol,k,kcomp))


             !Average size of insoluble core (average radius)
             averageRadiusCore(:ncol,k) = 0.5_r8*( (volumeCore(:ncol,k,kcomp)) &
                  / numberConcentration(:ncol,k,kcomp) * (6.0_r8/pi))**athird

             !Average size of total aerosol (average radius)
             averageRadiusTotal(:ncol,k) = 0.5_r8*((volumeConcentration(:ncol,k,kcomp)) &
                  / numberConcentration(:ncol,k,kcomp)*(6.0_r8/pi))**athird

             !do i=1,ncol
             !   if(numberConcentration(i,k,kcomp) .gt. 1.e6 .and. kcomp.eq.6 )then
             !      print*, "hygro_check",kcomp,numberConcentration(i,k,kcomp), averageRadiusTotal(i,k)*1.e6, averageRadiusCore(i,k)*1.e6 &
             !               , hygroscopicityCoat(i,k), hygroscopicityAvg(i,k), (averageRadiusTotal(i,k)-averageRadiusCore(i,k))*1.e9
             !   endif
             !end do

             ! use one or the other hygroscopicity based on coating
             where ( averageRadiusTotal(:ncol,k) - averageRadiusCore(:ncol,k)  .gt. coatingLimit )
                hygroscopicity(:ncol,k,kcomp) = hygroscopicityCoat(:ncol,k)
             elsewhere
                hygroscopicity(:ncol,k,kcomp) = hygroscopicityAvg(:ncol,k)
             endwhere

          elsewhere ! No aerosol

             hygroscopicity(:ncol,k,kcomp) = 1.e-10_r8

          end where

       end do !levels

    end do !kcomp /modes

  end subroutine calculateHygroscopicity

  !**************************************************************************************
  subroutine addModeHygroscopicity (ncol, hasAerosol, massConcentrationTracerInMode, &
       volumeCore, volumeCoat, hygroscopicityAvg, hygroscopicityCoat, tracerIndex)

    ! arguments
    integer  , intent(in)    :: ncol
    logical  , intent(in)    :: hasAerosol(pcols,pver)                    ![bool] true if we have any aerosol here
    real(r8) , intent(in)    :: massConcentrationTracerInMode(pcols,pver) ![kg/m3] mass concentration in
    integer  , intent(in)    :: tracerIndex                               !in physics space
    real(r8) , intent(inout) :: volumeCore(pcols, pver)                   !O [m3/m3] volume of insoluble core
    real(r8) , intent(inout) :: volumeCoat(pcols, pver)                   !O [m3/m3] volume of total aerosol
    real(r8) , intent(inout) :: hygroscopicityAvg(pcols, pver)            !O [-] average hygroscopicity
    real(r8) , intent(inout) :: hygroscopicityCoat(pcols, pver)           !O [-] average hygroscopicity

    ! local variables
    real(r8)                :: massFractionInCoating
    integer                 :: k                                        !counter for levels

    ! Only tracers more soluble than 20% can add to the coating volume
    if(solubleMassFraction(tracerIndex) .gt. solubleMassFractionCoatingLimit)then
       massFractionInCoating = 1.0_r8 !all volume goes to coating
    else
       massFractionInCoating = 0.0_r8 !zero volume goes to coating
    endif

    do k=1,pver

       where(hasAerosol(:ncol,k) .eqv. .true.)

          volumeCore(:ncol,k) = volumeCore(:ncol,k) &
               + massConcentrationTracerInMode(:ncol,k)/rhopart(tracerIndex)*(1.0_r8 - massFractionInCoating)

          volumeCoat(:ncol,k) = volumeCoat(:ncol,k) &
               + massConcentrationTracerInMode(:ncol,k)/rhopart(tracerIndex)*massFractionInCoating

          !sum up numerator in eqn 4 in Abdul-Razzak et al (average
          !hygrocopicity) Note that molecular weight is that of the
          !AEROSOL TYPE This is because of some conflict with mozart
          !which needs molecular weight of OC tracers to be 12 when
          !reading emissions So molecular weight is duplicated, and
          !the molecular weight of the TYPE is used here!

          hygroscopicityAvg(:ncol,k) = hygroscopicityAvg(:ncol,k) +  &
               massConcentrationTracerInMode(:ncol,k)*numberOfIons(tracerIndex)*osmoticCoefficient(tracerIndex) &
               *solubleMassFraction(tracerIndex)/aerosol_type_molecular_weight(aerosolType(tracerIndex))

          !Contribution to hygroscopicity of coating (only if goes to coating)
          !sum up numerator in eqn 4 in Abdul-Razzak et al (average hygrocopicity)
          !Note that molecular weight is that of the AEROSOL TYPE
          !This is because of some conflict with mozart which needs
          !molecular weight of OC tracers to be 12 when reading
          !emissions So molecular weight is duplicated, and the
          !molecular weight of the TYPE is used here!

         hygroscopicityCoat(:ncol,k) = hygroscopicityCoat(:ncol,k) +  &
               massConcentrationTracerInMode(:ncol,k)*numberOfIons(tracerIndex)*osmoticCoefficient(tracerIndex) &
               *solubleMassFraction(tracerIndex)/aerosol_type_molecular_weight(aerosolType(tracerIndex))           &
               *massFractionInCoating !Only add to this if mass goes to coating

       elsewhere

          hygroscopicityAvg(:ncol,k) = 1.0e-10_r8
          hygroscopicityCoat(:ncol,k)= 1.0e-10_r8

       end where

    end do

  end subroutine addModeHygroscopicity

  !****************************************************************
  subroutine doLognormalInterpolation(ncol, numberConcentration, hasAerosol, &
       cam, volumeConcentration, f_c, f_acm, f_bcm, f_aqm, f_aitbc, lnSigma)

    ! arguments
    integer  , intent(in)    :: ncol
    real(r8) , intent(in)    :: volumeConcentration(pcols,pver,nmodes)
    logical  , intent(in)    :: hasAerosol(pcols,pver,nmodes)
    real(r8) , intent(in)    :: cam(pcols,pver,nbmodes)                  ![kg/m3] total added mass per mode
    real(r8) , intent(in)    :: f_c(pcols,pver)                          ![frc] fraction of carbon in total add-on
    real(r8) , intent(in)    :: f_acm(pcols,pver,nbmodes)                ![frc] fraction of carbon per mode (in add-on)
    real(r8) , intent(in)    :: f_bcm(pcols,pver,nbmodes)                ![frc] fraction of bc in carbon per mode
    real(r8) , intent(in)    :: f_aqm(pcols,pver,nbmodes)                ![frc] fraction of aq in sulfate added
    real(r8) , intent(in)    :: f_aitbc(pcols,pver)                      ![frc] fraction of bc in coated bc/oc mode
    real(r8) , intent(inout) :: numberConcentration(pcols,pver,0:nmodes) ![#/m3] number concentration
    real(r8) , intent(out)   :: lnsigma(pcols,pver,nmodes)               ![-] log (base e) of std. dev

    ! local variables
    integer  :: kcomp
    integer  :: i,k
    real(r8) :: nconccm3(pcols,pver)
    real(r8) :: camUg(pcols,pver)
    real(r8) :: log10sig(pcols,pver)    ! [-] logarithm (base 10) of look up tables
    real(r8) :: f_ocm(pcols,pver,4)     ! [-] fraction of added mass which is either SOA condensate or OC coagulate
    real(r8) :: cxs(pcols,pver,nbmodes) ![ug/m3] NOTE NON-SI UNITS non-allocated mass
    real(r8) :: radius_tmp(pcols,pver)  ![m] radius in look up tables

    ! Initialize logn tables for interpolation
    if (.not. init_logn_tables) then
       call initlogn()
       init_logn_tables = .true.
    end if

    ! total mass not allocated to any mode
    ! this is non-zero if the look-up table can not cope with all the add-on mass
    ! cxstot(:,:) = 0.0_r8

    ! calculate fraction of added mass which is either SOA condensate or OC coagulate,
    ! which in AeroTab are both treated as condensate for kcomp=1-4
    do kcomp=1,4
       do k=1,pver
          do i=1,ncol
             f_ocm(i,k,kcomp) = f_acm(i,k,kcomp)*(1.0_r8-f_bcm(i,k,kcomp))
          enddo
       enddo
    enddo

    ! Go through all "background" size-modes (kcomp=1-10)
    do kcomp=1,nbmodes

       camUg(:,:) = cam(:,:,kcomp)*1.e9_r8
       nConccm3(:,:) = 1e-6_r8*numberConcentration(:,:,kcomp)

       ! Calculate growth from knowing added process specific internally mixed mass to each background mode
       ! (level sent but not needed, and kcomp not needed for intlog4_sub)

       if ( kcomp .ge. MODE_IDX_SO4SOA_AIT .and. kcomp .le. MODE_IDX_BC_AIT) then  ! kcomp=1,2

          do k=1,pver
             call intlog1to3_sub(   &
                  ncol,             & !I number of points
                  kcomp,            & !I [idx] mode index
                  camUg(:,k),       & !I [ug/m3] mass concentration
                  nConccm3(:,k),    & !I [#/cm3] number concentration
                  f_ocm(:,k,kcomp), & !I [frc] mass fraction which is SOA cond. or OC coag.
                  cxs(:,k,kcomp),   & !O [ug/m3] mass which did not fit the table
                  log10sig(:,k),    & !O [-]sigma, is later thrown away begause of volume balance
                  radius_tmp(:,k)   & !O [m] Number median radius
                  )
          end do  !loop on levels

       else if (kcomp .eq. MODE_IDX_OMBC_INTMIX_COAT_AIT) then ! kcomp=4

          do k=1,pver
             call intlog4_sub(      &
                  ncol,             & !I [nbr] number of points
                  kcomp,            & !I [idx] mode index
                  camUg(:,k),       & !I [ug/m3] mass concentration
                  nConccm3(:,k),    & !I [#/cm3] number concentration
                  f_ocm(:,k,kcomp), & !I [frc] mass fraction which is SOA cond. or OC coag.
                  f_aqm(:,k,kcomp), & !I [frc] fraction of sulfate which is aquous
                  cxs(:,k,kcomp),   & !O [ug/m3] mass which did not fit the table
                  log10sig(:,k),    & !O [-]sigma, is later thrown away begause of volume balance
                  radius_tmp(:,k)   & !O [m] Number median radius
                  )
          end do

       else if (kcomp .ge. MODE_IDX_SO4_AC .and. kcomp .le. MODE_IDX_SS_A3)then    ! kcomp=5-10

          do k=1,pver
             call intlog5to10_sub(  &
                  ncol,             & !I [nbr] number of points used
                  kcomp,            & !I [mode index]
                  camUg(:,k),       & !I [ug/m3] mass concentration
                  nConccm3(:,k),    & !I [#/cm3] number concentration
                  f_acm(:,k,kcomp), & !I [frc] fraction of aerosol which is carbon
                  f_bcm(:,k,kcomp), & !I [frc] fraction of carbon which is bc
                  f_aqm(:,k,kcomp), & !I [frc] fraction of sulfate which is aquous
                  cxs(:,k,kcomp),   & !O [ug/m3] mass which did not fit the table (not given to any mode)
                  log10sig(:,k),    & !O logarithm (base 10) sigma, is later thrown away begause of volume balance
                  radius_tmp(:,k)   & !O [m] Number median radius
                  )
          end do ! k

       endif

       !initialize
       lnsigma(:,:,kcomp) = log(2.0_r8)

       !The whole point of the interpolation routines is to get the new sigma ==> so trust the sigma

       !This means that in order to conserve the volume (which is known), we have to throw away
       !the number concentration. Should create a diagnostic or a warning if number concenration is very different
       !from the original number concentration since in principal, the number concentration is
       !also conserved!
       do k=1,pver
          !Don't change number concentration unless "hasAerosol" is true
          where(hasAerosol(:ncol,k,kcomp))

             lnsigma(:ncol,k,kcomp) = ln10*log10sig(:ncol,k)

             numberConcentration(:ncol,k,kcomp) = volumeConcentration(:ncol,k,kcomp)*6.0_r8/pi      &
                  /(2.0_r8*radius_tmp(:ncol,k))**3  &
                  *DEXP(-4.5_r8*lnsigma(:ncol,k,kcomp)*lnsigma(:ncol,k,kcomp))

             !==> Now we have a set of n, vol, sigma which is consistent and gives back whatever the
             !lookup tables told us! If the look up tables were conserving volume we didn't have to do
             !the step just above!!

             !Sum up all mass which was not added to any mode (mass exceeding the max limit in the look-up tables)
             !cxstot(:ncol,k) = cxstot(:ncol,k) + cxs(:ncol,k,kcomp)*1.e-9_r8 ! ug/m3 ==> kg/m3

          end where
       end do

    end do !kcomp

    !The modes which do not have any added aerosol:
    do kcomp=nbmodes+1,nmodes
       do k=1,pver
          lnsigma(:ncol,k,kcomp) = log(originalSigma(kcomp))
       end do
    end do

    !AK (fxm): "unactivated" code below...
    !Excessive internally mixed process mass added to the background modes (exceeding the max limit in the look-up tables)
    !is instead added to / lumped with the externally mixed non-background modes (kcomp=11,12,14)
    !numberConcentration(:,:,MODE_IDX_SO4_NUC) = numberConcentration(:,:,MODE_IDX_SO4_NUC) &
    !                                    + (volumeToNumber(MODE_IDX_SO4_NUC) &          !excess sulfate mass is moved to this mode
    !                                     *RESHAPE(cxstot,(/pcols,pver/)) &
    !                                     *(1.0_r8-f_c(:,:))/rhopart(l_so4_n))

    !numberConcentration(:,:,MODE_IDX_BC_NUC) = numberConcentration(:,:,MODE_IDX_BC_NUC) &
    !                                    + (volumeToNumber(MODE_IDX_BC_NUC)  &          !excess carbon mass is moved to this mode
    !                                    * RESHAPE(cxstot,(/pcols,pver/)) &
    !                                    * f_c(:,:)/rhopart(l_bc_n))

    !SKIP LUMPING OF OC-MODE TO MODE MODE_IDX_LUMPED ORGANICS SINCE THIS WILL MESS UP THE HASAEROSOL-MASK!
    !   modedefs(i)%Nnatk(MODE_IDX_LUMPED_ORGANICS) = efact_omn &   !excess OM mass is moved to this mode (originally kcomp=13)
    !            * (modedefs(i)%Nnatk(MODE_IDX_LUMPED_ORGANICS) + cxstot(i)*modedefs(i)%f_c*(1.0_r8-modedefs(i)%f_bc))

  end subroutine doLognormalInterpolation

  !********************************************************************************************
  subroutine modalapp2d(ncol,Nnatkbg,Ca,f_c,f_bc,f_aq,f_so4_cond,f_soa,Cam,fcm,fbcm,faqm,fso4condm,fsoam)

    !     Calculation of the apportionment of internally mixed SO4, BC and OC
    !     mass between the various background mineral and sea-salt modes.
    !     Now also Aitken-modes are subject to condensation of H2SO4, and both n and
    !     Aitken modes may coagulate onto the mineral/sea-salt background aerosol.
    ! SOA
    !     May 2013: The SO4(Ait) mode now takes into account condensed SOA in addition
    !     to H2SO4, but as long as SOA is not allowed to condense on more than one
    !     mode, no changes are necessary here. NB: to allow SOA to condense also on
    !     the BC(Ait) and/or other modes, change this code accordingly! Without any
    !     changes, Cam(pcols,1) = condensed SO4 onto the SO4(ait) mode still.
    ! SOA
    !     Alf Grini, february 2014 : Added info about units,
    !     used values calculated at initialization.
    !     changed in-out variables to components of derived data types (modedefs)
    !     defined in microphysics_oslo.F90, and corrected for mass balance error
    !     for SO4 due to lumping of coagulate and condensate.

    ! Arguments
    integer , intent(in)  :: ncol                          ! number of columns used
    real(r8), intent(in)  :: Nnatkbg(pcols,pver,nbmodes)   ! aerosol background mode number concentration     #/m3
    real(r8), intent(in)  :: Ca(pcols,pver)                ! internally mixed mass, tot=SO4+OC+BC
    real(r8), intent(in)  :: f_c(pcols,pver)               ! mass fraction (OC+BC)/tot
    real(r8), intent(in)  :: f_bc(pcols,pver)              ! mass fraction BC/(OC+BC)
    real(r8), intent(in)  :: f_aq(pcols,pver)              ! mass fraction SO4(aq)/SO4
    real(r8), intent(in)  :: f_soa(pcols,pver)             ! mass fraction SOA/(POM+SOA)
    real(r8), intent(in)  :: f_so4_cond(pcols,pver)        ! mass fraction SO4_COND/(COND+COAG)
    real(r8), intent(out) :: Cam(pcols,pver,nbmodes)       ! modal internal mass, tot=SO4+BC+OC
    real(r8), intent(out) :: fcm(pcols,pver,nbmodes)       ! modal mass fraction (OC+BC)/tot
    real(r8), intent(out) :: fbcm(pcols,pver,nbmodes)      ! modal mass fraction BC/(OC+BC)
    real(r8), intent(out) :: faqm(pcols,pver,nbmodes)      ! modal mass fraction SO4(aq)/SO4
    real(r8), intent(out) :: fso4condm(pcols,pver,nbmodes) ! modal mass fraction (SO4(cond)/SO4(cond+coag))
    real(r8), intent(out) :: fsoam(pcols,pver,nbmodes)     ! modal mass fraction SOA / (POM+SOA)

    !
    ! Local variables
    real(r8) condensationSinkSO4(pcols,pver,nbmodes) ![1/s] loss rate of cond. vap on any mode
    real(r8) condensationSinkOA(pcols,pver,nbmodes)  ![1/s] loss rate of cond. vap on any mode
    real(r8) coagulationSink(pcols,pver,nbmodes)     ![1/s] loss rate of BC through coagulation on any mode
    real(r8) aquousPhaseSink(pcols,pver,nbmodes)     ![-] fraction of particles available for aq. phase in any mode

    real(r8) sumCondensationSinkSO4(pcols,pver)      ![1/s] sum condensation sink to all modes
    real(r8) sumCondensationSinkOA(pcols,pver)       ![1/s] sum condensation sink to all modes
    real(r8) sumCoagulationSink(pcols,pver)          ![1/s] sum coagulation sink to all modes
    real(r8) sumAquousPhaseSink(pcols,pver)          ![1/s] sum aquous phase sink to all modes

    real(r8) fcondkSO4(pcols,pver,nbmodes)
    real(r8) fcondkOA(pcols,pver,nbmodes)
    real(r8) fcoagk(pcols,pver,nbmodes)
    real(r8) faqk(pcols,pver,nbmodes)

    real(r8) cabck(pcols,pver,nbmodes)               ![kg/m3] bc distributed to each mode
    real(r8) caock(pcols,pver,nbmodes)               ![kg/m3] pom coagulate distributed to each mode
    real(r8) csoacondsk(pcols,pver,nbmodes)
    real(r8) caqsk(pcols,pver,nbmodes)               ![kg/m3] aq phase sulfate distributed to each mode
    real(r8) cso4condsk(pcols,pver,nbmodes)          ![kg/m3] non-aq sulfate condensate distributed to each mode
    real(r8) cso4coagsk(pcols,pver,nbmodes)          ![kg/m3] non-aq sulfate coagulate distributed to each mode
    real(r8) cso4condcoagsk(pcols,pver,nbmodes)      ![kg/m3] non-aq sulfate condensate distributed to each mode
    real(r8) coccondcoagsk(pcols,pver,nbmodes)       ![kg/m3] non-aq sulfate coagulate distributed to each mode

    integer :: i !counter for modes
    integer :: k !counter for levels

    !Find the sink on any mode (0 is omitted here, WHY??, it does receive matter in oslo_aero_coag/condtend!!))
    !Should either remove it from there or add something to it here!
    do i=1,nbmodes
       do k=1,pver
          condensationSinkSO4(:ncol,k,i) = normalizedCondensationSink(i,COND_VAP_H2SO4)*Nnatkbg(:ncol,k,i)
          condensationSinkOA(:ncol,k,i) = normalizedCondensationSink(i,COND_VAP_ORG_SV)*Nnatkbg(:ncol,k,i)
          coagulationSink(:ncol,k,i)  = normalizedCoagulationSink(i,MODE_IDX_BC_NUC)*Nnatkbg(:ncol,k,i) !use a typical coagulator (BC_NUC)
          aquousPhaseSink(:ncol,k,i)  = numberFractionAvailableAqChem(i)*Nnatkbg(:ncol,k,i)             !aq phase sink to this mode
       end do
    enddo

    !Sum the sinks
    sumCondensationSinkSO4(:,:) = 0.0_r8
    sumCondensationSinkOA(:,:) = 0.0_r8
    sumCoagulationSink(:,:) = 0.0_r8
    sumAquousPhaseSink(:,:) = 0.0_r8
    do i=1,nbmodes
       do k=1,pver
          sumCondensationSinkSO4(:ncol,k) = sumCondensationSinkSO4(:ncol,k) + condensationSinkSO4(:ncol,k,i)
          sumCondensationSinkOA(:ncol,k) = sumCondensationSinkOA(:ncol,k) + condensationSinkOA(:ncol,k,i)
          sumCoagulationSink(:ncol,k) = sumCoagulationSink(:ncol,k) + coagulationSink(:ncol,k,i)
          sumAquousPhaseSink(:ncol,k) = sumAquousPhaseSink(:ncol,k) + aquousPhaseSink(:ncol,k,i)
       end do
    end do

    ! And finally the contribution from each mode relative to the totals are calculated,
    ! assuming that the apportionment of mass for the first iteration (in time) is representative
    ! for the whole apportionment process (which is ok for small and moderate masses added):
    do i=1,nbmodes
       do k=1,pver
          !Get the fraction of contribution per process per mode
          fcondkSO4(:ncol,k,i)=condensationSinkSO4(:ncol,k,i)/(sumCondensationSinkSO4(:ncol,k)+1.e-100_r8)  !fraction of condensation sink in this mode
          fcondkOA(:ncol,k,i)=condensationSinkOA(:ncol,k,i)/(sumCondensationSinkOA(:ncol,k)+1.e-100_r8)  !fraction of condensation sink in this mode
          fcoagk(:ncol,k,i)=coagulationSink(:ncol,k,i)/(sumCoagulationSink(:ncol,k)+1.e-100_r8)    !fraction of coagulation sink in this mode
          faqk(:ncol,k,i)=aquousPhaseSink(:ncol,k,i)/(sumAquousPhaseSink(:ncol,k)+1.e-100_r8)      !fraction of aquous phase sink in this mode

          !BC coagulate to this mode [kg/m3]
          cabck(:ncol,k,i)=fcoagk(:ncol,k,i)*f_c(:ncol,k)*f_bc(:ncol,k)*Ca(:ncol,k)

          !OC coagulate to this mode [kg/m3]
          caock(:ncol,k,i)=fcoagk(:ncol,k,i)*f_c(:ncol,k)*(1.0_r8-f_bc(:ncol,k))*(1.0_r8-f_soa(:ncol,k))*Ca(:ncol,k)

          !SOA condensate to this mode [kg/m3]
          csoacondsk(:ncol,k,i) = fcondkOA(:ncol,k,i)*f_c(:ncol,k)*(1.0_r8-f_bc(:ncol,k))*f_soa(:ncol,k)*Ca(:ncol,k)

          !Aquous phase SO4 to this mode [kg/m3]
          caqsk(:ncol,k,i)=faqk(:ncol,k,i)*f_aq(:ncol,k)*(1.0_r8-f_c(:ncol,k))*Ca(:ncol,k)

          !so4 condensate
          cso4condsk(:ncol,k,i)=fcondkSO4(:ncol,k,i)*(1.0_r8-f_aq(:ncol,k))*f_so4_cond(:ncol,k)*(1.0_r8-f_c(:ncol,k))*Ca(:ncol,k)

          !soa coagulate
          cso4coagsk(:ncol,k,i) = fcoagk(:ncol,k,i)*(1.0_r8-f_aq(:ncol,k))*(1.0_r8-f_so4_cond(:ncol,k))*(1.0_r8-f_c(:ncol,k))*Ca(:ncol,k) ![kg/m3] so4 coagulate
       end do
    enddo

    !The tables take as input the combined coagulate and condensate (both POM and SOA)
    !The activation needs them separately for mass balance!
    cso4condcoagsk(:ncol,:,:) = cso4condsk(:ncol,:,:) + cso4coagsk(:ncol,:,:)
    coccondcoagsk(:ncol,:,:) =  caock(:ncol,:,:) + csoacondsk(:ncol,:,:)

    do i=1,nbmodes
       do k=1,pver
          Cam(:ncol,k,i)=  cabck(:ncol,k,i)           &                               !BC
               + coccondcoagsk(:ncol,k,i)   &                               !OM
               + caqsk(:ncol,k,i) + cso4condcoagsk(:ncol,k,i)  + smallNumber!SO4 ==>   !total process mode mass to mode i

          fcm(:ncol,k,i)=(cabck(:ncol,k,i)+coccondcoagsk(:ncol,k,i))/(Cam(:ncol,k,i)+smallNumber)       !fraction of mass being carbon (oc or bc)
          fbcm(:ncol,k,i)=cabck(:ncol,k,i)/(cabck(:ncol,k,i)+coccondcoagsk(:ncol,k,i)+smallNumber)      !fraction of carbon mass being bc
          faqm(:ncol,k,i)=caqsk(:ncol,k,i)/(caqsk(:ncol,k,i)+cso4condcoagsk(:ncol,k,i)+smallNumber)     !fraction of sulfate being aq phase

          !Not  needed for tables, but for mass balances in activation
          fso4condm(:ncol,k,i) = cso4condsk(:ncol,k,i)/(cso4condcoagsk(:ncol,k,i) + smallNumber) !fraction of cond+coag which is coag
          fsoam(:ncol,k,i) = csoacondsk(:ncol,k,i)/(coccondcoagsk(:ncol,k,i) + smallNumber) !fraction of OC which is SOA
       end do
    enddo

  end subroutine modalapp2d

end module oslo_aero_conc
