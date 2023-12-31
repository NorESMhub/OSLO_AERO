module oslo_aero_dust_sediment

  !---------------------------------------------------------------------------------
  ! Routines to compute tendencies from sedimentation of dust
  ! Author: Phil Rasch
  !---------------------------------------------------------------------------------

  use shr_kind_mod,   only: r8=>shr_kind_r8
  use ppgrid,         only: pcols, pver, pverp
  use physconst,      only: gravit, rair
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none
  private

  ! public routines
  public :: oslo_aero_dust_sediment_vel
  public :: oslo_aero_dust_sediment_tend

  ! private routines
  private :: getflx
  private :: cfint2
  private :: cfdotmc_pro

  real (r8), parameter :: vland  = 2.8_r8            ! dust fall velocity over land  (cm/s)
  real (r8), parameter :: vocean = 1.5_r8            ! dust fall velocity over ocean (cm/s)
  real (r8), parameter :: mxsedfac = 0.99_r8       ! maximum sedimentation flux factor

!===============================================================================
contains
!===============================================================================

  subroutine oslo_aero_dust_sediment_vel(ncol, icefrac, landfrac, ocnfrac, pmid, pdel, t, dustmr, pvdust)

    ! Compute gravitational sedimentation velocities for dust
    ! note that pvel is at the interfaces (loss from cell is based on pvel(k+1))

    ! Arguments
    integer,  intent(in)  :: ncol                 ! number of colums to process
    real(r8), intent(in)  :: icefrac (pcols)      ! sea ice fraction (fraction)
    real(r8), intent(in)  :: landfrac(pcols)      ! land fraction (fraction)
    real(r8), intent(in)  :: ocnfrac (pcols)      ! ocean fraction (fraction)
    real(r8), intent(in)  :: pmid(pcols,pver)     ! pressure of midpoint levels (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)     ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t(pcols,pver)        ! temperature (K)
    real(r8), intent(in)  :: dustmr(pcols,pver)   ! dust (kg/kg)
    real(r8), intent(out) :: pvdust (pcols,pverp) ! vertical velocity of dust (Pa/s)

    ! Local variables
    real (r8) :: rho(pcols,pver)                    ! air density in kg/m3
    real (r8) :: vfall(pcols)                       ! settling velocity of dust particles (m/s)
    integer   :: i,k
    real (r8) :: lbound, ac, bc, cc

    ! dust fall velocity
    do k = 1,pver
       do i = 1,ncol
          ! merge the dust fall velocities for land and ocean (cm/s) SHOULD ALSO ACCOUNT FOR ICEFRAC
          vfall(i) = vland*landfrac(i) + vocean*(1._r8-landfrac(i))

          ! fall velocity (assume positive downward)
          pvdust(i,k+1) = vfall(i)     
       end do
    end do
  end subroutine oslo_aero_dust_sediment_vel

  !===============================================================================
  subroutine oslo_aero_dust_sediment_tend ( ncol,   dtime,  pint, pmid, pdel, t, &
       dustmr, pvdust, dusttend, sfdust, dusttend_to_ll_out )

    !----------------------------------------------------------------------
    !  Apply Particle Gravitational Sedimentation 
    ! -> note that pvel is at the interfaces (loss from cell is based on pvel(k+1))
    !----------------------------------------------------------------------

    ! Arguments
    integer,  intent(in)  :: ncol                      ! number of colums to process
    real(r8), intent(in)  :: dtime                     ! time step
    real(r8), intent(in)  :: pint(pcols,pverp)         ! interfaces pressure (Pa)
    real(r8), intent(in)  :: pmid(pcols,pver)          ! midpoint pressures (Pa)
    real(r8), intent(in)  :: pdel(pcols,pver)          ! pressure diff across layer (Pa)
    real(r8), intent(in)  :: t(pcols,pver)             ! temperature (K)
    real(r8), intent(in)  :: dustmr(pcols,pver)        ! dust (kg/kg)
    real(r8), intent(in)  :: pvdust (pcols,pverp)      ! vertical velocity of dust drops  (Pa/s)
    real(r8), intent(out) :: dusttend(pcols,pver)      ! dust tend
    real(r8), intent(out) :: sfdust(pcols)             ! surface flux of dust (rain, kg/m/s)
    real(r8),intent(out),optional  :: dusttend_to_ll_out(pcols) ! fluxes at the interfaces, dust (positive = down

    ! Local variables
    integer  :: i,k
    real(r8) :: fxdust(pcols,pverp)                     ! fluxes at the interfaces, dust (positive = down)
    !----------------------------------------------------------------------

    ! initialize variables
    fxdust  (:ncol,:) = 0._r8 ! flux at interfaces (dust)
    dusttend(:ncol,:) = 0._r8 ! tend (dust)
    sfdust(:ncol)     = 0._r8 ! sedimentation flux out bot of column (dust)

    ! fluxes at interior points
    call getflx(ncol, pint, dustmr, pvdust, dtime, fxdust)

    ! calculate fluxes at boundaries
    do i = 1,ncol
       fxdust(i,1) = 0
       ! surface flux by upstream scheme
       fxdust(i,pverp) = dustmr(i,pver) * pvdust(i,pverp) * dtime
    end do

    ! filter out any negative fluxes from the getflx routine
    do k = 2,pver
       fxdust(:ncol,k) = max(0._r8, fxdust(:ncol,k))
    end do

    ! Limit the flux out of the bottom of each cell to the water content in each phase.
    ! Apply mxsedfac to prevent generating very small negative cloud water/ice
    ! NOTE, REMOVED CLOUD FACTOR FROM AVAILABLE WATER. ALL CLOUD WATER IS IN CLOUDS.
    ! ***Should we include the flux in the top, to allow for thin surface layers?
    ! ***Requires simple treatment of cloud overlap, already included below.
    do k = 1,pver
       do i = 1,ncol
          fxdust(i,k+1) = min( fxdust(i,k+1), mxsedfac * dustmr(i,k) * pdel(i,k) )
       end do
    end do

    ! Now calculate the tendencies 
    do k = 1,pver
       do i = 1,ncol
          ! net flux into cloud changes cloud dust/ice (all flux is out of cloud)
          dusttend(i,k)  = (fxdust(i,k) - fxdust(i,k+1)) / (dtime * pdel(i,k))
       end do
    end do

    ! convert flux out the bottom to mass units Pa -> kg/m2/s
    sfdust(:ncol) = fxdust(:ncol,pverp) / (dtime*gravit)

    ! fluxes at the interface
    if(present(dusttend_to_ll_out))then
       dusttend_to_ll_out(1:ncol) = fxdust(:ncol,pver)/(dtime*pdel(:ncol,pver))
    end if

  end subroutine oslo_aero_dust_sediment_tend

  !===============================================================================
  subroutine getflx(ncol, xw, phi, vel, deltat, flux)

    !.....xw1.......xw2.......xw3.......xw4.......xw5.......xw6
    !....psiw1.....psiw2.....psiw3.....psiw4.....psiw5.....psiw6
    !....velw1.....velw2.....velw3.....velw4.....velw5.....velw6
    !.........phi1......phi2.......phi3.....phi4.......phi5.......

    ! arguments
    integer , intent(in)  :: ncol  ! number of colums to process
    real(r8), intent(out) :: flux(pcols,pverp)
    real(r8), intent(in)  :: xw(pcols,pverp)
    real(r8), intent(in)  :: vel(pcols,pverp)
    real(r8), intent(in)  :: deltat

    ! local variables
    integer   :: i
    integer   :: k
    real (r8) :: psi(pcols,pverp)
    real (r8) :: phi(pcols,pverp-1)
    real (r8) :: fdot(pcols,pverp)
    real (r8) :: xx(pcols)
    real (r8) :: fxdot(pcols)
    real (r8) :: fxdd(pcols)
    real (r8) :: psistar(pcols)
    real (r8) :: xxk(pcols,pver)

    do i = 1,ncol
       ! integral of phi
       psi(i,1) = 0._r8
       ! fluxes at boundaries
       flux(i,1) = 0
       flux(i,pverp) = 0._r8
    end do

    ! integral function
    do k = 2,pverp
       do i = 1,ncol
          psi(i,k) = phi(i,k-1)*(xw(i,k)-xw(i,k-1)) + psi(i,k-1)
       end do
    end do

    ! calculate the derivatives for the interpolating polynomial
    call cfdotmc_pro (ncol, xw, psi, fdot)

    ! calculate fluxes at interior pts
    do k = 2,pver
       do i = 1,ncol
          xxk(i,k) = xw(i,k)-vel(i,k)*deltat
       end do
    end do
    do k = 2,pver
       call cfint2(ncol, xw, psi, fdot, xxk(1,k), fxdot, fxdd, psistar)
       do i = 1,ncol
          flux(i,k) = (psi(i,k)-psistar(i))
       end do
    end do

  end subroutine getflx

  !===============================================================================
  subroutine cfint2 (ncol, x, f, fdot, xin, fxdot, fxdd, psistar)

    ! arguments
    integer   , intent(in)  :: ncol                      ! number of colums to process
    real (r8) , intent(in)  :: x(pcols, pverp)
    real (r8) , intent(in)  :: f(pcols, pverp)
    real (r8) , intent(out) :: fdot(pcols, pverp)
    real (r8) , intent(in)  :: xin(pcols)
    real (r8) , intent(out) :: fxdot(pcols)
    real (r8) , intent(out) :: fxdd(pcols)
    real (r8) , intent(out) :: psistar(pcols)

    ! local variables
    integer   :: i
    integer   :: k
    integer   :: intz(pcols)
    real (r8) :: dx
    real (r8) :: s
    real (r8) :: c2
    real (r8) :: c3
    real (r8) :: xx
    real (r8) :: xinf
    real (r8) :: psi1, psi2, psi3, psim
    real (r8) :: cfint
    real (r8) :: cfnew
    real (r8) :: xins(pcols)
    real (r8) :: a, b, c ! the minmod function 
    real (r8) :: minmod ! the minmod function 
    real (r8) :: medan ! the minmod function 

    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do i = 1,ncol
       xins(i) = medan(x(i,1), xin(i), x(i,pverp))
       intz(i) = 0
    end do

    ! first find the interval 
    do k =  1,pverp-1
       do i = 1,ncol
          if ((xins(i)-x(i,k))*(x(i,k+1)-xins(i)).ge.0._r8) then
             intz(i) = k
          endif
       end do
    end do

    do i = 1,ncol
       if (intz(i).eq.0) then
          write(iulog,*) ' interval was not found for col i ', i
          call endrun('DUST_SEDIMENT_MOD:cfint2 -- interval was not found ')
       endif
    end do

    ! now interpolate
    do i = 1,ncol
       k = intz(i)
       dx = (x(i,k+1)-x(i,k))
       s = (f(i,k+1)-f(i,k))/dx
       c2 = (3*s-2*fdot(i,k)-fdot(i,k+1))/dx
       c3 = (fdot(i,k)+fdot(i,k+1)-2*s)/dx**2
       xx = (xins(i)-x(i,k))
       fxdot(i) =  (3*c3*xx + 2*c2)*xx + fdot(i,k)
       fxdd(i) = 6*c3*xx + 2*c2
       cfint = ((c3*xx + c2)*xx + fdot(i,k))*xx + f(i,k)

       ! limit the interpolant
       psi1 = f(i,k)+(f(i,k+1)-f(i,k))*xx/dx
       if (k.eq.1) then
          psi2 = f(i,1)
       else
          psi2 = f(i,k) + (f(i,k)-f(i,k-1))*xx/(x(i,k)-x(i,k-1))
       endif
       if (k+1.eq.pverp) then
          psi3 = f(i,pverp)
       else
          psi3 = f(i,k+1) - (f(i,k+2)-f(i,k+1))*(dx-xx)/(x(i,k+2)-x(i,k+1))
       endif
       psim = medan(psi1, psi2, psi3)
       cfnew = medan(cfint, psi1, psim)
       if (abs(cfnew-cfint)/(abs(cfnew)+abs(cfint)+1.e-36_r8)  .gt..03_r8) then
       endif
       psistar(i) = cfnew
    end do

  end subroutine cfint2

  !===============================================================================
  subroutine cfdotmc_pro (ncol, x, f, fdot)

    ! prototype version; eventually replace with final SPITFIRE scheme
    ! calculate the derivative for the interpolating polynomial multi column version
    ! assumed variable distribution

    !     x1.......x2.......x3.......x4.......x5.......x6     1,pverp points
    !     f1.......f2.......f3.......f4.......f5.......f6     1,pverp points
    !     ...sh1.......sh2......sh3......sh4......sh5....     1,pver points
    !     .........d2.......d3.......d4.......d5.........     2,pver points
    !     .........s2.......s3.......s4.......s5.........     2,pver points
    !     .............dh2......dh3......dh4.............     2,pver-1 points
    !     .............eh2......eh3......eh4.............     2,pver-1 points
    !     ..................e3.......e4..................     3,pver-1 points
    !     .................ppl3......ppl4................     3,pver-1 points
    !     .................ppr3......ppr4................     3,pver-1 points
    !     .................t3........t4..................     3,pver-1 points
    !     ................fdot3.....fdot4................     3,pver-1 points


    ! arguments
    integer   , intent(in)  :: ncol               ! number of colums to process
    real (r8) , intent(in)  :: x(pcols, pverp)
    real (r8) , intent(in)  :: f(pcols, pverp)
    real (r8) , intent(out) :: fdot(pcols, pverp) ! derivative at nodes

    ! local variables
    integer  :: i,k
    real(r8) :: a,b,c            ! work vars
    real(r8) :: s(pcols,pverp)   ! first divided differences at nodes
    real(r8) :: sh(pcols,pverp)  ! first divided differences between nodes
    real(r8) :: d(pcols,pverp)   ! second divided differences at nodes
    real(r8) :: dh(pcols,pverp)  ! second divided differences between nodes
    real(r8) :: e(pcols,pverp)   ! third divided differences at nodes
    real(r8) :: eh(pcols,pverp)  ! third divided differences between nodes
    real(r8) :: pp               ! p prime
    real(r8) :: ppl(pcols,pverp) ! p prime on left
    real(r8) :: ppr(pcols,pverp) ! p prime on right
    real(r8) :: qpl
    real(r8) :: qpr
    real(r8) :: ttt
    real(r8) :: t
    real(r8) :: tmin
    real(r8) :: tmax
    real(r8) :: delxh(pcols,pverp)
    real(r8) :: minmod           ! the minmod function 
    real(r8) :: medan            ! the minmod function 

    minmod(a,b) = 0.5_r8*(sign(1._r8,a) + sign(1._r8,b))*min(abs(a),abs(b))
    medan(a,b,c) = a + minmod(b-a,c-a)

    do k = 1,pver
       ! first divided differences between nodes
       do i = 1, ncol
          delxh(i,k) = (x(i,k+1)-x(i,k))
          sh(i,k) = (f(i,k+1)-f(i,k))/delxh(i,k)
       end do

       ! first and second divided differences at nodes
       if (k.ge.2) then
          do i = 1,ncol
             d(i,k) = (sh(i,k)-sh(i,k-1))/(x(i,k+1)-x(i,k-1))
             s(i,k) = minmod(sh(i,k),sh(i,k-1))
          end do
       endif
    end do

    ! second and third divided diffs between nodes
    do k = 2,pver-1
       do i = 1, ncol
          eh(i,k) = (d(i,k+1)-d(i,k))/(x(i,k+2)-x(i,k-1))
          dh(i,k) = minmod(d(i,k),d(i,k+1))
       end do
    end do

    ! treat the boundaries
    do i = 1,ncol
       e(i,2) = eh(i,2)
       e(i,pver) = eh(i,pver-1)
       !  outside level
       fdot(i,1) = sh(i,1) - d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*(x(i,1)-x(i,3))
       fdot(i,1) = minmod(fdot(i,1),3*sh(i,1))
       fdot(i,pverp) = sh(i,pver) + d(i,pver)*delxh(i,pver) + eh(i,pver-1)*delxh(i,pver)*(x(i,pverp)-x(i,pver-1))
       fdot(i,pverp) = minmod(fdot(i,pverp),3*sh(i,pver))

       ! one in from boundary
       fdot(i,2) = sh(i,1) + d(i,2)*delxh(i,1) - eh(i,2)*delxh(i,1)*delxh(i,2)
       fdot(i,2) = minmod(fdot(i,2),3*s(i,2))
       fdot(i,pver) = sh(i,pver) - d(i,pver)*delxh(i,pver) - eh(i,pver-1)*delxh(i,pver)*delxh(i,pver-1)
       fdot(i,pver) = minmod(fdot(i,pver),3*s(i,pver))
    end do

    do k = 3,pver-1
       do i = 1,ncol
          e(i,k) = minmod(eh(i,k),eh(i,k-1))
       end do
    end do

    do k = 3,pver-1
       do i = 1,ncol
          ! p prime at k-0.5
          ppl(i,k)=sh(i,k-1) + dh(i,k-1)*delxh(i,k-1)  

          ! p prime at k+0.5
          ppr(i,k)=sh(i,k)   - dh(i,k)  *delxh(i,k)
          t = minmod(ppl(i,k),ppr(i,k))

          ! derivate from parabola thru f(i,k-1), f(i,k), and f(i,k+1)
          pp = sh(i,k-1) + d(i,k)*delxh(i,k-1) 

          ! quartic estimate of fdot
          fdot(i,k) = pp - delxh(i,k-1)*delxh(i,k)*(eh(i,k-1)*(x(i,k+2)-x(i,k)) &
               + eh(i,k  )*(x(i,k  )-x(i,k-2)))/(x(i,k+2)-x(i,k-2))

          ! now limit it
          qpl = sh(i,k-1) + delxh(i,k-1)*minmod(d(i,k-1)+ e(i,k-1)*(x(i,k)-x(i,k-2)), &
               d(i,k) - e(i,k)*delxh(i,k))
          qpr = sh(i,k) + delxh(i,k  )*minmod(d(i,k) + e(i,k)*delxh(i,k-1), &
               d(i,k+1)+e(i,k+1)*(x(i,k)-x(i,k+2)))

          fdot(i,k) = medan(fdot(i,k), qpl, qpr)

          ttt = minmod(qpl, qpr)
          tmin = min(0._r8,3*s(i,k),1.5_r8*t,ttt)
          tmax = max(0._r8,3*s(i,k),1.5_r8*t,ttt)
          fdot(i,k) = fdot(i,k) + minmod(tmin-fdot(i,k), tmax-fdot(i,k))
       end do
    end do

  end subroutine cfdotmc_pro

end module oslo_aero_dust_sediment
