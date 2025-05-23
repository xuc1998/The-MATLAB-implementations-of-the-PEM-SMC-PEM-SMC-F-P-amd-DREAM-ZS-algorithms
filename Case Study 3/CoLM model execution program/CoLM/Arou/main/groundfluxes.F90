
 subroutine groundfluxes (zlnd, zsno, hu, ht, hq, &
                          us, vs, tm, qm, rhoair, psrf, &
                          ur, thm, th, thv, tg, qg, dqgdT, htvp, &
                          fsno, sigf, cgrnd, cgrndl, cgrnds, & 
                          taux, tauy, fsena, fevpa, fseng, fevpg, tref, qref, &
                          z0ma, zol, rib, ustar, qstar, tstar, f10m, fm, fh, fq)

!=======================================================================
! this is the main subroutine to execute the calculation of thermal processes
! and surface fluxes
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002
!=======================================================================

  use precision
  use phycon_module, only : cpair,vonkar,grav
  implicit none
 
!----------------------- Dummy argument --------------------------------
  real(r8), INTENT(in) :: &
        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]

        ! atmospherical variables and observational height
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        tm,       &! temperature at agcm reference height [kelvin] [not used]
        qm,       &! specific humidity at agcm reference height [kg/kg]
        rhoair,   &! density air [kg/m3]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]

        fsno,     &! fraction of ground covered by snow
        sigf,     &! fraction of veg cover, excluding snow-covered veg [-]

        ur,       &! wind speed at reference height [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)

        tg,       &! ground surface temperature [K]
        qg,       &! ground specific humidity [kg/kg]
        dqgdT,    &! d(qg)/dT
        htvp       ! latent heat of vapor of water (or sublimation) [j/kg]

  real(r8), INTENT(out) :: &
        taux,     &! wind stress: E-W [kg/m/s**2]
        tauy,     &! wind stress: N-S [kg/m/s**2]
        fsena,    &! sensible heat from canopy height to atmosphere [W/m2]
        fevpa,    &! evapotranspiration from canopy height to atmosphere [mm/s]
        fseng,    &! sensible heat flux from ground [W/m2]
        fevpg,    &! evaporation heat flux from ground [mm/s]
        cgrnd,    &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,   &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,   &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,     &! 2 m height air temperature [kelvin]
        qref,     &! 2 m height air humidity

        z0ma,     &! effective roughness [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        f10m,     &! integral of profile function for momentum at 10m
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq         ! integral of profile function for moisture

!------------------------ LOCAL VARIABLES ------------------------------
  integer i
  integer niters, &! maximum number of iterations for surface temperature
       iter,      &! iteration index
       nmozsgn     ! number of times moz changes sign

  real(r8) :: &
       beta,      &! coefficient of conective velocity [-]
       displax,   &! zero-displacement height [m]
       dth,       &! diff of virtual temp. between ref. height and surface
       dqh,       &! diff of humidity between ref. height and surface
       dthv,      &! diff of vir. poten. temp. between ref. height and surface
       obu,       &! monin-obukhov length (m)
       obuold,    &! monin-obukhov length from previous iteration
       ram,       &! aerodynamical resistance [s/m]
       rah,       &! thermal resistance [s/m]
       raw,       &! moisture resistance [s/m]
       raih,      &! temporary variable [kg/m2/s]
       raiw,      &! temporary variable [kg/m2/s]
       temp1,     &! relation for potential temperature profile
       temp2,     &! relation for specific humidity profile
       temp12m,   &! relation for temperature at 2m
       temp22m,   &! relation for specific humidity at 2m
       thvstar,   &! virtual potential temperature scaling parameter
       um,        &! wind speed including the stablity effect [m/s]
       wc,        &! convective velocity [m/s]
       wc2,       &! wc**2
       zeta,      &! dimensionless height used in Monin-Obukhov theory
       zii,       &! convective boundary height [m]
       zldis,     &! reference height "minus" zero displacement heght [m]
       z0mg,      &! roughness length over ground, momentum [m]
       z0hg,      &! roughness length over ground, sensible heat [m]
       z0qg        ! roughness length over ground, latent heat [m]

!----------------------- Dummy argument --------------------------------
! initial roughness length
      if(fsno > 0.)then
         z0mg = zsno; z0hg = z0mg; z0qg = z0mg
      else
         z0mg = zlnd; z0hg = z0mg; z0qg = z0mg
      endif

! potential temperatur at the reference height
      beta = 1.      ! -  (in computing W_*)
      zii = 1000.    ! m  (pbl height)
      z0ma = z0mg

!-----------------------------------------------------------------------
!     Compute sensible and latent fluxes and their derivatives with respect 
!     to ground temperature using ground temperatures from previous time step.
!-----------------------------------------------------------------------
! Initialization variables
      nmozsgn = 0
      obuold = 0.

      dth   = thm-tg
      dqh   = qm-qg
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
      zldis = hu-0.

      call moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)
 
! Evaluated stability-dependent variables using moz from prior iteration
      niters=6

      !----------------------------------------------------------------
      do iter = 1, niters         ! begin stability iteration
      !----------------------------------------------------------------
         displax = 0.
         call moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                        ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq)

         tstar = temp1*dth
         qstar = temp2*dqh

         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         thvstar=tstar+0.61*th*qstar
         zeta=zldis*vonkar*grav*thvstar/(ustar**2*thv)
         if(zeta >= 0.) then     !stable
           zeta = min(2.,max(zeta,1.e-6))
         else                    !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
         endif
         obu = zldis/zeta

         if(zeta >= 0.)then
           um = max(ur,0.1)
         else
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
          wc2 = beta*beta*(wc*wc)
           um = sqrt(ur*ur+wc2)
         endif

         if (obuold*obu < 0.) nmozsgn = nmozsgn+1
         if(nmozsgn >= 4) EXIT

         obuold = obu

      !----------------------------------------------------------------
      enddo                       ! end stability iteration
      !----------------------------------------------------------------

! Get derivative of fluxes with repect to ground temperature
      ram    = 1./(ustar*ustar/um)
      rah    = 1./(temp1*ustar) 
      raw    = 1./(temp2*ustar) 

      raih   = (1.-sigf)*rhoair*cpair/rah
      raiw   = (1.-sigf)*rhoair/raw          
      cgrnds = raih
      cgrndl = raiw*dqgdT
      cgrnd  = cgrnds + htvp*cgrndl

      zol = zeta
      rib = min(5.,zeta*ustar**2/(vonkar*temp1*um**2))

! surface fluxes of momentum, sensible and latent 
! using ground temperatures from previous time step
      taux   = -(1.-sigf)*rhoair*us/ram        
      tauy   = -(1.-sigf)*rhoair*vs/ram
      fseng  = -raih*dth
      fevpg  = -raiw*dqh 

      fsena  = fseng
      fevpa  = fevpg

! 2 m height air temperature
      tref   = (1.-sigf)*(thm + temp1*dth * (1./temp12m - 1./temp1))
      qref   = (1.-sigf)*( qm + temp2*dqh * (1./temp22m - 1./temp2))

 end subroutine groundfluxes
