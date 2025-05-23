
 subroutine moninobuk(hu,ht,hq,displa,z0m,z0h,z0q,obu,um,&
		      ustar,temp1,temp2,temp12m,temp22m,f10m,fm,fh,fq)

! ======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! calculation of friction velocity, relation for potential temperatur
! and humidity profiles of surface boundary layer.
! the scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
! ======================================================================

  use precision
  use phycon_module, only : vonkar
  implicit none

! ---------------------- dummy argument --------------------------------

  real(r8), INTENT(in) :: hu       ! observational height of wind [m]
  real(r8), INTENT(in) :: ht       ! observational height of temperature [m]
  real(r8), INTENT(in) :: hq       ! observational height of humidity [m]
  real(r8), INTENT(in) :: displa   ! displacement height [m]
  real(r8), INTENT(in) :: z0m      ! roughness length, momentum [m]
  real(r8), INTENT(in) :: z0h      ! roughness length, sensible heat [m]
  real(r8), INTENT(in) :: z0q      ! roughness length, latent heat [m]
  real(r8), INTENT(in) :: obu      ! monin-obukhov length (m)
  real(r8), INTENT(in) :: um       ! wind speed including the stablity effect [m/s]

  real(r8), INTENT(out) :: ustar   ! friction velocity [m/s]
  real(r8), INTENT(out) :: temp1   ! relation for potential temperature profile
  real(r8), INTENT(out) :: temp2   ! relation for specific humidity profile
  real(r8), INTENT(out) :: temp12m ! relation for temperature at 2m
  real(r8), INTENT(out) :: temp22m ! relation for specific humidity at 2m
  real(r8), INTENT(out) :: f10m    ! integral of profile function for momentum at 10m
  real(r8), INTENT(out) :: fm      ! integral of profile function for momentum
  real(r8), INTENT(out) :: fh      ! integral of profile function for heat 
  real(r8), INTENT(out) :: fq      ! integral of profile function for moisture

!------------------------ local variables ------------------------------

  real(r8) zldis  ! reference height "minus" zero displacement heght [m]
  real(r8) psi    ! stability function for unstable case
  real(r8) zetam  ! transition point of flux-gradient relation (wind profile)
  real(r8) zetat  ! transition point of flux-gradient relation (temp. profile)
  real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory

!-----------------------------------------------------------------------
! adjustment factors for unstable (moz < 0) or stable (moz > 0) conditions.

! wind profile
        zldis=hu-displa
        zeta=zldis/obu
        zetam=1.574
        if(zeta < -zetam)then           ! zeta < -1
          fm    = log(-zetam*obu/z0m) - psi(1,-zetam) &
                + psi(1,z0m/obu) + 1.14*((-zeta)**0.333-(zetam)**0.333) 
          ustar = vonkar*um / fm 
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fm    = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu) 
          ustar = vonkar*um / fm
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fm    = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu 
          ustar = vonkar*um / fm
        else                            !  1 < zeta, phi=5+zeta
          fm    = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.) 
          ustar = vonkar*um / fm
        endif
        
        ! for 10 meter wind-velocity
        zldis=10.+z0m
        zeta=zldis/obu
        zetam=1.574
        if(zeta < -zetam)then           ! zeta < -1
          f10m  = log(-zetam*obu/z0m) - psi(1,-zetam) &
                + psi(1,z0m/obu) + 1.14*((-zeta)**0.333-(zetam)**0.333) 
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          f10m  = log(zldis/z0m) - psi(1,zeta) + psi(1,z0m/obu) 
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          f10m  = log(zldis/z0m) + 5.*zeta - 5.*z0m/obu 
        else                            !  1 < zeta, phi=5+zeta
          f10m  = log(obu/z0m) + 5. - 5.*z0m/obu + (5.*log(zeta)+zeta-1.) 
        endif

! temperature profile
        zldis=ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fh    = log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)) 
          temp1 = vonkar / fh 
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fh    = log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu) 
          temp1 = vonkar / fh 
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fh    = log(zldis/z0h) + 5.*zeta - 5.*z0h/obu 
          temp1 = vonkar / fh
        else                            !  1 < zeta, phi=5+zeta
          fh    = log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.) 
          temp1 = vonkar / fh
        endif

        ! for 2 meter screen temperature
        zldis=2.+z0h  ! ht-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          temp12m = vonkar / ( log(-zetat*obu/z0h)-psi(2,-zetat) &
                + psi(2,z0h/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)) )
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          temp12m = vonkar / ( log(zldis/z0h) - psi(2,zeta) + psi(2,z0h/obu) )
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          temp12m = vonkar / ( log(zldis/z0h) + 5.*zeta - 5.*z0h/obu )
        else                            !  1 < zeta, phi=5+zeta
          temp12m = vonkar / ( log(obu/z0h) + 5. - 5.*z0h/obu + (5.*log(zeta)+zeta-1.) )
        endif

! humidity profile
        zldis=hq-displa
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
          fq    = log(-zetat*obu/z0q) - psi(2,-zetat) &
                + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)) 
          temp2 = vonkar / fq
        else if(zeta < 0.)then          ! -1 <= zeta < 0
          fq    = log(zldis/z0q) - psi(2,zeta) + psi(2,z0q/obu) 
          temp2 = vonkar / fq
        else if(zeta <= 1.)then         !  0 <= ztea <= 1
          fq    = log(zldis/z0q) + 5.*zeta - 5.*z0q/obu 
          temp2 = vonkar / fq
        else                            !  1 < zeta, phi=5+zeta
          fq    = log(obu/z0q) + 5. - 5.*z0q/obu + (5.*log(zeta)+zeta-1.) 
          temp2 = vonkar / fq
        endif

      ! for 2 meter screen humidity
        zldis=2.+z0h
        zeta=zldis/obu
        zetat=0.465
        if(zeta < -zetat)then           ! zeta < -1
           temp22m = vonkar/(log(-zetat*obu/z0q)-psi(2,-zetat) &
                 + psi(2,z0q/obu) + 0.8*((zetat)**(-0.333)-(-zeta)**(-0.333)))
        elseif (zeta < 0.) then         ! -1 <= zeta < 0
           temp22m = vonkar/(log(zldis/z0q)-psi(2,zeta)+psi(2,z0q/obu))
        else if (zeta <= 1.) then       !  0 <= zeta <= 1
           temp22m = vonkar/(log(zldis/z0q)+5.*zeta-5.*z0q/obu)
        else                            ! 1 < zeta, phi=5+zeta
           temp22m = vonkar/(log(obu/z0q)+5.-5.*z0q/obu+(5.*log(zeta)+zeta-1.))
        endif

 end subroutine moninobuk




 subroutine moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0m,um,obu)

! ======================================================================
! Original author : Yongjiu Dai, September 15, 1999
!
! initialzation of Monin-Obukhov length,
! the scheme is based on the work of Zeng et al. (1998):
! Intercomparison of bulk aerodynamic algorithms for the computation
! of sea surface fluxes using TOGA CORE and TAO data. J. Climate, Vol. 11: 2628-2644
! ======================================================================

   use precision
   use phycon_module, only : grav, vonkar
   implicit none

! Dummy argument
  real(r8), INTENT(in) :: ur    ! wind speed at reference height [m/s]
  real(r8), INTENT(in) :: thm   ! intermediate variable (tm+0.0098*ht)
  real(r8), INTENT(in) :: th    ! potential temperature [kelvin]
  real(r8), INTENT(in) :: thv   ! virtual potential temperature (kelvin)
  real(r8), INTENT(in) :: dth   ! diff of virtual temp. between ref. height and surface
  real(r8), INTENT(in) :: dthv  ! diff of vir. poten. temp. between ref. height and surface
  real(r8), INTENT(in) :: dqh   ! diff of humidity between ref. height and surface
  real(r8), INTENT(in) :: zldis ! reference height "minus" zero displacement heght [m]
  real(r8), INTENT(in) :: z0m   ! roughness length, momentum [m]

  real(r8), INTENT(out) :: um   ! wind speed including the stablity effect [m/s]
  real(r8), INTENT(out) :: obu  ! monin-obukhov length (m)

! Local       
  real(r8) wc     ! convective velocity [m/s]
  real(r8) rib    ! bulk Richardson number
  real(r8) zeta   ! dimensionless height used in Monin-Obukhov theory
  real(r8) ustar  ! friction velocity [m/s]

!-----------------------------------------------------------------------
! Initial values of u* and convective velocity

        ustar=0.06
        wc=0.5
        if(dthv >= 0.)then
          um=max(ur,0.1)
        else
          um=sqrt(ur*ur+wc*wc)
        endif

        rib=grav*zldis*dthv/(thv*um*um)

        if(rib >= 0.)then      ! neutral or stable
          zeta = rib*log(zldis/z0m)/(1.-5.*min(rib,0.19))
          zeta = min(2.,max(zeta,1.e-6))
        else                   ! unstable
          zeta = rib*log(zldis/z0m)
          zeta = max(-100.,min(zeta,-1.e-6))
        endif
        obu=zldis/zeta

 end subroutine moninobukini



 function psi(k,zeta)

!=======================================================================
! stability function for rib < 0

  use precision
  implicit none
      
  integer k    
  real(r8) psi   ! stability function for unstable case
  real(r8) zeta  ! dimensionless height used in Monin-Obukhov theory
  real(r8) chik  ! 

      chik = (1.-16.*zeta)**0.25
      if(k == 1)then
        psi = 2.*log((1.+chik)*0.5)+log((1.+chik*chik)*0.5)-2.*atan(chik)+2.*atan(1.)
      else
        psi = 2.*log((1.+chik*chik)*0.5)
      endif

 end function psi
