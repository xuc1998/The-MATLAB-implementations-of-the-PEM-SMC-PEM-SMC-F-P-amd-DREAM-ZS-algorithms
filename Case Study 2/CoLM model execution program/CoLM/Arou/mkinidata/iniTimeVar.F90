
#include <define.h>

 subroutine iniTimeVar(nl_soil,maxsnl,itypwat&
	              ,porsl,albsol,z0m,chil,ref,tran&
                      ,z,dz,tss,wliq,wice&
                      ,tg,tlsun,tlsha,ldew,sag,scv&   
                      ,snowdp,fveg,fsno,sigf,green,lai,sai,coszen&
                      ,albg,albv,alb,ssun,ssha,thermk,extkb,extkd&
                      ,trad,tref,qref,rst,emis,z0ma,zol,rib&
                      ,ustar,qstar,tstar,fm,fh,fq&
#if(defined SOILINI)
                      ,nl_soil_ini,soil_z,soil_t,soil_w,snow_d)
#else
                      )
#endif

!=======================================================================
! Original author : Yongjiu Dai, 08/30/2002, revised February 2004
!=======================================================================

  use precision
  use phycon_module, only : tfrz
  implicit none 

  integer, INTENT(in) ::        &! 
        nl_soil,                &! soil layer number
        maxsnl,                 &! maximum snow layer number
        itypwat                  ! index for land cover type [-]

  real(r8), INTENT(in) ::       &!
        fveg,                   &! fraction of vegetation cover
        green,                  &! leaf greenness
        lai,                    &! leaf area index
        sai,                    &! stem area index
        coszen,                 &! cosine of solar zenith angle
	albsol,                 &! soil albedo for different coloured soils [-]
	z0m,                    &! aerodynamic roughness length [m]
        chil,                   &! leaf angle distribution factor
        ref (2,2),              &! leaf reflectance (iw=iband, il=life and dead)
        tran(2,2),              &! leaf transmittance (iw=iband, il=life and dead)
        porsl(1:nl_soil)         ! porosity of soil

#if(defined SOILINI)
  integer, INTENT(in) :: nl_soil_ini
  real(r8), INTENT(in) ::       &!
        soil_z(nl_soil_ini),    &! soil layer depth for initial (m)
        soil_t(nl_soil_ini),    &! soil temperature from initial file (K)
        soil_w(nl_soil_ini),    &! soil wetness from initial file (-)
        snow_d                   ! snow depth (m)
#endif

  real(r8), INTENT(inout) ::    &!
        z (maxsnl+1:nl_soil),   &! node depth [m]
        dz(maxsnl+1:nl_soil)     ! interface depth [m]

  real(r8), INTENT(out) ::      &!
        tss (maxsnl+1:nl_soil), &! soil temperature [K]
        wliq(maxsnl+1:nl_soil), &! liquid water in layers [kg/m2]
        wice(maxsnl+1:nl_soil), &! ice lens in layers [kg/m2]
        tg,                     &! ground surface temperature [K]
        tlsun,                  &! sunlit leaf temperature [K]
        tlsha,                  &! shaded leaf temperature [K]
        ldew,                   &! depth of water on foliage [mm]
        sag,                    &! non dimensional snow age [-]
        scv,                    &! snow cover, water equivalent [mm]
        snowdp,                 &! snow depth [meter]
        fsno,                   &! fraction of snow cover on ground
        sigf,                   &! fraction of veg cover, excluding snow-covered veg [-]

        albg(2,2),              &! albedo, ground [-]
        albv(2,2),              &! albedo, vegetation [-]
        alb (2,2),              &! averaged albedo [-]
        ssun(2,2),              &! sunlit canopy absorption for solar radiation
        ssha(2,2),              &! shaded canopy absorption for solar radiation
        thermk,                 &! canopy gap fraction for tir radiation
        extkb,                  &! (k, g(mu)/mu) direct solar extinction coefficient
        extkd,                  &! diffuse and scattered diffuse PAR extinction coefficient

                    ! Additional variables required by reginal model (WRF & RSM) 
        trad,                   &! radiative temperature of surface [K]
        tref,                   &! 2 m height air temperature [kelvin]
        qref,                   &! 2 m height air specific humidity
        rst,                    &! canopy stomatal resistance (s/m)
        emis,                   &! averaged bulk surface emissivity
        z0ma,                   &! effective roughness [m]
        zol,                    &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,                    &! bulk Richardson number in surface layer
        ustar,                  &! u* in similarity theory [m/s]
        qstar,                  &! q* in similarity theory [kg/kg]
        tstar,                  &! t* in similarity theory [K]
        fm,                     &! integral of profile function for momentum
        fh,                     &! integral of profile function for heat
        fq                       ! integral of profile function for moisture

        integer j, snl                      
        real(r8) wet(nl_soil), wt, ssw, oro, rhosno_ini, a
!-----------------------------------------------------------------------

  if(itypwat <= 5)then ! land grid
     rhosno_ini = 250.
#if(defined SOILINI)
     do j = 1, nl_soil
        call polint(soil_z,soil_t,nl_soil_ini,z(j),tss(j))
        call polint(soil_z,soil_w,nl_soil_ini,z(j),wet(j))
        a = min(soil_t(1),soil_t(2),soil_t(3))-5.
        tss(j) = max(tss(j), a)
        a = max(soil_t(1),soil_t(2),soil_t(3))+5.
        tss(j) = min(tss(j), a)

        a = min(soil_w(1),soil_w(2),soil_w(3))
	wet(j) = max(wet(j), a, 0.1)
        a = max(soil_w(1),soil_w(2),soil_w(3))
	wet(j) = min(wet(j), a, 0.5)

        if(tss(j).ge.tfrz)then
           wliq(j) = wet(j)*dz(j)*1000.
!          wliq(j) = porsl(j)*wet(j)*dz(j)*1000.
           wice(j) = 0.
        else
           wliq(j) = 0.
           wice(j) = wet(j)*dz(j)*1000.
!          wliq(j) = porsl(j)*wet(j)*dz(j)*1000.
        endif
     enddo

     snowdp = snow_d
     sag    = 0.
     scv    = snowdp*rhosno_ini

     call snowfraction (fveg,z0m,snowdp,wt,sigf,fsno)
     call snow_ini (itypwat,maxsnl,snowdp,snl,z,dz)

     if(snl.lt.0)then
        do j = snl+1, 0
           tss(j) = min(tfrz-1., tss(1))
           wliq(j) = 0.
           wice(j) = dz(j)*rhosno_ini         ! m * kg m-3 = kg m-2
        enddo
     endif

     if(snl>maxsnl)then
        tss (maxsnl+1:snl) = -999.
        wice(maxsnl+1:snl) = 0.
        wliq(maxsnl+1:snl) = 0.
        z   (maxsnl+1:snl) = 0.
        dz  (maxsnl+1:snl) = 0.
     endif

     ldew  = 0.
     tlsun = tss(1)
     tlsha = tss(1)
     tg    = tss(1)
#else
! soil temperature and water content
     do j = 1, nl_soil
        if(itypwat==3)then ! land ice 
           tss(j) = 253.
           wliq(j) = 0.
          wice(j) = dz(j)*1000.
        else
           tss(j) = 270
           wliq(j) = dz(j)*porsl(j)*1000.
           wice(j) = 0.
        endif
     enddo
   
! snow temperature and water content
     tss(maxsnl+1:0) = -999.
     wice(maxsnl+1:0) = 0.
     wliq(maxsnl+1:0) = 0.
     z (maxsnl+1:0) = 0.
     dz(maxsnl+1:0) = 0.

!     print*,'iniTimeVar.F90-test-myinitial-begin'
     tss(maxsnl+1:nl_soil)=(/0.00,0.00,0.00,0.00,0.00,&
                            276.678,284.10,285.17,283.62,&
                            282.94,279.823,279.234,275.61,275.042,273.97/)
     wice(maxsnl+1:nl_soil)=(/0.00,0.00,0.00,0.00,0.00,&
                             0.00,0.00,0.00,0.00,0.00,0.00,&
                             0.00,0.00,0.00,0.00/)
     wliq(maxsnl+1:nl_soil)=(/0.00,0.00,0.00,0.00,0.00,&
                             7.28,7.28,7.3,16.92,17.1,&
                             26.7,19.4,38.4,29.2,24.4/)
     z(maxsnl+1:nl_soil)=(/0.000,0.000,0.000,0.000,0.000,&
                          0.007,0.027,0.062,0.118,0.212,&
                          0.366,0.619,1.038,1.727,2.864/)
     dz(maxsnl+1:nl_soil)=(/0.00,0.00,0.00,0.0,0.00,&
                           0.017,0.027,0.045,0.074,0.123,&
                           0.203,0.335,0.553,0.913,1.136/)
!     print*,'iniTimeVar.F90-test-tss6',tss(6)


     sigf   = fveg
     fsno   = 0.160
     ldew   = 0.407
     scv    = 2.662
     sag    = 0.171
     snowdp = 0.019
     tlsun  = tss(1)
!     tlsun=272.378
     tlsha  = tss(1)
!     tlsha=262.378
     tg     = tss(1)
!      tg=272.5879
#endif

! surface albedo
     ssw = min(1.,1.e-3*wliq(1)/dz(1))
     call albland (itypwat,albsol,chil,ref,tran,&
                   fveg,green,lai,sai,coszen,wt,fsno,scv,sag,ssw,tg,&
                   alb,albg,albv,ssun,ssha,thermk,extkb,extkd)
  else                 ! ocean grid
     tss(:) = 300.
     wice(:) = 0.
     wliq(:) = 1000.
     z (maxsnl+1:0) = 0.
     dz(maxsnl+1:0) = 0.
     sigf   = 0.
     fsno   = 0.
     ldew   = 0.
     scv    = 0.
     sag    = 0.
     snowdp = 0.
     tlsun  = 300.
     tlsha  = 300.
     tg     = 300.

     oro = 0
     call albocean (oro,scv,coszen,alb)
     albg(:,:) = alb(:,:)
     albv(:,:) = 0.0
     ssun(:,:) = 0.0
     ssha(:,:) = 0.0
     thermk = 0.0
     extkb = 0.0
     extkd = 0.0
  endif

! Additional variables required by reginal model (WRF & RSM)
! totally arbitrarily assigned here
  trad  = tg      
  tref  = tg      
  qref  = 0.3     
  rst   = 1.e36   
  emis  = 1.0     
  z0ma  = 0.01    
  zol   = -1.0    
  rib   = -0.1    
  ustar = 0.25    
  qstar = 0.001   
  tstar = -1.5    
  fm    = alog(30.)  
  fh    = alog(30.)  
  fq    = alog(30.)  

 end subroutine iniTimeVar
!-----------------------------------------------------------------------
! EOP


  subroutine snow_ini(itypwat,maxsnl,snowdp,snl,z,dz)

! Snow spatial discretization initially

  use precision
  implicit none

  integer, intent(in) :: maxsnl  ! maximum of snow layers
  integer, intent(in) :: itypwat ! index for land cover type [-]
  real(r8), intent(in) :: snowdp ! snow depth [m]
  real(r8), intent(out) :: z (maxsnl+1:0) ! node depth [m]
  real(r8), intent(out) :: dz(maxsnl+1:0) ! layer thickness [m]
  integer, intent(out) :: snl ! number of snow layer
  real(r8) zi
  integer i
!-----------------------------------------------------------------------

  dz(:0) = 0.
  z(:0) = 0.
  snl = 0
  if(itypwat.le.3)then ! non water bodies

     if(snowdp.lt.0.01)then
        snl = 0
     else
        if(snowdp>=0.01 .and. snowdp<=0.03)then
           snl = -1
           dz(0)  = snowdp
        else if(snowdp>0.03 .and. snowdp<=0.04)then
           snl = -2
           dz(-1) = snowdp/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.04 .and. snowdp<=0.07)then
           snl = -2
           dz(-1) = 0.02
           dz( 0) = snowdp - dz(-1)
        else if(snowdp>0.07 .and. snowdp<=0.12)then
           snl = -3
           dz(-2) = 0.02
           dz(-1) = (snowdp - 0.02)/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.12 .and. snowdp<=0.18)then
           snl = -3
           dz(-2) = 0.02
           dz(-1) = 0.05
           dz( 0) = snowdp - dz(-2) - dz(-1)
        else if(snowdp>0.18 .and. snowdp<=0.29)then
           snl = -4
           dz(-3) = 0.02
           dz(-2) = 0.05
           dz(-1) = (snowdp - dz(-3) - dz(-2))/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.29 .and. snowdp<=0.41)then
           snl = -4
           dz(-3) = 0.02
           dz(-2) = 0.05
           dz(-1) = 0.11
           dz( 0) = snowdp - dz(-3) - dz(-2) - dz(-1)
        else if(snowdp>0.41 .and. snowdp<=0.64)then
           snl = -5
           dz(-4) = 0.02
           dz(-3) = 0.05
           dz(-2) = 0.11
           dz(-1) = (snowdp - dz(-4) - dz(-3) - dz(-2))/2.
           dz( 0) = dz(-1)
        else if(snowdp>0.64)then
           snl = -5
           dz(-4) = 0.02
           dz(-3) = 0.05
           dz(-2) = 0.11
           dz(-1) = 0.23
           dz( 0) = snowdp - dz(-4) - dz(-3) - dz(-2) - dz(-1)
        endif

        zi = 0.
        do i = 0, snl+1, -1
           z(i) = zi - dz(i)/2.
           zi = -zi-dz(i)
        enddo
     endif

  endif

  end subroutine snow_ini
!-----------------------------------------------------------------------
! EOP


  subroutine polint(xa,ya,n,x,y)

! Given arrays xa and ya, each of length n, and gi
! value y, and an error estimate dy. If P (x) is the p
! P (xa(i)) = ya(i), i = 1, . . . , n, then the returned value
! (from: "Numerical Recipes")

  use precision
  implicit none
  integer n,NMAX
  real(r8) dy,x,y,xa(n),ya(n)
  parameter (NMAX=10)      !Largest anticipated val
  integer i,m,ns
  real(r8) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)

  ns=1
  dif=abs(x-xa(1))

  do i=1,n       !Here we find the index ns of the closest table entry,
     dift=abs(x-xa(i))
     if(dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)  !and initialize the tableau of c's and d's.
     d(i)=ya(i)
  enddo

  y=ya(ns)       !This is the initial approximation to y.
  ns=ns-1

  do m=1,n-1  !For each column of the tableau,
     do i=1,n-m   !we loop over the current c's and d's and update them.
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) pause 'failure in polint'  !two input xa's are identical.
           den=w/den
           d(i)=hp*den    !Here the c's and d's are updated.
           c(i)=ho*den
     enddo
     if(2*ns.lt.n-m)then  !After each column in the tableau is completed, we decide
        dy=c(ns+1)        !which correction, c or d, we want to add to our accumulating
     else                 !value of y, i.e., which path to take through
        dy=d(ns)          !the tableau-forking up or down. We do this in such a
        ns=ns-1           !way as to take the most "straight line" route through the
     endif                !tableau to its apex, updating ns accordingly to keep track
     y=y+dy               !of where we are. This route keeps the partial approximations
  enddo                   !centered (insofar as possible) on the target x. T he
                          !last dy added is thus the error indication.
  end subroutine polint
!-----------------------------------------------------------------------
! EOP

