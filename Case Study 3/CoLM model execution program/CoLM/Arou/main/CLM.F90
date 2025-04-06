#include <define.h>

#if(defined offline)
  PROGRAM CLM
! ======================================================================
! The Common Land Model was developed in cooperation with
!     Beijing Normal University and IRI-Columbia University (Dai)
!     Georgia Institute of Technology                       (Dickinson)
!     National Center for Atmospheric Research              (Bonan, Oleson)
!     University of Arizona                                 (Zeng)
!     University of Texas at Austin                         (Yang)
!     GSFC/NASA                                             (Houser, Bosilovich)
!     COLA                                                  (Dirmeyer, Schlosser)
!     Colorado State University at Fort Collins             (Denning, Baker)
!
! Reference: 
!     [1] Dai et al., 2003: The Common Land Model (CLM). 
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. J. Climate
!
!     Author: Yongjiu Dai, January 2004
! ======================================================================

      use precision
      implicit none

#include <paramodel.h>

! ----------------local variables ---------------------------------
      integer, parameter :: nl_soil  = nl_soil_  ! number of soil layers
      integer, parameter :: maxsnl   = maxsnl_   ! max number of snow layers
      integer, parameter :: nftune   = nftune_   ! number of clm tunable constants
      integer, parameter :: nfcon    = nfcon_    ! number of time constant variables
      integer, parameter :: nforc    = nforc_    ! number of forcing variables
      integer, parameter :: nfvar    = nfvar_    ! number of time varying variables
      integer, parameter :: nfldv    = nfldv_    ! number of output variables
      integer, parameter :: nflai    = nflai_    ! number of leaf time varying variables
      integer, parameter :: maxpatch = maxpatch_ ! max number of patches in a grid

      real    :: deltim                          ! time step (senconds)
      integer :: mstep                           ! model step for simulation [-]
                                                 !
      integer :: lusrf                           ! logical unit number of surface data
      integer :: lulai                           ! logical unit number of LAI data
      integer :: lumet                           ! logical unit number of meteorological forcing
      integer :: lhistTimeConst                  ! logical unit number of restart time-invariant file
      integer :: lhistTimeVar                    ! logical unit number of restart time-varying file
      integer :: luout                           ! logical unit number of output
      integer :: lon_points                      ! number of longitude points on model grid
      integer :: lat_points                      ! number of latitude points on model grid
      integer :: numpatch                        ! total number of patches of grids
      integer :: lsao                       !my SA analysis file number
      integer :: lsao2
      integer :: lsao3
      integer :: lsao4
      character(LEN=256):: fsao             !my SA analysis file name
      character(LEN=256):: fsao2
      character(LEN=256):: fsao3
      character(LEN=256):: fsao4
                                                 !
      character(LEN=256) :: site                 ! site name
      character(LEN=256) :: fsurdat              ! file name of surface data
      character(LEN=256) :: flaidat              ! file name of time-varying vegetation data
      character(LEN=256) :: fmetdat              ! file name meteorological data
      character(LEN=256) :: fhistTimeConst       ! file name of restart time-invariant file
      character(LEN=256) :: fhistTimeVar         ! file name of restart time-varying file
      character(LEN=256) :: foutdat              ! file name of output file
                                                 !
#if(!defined EcoDynamics)
      real(r8), allocatable :: mflai(:,:,:,:)    ! time varying vegetation parameters
      real(r8), allocatable :: flai(:,:)         ! time varying vegetation parameters
#endif
      real(r8), allocatable :: forcxy(:,:,:)     ! atmospheric forcing
      real(r8), allocatable :: fldxy(:,:,:)      ! output fluxes in 2-dimension
                                                 !
      integer   idate(3)                         ! calendar (year, julian day, seconds)
      integer,  allocatable :: numpatch_lat(:)   ! number of patches of grids at lon. strip
      integer,  allocatable :: ixy_patch(:)      ! patch longitude index
      integer,  allocatable :: jxy_patch(:)      ! patch latitude index
      integer,  allocatable :: mxy_patch(:)      ! patch subgrid index of lnd point
      real(r8), allocatable :: wtxy_patch(:)     ! patch weight
      real(r8), allocatable :: fcon(:,:)         ! time constant variables
      real(r8), allocatable :: forc(:,:)         ! forcing variables
      real(r8), allocatable :: fvar(:,:)         ! time varying variables
      real(r8), allocatable :: fldv(:,:)         ! output fluxes
      real(r8)  ftune(nftune)                    ! clm tunable constants
      real(r8), allocatable :: rstfac(:,:)       !factor of water stress
      real(r8), allocatable :: rstfac_all(:,:)   ! all time factor of water      
      real(r8), allocatable :: oro(:)            ! ocean(0)/seaice(2)/ flag
      real(r8), allocatable :: a(:,:,:)          !
      real(r8), allocatable :: fldxy_r(:,:,:)    ! output fluxes in 2-dimension
      integer, allocatable :: itypwat(:)         ! land water type
      real(r8), allocatable :: SA_LE(:)      ! SA analysis output matrix
      real(r8), allocatable :: SA_NEE(:)
      real(r8), allocatable :: SA_H(:)
      real(r8), allocatable :: SA_SWC(:)
                                                 !
      logical :: lwrite                          ! true: write out frequency
      logical :: doalb                           ! true if time for surface albedo calculation
      logical :: dolai                           ! true if time for time-varying vegetation paramter
      logical :: dosst                           ! true if time for update sst/ice/snow
                                                 !
      integer :: istep                           ! looping step
      integer :: i,j,k,l,m                       ! looping indices
      integer :: nac                             ! number of accumulation
      integer :: idate_p(3)                      ! current model calendar 
      integer :: luout2
      integer :: lhistTimeVar2
      integer :: lmyvariance
      luout2=17
      lhistTimeVar2=16
      lmyvariance=11   
      lsao=22
      lsao2=23
!      lsao3=24
      lsao4=25
      fsao='../output_LE.txt'
      fsao2='../output_NEE.txt'
!      fsao3='../output_H.txt'
      fsao4='../output_RSM.txt'

      namelist /clmexp/ site,                   &!1
                        flaidat,                &!2
                        fmetdat,                &!3
                        fhistTimeConst,         &!4
                        fhistTimeVar,           &!5
                        foutdat,                &!6
                        lhistTimeConst,         &!7
                        lhistTimeVar,           &!8
                        lulai,                  &!9
                        lumet,                  &!10
                        luout,                  &!11
                        lon_points,             &!12
                        lat_points,             &!13
                        numpatch,               &!14
                        deltim,                 &!15
                        mstep                    !16
! ======================================================================

    ! define the run and open files (for off-line use)
      read(5,clmexp) 

      allocate (numpatch_lat(lat_points))
      allocate (ixy_patch(numpatch))
      allocate (jxy_patch(numpatch))
      allocate (mxy_patch(numpatch))
      allocate (wtxy_patch(numpatch))
      allocate (fcon(numpatch,nfcon))
      allocate (forc(numpatch,nforc))
      allocate (fvar(numpatch,nfvar))
      allocate (fldv(numpatch,nfldv))
      allocate (rstfac(numpatch,1))
      allocate (rstfac_all(mstep,numpatch))
#if(!defined EcoDynamics)
      allocate (mflai(lon_points,lat_points,maxpatch,nflai))
      allocate (flai(numpatch,4))
#endif
      allocate (forcxy(lon_points,lat_points,nforc))
      allocate (fldxy(lon_points,lat_points,nfldv))
      allocate (fldxy_r(lon_points,lat_points,nfldv))
      allocate (a(lon_points,lat_points,nfldv))
      allocate (oro(numpatch)) 
      allocate (itypwat(numpatch)) 
      allocate (SA_LE(mstep))
      allocate (SA_NEE(mstep))
      allocate (SA_H(mstep))
      allocate (SA_SWC(mstep))

#if(!defined EcoDynamics)
      OPEN(unit=lulai,file=flaidat,form='unformatted',&
                           status='old',action='read')
#endif
    ! Open for meteorological forcing data
      OPEN(unit=lumet,file=fmetdat,form='formatted',&
                           status='old',action='read')
     
    ! Open for model time invariant constant data
      OPEN(unit=lhistTimeConst,file=fhistTimeConst,form='unformatted',&
                           status='unknown',action='read')
    
    ! Open for model time varying data (model state variables)
      OPEN(unit=lhistTimeVar,file=fhistTimeVar,form='unformatted',&
                           status='unknown',action='read')
    ! Open my interest variance output file
!      OPEN(unit=lmyvariance,file='/home/u120220909911/RU-FY2/Final/psuade/NEE/LH-180/RU-FY2/output/rst/rstfac.txt',&
!                      form='formatted',status='new',action='write')     
      
!      print*,'numpatch:',numpatch
      CALL rstTimeConstRead (lat_points,lhistTimeConst,nfcon,nftune,&
                             numpatch,numpatch_lat,ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                             fcon,ftune)
      
      CALL rstTimeVarRead (lhistTimeVar,nfvar,numpatch,idate,fvar)
      
! ======================================================================
! begin time stepping loop
! ======================================================================

      fldxy(:,:,:)=0.
      nac = 0

      do istep = 1, mstep
         idate_p(:) = idate(:)
!         print*,'CLM.F90-test:istep:',istep
       ! Read in the atmospheric forcing
         CALL GETMET (lumet,lon_points,lat_points,nforc,forcxy,idate(2),idate(3))
        
       ! Mapping atmospheric fields to force clm: [lon_points]x[lat_points] grid
       !     -> [numpatch] vector of subgrid points
         do k = 1, numpatch   !clm vector index
            i = ixy_patch(k)  !longitude index
            j = jxy_patch(k)  !latitude index
            forc(k,1:) = forcxy(i,j,1:) 
         end do

       ! doalb is true when the next time step is a radiation time step
         doalb = .true. 
         dolai = .true.
         dosst = .false.
         oro(:) = 1.

#if(!defined EcoDynamics)
         if(dolai)then
       ! read(lulai) 
         do j = 1, lat_points
            do i = 1, lon_points
               read(lulai) mflai(i,j,1:maxpatch,1:nflai)
            enddo
         enddo

         do k = 1, numpatch
            i = ixy_patch(k)
            j = jxy_patch(k)
            m = mxy_patch(k)
            flai(k,:) = mflai(i,j,m,:) !leaf area index
         enddo
         endif
#endif
       ! Calendar for NEXT time step
         CALL TICKTIME (deltim,idate)

       ! Call clm driver
         CALL CLMDRIVER (nl_soil,maxsnl,numpatch,idate,deltim,&
                         nftune,nfcon,nforc,nfvar,nfldv,&
                         ftune,fcon,forc,fvar,fldv,rstfac,&
#if(!defined EcoDynamics)
                         flai,&
#endif
                         dolai,doalb,dosst,oro)
!        print*,'CLM.F90-test-rstfac',rstfac
!        rstfac_all(istep,1)=rstfac(1,1)
!        rstfac_all(istep,2)=rstfac(2,1)
!        print*,'CLM.F90-test-rstfac_all',rstfac_all(istep,:)
       ! Mapping subgrid patch [numpatch] vector of subgrid points to 
       !     -> [lon_points]x[lat_points] grid average
         itypwat(:) = nint(fcon(:,3))
         CALL vec2xy(lat_points,lon_points,numpatch,&
                     ixy_patch,jxy_patch,wtxy_patch,itypwat,nfcon,nforc,nfldv,&
                     fcon,forcxy,fldv,fldxy_r)

         nac = nac + 1
         fldxy(:,:,:) = fldxy(:,:,:) + fldxy_r(:,:,:)

       ! Logical idenfication for writing output and open file
!         lwrite = .false.
!         CALL lpwrite(idate_p,idate,lhistTimeVar,lhistTimeVar2,&
!                       luout,luout2,foutdat,lwrite)

       ! Write out the model variables for restart run [histTimeVar] and the histroy file
!      my test for do not output the file to save running time
       lwrite=.true.        
        if(lwrite)then

!            CALL rstTimeVarWrite (lhistTimeVar,nfvar,numpatch,idate,fvar)
            
!            CALL rstTimeVarWrite2(lhistTimeVar2,nfvar,numpatch,idate,fvar)

            a(:,:,:) = fldxy(:,:,:) / float(nac)
!            CALL flxwrite (luout,lon_points,lat_points,nfldv,a)
            
!            CALL flxwrite2(luout2,lon_points,lat_points,nfldv,a)
!output interest var to my SA output matrix:eg H--3
            SA_LE(istep)=a(1,1,4)
            SA_NEE(istep)=a(1,1,22)*1000000-a(1,1,21)*1000000
!            SA_H(istep)=a(1,1,62)
            SA_SWC(istep)=a(1,1,36)

            fldxy(:,:,:) = 0.0
            lwrite = .false.
            nac = 0

         endif

      end do
!      print*,'begin output my interest variance'
!      do istep=1,mstep
!         write(lmyvariance,101) rstfac_all(istep,:)
!101      format(1x,2f10.7)
!      end do
!     print*,'end output my interest variance'
!     print*,'begin output my SA output matrix'
!open the SA output file
!     print*,'CLM.F90-test-SA_LE',SA_LE(1)
     open(unit=lsao,file=fsao,form='formatted',&
               status='unknown',action='write')
!write the output matrix to certain file:output_step.txt
     CALL Myoutputwrite_vector(lsao,mstep,SA_LE)

!     print*,'CLM.F90-TEST-SA_NEE',SA_NEE(1)    
     open(unit=lsao2,file=fsao2,form='formatted',&
               status='unknown',action='write')
     CALL Myoutputwrite_vector(lsao2,mstep,SA_NEE)

!     print*,'CLM.F90-TEST-SA_H',SA_H(1)
!     open(unit=lsao3,file=fsao3,form='formatted',&
!               status='unknown',action='write')
!     CALL Myoutputwrite_vector(lsao3,mstep,SA_H)
!     CALL Myoutputwrite_matrix(lsao3,mstep,3,SA_H)

!     print*,'CLM.F90-TEST-SA-SWC',SA_SWC(1)
     open(unit=lsao4,file=fsao4,form='formatted',&
               status='unknown',action='write')
     CALL Myoutputwrite_vector(lsao4,mstep,SA_SWC)
!      write(6,*) 'CLM Execution Completed'

  END PROGRAM CLM
! ----------------------------------------------------------------------
! EOP
#endif



#if(defined coup_atmosmodel)

  SUBROUTINE CLM(numpatch,lon_points,lat_points,& 
                 nl_soil,maxsnl,nftune,nfcon,nforc,&
                 nfvar,nfldv,maxpatch,& 
                 ixy_patch,jxy_patch,mxy_patch,wtxy_patch,&
                 idate,deltim,ftune,fcon,fvar,&
#if(!defined EcoDynamics)
                 mflai,&
#endif
                 dolai,doalb,dosst,oro,forcxy,fldxy_r)

! ======================================================================
! The Common Land Model was developed in cooperation with
!     Beijing Normal University and IRI-Columbia University (Dai)
!     Georgia Institute of Technology                       (Dickinson)
!     National Center for Atmospheric Research              (Bonan, Oleson)
!     University of Arizona                                 (Zeng)
!     University of Texas at Austin                         (Yang)
!     GSFC/NASA                                             (Houser, Bosilovich)
!     COLA                                                  (Dirmeyer, Schlosser)
!     Colorado State University at Fort Collins             (Denning, Baker)
!
! Reference: 
!     [1] Dai et al., 2003: The Common Land Model (CLM). 
!         Bull. of Amer. Meter. Soc., 84: 1013-1023
!     [2] Dai et al., 2004: A two-big-leaf model for canopy temperature,
!         photosynthesis and stomatal conductance. Journal of Climate
!
!     Author: Yongjiu Dai, January 2004
! ======================================================================

      use precision
      implicit none

! ----------------local variables ---------------------------------

      integer, INTENT(in) :: numpatch   ! total number of patches of grids
      integer, INTENT(in) :: lon_points ! number of longitude points on model grid
      integer, INTENT(in) :: lat_points ! number of latitude points on model grid
      integer, INTENT(in) :: nl_soil    ! number of soil layers
      integer, INTENT(in) :: maxsnl     ! max number of snow layers
      integer, INTENT(in) :: nftune     ! number of clm tunable constants
      integer, INTENT(in) :: nfcon      ! number of time constant variables
      integer, INTENT(in) :: nforc      ! number of forcing variables
      integer, INTENT(in) :: nfvar      ! number of time varying variables
      integer, INTENT(in) :: nfldv      ! number of output variables
      integer, INTENT(in) :: maxpatch   ! max number of patches in a grid

      real   , INTENT(in) :: deltim     ! time step (senconds)
                                                 
#if(!defined EcoDynamics) 
      real(r8), INTENT(in) :: mflai(lon_points,lat_points,maxpatch,4)
      real(r8) flai(numpatch,4)         ! time varying vegetation parameters
#endif
      real(r8), INTENT(in) :: forcxy(lon_points,lat_points,nforc)   ! atmospheric forcing
      real(r8), INTENT(out) :: fldxy_r(lon_points,lat_points,nfldv) ! output fluxes in 2-dimension
                                                     
      integer,  INTENT(in) :: idate(3)               ! calendar (year, julian day, seconds)
      integer,  INTENT(in) :: ixy_patch(numpatch)    ! patch longitude index
      integer,  INTENT(in) :: jxy_patch(numpatch)    ! patch latitude index
      integer,  INTENT(in) :: mxy_patch(numpatch)    ! patch subgrid index of lnd point
      real(r8), INTENT(in) :: wtxy_patch(numpatch)   ! patch weight
      real(r8), INTENT(in) :: fcon(numpatch,nfcon)   ! time constant variables
      real(r8), INTENT(inout) :: fvar(numpatch,nfvar)! time varying variables
      real(r8), INTENT(in) :: ftune(nftune)          ! clm tunable constants
      real(r8), INTENT(inout) :: oro(numpatch)          ! ocean(0)/seaice(2)/ flag
                                               
      logical, INTENT(in) :: doalb  ! true if time for surface albedo calculation
      logical, INTENT(in) :: dolai  ! true if time for time-varying vegetation paramter
      logical, INTENT(in) :: dosst  ! true if time for update sst/ice/snow
                                              
      integer :: i,j,k,l,m          ! looping indices
      real(r8) forc(numpatch,nforc) ! forcing variables
      real(r8) fldv(numpatch,nfldv) ! output fluxes
      integer itypwat(numpatch)     ! land water type

! ----------------------------------------------------------------------
! begin the land surface model
! ----------------------------------------------------------------------

       ! mapping atmospheric fields to force clm: [lon_points]x[lat_points] grid
       !     -> [numpatch] vector of subgrid points
         do k = 1, numpatch   !clm vector index
            i = ixy_patch(k)  !longitude index
            j = jxy_patch(k)  !latitude index
            forc(k,1:) = forcxy(i,j,1:) 
         end do

#if(!defined EcoDynamics)
         if(dolai)then
         do k = 1, numpatch
            i = ixy_patch(k)
            j = jxy_patch(k)
            m = mxy_patch(k)
            flai(k,1) = mflai(i,j,m,1) !leaf area index
            flai(k,2) = mflai(i,j,m,2) !stem are index
            flai(k,3) = mflai(i,j,m,3) !fraction of green leaf
            flai(k,4) = mflai(i,j,m,4) !fraction of veg. cover
         enddo
         endif
#endif

       ! call clm driver
         CALL CLMDRIVER (nl_soil,maxsnl,numpatch,idate,deltim,&
                         nftune,nfcon,nforc,nfvar,nfldv,&
                         ftune,fcon,forc,fvar,fldv,&
#if(!defined EcoDynamics)
                         flai,&
#endif
                         dolai,doalb,dosst,oro)

       ! Mapping subgrid patch [numpatch] vector of subgrid points to 
       !     -> [lon_points]x[lat_points] grid average
         itypwat(:) = nint(fcon(:,3))
         CALL vec2xy(lat_points,lon_points,numpatch,&
                     ixy_patch,jxy_patch,wtxy_patch,itypwat,nfcon,nforc,nfldv,&
                     fcon,forcxy,fldv,fldxy_r)

  END SUBROUTINE CLM
! ----------------------------------------------------------------------
! EOP
#endif
