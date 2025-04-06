
#include <define.h>
 
  subroutine lpwrite(idate_p,idate,lhistTimeVar,lhistTimeVar2,luout,luout2,foutdat,lwrite)

  implicit none
  integer, intent(in) :: idate_p(3)
  integer, intent(in) :: idate(3)
  integer, intent(in) :: lhistTimeVar
  integer, intent(in) :: lhistTimeVar2
  integer, intent(in) :: luout
  integer, intent(in) :: luout2
  character(LEN=256), intent(in) :: foutdat 
  logical, intent(inout) :: lwrite
  integer :: months(12)
  integer :: j_month_e
  logical :: leapyear
  integer i, l
  character(LEN=256) fhistTimeVar_name
  character(LEN=256) fhistTimeVar_name2
  character(LEN=256) fout
  character(LEN=256) fout2
  character(LEN=256) cdate
 
         
#if(defined WR_DAILY)
  if(idate_p(2).ne.idate(2))then
     lwrite = .true.
  endif
#elif(defined WR_MONTHLY)
  leapyear = (mod(idate_p(1),4)==0.and.mod(idate_p(1),100)/=0).or.&
                                       mod(idate_p(1),400)==0
  if(leapyear)then
     months = (/31,60,91,121,152,182,213,244,274,305,335,366/)
  else
     months = (/31,59,90,120,151,181,212,243,273,304,334,365/)
  endif
  j_month_e = -999
  do l = 1, 12
     if(idate_p(2).eq.months(l))then
        j_month_e = 1
        exit
     endif
  enddo
  if((j_month_e.gt.0) .and. (idate_p(2).ne.idate(2)))then
     lwrite = .true.
  endif
#elif(defined WR_HOURLY)
     lwrite = .true.
#endif

! Open for model time varying data (model state variables) and history filed
  if(lwrite)then
     write(cdate,'(i4.4,"-",i3.3,"-",i5.5)') idate_p(1),idate_p(2),idate_p(3)
     fhistTimeVar_name = trim(foutdat)//'-rstTimeVar'//'-'//trim(cdate)
     fhistTimeVar_name2= '/home/u120220909911/RU-FY2/Third/psuade/NEE/LH/RU-FY2/output/rst/Valdai'//&
                           '-rstTimeVar'//'-'//trim(cdate)//'.txt'
     fout=trim(foutdat)//'-'//trim(cdate)
     fout2=trim(foutdat)//'-'//trim(cdate)//'.txt'

     open(unit=lhistTimeVar,file=fhistTimeVar_name,form='unformatted',&
                            status='new',action='write')
     open(unit=lhistTimeVar2,file=fhistTimeVar_name2,form='formatted',&
                            status='new',action='write')
     open(unit=luout,file=fout,access='sequential',form='unformatted',&
                     status='new',action='write')
   
     open(unit=luout2,file=fout2,access='sequential',form='formatted',&
                     status='new',action='write')
  endif

  end subroutine lpwrite
! ------------------------------------------------------------------------
! EOP
