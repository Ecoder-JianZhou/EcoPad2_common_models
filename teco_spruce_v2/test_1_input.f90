
program test
    character(len=100) parsFile = "SPRUCE_pars.txt"
    character(len=100) forcingFile = "SPRUCE_forcing.txt"
    real, dimension (:,:), forcingData :: darray
    integer :: nLines 
    integer, parameter :: ilines=150000
    integer, parameter :: iiterms=10
    real forcing_data_raw(iiterms,ilines)
    call Getclimate(forcingFile, ilines, iiterms, nLines, forcingData)
end 



! subroutine Getclimate(year_seq,doy_seq,hour_seq,          &
!     &   forcing_data,climatefile,lines,yr_length)
subroutine Getclimate(forcingFile, ilines, iiterms, nLines, forcingData)
    implicit none
    integer, parameter :: ilines
    integer, parameter :: iiterms
    integer nLines
    ! integer,dimension(ilines):: year_seq,doy_seq,hour_seq
    real forcingData(iiterms,ilines)
    character(len=150) forcingFile,commts
    integer m,n,istat1,lines,yr_length

    open(11,file=forcingFile,status='old',ACTION='read',     &
    &     IOSTAT=istat1)
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    m=0  ! to record the lines in a file
    yr_length=0 ! to record years of a dataset
    do    ! read forcing files
        m=m+1
        read(11,*,IOSTAT=istat1)year_seq(m),      &
        &       doy_seq(m),hour_seq(m),           &
        &       (forcing_data(n,m),n=1,iiterms)
        if(istat1<0)exit
    enddo ! end of reading the forcing file
    lines=m-1
    yr_length=(year_seq(lines)-year_seq(1))+1
    close(11)    ! close forcing file
    return
end