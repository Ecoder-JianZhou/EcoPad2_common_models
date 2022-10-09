
program test
    character(len=100) parsFile
    character(len=100) forcingFile
    real, dimension (:,:), allocatable :: forcingData
    integer :: nLines 
    integer, parameter :: ilines=20
    integer, parameter :: iiterms=10
    real forcingData_raw(iiterms,ilines)
    parsFile = "SPRUCE_pars.txt"
    forcingFile = "SPRUCE_forcing.txt"
    write(*,*) "test1"
    call Getclimate(forcingFile, ilines, iiterms, nLines, forcingData_raw)
    allocate (forcingData(iiterms,nLines))
    forcingData = forcingData_raw(:,1:nLines)
    write(*,*)forcingData 
end 




subroutine Getclimate(forcingFile, ilines, iiterms, nLines, forcingData_raw)
    ! arg: fileName, row, column, tot_line, return_data
    implicit none
    integer ilines, iiterms
    integer nLines
    real forcingData_raw(iiterms,ilines)
    character(len=100) forcingFile,commts
    integer nRow, mCol, istat1

    open(11,file=forcingFile,status='old',ACTION='read',IOSTAT=istat1)
    ! skip 2 lines of input met data file
    read(11,'(a160)') commts
    nRow = 0  ! to record the lines in a file
    do while(.true.)   ! read forcing files
        nRow = nRow + 1
        read(11,*,IOSTAT=istat1)(forcingData_raw(mCol, nRow),mCol=1,iiterms)
        if(istat1<0)exit
    enddo ! end of reading the forcing file
    nLines = nRow - 1
    close(11)    ! close forcing file
    return
end