! This version is used to run forecasting in SPRUCE or other site.  -- Jian Zhou Oct. 05 2022
! 
! ==============================================================================================================
! this version, Int_CWE is a combined version of improved soil thermal (from Yuanyuan) and methane (from Shuang) 
! that adapted to SPRUCE site, a boreal forest wetland ecosystem
! The Int_CWE is adapted to be ready for ECOPAD webpage run 
!                                                   -- Shuang Ma Mar. 28 2018
! NOTES: edit on March 28, convert format from JJ Github(webpage ready, but no switches) to a full version of CWE with switches
! carbon and soil water function well before switching on the soil thermal module
! ice(i) = ice actual thickness of ice unit m
! liq_water(i) = water actual thickness unit m
! wcl(i) = water + ice actual thickness unit m
! *********************************************************
! STATEMENT OF CHANGE IN THIS NEW VERSION 1.2 (v1.2)
! 1. CHANGED CH4 UNIT
! 2. changed /wsc(i) to /thicksl(i)
! 3. equation for diffusion flux difu(1)
!
! need to change together:      1) do i=1,5114         !number of days of simulated and projected years
!                               2) forcing data
! two switch needed in ebullition part to exclude/include ebullition into chamber flux

!                if(yr.gt.yrs_eq)then  !!!!!!!!!!!!!   commented for spin up
!            stopped write out file 82 83 in spin up
! ===============================================================================================================


program TECO_main
    write(*,*) "This is TECO_main ..."
    call a1()
    call a2()
    call a3()
    call a4()
    call a5()
    call a6()
    call a7()
    call a8()
    call a9()
    call a10()
    call a11()
    call a12()
    call a13()
    call a14()
    call a15()
    call a16()
    call a17()
    call a18()
    call a19()
    call a20()
    call a21()
    call a22()
    call a23()
    call a24()
    call a25()
    call a26()
    call a27()
    call a28()
    call a29()
    call a30()
    call a31()
end

subroutine a1()
    write(*,*) "this is a1 ..."
end subroutine a1

subroutine a2()
    write(*,*) "this is a2 ..."
end subroutine a2

subroutine a3()
    write(*,*) "this is a3 ..."
end subroutine a3

subroutine a4()
    write(*,*) "this is a4 ..."
end subroutine a4

subroutine a5()
    write(*,*) "this is a5 ..."
end subroutine a5

subroutine a6()
    write(*,*) "this is a6 ..."
end subroutine a6

subroutine a7()
    write(*,*) "this is a7 ..."
end subroutine a7

subroutine a8()
    write(*,*) "this is a8 ..."
end subroutine a8

subroutine a9()
    write(*,*) "this is a9 ..."
end subroutine a9

subroutine a10()
    write(*,*) "this is a10 ..."
end subroutine a10

subroutine a11()
    write(*,*) "this is a11 ..."
end subroutine a11

subroutine a12()
    write(*,*) "this is a12 ..."
end subroutine a12

subroutine a13()
    write(*,*) "this is a13 ..."
end subroutine a13

subroutine a14()
    write(*,*) "this is a14 ..."
end subroutine a14

subroutine a15()
    write(*,*) "this is a15 ..."
end subroutine a15

subroutine a16()
    write(*,*) "this is a16 ..."
end subroutine a16

subroutine a17()
    write(*,*) "this is a17 ..."
end subroutine a17

subroutine a18()
    write(*,*) "this is a18 ..."
end subroutine a18

subroutine a19()
    write(*,*) "this is a19 ..."
end subroutine a19

subroutine a20()
    write(*,*) "this is a20 ..."
end subroutine a20

subroutine a21()
    write(*,*) "this is a21 ..."
end subroutine a21

subroutine a22()
    write(*,*) "this is a22 ..."
end subroutine a22

subroutine a23()
    write(*,*) "this is a23 ..."
end subroutine a23

subroutine a24()
    write(*,*) "this is a24 ..."
end subroutine a24

subroutine a25()
    write(*,*) "this is a25 ..."
end subroutine a25

subroutine a26()
    write(*,*) "this is a26 ..."
end subroutine a26

subroutine a27()
    write(*,*) "this is a27 ..."
end subroutine a27

subroutine a28()
    write(*,*) "this is a28 ..."
end subroutine a28

subroutine a29()
    write(*,*) "this is a29 ..."
end subroutine a29

subroutine a30()
    write(*,*) "this is a30 ..."
end subroutine a30

subroutine a31()
    write(*,*) "this is a31 ..."
end subroutine a31
