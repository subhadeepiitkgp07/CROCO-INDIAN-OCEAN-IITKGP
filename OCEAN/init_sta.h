!======================================================================
! CROCO is a branch of ROMS developped at IRD, INRIA, 
! Ifremer, CNRS and Univ. Toulouse III  in France
! The two other branches from UCLA (Shchepetkin et al)
! and Rutgers University (Arango et al) are under MIT/X style license.
! CROCO specific routines (nesting) are under CeCILL-C license.
!
! CROCO website : http://www.croco-ocean.org
!======================================================================
!
! Include file "init_sta.h".
! ==============================
!
! statitle         Stations application title.
! stat0,stax0, ... input start time and positions from sta.in
! stacoor          type of coordinates in input(lat,lon or x,y)a
! STgrd            input grid level station in nested applications

      real stat0(Msta), stax0(Msta), stay0(Msta), staz0(Msta)
      common /ncrealsta/ stat0, stax0, stay0, staz0

      integer  stacoor(Msta), STgrd(Msta)
      common /ncintsta/ stacoor,STgrd

! to be tested if STgrd is really necessary ALVARO

      character*80 statitle
      common /nccharsta/ statitle
