!
!----------------------------------------------------------
! Compute vertical advective fluxes 
! using 5th-order WENO scheme
!----------------------------------------------------------
!
#  ifdef EW_PERIODIC
          imin=1
          imax=LOCALLM+1
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
#   else
          imin=3
          imax=Lm-1
#   endif
#  endif
!
!----------------------------------------------------------------------
!  k loop: FC
!----------------------------------------------------------------------
!

!$acc loop independent
        do k=3,N-3
          do i=IstrU,Iend
            if ( i.ge.imin .and. i.le.imax ) then
              FC(i,k)=flux6(We(i-3,j,k),We(i-2,j,k),We(i-1,j,k),
     &                  We(i  ,j,k),We(i+1,j,k),We(i+2,j,k),1.)
            else
              FC(i,k)=0.5*(We(i-1,j,k)+We(i,j,k))
            endif
            FC(i,k)=FC(i,k)*FLUX5(
     &           u(i,j,k-2,nrhs), u(i,j,k-1,nrhs), 
     &           u(i,j,k  ,nrhs), u(i,j,k+1,nrhs),
     &           u(i,j,k+2,nrhs), u(i,j,k+3,nrhs), FC(i,k))
          enddo
        enddo

        do i=IstrU,Iend
          if ( i.ge.imin .and. i.le.imax ) then
            FC(i,2)=flux6(We(i-3,j,2),We(i-2,j,2),We(i-1,j,2),
     &                We(i  ,j,2),We(i+1,j,2),We(i+2,j,2),1.)
          else
            FC(i,2)=0.5*(We(i-1,j,2)+We(i,j,2))
          endif
          FC(i,2)=FC(i,2)*FLUX3(
     &         u(i,j,1,nrhs), u(i,j,2,nrhs), 
     &         u(i,j,3,nrhs), u(i,j,4,nrhs), FC(i,2))

          if ( i.ge.imin .and. i.le.imax ) then
            FC(i,N-2)=flux6(We(i-3,j,N-2),We(i-2,j,N-2),We(i-1,j,N-2),
     &                We(i  ,j,N-2),We(i+1,j,N-2),We(i+2,j,N-2),1.)
          else   
            FC(i,N-2)=0.5*(We(i-1,j,N-2)+We(i,j,N-2))
          endif
          FC(i,N-2)=FC(i,N-2)*FLUX3(
     &         u(i,j,N-3,nrhs), u(i,j,N-2,nrhs), 
     &         u(i,j,N-1,nrhs), u(i,j,N  ,nrhs), FC(i,N-2))

          if ( i.ge.imin .and. i.le.imax ) then
            FC(i,1)=flux6(We(i-3,j,1),We(i-2,j,1),We(i-1,j,1),
     &                We(i  ,j,1),We(i+1,j,1),We(i+2,j,1),1.)
          else
            FC(i,1)=0.5*(We(i-1,j,1)+We(i,j,1))
          endif
          FC(i,1)=FC(i,1)*FLUX2(
     &         u(i,j,1,nrhs), u(i,j,2,nrhs), FC(i,1), cdif)

          if ( i.ge.imin .and. i.le.imax ) then
            FC(i,N-1)=flux6(We(i-3,j,N-1),We(i-2,j,N-1),We(i-1,j,N-1),
     &                We(i  ,j,N-1),We(i+1,j,N-1),We(i+2,j,N-1),1.)
          else
            FC(i,N-1)=0.5*(We(i-1,j,N-1)+We(i,j,N-1))
          endif
          FC(i,N-1)=FC(i,N-1)*FLUX2(
     &         u(i,j,N-1,nrhs), u(i,j,N,nrhs), FC(i,N-1), cdif)
	    
          FC(i,0)=0.
          FC(i,N)=0.
        enddo
