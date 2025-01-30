!
!===============================================================
!
! Compute 5th order horizontal advection
!
!===============================================================
!
#  ifdef NS_PERIODIC
          jmin=1
          jmax=LOCALMM+1
#  else
#   ifdef MPI
          if (SOUTH_INTER) then
            jmin=1
          else
            jmin=3
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-1
          endif
#   else
          jmin=3
          jmax=Mm-1
#   endif
#  endif
#  ifdef EW_PERIODIC
          imin=0
          imax=LOCALLM
#  else
#   ifdef MPI
          if (WEST_INTER) then
            imin=0
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi
          else
            imax=Lmmpi-2
          endif
#   else
          imin=3
          imax=Lm-2
#   endif
#  endif
!
!----------------------------------------------------------------------
!  j loop: UFe
!----------------------------------------------------------------------
!

!$acc loop independent
          DO j = max(Jstr,jmin),min(Jend+1,jmax)  !j_loop_y_flux_5
                                                  ! use full stencil
                                                  !
              DO i = IstrU,Iend
                if ( i.ge.imin+1 .and. i.le.imax ) then
           UFe(i,j) = flux6(Hvom(i-3,j,k),Hvom(i-2,j,k),Hvom(i-1,j,k),
     &                        Hvom(i  ,j,k),Hvom(i+1,j,k),Hvom(i+2,j,k),1.)
                else
           UFe(i,j) = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                endif
                flx5 = UFe(i,j)*FLUX5(
     &             u(i,j-3,k,nrhs), u(i,j-2,k,nrhs), 
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs),
     &             u(i,j+1,k,nrhs), u(i,j+2,k,nrhs), UFe(i,j) )
#  ifdef MASKING
!                UFe(i,j) = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                      Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                flx3 = UFe(i,j)*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), UFe(i,j) ) 
!                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                flx2 = UFe(i,j)*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), UFe(i,j), cdif)
#   ifdef UP5_MASKING
                mask2=umask(i,j-2)*umask(i,j+1)
                IF (UFe(i,j).gt.0) THEN
                  mask1=umask(i,j-2)
                  mask0=umask(i,j-3)*mask2          
                ELSE
                  mask1=umask(i,j+1)
                  mask0=umask(i,j+2)*mask2
                ENDIF
#   else
                mask1=umask(i,j-2)*umask(i,j+1)
                mask2=umask(i,j-3)*umask(i,j+2)
                mask0=mask1*mask2
#   endif
                UFe(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2

#  else
                UFe(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO ! j_loop_y_flux_5

!$acc loop independent
          DO j = Jstr,Jend+1  !j_loop_y_flux_5
              DO i = IstrU,Iend
                IF ( j.eq.jmin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                UFe(i,j) = UFe(i,j)*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), UFe(i,j), cdif)
                                                             !
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
               UFe(i,j) = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
!                UFe(i,j) = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                      Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                flx3 = UFe(i,j)*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), UFe(i,j) )
#  ifdef MASKING
!                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
                flx2 = UFe(i,j)*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), UFe(i,j), cdif)
                mask1=umask(i,j-2)*umask(i,j+1)
                UFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFe(i,j)=flx3
#  endif

            ELSE IF ( j.eq.jmax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+Hvom(i,j,k))
                UFe(i,j) = UFe(i,j)*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), UFe(i,j), cdif)
                                          !
            ELSE IF ( j.eq.jmax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
!               UFe(i,j) = flux4(Hvom(i-2,j,k),Hvom(i-1,j,k),
!     &                     Hvom(i  ,j,k),Hvom(i+1,j,k),1.)
                flx3 = UFe(i,j)*FLUX3(
     &             u(i,j-2,k,nrhs), u(i,j-1,k,nrhs),
     &             u(i,j  ,k,nrhs), u(i,j+1,k,nrhs), UFe(i,j) )
#  ifdef MASKING
!                UFe(i,j) = 0.5*(Hvom(i-1,j,k)+ Hvom(i,j,k))
                flx2 = UFe(i,j)*FLUX2(
     &             u(i,j-1,k,nrhs), u(i,j  ,k,nrhs), UFe(i,j), cdif)
                mask1=umask(i,j-2)*umask(i,j+1)
                UFe(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFe(i,j)=flx3
#  endif
              ENDIF
            ENDDO
          ENDDO ! j_loop_y_flux_5
!
!----------------------------------------------------------------------
!  i loop: UFx
!----------------------------------------------------------------------
!
          DO j = Jstr,Jend
!$acc loop independent
            DO i = max(IstrU-1,imin),min(Iend,imax)  !i_loop_x_flux_5
                                                  ! use full stencil
              
            UFx(i,j) = flux6(Huon(i-2,j,k),Huon(i-1,j,k),Huon(i  ,j,k),
     &                   Huon(i+1,j,k),Huon(i+2,j,k),Huon(i+3,j,k),1.)
!           UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx5 = UFx(i,j)*FLUX5(
     &             u(i-2,j,k,nrhs), u(i-1,j,k,nrhs),
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs),
     &             u(i+2,j,k,nrhs), u(i+3,j,k,nrhs), UFx(i,j) )
#  ifdef MASKING
            UFx(i,j) = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)
                flx3 = UFx(i,j)*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), UFx(i,j) )       
            UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = UFx(i,j)*FLUX2(
     &        u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), UFx(i,j), cdif)
#   ifdef UP5_MASKING
                mask2=umask(i-1,j)*umask(i+2,j)
                IF (UFx(i,j).gt.0) THEN
                  mask1=umask(i-1,j)
                  mask0=umask(i-2,j)*mask2          
                ELSE
                  mask1=umask(i+2,j)
                  mask0=umask(i+3,j)*mask2
                ENDIF
#   else
                mask1=umask(i-1,j)*umask(i+2,j)
                mask2=umask(i-2,j)*umask(i+3,j)
                mask0=mask1*mask2
#   endif
                UFx(i,j)=mask0*flx5+(1-mask0)*mask1*flx3+
     &                              (1-mask0)*(1-mask1)*flx2
#  else
                UFx(i,j)=flx5
#  endif /* MASKING */
            ENDDO
          ENDDO ! i_loop_x_flux_5


          DO j = Jstr,Jend
!$acc loop independent
            DO i = IstrU-1,Iend  !i_loop_x_flux_5                                              !
            IF ( i.eq.imin-2 ) THEN   ! 2nd order flux next to south
                                           ! boundary
                UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx(i,j) = UFx(i,j)*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), UFx(i,j), cdif)
                                                             !
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN  ! 3rd of 4th order flux 2 in
                                                             ! from south boundary
!               UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx(i,j) = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)    
                flx3 = UFx(i,j)*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs), UFx(i,j) )
#  ifdef MASKING
                UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = UFx(i,j)*FLUX2(
     &             u(i  ,j,k,nrhs), u(i+1,j,k,nrhs), UFx(i,j), cdif)
                mask1=umask(i-1,j)*umask(i+2,j)
                UFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFx(i,j)=flx3
#  endif
                                          !
            ELSE IF ( i.eq.imax+2 ) THEN  ! 2nd order flux next to north
                                          ! boundary
                UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx(i,j) = UFx(i,j)*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), UFx(i,j), cdif)
                                          !
            ELSE IF ( i.eq.imax+1 ) THEN  ! 3rd or 4th order flux 2 in from
                                          ! north boundary
!               UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                UFx(i,j) = flux4(Huon(i-1,j,k),Huon(i  ,j,k),
     &                      Huon(i+1,j,k),Huon(i+2,j,k),1.)  
                flx3 = UFx(i,j)*FLUX3(
     &             u(i-1,j,k,nrhs), u(i  ,j,k,nrhs),
     &             u(i+1,j,k,nrhs), u(i+2,j,k,nrhs),  UFx(i,j) )
#  ifdef MASKING
                UFx(i,j) = 0.5*(Huon(i,j,k)+Huon(i+1,j,k))
                flx2 = UFx(i,j)*FLUX2(
     &             u(i,j,k,nrhs), u(i+1,j,k,nrhs), UFx(i,j), cdif)
                mask1=umask(i-1,j)*umask(i+2,j)
                UFx(i,j)=mask1*flx3+(1-mask1)*flx2
#  else
                UFx(i,j)=flx3
#  endif
            ENDIF
          ENDDO
        ENDDO ! i_loop_x_flux_5

