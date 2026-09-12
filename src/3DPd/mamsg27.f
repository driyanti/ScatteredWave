      subroutine mamsg27(acor,w0,w1,nallx2,nally2, 
     &                   ixs,iys,
     &                   nx,ny, 
     &                   nz,
     &                   neigz,neigx,neigy)
c
      implicit    real*8 (a-h,o-z)
      complex*16  acor(0:nz+1,nallx2,nally2,27) 
      complex*16  w0(0:nz+1,nallx2,nally2),w1(0:nz+1,nallx2,nally2) 
      integer     neigz(2),neigx(2),neigy(2) 
c
      do j = iys+1,iys+ny
         do i = ixs+1,ixs+nx
            do k = 1, nz
               w1(k, i,j) = acor(k,i,j,1 )*w0(k-1,i-1,j-1)
     &                    + acor(k,i,j,2 )*w0(k  ,i-1,j-1)
     &                    + acor(k,i,j,3 )*w0(k+1,i-1,j-1)
     &                    + acor(k,i,j,4 )*w0(k-1,i  ,j-1)
     &                    + acor(k,i,j,5 )*w0(k  ,i  ,j-1)
     &                    + acor(k,i,j,6 )*w0(k+1,i  ,j-1)
     &                    + acor(k,i,j,7 )*w0(k-1,i+1,j-1)
     &                    + acor(k,i,j,8 )*w0(k  ,i+1,j-1)
     &                    + acor(k,i,j,9 )*w0(k+1,i+1,j-1)
     &                    + acor(k,i,j,10)*w0(k-1,i-1,j  )
     &                    + acor(k,i,j,11)*w0(k  ,i-1,j  )
     &                    + acor(k,i,j,12)*w0(k+1,i-1,j  )
     &                    + acor(k,i,j,13)*w0(k-1,i  ,j  )
     &                    + acor(k,i,j,14)*w0(k  ,i  ,j  )
     &                    + acor(k,i,j,15)*w0(k+1,i  ,j  )
     &                    + acor(k,i,j,16)*w0(k-1,i+1,j  )
     &                    + acor(k,i,j,17)*w0(k  ,i+1,j  )
     &                    + acor(k,i,j,18)*w0(k+1,i+1,j  )
     &                    + acor(k,i,j,19)*w0(k-1,i-1,j+1)
     &                    + acor(k,i,j,20)*w0(k  ,i-1,j+1)
     &                    + acor(k,i,j,21)*w0(k+1,i-1,j+1)
     &                    + acor(k,i,j,22)*w0(k-1,i  ,j+1)
     &                    + acor(k,i,j,23)*w0(k  ,i  ,j+1)
     &                    + acor(k,i,j,24)*w0(k+1,i  ,j+1)
     &                    + acor(k,i,j,25)*w0(k-1,i+1,j+1)
     &                    + acor(k,i,j,26)*w0(k  ,i+1,j+1)
     &                    + acor(k,i,j,27)*w0(k+1,i+1,j+1)
            enddo
         enddo
      enddo
c
      return
      end
