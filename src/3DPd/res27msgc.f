      subroutine res27msgc(acor,q,p,r,nallx2,nally2, 
     &                     ixs,iys,nx,ny,nz, 
     &                     neigz,neigx,neigy)
c
      implicit real*8 (a-h,o-z)
      complex*16  acor(0:nz+1,nallx2,nally2,27) 
      complex*16  q(0:nz+1,nallx2,nally2),p(0:nz+1,nallx2,nally2) 
      complex*16  r(0:nz+1,nallx2,nally2) 
c
      do j = iys+1,iys+ny
         do i = ixs+1,ixs+nx
            do k = 1,nz
               r(k,i,j) = p(k,i,j) 
     &                  - acor(k,i,j,1 )*q(k-1,i-1,j-1)
     &                  - acor(k,i,j,2 )*q(k  ,i-1,j-1)
     &                  - acor(k,i,j,3 )*q(k+1,i-1,j-1)
     &                  - acor(k,i,j,4 )*q(k-1,i  ,j-1)
     &                  - acor(k,i,j,5 )*q(k  ,i  ,j-1)
     &                  - acor(k,i,j,6 )*q(k+1,i  ,j-1)
     &                  - acor(k,i,j,7 )*q(k-1,i+1,j-1)
     &                  - acor(k,i,j,8 )*q(k  ,i+1,j-1)
     &                  - acor(k,i,j,9 )*q(k+1,i+1,j-1)
     &                  - acor(k,i,j,10)*q(k-1,i-1,j  )
     &                  - acor(k,i,j,11)*q(k  ,i-1,j  )
     &                  - acor(k,i,j,12)*q(k+1,i-1,j  )
     &                  - acor(k,i,j,13)*q(k-1,i  ,j  )
     &                  - acor(k,i,j,14)*q(k  ,i  ,j  )
     &                  - acor(k,i,j,15)*q(k+1,i  ,j  )
     &                  - acor(k,i,j,16)*q(k-1,i+1,j  )
     &                  - acor(k,i,j,17)*q(k  ,i+1,j  )
     &                  - acor(k,i,j,18)*q(k+1,i+1,j  )
     &                  - acor(k,i,j,19)*q(k-1,i-1,j+1)
     &                  - acor(k,i,j,20)*q(k  ,i-1,j+1)
     &                  - acor(k,i,j,21)*q(k+1,i-1,j+1)
     &                  - acor(k,i,j,22)*q(k-1,i  ,j+1)
     &                  - acor(k,i,j,23)*q(k  ,i  ,j+1)
     &                  - acor(k,i,j,24)*q(k+1,i  ,j+1)
     &                  - acor(k,i,j,25)*q(k-1,i+1,j+1)
     &                  - acor(k,i,j,26)*q(k  ,i+1,j+1)
     &                  - acor(k,i,j,27)*q(k+1,i+1,j+1)
            enddo
         enddo
      enddo
c
      return
      end
