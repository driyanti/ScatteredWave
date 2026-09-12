      subroutine ext3a(w0,w1,witp,nallx2,nally2, 
     &                 nx,nxf,ixs,ixsf,ix0, 
     &                 ny,nyf,iys,iysf,iy0, 
     &                 nz, 
     &                 neigx,neigy) 
c
c       w1 is the extension of w0 along diagonal.
c       w1 <- exten_DA(w0)
c
      implicit real*8 (a-h,o-z)
      complex*16  w0(0:nz+1,nallx2,nally2),w1(0:nz+1,nallx2,nally2) 
      complex*16  witp(0:nz+1,nallx2,nally2,4) 
c
      if (nxf.le.0 .or. nyf.le.0 ) return
      nx1 = nxf
      ny1 = nyf
c
      ix1 = mod(ix0,2)+1
      iy1 = iy0
      ix2 = ix0
      iy2 = mod(iy0,2)+1
      ix3 = mod(ix0,2)+1
      iy3 = mod(iy0,2)+1
c
      do iy = 1,ny
         j = iys + iy
         jf = iysf + iy0 + iy*2 - 2
         do ix = 1,nx
            i = ixs + ix
            if = ixsf + ix0 + ix*2 - 2
            do k = 1,nz
               w1(k,if,jf) = w0(k,i,j)
            enddo
         enddo
      enddo
c
c     compute on grid points 1
c
      do j = iysf+iy1,iysf+nyf,2
         do i = ixsf+ix1,ixsf+nxf,2
            do k = 1,nz
               w1(k,i,j) = witp(k,i,j,1)*w1(k,i-1,j)
     &                   + witp(k,i,j,2)*w1(k,i+1,j)
            enddo
         enddo
      enddo
c
c     compute on grid points 2
c
      do j = iysf+iy2,iysf+nyf,2
         do i = ixsf+ix2,ixsf+nxf,2
            do k = 1,nz
               w1(k,i,j) = witp(k,i,j,1)*w1(k,i,j-1)
     &                   + witp(k,i,j,2)*w1(k,i,j+1)
            enddo
         enddo
      enddo
c
c     compute on grid points 3
c
      do j = iysf+iy3,iysf+nyf,2
         do i = ixsf+ix3,ixsf+nxf,2
            do k = 1,nz
               w1(k,i,j) = witp(k,i,j,1)*w1(k,i-1,j-1)
     &                   + witp(k,i,j,2)*w1(k,i+1,j-1)
     &                   + witp(k,i,j,3)*w1(k,i-1,j+1)
     &                   + witp(k,i,j,4)*w1(k,i+1,j+1)
            enddo
         enddo
      enddo
c
      return
      end
