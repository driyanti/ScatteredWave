      subroutine extup3a(w,work,witp,nallx2,nally2, 
     &                   nx,nxf,ixs,ixsf,ix0f, 
     &                   ny,nyf,iys,iysf,iy0f, 
     &                   nz, 
     &                   neigx,neigy)
c
      implicit real*8 (a-h,o-z)
      complex*16  w(0:nz+1,nallx2,nally2) 
      complex*16  work(0:nz+1,nallx2,nally2) 
      complex*16  witp(0:nz+1,nallx2,nally2,4) 
      integer     neigx(2), neigy(2) 
c
      if (nxf.le.0 .or. nyf.le.0 ) return
      nx1 = nxf
      ny1 = nyf
c
      do iy = 1, ny
         j = iys + iy
         jf = iysf + iy0f + iy*2 - 2
         do ix = 1, nx
            i = ixs + ix
            if = ixsf + ix0f + ix*2 - 2
            do k = 1, nz
               work(k,if,jf) = w(k,i,j)
               w(k,if,jf) = w(k,if,jf) + work(k,if,jf)
            enddo
         enddo
      enddo
c
      ix1f = mod(ix0f,2)+1
      iy1f = mod(iy0f,2)+1
c
      do jf = iysf + iy0f, iysf + nyf, 2
         do if = ixsf + ix1f, ixsf + nxf, 2
            do k = 1, nz
               w(k,if,jf) = w(k,if,jf) 
     &                    + witp(k,if,jf,1)*work(k,if-1,jf)
     &                    + witp(k,if,jf,2)*work(k,if+1,jf)
            enddo
         enddo
      enddo
c
      do jf = iysf + iy1f, iysf + nyf, 2
         do if = ixsf + ix0f, ixsf + nxf, 2
            do k = 1, nz
               w(k,if,jf) = w(k,if,jf) 
     &                    + witp(k,if,jf,1)*work(k,if,jf-1)
     &                    + witp(k,if,jf,2)*work(k,if,jf+1)
            enddo
         enddo
      enddo
c
      do jf = iysf + iy1f, iysf + nyf, 2
         do if = ixsf + ix1f, ixsf + nxf, 2
            do k = 1, nz
               w(k,if,jf) = w(k,if,jf) 
     &                    + witp(k,if,jf,1)*work(k,if-1,jf-1)
     &                    + witp(k,if,jf,2)*work(k,if+1,jf-1)
     &                    + witp(k,if,jf,3)*work(k,if-1,jf+1)
     &                    + witp(k,if,jf,4)*work(k,if+1,jf+1)
            enddo
         enddo
      enddo
c
      return
      end
