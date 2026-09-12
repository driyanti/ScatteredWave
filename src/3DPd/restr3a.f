      subroutine restr3a(r,f,wrst,nallx2,nally2, 
     &                   nx,nxf,ixs,ixsf,ix0f, 
     &                   ny,nyf,iys,iysf,iy0f, 
     &                   nz, 
     &                   neigx, neigy)
c
c     restriction of r to f
c
      implicit    real*8 (a-h,o-z)
      complex*16  r(0:nz+1,nallx2,nally2),f(0:nz+1,nallx2,nally2) 
      complex*16  wrst(0:nz+1,nallx2,nally2,-1:1,-1:1) 
c
      if (nxf.le.0 .or. nyf.le.0 ) return
      nx1 = nxf
      ny1 = nyf
c
      do iy = 1,ny
         j = iy + iys
         jf = iysf+iy0f+2*iy-2
         do ix = 1,nx
            i = ix + ixs
            if = ixsf+ix0f+2*ix-2
            do k = 1,nz
               f(k,i,j) = wrst(k,i,j,-1,-1)*r(k,if-1,jf-1)
     &                  + wrst(k,i,j, 0,-1)*r(k,if,jf-1)
     &                  + wrst(k,i,j, 1,-1)*r(k,if+1,jf-1)
     &                  + wrst(k,i,j,-1, 0)*r(k,if-1,jf)
     &                  + wrst(k,i,j, 0, 0)*r(k,if,jf)
     &                  + wrst(k,i,j, 1, 0)*r(k,if+1,jf)
     &                  + wrst(k,i,j,-1, 1)*r(k,if-1,jf+1)
     &                  + wrst(k,i,j, 0, 1)*r(k,if,jf+1)
     &                  + wrst(k,i,j, 1, 1)*r(k,if+1,jf+1)
            enddo
         enddo
      enddo
c
      return
      end
