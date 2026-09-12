      subroutine mkacor(a,acor,ngx,ngy,nz,nallx2,nally2,ixsf,iysf,
     &                   idata,omg)
c
      implicit none
      integer     nz,ngx,ngy,nallx2,nally2,ngx2,ngy2,nz2 
      integer     k,ix,iy,if,jf,kf,ixsf,iysf,idata 
      real*8      hx,hy,hz,hx2,hy2,hz2,hh,hh2,omg,wdata,wnum  
      complex*16  a(nz,ngx,ngy,7) 
      complex*16  acor(0:nz+1,nallx2,nally2,27) 
      complex*16  cidd 
      parameter   (cidd = (0.d0,1.d0))
c
      do iy = 1,ngy
         do ix = 1,ngx
            do k = 1,nz
               acor(k,ix+ixsf,iy+iysf, 5) = a(k,ix,iy,1)
               acor(k,ix+ixsf,iy+iysf,11) = a(k,ix,iy,2)
               acor(k,ix+ixsf,iy+iysf,13) = a(k,ix,iy,3)
               acor(k,ix+ixsf,iy+iysf,14) = a(k,ix,iy,4)
               acor(k,ix+ixsf,iy+iysf,15) = a(k,ix,iy,5)
               acor(k,ix+ixsf,iy+iysf,17) = a(k,ix,iy,6)
               acor(k,ix+ixsf,iy+iysf,23) = a(k,ix,iy,7)
            enddo
         enddo
      enddo
c
      open(10,file='matrixa.dat')
      do iy = 1,ngy
         do ix = 1,ngx
            do k = 1,nz
               write(10,*) iy,ix,k
               write(10,*) acor(k,ix+ixsf,iy+iysf, 5)
               write(10,*) acor(k,ix+ixsf,iy+iysf,11)
               write(10,*) acor(k,ix+ixsf,iy+iysf,13)
               write(10,*) acor(k,ix+ixsf,iy+iysf,14)
               write(10,*) acor(k,ix+ixsf,iy+iysf,15)
               write(10,*) acor(k,ix+ixsf,iy+iysf,17)
               write(10,*) acor(k,ix+ixsf,iy+iysf,23)
               write(10,*)
            enddo
         enddo
      enddo 
      close(10)
c     
      return
      end
      
