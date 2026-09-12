      subroutine makep(a,ngx,ngy,ngz,omg,idata)
     
      implicit none
      integer     ix,iy,iz,k,idata 
      integer     ngx,ngy,ngz 
      complex*16  a(ngx,ngy,ngz,7) 
      complex*16  cidd,cone,czero 
      real*8      hx,hy,hz,hx2,hy2,hz2,omg,wdata,wnum  
      parameter   (cidd  = (0.d0,1.d0), cone = (1.d0,0.d0),
     &             czero = (0.d0,0.d0))
      
      do k = 1,7
         do iz = 1,ngz
            do iy = 1,ngy
               do ix = 1,ngx
                  a(ix,iy,iz,k) = (0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      
      hx = 1.d0/dble(float(ngx+1))
      hy = 1.d0/dble(float(ngy+1))  
      hz = 1.d0/dble(float(ngz+1))  
      hx2 = 1.d0/hx**2
      hy2 = 1.d0/hy**2
      hz2 = 1.d0/hz**2 
c 
      do iz = 1,ngz
         do iy = 1,ngy
            do ix = 1,ngx
               wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
               a(ix,iy,iz,1) = cmplx(-hz2,0.d0)
               a(ix,iy,iz,2) = cmplx(-hy2,0.d0)
               a(ix,iy,iz,3) = cmplx(-hx2,0.d0)
               a(ix,iy,iz,4) = 2.d0*(hx2+hy2+hz2)-(1.d0,-0.5d0)*wnum**2
               a(ix,iy,iz,5) = cmplx(-hx2,0.d0)
               a(ix,iy,iz,6) = cmplx(-hy2,0.d0)
               a(ix,iy,iz,7) = cmplx(-hz2,0.d0)
            enddo
         enddo
      enddo
c
c     Boundary conditions ::
c
c     N : North
      do iz = 1,ngz 
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,ngy+1,iz)
            a(ix,ngy,iz,4) = a(ix,ngy,iz,4) - 1.d0/((1.d0,0.d0) 
     &                     + 1.d0*cidd*wnum*hy)/hy**2.d0
            a(ix,ngy,iz,6) = (0.d0,0.d0)
         enddo
      enddo
c     S : South
      do iz = 1,ngz
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,0,iz)
            a(ix,1,iz,4) = a(ix,1,iz,4) - 1.d0/((1.d0,0.d0) 
     &                     + 1.d0*cidd*wnum*hy)/hy**2.d0
            a(ix,1,iz,2) = (0.d0,0.d0)
         enddo
      enddo
c     E : East
      do iz = 1,ngz
         do iy = 1,ngy
            wnum = wdata(omg,ngx,ngy,ngz,idata,ngx+1,iy,iz)
            a(ngx,iy,iz,4) = a(ngx,iy,iz,4) - 1.d0/((1.d0,0.d0)  
     &                     + 1.d0*cidd*wnum*hx)/hx**2.d0
            a(ngx,iy,iz,5) = (0.d0,0.d0)
         enddo
      enddo
c     W : West
      do iz = 1,ngz
         do iy = 1,ngy
            wnum = wdata(omg,ngx,ngy,ngz,idata,0,iy,iz)
            a(1,iy,iz,4) = a(1,iy,iz,4) - 1.d0/((1.d0,0.d0)  
     &                   + 1.d0*cidd*wnum*hx)/hx**2.d0
            a(1,iy,iz,3) = (0.d0,0.d0)
         enddo
      enddo
c     R : Right
      do iy = 1,ngy
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,ngz+1)
            a(ix,iy,ngz,4) = a(ix,iy,ngz,4) - 1.d0/((1.d0,0.d0)  
     &                     + 1.d0*cidd*wnum*hz)/hz**2.d0
            a(ix,iy,ngz,7) = (0.d0,0.d0)
         enddo
      enddo
c     L : Left 
      do iy = 1,ngy
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,0)
            a(ix,iy,1,4) = a(ix,iy,1,4) - 1.d0/((1.d0,0.d0)  
     &                   + 1.d0*cidd*wnum*hz)/hz**2.d0
            a(ix,iy,1,1) = (0.d0,0.d0)
         enddo
      enddo     
c
      return
      end
      
