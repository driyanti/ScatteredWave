      subroutine makeah2(a,ngx,ngy,ngz,omg,idata)
c
c     This subroutine returns the preconditioning matrix M with
c     second order BC of Engquist and Majda discretized using second
c     order finite difference scheme.
c
c     Remarks: Virtual points required. NGX,NGY,NGZ are number of grid 
c              points in the INTERIOR.
c              "The sign" may not be correct!
c     
      implicit none
      integer     ix,iy,iz,k,idata 
      integer     ngx,ngy,ngz 
      complex*16  a(0:ngx+1,0:ngy+1,0:ngz+1,7) 
      complex*16  cidd,cone,czero 
      real*8      hx,hy,hz,hx2,hy2,hz2,hh,hh2,omg,wdata,wnum  
      parameter   (cidd  = (0.d0,1.d0), cone = (1.d0,0.d0),
     &             czero = (0.d0,0.d0))
      
      write(*,*) ngx,ngy,ngz 
      do k = 1,7
         do iz = 0,ngz+1
            do iy = 0,ngy+1
               do ix = 0,ngx+1
                  a(ix,iy,iz,k) = (0.d0,0.d0)
               enddo
            enddo
         enddo
      enddo
      
      hx  =  1.d0/dble(float(ngx+1))
      hy  =  1.d0/dble(float(ngy+1))  
      hz  =  1.d0/dble(float(ngz+1))  
      hx2 =  1.d0/hx**2
      hy2 =  1.d0/hy**2
      hz2 =  1.d0/hz**2
      hh  =  hx
      hh2 =  hx2
c 
      do iz = 0,ngz+1
         do iy = 0,ngy+1
            do ix = 0,ngx+1
               wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
               a(ix,iy,iz,1) = -cmplx(hz2,0.d0)
               a(ix,iy,iz,2) = -cmplx(hy2,0.d0)
               a(ix,iy,iz,3) = -cmplx(hx2,0.d0)
               a(ix,iy,iz,4) =  2.d0*(hx2+hy2+hz2)-wnum**2
               a(ix,iy,iz,5) = -cmplx(hx2,0.d0)
               a(ix,iy,iz,6) = -cmplx(hy2,0.d0)
               a(ix,iy,iz,7) = -cmplx(hz2,0.d0)
            enddo
         enddo
      enddo
c
c     Boundary conditions ::
c
c     FACES::
c     N : North
      iy = ngy+1
      do iz = 1,ngz 
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = a(ix,iy,iz,1) + cidd/(wnum*hz)*hz2
            a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
            a(ix,iy,iz,3) = a(ix,iy,iz,3) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2
            a(ix,iy,iz,5) = a(ix,iy,iz,5) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,6) = (0.d0,0.d0)
            a(ix,iy,iz,7) = a(ix,iy,iz,7) + cidd/(wnum*hz)*hz2
         enddo
      enddo
c     S : South
      iy =  0
      do iz = 1,ngz
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = a(ix,iy,iz,1) + cidd/(wnum*hz)*hz2
            a(ix,iy,iz,2) = (0.d0,0.d0)
            a(ix,iy,iz,3) = a(ix,iy,iz,3) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2
            a(ix,iy,iz,5) = a(ix,iy,iz,5) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
            a(ix,iy,iz,7) = a(ix,iy,iz,7) + cidd/(wnum*hz)*hz2
         enddo
      enddo
c     E : East
      ix = ngx+1
      do iz = 1,ngz
         do iy = 1,ngy
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = a(ix,iy,iz,1) + cidd/(wnum*hz)*hz2
            a(ix,iy,iz,2) = a(ix,iy,iz,2) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2     
            a(ix,iy,iz,5) = (0.d0,0.d0)
            a(ix,iy,iz,6) = a(ix,iy,iz,6) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,7) = a(ix,iy,iz,7) + cidd/(wnum*hz)*hz2
         enddo
      enddo
c     W : West
      ix = 0
      do iz = 1,ngz
         do iy = 1,ngy
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = a(ix,iy,iz,1) + cidd/(wnum*hz)*hz2
            a(ix,iy,iz,2) = a(ix,iy,iz,2) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,3) = (0.d0,0.d0)
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
            a(ix,iy,iz,6) = a(ix,iy,iz,6) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,7) = a(ix,iy,iz,7) + cidd/(wnum*hz)*hz2
         enddo
      enddo
c     R : Right
      iz = ngz+1
      do iy = 1,ngy
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
            a(ix,iy,iz,2) = a(ix,iy,iz,2) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,3) = a(ix,iy,iz,3) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hz)*((wnum*hz)**2-2.d0)*hz2
            a(ix,iy,iz,5) = a(ix,iy,iz,5) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,6) = a(ix,iy,iz,6) + cidd/(wnum*hy)*hz2
            a(ix,iy,iz,7) = (0.d0,0.d0)
         enddo
      enddo
c     L : Left 
      iz = 0
      do iy = 1,ngy
         do ix = 1,ngx
            wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
            a(ix,iy,iz,1) = (0.d0,0.d0)
            a(ix,iy,iz,2) = a(ix,iy,iz,2) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,3) = a(ix,iy,iz,3) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,4) = a(ix,iy,iz,4)
     &                    + 2.d0*cidd/(wnum*hz)*((wnum*hz)**2-2.d0)*hz2
            a(ix,iy,iz,5) = a(ix,iy,iz,5) + cidd/(wnum*hx)*hx2
            a(ix,iy,iz,6) = a(ix,iy,iz,6) + cidd/(wnum*hy)*hy2
            a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
         enddo
      enddo     
c
c     EDGES::
c     NORTH FACE:
c     (1,NGY+1,0) --> (NGX,NGY+1,0)
      iy = ngy+1
      iz = 0
      do ix = 1,ngx
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = (0.d0,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,6) = (0.d0,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
      enddo
c
c     (1,NGY+1,NGZ+1) --> (NGX,NGY+1,NGZ+1)
      iy = ngy+1
      iz = ngz+1
      do ix = 1,ngx
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,6) = (0.d0,0.d0)
         a(ix,iy,iz,7) = (0.d0,0.d0)
      enddo
c
c     (0,NGY+1,1) --> (0,NGY+1,NGZ)
      ix = 0
      iy = ngy+1
      do iz = 1,ngz
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) + 1.d0*cidd/(wnum*hz)*hz2
         a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
         a(ix,iy,iz,3) = (0.d0,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
         a(ix,iy,iz,6) = (0.d0,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) + 1.d0*cidd/(wnum*hz)*hz2
      enddo
c
c     (NGX+1,NGY+1,1) --> (NGX+1,NGY+1,NGZ)
      ix = ngx+1
      iy = ngy+1
      do iz = 1,ngz
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) + 1.d0*cidd/(wnum*hz)*hz2
         a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a(ix,iy,iz,5) = (0.d0,0.d0)
         a(ix,iy,iz,6) = (0.d0,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) + 1.d0*cidd/(wnum*hz)*hz2
      enddo
c
c     SOUTH FACE:
c     (1,0,0) --> (NGX,0,0)
      iy = 0
      iz = 0
      do ix = 1,ngx
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = (0.d0,0.d0)
         a(ix,iy,iz,2) = (0.d0,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
      enddo
c
c     (1,0,NGZ+1) --> (NGX,0,NGZ+1)
      iy = 0
      iz = ngz+1
      do ix = 1,ngx
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
         a(ix,iy,iz,2) = (0.d0,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) + 1.d0*cidd/(wnum*hx)*hx2
         a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
         a(ix,iy,iz,7) = (0.d0,0.d0)
      enddo
c
c     (0,0,1) --> (0,0,NGZ)
      ix = 0
      iy = 0
      do iz = 1,ngz
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) + 1.d0*cidd/(wnum*hz)*hz2
         a(ix,iy,iz,2) = (0.d0,0.d0)
         a(ix,iy,iz,3) = (0.d0,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
         a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) + 1.d0*cidd/(wnum*hz)*hz2
      enddo
c
c     (NGX+1,0,1) --> (NGX+1,0,NGZ)
      ix = ngx+1
      iy = 0
      do iz = 1,ngz
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) + 1.d0*cidd/(wnum*hz)*hz2
         a(ix,iy,iz,2) = (0.d0,0.d0)
         a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a(ix,iy,iz,5) = (0.d0,0.d0)
         a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
         a(ix,iy,iz,7) = a(ix,iy,iz,7) + 1.d0*cidd/(wnum*hz)*hz2
      enddo

c     LEFT FACE:
c     P(0,1,0) --> (0,NGY,0)
      ix = 0
      iz = 0
      do iy = 1,ngy
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = (0.d0,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,3) = (0.d0,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
         a(ix,iy,iz,6) = a(ix,iy,iz,6) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
      enddo
c
c     P(NGX+1,1,0) --> (NGX+1,NGY,0)
      ix = ngx+1
      iz = 0
      do iy = 1,ngy
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = (0.d0,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &                - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a(ix,iy,iz,5) = (0.d0,0.d0)
         a(ix,iy,iz,6) = a(ix,iy,iz,6) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
      enddo
c
c     RIGHT FACE:
c     P(0,1,NGZ+1) --> (0,NGY,NGZ+1)
      ix = 0
      iz = ngz+1
      do iy = 1,ngy
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,3) = (0.d0,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
         a(ix,iy,iz,6) = a(ix,iy,iz,6) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,7) = (0.d0,0.d0)
      enddo
c
c     P(NGX+1,1,NGZ+1) --> (NGX+1,NGY,NGZ+1)
      ix = ngx+1
      iz = ngz+1
      do iy = 1,ngy
         wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
         a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
         a(ix,iy,iz,2) = a(ix,iy,iz,2) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
         a(ix,iy,iz,4) = a(ix,iy,iz,4) 
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a(ix,iy,iz,5) = 0.d0
         a(ix,iy,iz,6) = a(ix,iy,iz,6) + 1.d0*cidd/(wnum*hy)*hy2
         a(ix,iy,iz,7) = 0.d0
      enddo

c     CORNER:
c     P(0,0,0)
      ix = 0
      iy = 0
      iz = 0
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = 0.d0
      a(ix,iy,iz,2) = 0.d0
      a(ix,iy,iz,3) = 0.d0
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
      a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
      a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
c
c     P(NGX+1,0,0)
      ix = ngx+1
      iy = 0
      iz = 0
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = 0.d0
      a(ix,iy,iz,2) = 0.d0
      a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = 0.d0
      a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
      a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
c
c     P(0,NGY+1,0)
      ix = 0
      iy = ngy+1
      iz = 0
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = 0.d0
      a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
      a(ix,iy,iz,3) = 0.d0
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
      a(ix,iy,iz,6) = 0.d0
      a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
c
c     P(NGX+1,NGY+1,0)
      ix = ngx+1
      iy = ngy+1
      iz = 0
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = 0.d0
      a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
      a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = 0.d0
      a(ix,iy,iz,6) = 0.d0
      a(ix,iy,iz,7) = a(ix,iy,iz,7) - cmplx(hz2,0.d0)
c
c     P(0,0,NGZ+1)
      ix = 0
      iy = 0
      iz = ngz+1
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
      a(ix,iy,iz,2) = 0.d0
      a(ix,iy,iz,3) = 0.d0
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
      a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
      a(ix,iy,iz,7) = 0.d0
c
c     P(NGX+1,0,NGZ+1)
      ix = ngx+1
      iy = 0
      iz = ngz+1
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
      a(ix,iy,iz,2) = 0.d0
      a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = 0.d0
      a(ix,iy,iz,6) = a(ix,iy,iz,6) - cmplx(hy2,0.d0)
      a(ix,iy,iz,7) = 0.d0
c
c     P(0,NGY+1,NGZ+1)
      ix = 0
      iy = ngy+1
      iz = ngz+1
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
      a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
      a(ix,iy,iz,3) = 0.d0
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = a(ix,iy,iz,5) - cmplx(hx2,0.d0)
      a(ix,iy,iz,6) = 0.d0
      a(ix,iy,iz,7) = 0.d0
c
c     P(NGX+1,NGY+1,NGZ+1)
      ix = ngx+1
      iy = ngy+1
      iz = ngz+1
      wnum = wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
      a(ix,iy,iz,1) = a(ix,iy,iz,1) - cmplx(hz2,0.d0)
      a(ix,iy,iz,2) = a(ix,iy,iz,2) - cmplx(hy2,0.d0)
      a(ix,iy,iz,3) = a(ix,iy,iz,3) - cmplx(hx2,0.d0)
      a(ix,iy,iz,4) = a(ix,iy,iz,4) + 4.d0*cidd*wnum*hh*hh2
      a(ix,iy,iz,5) = 0.d0
      a(ix,iy,iz,6) = 0.d0
      a(ix,iy,iz,7) = 0.d0
c
      return
      end
      
