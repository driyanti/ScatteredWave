      subroutine mkacorf(acor,ngx,ngy,nz,nallx2,nally2,ixsf,iysf,
     &                   idata,omg,bet1,bet2)
c
      implicit none
      integer     nz,ngx,ngy,nallx2,nally2,ngx2,ngy2,nz2 
      integer     k,ix,iy,ixo,iyo,izo,ixsf,iysf,idata 
      real*8      hx,hy,hz,hx2,hy2,hz2,hh,hh2,omg,wdata,wnum  
      real*8      bet1,bet2 
c      complex*16  a(nz,ngx,ngy,7) 
      complex*16  acor(0:nz+1,nallx2,nally2,27) 
      complex*16  cidd 
      parameter   (cidd = (0.d0,1.d0))
c
c      do iy = 1,ngy
c         do ix = 1,ngx
c            do k = 1,nz
c               acor(k,ix+ixsf,iy+iysf, 5) = a(k,ix,iy,1)
c               acor(k,ix+ixsf,iy+iysf,11) = a(k,ix,iy,2)
c               acor(k,ix+ixsf,iy+iysf,13) = a(k,ix,iy,3)
c               acor(k,ix+ixsf,iy+iysf,14) = a(k,ix,iy,4)
c               acor(k,ix+ixsf,iy+iysf,15) = a(k,ix,iy,5)
c               acor(k,ix+ixsf,iy+iysf,17) = a(k,ix,iy,6)
c               acor(k,ix+ixsf,iy+iysf,23) = a(k,ix,iy,7)
c            enddo
c         enddo
c      enddo
c
c      write(*,*) ngx,ngy,nz,omg,idata,ixsf,iysf 
c      do iy = 1,ngy
c         do ix = 1,ngx
c            do k = 1,nz
c               acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,14) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
c               acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
c            enddo
c         enddo
c      enddo
c
      ngx2 = ngx-2
      ngy2 = ngy-2
      nz2  = nz -2
c      write(*,*) ngx2,ngy2,nz2 
c      
      hx  =  1.d0/dble(float(ngx2+1))
      hy  =  1.d0/dble(float(ngy2+1))  
      hz  =  1.d0/dble(float(nz2+1))  
      hx2 =  1.d0/hx**2
      hy2 =  1.d0/hy**2
      hz2 =  1.d0/hz**2
      hh  =  hx
      hh2 =  hx2
c      write(*,*) hx,hy,hz 
c 
      do iy = 2,ngy-1
         iyo = iy-1
         do ix = 2,ngx-1
            ixo = ix-1
            do k = 2,nz-1
               izo = k-1
               wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
               acor(k,ix+ixsf,iy+iysf, 5) = -cmplx(hy2,0.d0)
               acor(k,ix+ixsf,iy+iysf,11) = -cmplx(hx2,0.d0)
               acor(k,ix+ixsf,iy+iysf,13) = -cmplx(hz2,0.d0)
               acor(k,ix+ixsf,iy+iysf,14) = 2.d0*(hx2+hy2+hz2)
     &                                    - cmplx(bet1,-bet2)*wnum**2
               acor(k,ix+ixsf,iy+iysf,15) = -cmplx(hz2,0.d0)
               acor(k,ix+ixsf,iy+iysf,17) = -cmplx(hx2,0.d0)
               acor(k,ix+ixsf,iy+iysf,23) = -cmplx(hy2,0.d0)
            enddo
         enddo
      enddo
c
c     Boundary conditions ::
c
c     FACES::
c     N : North
      iy = ngy
      iyo = iy-1
      do ix = 2,ngx-1
         ixo = ix-1
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                   - cmplx(hy2,0.d0)
            acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,14) =   
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &                  + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2  
            acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)     
         enddo
      enddo
c     S : South
      iy = 1
      iyo = iy-1
      do ix = 2,ngx-1
         ixo = ix-1
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
            acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0) 
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,14) =
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &                  + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2
            acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0) 
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                   - cmplx(hy2,0.d0)
         enddo
      enddo
c     E : East
      ix = ngx
      ixo = ix-1
      do iy = 2,ngy-1
         iyo = iy-1
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
            acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                   - cmplx(hx2,0.d0)
            acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,14) = 
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2    
     &                  + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0) 
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
            acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
        enddo
      enddo
c     W : West
      ix = 1
      ixo = ix-1
      do iy = 2,ngy-1
         iyo = iy-1
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
            acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
            acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,14) = 
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2 
     &                  + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0) 
     &                                   + cidd/(wnum*hz)*hz2
            acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                                   - cmplx(hx2,0.d0)
            acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
         enddo
      enddo
c     R : Right
      k  = nz
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         do ix = 2,ngx-1
            ixo = ix-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
            acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0) 
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                   - cmplx(hz2,0.d0)
            acor(k,ix+ixsf,iy+iysf,14) =
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2 
     &                  + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
            acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                                 + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                 + cidd/(wnum*hy)*hy2
         enddo
      enddo
c     L : Left 
      k  = 1
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         do ix = 2,ngx-1
            ixo = ix-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                   + cidd/(wnum*hy)*hy2
            acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0) 
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
            acor(k,ix+ixsf,iy+iysf,14) =
     &                    2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &                  + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                   - cmplx(hz2,0.d0)
            acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                                   + cidd/(wnum*hx)*hx2
            acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                 + cidd/(wnum*hy)*hy2
         enddo
      enddo     
c
c     EDGES::
c     NORTH FACE:
c     (2,NGY,1) --> (NGX-1,NGY,1)
      iy = ngy
      iyo = iy-1
      k  = 1
      izo = k-1
      do ix = 2,ngx-1
         ixo = ix-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
      enddo
c
c     (2,NGY,NZ) --> (NGX-1,NGY,NZ)
      iy = ngy
      iyo = iy-1
      k  = nz
      izo = k-1
      do ix = 2,ngx-1
         ixo = ix-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0) 
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
      enddo
c
c     (1,NGY,2) --> (1,NGY,NZ-1)
      ix = 1
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)         
      enddo
c
c     (NGX,NGY,2) --> (NGX,NGY,NZ-1)
      ix = ngx
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
      enddo
c
c     SOUTH FACE:
c     (2,1,1) --> (NGX-1,1,1)
      iy = 1
      iyo = iy-1
      k  = 1
      izo = k-1
      do ix = 2,ngx-1
         ixo = ix-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                              + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) =
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
      enddo
c
c     (2,1,NZ) --> (NGX-1,1,NZ)
      iy = 1
      iyo = iy-1
      k  = nz
      izo = k-1
      do ix = 2,ngx-1
         ixo = ix-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                + 1.d0*cidd/(wnum*hx)*hx2
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
      enddo
c
c     (1,1,2) --> (1,1,NZ-1)
      ix = 1
      ixo = ix-1
      iy = 1
      iyo = iy-1
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,14) =
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
      enddo
c
c     (NGX,1,2) --> (NGX,1,NZ-1)
      ix = ngx
      ixo = ix-1
      iy = 1
      iyo = iy-1
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                + 1.d0*cidd/(wnum*hz)*hz2
         acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                - cmplx(hy2,0.d0)
      enddo
c
c     LEFT FACE:
c     P(1,2,1) --> (1,NGY-1,1)
      ix = 1
      ixo = ix-1
      k  = 1
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
         acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
      enddo
c
c     P(NGX,2,1) --> (NGX,NGY-1,1)
      ix = ngx
      ixo = ix-1
      k  = 1
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
      enddo
c
c     RIGHT FACE:
c     P(1,2,NZ) --> (1,NGY-1,NZ)
      ix = 1
      ixo = ix-1
      k  = nz
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
         acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
      enddo
c
c     P(NGX,2,NZ) --> (NGX,NGY-1,NZ)
      ix = ngx
      ixo = ix-1
      k  = nz
      izo = k-1
      do iy = 2,ngy-1
         iyo = iy-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
         acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                                - cmplx(hx2,0.d0)
         acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                                - cmplx(hz2,0.d0)
         acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
         acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                                + 1.d0*cidd/(wnum*hy)*hy2
      enddo
c
c     CORNER:
c     P(1,1,1)
      ix = 1
      ixo = ix-1
      iy = 1
      iyo = iy-1
      k  = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
c      write(*,*) ixo,iyo,izo,wnum
      acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                             - cmplx(hx2,0.d0) 
      acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
c
c     P(NGX,1,1)
      ix  = ngx
      ixo = ix-1
      iy  = 1
      iyo = iy-1
      k   = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
c      write(*,*) wnum
      acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                             - cmplx(hx2,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = !wnum
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2.d0
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)     
c
c     P(1,NGY,1)
      ix = 1
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      k  = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                             - cmplx(hx2,0.d0) 
      acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
c      
c     P(NGX,NGY,1)
      ix = ngx
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      k  = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0) 
     &                             - cmplx(hx2,0.d0) 
      acor(k,ix+ixsf,iy+iysf,13) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
c
c     P(1,1,NZ)
      ix = 1
      ixo = ix-1
      iy = 1
      iyo = iy-1
      k  = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0) 
     &                             - cmplx(hx2,0.d0) 
      acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
c
c     P(NGX,1,NZ)
      ix = ngx
      ixo = ix-1
      iy = 1
      iyo = iy-1
      k  = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                             - cmplx(hx2,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,23) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
c
c     P(1,NGY,NZ)
      ix = 1
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      k  = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = - cmplx(hx2,0.d0)
     &                             - cmplx(hx2,0.d0)
      acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
c
c     P(NGX,NGY,NZ)
      ix = ngx
      ixo = ix-1
      iy = ngy
      iyo = iy-1
      k  = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      acor(k,ix+ixsf,iy+iysf, 5) = - cmplx(hy2,0.d0)
     &                             - cmplx(hy2,0.d0)
      acor(k,ix+ixsf,iy+iysf,11) = - cmplx(hx2,0.d0)
     &                             - cmplx(hx2,0.d0)
      acor(k,ix+ixsf,iy+iysf,13) = - cmplx(hz2,0.d0)
     &                             - cmplx(hz2,0.d0)
      acor(k,ix+ixsf,iy+iysf,14) = 
     &                 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &               + 4.d0*cidd*wnum*hh*hh2
      acor(k,ix+ixsf,iy+iysf,15) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,17) = (0.d0,0.d0)
      acor(k,ix+ixsf,iy+iysf,23) = (0.d0,0.d0)
c
c      open(10,file='matrixaf.dat')
c      do iy = 1,ngy
c         do ix = 1,ngx
c            do k = 1,nz
c               write(10,*) iy,ix,k
c               write(10,*) acor(k,ix+ixsf,iy+iysf, 5)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,11)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,13)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,14)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,15)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,17)
c               write(10,*) acor(k,ix+ixsf,iy+iysf,23)
c               write(10,*)
c            enddo
c         enddo
c      enddo 
c      close(10)
c     
      return
      end
      
