      subroutine mat3d7pf(ngx,ngy,ngz,omg,idata,x,y,work)
c
      implicit none
      integer     ngx,ngy,ngz,nx,ny,nz 
      integer     ix,iy,iz,if,jf,kf,idata 
      complex*16  x(ngx,ngy,ngz),y(ngx,ngy,ngz) 
      complex*16  work(0:ngx+1,0:ngy+1,0:ngz+1) 
      real*8      hx,hy,hz,hh,hx2,hy2,hz2,hh2 
      real*8      a1,a2,a3,a4,a5,a6,a7 
      real*8      omg,wnum,wdata 
      complex*16  cidd 
      parameter   (cidd = (0.d0,1.d0))

      nx  = ngx - 2
      ny  = ngy - 2
      nz  = ngz - 2
      hx  = 1.d0/dble(float(nx+1))
      hy  = 1.d0/dble(float(ny+1))
      hz  = 1.d0/dble(float(nz+1))
      hx2 = 1.d0/hx**2
      hy2 = 1.d0/hy**2
      hz2 = 1.d0/hz**2
      hh  = hx
      hh2 = hx2
c
      a1 = -hz2
      a2 = -hy2
      a3 = -hx2
      a4 = 2.d0*(hx2+hy2+hz2)
      a5 = -hx2
      a6 = -hy2
      a7 = -hz2
      
c      write(*,*) ngx,ngy,ngz
c      write(*,*) nx,ny,nz 
c      write(*,*)'cek_mat',ngx,ngy,ngz,hx,hy,hz
c
c     Interior points (1:NGX,1:NGY,1:NGZ)
c
      do iz = 0,ngz+1
         do iy = 0,ngy+1
            do ix = 0,ngx+1
               work(ix,iy,iz) = 0.d0
            enddo
         enddo
      enddo
      do iz = 1,ngz 
         do iy = 1,ngy 
            do ix = 1,ngx 
               work(ix,iy,iz) = x(ix,iy,iz)
            enddo
         enddo
      enddo
c 
c      do iz=1,ngz
c       do iy=1,ngy
c        do ix=1,ngx
c	 write(22,*)iz,iy,ix,x(ix,iy,iz)	 
c	end do
c       end do
c      end do  	
    
      do iz = 2,ngz-1
         kf = iz-1
         do iy = 2,ngy-1
            jf = iy-1
            do ix = 2,ngx-1
               if = ix-1
               wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
               y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*a1
     &                     + x(ix  ,iy-1,iz  )*a2
     &                     + x(ix-1,iy  ,iz  )*a3
     &                     + x(ix  ,iy  ,iz  )*(a4 - wnum**2)
     &                     + x(ix+1,iy  ,iz  )*a5
     &                     + x(ix  ,iy+1,iz  )*a6
     &                     + x(ix  ,iy  ,iz+1)*a7
            enddo
         enddo
      enddo
c
c     FACES::
c     N : North
      iy = ngy
      jf = iy-1
      do iz = 2,ngz-1
         kf = iz-1
         do ix = 2,ngx-1
            if = ix-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
            y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &                  + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &                  + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 2.d0*cidd/
     &                               (wnum*hy)*((wnum*hy)**2-2.d0)*hy2)
     &                  + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2) 
c       print*,'iif',ix,iy,iz	 
         enddo
      enddo
c       stop	 

c
c     S : South
      iy = 1
      jf = iy-1
      do iz = 2,ngz-1
         kf = iz-1
         do ix = 2,ngx-1
            if = ix-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)     
            y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &                  + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 +2.d0*cidd/
     &                                (wnum*hy)*((wnum*hy)**2-2.d0)*hy2)
     &                  + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &                  + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)  
         enddo
      enddo
c
c     E : East
      ix = ngx
      if = ix-1
      do iz = 2,ngz-1
         kf = iz-1
         do iy = 2,ngy-1
            jf = iy-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)   
            y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &                  + x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &                  + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 2.d0*cidd/
     &                              (wnum*hx)*((wnum*hx)**2-2.d0)*hx2)
     &                  + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
     &                  + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)  
         enddo
      enddo
c
c     W : West
      ix = 1
      if = ix-1
      do iz = 2,ngz-1
         kf = iz-1
         do iy = 2,ngy-1
            jf = iy-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)                        
            y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &                  + x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 2.d0*cidd/
     &                              (wnum*hx)*((wnum*hx)**2-2.d0)*hx2)
     &                  + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &                  + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
     &                  + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)   
         enddo
      enddo
c
c     R : Right
      iz = ngz
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         do ix = 2,ngx-1
            if = ix-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)      
            y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &                  + x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &                  + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 2.d0*cidd/
     &                              (wnum*hz)*((wnum*hz)**2-2.d0)*hz2)
     &                  + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)  
         enddo
      enddo
c
c     L : Left
      iz = 1
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         do ix = 2,ngx-1
            if = ix-1
            wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)         
            y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &                  + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 2.d0*cidd/
     &                              (wnum*hz)*((wnum*hz)**2-2.d0)*hz2)
     &                  + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &                  + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
     &                  + x(ix  ,iy  ,iz+1)*(a7 - hz2)   
         enddo
      enddo
c
c     EDGES::
c     North face::
c     (2,NGY,1) ---> (NGX-1,NGY,1)
      iy = ngy
      jf = iy-1
      iz = 1
      kf = iz-1
      do ix = 2,ngx-1
         if = ix-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 - hy2)
     &               + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hx)
     &                              *(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2)
     &               + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy  ,iz+1)*(a7 - hz2) 
      enddo
c
c     (2,NGY,NGZ) ---> (NGX-1,NGY,NGZ)
      iy = ngy
      jf = iy-1
      iz = ngz
      kf = iz-1
      do ix = 2,ngx-1
         if = ix-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &               + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &               + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hx)
     &                               *(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2)
     &               + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
      enddo
c
c     (1,NGY,2) ---> (1,NGY,NGZ-1)
      ix = 1
      if = ix-1
      iy = ngy
      jf = iy-1
      do iz = 2,ngz-1
         kf = iz-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &               + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hz)
     &                               *(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2)
     &               + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &               + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2) 
      enddo        
c
c     (NGX,NGY,2) ---> (NGX,NGY,NGZ-1)
      ix = ngx
      if = ix-1
      iy = ngy
      jf = iy-1
      do iz = 2,ngz-1
         kf = iz-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &               + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &               + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hz)
     &                              *(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2)
     &               + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)     
      enddo
c
c     South Face::  
c     (2,1,1) ---> (NGX-1,1,1)
      iy = 1
      jf = iy-1
      iz = 1
      kf = iz-1
      do ix = 2,ngx-1
         if = ix-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)  
         y(ix,iy,iz) = x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hx)
     &                              *(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2)
     &               + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &               + x(ix  ,iy  ,iz+1)*(a7 - hz2)
      enddo
c
c     (2,1,NGZ) ---> (NGX-1,1,NGZ)
      iy = 1
      jf = iy-1
      iz = ngz
      kf = iz-1
      do ix = 2,ngx-1
         if = ix-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &               + x(ix-1,iy  ,iz  )*(a3 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hx)
     &                              *(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2)
     &               + x(ix+1,iy  ,iz  )*(a5 + cidd/(wnum*hx)*hx2)
     &               + x(ix  ,iy+1,iz  )*(a6 - hy2)
      enddo
c
c     (1,1,2) ---> (1,1,NGZ-1)
      ix = 1
      if = ix-1
      iy = 1
      jf = iy-1   
      do iz = 2,ngz-1
         kf = iz-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)    
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hz)
     &                              *(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2)
     &               + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &               + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &               + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)  
      enddo 
c
c     (NGX,1,2) ---> (NGX,1,NGZ-1)
      ix = ngx
      if = ix-1
      iy = 1
      jf = iy-1
      do iz = 2,ngz-1
         kf = iz-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 + cidd/(wnum*hz)*hz2)
     &               + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hz)
     &                              *(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2)
     &               + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &               + x(ix  ,iy  ,iz+1)*(a7 + cidd/(wnum*hz)*hz2)         
      enddo
c
c     Left face::
c     (1,2,1) ---> (1,NGY-1,1)
      ix = 1
      if = ix-1
      iz = 1
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hy)
     &                              *(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2)
     &               + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &               + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
     &               + x(ix  ,iy  ,iz+1)*(a7 - hz2)
      enddo
c
c     (NGX,2,1) ---> (NGX,NGY-1,1)
      ix = ngx
      if = ix-1
      iz = 1
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &               + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hy)
     &                              *(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2)
     &               + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
     &               + x(ix  ,iy  ,iz+1)*(a7 - hz2)
      enddo
c
c     Right face:
c     (1,2,NGZ) ---> (1,NGY-1,NGZ)
      ix = 1
      if = ix-1
      iz = ngz
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf) 
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &               + x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hy)
     &                               *(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2)
     &               + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &               + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2) 
      enddo
c
c     (NGX,2,NGZ) ---> (NGX,NGY-1,NGZ)
      ix = ngx
      if = ix-1
      iz = ngz
      kf = iz-1
      do iy = 2,ngy-1
         jf = iy-1
         wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)      
         y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &               + x(ix  ,iy-1,iz  )*(a2 + cidd/(wnum*hy)*hy2)
     &               + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &               + x(ix  ,iy  ,iz  )*(a4 - wnum**2 - 1.d0/(wnum*hy)
     &                              *(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2)
     &               + x(ix  ,iy+1,iz  )*(a6 + cidd/(wnum*hy)*hy2)
      enddo
c
c     CORNERS:
c
c     (1,1,1)
      ix = 1
      if = ix-1
      iy = 1
      jf = iy-1
      iz = 1
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &            + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &            + x(ix  ,iy  ,iz+1)*(a7 - hz2)
c
c     (NGX,1,1)
      ix = ngx
      if = ix-1
      iy = 1
      jf = iy-1
      iz = 1
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix-1,iy  ,iz  )*(a3 - hx2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix  ,iy+1,iz  )*(a6 - hy2)
     &            + x(ix  ,iy  ,iz+1)*(a7 - hz2)   
c
c     (1,NGY,1)
      ix = 1
      if = ix-1
      iy = ngy
      jf = iy-1
      iz = 1
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)   
      y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 - hy2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &            + x(ix  ,iy  ,iz+1)*(a7 - hz2)
c
c     (NGX,NGY,1)
      ix = ngx
      if = ix-1
      iy = ngy
      jf = iy-1
      iz = 1
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix  ,iy-1,iz  )*(a2 - hy2)
     &            + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix  ,iy  ,iz+1)*(a7 - hz2)
c
c     (1,1,NGZ)      
      ix = 1
      if = ix-1
      iy = 1
      jf = iy-1
      iz = ngz
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix+1,iy  ,iz  )*(a5 - hx2)
     &            + x(ix  ,iy+1,iz  )*(a6 - hy2)
c
c     (NGX,1,NGZ)
      ix = ngx
      if = ix-1
      iy = 1
      jf = iy-1
      iz = ngz
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)     
      y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &            + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix  ,iy+1,iz  )*(a6 - hy2)
c
c     (1,NGY,NGZ)
      ix = 1
      if = ix-1
      iy = ngy
      jf = iy-1
      iz = ngz
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &            + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2)
     &            + x(ix+1,iy  ,iz  )*(a5 - hx2)
c
c     (NGX,NGY,NGZ)
      ix = ngx
      if = ix-1
      iy = ngy
      jf = iy-1
      iz = ngz
      kf = iz-1
      wnum = wdata(omg,nx,ny,nz,idata,if,jf,kf)
      y(ix,iy,iz) = x(ix  ,iy  ,iz-1)*(a1 - hz2)
     &            + x(ix  ,iy-1,iz  )*(a2 - hy2)
     &            + x(ix-1,iy  ,iz  )*(a3 - hx2)
     &            + x(ix  ,iy  ,iz  )*(a4 - wnum**2 + 
     &                                 4.d0*cidd*wnum*hh*hh2) 
c   
c      print*,'finish mat3dmf',y(ix,iy,iz)
   
c      do iz=1,ngz
c       do iy=1,ngy
c        do ix=1,ngx
c	 write(19,*)ix,iy,iz,y(ix,iy,iz)
c	end do
c       end do
c      end do  	
            
      return
      end
      
     
     
     
      
