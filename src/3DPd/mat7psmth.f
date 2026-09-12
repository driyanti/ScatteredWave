      subroutine mat7pda(y,x,idata,omg,bet1,bet2,
     &                   ix,iy,iz,w1,w2,w6,w7,ngx,ngy,ngz)
c
      implicit    none 
      integer     ix,iy,iz,ixo,iyo,izo,ngx,ngy,ngz,ngx2,ngy2,ngz2,idata 
      real*8      omg,hx,hy,hz,hx2,hy2,hz2,wnum,wdata,bet1,bet2 
      complex*16  w1,w2,w6,w7,a1,a2,a6,a7,cidd   
      complex*16  x,y 
      parameter   (cidd = (0.d0,1.d0))
c      
c     receive i,j,k. Sweep subsequently in z-,x-, and y-direction (or
c     k-,i-, and j-direction, respectively).
c
c     compute gridsize paramater
c          
      ngx2 = ngx-2
      ngy2 = ngy-2
      ngz2 = ngz-2
c      write(*,*) 
c      
      hx = 1.d0/dble(float(ngx2+1))
      hy = 1.d0/dble(float(ngy2+1))
      hz = 1.d0/dble(float(ngz2+1))
      hx2 = 1.d0/hx**2
      hy2 = 1.d0/hy**2
      hz2 = 1.d0/hz**2
c
      ixo = ix-1
      iyo = iy-1
      izo = iz-1
c      
      if (iy .eq. 1) then
         if (ix .eq. 1) then
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)  
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a6*w6 - a7*w7  
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a6*w6 - a7*w7  
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a6*w6 - a7*w7 
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then 
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a7*w7
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a7*w7
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a7*w7
            endif
         else
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a6*w6 - a7*w7
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a6*w6 - a7*w7
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               y  = x - a2*w2 - a6*w6 - a7*w7
            endif
         endif
      elseif (iy .eq. ngy) then
         if (ix .eq. 1) then 
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a6*w6
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a6*w6
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a6*w6
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a2*w2
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a2*w2
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
               y  = x - a1*w1 - a2*w2
            endif
         else
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               y  = x - a1*w1 - a2*w2 - a6*w6
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo) 
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
               y  = x - a1*w1 - a2*w2 - a6*w6
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
               a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               y  = x - a1*w1 - a2*w2 - a6*w6
            endif
         endif
      else
         if (ix .eq. 1) then
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a6*w6 - a7*w7
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a6*w6 - a7*w7
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a6*w6 - a7*w7
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a2*w2 - a7*w7
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a2*w2 - a7*w7
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a2*w2 - a7*w7
            endif
         else
            if (iz .eq. 1) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a2*w2 - a6*w6 - a7*w7
            elseif (iz .eq. ngz) then
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
               a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
               y  = x - a1*w1 - a2*w2 - a6*w6 - a7*w7
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               a1 = - cmplx(hy2,0.d0) 
               a2 = - cmplx(hx2,0.d0)
               a6 = - cmplx(hx2,0.d0)
               a7 = - cmplx(hy2,0.d0)
               y  = x - a1*w1 - a2*w2 - a6*w6 - a7*w7
            endif
         endif
      endif
c
      return
      end
c      
c***********************************************************************
c
      subroutine mat7pdd(ad,idata,omg,bet1,bet2,
     &                   ix,iy,iz,idz,ngx,ngy,ngz)
      
      implicit none
      integer     ix,iy,iz,ixo,iyo,izo,idz,ngx,ngy,ngz,idata 
      integer     ngx2,ngy2,ngz2 
      complex*16  ad,cidd  
      real*8      hx,hy,hz,hx2,hy2,hz2,omg,wnum,wdata,bet1,bet2 
      parameter   (cidd = (0.d0,1.d0))
      
      ngx2 = ngx-2
      ngy2 = ngy-2
      ngz2 = ngz-2
      hx   = 1.d0/dble(float(ngx2+1))
      hy   = 1.d0/dble(float(ngy2+1))
      hz   = 1.d0/dble(float(ngz2+1))
      hx2  = 1.d0/hx**2
      hy2  = 1.d0/hy**2
      hz2  = 1.d0/hz**2
      ixo  = ix-1
      iyo  = iy-1
      izo  = iz-1
c      
      if (iy .eq. 1) then
         if (ix .eq. 1) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = (0.d0,0.d0)
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                 ad = (0.d0,0.d0)
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         else
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = (0.d0,0.d0)
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                  ad = (0.d0,0.d0)
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2   
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2   
               endif
            endif
         endif
      elseif (iy .eq. ngy) then
         if (ix .eq. 1) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               endif 
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         else
            if (iz .eq. 1) then
               if (idz .eq. -1) then 
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0) 
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         endif
      else
         if (ix .eq. 1) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         elseif (ix .eq. ngx) then
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         else
            if (iz .eq. 1) then
               if (idz .eq. -1) then
                  ad = 0.d0
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               endif
            elseif (iz .eq. ngz) then
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
               elseif (idz .eq. 1) then
                  ad = 0.d0
               endif
            else
               wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
               if (idz .eq. -1) then
                  ad = - cmplx(hz2,0.d0) !+ 1.d0*cidd/(wnum*hz)*hz2
               elseif (idz .eq. 1) then
                  ad = - cmplx(hz2,0.d0) !+ 1.d0*cidd/(wnum*hz)*hz2
               endif
            endif
         endif
      endif

c     wnum = wdata(omg,ngx2,ngy2,ngz2,idata,izo,ixo,iyo)
c      if (iz .eq. 1) then
c         if (idz .eq. -1) then
c            ad = 0.d0
c         elseif (idz .eq. 1) then
c            ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
c         endif
c      elseif (iz .eq. ngz) then
c         if (idz .eq. -1) then
c            ad = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
c         elseif (idz .eq. 1) then
c            ad = 0.d0
c         endif
c      else
c         ad = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
c      endif
      return
      end
c      
        
