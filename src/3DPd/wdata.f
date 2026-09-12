      function wdata(omg,ngx,ngy,ngz,idata,ix,iy,iz)
c      
      implicit none
      integer   ix,iy,iz,idata,ngx,ngy,ngz 
      real*8    hx,hy,hz,omg,xh,yh,zh,wdata 
      real*8    ya,yb 
c      
      hx = 1.d0/dble(float(ngx+1))
      hy = 1.d0/dble(float(ngy+1))
      hz = 1.d0/dble(float(ngz+1))
      if (idata .eq. 0) then
         wdata = omg
      elseif (idata .eq. 1) then
         yh = hy*dble(float(iy))
         if (yh .le. 1.d0/3.d0) then
            wdata = 1.2d0*omg
         elseif (yh. gt. 1.d0/3.d0 .and. yh .le. 2.d0/3.d0) then
            wdata = omg
         elseif (yh .gt. 2.d0/3.d0) then
            wdata = 1.5d0*omg
         endif
      elseif (idata .eq. 2) then
         zh = hz*dble(float(iz))
         if (zh .le. 1.d0/3.d0) then
            wdata = 1.2d0*omg
         elseif (zh. gt. 1.d0/3.d0 .and. zh .le. 2.d0/3.d0) then
            wdata = omg
         elseif (zh .gt. 2.d0/3.d0) then
            wdata = 1.5d0*omg
         endif
      elseif (idata .eq. 3) then
         xh = hx*dble(float(ix))
         yh = hy*dble(float(iy))
         zh = hz*dble(float(iz))
         ya = (1.d0 - 0.5d0*xh - 0.375d0*zh)/2.5d0
         yb = (1.d0 + 1.d0/6.d0*xh + 1.d0/3.d0*zh)/(5.d0/3.d0)
         if (yh .lt. ya) then
            wdata = 1.25d0*omg
         elseif (yh .ge. ya .and. yh .le. yb) then
            wdata = omg
         elseif (yh .gt. yb) then
            wdata = 2.0d0*omg
         endif      
      endif
c      
      return
      end
      
      
