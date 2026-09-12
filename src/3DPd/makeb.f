      subroutine makeb(b,ngx,ngy,ngz,irhs) 
      
      implicit    none 
      integer     ix,iy,iz,ngx,ngy,ngz,irhs  
      complex*16  b(ngx,ngy,ngz) 
      real*8      hx,hy,hz 
c      
      hx = 1.d0/dble(float(ngx+1))
      hy = 1.d0/dble(float(ngy+1))
      hz = 1.d0/dble(float(ngz+1))
c
      if (irhs .eq. 0) then
         do iz = 1,ngz
            do iy = 1,ngy
               do ix = 1,ngx
                  if (ix .eq. (ngx+1)/2 .and.
     &                iy .eq. (ngy+1)/2 .and.
     &                iz .eq. (ngz+1)/2) then
                     b(ix,iy,iz) = (1.d0,0.d0)/(hx*hy*hz)
                  else
                     b(ix,iy,iz) = (0.d0,0.d0)
                  endif
               enddo
            enddo
         enddo
      elseif (irhs .eq. 1) then
         do iz = 1,ngz
            do iy = 1,ngy
               do ix = 1,ngx
                  if (ix .eq. (ngx+1)/2 .and.
     &                iy .eq. 2 .and.
     &                iz .eq. (ngz+1)/2) then
                     b(ix,iy,iz) = (1.d0,0.d0)/(hx*hy*hz)
                  else
                     b(ix,iy,iz) = (0.d0,0.d0)
                  endif
               enddo
            enddo
         enddo
      endif
      
      return
      end
