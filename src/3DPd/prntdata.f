      subroutine prntdata(x,ngx,ngy,ngz)
      
      
      implicit none
      integer     ix,iy,iz,n,ngx,ngy,ngz 
      complex*16  x(ngx,ngy,ngz) 
      real*8      hx,hy,hz,dx,dy,dz 
      
      hx = 1.d0/dble(float(ngx+1))
      hy = 1.d0/dble(float(ngy+1))
      hz = 1.d0/dble(float(ngz+1))
c
      open(10,file='solution.dat')
      iz = (ngz+1)/2      
      do iy = 1,ngy
         do ix = 1,ngx
            dx = ix*hx
            dy = iy*hy
            write(10,10) dx,dy,realpart(x(ix,iy,iz)),
     &                   imagpart(x(ix,iy,iz))
         enddo
      enddo
      close(10)
  
   10 format(d12.5,' ',d12.5,' ',d12.5,' ',d12.5)   
      return
      end
