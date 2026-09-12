      subroutine mat3d7p(ngx,ngy,ngz,x,y,a,work)
c
c     The subroutine returns the Y = A*X, where A is a matrix related to
c     a nine-point Finite Difference stencil.
c     
      implicit    none
c      
      integer     ix,iy,iz,n,ngx,ngy,ngz 
      complex*16  a(ngx,ngy,ngz,7) 
      complex*16  work(0:ngx+1,0:ngy+1,0:ngz+1) 
      complex*16  x(ngx,ngy,ngz),y(ngx,ngy,ngz)  
c      real*8      ww,xx,yy,dnorm2  
c
c      write(*,*) ngx,ngy,ngz 
      do iz = 0,ngz+1
         do iy = 0,ngy+1
            do ix = 0,ngx+1
               work(ix,iy,iz) = (0.d0,0.d0)
            enddo
         enddo
      enddo
c      xx = dnorm2(x,n)
c      write(*,*) 'xx = ',xx
c      
c      n = ngx*ngy*ngz 
      do iz = 1,ngz
         do iy = 1,ngy
            do ix = 1,ngx
               work(ix,iy,iz) = x(ix,iy,iz)
            enddo
         enddo
      enddo
      
c      ww = dnorm2(work,n)
c      write(*,*) 'ww = ',ww
C 
      do iz = 1,ngz
         do iy = 1,ngy
            do ix = 1,ngx
               y(ix,iy,iz) = a(ix,iy,iz,1)*work(ix  ,iy  ,iz-1)
     &                     + a(ix,iy,iz,2)*work(ix  ,iy-1,iz  )
     &                     + a(ix,iy,iz,3)*work(ix-1,iy  ,iz  )
     &                     + a(ix,iy,iz,4)*work(ix  ,iy  ,iz  )
     &                     + a(ix,iy,iz,5)*work(ix+1,iy  ,iz  )
     &                     + a(ix,iy,iz,6)*work(ix  ,iy+1,iz  )
     &                     + a(ix,iy,iz,7)*work(ix  ,iy  ,iz+1)
            enddo
         enddo
      enddo
c      read(*,*)
      
c      yy = dnorm2(y,n)
c      write(*,*) 'yy = ',yy
c      
      return
      end
