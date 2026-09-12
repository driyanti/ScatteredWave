      subroutine bicgstab(b,x,n,bet1,bet2,wreff,idata,
     &                    ia,ja,ipm,rpm,omg,res,w,
     &                    temp,rpred,ipred,ip)
      
      implicit none 
      integer     ip(*),ipm(4),ia(3),ja(*),ipred(*),img 
      integer     j,n,ngx,ngy,ngz,idata 
      integer     iter,maxiter 
      real*8      rpm(2),res(*) 
c      complex*16  a(*) 
      complex*16  b(n),x(n),w(n,7),temp(n+2,3) 
      complex*16  rpred(*)       
      complex*16  rho0,rho1,alpha,beta,omega,r0v,ts,tt 
      complex*16  herpr      
c
      real*8      dnorm2,bb,rr,rr0,rel,tol,omg,wreff,bet1,bet2   
c
      img  = 1
c      img  = 0
      iter = 0
      tol  = 1.0d-6
      maxiter = 10000
c      maxiter = 2
      open(1,file='residual.dat') 
      ngx = ia(1)
      ngy = ia(2)
      ngz = ia(3) 
      
      write(*,*) n,ngx,ngy,ngz,ngx*ngy*ngz  
c
      bb = dnorm2(b,n)
c      
      call mat3d7pf(ngx,ngy,ngz,wreff,idata,x(1),w(1,5),temp(1,1))
c      
      do j = 1,n
         w(j,5) = b(j) - w(j,5)
         w(j,1) = w(j,5)
         w(j,2) = (0.d0,0.d0)
         w(j,3) = (0.d0,0.d0)
      enddo
c      
      rr = dnorm2(w(1,1),n)
      rr0 = rr
      rel = rr/rr0
      write(*,100) iter,rr,rel
      write(1,100) iter,rr,rel
c
      rho0  = (1.d0,0.d0)
      alpha = (1.d0,0.d0)
      omega = (1.d0,0.d0)
c
      do while (iter .lt. maxiter .and. rel .gt. tol)
c
         iter = iter + 1
c
         rho1 = herpr(w(1,5),w(1,1),n)
         beta = (rho1/rho0)*(alpha/omega)
         rho0 = rho1 
c         
         do j = 1,n
            w(j,3) = w(j,1) + beta*(w(j,3) - omega*w(j,2))
         enddo
	 
c      
         if (img .eq. 0) then
            do j = 1,n
               w(j,4) = w(j,3)
            enddo 
         else
            call pmsa3d19(n,w(1,3),w(1,4),bet1,bet2,idata,wreff,!Ap,
     &                    ia,ja,
     &                    temp,rpred,ipred,ip,omg)
         endif
c

c         call mat3d7p(ngx,ngy,ngz,w(1,4),w(1,2),a,temp(1,1))
         call mat3d7pf(ngx,ngy,ngz,wreff,idata,w(1,4),w(1,2),temp(1,1))
c
	 
         r0v = herpr(w(1,5),w(1,2),n)
         alpha = rho0/r0v

c
         do j = 1,n
            w(j,6) = w(j,1) - alpha*w(j,2)
	 write(17,*)j,w(j,4),w(j,3)	 	    
         enddo
c
         if (img .eq. 0) then
            do j = 1,n
               w(j,1) = w(j,6)
            enddo
         else
            call pmsa3d19(n,w(1,6),w(1,1),bet1,bet2,idata,wreff,!Ap,
     &                    ia,ja,
     &                    temp,rpred,ipred,ip,omg)
         endif
c
c         call mat3d7p(ngx,ngy,ngz,w(1,1),w(1,7),a,temp(1,1)) 
         call mat3d7pf(ngx,ngy,ngz,wreff,idata,w(1,1),w(1,7),temp(1,1))
c
         ts = herpr(w(1,7),w(1,6),n)
         tt = herpr(w(1,7),w(1,7),n)
         omega = ts/tt
c
         do j=1,n
	 write(18,*)j,w(j,1)	 	    
            x(j)   = x(j) + alpha*w(j,4) + omega*w(j,1)
            w(j,1) = w(j,6) - omega*w(j,7)
	 write(19,*)j,w(j,6),w(j,7)	 	    
         enddo
	 
c
         rr  = dnorm2(w(1,1),n)
         rel = rr/rr0
c
         write(*,100) iter,rr,rel
         write(1,100) iter,rr,rel
	 stop	 

c
      enddo

c      
      close(1)      
c      
  100 format(' ',i7,'   ',d14.8,'   ',d14.8)
c
      return
      end
             
      
