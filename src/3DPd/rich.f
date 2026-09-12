c      subroutine mgrid(ia,n,ngx,ngy,nxpml,nypml,xl,yl,
c     &                 a,b,x,ap,inump,itdata,
c     &                 imgrid,rmgrid,mmax,ngall0,ngall1,nipre,nrpre,
c     &                 ipm,ipre,rpre,wpre,w,work)
     
      subroutine rich(b,x,n,bet1,bet2,wreff,idata,
     &                ia,ja,ipm,rpm,omg,res,w,
     &                temp,rpred,ipred,ip)
c
      implicit    none 
      integer     ip(*),ipm(4),ia(3),ja(*),ipred(*),img,idata  
      integer     i,n,ngx,ngy,ngz 
      integer     iter,maxiter 
      real*8      rpm(2),res(*) 
c      complex*16  a(*) 
      complex*16  b(n),x(n),w(n,2),temp(n+2,3) 
      complex*16  rpred(*)       
c      complex*16  rho0,rho1,alpha,beta,omega,r0v,ts,tt 
c      complex*16  herpr      
c
      real*8      bnorm,rnorm,rnorm0 
      real*8      dnorm2,bb,rr,rr0,rel,tol,omg,wreff 
      real*8      bet1,bet2   
c

      iter    = 0
      ngx     = ia(1)
      ngy     = ia(2)
      ngz     = ia(3)
      maxiter = 300
      tol     = 1.0d-6

      open(1,file='residual.out')
c
      bnorm = dnorm2(b,n)
      if (bnorm .lt. 1.d-32) return
c      
c      call mat3d7p(ngx,ngy,ngz,x,w(1,1),a,temp(1,1))
      call mat3d7pfp(ngx,ngy,ngz,bet1,bet2,wreff,idata,
     &               x,w(1,1),temp(1,1))
c
      do i = 1,n
         w(i,1) = b(i) - w(i,1)
      enddo
c
      rnorm  = dnorm2(w(1,1),n)
      rnorm0 = rnorm
      rel = rnorm/bnorm
      write(*,100) iter, rnorm, rel,rnorm/rnorm0
c      
c      call pmgrid(n,ngx,ngy,ap,ia,xl,yl,
c     &            imgrid,rmgrid,mmax,ngall0,ngall1,nipre,nrpre,    
c     &            wpre,rpre,ipre,ip,work)
c     
c      call cputime(tend)
c      time(1) = time(1) + (tend-tinit)
c      time(5) = time(5) + (tend-tinit)
c
c      call cputime(tinit)
      do while (iter .lt. maxiter .and. rel .gt. tol)     
c         
         iter = iter + 1
c         
c         call cputime(t3)
c         call pmg2d(n,xl,yl,w(1,1),w(1,2),ap,wrelax,ia,
c     &              wpre,rpre,ipre,work)
         call pmsa3d19(n,w(1,1),w(1,2),!Ap,
     &                 ia,ja,
     &                 temp,rpred,ipred,ip,omg)
c         call cputime(t4)
c         time(3) = time(3) + (t4-t3) 
c     
         do i = 1,n
            x(i) = x(i) + w(i,2)
         enddo
c      
c         call cputime(t3)
c         call mat3d7p(ngx,ngy,ngz,x,w(1,1),a,temp(1,1))      
         call mat3d7pfp(ngx,ngy,ngz,bet1,bet2,wreff,idata,
     &                  x,w(1,1),temp(1,1))
c         call cputime(t4)
c         time(2) = time(2) + (t2-t1)
c         nmatvec = nmatvec + 1
c
         do i = 1,n
            w(i,1) = b(i) - w(i,1)
         enddo
c      
         rnorm = dnorm2(w(1,1),n)
         rel =  rnorm/bnorm
c         
         write(*,100)  iter, rnorm, rel, rnorm/rnorm0
         write(1,100) iter, rnorm, rel, rnorm/rnorm0
         rnorm0 = rnorm
c
      enddo
      close(1)
c
c      call cputime(tend)
c      time(4) = time(4) + (tend-tinit)
c      time(5) = time(5) + (tend-tinit)
c
c      itdata(1) = iter
c      itdata(2) = nmatvec
c      resdata(1) = rnorm
c      resdata(2) = rel
c
c      call prntdata(ia,n,nxpml,nypml,xl,yl,inump,itdata,time,x,resdata)
c      
 100  format(' ',i7,'   ',d14.8,'   ',d14.8,'   ',d14.8)  
    
      return
      end
