      subroutine msg3da(nz,ngx,ngy,nallx2,nally2,p,q,
     &                  acor,dinv,dc, 
     &                  mx,my,nx,ny,idata,omg,bet1,bet2, 
     &                  omega,ids,nr1,nr2,nr3,nr4,niter, 
     &                  icyc,md, 
     &                  wrst,
     &                  witpd,r,w,f, 
     &                  ixs,iys,ix0,iy0,
     &                  iswp, 
     &                  ipz,ipx,ipy,npz,npx,npy, 
     &                  neigz, 
     &                  neigx,neigy) 
c
      implicit   real*8 (a-h,o-z)
      parameter (eps = 1.d-8, izero = 0, ione = 1,
     &           dzero = 0.d0, done = 1.d0, dhalf = 0.5d0)
      complex*16  p(nz,ngx,ngy),q(nz,ngx,ngy) 
      complex*16  acor(nallx2,nally2,27),dinv(nallx2,nally2) 
      integer     nx(mx),ny(my) 
      integer     ixs(mx),iys(my),ix0(mx),iy0(my) 
      complex*16  wrst(0:nz+1,nallx2,nally2,-1:1,-1:1)  
      complex*16  witpd(0:nz+1,nallx2,nally2,4) 
      complex*16  r(0:nz+1,nallx2,nally2) 
      complex*16  w(0:nz+1,nallx2,nally2) 
      complex*16  f(0:nz+1,nallx2,nally2) 
      integer     neigz(2),neigx(2,mx),neigy(2,my) 
c
c     Initiallization 
c
      do iy = 1, ngy
         do ix = 1, ngx
            do iz = 1, nz
               f(iz,ixs(mx)+ix,iys(my)+iy) = p(iz,ix,iy)
               w(iz,ixs(mx)+ix,iys(my)+iy) = q(iz,ix,iy)
            enddo
         enddo
      enddo
c     
      m = min(mx,my)
      mmax = m
      kmax = m
      ipe = ipz + ipx*npz+ipy*npz*npx
c
      do iter = 1, niter
         if (icyc.eq.0) then
            call vcycle3(f,w,acor,dinv,dc, 
     &                   mmax,m,kmax,mx,my,md,
     &                   nz,nallx2,nally2,nx,ny,idata,omg,bet1,bet2,
     &                   ixs,iys,ix0,iy0, 
     &                   omega,iswp,icyc, 
     &                   nr1,nr2,nr3,nr4, 
     &                   wrst,witpd,r, 
     &                   ip, 
     &                   neigz,neigx,neigy, 
     &                   ipe,ipz,npz) 
         elseif (icyc.eq.1) then
            call fcycle3(f,w,acor,dinv,dc, 
     &                   mmax,m,kmax,mx,my,md, 
     &                   nz,nallx2,nally2,nx,ny,idata,omg,bet1,bet2,
     &                   ixs,iys,ix0,iy0, 
     &                   omega,iswp,icyc, 
     &                   nr1,nr2,nr3,nr4, 
     &                   wrst,witpd,r, 
     &                   ip, 
     &                   neigz,neigx,neigy, 
     &                   ipe,ipz,npz)
         endif
      enddo
c
      do iy = 1, ngy
         do ix = 1, ngx
            do iz = 1, nz
               q(iz,ix,iy) = w(iz,ixs(mx)+ix,iys(my)+iy)
               write(19,*),iz,ix,iy,q(iz,ix,iy)	                   
	    enddo
         enddo
      enddo
c
      return
      end
