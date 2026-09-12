      subroutine vcycle3(f,w,acor,dinv,dc,mmax,kmg,kmax,mx,my,md,
     &                   nz,nallx2,nally2,
     &                   nx,ny,idata,omg,bet1,bet2,
     &                   ixs,iys,
     &                   ix0,iy0,
     &                   omega,iswp,icyc,
     &                   nr1,nr2,nr3,nr4,
     &                   wrst,witpd,r,
     &                   ip,neigz,neigx,neigy,ipe,ipz,npz) 
c     
      implicit  real*8 (a-h,o-z)
c      
      integer     nz,nallx2,nally2,mx,my,kmg 
      integer     nr1,nr2,nr3,nr4 
c      complex*16  f(*),w((nz+2)*nallx2*nally2),r(*) 
      complex*16  f(*),w(*),r(*) 
      complex*16  acor(*),dinv(*),dc(*) 
      complex*16  wrst(*),witpd(*) 
      integer     nx(mx),ny(my),ixs(mx),iys(my),ix0(mx),iy0(my) 
      integer     ip(*),neigz(2),neigy(2*mx),neigx(2*my) 
c       
      mxy = min(mx,my)
      kx  = mx - mxy + kmg  
      ky  = my - mxy + kmg  
      nz2 = nz + 2 
c
      nxf = nx(kx)
      nyf = ny(ky)
c
c     Initialize solution vector
c      
      do j = iys(ky)+1,iys(ky)+nyf
         do i = ixs(kx)+1,ixs(kx)+nxf
            do k = 1,nz
	             w(k+(i+(j-1)*(nallx2)-1)*nz2+1) = (0.d0,0.d0)
            enddo
         enddo
      enddo
c
      nr = max(nr1,nr2)
c
      if (kmg .eq. 1) then
         if (kmg .eq. mmax) then
            if (iswp .eq. 0) then        
               call sor27msgf(!acor,
     &                        dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,idata,omg,bet1,bet2,
     &                        neigz,neigx,neigy)
            elseif (iswp .eq. 1) then
               call jac27msgf(!acor,
     &                        dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,idata,omg,bet1,bet2,
     &                        neigz,neigx,neigy)
            endif
         else
            if (iswp .eq. 0) then        
               call sor27msgc(acor,dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,
     &                        neigz,neigx,neigy)
            elseif (iswp .eq. 1) then
               call jac27msgc(acor,dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,
     &                        neigz,neigx,neigy)
            endif    
         endif
         return
      endif
c      
      if (nr1 .gt. 0) then
c
c        Pre-smoothing
c
         if (kmg .eq. mmax) then
            if (iswp .eq. 0) then
               call sor27msgf(!acor,
     &                        dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr1, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,idata,omg,bet1,bet2,
     &                        neigz,neigx,neigy)
            elseif (iswp .eq. 1) then
               call jac27msgf(!acor,
     &                        dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr1, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,idata,omg,bet1,bet2,
     &                        neigz,neigx,neigy)
            endif
            call res27msgf(!acor,
     &                     idata,omg,bet1,bet2,w,f,r,nallx2,nally2,
     &                     ixs(kx),iys(ky),nx(kx),ny(ky),nz, 
     &                     neigz,neigx(2*kx), neigy(2*ky))
         else
            if (iswp .eq. 0) then
               call sor27msgc(acor,dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr1, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,
     &                        neigz,neigx,neigy)
            elseif (iswp .eq. 1) then
               call jac27msgc(acor,dinv,dc,w,f,r, 
     &                        nz,nallx2,nally2, 
     &                        omega,nr1, 
     &                        ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                        ix0(kx),iy0(ky), 
     &                        ipe,ipz,npz,
     &                        neigz,neigx,neigy)
            endif
            call res27msgc(acor,w,f,r,nallx2,nally2,
     &                     ixs(kx),iys(ky),nx(kx),ny(ky),nz, 
     &                     neigz,neigx(2*kx), neigy(2*ky))
         endif 
c
c        Compute residual
c
c         call res27msg(acor,w,f,r,nallx2,nally2,
c     &                 ixs(kx),iys(ky),nx(kx),ny(ky),nz, 
c     &                 neigz,neigx(2*kx),neigy(2*ky))
      else
         do j = iys(ky)+1,iys(ky)+nyf
            do i = ixs(kx)+1,ixs(kx)+nxf
	             do k = 1,nz
	                r(k+(i+(j-1)*(nallx2)-1)*nz2+1) =
     &                      f(k+(i+(j-1)*(nallx2)-1)*nz2+1)
               enddo
            enddo
         enddo
      endif
c
c     Restrict the residual to the coarse grid
c    
      call restr3a(r,f,wrst,nallx2,nally2,
     &             nx(kx-1),nx(kx),ixs(kx-1),ixs(kx),ix0(kx),
     &             ny(ky-1),ny(ky),iys(ky-1),iys(ky),iy0(ky),
     &             nz,neigx(2*kx),neigy(2*ky))
c
c     call fcycmgs recursively
c
      kmgc = kmg-1
c
      call msg3recs(f,w,acor,dinv,dc,
     &              mmax,kmgc,kmax,mx,my,md,
     &              nz,nallx2,nally2,nx,ny,idata,omg,bet1,bet2,
     &              ixs,iys,ix0,iy0,
     &              omega,iswp,0,
     &              nr1,nr2,nr3,nr4,
     &              wrst,
     &        	    witpd,r,
     &              ip,
     &              neigz,
     &              neigx,neigy,
     &              ipe,ipz,npz)
c
c     Interpolation from coarse grids w_(kmg-1) to fine grids w_kmg
c
      call extup3a(w,r,witpd,nallx2,nally2, 
     &             nx(kx-1),nx(kx),ixs(kx-1),ixs(kx),ix0(kx),
     &	           ny(ky-1),ny(ky),iys(ky-1),iys(ky),iy0(ky),
     &          	 nz,neigx(2*kx),neigy(2*ky))
c
c     Post-smoother
c
      if (kmg .eq. mmax) then
         if (iswp .eq. 0) then
            call sor27msgf(!acor,
     &                     dinv,dc,w,f,r, 
     &                     nz,nallx2,nally2, 
     &                     omega,nr2, 
     &                     ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                     ix0(kx),iy0(ky), 
     &                     ipe,ipz,npz,idata,omg,bet1,bet2,
     &                     neigz,neigx,neigy)
         elseif (iswp .eq. 1) then
            call jac27msgf(!acor,
     &                     dinv,dc,w,f,r, 
     &                     nz,nallx2,nally2, 
     &                     omega,nr2, 
     &                     ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                     ix0(kx),iy0(ky), 
     &                     ipe,ipz,npz,idata,omg,bet1,bet2,
     &                     neigz,neigx,neigy)
         endif
      else
         if (iswp .eq. 0) then
            call sor27msgc(acor,dinv,dc,w,f,r, 
     &                     nz,nallx2,nally2, 
     &                     omega,nr2, 
     &                     ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                     ix0(kx),iy0(ky), 
     &                     ipe,ipz,npz,
     &                     neigz,neigx,neigy)
         elseif (iswp .eq. 1) then
            call jac27msgc(acor,dinv,dc,w,f,r, 
     &                     nz,nallx2,nally2, 
     &                     omega,nr2, 
     &                     ixs(kx),iys(ky),nx(kx),ny(ky), 
     &                     ix0(kx),iy0(ky), 
     &                     ipe,ipz,npz,
     &                     neigz,neigx,neigy)
         endif
      endif
c
      return
      end
