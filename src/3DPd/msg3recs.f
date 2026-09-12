      subroutine msg3recs(f,w,acor,dinv,dc,
     &              mmax,kmg,kmax,mx, my,md,
     &              nz, nallx2, nally2, nx, ny,idata,omg,bet1,bet2,
     &              ixs, iys, ix0, iy0,
     &              omega, iswp, icyc,
     &              nr1, nr2, nr3, nr4,
     &              wrst,
     &        	    witpd, r,
     &              ip,
     &              neigz,
     &              neigx, neigy,
     &              ipe, ipz, npz)
c     
      complex*16  f(*),w(*),r(*) 
      complex*16  acor(*),dinv(*),dc(*) 
      complex*16  wrst(*),witpd(*) 
      real*8      omega 
      integer     nx(*),ny(*),ixs(*),iys(*),ix0(*),iy0(*) 
      integer     ip(*),neigz(*),neigy(*),neigx(*) 
c      
      if (icyc .eq. 0) then
         call vcycle3(f,w,acor,dinv,dc,mmax,kmg,kmax,mx, my,md,
     &                nz, nallx2, nally2, nx, ny,idata,omg,bet1,bet2,
     &                ixs, iys, ix0, iy0,
     &                omega, iswp,icyc,
     &                nr1, nr2, nr3, nr4,
     &                wrst,witpd, r,
     &                ip,neigz,neigx, neigy,
     &                ipe, ipz, npz)
      elseif (icyc .eq. 1) then
          call fcycle3(f,w,acor,dinv,dc,mmax,kmg,kmax,mx, my,md,
     &                nz, nallx2, nally2, nx, ny,idata,omg,bet1,bet2,
     &                ixs, iys, ix0, iy0,
     &                omega, iswp, icyc,
     &                nr1, nr2, nr3, nr4,
     &                wrst,witpd, r,
     &                ip,neigz,neigx, neigy,
     &                ipe, ipz, npz)
c
      endif
     
      return
      end
