      subroutine defmsg3d(nz,ngx,ngy,mx,my, 
     &                    nx,ny,nallx2,nally2, 
     &                    npz,npx,npy,ipz,ipx,ipy, 
     &                    ixs,iys,ix0,iy0,ixg,iyg,
     &                    izg, 
     &                    neigz,neigx,neigy) 
c
      implicit real*8 (a-h,o-z)
      integer     nx(mx),ny(my),ixs(mx),iys(my) 
      integer     ix0(mx),iy0(my),ixg(mx),iyg(my) 
      integer     neigz(2),neigx(2,mx),neigy(2,my) 
      integer     nza(0:255),nxa(0:255),nya(0:255) 
c
      ipe  = ipz+ipx*npz+ipy*npz*npx
      ierr  = 0
      icomm = 0
c
c     Make neigz
c
      if (ipz.ne.0) then
         neigz(1) = ipe - 1
      else
         neigz(1) = -1
      endif
      if (ipz.ne.npz-1) then
         neigz(2) = ipe + 1
      else
         neigz(2) = -1
      endif
      nz0 = nz
      nx0 = ngx
      ny0 = ngy
      call findis3d(npz,npx,npy,ipz,ipx,ipy, 
     &              nz0,nx0, ny0,izg0,ixg0,iyg0,icomm)
      izg = izg0
c
c     Make x y directions
c
      do kx = mx,1,-1
c         if (kx.eq.mx) then
c            nx(kx) = ngx
c            ixs(kx) = 1
c         else
c            nx(kx) = (nx(kx+1) - ix0(kx+1) + 2)/2
c            ixs(kx) = ixs(kx+1) + nx(kx+1) + 1
c         endif
c         
         if (kx.eq.mx) then
            nx(kx)  = ngx
            ixs(kx) = 1             
         elseif (kx.eq.(mx-1)) then
            nx(kx)  = (nx(kx+1) - ix0(kx+1) + 2)/2
            ixs(kx) = ixs(kx+1) + nx(kx+1) + 1
            ixsd = ixs(kx)
            nxd  = nx(kx) + 1
         elseif (kx.eq.(mx-2)) then
            nx(kx)  = (nx(kx+1) - ix0(kx+1) + 2)/2
            ixs(kx) = ixs(kx+1) + nx(kx+1) + 1 - nxd
         else
            nx(kx)  = (nx(kx+1) - ix0(kx+1) + 2)/2
            ixs(kx) = ixs(kx+1) + nx(kx+1) + 1
         endif
c         
         nz0 = nz
         nx0 = nx(kx)
         ny0 = 0
         call findis3d(npz,npx,npy,ipz,ipx,ipy, 
     &                 nz0,nx0,ny0,izg0,ixg0,iyg0,icomm)
         ixg(kx) = ixg0
         ix0(kx) = mod(ixg(kx)+1,2) + 1
c
         neigx(1,kx) = -1
         neigx(2,kx) = -1
c
      enddo
c
      do ky = my,1,-1
c
c         if (ky.eq.my) then
c            ny(ky) = ngy
c            iys(ky) = 1
c         else
c            ny(ky) = (ny(ky+1) - iy0(ky+1) + 2)/2
c            iys(ky) = iys(ky+1) + ny(ky+1) + 1
c         endif
c
         if (ky.eq.my) then
            ny(ky) = ngy
            iys(ky) = 1
         elseif (ky.eq.(my-1)) then
            ny(ky) = (ny(ky+1) - iy0(ky+1) + 2)/2
            iys(ky) = 1
         else
            ny(ky) = (ny(ky+1) - iy0(ky+1) + 2)/2
            iys(ky) = ny(my-1)+2
         endif
c
         nz0 = nz
         nx0 = 0
         ny0 = ny(ky)
         call findis3d(npz,npx,npy,ipz,ipx,ipy, 
     &                 nz0,nx0,ny0,izg0,ixg0,iyg0,icomm)
         iyg(ky) = iyg0
         iy0(ky) = mod(iyg(ky)+1,2) + 1
c
         neigy(1,ky) = -1
         neigy(2,ky) = -1
c
      enddo
c
c      nallx2 = ixs(1) + nx(1) + 1
c      nally2 = iys(1) + ny(1) + 1
      nxs1 = ixs(1) + nx(1) + 1
      nxs2 = ixs(mx-1) + nx(mx-1) + 1
      
      if (nxs1 .ge. nxs2) then
         nallx2 = nxs1
      else
         nallx2 = nxs2
      endif
c
c      nallx2 = ixs(1) + nx(1) + 1
c      nally2 = iys(1) + ny(1) + 1
      nally2 = ny(my) + 1
c      
      return
      end
