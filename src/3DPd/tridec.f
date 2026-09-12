      subroutine tridec(a,dinv,dc,w1,w2, 
     &                  nz,nallx2,nally2,
     &                  ixs,iys,nx,ny,
     &                  ierr, 
     &                  ipe, ipz, npz)
c
      implicit real*8 (a-h,o-z)
      parameter   (eps=1.d-20)
      complex*16  a(0:nz+1,nallx2,nally2,27) 
      complex*16  dinv(0:nz+1,nallx2,nally2),dc(0:nz+1,nallx2,nally2) 
      complex*16  w1(0:nz+1,nallx2,nally2),w2(0:nz+1,nallx2,nally2) 
      complex*16  d1 
      integer     neigz(2) 
c
      ierr = 0 
      nxy = nx*ny
      if (nxy.eq.0) return
      if (ipz.ne.0) then
         neigz(1) = ipe - 1
      else
         neigz(1) = ipe + npz - 1
      endif
      if (ipz.ne.npz-1) then
         neigz(2) = ipe + 1
      else
         neigz(2) = ipe - npz + 1
      endif
c
c     Compute dinv, dc
c
      do j = iys, iys+ny+1
         do i = ixs, ixs+nx+1
            do k = 0, nz+1
               w1(k,i,j) = (0.d0,0.d0)
               w2(k,i,j) = (0.d0,0.d0)
            enddo
         enddo
      enddo
cwas-s
c      write(6,*) "trzp-str" ,ipe
cwas-e
c
c      call trzp(a(0,1,1,15), a(0,1,1,15), nz, nallx2, nally2, 
c     &     ixs, iys, nx, ny, 
c     &     1, 1, nx, ny, 1,1, nx, ny, 
c     &     neigz, icomm, sb, rb)
cwas-s
c      write(6,*) "trzp-end" ,ipe
cwas-e
c
c     Phase 0
c
      iph0 = 0
      ipzh0 = mod(ipz - iph0 + npz , npz)
      ixys0 = (ipzh0*nxy)/npz+1
      ixye0 = ((ipzh0+1)*nxy)/npz
c
      iys0 = (ixys0-1)/nx + 1
      ixs0 = ixys0 - (iys0-1)*nx
      iye0 = (ixye0-1)/nx + 1
      ixe0 = ixye0 - (iye0-1)*nx
c
      iph1  = mod(iph0 + 1,npz)
      ipzh1 = mod(ipz - iph1 + npz, npz)
      ixys1 = (ipzh1*nxy)/npz+1
      ixye1 = ((ipzh1+1)*nxy)/npz
c
      iys1 = (ixys1-1)/nx + 1
      ixs1 = ixys1 - (iys1-1)*nx
      iye1 = (ixye1-1)/nx + 1
      ixe1 = ixye1 - (iye1-1)*nx
c
      do iy = iys0, iye0
         j = iys + iy
         if (iys0.eq.iye0) then
            ix0 = ixs0
            ix1 = ixe0
         elseif (iy.eq.iys0) then
            ix0 = ixs0
            ix1 = nx
         elseif (iy.eq.iye0) then
            ix0 = 1
            ix1 = ixe0
         else
            ix0 = 1
            ix1 = nx
         endif
         do i = ixs+ix0, ixs+ix1
            w1(1,i,j) = 1.d0/a(1,i,j,14)
            do k = 2, nz
               d1 = a(k,i,j,14) 
     &            - a(k,i,j,13)*w1(k-1,i,j)*a(k-1,i,j,15)
               if (abs(d1) .lt. 1.d-20) then
                  ierr = 10
                  d1 = (1.d-20,0.d0)
               endif
               w1(k,i,j) = 1.d0/d1
            enddo
         enddo
      enddo
c
c      call trzp(w1, w1, nz, nallx2, nally2, 
c     &     ixs, iys, nx, ny, 
c     &     ixs0, iys0, ixe0, iye0, ixs1, iys1, ixe1, iye1, 
c     &     neigz, icomm, sb, rb)
c
c     Phase 1 to npz-2
c
      do iph0 = 1, npz-2
         ipzh0 = mod(ipz - iph0 + npz , npz)
         ixys0 = (ipzh0*nxy)/npz+1
         ixye0 = ((ipzh0+1)*nxy)/npz
c
         iys0 = (ixys0-1)/nx + 1
         ixs0 = ixys0 - (iys0-1)*nx
         iye0 = (ixye0-1)/nx + 1
         ixe0 = ixye0 - (iye0-1)*nx
c
         iph1  = mod(iph0 + 1,npz)
         ipzh1 = mod(ipz - iph1 + npz, npz)
         ixys1 = (ipzh1*nxy)/npz+1
         ixye1 = ((ipzh1+1)*nxy)/npz
c
         iys1 = (ixys1-1)/nx + 1
         ixs1 = ixys1 - (iys1-1)*nx
         iye1 = (ixye1-1)/nx + 1
         ixe1 = ixye1 - (iye1-1)*nx
c
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               do k = 1,nz
                  d1 = a(k,i,j,14) 
     &               - a(k,i,j,13)*w1(k-1,i,j)*a(k-1,i,j,15)
                  if (abs(d1).lt.1.d-20) then
                     ierr = 11
                     d1 = (1.d-20,0.d0)
                  endif
                  w1(k,i,j) = 1.d0/d1
               enddo
            enddo
         enddo
c
c         call trzp(w1, w1, nz, nallx2, nally2, 
c     &        ixs, iys, nx, ny, 
c     &        ixs0, iys0, ixe0, iye0, ixs1, iys1, ixe1, iye1, 
c     &        neigz, icomm, sb, rb)
      enddo
c
c     Phase npz-1
c
      if (npz.ne.1) then
         iph0  = npz - 1
         ipzh0 = mod(ipz - iph0 + npz , npz)
         ixys0 = (ipzh0*nxy)/npz+1
         ixye0 = ((ipzh0+1)*nxy)/npz
c
         iys0 = (ixys0-1)/nx + 1
         ixs0 = ixys0 - (iys0-1)*nx
         iye0 = (ixye0-1)/nx + 1
         ixe0 = ixye0 - (iye0-1)*nx
c
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               do k = 1,nz-1
                  d1 = a(k,i,j,14) 
     &               - a(k,i,j,13)*w1(k-1,i,j)*a(k-1,i,j,15)
                  if (abs(d1).lt.1.d-20) then
                     ierr = 12
                     d1 = (1.d-20,0.d0)
                  endif
                  w1(k,i,j) = 1.d0/d1
               enddo
            enddo
         enddo
      endif
c
      do j = iys+1, iys+ny
         do i = ixs+1, ixs+nx
            do k = 1, nz
               dinv(k,i,j) = w1(k,i,j)
            enddo
         enddo
      enddo
c     
      do j = iys,iys+ny+1
         do i = ixs,ixs+nx+1
            do k = 0,nz+1
               w1(k,i,j) = (0.d0,0.d0)
               w2(k,i,j) = (0.d0,0.d0)
            enddo
         enddo
      enddo
c
c     Phase 0 (F)
c
      iph0 = 0
      ipzh0 = mod(ipz - iph0 + npz , npz)
      ixys0 = (ipzh0*nxy)/npz+1
      ixye0 = ((ipzh0+1)*nxy)/npz
c
      iys0 = (ixys0-1)/nx + 1
      ixs0 = ixys0 - (iys0-1)*nx
      iye0 = (ixye0-1)/nx + 1
      ixe0 = ixye0 - (iye0-1)*nx
c
      iph1  = mod(iph0 + 1,npz)
      ipzh1 = mod(ipz - iph1 + npz, npz)
      ixys1 = (ipzh1*nxy)/npz+1
      ixye1 = ((ipzh1+1)*nxy)/npz
c
      iys1 = (ixys1-1)/nx + 1
      ixs1 = ixys1 - (iys1-1)*nx
      iye1 = (ixye1-1)/nx + 1
      ixe1 = ixye1 - (iye1-1)*nx
      do iy = iys0, iye0
         j = iys + iy
         if (iys0.eq.iye0) then
            ix0 = ixs0
            ix1 = ixe0
         elseif (iy.eq.iys0) then
            ix0 = ixs0
            ix1 = nx
         elseif (iy.eq.iye0) then
            ix0 = 1
            ix1 = ixe0
         else
            ix0 = 1
            ix1 = nx
         endif
         do i = ixs+ix0, ixs+ix1
            w1(1,i,j) = dinv(1,i,j)*a(1,i,j,13)
            do k = 2, nz
               w1(k,i,j) = 
     &              -dinv(k,i,j)*a(k,i,j,13)*w1(k-1,i,j)
            enddo
         enddo
      enddo
c      call trzp(w1, w1, nz, nallx2, nally2, 
c     &     ixs, iys, nx, ny, 
c     &     ixs0, iys0, ixe0, iye0, ixs1, iys1, ixe1, iye1, 
c     &     neigz, icomm, sb, rb)
c
c     Phase 1 to npz-2 (F)
c
      do iph0 = 1, npz-2
         ipzh0 = mod(ipz - iph0 + npz , npz)
         ixys0 = (ipzh0*nxy)/npz+1
         ixye0 = ((ipzh0+1)*nxy)/npz
c
         iys0 = (ixys0-1)/nx + 1
         ixs0 = ixys0 - (iys0-1)*nx
         iye0 = (ixye0-1)/nx + 1
         ixe0 = ixye0 - (iye0-1)*nx
c
         iph1 = mod(iph0 + 1,npz)
         ipzh1 = mod(ipz - iph1 + npz, npz)
         ixys1 = (ipzh1*nxy)/npz+1
         ixye1 = ((ipzh1+1)*nxy)/npz
c
         iys1 = (ixys1-1)/nx + 1
         ixs1 = ixys1 - (iys1-1)*nx
         iye1 = (ixye1-1)/nx + 1
         ixe1 = ixye1 - (iye1-1)*nx
c
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               do k = 1, nz
                  w1(k,i,j) = 
     &                 -dinv(k,i,j)*a(k,i,j,13)*w1(k-1,i,j)
               enddo
            enddo
         enddo
c         call trzp(w1, w1, nz, nallx2, nally2, 
c     &        ixs, iys, nx, ny, 
c     &        ixs0, iys0, ixe0, iye0, ixs1, iys1, ixe1, iye1, 
c     &        neigz, icomm, sb, rb)
      enddo
c
c     Phase npz-1 (F)
c
      if (npz.ne.1) then
         iph0  = npz - 1
         ipzh0 = mod(ipz - iph0 + npz , npz)
         ixys0 = (ipzh0*nxy)/npz+1
         ixye0 = ((ipzh0+1)*nxy)/npz
c
         iys0 = (ixys0-1)/nx + 1
         ixs0 = ixys0 - (iys0-1)*nx
         iye0 = (ixye0-1)/nx + 1
         ixe0 = ixye0 - (iye0-1)*nx
c
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               do k = 1, nz-2
                  w1(k,i,j) = 
     &                 -dinv(k,i,j)*a(k,i,j,13)*w1(k-1,i,j)
               enddo
            enddo
         enddo
      endif
c
      iph0 = npz-1
      ipzh0 = mod(ipz - iph0 + npz , npz)
      ixys0 = (ipzh0*nxy)/npz+1
      ixye0 = ((ipzh0+1)*nxy)/npz
c
      iys0 = (ixys0-1)/nx + 1
      ixs0 =ixys0 - (iys0-1)*nx
      iye0 = (ixye0-1)/nx + 1
      ixe0 =ixye0 - (iye0-1)*nx
c
      if (nz.ne.1) then
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               w1(nz-1,i,j) = dinv(nz-1,i,j)*(a(nz-1,i,j,15)
     &              - a(nz-1,i,j,13)*w1(nz-2,i,j))
            enddo
         enddo
      else
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               w1(1,i,j) = (-1.d0,0.d0)
            enddo
         enddo
      endif
c
c     Phase npz-1 (B)
c
      iph0  = npz-1
      ipzh0 = mod(ipz - iph0 + npz , npz)
      ixys0 = (ipzh0*nxy)/npz+1
      ixye0 = ((ipzh0+1)*nxy)/npz
c
      iys0 = (ixys0-1)/nx + 1
      ixs0 = ixys0 - (iys0-1)*nx
      iye0 = (ixye0-1)/nx + 1
      ixe0 = ixye0 - (iye0-1)*nx
      do iy = iys0, iye0
         j = iys + iy
         if (iys0.eq.iye0) then
            ix0 = ixs0
            ix1 = ixe0
         elseif (iy.eq.iys0) then
            ix0 = ixs0
            ix1 = nx
         elseif (iy.eq.iye0) then
            ix0 = 1
            ix1 = ixe0
         else
            ix0 = 1
            ix1 = nx
         endif
         do i = ixs+ix0, ixs+ix1
            do k = nz-2,1,-1
               w1(k,i,j) = w1(k,i,j) 
     &              - dinv(k,i,j)*a(k,i,j,15)*w1(k+1,i,j)
            enddo
         enddo
      enddo
c
c     Phase npz-2 to 0 (B)
c
      do iph0 = npz-2, 0, -1
         ipzh0 = mod(ipz - iph0 + npz , npz)
         ixys0 = (ipzh0*nxy)/npz+1
         ixye0 = ((ipzh0+1)*nxy)/npz
c
         iys0 = (ixys0-1)/nx + 1
         ixs0 = ixys0 - (iys0-1)*nx
         iye0 = (ixye0-1)/nx + 1
         ixe0 = ixye0 - (iye0-1)*nx
c
         iph1 = mod(iph0 + 1,npz)
         ipzh1 = mod(ipz - iph1 + npz, npz)
         ixys1 = (ipzh1*nxy)/npz+1
         ixye1 = ((ipzh1+1)*nxy)/npz
c
         iys1 = (ixys1-1)/nx + 1
         ixs1 = ixys1 - (iys1-1)*nx
         iye1 = (ixye1-1)/nx + 1
         ixe1 = ixye1 - (iye1-1)*nx
c         call trzm(w1, w1, nz, nallx2, nally2, 
c     &        ixs, iys, nx, ny, 
c     &        ixs1, iys1, ixe1, iye1, ixs0, iys0, ixe0, iye0, 
c     &        neigz, icomm, sb, rb)
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               ix0 = ixs0
               ix1 = ixe0
            elseif (iy.eq.iys0) then
               ix0 = ixs0
               ix1 = nx
            elseif (iy.eq.iye0) then
               ix0 = 1
               ix1 = ixe0
            else
               ix0 = 1
               ix1 = nx
            endif
            do i = ixs+ix0, ixs+ix1
               do k = nz, 1, -1
                  w1(k,i,j) = w1(k,i,j) 
     &                 - dinv(k,i,j)*a(k,i,j,15)*w1(k+1,i,j)
               enddo
            enddo
         enddo
      enddo
      do j = iys+1, iys+ny
         do i = ixs+1, ixs+nx
            do k = 1, nz
               if (abs(w1(k,i,j)).gt.eps) then
                  dc(k,i,j) = w1(k,i,j)
               else
                  dc(k,i,j) = (0.d0,0.d0)
               endif
            enddo
         enddo
      enddo
c
      iph0 = npz-1
      ipzh0 = mod(ipz - iph0 + npz , npz)
      ixys0 = (ipzh0*nxy)/npz+1
      ixye0 = ((ipzh0+1)*nxy)/npz
c
      iys0 = (ixys0-1)/nx + 1
      ixs0 =ixys0 - (iys0-1)*nx
      iye0 = (ixye0-1)/nx + 1
      ixe0 =ixye0 - (iye0-1)*nx
c
      iph1 = mod(iph0 + 1,npz)
      ipzh1 = mod(ipz - iph1 + npz, npz)
      ixys1 = (ipzh1*nxy)/npz+1
      ixye1 = ((ipzh1+1)*nxy)/npz
c
      iys1 = (ixys1-1)/nx + 1
      ixs1 =ixys1 - (iys1-1)*nx
      iye1 = (ixye1-1)/nx + 1
      ixe1 =ixye1 - (iye1-1)*nx
c      call trzm(w1, w1, nz, nallx2, nally2, 
c     &     ixs, iys, nx, ny, 
c     &     ixs1, iys1, ixe1, iye1, ixs0, iys0, ixe0, iye0, 
c     &     neigz, icomm, sb, rb)
      iphm = mod(iph0 - 1 + npz,npz)
      ipzhm = mod(ipz - iphm + npz , npz)
      ixysm = (ipzhm*nxy)/npz+1
      ixyem = ((ipzhm+1)*nxy)/npz
c
      iysm = (ixysm-1)/nx + 1
      ixsm =ixysm - (iysm-1)*nx
      iyem = (ixyem-1)/nx + 1
      ixem =ixyem - (iyem-1)*nx
c      call trzp(w1, w1, nz, nallx2, nally2, 
c     &     ixs, iys, nx, ny, 
c     &     ixsm, iysm, ixem, iyem, ixs0, iys0, ixe0, iye0, 
c     &     neigz, icomm, sb, rb)
      do iy = iys0, iye0
         j = iys + iy
         if (iys0.eq.iye0) then
            ix0 = ixs0
            ix1 = ixe0
         elseif (iy.eq.iys0) then
            ix0 = ixs0
            ix1 = nx
         elseif (iy.eq.iye0) then
            ix0 = 1
            ix1 = ixe0
         else
            ix0 = 1
            ix1 = nx
         endif
         do i = ixs+ix0, ixs+ix1
            d1 = a(nz,i,j,14) 
     &           - a(nz,i,j,13)*w1(nz-1,i,j) 
     &           - a(nz,i,j,15)*w1(nz+1,i,j)
            if (abs(d1).lt.1.d-20) then
               ierr = 13
               d1 = (1.d-20,0.d0)
            endif
            dinv(nz,i,j) = 1.d0/d1
         enddo
      enddo
      return
      end
