      subroutine jac27msg(a,dinv,dc,w,f,work, 
     &                    nz,nallx2,nally2,
     &                    omega,nr, 
     &                    ixs,iys,nx,ny, 
     &                    ix0,iy0,
     &                    ipe,ipz,npz, 
     &                    neigz,neigx,neigy)
c
      implicit real*8 (a-h,o-z)
      complex*16  a(0:nz+1,nallx2,nally2,27) 
      complex*16  dinv(0:nz+1,nallx2,nally2),dc(0:nz+1,nallx2,nally2) 
      complex*16  w(0:nz+1,nallx2,nally2),f(0:nz+1,nallx2,nally2) 
      complex*16  work(0:nz+1,nallx2,nally2) 
      integer     neigz(2),neigx(2),neigy(2),neigzc(2) 
c
      nxy = nx*ny
      if (nxy.eq.0) return
      if (ipz.ne.0) then
         neigzc(1) = ipe - 1
      else
         neigzc(1) = ipe + npz - 1
      endif
      if (ipz.ne.npz-1) then
         neigzc(2) = ipe + 1
      else
         neigzc(2) = ipe - npz + 1
      endif
c     iteration
      do 100 ir = 1, nr
c
         do j = iys+1, iys+ny
            do i = ixs+1, ixs+nx
               do k = 1, nz
                  work(k,i,j) = f(k,i,j) 
     &                       - a(k,i,j,1 )*w(k-1,i-1,j-1)
     &                       - a(k,i,j,2 )*w(k  ,i-1,j-1)
     &                       - a(k,i,j,3 )*w(k+1,i-1,j-1)
     &                       - a(k,i,j,4 )*w(k-1,i  ,j-1)
     &                       - a(k,i,j,5 )*w(k  ,i  ,j-1)
     &                       - a(k,i,j,6 )*w(k+1,i  ,j-1)
     &                       - a(k,i,j,7 )*w(k-1,i+1,j-1)
     &                       - a(k,i,j,8 )*w(k  ,i+1,j-1)
     &                       - a(k,i,j,9 )*w(k+1,i+1,j-1)
     &                       - a(k,i,j,10)*w(k-1,i-1,j  )
     &                       - a(k,i,j,11)*w(k  ,i-1,j  )
     &                       - a(k,i,j,12)*w(k+1,i-1,j  )
     &                       - a(k,i,j,16)*w(k-1,i+1,j  )
     &                       - a(k,i,j,17)*w(k  ,i+1,j  )
     &                       - a(k,i,j,18)*w(k+1,i+1,j  )
     &                       - a(k,i,j,19)*w(k-1,i-1,j+1)
     &                       - a(k,i,j,20)*w(k  ,i-1,j+1)
     &                       - a(k,i,j,21)*w(k+1,i-1,j+1)
     &                       - a(k,i,j,22)*w(k-1,i  ,j+1)
     &                       - a(k,i,j,23)*w(k  ,i  ,j+1)
     &                       - a(k,i,j,24)*w(k+1,i  ,j+1)
     &                       - a(k,i,j,25)*w(k-1,i+1,j+1)
     &                       - a(k,i,j,26)*w(k  ,i+1,j+1)
     &                       - a(k,i,j,27)*w(k+1,i+1,j+1)
               enddo
            enddo
         enddo
c
c     
c     Phase 0 (F)
c
         iph0  = 0
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
               i0 = ixs0
               i1 = ixe0
            elseif (iy.eq.iys0) then
               i0 = ixs0
               i1 = nx
            elseif (iy.eq.iye0) then
               i0 = 1
               i1 = ixe0
            else
               i0 = 1
               i1 = nx
            endif
c
            do i = ixs+i0, ixs+i1
               work(1,i,j) = dinv(1,i,j)*work(1,i,j)
               if (npz.eq.1) then
                  do k = 2, nz-1
                     work(k,i,j) = dinv(k,i,j)*
     &                    (work(k,i,j)-a(k,i,j,13)*work(k-1,i,j))
                  enddo
               else
                  do k = 2, nz
                     work(k,i,j) = dinv(k,i,j)*
     &                    (work(k,i,j)-a(k,i,j,13)*work(k-1,i,j))
                  enddo
               endif
            enddo
         enddo
c
         if ((npz.eq.1).and.(nz.eq.1)) goto 500
c     
c     Phase 1 to npz-2 (F)
c
         do iph0 = 1, npz-2
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
c
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
c
               do i = ixs+i0, ixs+i1
                  do k = 1, nz
                     work(k,i,j) = dinv(k,i,j)*
     &                    (work(k,i,j)-a(k,i,j,13)*work(k-1,i,j))
                  enddo
               enddo
            enddo
         enddo
c
c     Phase npz-1 (F)
c
         if (npz.ne.1) then
            iph0 = npz - 1
            ipzh0 = mod(ipz - iph0 + npz , npz)
            ixys0 = (ipzh0*nxy)/npz+1
            ixye0 = ((ipzh0+1)*nxy)/npz
c
            iys0 = (ixys0-1)/nx + 1
            ixs0 =ixys0 - (iys0-1)*nx
            iye0 = (ixye0-1)/nx + 1
            ixe0 =ixye0 - (iye0-1)*nx
c     
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
               do i = ixs+i0, ixs+i1
                  do k = 1, nz-1
                     work(k,i,j) = dinv(k,i,j)*
     &                    (work(k,i,j)-a(k,i,j,13)*work(k-1,i,j))
                  enddo
               enddo
            enddo
         endif
c     
c     Phase npz-1 (B)
c     
         if (nz.ne.1) then
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
            iphm = mod(iph0 - 1 + npz,npz)
            ipzhm = mod(ipz - iphm + npz , npz)
            ixysm = (ipzhm*nxy)/npz+1
            ixyem = ((ipzhm+1)*nxy)/npz
c
            iysm = (ixysm-1)/nx + 1
            ixsm =ixysm - (iysm-1)*nx
            iyem = (ixyem-1)/nx + 1
            ixem =ixyem - (iyem-1)*nx
c
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
c
               do i = ixs+i0, ixs+i1
                  do k = nz-2, 1, -1
                     work(k,i,j) = work(k,i,j)
     &                    -dinv(k,i,j)*a(k,i,j,15)*work(k+1,i,j)
                  enddo
               enddo
            enddo
         elseif (npz.gt.1) then
            iph0 = npz-2
            ipzh0 = mod(ipz - iph0 + npz , npz)
            ixys0 = (ipzh0*nxy)/npz+1
            ixye0 = ((ipzh0+1)*nxy)/npz
c
            iys0 = (ixys0-1)/nx + 1
            ixs0 =ixys0 - (iys0-1)*nx
            iye0 = (ixye0-1)/nx + 1
            ixe0 =ixye0 - (iye0-1)*nx
c
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
c
               do i = ixs+i0, ixs+i1
                  work(nz+1,i,j) = 0.d0
               enddo
            enddo
         endif
c     
c     Phase npz-2 to 0 (B)
c     
         do iph0 = npz-2, 0 , -1
            ipzh0 = mod(ipz - iph0 + npz , npz)
            ixys0 = (ipzh0*nxy)/npz+1
            ixye0 = ((ipzh0+1)*nxy)/npz
c
            iys0 = (ixys0-1)/nx + 1
            ixs0 =ixys0 - (iys0-1)*nx
            iye0 = (ixye0-1)/nx + 1
            ixe0 =ixye0 - (iye0-1)*nx
c
            iphm = mod(iph0 - 1 + npz,npz)
            ipzhm = mod(ipz - iphm + npz , npz)
            ixysm = (ipzhm*nxy)/npz+1
            ixyem = ((ipzhm+1)*nxy)/npz
c
            iysm = (ixysm-1)/nx + 1
            ixsm =ixysm - (iysm-1)*nx
            iyem = (ixyem-1)/nx + 1
            ixem =ixyem - (iyem-1)*nx
c     
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
c
               do i = ixs+i0, ixs+i1
                  do k = nz, 1, -1
                     work(k,i,j) = work(k,i,j)
     &                    -dinv(k,i,j)*a(k,i,j,15)*work(k+1,i,j)
                  enddo
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
         iphm = mod(iph0 - 1 + npz,npz)
         ipzhm = mod(ipz - iphm + npz , npz)
         ixysm = (ipzhm*nxy)/npz+1
         ixyem = ((ipzhm+1)*nxy)/npz
c     
         iysm = (ixysm-1)/nx + 1
         ixsm =ixysm - (iysm-1)*nx
         iyem = (ixyem-1)/nx + 1
         ixem =ixyem - (iyem-1)*nx
c     
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               i0 = ixs0
               i1 = ixe0
            elseif (iy.eq.iys0) then
               i0 = ixs0
               i1 = nx
            elseif (iy.eq.iye0) then
               i0 = 1
               i1 = ixe0
            else
               i0 = 1
               i1 = nx
            endif
c     
            do i = ixs+i0, ixs+i1
               work(nz,i,j) = dinv(nz,i,j)*
     &              (work(nz,i,j)-a(nz,i,j,13)*work(nz-1,i,j)
     &              -a(nz,i,j,15)*work(nz+1,i,j))
               work(nz+1,i,j) = work(nz,i,j)
            enddo
         enddo
c     
c     Propagation of work_n
c     
         do iph0 = 0, npz-2
            ipzh0 = mod(ipz - iph0 + npz , npz)
            ixys0 = (ipzh0*nxy)/npz+1
            ixye0 = ((ipzh0+1)*nxy)/npz
c
            iys0 = (ixys0-1)/nx + 1
            ixs0 =ixys0 - (iys0-1)*nx
            iye0 = (ixye0-1)/nx + 1
            ixe0 =ixye0 - (iye0-1)*nx
c     
            iphm = mod(iph0 - 1 + npz,npz)
            ipzhm = mod(ipz - iphm + npz , npz)
            ixysm = (ipzhm*nxy)/npz+1
            ixyem = ((ipzhm+1)*nxy)/npz
c     
            iysm = (ixysm-1)/nx + 1
            ixsm =ixysm - (iysm-1)*nx
            iyem = (ixyem-1)/nx + 1
            ixem =ixyem - (iyem-1)*nx
         enddo
c     
c     Correction by dc
c
         do iph0 = 0, npz-2
            ipzh0 = mod(ipz - iph0 + npz , npz)
            ixys0 = (ipzh0*nxy)/npz+1
            ixye0 = ((ipzh0+1)*nxy)/npz
c
            iys0 = (ixys0-1)/nx + 1
            ixs0 =ixys0 - (iys0-1)*nx
            iye0 = (ixye0-1)/nx + 1
            ixe0 =ixye0 - (iye0-1)*nx
c     
            do iy = iys0, iye0
               j = iys + iy
               if (iys0.eq.iye0) then
                  i0 = ixs0
                  i1 = ixe0
               elseif (iy.eq.iys0) then
                  i0 = ixs0
                  i1 = nx
               elseif (iy.eq.iye0) then
                  i0 = 1
                  i1 = ixe0
               else
                  i0 = 1
                  i1 = nx
               endif
c     
               do i = ixs+i0, ixs+i1
                  do k = 1, nz
                     work(k,i,j) = work(k,i,j) 
     &                    - work(nz+1,i,j)*dc(k,i,j)
                  enddo
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
         do iy = iys0, iye0
            j = iys + iy
            if (iys0.eq.iye0) then
               i0 = ixs0
               i1 = ixe0
            elseif (iy.eq.iys0) then
               i0 = ixs0
               i1 = nx
            elseif (iy.eq.iye0) then
               i0 = 1
               i1 = ixe0
            else
               i0 = 1
               i1 = nx
            endif
c     
            do i = ixs+i0, ixs+i1
               do k = 1, nz-1
                  work(k,i,j) = work(k,i,j) 
     &                 - work(nz+1,i,j)*dc(k,i,j)
               enddo
            enddo
         enddo
c     
c     Update
c
 500     continue
         do j = iys+1, iys+ny
            do i = ixs+1, ixs+nx
               do k = 1, nz
                  w(k,i,j) = (1.d0 - omega)*w(k,i,j)
     &                    + omega*work(k,i,j)
               enddo
            enddo
         enddo
 100  continue
c
      return
      end
