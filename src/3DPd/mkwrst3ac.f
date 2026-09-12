      subroutine mkwrst3ac(wrst,witp,
     &                     nallx2,nally2, 
     &                     nx,nxf,ixs,ixsf,ix0f, 
     &                     ny,nyf,iys,iysf,iy0f, 
     &                     nz,
     &                     neigx, neigy)
c
c     This subroutine makes weights of the restriction mapping from
c     G(kx+1,ky+1) to G(kx,ky). Here it is defined as the transpose 
c     of the prolongation mapping.
c
      implicit real*8 (a-h,o-z)
      complex*16 czero,cone,chalf 
      parameter (eps = 1.d-8, izero = 0, ione = 1,epsvs = 1.d-20, 
     &           czero = (0.d0,0.d0),cone = (1.d0,0.d0),
     &           chalf = (0.5d0,0.d0))
      complex*16  wrst(0:nz+1,nallx2,nally2,-1:1,-1:1)  
      complex*16  witp(0:nz+1,nallx2,nally2,4) 
      complex*16  d 
c
      if (nxf.le.0 .or. nyf.le.0) return
c
c     initialize
c
      nx1 = nxf
      ny1 = nyf
c
      do iw = -1,1
         do jw = -1,1
            do j = iys+1, iys+ny
               do i = ixs+1, ixs+nx
                  do k = 1, nz
                     wrst(k,i,j,iw,jw) = czero
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      do j = iysf, iysf+nyf+1
         do k = 0, nz+1
            witp(k,ixsf,j,1)       = czero
            witp(k,ixsf,j,2)       = czero
            witp(k,ixsf,j,3)       = czero
            witp(k,ixsf,j,4)       = czero
            witp(k,ixsf+nxf+1,j,1) = czero
            witp(k,ixsf+nxf+1,j,2) = czero
            witp(k,ixsf+nxf+1,j,3) = czero
            witp(k,ixsf+nxf+1,j,4) = czero
         enddo
      enddo
      do i = ixsf, ixsf+nxf+1
         do k = 0, nz+1
            witp(k,i,iysf,1)       = czero
            witp(k,i,iysf,2)       = czero
            witp(k,i,iysf,3)       = czero
            witp(k,i,iysf,4)       = czero
            witp(k,i,iysf+nyf+1,1) = czero
            witp(k,i,iysf+nyf+1,2) = czero
            witp(k,i,iysf+nyf+1,3) = czero
            witp(k,i,iysf+nyf+1,4) = czero
         enddo
      enddo

c
      irstr = 0
      if (irstr .eq. 0) then
c
c        (-1,-1)-component
c
         do iy = 1,ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1,nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1,nz
                  wrst(k,i,j,-1,-1) = witp(k,if-1,jf-1,4)
c                  write(*,*) wrst(k,i,j,-1,-1)
               enddo
            enddo
         enddo
c        ( 0,-1)-component
         do iy = 1,ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1,nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1,nz
                  wrst(k,i,j,0,-1) = witp(k,if,jf-1,2)
               enddo
            enddo
         enddo
c        ( 1,-1)-component
         do iy = 1,ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1,nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1,nz
                  wrst(k,i,j,1,-1) = witp(k,if+1,jf-1,3)
               enddo
            enddo
         enddo
c        (-1, 0)-component
         do iy = 1,ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1,nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1,nz
                  wrst(k,i,j,-1, 0) = witp(k,if-1,jf,2)
               enddo
            enddo
         enddo
c        ( 0, 0)-component
         do iy = 1, ny
            j = iys+iy
            do ix = 1, nx
               i = ixs+ix
               do k = 1, nz
                  wrst(k,i,j,0,0) = cone
               enddo
            enddo
         enddo
c        ( 1, 0)-component
         do iy = 1, ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1, nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1, nz
                  wrst(k,i,j, 1, 0) = witp(k,if+1,jf,1)
               enddo
            enddo
         enddo
c        (-1, 1)-component
         do iy = 1, ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1, nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1, nz
                  wrst(k,i,j,-1, 1) = witp(k,if-1,jf+1,2)
               enddo
            enddo
         enddo
c        ( 0, 1)-component
         do iy = 1, ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1, nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1, nz
                  wrst(k,i,j, 0, 1) = witp(k,if,jf+1,1)
               enddo
            enddo
         enddo
c        ( 1, 1)-component
         do iy = 1, ny
            j = iys+iy
            jf = iysf+iy0f+iy*2-2
            do ix = 1, nx
               i = ixs+ix
               if = ixsf+ix0f+ix*2-2
               do k = 1, nz
                  wrst(k,i,j, 1, 1) = witp(k,if+1,jf+1,1)
               enddo
            enddo
         enddo
c
      elseif (irstr .eq. 1) then
c
         do iy = 1,ny
            j = iys+iy
            do ix = 1,nx
               i = ixs+ix
               do k = 1,nz
                  wrst(k,i,j,-1,-1) = czero
                  wrst(k,i,j, 0,-1) = czero
                  wrst(k,i,j, 1,-1) = czero
                  wrst(k,i,j,-1, 0) = czero
                  wrst(k,i,j, 0, 0) = 0.5d0*cone
                  wrst(k,i,j, 1, 0) = czero
                  wrst(k,i,j,-1, 1) = czero
                  wrst(k,i,j, 0, 1) = czero
                  wrst(k,i,j, 1, 1) = czero
               enddo
            enddo
         enddo
c
      elseif (irstr .eq. 2) then
c
         d = 4.d0/16.d0*cone
         do iy = 1,ny
            j = iys+iy
            if (iy .eq. 1) then
               do ix = 1,nx
                  i = ixs+ix
                  if (ix .eq. 1) then
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = 0.d0
                        wrst(k,i,j, 0,-1) = 0.d0
                        wrst(k,i,j, 1,-1) = 0.d0
                        wrst(k,i,j,-1, 0) = 0.d0
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = 0.d0
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = d
                     enddo
                  elseif (ix .eq. nx) then
                     do k=1,nz
                        wrst(k,i,j,-1,-1) = 0.d0
                        wrst(k,i,j, 0,-1) = 0.d0
                        wrst(k,i,j, 1,-1) = 0.d0
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 0.d0
                        wrst(k,i,j,-1, 1) = d
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = 0.d0
                     enddo
                  else
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = 0.d0
                        wrst(k,i,j, 0,-1) = 0.d0
                        wrst(k,i,j, 1,-1) = 0.d0
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = d
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = d
                     enddo
                  endif
               enddo
            elseif (iy .eq. ny) then
               do ix = 1,nx
                  i = ixs+ix
                  if (ix .eq. 1) then
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = 0.d0
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = d
                        wrst(k,i,j,-1, 0) = 0.d0
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = 0.d0
                        wrst(k,i,j, 0, 1) = 0.d0
                        wrst(k,i,j, 1, 1) = 0.d0
                     enddo
                  elseif (ix .eq. nx) then
                     do k =1,nz
                        wrst(k,i,j,-1,-1) = d
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = 0.d0
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 0.d0
                        wrst(k,i,j,-1, 1) = 0.d0
                        wrst(k,i,j, 0, 1) = 0.d0
                        wrst(k,i,j, 1, 1) = 0.d0
                     enddo
                  else
                     do k =1,nz
                        wrst(k,i,j,-1,-1) = d
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = d
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = 0.d0
                        wrst(k,i,j, 0, 1) = 0.d0
                        wrst(k,i,j, 1, 1) = 0.d0
                     enddo
                  endif
               enddo
            else
              do ix = 1,nx
                  i = ixs+ix
                  if (ix .eq. 1) then
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = d
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = d
                        wrst(k,i,j,-1, 0) = 0.d0
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = d
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = d
                     enddo
                  elseif (ix .eq. nx) then
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = d
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = d
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 0.d0
                        wrst(k,i,j,-1, 1) = d
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = d
                     enddo
                  else 
                     do k = 1,nz
                        wrst(k,i,j,-1,-1) = d
                        wrst(k,i,j, 0,-1) = 2.d0*d
                        wrst(k,i,j, 1,-1) = d
                        wrst(k,i,j,-1, 0) = 2.d0*d
                        wrst(k,i,j, 0, 0) = 4.d0*d
                        wrst(k,i,j, 1, 0) = 2.d0*d
                        wrst(k,i,j,-1, 1) = d
                        wrst(k,i,j, 0, 1) = 2.d0*d
                        wrst(k,i,j, 1, 1) = d
                     enddo
                  endif
               enddo
            endif
         enddo 
c
      elseif (irstr .eq. 3) then
c
         d = 1.d0/8.d0*cone
         do iy = 1,ny
            j = iys+iy
            do ix = 1,nx
               i = ixs+ix
               do k = 1,nz
                  wrst(k,i,j,-1,-1) = 0.d0
                  wrst(k,i,j, 0,-1) = d
                  wrst(k,i,j, 1,-1) = 0.d0
                  wrst(k,i,j,-1, 0) = d
                  wrst(k,i,j, 0, 0) = 4.d0*d
                  wrst(k,i,j, 1, 0) = d
                  wrst(k,i,j,-1, 1) = 0.d0
                  wrst(k,i,j, 0, 1) = d
                  wrst(k,i,j, 1, 1) = 0.d0
               enddo
            enddo
         enddo
c
      endif  
c
      do iy = 1,ny
         j = iys+iy
         do ix = 1,nx
            i = ixs+ix
            do k = 1,nz
               wrst(k,i,j,-1,-1) = conjg(wrst(k,i,j,-1,-1))
               wrst(k,i,j, 0,-1) = conjg(wrst(k,i,j, 0,-1))
               wrst(k,i,j, 1,-1) = conjg(wrst(k,i,j, 1,-1))
               wrst(k,i,j,-1, 0) = conjg(wrst(k,i,j,-1, 0))
               wrst(k,i,j, 0, 0) = conjg(wrst(k,i,j, 0, 0))
               wrst(k,i,j, 1, 0) = conjg(wrst(k,i,j, 1, 0))
               wrst(k,i,j,-1, 1) = conjg(wrst(k,i,j,-1, 1))
               wrst(k,i,j, 0, 1) = conjg(wrst(k,i,j, 1, 0))
               wrst(k,i,j, 1, 1) = conjg(wrst(k,i,j, 1, 1))
c               write(*,*) wrst(k,i,j,-1,-1)
c               write(*,*) wrst(k,i,j, 0,-1)
c               write(*,*) wrst(k,i,j, 1,-1)
c               write(*,*) wrst(k,i,j,-1, 0)
c               write(*,*) wrst(k,i,j, 0, 0)
c               write(*,*) wrst(k,i,j, 1, 0)
c               write(*,*) wrst(k,i,j,-1, 1)
c               write(*,*) wrst(k,i,j, 0, 1)
c               write(*,*) wrst(k,i,j, 1, 1)
            enddo
         enddo
      enddo
c
      return
      end
