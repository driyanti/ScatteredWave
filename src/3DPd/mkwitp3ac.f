      subroutine mkwitp3ac(witp,
     &                     nallx2,nally2,
     &                     nxf,ixsf,ixs0,
     &                     nyf,iysf,iys0,
     &                     nz, 
     &                     a,ids, 
     &                     work1,a2,a2t,
     &                     neigx,neigy)
c
c     This subroutine makes the weights of the prolongation mapping
c     from G(kx,ky) to G(k+1,ky+1). wobst(i,j) indicates where are 
c     obstacles. If wobst(i,j) gt 0.5, then (i,j) is included in an 
c     obstacle.
c
      implicit real*8 (a-h,o-z)
      complex*16  cone,czero,chalf 
      parameter  (eps = 1.d-8,izero = 0,ione = 1,epsvs = 1.d-20, 
     &            czero = (0.d0,0.d0),cone = (1.d0,0.d0),
     &            chalf = (0.5d0,0.d0), done = 1.d0, dzero = 0.d0, 
     &            dhalf = 0.5d0)
      complex*16  witp(0:nz+1,nallx2,nally2,4) 
      complex*16  a(0:nz+1,nallx2,nally2,-1:1,-1:1,-1:1)  
c
      complex*16  a2(nz,nxf,nyf,-1:1,-1:1),a2t(nz,nxf,nyf,-1:1,-1:1)  
      complex*16  work1(0:nz+1,0:nxf+1,0:nyf+1) 
      complex*16  c1,c2,csum,d,d1,d2,d3,d4,w1w,w1e,w1s,w2s,wln 
      real*8      da,db,dc,dw,de,dn,ds,w1,w2 
c
      integer     neigx(2),neigy(2) 
c
      nx1 = nxf
      ny1 = nyf
c     
c     initialize
c
      do j =  iysf, iysf+nyf+1
         do i = ixsf, ixsf+nxf+1
            do k = 1, nz
               witp(k,i,j,1) = (0.d0,0.d0)
               witp(k,i,j,2) = (0.d0,0.d0)
               witp(k,i,j,3) = (0.d0,0.d0)
               witp(k,i,j,4) = (0.d0,0.d0)
            enddo
         enddo
      enddo
c
      ixs1 = mod(ixs0,2)+1
      iys1 = iys0
      ixs2 = ixs0
      iys2 = mod(iys0,2)+1
      ixs3 = mod(ixs0,2)+1
      iys3 = mod(iys0,2)+1
c
      if (ids.eq.0) then
c
c     on grid points 1
c
         do j = iysf+iys1,iysf+nyf,2
            do i = ixsf+ixs1,ixsf+nxf,2
               do k = 1, nz
                  d1 = czero
                  d2 = czero
                  if (i .ne. ixsf+1)   d1 = cone
                  if (i .ne. ixsf+nxf) d2 = cone
                  d = d1 + d2
                  if (abs(d) .gt. dhalf) then
                     witp(k,i,j,1) = d1/d
                     witp(k,i,j,2) = d2/d
                  endif
                  witp(k,i,j,3) = czero
                  witp(k,i,j,4) = czero
               enddo
            enddo
         enddo
c
c     on grid points 2
c
         do j = iysf+iys2, iysf+nyf, 2
            do i = ixsf+ixs2, ixsf+nxf, 2
               do k = 1, nz
                  d1 = czero
                  d2 = czero
                  if (j .ne. iysf+1)   d1 = cone
                  if (j .ne. iysf+nyf) d2 = cone
                  d = d1 + d2
                  if (abs(d) .gt. dhalf) then
                     witp(k,i,j,1) = d1/d
                     witp(k,i,j,2) = d2/d
                  endif
                  witp(k,i,j,3) = czero
                  witp(k,i,j,4) = czero
               enddo
            enddo
         enddo
c
c     on grid points 3
c
         do j = iysf+iys3, iysf+nyf, 2
            do i = ixsf+ixs3, ixsf+nxf, 2
               do k = 1, nz
                  d1 = czero
                  d2 = czero
                  d3 = czero
                  d4 = czero
                  if (i .ne. ixsf+1   .and. j .ne. iysf+1  ) d1 = cone
                  if (i .ne. ixsf+nxf .and. j .ne. iysf+1  ) d2 = cone 
                  if (i .ne. ixsf+1   .and. j .ne. iysf+nyf) d3 = cone
                  if (i .ne. ixsf+nxf .and. j .ne. iysf+nyf) d4 = cone
                  d = d1 + d2 + d3 + d4
                  if (abs(d) .gt. dhalf) then
                     witp(k,i,j,1) = d1/d
                     witp(k,i,j,2) = d2/d
                     witp(k,i,j,3) = d3/d
                     witp(k,i,j,4) = d4/d
                  else
                     witp(k,i,j,1) = czero
                     witp(k,i,j,2) = czero
                     witp(k,i,j,3) = czero
                     witp(k,i,j,4) = czero
                  endif
               enddo
            enddo
         enddo
c
      elseif (ids.eq.1) then
c
c     Matrix dependent
c
c     lumping to 2d
c         
         do i0 = -1,1
            do j0 = -1,1
               do iy1 = 1,nyf
                  do ix1 = 1,nxf
                     j = iysf+iy1
                     i = ixsf+ix1
                     do k = 1,nz
                        a2(k,ix1,iy1,i0,j0) = a(k,i,j,-1,i0,j0)
     &                                      + a(k,i,j, 0,i0,j0)
     &                                      + a(k,i,j, 1,i0,j0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
c     on grid points 1
c
         do iy1 = iys1, nyf, 2
            do ix1 = ixs1, nxf, 2
               j = iysf+iy1
               i = ixsf+ix1
               do k = 1,nz
                  d1 = a2(k,ix1,iy1,-1,-1)
     &               + a2(k,ix1,iy1,-1,0)
     &               + a2(k,ix1,iy1,-1,1)
                  d2 = a2(k,ix1,iy1,1,-1) 
     &               + a2(k,ix1,iy1,1,0)
     &               + a2(k,ix1,iy1,1,1)
                  if (i .eq. ixsf+1)   d1 = czero
                  if (i .eq. ixsf+nxf) d2 = czero
                  d = -(  a2(k,ix1,iy1, 0, 0)
     &                  + a2(k,ix1,iy1, 0,-1)
     &                  + a2(k,ix1,iy1, 0, 1))
                  iw = 0
                  if (abs(d) .gt. epsvs) then 
                     w1 = d1/d
                     w2 = d2/d
                     if (abs(w1).ge.0.d0 .and. abs(w2).ge.0.d0) then
                        witp(k,i,j,1) = w1
                        witp(k,i,j,2) = w2
                        iw = 1
                     endif
                  endif
                  if (iw .eq. 0) then
                     do i0 = -1,1,2
                        do j0 = -1, 1
                          if (abs(a2(k,ix1,iy1,i0,j0)/a2(k,ix1,iy1,0,0)) 
     &                         .gt. 0.d0) then
                               a2(k,ix1,iy1,0,0) = a2(k,ix1,iy1,0,0)
     &                                           + a2(k,ix1,iy1,i0,j0)
                               a2(k,ix1,iy1,i0,j0) = (0.d0,0.d0)
                          endif
                        enddo
                     enddo
                     d1 = a2(k,ix1,iy1,-1,-1)
     &                  + a2(k,ix1,iy1,-1,0)
     &                  + a2(k,ix1,iy1,-1,1)
                     d2 = a2(k,ix1,iy1,1,-1)
     &                  + a2(k,ix1,iy1,1,0)
     &                  + a2(k,ix1,iy1,1,1)
                     if (i .eq. ixsf+1)   d1 = czero
                     if (i .eq. ixsf+nxf) d2 = czero
                     d = -(  a2(k,ix1,iy1,0,0)+a2(k,ix1,iy1, 0,-1)
     &                    +a2(k,ix1,iy1,0,1))
                     if (abs(d) .gt. epsvs) then
                        witp(k,i,j,1) = d1/d
                        witp(k,i,j,2) = d2/d
                     endif
                  endif
               enddo
            enddo
         enddo
c
c     on grid points 2
c
         do iy1 = iys2,nyf,2
            do ix1 = ixs2,nxf,2
               j = iysf+iy1
               i = ixsf+ix1
               do k = 1, nz
                  d1 = a2(k,ix1,iy1,-1,-1)
     &               + a2(k,ix1,iy1, 0,-1)
     &               + a2(k,ix1,iy1, 1,-1)
                  d2 = a2(k,ix1,iy1,-1,1)+a2(k,ix1,iy1,0,1)
     &               +a2(k,ix1,iy1,1,1)
                  if (j .eq. iysf+1)   d1 = czero
                  if (j .eq. iysf+nyf) d2 = czero
                  d  = -(a2(k,ix1,iy1,0,0)+a2(k,ix1,iy1,-1, 0)
     &               + a2(k,ix1,iy1,1,0))
                  if (abs(d).gt.epsvs) then
                     w1 = d1/d
                     w2 = d2/d
                     if (abs(w1) .ge. 0.d0 .and. abs(w2) .ge. 0.d0) then
                        witp(k,i,j,1) = w1
                        witp(k,i,j,2) = w2
                        iw = 1
                     endif
                  endif
                  if (iw .eq. 0) then
                     do i0 = -1,1
                        do j0 = -1,1,2
                          if (abs(a2(k,ix1,iy1,i0,j0)/a2(k,ix1,iy1,0,0)) 
     &                         .gt. 0.d0) then
                               a2(k,ix1,iy1,0,0) = a2(k,ix1,iy1,0,0)
     &                                           + a2(k,ix1,iy1,i0,j0)
                               a2(k,ix1,iy1,i0,j0) = (0.d0,0.d0)
                          endif
                        enddo
                     enddo
                     d1 = a2(k,ix1,iy1,-1,-1)
     &                  + a2(k,ix1,iy1, 0,-1)
     &                  + a2(k,ix1,iy1, 1,-1)
                     d2 = a2(k,ix1,iy1,-1, 1)
     &                  + a2(k,ix1,iy1, 0, 1)
     &                  + a2(k,ix1,iy1, 1, 1)
                     if (j .eq. iysf+1)   d1 = czero
                     if (j .eq. iysf+nyf) d2 = czero 
                     d  = -(  a2(k,ix1,iy1, 0, 0)
     &                      + a2(k,ix1,iy1,-1, 0)
     &                      + a2(k,ix1,iy1, 1, 0))
                     if (abs(d).gt.epsvs) then
                        witp(k,i,j,1) = d1/d
                        witp(k,i,j,2) = d2/d
                     endif
                  endif
               enddo
            enddo
         enddo
c
c     on grid points 3
c
         do iy1 = iys3,nyf,2
            do ix1 = ixs3,nxf,2
               j = iysf+iy1
               i = ixsf+ix1
               do k = 1,nz
                  d1 = a2(k,ix1,iy1,-1,-1)
     &               + a2(k,ix1,iy1,0,-1)*witp(k,i,j-1,1)
     &               + a2(k,ix1,iy1,-1,0)*witp(k,i-1,j,1)
                  d2 = a2(k,ix1,iy1,1,-1)
     &               + a2(k,ix1,iy1,0,-1)*witp(k,i,j-1,2)
     &               + a2(k,ix1,iy1,1,0)*witp(k,i+1,j,1)
                  d3 = a2(k,ix1,iy1,-1,1)
     &               + a2(k,ix1,iy1,0,1)*witp(k,i,j+1,1)
     &               + a2(k,ix1,iy1,-1,0)*witp(k,i-1,j,2)
                  d4 = a2(k,ix1,iy1,1,1)
     &               + a2(k,ix1,iy1,0,1)*witp(k,i,j+1,2)
     &               + a2(k,ix1,iy1,1,0)*witp(k,i+1,j,2)
                  if (i .eq. ixsf+1   .and. j .eq. iysf+1  ) d1 = czero
                  if (i .eq. ixsf+nxf .and. j .eq. iysf+1  ) d2 = czero 
                  if (i .eq. ixsf+1   .and. j .eq. iysf+nyf) d3 = czero
                  if (i .eq. ixsf+nxf .and. j .eq. iysf+nyf) d4 = czero
c
                  d = -a2(k,ix1,iy1,0,0)
                  iw = 0
                  if (abs(d) .gt. epsvs) then
                     witp(k,i,j,1) = d1/d
                     witp(k,i,j,2) = d2/d
                     witp(k,i,j,3) = d3/d
                     witp(k,i,j,4) = d4/d
                  endif
               enddo
            enddo
         enddo
c
      elseif (ids.eq.2 .or. ids.eq.3) then
c
c     lumping to 2d
c         
         do i0 = -1,1
            do j0 = -1,1
               do iy1 = 1,nyf
                  do ix1 = 1,nxf
                     j = iysf+iy1
                     i = ixsf+ix1
                     do k = 1,nz
                        a2(k,ix1,iy1,i0,j0) = a(k,i,j,-1,i0,j0)
     &                                      + a(k,i,j,0,i0,j0)
     &                                      + a(k,i,j,1,i0,j0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
         nx1 = nxf
         ny1 = nyf
c
         do j0 = -1,1
            do i0 = -1,1
               do iy1 = 0, ny1+1
                  do ix1 = 0, nx1+1
                     do k = 0, nz+1
                        work1(k,ix1,iy1) = czero
                     enddo
                  enddo
               enddo
               do iy1 = 1, ny1
                  do ix1 = 1, nx1
                     do k = 1, nz
                        work1(k,ix1,iy1) = a2(k,ix1,iy1,-i0,-j0)
                     enddo
                  enddo
               enddo
c
               do iy1 = 1, ny1
                  do ix1 = 1, nx1
                     do iz = 1, nz
                        a2t(iz,ix1,iy1,i0,j0) = work1(iz,ix1+i0,iy1+j0)
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
         do l = -1,1
            do k = -1,1
               do iy1 = 1, ny1
                  do ix1 = 1, nx1
                     do iz = 1, nz
                        i = ixsf+ix1
                        j = iysf+iy1
                        a2t(iz,ix1,iy1,k,l) = 0.5d0*(a2t(iz,ix1,iy1,k,l) 
     &                                      + a2(iz,ix1,iy1,k,l))
                     enddo
                  enddo
               enddo
            enddo
         enddo
c
c     on grid points 1
c
         do iy1 = iys1, ny1, 2
            do ix1 = ixs1, nx1, 2
               do iz = 1, nz
c
                  i = ixsf+ix1
                  j = iysf+iy1
                  da = abs(  a2t(iz,ix1,iy1,-1,-1)
     &                     + a2t(iz,ix1,iy1,-1, 0)
     &                     + a2t(iz,ix1,iy1,-1, 1))
                  db = abs(a2t(iz,ix1,iy1,-1,-1))
                  dc = abs(a2t(iz,ix1,iy1,-1, 1))
                  dw = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1, 1,-1) 
     &                     + a2t(iz,ix1,iy1, 1, 0)
     &                     + a2t(iz,ix1,iy1, 1, 1))
                  db = abs(a2t(iz,ix1,iy1, 1,-1))
                  dc = abs(a2t(iz,ix1,iy1, 1, 1))
                  de = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1,-1, 1) 
     &                     + a2t(iz,ix1,iy1, 0, 1) 
     &                     + a2t(iz,ix1,iy1, 1, 1))
                  db = abs(a2t(iz,ix1,iy1,-1, 1))
                  dc = abs(a2t(iz,ix1,iy1, 1, 1))
                  dn = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1,-1,-1) 
     &                     + a2t(iz,ix1,iy1, 0,-1) 
     &                     + a2t(iz,ix1,iy1, 1,-1))
                  db = abs(a2t(iz,ix1,iy1,-1,-1))
                  dc = abs(a2t(iz,ix1,iy1, 1,-1))
                  ds = max(da,db,dc)
c
                  csum = (0.d0,0.d0)
                  do l = -1,1
                     do k = -1,1
                        csum = csum + a2(iz,ix1,iy1,k,l) ! modified abs
                     enddo
                  enddo
                  sigm = 0.5d0*
     &                   min(1.d0, abs(1.d0-csum/a2(iz,ix1,iy1,0,0)))
                  c1 = a2 (iz,ix1,iy1, 1,-1) + a2 (iz,ix1,iy1, 1, 0) + 
     &                 a2 (iz,ix1,iy1, 1, 1) - a2t(iz,ix1,iy1, 1,-1) -
     &                 a2t(iz,ix1,iy1, 1, 0) - a2t(iz,ix1,iy1, 1, 1) -
     &                 a2 (iz,ix1,iy1,-1,-1) - a2 (iz,ix1,iy1,-1, 0) -
     &                 a2 (iz,ix1,iy1,-1, 1) + a2t(iz,ix1,iy1,-1,-1) +
     &                 a2t(iz,ix1,iy1,-1, 0) + a2t(iz,ix1,iy1,-1, 1)
c
                  if (dw+de .gt. epsvs) then
                     w1w = sigm*(1.d0+(dw-de)/(dw+de)+c1/(dw+de+ds+dn))
                     w1e = 2.d0*sigm - w1w
c
                     w1w = max(0.d0,abs(w1w))*cone ! wlw = max(0.d0,abs(w1w))
                     w1e = max(0.d0,abs(w1e))*cone ! wlw = max(0.d0,abs(w1e))
                     w1  = min(2.d0*sigm,abs(w1w)) ! w1 = min(2.d0*sigm,w1w)
                     w2  = min(2.d0*sigm,abs(w1e)) ! w2 = min(2.d0*sigm,w1e)
c
c    Obstacle
c
                     if (i .eq. ixsf+1)   w1 = czero
                     if (i .eq. ixsf+nxf) w2 = czero
                     witp(iz,i,j,1) = w1*cone
                     witp(iz,i,j,2) = w2*cone
                  else
                     witp(iz,i,j,1) = czero
                     witp(iz,i,j,2) = czero
                  endif
               enddo
            enddo
         enddo
c
c     on grid points 2
c
         do iy1 = iys2, ny1, 2
            do ix1 = ixs2, nx1, 2
               do iz = 1, nz
c
                  i = ixsf+ix1
                  j = iysf+iy1
c
                  da = abs(  a2t(iz,ix1,iy1,-1,-1) 
     &                     + a2t(iz,ix1,iy1,-1, 0) 
     &                     + a2t(iz,ix1,iy1,-1, 1))
                  db = abs(a2t(iz,ix1,iy1,-1,-1))
                  dc = abs(a2t(iz,ix1,iy1,-1, 1))
                  dw = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1, 1,-1) 
     &                     + a2t(iz,ix1,iy1, 1, 0)
     &                     + a2t(iz,ix1,iy1, 1, 1))
                  db = abs(a2t(iz,ix1,iy1, 1,-1))
                  dc = abs(a2t(iz,ix1,iy1, 1, 1))
                  de = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1,-1, 1) 
     &                     + a2t(iz,ix1,iy1, 0, 1)
     &                     + a2t(iz,ix1,iy1, 1, 1))
                  db = abs(a2t(iz,ix1,iy1,-1, 1))
                  dc = abs(a2t(iz,ix1,iy1, 1, 1))
                  dn = max(da,db,dc)
c
                  da = abs(  a2t(iz,ix1,iy1,-1,-1) 
     &                     + a2t(iz,ix1,iy1, 0,-1)
     &                     + a2t(iz,ix1,iy1, 1,-1))
                  db = abs(a2t(iz,ix1,iy1,-1,-1))
                  dc = abs(a2t(iz,ix1,iy1, 1,-1))
                  ds = max(da,db,dc)
c
                  csum = (0.d0,0.d0)
                  do l = -1, 1
                     do k = -1, 1
                       csum = csum + a2(iz,ix1,iy1, k, l)
                     enddo
                  enddo
                  sigm = 0.5d0*
     &                   min(1.d0,abs(1.d0-csum/a2(iz,ix1,iy1,0,0)))
                  c2 = a2 (iz,ix1,iy1,-1, 1) + a2 (iz,ix1,iy1, 0, 1) +
     &                 a2 (iz,ix1,iy1, 1, 1) - a2t(iz,ix1,iy1,-1, 1) -
     &                 a2t(iz,ix1,iy1, 0, 1) - a2t(iz,ix1,iy1, 1, 1) -
     &                 a2 (iz,ix1,iy1,-1,-1) - a2 (iz,ix1,iy1, 0,-1) -
     &                 a2 (iz,ix1,iy1, 1,-1) + a2t(iz,ix1,iy1,-1,-1) +
     &                 a2t(iz,ix1,iy1, 0,-1) + a2t(iz,ix1,iy1, 1,-1)
c
                  if (ds+dn .gt. epsvs) then
                     w1s = sigm*(1.d0+(ds-dn)/(ds+dn)+c2/(dw+de+ds+dn))
                     w1n = 2.d0*sigm - w1s
c
                     w1s = max(0.d0,abs(w1s))*cone  ! w1s = max(0.d0,w1s)
                     w2s = max(0.d0,abs(w1n))*cone  ! w2s = max(0.d0,w1n) 
                     w1  = min(2.d0*sigm,abs(w1s)) ! w1  = min(2.d0*sigm,w1s)
                     w2  = min(2.d0*sigm,abs(w1n)) ! w2  = min(2.d0*sigm,w1n)
                     if (j .eq. iysf+1)   w1 = czero
                     if (j .eq. iysf+nyf) w2 = czero
                     witp(iz,i,j,1) = w1*cone
                     witp(iz,i,j,2) = w2*cone
                  else
                     witp(iz,i,j,1) = czero
                     witp(iz,i,j,2) = czero
                  endif
               enddo
            enddo
         enddo
c
c     on grid points 3
c
         do iy1 = iys3, nyf, 2
            do ix1 = ixs3, nxf, 2
               j = iysf+iy1
               i = ixsf+ix1
               do k = 1, nz
                  d1 = a2(k,ix1,iy1,-1,-1)
     &               + a2(k,ix1,iy1,0,-1)*witp(k,i,j-1,1)
     &               + a2(k,ix1,iy1,-1,0)*witp(k,i-1,j,1)
                  d2 = a2(k,ix1,iy1,1,-1)
     &               + a2(k,ix1,iy1,0,-1)*witp(k,i,j-1,2)
     &               + a2(k,ix1,iy1,1, 0)*witp(k,i+1,j,1)
                  d3 = a2(k,ix1,iy1,-1,1)
     &               + a2(k,ix1,iy1, 0,1)*witp(k,i,j+1,1)
     &               + a2(k,ix1,iy1,-1,0)*witp(k,i-1,j,2)
                  d4 = a2(k,ix1,iy1,1,1)
     &               + a2(k,ix1,iy1,0,1)*witp(k,i,j+1,2)
     &               + a2(k,ix1,iy1,1,0)*witp(k,i+1,j,2)
c
                  if (i .eq. ixsf+1   .and. j .eq. iysf+1  ) d1 = czero
                  if (i .eq. ixsf+nxf .and. j .eq. iysf+1  ) d2 = czero 
                  if (i .eq. ixsf+1   .and. j .eq. iysf+nyf) d3 = czero
                  if (i .eq. ixsf+nxf .and. j .eq. iysf+nyf) d4 = czero
c
                  d = -a2(k,ix1,iy1,0,0)
                  if (abs(d) .gt. epsvs) then
                     witp(k,i,j,1) = d1/d
                     witp(k,i,j,2) = d2/d
                     witp(k,i,j,3) = d3/d
                     witp(k,i,j,4) = d4/d
                  endif
               enddo
            enddo
         enddo
      endif
c
      return
      end
