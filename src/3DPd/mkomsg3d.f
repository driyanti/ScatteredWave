      subroutine mkomsg3d(nz,ngx,ngy, 
     &                    nallx2,nally2, 
     &                    mx,my,ids, 
     &                    nx,ny,  
     &                    wrst, 
     &                    witpd, 
     &                    a2,a2t,
     &                    acor,dinv,dc, 
     &                    work1,work2, 
     &                    ixs,iys,ix0,iy0,ixg,iyg, 
     &                    izg, 
     &                    ipz,ipx,ipy,npz,npx,npy, 
     &                    neigz, 
     &                    neigx,neigy, 
     &                    ierr,
     &                    idata,omg,bet1,bet2)
c
c     make multi semi-coarsening grid environment
c     
      implicit   real*8 (a-h,o-z)
      complex*16  czero,cone,chalf 
      parameter (eps = 1.d-8, izero = 0, ione = 1, epsvs=1.d-20, 
     &           czero = (0.d0,0.d0), cone = (1.d0,0.d0),
     &           chalf = (0.5d0,0.d0), done = 1.d0,dzero = 0.d0)
c
      integer     idata 
      integer     nx(mx),ny(my) 
      integer     ixs(mx),iys(my),ix0(mx),iy0(my),ixg(mx),iyg(my) 
      integer     neigx(2,mx),neigy(2,my),neigz(2) 
      complex*16  wrst(0:nz+1,nallx2,nally2,9) 
      complex*16  witpd(0:nz+1,nallx2,nally2,4) 
      complex*16  a2(nz,ngx,ngy,9),a2t(nz,ngx,ngy,9) 
      complex*16  acor(0:nz+1,nallx2,nally2,27)  
      complex*16  dinv(0:nz+1,nallx2,nally2)
      complex*16  dc(0:nz+1,nallx2,nally2) 
      complex*16  work1(0:nz+1,nallx2,nally2) 
      complex*16  work2(0:nz+1,nallx2,nally2)  
      real*8      omg,bet1,bet2 
 
c
c     initialization 
c
      do j = 1, nally2
         do i = 1, nallx2
            do k = 0, nz+1
               wrst(k,i,j,1)  = czero
               wrst(k,i,j,2)  = czero
               wrst(k,i,j,3)  = czero
               wrst(k,i,j,4)  = czero
               wrst(k,i,j,5)  = czero
               wrst(k,i,j,6)  = czero
               wrst(k,i,j,7)  = czero
               wrst(k,i,j,8)  = czero
               wrst(k,i,j,9)  = czero
               witpd(k,i,j,1) = czero
               witpd(k,i,j,2) = czero
               witpd(k,i,j,3) = czero
               witpd(k,i,j,4) = czero
               acor(k,i,j,1) = czero
               acor(k,i,j,2) = czero
               acor(k,i,j,3) = czero
               acor(k,i,j,4) = czero
               acor(k,i,j,5) = czero
               acor(k,i,j,6) = czero
               acor(k,i,j,7) = czero
               acor(k,i,j,8) = czero
               acor(k,i,j,9) = czero
               acor(k,i,j,10) = czero
               acor(k,i,j,11) = czero
               acor(k,i,j,12) = czero
               acor(k,i,j,13) = czero
               acor(k,i,j,14) = czero ! cone
               acor(k,i,j,15) = czero
               acor(k,i,j,16) = czero
               acor(k,i,j,17) = czero
               acor(k,i,j,18) = czero
               acor(k,i,j,19) = czero
               acor(k,i,j,20) = czero
               acor(k,i,j,21) = czero
               acor(k,i,j,22) = czero
               acor(k,i,j,23) = czero
               acor(k,i,j,24) = czero
               acor(k,i,j,25) = czero
               acor(k,i,j,26) = czero
               acor(k,i,j,27) = czero
               dinv(k,i,j)    = czero
               dc(k,i,j)      = czero
               work1(k,i,j)   = czero
               work2(k,i,j)   = czero
            enddo
         enddo
      enddo
c
c     copy information on the finest grid
c
c      do iy = 1, ngy
c         do ix = 1, ngx
c            do k = 1, nz
c               acor(k,ix+ixs(mx),iy+iys(my),2 ) = a(k,ix,iy,1 )
c               acor(k,ix+ixs(mx),iy+iys(my),4 ) = a(k,ix,iy,2 )
c               acor(k,ix+ixs(mx),iy+iys(my),5 ) = a(k,ix,iy,3 )
c               acor(k,ix+ixs(mx),iy+iys(my),6 ) = a(k,ix,iy,4 )
c               acor(k,ix+ixs(mx),iy+iys(my),8 ) = a(k,ix,iy,5 )
c               acor(k,ix+ixs(mx),iy+iys(my),10) = a(k,ix,iy,6 )
c               acor(k,ix+ixs(mx),iy+iys(my),11) = a(k,ix,iy,7 )
c               acor(k,ix+ixs(mx),iy+iys(my),12) = a(k,ix,iy,8 )
c               acor(k,ix+ixs(mx),iy+iys(my),13) = a(k,ix,iy,9 )
c               acor(k,ix+ixs(mx),iy+iys(my),14) = a(k,ix,iy,10)
c               acor(k,ix+ixs(mx),iy+iys(my),15) = a(k,ix,iy,11)
c               acor(k,ix+ixs(mx),iy+iys(my),16) = a(k,ix,iy,12)
c               acor(k,ix+ixs(mx),iy+iys(my),17) = a(k,ix,iy,13)
c               acor(k,ix+ixs(mx),iy+iys(my),18) = a(k,ix,iy,14)
c               acor(k,ix+ixs(mx),iy+iys(my),20) = a(k,ix,iy,15)
c               acor(k,ix+ixs(mx),iy+iys(my),22) = a(k,ix,iy,16)
c               acor(k,ix+ixs(mx),iy+iys(my),23) = a(k,ix,iy,17)
c               acor(k,ix+ixs(mx),iy+iys(my),24) = a(k,ix,iy,18)
c               acor(k,ix+ixs(mx),iy+iys(my),26) = a(k,ix,iy,19)
c            enddo
c         enddo
c      enddo
c      
c      do iy = 1, ngy
c         do ix = 1, ngx
c            do k = 1, nz
c               acor(k,ix+ixs(mx),iy+iys(my), 5) = a(k,ix,iy,1)
c               acor(k,ix+ixs(mx),iy+iys(my),11) = a(k,ix,iy,2)
c               acor(k,ix+ixs(mx),iy+iys(my),13) = a(k,ix,iy,3)
c               acor(k,ix+ixs(mx),iy+iys(my),14) = a(k,ix,iy,4)
c               acor(k,ix+ixs(mx),iy+iys(my),15) = a(k,ix,iy,5)
c               acor(k,ix+ixs(mx),iy+iys(my),17) = a(k,ix,iy,6)
c               acor(k,ix+ixs(mx),iy+iys(my),23) = a(k,ix,iy,7)
c            enddo
c         enddo
c      enddo
c
      call mkacorf(acor,ngx,ngy,nz,nallx2,nally2,ixs(mx),iys(my),
     &             idata,omg,bet1,bet2)
c
c     Coarsening to the diagonal direction
c
      m = min(mx,my)
c
      do k = m-1,1,-1
c
         kx = mx-m+k
         ky = my-m+k
c
c     make witpd along the diagonal direction. 
c
         irst = 1
         if (irst .eq. 0) then
            if (k .eq. m-1) then
               call mkwitp3af(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,ids, ! ids = 0
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            else
               call mkwitp3ac(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,ids, ! ids = 0
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            endif
            call mkwrst3a(wrst, witpd,
     &                    nallx2,nally2, 
     &                    nx(kx),nx(kx+1),ixs(kx),ixs(kx+1),ix0(kx+1), 
     &                    ny(ky),ny(ky+1),iys(ky),iys(ky+1),iy0(ky+1), 
     &                    nz, 
     &                    neigx(1,kx+1),neigy(1,ky+1))         
         elseif (irst .eq. 1) then
            if (k .eq. m-1) then
               call mkwitp3af(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,0, ! ids = 0
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            else
               call mkwitp3ac(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,0, ! ids = 0
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            endif
c     
c     make wrst along the horizontal direction.
c
            call mkwrst3a(wrst, witpd,
     &                    nallx2,nally2, 
     &                    nx(kx),nx(kx+1),ixs(kx),ixs(kx+1),ix0(kx+1), 
     &                    ny(ky),ny(ky+1),iys(ky),iys(ky+1),iy0(ky+1), 
     &                    nz, 
     &                    neigx(1,kx+1), neigy(1,ky+1))
     
            if (k .eq. m-1) then
               call mkwitp3af(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,ids, ! ids = 2
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            else
               call mkwitp3ac(witpd,
     &                        nallx2,nally2, 
     &                        nx(kx+1),ixs(kx+1),ix0(kx+1), 
     &                        ny(ky+1),iys(ky+1),iy0(ky+1), 
     &                        nz, 
     &                        acor,ids, ! ids = 2
     &                        work1,a2,a2t, 
     &                        neigx(1,kx+1),neigy(1,ky+1))
            endif
         endif
         
c
         do jv = 1, 3
            do iv = 1, 3
               do kv = 1, 3
c
                  call prep3msg(work1, nz, nallx2, nally2, 
     &                 kv, iv, jv, izg, ixg(kx), iyg(ky), 
     &                 ixs(kx), iys(ky), nx(kx), ny(ky))
c
                  call ext3a(work1, work1, witpd, nallx2, nally2, 
     &                 nx(kx), nx(kx+1), ixs(kx), ixs(kx+1), ix0(kx+1), 
     &                 ny(ky), ny(ky+1), iys(ky), iys(ky+1), iy0(ky+1), 
     &                 nz,
     &                 neigx(1,kx+1), neigy(1,ky+1))
c
                  call mamsg27(acor,work1,work2,nallx2,nally2, 
     &                 ixs(kx+1),iys(ky+1),
     &                 nx(kx+1),ny(ky+1),
     &                 nz,
     &                 neigz,neigx(1,kx+1),neigy(1,ky+1))
c
                  call restr3a(work2,work2,wrst,nallx2,nally2, 
     &                 nx(kx),nx(kx+1),ixs(kx),ixs(kx+1),ix0(kx+1), 
     &                 ny(ky),ny(ky+1),iys(ky),iys(ky+1),iy0(ky+1), 
     &                 nz, 
     &                 neigx(1,kx+1), neigy(1,ky+1))
c
                  call post3msg(work2, acor, nz, nallx2, nally2, 
     &                 kv, iv, jv, izg, ixg(kx), iyg(ky), 
     &                 ixs(kx), iys(ky), nx(kx), ny(ky), 
     &                 neigz, neigx(1,kx), neigy(1,ky)) 
               enddo
            enddo
         enddo
c
         call trsp3msg(acor, work1, work2, nallx2, nally2, 
     &        ixs(kx), iys(ky), nx(kx), ny(ky),
     &        nz, 
     &        neigz, neigx(1,kx), neigy(1,ky)) 
c
      enddo
c
c     Compute dinv
c
      ipe = ipz + ipx*npz+ipy*npz*npx
      do ky = my,1,-1
         do kx = mx,1,-1
            ierr0 = 0
            call tridec(acor, dinv, dc, work1, work2, 
     &           nz, nallx2, nally2,
     &           ixs(kx), iys(ky), nx(kx), ny(ky),
     &           ierr, 
     &           ipe, ipz, npz)
         enddo
      enddo
c
      return
      end
