      subroutine post3msg(y,acor,nz,nallx2,nally2, 
     &                    kst,ist,jst,izg,ixg,iyg,
     &                    ixs,iys,nx,ny,
     &                    neigz,neigx,neigy)
c
      implicit    real*8 (a-h,o-z)
c      
      complex*16  czero,cone,chalf 
      parameter  (eps = 1.d-8, izero = 0, ione = 1,
     &            czero = (0.d0,0.d0),cone = (1.d0,0.d0),
     &            chalf = (0.5d0,0.d0))
c      
      complex*16  y(0:nz+1,nallx2,nally2) 
      complex*16  acor(0:nz+1,nallx2,nally2,-1:1,-1:1,-1:1) 
      integer     neigz(2),neigx(2),neigy(2) 
c
      do i = ixs,ixs+nx+1
         do k = 0,nz+1
            y(k,i,iys) = czero
            y(k,i,iys+ny+1) = czero
         enddo
      enddo
      do j = iys,iys+ny+1
         do k = 0,nz+1
            y(k,ixs,j) = czero
            y(k,ixs+nx+1,j) = czero
         enddo
      enddo
      do j = iys,iys+ny+1
         do i = ixs,ixs+nx+1
            y(0,i,j) = czero
            y(nz+1,i,j) = czero
         enddo
      enddo
c
      iz03 = mod(izg-1,3) + 1
      ix03 = mod(ixg-1,3) + 1
      iy03 = mod(iyg-1,3) + 1
c
      do jo = -1,1
         do io = -1,1
            do ko = -1,1
               do j = iys + mod(jst-iy03+3,3)+1,iys+ny,3
                  do i = ixs + mod(ist-ix03+3,3)+1,ixs+nx,3
                     do k = mod(kst-iz03+3,3)+1,nz,3
                        acor(k,i,j,ko,io,jo) = y(k+ko,i+io,j+jo)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
      
