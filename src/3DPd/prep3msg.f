      subroutine prep3msg(p,nz,nallx2,nally2, 
     &                    kst,ist,jst,izg,ixg,iyg,ixs,iys,nx,ny)
c
      implicit real*8 (a-h,o-z)
c
      complex*16  cone,zero,chalf       
      parameter  (eps = 1.d-8, izero = 0, ione = 1,
     &            czero = (0.d0,0.d0),cone = (1.d0,0.d0),
     &            chalf = (0.5d0,0.d0))
      complex*16  p(0:nz+1,nallx2,nally2) 
c
      do j = iys,iys+ny+1
         do i = ixs,ixs+nx+1
            do k = 0,nz+1
               p(k,i,j) = czero
            enddo
         enddo
      enddo
c
      iz03 = mod(izg-1,3) + 1
      ix03 = mod(ixg-1,3) + 1
      iy03 = mod(iyg-1,3) + 1
      do j = iys + mod(jst-iy03+3,3)+1,iys+ny,3
         do i = ixs + mod(ist-ix03+3,3)+1,ixs+nx,3
            do k = mod(kst-iz03+3,3)+1,nz,3
               p(k,i,j) = cone
            enddo
         enddo
      enddo
c
      return
      end
