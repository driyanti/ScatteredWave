      subroutine trsp3msg(acor,work1,work2,nallx2,nally2, 
     &                    ixs,iys,nx,ny,
     &                    nz,
     &                    neigz,neigx,neigy)
c
      implicit real*8 (a-h,o-z)
c      
      complex*16  cone,czero,chalf 
      parameter (eps = 1.d-8, izero = 0, ione = 1,
     &           czero = (0.d0,0.d0),cone = (1.d0,0.d0), 
     &           chalf = (0.5d0,0.d0))
      complex*16  acor(0:nz+1,nallx2,nally2,-1:1,-1:1,-1:1) 
      complex*16  work1(0:nz+1,nallx2,nally2)  
      complex*16  work2(0:nz+1,nallx2,nally2) 
c
      do icy = -1, 0
         if (icy.eq.-1) iex = 1
         if (icy.eq.0)  iex = 0
         do icx = -1, iex
            if (icy.eq.-1) then
               iez = 1
            elseif (icx.eq.-1) then
               iez = 1
            else
               iez = -1
            endif 
            do icz = -1, iez
               do j = iys, iys+ny+1
                  do i = ixs, ixs+nx+1
                     do k = 0 , nz+1
                        work1(k,i,j) = 0.d0
                        work2(k,i,j) = 0.d0
                     enddo
                  enddo
               enddo
               do j = iys+1, iys+ny
                  do i = ixs+1, ixs+nx
                     do k = 1 , nz
                        work1(k,i,j) = acor(k,i,j,-icz,-icx,-icy)
                        work2(k,i,j) = acor(k,i,j, icz, icx, icy)
                     enddo
                  enddo
               enddo
c
               do j = iys+1, iys+ny
                  do i = ixs+1, ixs+nx
                     do k = 1 , nz
                           acor(k,i,j,-icz,-icx,-icy) 
     &                    = work2(k-icz,i-icx,j-icy)
                        acor(k,i,j, icz, icx, icy) 
     &                    = work1(k+icz,i+icx,j+icy)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
c
      return
      end
