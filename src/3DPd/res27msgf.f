      subroutine res27msgf(!acor,
     &                     idata,omg,bet1,bet2,
     &                     q,p,r,nallx2,nally2, 
     &                     ixs,iys,ngx,ngy,nz, 
     &                     neigz,neigx,neigy)
c
      implicit real*8 (a-h,o-z)
c      complex*16  acor(0:nz+1,nallx2,nally2,27) 
      complex*16  q(0:nz+1,nallx2,nally2),p(0:nz+1,nallx2,nally2) 
      complex*16  r(0:nz+1,nallx2,nally2) 
      complex*16  a1,a2,a3,a4,a5,a6,a7,cidd 
      parameter   (cidd = (0.d0,1.d0))
c
c      do j = iys+1,iys+ny
c         do i = ixs+1,ixs+nx
c            do k = 1,nz
c               r(k,i,j) = p(k,i,j) 
!     &                  - acor(k,i,j,1 )*q(k-1,i-1,j-1)
!     &                  - acor(k,i,j,2 )*q(k  ,i-1,j-1)
!     &                  - acor(k,i,j,3 )*q(k+1,i-1,j-1)
!     &                  - acor(k,i,j,4 )*q(k-1,i  ,j-1)
c     &                  - acor(k,i,j,5 )*q(k  ,i  ,j-1)
!     &                  - acor(k,i,j,6 )*q(k+1,i  ,j-1)
!     &                  - acor(k,i,j,7 )*q(k-1,i+1,j-1)
!     &                  - acor(k,i,j,8 )*q(k  ,i+1,j-1)
!     &                  - acor(k,i,j,9 )*q(k+1,i+1,j-1)
!     &                  - acor(k,i,j,10)*q(k-1,i-1,j  )
c     &                  - acor(k,i,j,11)*q(k  ,i-1,j  )
!     &                  - acor(k,i,j,12)*q(k+1,i-1,j  )
c     &                  - acor(k,i,j,13)*q(k-1,i  ,j  )
c     &                  - acor(k,i,j,14)*q(k  ,i  ,j  )
c     &                  - acor(k,i,j,15)*q(k+1,i  ,j  )
!     &                  - acor(k,i,j,16)*q(k-1,i+1,j  )
c     &                  - acor(k,i,j,17)*q(k  ,i+1,j  )
!     &                  - acor(k,i,j,18)*q(k+1,i+1,j  )
!     &                  - acor(k,i,j,19)*q(k-1,i-1,j+1)
!     &                  - acor(k,i,j,20)*q(k  ,i-1,j+1)
!     &                  - acor(k,i,j,21)*q(k+1,i-1,j+1)
!     &                  - acor(k,i,j,22)*q(k-1,i  ,j+1)
c     &                  - acor(k,i,j,23)*q(k  ,i  ,j+1)
!     &                  - acor(k,i,j,24)*q(k+1,i  ,j+1)
!     &                  - acor(k,i,j,25)*q(k-1,i+1,j+1)
!     &                  - acor(k,i,j,26)*q(k  ,i+1,j+1)
!     &                  - acor(k,i,j,27)*q(k+1,i+1,j+1)
c            enddo
c         enddo
c      enddo
c
      ngx2 = ngx-2
      ngy2 = ngy-2
      nz2  = nz -2
c      
      hx  =  1.d0/dble(float(ngx2+1))
      hy  =  1.d0/dble(float(ngy2+1))  
      hz  =  1.d0/dble(float(nz2+1))  
      hx2 =  1.d0/hx**2
      hy2 =  1.d0/hy**2
      hz2 =  1.d0/hz**2
      hh  =  hx
      hh2 =  hx2
c
c     Interior points
c                
      do j = iys+2,iys+ngy-1
         iyo = j-2
         do i = ixs+2,ixs+ngx-1
            ixo = i-2
            do k = 2,nz-1
               izo = k-1
               wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)  
               a1   = -cmplx(hy2,0.d0)
               a2   = -cmplx(hx2,0.d0)
               a3   = -cmplx(hz2,0.d0)
               a4   = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2 
               a5   = -cmplx(hz2,0.d0)
               a6   = -cmplx(hx2,0.d0)
               a7   = -cmplx(hy2,0.d0)
               r(k,i,j) = p(k,i,j) 
     &                  - a1*q(k  ,i  ,j-1)
     &                  - a2*q(k  ,i-1,j  )
     &                  - a3*q(k-1,i  ,j  )
     &                  - a4*q(k  ,i  ,j  )
     &                  - a5*q(k+1,i  ,j  ) 
     &                  - a6*q(k  ,i+1,j  )
     &                  - a7*q(k  ,i  ,j+1)
            enddo
         enddo
      enddo
c
c     Boundary points ::
c
c     Faces ::
c     North
      j   = iys+ngy
      iyo = j-2
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)      
            a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
            a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a3 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &           + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2  
            a5 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            r(k,i,j) = p(k,i,j) 
     &               - a1*q(k  ,i  ,j-1)
     &               - a2*q(k  ,i-1,j  )
     &               - a3*q(k-1,i  ,j  )
     &               - a4*q(k  ,i  ,j  )
     &               - a5*q(k+1,i  ,j  ) 
     &               - a6*q(k  ,i+1,j  )
         enddo
      enddo
c
c     South
      j   = iys+1
      iyo = j-2
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)  
            a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a3 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &           + 2.d0*cidd/(wnum*hy)*((wnum*hy)**2-2.d0)*hy2
            a5 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
            r(k,i,j) = p(k,i,j) 
     &               - a2*q(k  ,i-1,j  )
     &               - a3*q(k-1,i  ,j  )
     &               - a4*q(k  ,i  ,j  )
     &               - a5*q(k+1,i  ,j  ) 
     &               - a6*q(k  ,i+1,j  )
     &               - a7*q(k  ,i  ,j+1)   
         enddo
      enddo
c
c     East
      i = ixs+ngx
      ixo = i-2
      do j = iys+2,iys+ngy-1
         iyo = j-2
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo) 
            a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
            a3 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2       
     &           + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            a5 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            r(k,i,j) = p(k,i,j) 
     &               - a1*q(k  ,i  ,j-1)
     &               - a2*q(k  ,i-1,j  )
     &               - a3*q(k-1,i  ,j  )
     &               - a4*q(k  ,i  ,j  )
     &               - a5*q(k+1,i  ,j  ) 
     &               - a7*q(k  ,i  ,j+1)
         enddo
      enddo
c
c     West
      i   = ixs+1
      ixo = i-2
      do j = iys+2,iys+ngy-1
         iyo = j-2
         do k = 2,nz-1
            izo = k-1
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo) 
            a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            a3 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2 
     &           + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            a5 = - cmplx(hz2,0.d0) + cidd/(wnum*hz)*hz2
            a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
            a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            r(k,i,j) = p(k,i,j) 
     &               - a1*q(k  ,i  ,j-1)
     &               - a3*q(k-1,i  ,j  )
     &               - a4*q(k  ,i  ,j  )
     &               - a5*q(k+1,i  ,j  ) 
     &               - a6*q(k  ,i+1,j  )
     &               - a7*q(k  ,i  ,j+1)
         enddo
      enddo
c
c     Right
      k   = nz
      izo = k-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         do i = ixs+2,ixs+ngx-1
            ixo = i-2
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2 
     &           + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            r(k,i,j) = p(k,i,j) 
     &               - a1*q(k  ,i  ,j-1)
     &               - a2*q(k  ,i-1,j  )
     &               - a3*q(k-1,i  ,j  )
     &               - a4*q(k  ,i  ,j  ) 
     &               - a6*q(k  ,i+1,j  )
     &               - a7*q(k  ,i  ,j+1)
         enddo
      enddo
c
c     Left
      k   = 1
      izo = k-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         do i = ixs+2,ixs+ngx-1
            ixo = i-2
            wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
            a1 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            a2 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &           + 2.d0*cidd/(wnum*hx)*((wnum*hx)**2-2.d0)*hx2
            a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
            a6 = - cmplx(hx2,0.d0) + cidd/(wnum*hx)*hx2
            a7 = - cmplx(hy2,0.d0) + cidd/(wnum*hy)*hy2
            r(k,i,j) = p(k,i,j) 
     &               - a1*q(k  ,i  ,j-1)
     &               - a2*q(k  ,i-1,j  )
     &               - a4*q(k  ,i  ,j  )
     &               - a5*q(k+1,i  ,j  ) 
     &               - a6*q(k  ,i+1,j  )
     &               - a7*q(k  ,i  ,j+1)
         enddo
      enddo 
c
c     Edges ::
c     North face ::
c     (2,ngy,1) --> (ngx-1,ngy,1)
      j   = iys+ngy
      iyo = j-2
      k   = 1
      izo = k-1
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a2*q(k  ,i-1,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )
      enddo
c
c     (2,ngy,nz) --> (ngx-1,ngy,nz)
      j = iys+ngy
      iyo = j-2
      k = nz
      izo = k-1
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a2*q(k  ,i-1,j  )
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a6*q(k  ,i+1,j  )
       
      enddo
c     (1,ngy,2) --> (1,ngy,nz-1)
      i   = ixs+1
      ixo = i-2
      j   = iys+ngy
      iyo = j-2
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         a3 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a5 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)   
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )   
      enddo 
c
c     (ngx,ngy,2) --> (ngx,ngy,nz-1)
      i = ixs+ngx
      ixo = i-2
      j = iys+ngy
      iyo = j-2
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a3 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a5 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a2*q(k  ,i-1,j  )
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
      enddo
c
c     South face ::
c     (2,1,1) --> (ngx-1,1,1)
      j = iys+1
      iyo = j-2
      k = 1
      izo = k-1
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         r(k,i,j) = p(k,i,j) 
     &            - a2*q(k  ,i-1,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo 
c
c     (2,1,nz) --> (ngx-1,1,nz)
      j = iys+1
      iyo = j-2  
      k = nz
      izo = k-1
      do i = ixs+2,ixs+ngx-1
         ixo = i-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)  
         a2 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hx)*(2.d0*cidd-3.d0*(wnum*hx)**2)*hx2
         a6 = - cmplx(hx2,0.d0) + 1.d0*cidd/(wnum*hx)*hx2
         a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         r(k,i,j) = p(k,i,j) 
     &            - a2*q(k  ,i-1,j  )
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a6*q(k  ,i+1,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     (1,1,2) --> (1,1,nz-1)
      i = ixs+1
      ixo = i-2
      j = iys+1
      iyo = j-2
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a3 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a5 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         r(k,i,j) = p(k,i,j) 
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     (ngx,1,2) --> (ngx,1,nz-1)
      i = ixs+ngx
      ixo = i-2
      j = iys+1
      iyo = j-2
      do k = 2,nz-1
         izo = k-1
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a3 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hz)*(2.d0*cidd-3.d0*(wnum*hz)**2)*hz2
         a5 = - cmplx(hz2,0.d0) + 1.d0*cidd/(wnum*hz)*hz2
         a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
         r(k,i,j) = p(k,i,j) 
     &            - a2*q(k  ,i-1,j  )
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     Left face ::
c     (1,2,1) --> (1,ngy-1,1)
      i = ixs+1
      ixo = i-2
      k = 1
      izo = i-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     (ngx,2,1) --> (ngx,ngy-1,1)
      i = ixs+ngx
      ixo = i-2
      k = 1
      izo = k-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a2*q(k  ,i-1,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a5*q(k+1,i  ,j  ) 
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     Right face ::
c     (1,2,nz) --> (1,ngy-1,nz)
      i = ixs+1
      ixo = i-2
      k = nz
      izo = k-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  ) 
     &            - a6*q(k  ,i+1,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     (ngx,2,nz) --> (ngx,ngy-1,nz)
      i = ixs+ngx
      ixo = i-2
      k = nz
      izo = k-1
      do j = iys+2,iys+ngy-1
         iyo = j-2
         wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
         a1 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
         a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
         a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &        - 1.d0/(wnum*hy)*(2.d0*cidd-3.d0*(wnum*hy)**2)*hy2
         a7 = - cmplx(hy2,0.d0) + 1.d0*cidd/(wnum*hy)*hy2
         r(k,i,j) = p(k,i,j) 
     &            - a1*q(k  ,i  ,j-1)
     &            - a2*q(k  ,i-1,j  )
     &            - a3*q(k-1,i  ,j  )
     &            - a4*q(k  ,i  ,j  )
     &            - a7*q(k  ,i  ,j+1)     
      enddo
c
c     Corner points ::
c
c     (1,1,1)
      i = ixs+1
      ixo = i-2
      j = iys+1
      iyo = j-2
      k = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
      a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      r(k,i,j) = p(k,i,j) 
     &         - a4*q(k  ,i  ,j  )
     &         - a5*q(k+1,i  ,j  ) 
     &         - a6*q(k  ,i+1,j  )
     &         - a7*q(k  ,i  ,j+1)  
c
c     (ngx,1,1)
      i = ixs+ngx
      ixo = i-2
      j = iys+1
      iyo = j-2
      k = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2.d0
     &     + 4.d0*cidd*wnum*hh*hh2
      a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)      
      r(k,i,j) = p(k,i,j) 
     &         - a2*q(k  ,i-1,j  )
     &         - a4*q(k  ,i  ,j  )
     &         - a5*q(k+1,i  ,j  ) 
     &         - a7*q(k  ,i  ,j+1)   
c
c     (1,ngy,1)
      i = ixs+1
      ixo = i-2
      j = iys+ngy
      iyo = j-2
      k = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)   
      r(k,i,j) = p(k,i,j) 
     &         - a1*q(k  ,i  ,j-1)
     &         - a4*q(k  ,i  ,j  )
     &         - a5*q(k+1,i  ,j  ) 
     &         - a6*q(k  ,i+1,j  )
c
c     (ngx,ngy,1)
      i = ixs+ngx
      ixo = i-2
      j = iys+ngy
      iyo = j-2
      k = 1
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a5 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      r(k,i,j) = p(k,i,j) 
     &         - a1*q(k  ,i  ,j-1)
     &         - a2*q(k  ,i-1,j  )
     &         - a4*q(k  ,i  ,j  )
     &         - a5*q(k+1,i  ,j  ) 
c
c     (1,1,nz)
      i = ixs+1
      ixo = i-2
      j = iys+1
      iyo = j-2
      k = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0) 
      a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      r(k,i,j) = p(k,i,j) 
     &         - a3*q(k-1,i  ,j  )
     &         - a4*q(k  ,i  ,j  ) 
     &         - a6*q(k  ,i+1,j  )
     &         - a7*q(k  ,i  ,j+1)
c
c     (ngx,1,nz)
      i = ixs+ngx
      ixo = i-2
      j = iys+1
      iyo = j-2
      k = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
      a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a7 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      r(k,i,j) = p(k,i,j) 
     &         - a2*q(k  ,i-1,j  )
     &         - a3*q(k-1,i  ,j  )
     &         - a4*q(k  ,i  ,j  )
     &         - a7*q(k  ,i  ,j+1)
c
c     (1,ngy,nz)
      i = ixs+1
      ixo = i-2
      j = iys+ngy
      iyo = j-2
      k = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      a6 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
      r(k,i,j) = p(k,i,j) 
     &         - a1*q(k  ,i  ,j-1)
     &         - a3*q(k-1,i  ,j  )
     &         - a4*q(k  ,i  ,j  )
     &         - a6*q(k  ,i+1,j  )
c
c     (ngx,ngy,nz)
      i = ixs+ngx
      ixo = i-2
      j = iys+ngy
      iyo = j-2
      k = nz
      izo = k-1
      wnum = wdata(omg,ngx2,ngy2,nz2,idata,izo,ixo,iyo)
      a1 = - cmplx(hy2,0.d0) - cmplx(hy2,0.d0)
      a2 = - cmplx(hx2,0.d0) - cmplx(hx2,0.d0)
      a3 = - cmplx(hz2,0.d0) - cmplx(hz2,0.d0)
      a4 = 2.d0*(hx2+hy2+hz2) - cmplx(bet1,-bet2)*wnum**2
     &     + 4.d0*cidd*wnum*hh*hh2
      r(k,i,j) = p(k,i,j) 
     &         - a1*q(k  ,i  ,j-1)
     &         - a2*q(k  ,i-1,j  )
     &         - a3*q(k-1,i  ,j  )
     &         - a4*q(k  ,i  ,j  )
c
c      print*,'res27',k,i,j,r(k,i,j),p(k,i,j),q(k,i,j-1),q(k-1,i,j)

c       do k=1,nz
c        do j= iys+1,iys+ngy
c	 do i= ixs+1,ixs+ngx
c              write(40,*)k,i,j,r(k,i,j),p(k,i,j),q(k,i,j)
c          enddo
c	 enddo
c	enddo   

      return
      end
