      subroutine pmsa3d19(n,p,q,bet1,bet2,idata,omg,!a,
     &                    ia,ja,tempw, 
     &                    rpre,ipre,ip,omega)
c
      implicit real*8 (a-h,o-z)
      integer     ia(3),ipre(*),ip(7) 
c      complex*16  a(*),
      complex*16  p(n),q(n),tempw(*),rpre(*) 
      real*8      bet1,bet2 
c
c     Resolve ia
c
      nz = ia(1)
      nz2 = nz+2
      ngx = ia(2)
      ngy = ia(3)
c      write(*,*) nz,ngx,ngy 
c
c      icomm = ip(1)
      npz = ip(2)
      npx = ip(3)
      npy = ip(4)
      ipz = ip(5)
      ipx = ip(6)
      ipy = ip(7)
c
c     pointers
c
      mx = ipre(1)
      my = ipre(2)
c
      mx_ipre = 1
      my_ipre = 2
      nallx2_ipre = 3
      nally2_ipre = 4
      ids_ipre =5
      nit_ipre = 6
      nr1_ipre = 7
      nr2_ipre = 8
      nr3_ipre = 9
      nr4_ipre = 10
c
c
      iswp_ipre = 13
      icyc_ipre = 14
      md_ipre = 15
c
      neigz_ipre = 16
      izg_ipre = 18
c
      nx_ipre = 21
      ny_ipre = mx+21
      ixs_ipre = mx+my+21
      iys_ipre = 2*mx+my+21
      ix0_ipre = 2*mx+2*my+21
      iy0_ipre = 3*mx+2*my+21
      ixg_ipre = 3*mx+3*my+21
      iyg_ipre = 4*mx+3*my+21
      neigx_ipre = 4*mx+4*my+21
      neigy_ipre = 6*mx+4*my+21
c
      ids = ipre(ids_ipre) 
      nit = ipre(nit_ipre) 
      nr1 = ipre(nr1_ipre) 
      nr2 = ipre(nr2_ipre) 
      nr3 = ipre(nr3_ipre) 
      nr4 = ipre(nr4_ipre) 
      nallx2 = ipre(nallx2_ipre)
      nally2 = ipre(nally2_ipre)
c
c
      iswp = ipre(iswp_ipre)
      icyc = ipre(icyc_ipre)
      md = ipre(md_ipre)
c
c     Address of tempw. Length(tempw) = 3*nz2*nallx2*nally2
c
      nb = max((ngx+2)*(ngy+2),(ngx+2)*(nz+2))
      nb = max(nb,(ngy+2)*(nz+2))
      itw_r = 1
      itw_w = itw_r + nz2*nallx2*nally2
      itw_f = itw_w + nz2*nallx2*nally2
c
c     Address of rpre
c
c      irp_acor  = 1
c      irp_wrst  = irp_acor + 27*nz2*nallx2*nally2
c      irp_witph = irp_wrst +  9*nz2*nallx2*nally2
c      irp_witpv = irp_witph + 2*nz2*nallx2*nally2
c      irp_witpd = irp_witpv + 2*nz2*nallx2*nally2
c      irp_dinv  = irp_witpd + 4*nz2*nallx2*nally2
c      irp_dc    = irp_dinv + nz2*nallx2*nally2
c      irp_omega = irp_dc + nz2*nallx2*nally2
c      
      irp_acor  = 1
      irp_wrst  = irp_acor  + 27*nz2*nallx2*nally2
      irp_witpd = irp_wrst  +  9*nz2*nallx2*nally2
      irp_dinv  = irp_witpd + 4*nz2*nallx2*nally2
      irp_dc    = irp_dinv + nz2*nallx2*nally2
      irp_omega = irp_dc + nz2*nallx2*nally2
c
c      omega = rpre(irp_omega)
c
c     initilize of q
c
      do i = 1,nz*ngx*ngy
         q(i) = (0.d0,0.d0)
      enddo
c     
      call msg3da(nz,ngx,ngy,nallx2,nally2,p,q,
     &            rpre(irp_acor),rpre(irp_dinv),rpre(irp_dc), 
     &            mx,my,ipre(nx_ipre),ipre(ny_ipre),idata,omg,bet1,bet2, 
     &            omega,ids,nr1,nr2,nr3,nr4,nit, 
     &            icyc,md, 
     &            rpre(irp_wrst), 
     &            rpre(irp_witpd), 
     &            tempw(itw_r),tempw(itw_w), 
     &            tempw(itw_f),
     &            ipre(ixs_ipre),ipre(iys_ipre), 
     &            ipre(ix0_ipre),ipre(iy0_ipre), 
     &            iswp, 
     &            ipz,ipx,ipy,npz,npx,npy, 
     &            ipre(neigz_ipre), 
     &            ipre(neigx_ipre),ipre(neigy_ipre))
c
      return
      end
