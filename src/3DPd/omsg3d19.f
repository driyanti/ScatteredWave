      subroutine omsg3d19(ia,
     &                    tempw,ntempw,rpre,nrpre,ipre,nipre,
     &                    a2,a2t,
     &                    icyc,md, 
     &                    ip,
     &                    ids,iswp,mx,my,nit,nr1,nr2,nr3,nr4, 
     &                    omega,
     &                    idata,omg,bet1,bet2)
c
c     Prepare the multi semi-coarsening grid preconditioner. 
c
c     nallx2 = nx(mx)+nx(mx-1) + ... + nx(1) + 2*(mx+1)
c     nally2 = ny(my)+ny(my-1) + ... + ny(1) + 2*(my+1)
c     
c     ipre(1:nipre):  nipre >= 6*mx+6*my+40
c
c
c     tempw(1:ntemw): temporary working area for the preconditioning
c                     ntempw >=(nz+2)*nallx2*nally2*3 
c                             + 2*nb
c
c        contents:    r(0:nz+1,nallx2,nally2)
c                     w(0:nz+1,nallx2,nally2)
c                     f(0:nz+1,nallx2,nally2)
c                     sb and rb (send and recieve buffers)
c
c     rpre(1:nrpre): real data for the preconditioning
c                    nrpre >= (nz+2)*nallx2*nally2*46+1
c
c         contents:   acor(0:nz+1,nallx2,nally2,27)
c                     wrst(0:nz+1,nallx2,nally2,9)
c                     witph(0:nz+1,nallx2,nally2,2)
c                     witpv(0:nz+1,nallx2,nally2,2)
c                     witpd(0:nz+1,nallx2,nally2,4)
c                     dinv(0:nz+1,nallx2,nally2)
c                     dc(0:nz+1,nallx2,nally2)
c                     omega
c
      implicit real*8 (a-h,o-z)
c
      complex*16  a2(*),a2t(*) 
      integer     ia(3) 
      integer     ipre(nipre),ip(7) 
      complex*16  tempw(ntempw),rpre(nrpre) 
c
      ierr = 0
c
c     Resolve ia
c
      nz = ia(1)
      nz2 = nz+2
      ngx = ia(2)
      ngy = ia(3)
c
      icomm = ip(1)
      npz = ip(2)
      npx = ip(3)
      npy = ip(4)
      ipz = ip(5)
      ipx = ip(6)
      ipy = ip(7)
c
c     Check multigrid parameters
c
c      if ((ids.lt.0).or.(ids.gt.3)) ierr = 20
c      call imax(ierr, 1, icomm)
c      if (ierr.ne.0) then
c         return
c      endif
c      if ((mx.lt.1).or.(my.lt.1)) ierr = 21
c      call imax(ierr, 1, icomm)
c      if (ierr.ne.0) then
c         return
c      endif
c      if (nit.le.0) ierr = 22
c      call imax(ierr, 1, icomm)
c      if (ierr.ne.0) then
c         return
c      endif
c      if ((icyc.lt.0).or.(icyc.gt.2)) ierr = 24
c      call imax(ierr, 1, icomm)
c      if (ierr.ne.0) then
c         return
c      endif
c
c
c     Setup of ipre
c
c
c      if (nipre.lt.6*mx+6*my+40) ierr = 5
c      call imax(ierr, 1, icomm)
c      if (ierr.ne.0) then
c         return
c      endif
c
c     pointers
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
      ipre(mx_ipre) = mx
      ipre(my_ipre) = my
      ipre(ids_ipre) = ids
      ipre(nit_ipre) = nit
      ipre(nr1_ipre) = nr1
      ipre(nr2_ipre) = nr2
      ipre(nr3_ipre) = nr3
      ipre(nr4_ipre) = nr4
c
c
      ipre(iswp_ipre) = iswp
      ipre(icyc_ipre) = icyc
      ipre(md_ipre) = md
c
c      write(*,*) 'OMSG3D'
c      write(*,*) ngx,ngy,nz 
      call defmsg3d(nz,ngx,ngy,mx,my,
     &              ipre(nx_ipre),ipre(ny_ipre), 
     &              nallx2,nally2, 
     &              npz,npx,npy,ipz,ipx,ipy, 
     &              ipre(ixs_ipre),ipre(iys_ipre), 
     &              ipre(ix0_ipre),ipre(iy0_ipre), 
     &              ipre(ixg_ipre),ipre(iyg_ipre), 
     &              ipre(izg_ipre), 
     &              ipre(neigz_ipre), 
     &              ipre(neigx_ipre),ipre(neigy_ipre))
c
      ipre(nallx2_ipre) = nallx2
      ipre(nally2_ipre) = nally2
c
c     check ntempw and set pointers to tempw
c
c
      nb = max((ngx+2)*(ngy+2),(ngx+2)*(nz+2))
      nb = max(nb,(ngy+2)*(nz+2))
c
      nrw = nz2*nallx2*nally2*3 + 2+nb
c
      do j = 1,nrw 
         tempw(j) = (0.d0,0.d0)
      enddo
c
c     Address of tempw. Length(tempw) = 3*nz2*nallx2*nally2
c
      itw_r = 1
      itw_w = itw_r + nz2*nallx2*nally2
      itw_f = itw_w + nz2*nallx2*nally2
c
c     check nrpre and set pointers to rpre
c
      nrw = nz2*nallx2*nally2*46 + 1
c
      irp_acor  = 1
      irp_wrst  = irp_acor  + 27*nz2*nallx2*nally2
      irp_witpd = irp_wrst  + 9*nz2*nallx2*nally2
      irp_dinv  = irp_witpd + 4*nz2*nallx2*nally2
      irp_dc    = irp_dinv + nz2*nallx2*nally2
      irp_omega = irp_dc + nz2*nallx2*nally2
c
c     check natw
c
      nrw = 2*nz*ngx*ngy*9
c      do j = 1, nrw
c         atw(j) = (0.d0,0.d0)
c      enddo
c
      call mkomsg3d(nz,ngx,ngy, 
     &              nallx2,nally2, 
     &              mx,my,ids, 
     &              ipre(nx_ipre),ipre(ny_ipre), 
     &              rpre(irp_wrst),
     &              rpre(irp_witpd), 
     &              a2,a2t,
     &              rpre(irp_acor),rpre(irp_dinv),rpre(irp_dc), 
     &              tempw(itw_r),tempw(itw_f), 
     &              ipre(ixs_ipre),ipre(iys_ipre), 
     &              ipre(ix0_ipre),ipre(iy0_ipre),
     &              ipre(ixg_ipre),ipre(iyg_ipre),
     &              ipre(izg_ipre), 
     &              ipz,ipx,ipy,npz,npx,npy, 
     &              ipre(neigz_ipre), 
     &              ipre(neigx_ipre),ipre(neigy_ipre), 
     &              ierr,idata,omg,bet1,bet2)
c
      return
      end
