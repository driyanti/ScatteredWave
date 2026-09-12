      subroutine das3dx(ia,nia,
     &                  ip,nip,
     &                  ngx,ngy,ngz, 
     &                  npx,npy,npz, 
     &                  ipx,ipy,ipz)
c
c     This subroutine prepares data to perform the 2D 9
c     matrix vector product. 
c 
      implicit none
      integer nia,nip 
      integer ia(nia),ip(nip) 
      integer ngx,ngy,ngz,npx,npy,npz,ipx,ipy,ipz  
      integer ierr,icomm  
c 
      ierr = 0
      icomm = 1
c
c     prepare ia
c 
      if (nia.lt.3) ierr = 1
c
      ia(1) = ngx
      ia(2) = ngy
      ia(3) = ngz
c
c     prepare the data for parallel execution
c
      if (nip.lt.7) ierr = 10
c
      ip(1) = icomm
      ip(2) = npx
      ip(3) = npy
      ip(4) = npz
      ip(5) = ipx
      ip(6) = ipy
      ip(7) = ipz
c
      return
      end
