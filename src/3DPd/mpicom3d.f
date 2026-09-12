      subroutine comm3da1(w, nx, ny, nz, neigx, neigy, neigz, 
     $     mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z) 
c
      return
      end
c
      subroutine commmsgx(work, nz, nallx2, nally2,
     $     ixs, iys, nx, ny,
     $     neigx, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z)
c
      return
      end
c
      subroutine commmsgy(work, nz, nallx2, nally2,
     $     ixs, iys, nx, ny,
     $     neigy, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z)
c
      return
      end
c
      subroutine commmsga(work, nz, nallx2, nally2, 
     $     ixs, iys, nx, ny,
     $     neigx, neigy, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z)
c
      return
      end
c
      subroutine commxyz(work, nz, nallx2, nally2, 
     $     ixs, iys, nx, ny,
     $     neigz, neigx, neigy, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z)
c
      return
      end
c
      subroutine findis3d(npx, npy, npz, ipx, ipy, ipz, 
     $           nx0, ny0, nz0, ixg, iyg, izg, 
     $           mpi_comm)
      implicit real*8 (a-h,o-z)
c
      ixg = 1
      iyg = 1
      izg = 1
      return
      end
c
      subroutine allgrid3d(npx, npy, npz, ipx, ipy, ipz, 
     $     nx0, ny0, nz0, nxa, nya, nza, 
     $     mpi_comm)
      implicit real*8 (a-h,o-z)
      nxa = nx0
      nya = ny0
      nza = nz0
      return
      end
c
      subroutine trzm(ws, wr, nz, nallx2, nally2, 
     $     ixs, iys, nx, ny, 
     $     ixss, iyss, ixes, iyes, ixsr, iysr, ixer, iyer, 
     $     neigz, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z) 
c
      return
      end
c
      subroutine trzp(ws, wr, nz, nallx2, nally2, 
     $     ixs, iys, nx, ny, 
     $     ixss, iyss, ixes, iyes, ixsr, iysr, ixer, iyer, 
     $     neigz, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z) 
c
      return
      end
c
      subroutine trzpr(ws, wr, nz, nallx2, nally2, 
     $     ixs, iys, nx, ny, 
     $     ixss, iyss, ixes, iyes, ixsr, iysr, ixer, iyer, 
     $     neigz, mpi_comm, sb, rb)
      implicit real*8 (a-h,o-z) 
c
      return
      end
c
      subroutine imax(ia,n,mpi_comm)
      implicit real*8 (a-h,o-z)
      return
      end
c
      subroutine isum(ia,n,mpi_comm)
      implicit real*8 (a-h,o-z)
      return
      end
c
      subroutine mpi_init(ierr)
      ierr = 0
      return
      end
c
      subroutine mpi_top_communicator(icomm,ierr)
      ierr = 0
      return
      end
c
      subroutine mpi_comm_rank(icomm, ipe, ierr)
      ipe = 0
      ierr = 0
      return
      end
c
      subroutine mpi_finalize(ierr)
      ierr = 0
      return
      end
c
      subroutine cjgettmrf(it)
      implicit real*8 (a-h,o-z)
c      call clock(t)
      it = t*1000
      return
      end

