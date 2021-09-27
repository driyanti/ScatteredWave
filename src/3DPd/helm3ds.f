C###{__Helmholtz_3D__}
C###{__Here_we_read_parameters__}
#define _YCUT_ 1
#define _LOGS_ 1
C###{}
      PROGRAM helm3ds

      IMPLICIT NONE
      include 'mpif.h'
      include 'prms3d.h'
      include 'fmpiu.h'
      include 'grids.h'
C###{__MPI_stuff__}
      integer GDIM,pdims(0:1),plocal(0:1)
      parameter(GDIM=2)
      logical prd_flg(0:1),remain_dims(0:1),reordr
C###{__rmp_holder_of_some_often_used_data__}
      integer nia,nip,irpm
      parameter(nia=3,nip=7,irpm=8)
      character logfl*14

      integer ngx,ngy,ngz,npx,npy,npz,pid
      integer idata,irhs,m1,m2,nr1,nr2,nr3,nr4
      integer nit,icyc,md,ids,iswp,itsol,itmax,itol
      integer row_bsz,col_bsz,zlr_bsz,nxall,nyall,mlev,nbf
      integer ipm(4),ia(nia),ja(nia),ip(nip)
      real*8 bet1,bet2,omega,omg,Lx,Ly,Lz,timing,tt,rpm(irpm)

C###{__GO__}
C###{__Init_MPI-enviroment__}
      call MPI_Init(ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD,pe_id,ierr)
      call MPI_Comm_size(MPI_COMM_WORLD,pe_size,ierr)

C###{__[TIMER:] __INIT__ }
      tt = 0.d0
      tt = timing(tt)

C###{__Read_data_file__}
      open(10,file='3dmsg.input',status='old')
      read(10,*) omg
      read(10,*) Lx,Ly,Lz
      read(10,*) npx,npy,npz
      read(10,*) ngx,ngy,ngz
      read(10,*) nax,nay,naz
      read(10,*) bet1,bet2
      read(10,*) idata,irhs
      read(10,*) m1,m2
      read(10,*) nr1,nr2,nr3,nr4
      read(10,*) nit
      read(10,*) icyc,md
      read(10,*) omega
      read(10,*) ids,iswp
      read(10,*) itsol,itmax,itol 
      close(10)

C###{__Make_boundary_offsets__}
      nglx = ngx +2
      ngly = ngy +2
      nglz = ngz +2

C###{__Store_Grid_Sizes_}
      rpm(ihx) = Lx/dble(float(nglx-1))
      rpm(ihy) = Ly/dble(float(ngly-1))
      rpm(ihz) = Lz/dble(float(nglz-1))
      rpm(ihx2)= 1.d0/rpm(ihx)**2
      rpm(ihy2)= 1.d0/rpm(ihy)**2
      rpm(ihz2)= 1.d0/rpm(ihz)**2

C###{__Take_into_account_[AbsLayer]__}
      nglx = ngx + 2*nax + 2
      ngly = ngy + 2*nay + 2 
      nglz = ngz + 2*naz + 2
      ndx  = ngx
      ndy  = ngy
      ndz  = ngz
C###{__Find_about_processor_2D_GRID__}
C###{__if_number_of_procs_le_3_than_slice_the_domain__}
      if(pe_size.le.3)then
#ifdef _YCUT_
        pdims(0)= pe_size
        pdims(1)= 1
#else
        pdims(0)= 1
        pdims(1)= pe_size
#endif
      else
        if(npx*npy.ne.pe_size)then
          print*,"Wrong_dimensions_of_PROC_GRID"
          call pehalt(pe_id,"ERR: npx*npy.ne.pe_size")
          !!call MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER,ierr)
        else
          pdims(0)=npy !!procs_rows<=[Y]
          pdims(1)=npx !!procs_cols<=[X] 
        endif
      endif
      prd_flg(0)=.false.
      prd_flg(1)=.false.
      reordr=.true.

C###{__Get_balanced_distribution_of_processes__}
      call MPI_Dims_create(pe_size,GDIM,pdims,ierr)
      if(pe_id.eq.0)write(*,100)pe_size,pdims(0),pdims(1)

      call pe_descr_alloc(pe_descr)
      pe_descr(PE,ID) = pe_id
      pe_descr(PE,SZ) = pe_size

      pe_descr(PGR,RSZ)=pdims(0)
      pe_descr(PGR,CSZ)=pdims(1)

      pe_descr(PE,ZCUT) = NO  !!No Z-direction-cut
      pe_descr(PE,RWWS) = YES !!RowWise block-str. decm.
      
      if(pe_descr(PGR,CSZ).ge.1)then
        pe_descr(PE,CLWS) = YES
      else
        pe_descr(PE,CLWS) = NO !ColumnWise block-str. decm.
      endif

C###{__Create_Cartesian_2D_Communicator__}
      call MPI_Cart_create(MPI_COMM_WORLD,GDIM,pdims,
     &  prd_flg,reordr,pe_descr(PE,COMM_2D),ierr)
      if(ierr.gt.0)then
        if(pe_id.eq.0)print*,"FAILED->MPI_Cart_create()"
        call MPI_Abort(MPI_COMM_WORLD,MPI_ERR_OTHER,ierr)
      endif

C###{__Determine the PS-position in the grid and split ...}
C###{__communicator_into_row_comm_and_col_comm__}
      call MPI_Cart_coords(pe_descr(PE,COMM_2D),pe_id,
     &  GDIM,plocal,ierr)
      pe_descr(PGR,YPS)=plocal(0)
      pe_descr(PGR,XPS)=plocal(1)

C###{__Get row and column communicators using ...}
C###{__cartesian sub-topology__}
      remain_dims(0) = .FALSE.
      remain_dims(1) = .TRUE.
      call MPI_Cart_sub(pe_descr(PE,COMM_2D),remain_dims,
     &  pe_descr(PGR,ROW_COMM),ierr)

      remain_dims(0) = .TRUE.
      remain_dims(1) = .FALSE.
      call MPI_Cart_sub(pe_descr(PE,COMM_2D),remain_dims,
     &  pe_descr(PGR,COL_COMM),ierr)

C###{__SetUp_the_Neighborhood__NEW_(!!!)}
C###{__[AND]_Boundaries_FLAGS__}
      call set_neighbrsN(pe_descr)
      
C###{__Setting up PE - enviroment}
      pid=pe_descr(PGR,YPS)
      pe_descr(ROW,LOW) = blck_low_bnd(pid,pdims(0),ngly)
      pe_descr(ROW,HIGH) = blck_high_bnd(pid,pdims(0),ngly)
      pe_descr(ROW,BSIZE) = blck_size(pid,pdims(0),ngly)
      pe_descr(ROW,GSIZE) = ngly

      pid=pe_descr(PGR,XPS)
      pe_descr(COL,LOW) = blck_low_bnd(pid,pdims(1),nglx)
      pe_descr(COL,HIGH) = blck_high_bnd(pid,pdims(1),nglx)
      pe_descr(COL,BSIZE) = blck_size(pid,pdims(1),nglx)
      pe_descr(COL,GSIZE) = nglx

      pid=0 !!__NO_Z_cut__
      pe_descr(ZLR,LOW) = blck_low_bnd(pid,1,nglz)
      pe_descr(ZLR,HIGH) = blck_high_bnd(pid,1,nglz)
      pe_descr(ZLR,BSIZE) = blck_size(pid,1,nglz)
      pe_descr(ZLR,GSIZE) = nglz

C###{__Now_Set_Parallel_Block_Sizes__2D_PROC_GRID !!!}
C###{__Vector_Signature :: {Z,X,Y}__}
      pe_descr(ZLR,PLSZ) = pe_descr(ZLR,BSIZE)
      pe_descr(COL,PLSZ) = pe_descr(COL,BSIZE)
      pe_descr(ROW,PLSZ) = pe_descr(ROW,BSIZE)

      if(pe_descr(BND,NORTH).eq.YES)then
        pe_descr(ROW,PLSZ) = pe_descr(ROW,PLSZ)-1
        pe_descr(ROW,PLHI) = YES
      endif
      if(pe_descr(BND,SOUTH).eq.YES)then
        pe_descr(ROW,PLSZ) = pe_descr(ROW,PLSZ)-1
      endif
      if(pe_descr(BND,EAST).eq.YES)then
        pe_descr(COL,PLSZ) = pe_descr(COL,PLSZ)-1
        pe_descr(COL,PLHI) = YES
      endif
      if(pe_descr(BND,WEST).eq.YES)then
        pe_descr(COL,PLSZ) = pe_descr(COL,PLSZ)-1
      endif
      if(pe_descr(BND,ZTOP).eq.YES)then
        pe_descr(ZLR,PLSZ) = pe_descr(ZLR,PLSZ)-1
        pe_descr(ZLR,PLHI) = YES
      endif
      if(pe_descr(BND,ZBTM).eq.YES)then
        pe_descr(ZLR,PLSZ) = pe_descr(ZLR,PLSZ)-1
      endif
C###{__Set_[Z,X,Y]_vector_limits__}
C###{__Limits_[Z]__}
      if(pe_descr(BND,ZTOP).eq.YES)then
        pe_descr(XYZ,ZHI) = pe_descr(ZLR,PLSZ)+1
      else
        pe_descr(XYZ,ZHI) = pe_descr(ZLR,PLSZ)
      endif
      if(pe_descr(BND,ZBTM).eq.YES)then
        pe_descr(XYZ,ZLO) = 0
      else
        pe_descr(XYZ,ZLO) = 1
      endif
C###{__Limits_[X]__}
      if(pe_descr(BND,EAST).eq.YES)then
        pe_descr(XYZ,XHI) = pe_descr(COL,PLSZ)+1
      else
        pe_descr(XYZ,XHI) = pe_descr(COL,PLSZ)
      endif
      if(pe_descr(BND,WEST).eq.YES)then
        pe_descr(XYZ,XLO) = 0
      else
        pe_descr(XYZ,XLO) = 1
      endif
C###{__Limits_[Y]__}
      if(pe_descr(BND,NORTH).eq.YES)then
        pe_descr(XYZ,YHI) = pe_descr(ROW,PLSZ)+1
      else
        pe_descr(XYZ,YHI) = pe_descr(ROW,PLSZ)
      endif
      if(pe_descr(BND,SOUTH).eq.YES)then
        pe_descr(XYZ,YLO) = 0
      else
        pe_descr(XYZ,YLO) = 1
      endif
C###{__done__}

C###{__Multigrid_parameters__}
C###{__by_calling_mglev()_we_good_estimation_for__}
C###{__the_total_number_of_levels__}
      call mglev(nglx,ngly,mlev)
      mlev = max(m2,max(mlev,m1))
      m1 = mlev
      m2 = mlev
      pe_descr(MGD,MLLT) = mlev
      pe_descr(MGD,PLLE) = mlev
      pe_descr(MGD,LCUR) = mlev

C###{__Estimate_PLLE_and_[nall2]_vector_sizes__}
      call est_arrs_szs(pe_descr,mlev,nxall,nyall,nbf)
      if(pe_id.eq.ROOT)then
        write(*,200)pe_descr(MGD,MLLT),pe_descr(MGD,PLLE)
      endif

C###{__Set_Up_The_logger_}
      logfl = log_file_name(pe_id)
      call initsink(logfl)

C###{__Show_some_info_}
#ifdef _LOGS_
      call pe_descr_show(pe_descr)
#endif
      row_bsz = 2 + pe_descr(ROW,PLSZ) !!->ngy
      col_bsz = 2 + pe_descr(COL,PLSZ) !!->ngy
      zlr_bsz = 2 + pe_descr(ZLR,PLSZ) !!->ngz

C###{__Call_memory_allocation_C_routine__}
C###{__Note :: C is cool }
      call memsetf(pe_descr,pe_id,pe_size,nbf,
     &  row_bsz,col_bsz,zlr_bsz,
     &  ierr,
     &  nxall,nyall,
     &  omg,ngx,ngy,ngz,
     &  npx,npy,npz,
     &  bet1,bet2,
     &  idata,irhs,m1,m2,nr1,nr2,nr3,nr4,
     &  nit,icyc,md,
     &  omega,rpm,
     &  ids,iswp,itsol,itmax,itol,
     &  ipm,ia,ja,ip)
      
C###{ I am BACK ...}
C###{__Stop_timer__}
      tt = timing(tt)
      if(pe_id.eq.ROOT)then
        write(*,'("[Total_Time]:=",D12.6)')tt
        if(ierr.eq.0)print*,":safe landing at:"
        call tmpclean()
      endif
C###{__This is the end, my friend, the end ...}
      call MPI_Finalize(ierr)
C###{__done__}
 100  format("[",I3,"]=>procs_in_use: [npy=",I2,"; npx=",I2,"]")
 200  format("[MG_Levels]: tot=",I2," pllst=",I2)
      END PROGRAM
C###{__EOF__}
