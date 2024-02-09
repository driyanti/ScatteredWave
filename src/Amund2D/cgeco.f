      subroutine  cscal(n,ca,cx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, linpack,  3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex ca,cx(*)
      integer i,incx,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
cccccccccccccccccccccccccccccccccccccc
      subroutine  csscal(n,sa,cx,incx)
c
c     scales a complex vector by a real constant.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex cx(*)
      real sa
      integer i,incx,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   30 continue
      return
      end
cccccccccccccccccccccccccccccc
      integer function icamax(n,cx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex cx(*)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end
ccccccccccccccccccccccc
      real function scasum(n,cx,incx)
c
c     takes the sum of the absolute values of a complex vector and
c     returns a single precision result.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      complex cx(*)
      real stemp
      integer i,incx,n,nincx
c
      scasum = 0.0e0
      stemp = 0.0e0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   10 continue
      scasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   30 continue
      scasum = stemp
      return
      end
ccccccccccccccccccccccccccc
      subroutine cgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      complex a(lda,1),z(1)
      real rcond
c
c     cgeco factors a complex matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, cgefa is slightly faster.
c     to solve  a*x = b , follow cgeco by cgesl.
c     to compute  inverse(a)*c , follow cgeco by cgesl.
c     to compute  determinant(a) , follow cgeco by cgedi.
c     to compute  inverse(a) , follow cgeco by cgedi.
c
c     on entry
c
c        a       complex(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   real
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       complex(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack cgefa
c     blas caxpy,cdotc,csscal,scasum
c     fortran abs,aimag,amax1,cmplx,conjg,real
c
c     internal variables
c
      complex cdotc,ek,t,wk,wkm
      real anorm,s,scasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
      complex zdum,zdum1,zdum2,csign1
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
      csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
      anorm = 0.0e0
      do 10 j = 1, n
         anorm = amax1(anorm,scasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call cgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e .
c     ctrans(a)  is the conjugate transpose of a .
c     the components of  e  are chosen to cause maximum local
c     growth in the elements of w  where  ctrans(u)*w = e .
c     the vectors are frequently rescaled to avoid overflow.
c
c     solve ctrans(u)*w = e
c
      ek = (1.0e0,0.0e0)
      do 20 j = 1, n
         z(j) = (0.0e0,0.0e0)
   20 continue
      do 100 k = 1, n
         if (cabs1(z(k)) .ne. 0.0e0) ek = csign1(ek,-z(k))
         if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
            s = cabs1(a(k,k))/cabs1(ek-z(k))
            call csscal(n,s,z,1)
            ek = cmplx(s,0.0e0)*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = cabs1(wk)
         sm = cabs1(wkm)
         if (cabs1(a(k,k)) .eq. 0.0e0) go to 40
            wk = wk/conjg(a(k,k))
            wkm = wkm/conjg(a(k,k))
         go to 50
   40    continue
            wk = (1.0e0,0.0e0)
            wkm = (1.0e0,0.0e0)
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + cabs1(z(j)+wkm*conjg(a(k,j)))
               z(j) = z(j) + wk*conjg(a(k,j))
               s = s + cabs1(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*conjg(a(k,j))
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + cdotc(n-k,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 110
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
c
      ynorm = 1.0e0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call caxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (cabs1(z(k)) .le. 1.0e0) go to 130
            s = 1.0e0/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
            s = cabs1(a(k,k))/cabs1(z(k))
            call csscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (cabs1(a(k,k)) .ne. 0.0e0) z(k) = z(k)/a(k,k)
         if (cabs1(a(k,k)) .eq. 0.0e0) z(k) = (1.0e0,0.0e0)
         t = -z(k)
         call caxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0e0/scasum(n,z,1)
      call csscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
      if (anorm .eq. 0.0e0) rcond = 0.0e0
      return
      end
cccccccccccccccc
      subroutine cgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex a(lda,1)
c
c     cgefa factors a complex matrix by gaussian elimination.
c
c     cgefa is usually called by cgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for cgeco) = (1 + 9/n)*(time for cgefa) .
c
c     on entry
c
c        a       complex(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that cgesl or cgedi will divide by zero
c                     if called.  use  rcond  in cgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas caxpy,cscal,icamax
c     fortran abs,aimag,real
c
c     internal variables
c
      complex t
      integer icamax,j,k,kp1,l,nm1
c
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = icamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0e0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0e0,0.0e0)/a(k,k)
            call cscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call caxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0e0) info = n
      return
      end


