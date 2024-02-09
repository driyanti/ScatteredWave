	double precision function dzasum(n,zx,incx)
c
c       takes the sum of the absolute values.
c       jack dongarra, 3/11/78.
c       modified 3/93 to return if incx .le. 0.
c       modified 12/3/93, array(1) declarations changed to array(*)
c
	double complex zx(*)
	
	double precision stemp,dcabs1
	integer i,incx,ix,n
c
	dzasum = 0.0d0
	stemp = 0.0d0
	if( n.le.0 .or. incx.le.0 )return
	if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
	ix = 1
	do 10 i = 1,n
	stemp = stemp + dcabs1(zx(ix))
	ix = ix + incx
10 	continue
	dzasum = stemp
	return
c
c        code for increment equal to 1
c
20 	do 30 i = 1,n
	stemp = stemp + dcabs1(zx(i))
30 	continue
	dzasum = stemp
	return
	end


	integer function izamax(n,zx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, 1/15/85.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
	double complex zx(*)
	double precision smax
	integer i,incx,ix,n
	double precision dcabs1
c
	izamax = 0
	if( n.lt.1 .or. incx.le.0 )return
	izamax = 1
	if(n.eq.1)return
	if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
	ix = 1
	smax = dcabs1(zx(1))
	ix = ix + incx
	do 10 i = 2,n
	if(dcabs1(zx(ix)).le.smax) go to 5
	izamax = i
	smax = dcabs1(zx(ix))
5    	ix = ix + incx
10 	continue
	return
c
c        code for increment equal to 1
c
20 	smax = dcabs1(zx(1))
	do 30 i = 2,n
	if(dcabs1(zx(i)).le.smax) go to 30
	izamax = i
	smax = dcabs1(zx(i))
30 	continue
	return
	end


	subroutine  zdscal(n,da,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
	double complex zx(*)
	double precision da
	integer i,incx,ix,n
c
	if( n.le.0 .or. incx.le.0 )return
	if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
	ix = 1
	do 10 i = 1,n
	zx(ix) = dcmplx(da,0.0d0)*zx(ix)
	ix = ix + incx
10 	continue
	return
c
c        code for increment equal to 1
c
20 	do 30 i = 1,n
	zx(i) = dcmplx(da,0.0d0)*zx(i)
30 	continue
	return
	end

	subroutine  zscal(n,za,zx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
	double complex za,zx(*)
	integer i,incx,ix,n
c
	if( n.le.0 .or. incx.le.0 )return
	if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
	ix = 1
	do 10 i = 1,n
	zx(ix) = za*zx(ix)
	ix = ix + incx
10 	continue
	return
c
c        code for increment equal to 1
c
20 	do 30 i = 1,n
	zx(i) = za*zx(i)
30 	continue
	return
	end

	subroutine zgeco(a,lda,n,ipvt,rcond,z)

	integer lda,n,ipvt(1)
	complex*16 a(lda,1),z(1)
	double precision rcond
c
c     zgeco factors a complex*16 matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, zgefa is slightly faster.
c     to solve  a*x = b , follow zgeco by zgesl.
c     to compute  inverse(a)*c , follow zgeco by zgesl.
c     to compute  determinant(a) , follow zgeco by zgedi.
c     to compute  inverse(a) , follow zgeco by zgedi.
c
c     on entry
c
c        a       complex*16(lda, n)
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
c        rcond   double precision
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
c        z       complex*16(n)
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
c     linpack zgefa
c     blas zaxpy,zdotc,zdscal,dzasum
c     fortran dabs,dmax1,dcmplx,dconjg
c
c     internal variables
c
	complex*16 zdotc,ek,t,wk,wkm
	double precision anorm,s,dzasum,sm,ynorm
	integer info,j,k,kb,kp1,l
c
	complex*16 zdum,zdum1,zdum2,csign1
	double precision cabs1
	double precision dreal,dimag
	complex*16 zdumr,zdumi
	dreal(zdumr) = zdumr
	dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
	cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
	csign1(zdum1,zdum2) = cabs1(zdum1)*(zdum2/cabs1(zdum2))
c
c     compute 1-norm of a
c
	anorm = 0.0d0
	do 10 j = 1, n
	anorm = dmax1(anorm,dzasum(n,a(1,j),1))
10 	continue
c
c     factor
c
	call zgefa(a,lda,n,ipvt,info)
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
	ek = (1.0d0,0.0d0)
	do 20 j = 1, n
	z(j) = (0.0d0,0.0d0)
20 	continue
	do 100 k = 1, n
	if (cabs1(z(k)) .ne. 0.0d0) ek = csign1(ek,-z(k))
	if (cabs1(ek-z(k)) .le. cabs1(a(k,k))) go to 30
	s = cabs1(a(k,k))/cabs1(ek-z(k))
	call zdscal(n,s,z,1)
	ek = dcmplx(s,0.0d0)*ek
30    continue
	wk = ek - z(k)
	wkm = -ek - z(k)
	s = cabs1(wk)
	sm = cabs1(wkm)
	if (cabs1(a(k,k)) .eq. 0.0d0) go to 40
	wk = wk/dconjg(a(k,k))
	wkm = wkm/dconjg(a(k,k))
	go to 50
40    continue
	wk = (1.0d0,0.0d0)
	wkm = (1.0d0,0.0d0)
50    continue
	kp1 = k + 1
	if (kp1 .gt. n) go to 90
	do 60 j = kp1, n
	sm = sm + cabs1(z(j)+wkm*dconjg(a(k,j)))
	z(j) = z(j) + wk*dconjg(a(k,j))
	s = s + cabs1(z(j))
60       continue
	if (s .ge. sm) go to 80
	t = wkm - wk
	wk = wkm
	do 70 j = kp1, n
	z(j) = z(j) + t*dconjg(a(k,j))
70          continue
80       continue
90    continue
	z(k) = wk
100 	continue
	s = 1.0d0/dzasum(n,z,1)
	call zdscal(n,s,z,1)
c
c     solve ctrans(l)*y = w
c
	do 120 kb = 1, n
	k = n + 1 - kb
	if (k .lt. n) z(k) = z(k) + zdotc(n-k,a(k+1,k),1,z(k+1),1)
	if (cabs1(z(k)) .le. 1.0d0) go to 110
	s = 1.0d0/cabs1(z(k))
	call zdscal(n,s,z,1)
110    continue
	l = ipvt(k)
	t = z(l)
	z(l) = z(k)
	z(k) = t
120 	continue
	s = 1.0d0/dzasum(n,z,1)
	call zdscal(n,s,z,1)
c
	ynorm = 1.0d0
c
c     solve l*v = y
c
	do 140 k = 1, n
	l = ipvt(k)
	t = z(l)
	z(l) = z(k)
	z(k) = t
	if (k .lt. n) call zaxpy(n-k,t,a(k+1,k),1,z(k+1),1)
	if (cabs1(z(k)) .le. 1.0d0) go to 130
	s = 1.0d0/cabs1(z(k))
	call zdscal(n,s,z,1)
	ynorm = s*ynorm
130    continue
140 	continue
	s = 1.0d0/dzasum(n,z,1)
	call zdscal(n,s,z,1)
	ynorm = s*ynorm
c
c     solve  u*z = v
c
	do 160 kb = 1, n
	k = n + 1 - kb
	if (cabs1(z(k)) .le. cabs1(a(k,k))) go to 150
	s = cabs1(a(k,k))/cabs1(z(k))
	call zdscal(n,s,z,1)
	ynorm = s*ynorm
150    continue
	if (cabs1(a(k,k)) .ne. 0.0d0) z(k) = z(k)/a(k,k)
	if (cabs1(a(k,k)) .eq. 0.0d0) z(k) = (1.0d0,0.0d0)
	t = -z(k)
	call zaxpy(k-1,t,a(1,k),1,z(1),1)
160	 continue
c     make znorm = 1.0
	s = 1.0d0/dzasum(n,z,1)
	call zdscal(n,s,z,1)
	ynorm = s*ynorm
c
	if (anorm .ne. 0.0d0) rcond = ynorm/anorm
	if (anorm .eq. 0.0d0) rcond = 0.0d0
	return
	end

	subroutine zgefa(a,lda,n,ipvt,info)
	integer lda,n,ipvt(1),info
	complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
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
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
	complex*16 t
	integer izamax,j,k,kp1,l,nm1
c
	complex*16 zdum
	double precision cabs1
	double precision dreal,dimag
	complex*16 zdumr,zdumi
	dreal(zdumr) = zdumr
	dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
	cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
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
	l = izamax(n-k+1,a(k,k),1) + k - 1
	ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
	if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
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
	t = -(1.0d0,0.0d0)/a(k,k)
	call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
	do 30 j = kp1, n
	t = a(l,j)
	if (l .eq. k) go to 20
	a(l,j) = a(k,j)
	a(k,j) = t
20          continue
	call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
30       continue
	go to 50
40    continue
	info = k
50    continue
60	 continue
70 	continue
	ipvt(n) = n
	if (cabs1(a(n,n)) .eq. 0.0d0) info = n
	return
	end



