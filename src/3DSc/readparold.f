      SUBROUTINE readpar
c
c author: Auke Ditzel
c date:   29-08-2000
c last change 12-08-2001
c
c to read the right parameters from file
c
      implicit none
      
      INCLUDE 'bounds.h'
      INCLUDE 'num.h'
      INCLUDE 'param.h'
      INCLUDE 'record.h'
      INCLUDE 'path.h'
      
      integer   ii,writechar
      real      rdummy
      integer   idummy
      real*8    oscildummy,dummy
      
      common/schrijf/writechar
      
      pi = 4.d0*atan(1.d0)
c
c time sampling papameters
c
      nsamp= 2048
      dt   = 3.e-3
      
      idummy = nsamp
      call getparint(iopath//'params.in','nsamp',idummy)
      nsamp = idummy
      
      rdummy = dt
      call getparfloat(iopath//'params.in','dt',rdummy)
      dt = dble(rdummy)
c
c spatial samplig parameters
c
      
      tot_inc_pmin   = 0.4d0
      tot_inc_pmax   = 0.04d0
      tot_inc_deltap = 4e-3

      rdummy = tot_inc_pmin
      call getparfloat
     .	    	 (iopath//'params.in','tot_inc_pmin',rdummy)
      tot_inc_pmin=dble(rdummy)

      rdummy = tot_inc_pmax
      call getparfloat
     .	    	 (iopath//'params.in','tot_inc_pmax',rdummy)
      tot_inc_pmax=dble(rdummy)
      
      rdummy = tot_inc_deltap
      call getparfloat
     .	    	 (iopath//'params.in','tot_inc_deltap',rdummy)
      tot_inc_deltap=dble(rdummy)
c
c train characteristics
c
      nroscils     = 1
      oscilfreq(1) = 0.d0
      OPEN (unit=12,file=iopath//'oscilfreqs.in')
      READ(12,*) nroscils
      IF (nroscils.gt.maxoscils) STOP 'number of oscils too big'
      DO ii = 1,nroscils      
        read(12,*) oscildummy
	  oscilfreq(ii) = oscildummy*2.*pi
      END DO
      CLOSE(12)
      write(*,*) 'variable nroscils, value ', nroscils
      write(*,*) 'variable oscilfreq, values ',(oscilfreq(ii)/(2.*pi),
     .						ii=1,nroscils), 'Hz'
c
c soil parameters
c      
      nlay  = 1
      rho(1)= 1500.d0
      D(1)  = 2.00d0
      G(1)  = cmplx(5.40e6,0.d0)
      K(1)  = cmplx(58.25e6,0.d0)
      h(1)  = 5.d0

      rho(2) = 1500.d0
      D(2)   = 2.00d0
      G(2)   = cmplx(5.40e6,0.d0)
      K(2)   = cmplx(58.25e6,0.d0)
      
      OPEN (unit=12,file=iopath//'soil.in')
      READ(12,*) nlay
      IF (nlay+1.gt.maxlay) STOP 'Too many layers'
      DO 100 ii=1,nlay
        READ(12,*)
        READ(12,*) dummy
        G(ii) = cmplx(dummy,0.d0)
        READ(12,*) dummy
        K(ii) = cmplx(dummy,0.d0)
        READ(12,*) D(ii)
        READ(12,*) rho(ii)
        READ(12,*) h(ii)
100   CONTINUE
      READ(12,*)
      READ(12,*) dummy
      G(nlay+1) = cmplx(dble(dummy),0.)
      READ(12,*) dummy
      K(nlay+1) = cmplx(dble(dummy),0.)
      READ(12,*) D(nlay+1)
      READ(12,*) rho(nlay+1)
      CLOSE(12)
      
      velo_src = 25.d0

      rdummy   = velo_src
      call getparfloat(iopath//'params.in','velo_src',rdummy)
      velo_src=dble(rdummy)
c
c receiver characteristics
c
      fmin      = 5.d0
      freqstart = 20.d0
      freqend   = 40.d0
      fmax      = 55.d0

      rdummy=fmin
      call getparfloat(iopath//'params.in','freqmin',rdummy)
      fmin=dble(rdummy)
      
      rdummy=freqstart
      call getparfloat(iopath//'params.in','freqstart',rdummy)
      freqstart=dble(rdummy)
      
      rdummy=freqend
      call getparfloat(iopath//'params.in','freqend',rdummy)
      freqend=dble(rdummy)
      
      rdummy=fmax
      call getparfloat(iopath//'params.in','freqmax',rdummy)
      fmax=dble(rdummy)
c
c source position
c
      srcpos(1) = 0.d0
      srcpos(2) = 0.d0
      srcpos(3) = 0.d0      

        rdummy=srcpos(1)
      call getparfloat(iopath//'params.in','xpos_src',rdummy)
	srcpos(1)=dble(rdummy)

        rdummy=srcpos(2)
      call getparfloat(iopath//'params.in','ypos_src',rdummy)
	srcpos(2)=dble(rdummy)

        rdummy=srcpos(3)
      call getparfloat(iopath//'params.in','zpos_src',rdummy)
	srcpos(3)=dble(rdummy)
c
c direction line receivers (perpendicular or parallel)
c
      grid = .false.
      perpendicular = .true.

      call getparbool(iopath//'params.in','grid',grid)
      if (.not.grid) then
        perpendicular = .true.
        call getparbool(iopath//'params.in',
     .			'perpendicular',perpendicular)
        if (.not.perpendicular) then
          perpendicular = .false.
        end if
      end if
c
c number of receivers, position first receiver and spacing
c      
      rcvpos(1) = 5.d0
      rcvpos(2) = 5.d0
      rcvpos(3) = 2.d0
      rcv_spac  = 2.5d0
      nr_rcv    = 10

      rdummy=rcvpos(1)
      call getparfloat(iopath//'params.in','xrcvpos',rdummy)
      rcvpos(1)=dble(rdummy)

      rdummy=rcvpos(2)
      call getparfloat(iopath//'params.in','yrcvpos',rdummy)
      rcvpos(2)=dble(rdummy)

      rdummy=rcvpos(3)
      call getparfloat(iopath//'params.in','zrcvpos',rdummy)
      rcvpos(3)=dble(rdummy)

      rdummy=rcv_spac
      call getparfloat(iopath//'params.in','rcv_spac',rdummy)
      rcv_spac=dble(rdummy)

      idummy=nr_rcv
      call getparint(iopath//'params.in','nr_rcv',idummy)
      nr_rcv=int(idummy)

      test = .false.
      call getparbool(iopath//'params.in','test',test)
      
      freqtest = .false.
      call getparbool(iopath//'params.in','freqtest',freqtest)
      
      logicread = .false.
      call getparbool(iopath//'params.in','logicread',logicread)

      IF (nsamp.gt.maxtsamp) 
     1		STOP 'Too many samples. Adjust include or input file.'
      IF (nr_rcv.gt.max_rcv) 
     1		STOP 'Too many receivers. Adjust include or input file.'

      END SUBROUTINE
