      integer  nk_total,nk_totpos,srclayer,rcvlayer,nlay
      logical    freqtest,test,logicread
  
      real*8   deltap,
     .	       tot_inc_pmax,tot_inc_pmin,tot_inc_deltap

      real*8   rho(maxlay),D(maxlay),z(0:maxlay),h(maxlay)
            
      complex*16  G(maxlay),K(maxlay),alpha(maxlay),beta(maxlay)
      complex*16   qalpha(maxlay),qbeta(maxlay)
     
     
      complex*16 matRd(3,3,0:maxlay),matTd(3,3,0:maxlay),
     .		matRu(3,3,0:maxlay),matTu(3,3,0:maxlay)
      

      common/soil/D,rho,G,K,z,h,alpha,beta
      common/soilextra/qalpha,qbeta
      common/spacesmpl/deltap,nk_total,nk_totpos,
     .	    	  	tot_inc_pmax,tot_inc_pmin,tot_inc_deltap
      common/layers/nlay,srclayer,rcvlayer
      common/transmat/matRd,matTd,matRu,matTu
      common/testslowness/freqtest,test,logicread
