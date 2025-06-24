C=======================================================================
                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, xvalc, F1, F2, FL)

CCC   Returns Fermi smeared structure functions.  Smearing function is a Gaussian.  CCC
CCC   Note:  not tested for nuclei with A < 12 or A > 64                            CCC      
CCC   July 20, 2020                                                                 CCC                      
!--------------------------------------------------------------------

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,xvalc(45),w2t,kappa2
      real*8 nu,x,mp,mp2,mpi,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta
      real*8 pf,pf2,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,xp,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at,Lcor,frac,fract,t1
      real*8 xxp(500),fytot,fytot2,norm,bw,nwid,ncor
      
      real*8 emcfac,emcfacL,qvt,exmin,ex,nuel,cof
      logical goodfit
      INTEGER ISM,j,nbins

      nbins = 98
      nwid = 3.3
      bw = 2.0*nwid/float(nbins)
           
      exmin = 0.0165
      mp = 0.938272
      mp2 = mp*mp
      mpi = 0.135
      pi = 3.141593
      x = q2/(q2+w2-mp2)
      nu = (w2+q2-mp2)/2./mp      
      qv = sqrt(nu**2 + q2)
      kappa2 = 1.0+4.0*mp2*x*x/q2
      nuel = q2/2./(0.931494*A)
      ex = nu-nuel

      
      if(A.GE.3.) then
      !!! energy shift !!!
        Es = 0.008
      !!! fermi momentum  !!!
        kf = xvalc(37)
        qvt = qv
        if(qv.GE.1.0) qvt = 1.0
        Es = xvalc(38)
        Es = Es-xvalc(39)*(1.0-qvt)   !!!  test  !!!
      endif

      norm = sqrt(pi)
      ncor = 1.000
c      ncor = 1.0027
c      ncor = 1.0663+4.7*(kf-0.225)
      norm = norm/ncor   !!! account for missing part of distribution past nwid for kf = 225 MeV


      f1p = 0.0D0
      f1n = 0.0D0
      f2p = 0.0D0
      f2n = 0.0D0
      fLp = 0.0D0
      fLn = 0.0D0


c        sigt = 0.0D0
c        sigL = 0.0D0
      fytot = 0.0D0
      fytot2 = 0.0D0


! adjust pf to give right width based on kf
      pf = 0.5 * kf
      pf2 = pf*1.5
! assume this is 2 * pf * qv
      DW2DPF = 2. * qv
c       DW2DPF = qv    
      dw2des = 2. * (nu + mp) 

      DO ism = 1,nbins

CCC   
        xxp(ism) = -nwid+bw*(float(ism-1))

        
        fyuse = bw/sqrt(2.0)/norm*exp(-0.5*xxp(ism)*xxp(ism)) !!! Gaussian !!!       
        
CCC  Next is from f1f209 CCC

        WSQP = W2 + XXp(ISM) * PF * DW2DPF - es * dw2des
        WSQP2 = W2 + XXp(ISM) * PF2 * DW2DPF - es * dw2des
        
CCC

        fytot = fytot+fyuse
        fytot2 = fytot2+fyuse

c           write(6,2000) w2,q2,ism,xxp(ism),fyuse, wsqp, fytot

        F1pp = 0.0D0
        F1nn = 0.0D0
        F2pp = 0.0D0
        F2nn = 0.0D0
        FLpp = 0.0D0
        FLnn = 0.0D0

        frac = 0.0
        
        do j=1,1
           if(j.EQ.1) then
              fract = 1.0D0-frac
              w2t = WSQP
           else
              fract = frac
              w2t = WSQP2
           endif   
           IF(w2t.GT. 1.159) THEN
             xt = q2/(q2+W2t-mp2)
             xp = 1.0D0+(w2t-mp)/(q2+xvalc(34))
             xp = 1.0D0/xp
          
             offshell = 1.0D0      !!!  test


CCC   Next is medium modification factor  CCC

            emcfac = (xvalc(26)+xvalc(27)*xp*xp)/
     &           (1.0+xvalc(28)*xp+xvalc(29)*xp*xp)

            emcfacL = 1.0     
            emcfacL = xvalc(30)*(1.0D0+xvalc(31)*xp*xp)*
     &       (1.+xvalc(32)*xp*xp)*exp(-1.0*xvalc(33)*xp)


            call sf(w2t,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)

c            t1 = (1.0+4.*xt*xt*mp2/q2)*f2pp-(2*xt*f1pp)

c            write(6,*) xt,q2,fLpp,t1         !!!  TEST
            
            f1pp = f1pp*emcfac*offshell
            f1nn = f1nn*emcfac*offshell
            fLpp = fLpp*emcfac*emcfacL*offshell
            fLnn = fLnn*emcfac*emcfacL*offshell
            f2pp = (2.*xt*f1pp+fLpp)/(1.+4.*xt*xt*mp2/q2)
            f2nn = (2.*xt*f1nn+fLnn)/(1.+4.*xt*xt*mp2/q2)

            F1p = F1p + F1pp * Fyuse * fract
            F1n = F1n + F1nn * Fyuse * fract
            F2p = F2p + F2pp * Fyuse * fract
            F2n = F2n + F2nn * Fyuse * fract      
            FLp = FLp + FLpp * Fyuse * fract
            FLn = FLn + FLnn * Fyuse * fract

         ENDIF
       ENDDO
      ENDDO

c      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

      cof = (ex-exmin)**0.5/(0.025-exmin)**0.5
      cof = min(1.0,cof)
      cof = max(cof,0.0)
      if(ex.LE.exmin) cof = 0.0
      F2 = F2*cof
      fL = FL*cof
      
      F1 = (kappa2*F2-FL)/2.0/x  !!!  Calculate to keep internal consistency  !!!


      
      if(F1.LT.0.0) F1 = 0.0
      if(F2.LT.0.0) F2 = 0.0
      if(FL.LT.0.0) FL = 0.0

      
c      write(6,*) w2,f1p,f1n
      
c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          
