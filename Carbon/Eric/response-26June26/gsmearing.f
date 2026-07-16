C=======================================================================
                                                                        
      SUBROUTINE GSMEARING(Z, A, W2, Q2, xvalc, F1, F2, FL)

CCC   Returns Fermi smeared structure functions.  Smearing function is a Gaussian.  CCC
CCC   Note:  not tested for nuclei with A < 12 or A > 64                            CCC      
CCC   Oct. 4, 2025                                                                  CCC                      
!--------------------------------------------------------------------

      implicit none
      real*8 Z,A,q2,w2,f1,f2,fL,xvalc(100),w2t,kappa2
      real*8 nu,x,mp,mp2,mpi,pi,f1p,f1pp,f1dp,f2p,f2pp,fLp,fLpp
      real*8 f1n,f1nn,f2n,f2nn,fLn,fLnn,f1d,offshell,delta,qvmax
      real*8 pf,pf2,kf,qv,es,dw2des,fyuse,fyuse2,epr,kbar2,fcor,deltae
      real*8 epr2,wsqp,wsqp2,frac2b,fracs,xt,xt2,xp,rc,emct,off_mKP_fit
      real*8 dw2dpf,r,zt,at,Lcor,frac,fract,t1,t2,q20,C1,C2,dES,a0
      real*8 xxp(500),fytot,fytot2,norm,bw,nwid,ncor,D1,p2
      
      real*8 emcfac,emcfacL,qvt,exmin,ex,nuel,cof
      logical goodfit
      INTEGER ISM,j,nbins

      nbins = 102
      nwid = 3.75
      bw = 2.0*nwid/float(nbins)

      C1 = 0.3
      q20 = 0.01
      
      C1 = xvalc(14)*(A/12.0)**0.25
      C1 = xvalc(14)
c      C1 = 0.0
c      q20 = xvalc(16)
c      C2 = (1.0+C1*q20/(q2+q20))**xvalc(39)
c      t2 = 1.0-xvalc(15)*(1.0+0.02*A**0.25)/(1.0+q2/0.8**2)**2
c     a0 = (1.0+xvalc(39)*(A/12.0-1.0)**0.25)
c      a0 = 1.0+C1
c      a0 = 1.0+xvalc(39)*(sqrt(A/12.0)-1.0)     

      a0 = 1.0D0
      
c      write(6,*) A, a0
      
      if(A.EQ.12) exmin = 0.0165
      If(A.EQ.27) exmin = 0.0085
      If(A.EQ.40) exmin = 0.0085
      
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

      qvmax = 1.5
      if(A.GE.3.) then
c      !!! energy shift !!!
c        Es = 0.0085
      !!! fermi momentum  !!!
        kf = xvalc(25)+0.00055*(A/12.0*sqrt(A/12.0)-1.0)

c        if(A.EQ.12) then
c           kf = 0.2154
c        elseif(A.GT.12) then
c           kf = 0.2375emacs
c        endif
        
c     kf = xvalc(25)*(1.0+0.001*A**0.25)
c        kf = kf*(1.0+0.00015*(A/12.0-1.0)**0.25)

c        qvt = qv
c     if(qv.GE.qvmax) qvt = qvmax
        Es = exmin
c        if(qv.LT.qvmax) then
c     Es = exmin + 2./3.*xvalc(24)  !!! Maximum Es
c          Es = Es-2./3.*xvalc(24)*(1.0-(qv/qvmax)**2)**0.5 !!!  test  !!!           
c          Es = Es-2./3.*xvalc(24)*(1.0-qv/qvmax)**xvalc(26) !!!  test  !!!
c        endif
      endif

      norm = sqrt(pi)
      ncor = 1.000
      ncor = 1.0027
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
        fyuse = fyuse*(1.0-xvalc(22)+xvalc(22)*xxp(ism)*xxp(ism))
        
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
        
c        t1 =   (1.0+xvalc(37)*A*A)**xvalc(38)
c        t1 = 1.0+xvalc(37)*(A/100.0)**xvalc(38)
c        t1 = xvalc(37)+xvalc(38)*(sqrt(A/12.0)-1.0)
c        t1 = sqrt(t1) 

        
        do j=1,1
           if(j.EQ.1) then
              fract = 1.0D0-frac
              w2t = WSQP
           else
              fract = frac
              w2t = WSQP2
           endif   
           IF(w2t.GT.1.159) THEN
             xt = q2/(q2+W2t-mp2)
             xp = 1.0D0+(w2t-mp2)/(q2+xvalc(36)/
     &            (1.0+xvalc(39)*q2**2)) !!! fixed on 11/02/23, mp -> mp2
c             xp = 1.0D0+(w2t-mp2)/(q2+xvalc(36)) !!! fixed on 11/02/23, mp -> mp2
             
             xp = 1.0D0/xp
c     xp = w2t-mp2+q2/(q2+xvalc(34))
c             xp = (q2+xvalc(36))/(q2+w2t-mp2)
c             xp = 1.0D0/xp
          
             offshell = 1.0D0      !!!  test

             xt2 = xt
             xt2 = xt
             xt2 = max(xt2,0.0001)

CCC   Next is medium modification factor  CCC
             
c            emcfac = 1.0-xvalc(12)/a0*abs(1.-xp)**xvalc(15)
c             emcfac = 1.0
c             emcfac = xvalc(11)*emcfac/(1.0-xvalc(12)/a0*
c     &            abs(1.-xvalc(13))**3)/
c     &            (1.0+xvalc(14)*abs(xp-xvalc(13))**1)

c             emcfac = emcfac*(1.0+xvalc(16)*(xt-xvalc(13)))

c             emcfac = emcfac*(1.0+xvalc(37)*(1.0-xt)**1.5+
c    &                     xvalc(38)*(1.0-xt)**2.5)

c             emcfac = 1.01-(0.65*(xt-0.25)**2.0)**0.75
c             emcfac = emcfac/(1.5-xt)**0.1
c             emcfac = emcfac*1.00*(1.0+0.15*xt**0.5)

             p2 = xvalc(21)-1.0
c             -xvalc(14)/(1.0-xvalc(12)**2)**xvalc(15)
c             write(6,*) xvalc(21),xvalc(12),p2
             emcfac = xvalc(21)
     &            -xvalc(11)*((xt2-xvalc(12))**2.0)**xvalc(13)
     &            +xvalc(14)*(1.0-xp*xp)**xvalc(15)


             emcfac = emcfac/(1.0+0.05*xp)**xvalc(40)
             
c            if(abs(xt2-xvalc(12)).LT.0.02) write(6,*) emcfac, p2 


c            emcfac = emcfac*(1.0+xvalc(14)*xt2**xvalc(15))**xvalc(16)

c              if(q2.EQ.0) emcfac = 0.89
c              emcfac = 1.01-(1.4*(xt-0.23)**2.0)**0.75 
c              emcfac = emcfac*(1.0+3.5*xt**5.0)**0.75
             
c             emcfac = max(emcfac,0.85)
             
c            emcfacL = xvalc(17)*(1.0D0+xvalc(18)*xp*xp)*
c     &       (1.+xvalc(19)*xp*xp*xp*xp)*exp(-1.0*xvalc(20)*xp)
            
c            emcfacL = sqrt(C2)*emcfacL
            emcfacL = 1.0
            
            call sf(w2t,q2,f1pp,fLpp,f2pp,f1nn,fLnn,f2nn)

c            t1 = (1.0+4.*xt*xt*mp2/q2)*f2pp-(2*xt*f1pp)

c            write(6,*) xt,q2,fLpp,t1         !!!  TEST

c            if(q2.LT.0.001) write(6,*) "TEST:  ",w2t,q2,xp,xt
            
            f1pp = f1pp*emcfac*offshell
            f1nn = f1nn*emcfac*offshell
            fLpp = fLpp*emcfac*emcfacL*offshell
            fLnn = fLnn*emcfac*emcfacL*offshell
c            fLpp = fLpp*emcfacL*offshell
c            fLnn = fLnn*emcfacL*offshell

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

      F1 = (Z*F1p+(A-Z)*F1n)
      F2 = (Z*F2p+(A-Z)*F2n)
      FL = (Z*FLp+(A-Z)*FLn)

      cof = (1.0-(es/nu)**4)
      cof = min(1.0,cof)
      cof = max(cof,0.0)
      if(nu.LE.es) cof = 0.0
      F1 = F1*cof
      F2 = F2*cof
      fL = FL*cof

      F2 = (2.0*x*F1+FL)/kappa2    !!! *** Check for consistency !!!
      
c      if(q2.LT.0.05)  write(6,*) "gsmearing:  ", w2,x,q2,F1   
      
c      if(q2.GT.0) F1 = (kappa2*F2-FL)/2.0/x  !!!  Calculate to keep internal consistency  !!!

c      if(q2.LT.0.05) write(6,*) "gsmearing2:  ", w2,x,q2,F1      

      
      if(F1.LT.0.0) F1 = 0.0
      if(F2.LT.0.0) F2 = 0.0
      if(FL.LT.0.0) FL = 0.0

      
c      write(6,*) w2,f1p,f1n
      
c      write(6,*) fytot,fytot2

 2000 format(2f7.3,1i4,4f10.4)



      RETURN                                                            
      END                                          
