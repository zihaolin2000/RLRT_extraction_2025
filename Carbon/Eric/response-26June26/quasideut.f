      SUBROUTINE QUASIDEUT(z,a,w2,q2,xvalm,f1qd)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCC   Subroutine for Transverse Enhancement in the QE and Delta region.         CCC
CCC   exchange currents and isobar excitations in the medium.  This is assumed  CCC
CCC   to be due to quasi-deuteron 2-body currents.  Shape is a distorted        CCC
CCC   Gaussian in W^2 with a cut-off at the single nucleon removal energy.      CCC
CCC                                                                             CCC      
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! fit to low q2 dip region: purefly empirical
! assume contribution is purely transverse
      implicit none
      real*8 Z,A,N,q2,w2,mp/0.938272/,mp2,mn/0.93957/,w,qv2,nu,numin
      real*8 Y,Y2,a1,a2,b1,b2,c1,c2,t1,t2,dw2,xmax,q20,w2min,sigqd,psid
      real*8 egam,ex,nuel,cof,x,f1qd, xvalm(100),alpha,pi2/9.86959/
      real*8 gep,gen,gmp,gmn
      
      N = A-Z
      alpha = 1/137.036
      xmax = 50.0
      q20 = 0.01
c      q20 = xvalm(26)
      mp2 = mp*mp

      f1qd = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      egam = nu*1000.0   !!! put in MeV
      x  = q2/(2.0*mp*nu)
      qv2 = q2+nu**2.0
      if(A.EQ.12) numin = 0.0165
      if(A.EQ.27) numin = 0.0085
      if(A.EQ.40) numin = 0.0085
      w2min = mp2+2.0*mp*numin-q2
      xmax = q2/2.0/mp/numin

      nuel = q2/2./(0.931494*A)
      ex = nu-nuel

      
      if(A.lt.2.5) return


      psid = 0.0
      if(egam.LT.20.0) then
        psid = exp(-73.3/egam)
      elseif(egam.GE.20.0.AND.egam.LT.140.0) then
        psid = 8.3714E-2-9.8343E-3*egam+4.1222E-4*egam*egam
     &        -3.4762E-6*egam**3+9.3537E-9*egam**4
      elseif(egam.GE.140.0) then
         psid = exp(-24.2/egam)
c         psid = exp(-24.2/egam)/(1.0+0.00003*egam*egam**0.5)/
c     &          exp(-24.2/140.0)
      endif

      if(egam.LE.20.0) psid = 0.0
      sigqd = 397.8*n*z*(egam-2.224)**1.5/egam**3*psid !!! in microbarns
c      sigqd = sigqd/(1.0+q2/xvalm(26))**xvalm(27)
      sigqd = sigqd/(1.0+Q2/0.1)**5.0

      sigqd = 0.0
      
      f1qd = sigqd/8.0/pi2/alpha/3.894e3*abs(w2-mp2)

      f1qd = f1qd*1000.0

      if(q2.GT.0.5) f1qd = 0.0
      
c      dw2 = w2-w2min

c      if(dw2.LT.0.0) dw2 = 0.0
      
      
c      if(nu.LT.numin) sigd = 0.0
c      cof = (nu-numin)**0.5/(0.020-numin)**0.5
c      cof = min(1.0,cof)
c      cof = max(cof,0.0)
c      if(ex.LE.numin) cof = 0.0
       
c       F1mec = F1mec*cof

      
c      if(dw2.LE.0.0.OR.x.GT.xmax) f1mec = 0.0
c      if(f1mec.LE.1.0E-9.OR.x.GE.xmax) f1mec=0.0


c      if(egam.LT.140.0.AND.q2.LT.0.02) write(6,*) egam,sigqd


c      sigqd = sigqd*1000.0
c       write(6,*) egam,sigqd
      
      return
      end


CCC-----------------    

      
