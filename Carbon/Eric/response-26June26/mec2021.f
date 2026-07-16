      SUBROUTINE MEC2021(z,a,w2,q2,xvalm,f1mec)

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
      real*8 Z,A,q2,w2,mp/0.938272/,mp2,mn/0.93957/,w,qv2,nu,numin
      real*8 Y,Y2,a1,a2,b1,b2,c1,c2,t1,t2,dw2,xmax,q20,w2min,w2min2
      real*8 dw22,ex,nuel,cof,x,f1mec, f1mec2, xvalm(100),q2t


      q2t = q2
c      if(q2.LE.0.01) q2t = 0.01
      xmax = 50.0
      q20 = 0.01
      q20 = xvalm(26)
      mp2 = mp*mp

      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu)
      qv2 = q2+nu**2.0
      if(A.EQ.12) numin = 0.0165
      if(A.EQ.27) numin = 0.0085
      if(A.EQ.40) numin = 0.0085
      w2min = mp2+2.0*mp*numin-q2
      w2min2 = 1.059
      xmax = q2/2.0/mp/numin              !!!  adjusted to improve fit  !!!

      nuel = q2/2./(0.931494*A)
      ex = nu-nuel

      
      if(A.lt.2.5) return

      a1 = xvalm(1)
      a2 = xvalm(35)

c      a2 = 0.0

      
       Y = (q2+q20)**xvalm(2)
     &     /(xvalm(3)+q2)**xvalm(4) 

c       Y = (q20+q2**xvalm(2))
c     &     /(1.0+q2/xvalm(3))**xvalm(4)
c       Y2 = Y
c       Y2 = Y*(q2+q20)**xvalm(19)/(q2+q20)**2
      Y2 =  Y*(q2+xvalm(18))**xvalm(19)
     &     /(xvalm(20)+q2)**2
c      Y2 = Y2/(1.0+0.5*sqrt(q2))

c      Y2 = 0.0
      
c      Y2 = xvalm(35)*Y*(xvalm(20)+q2)**0.5
      
      b1 = xvalm(5)
c      b1 = xvalm(5)-0.006*sqrt(A)
c      b2 = 1.275
      b2 = 1.3506
      b2 = xvalm(6)
      
      c1 = xvalm(32)+xvalm(33)*q2*sqrt(q2)

      c2 = 0.160
      c2 = xvalm(17)
c     &   +xvalm(33)*q2*sqrt(q2)
      
      t1 = (w2-b1)**2/(2.*c1**2)
      t2 = (w2-b2)**2/(2.*c2**2)
      
      dw2 = w2-w2min
      dw22 = w2-w2min2
      
      if(dw2.LT.0.0) dw2 = 0.0
      
      f1mec = Y*exp(-1.0*t1)
      f1mec = f1mec*(dw2)**1.00

      if(dw22.LT.0.0) dw22 = 0.0      
c      f1mec2 = Y2*exp(-1.*t2)

      f1mec2 = Y2/((w2-b2)**2+b2**2*c2**2)
      
c      f1mec2 = f1mec2*(dw22)**2.5
      f1mec2 = f1mec2*(dw22)**1.0 
      
c      f1mec2 = f1mec2*(q2+xvalm(18))/(q2+q20)**0.5
      
      f1mec = A*(a1*f1mec+a2*f1mec2)
      
      if(nu.LT.numin) f1mec = 0.0
      cof = (1.0-(numin/nu)**4)
      cof = min(1.0,cof)
      cof = max(cof,0.0)
      if(ex.LE.numin) cof = 0.0
       
       F1mec = F1mec*cof

      
      if(dw2.LT.0.0.OR.dw22.LT.0.0) f1mec = 0.0
      if(f1mec.LE.1.0E-9) f1mec=0.0

      

      return
      end


CCC-----------------    

      
