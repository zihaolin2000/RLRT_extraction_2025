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
      real*8 Y,a1,a2,b1,b2,c1,c2,t1,t2,dw2,xmax,q20,w2min
      real*8 ex,nuel,cof,x,f1mec, f1mec2, xvalm(45)

      xmax = 50.0
      q20 = 0.00001
c      q20 = xvalm(6)
      mp2 = mp*mp

      f1mec = 0.0
      if(w2.le.0.0) return
      w  = sqrt(w2)
      nu = (w2 - mp2 + q2)/2./mp
      x  = q2/(2.0*mp*nu)
      qv2 = q2+nu**2.0
      numin = 0.0165
      w2min = mp2+2.0*mp*numin-q2
      xmax = q2/2.0/mp/numin

      nuel = q2/2./(0.931494*A)
      ex = nu-nuel

      
      if(A.lt.2.5) return

      a1 = xvalm(1)
      a2 = xvalm(45)

c      a2 = 0.01045

      
      Y = A*exp(-1.0*q2*q2/xvalm(2))*(q2+q20)**2
     &     /(xvalm(3)+q2)**xvalm(4) 
      
      b1 = xvalm(5)
      b2 = 1.275
  
      c1 = xvalm(42)+xvalm(43)*q2
      c2 = 0.265
      
      t1 = (w2-b1)**2/(2.*c1**2)
      t2 = (w2-b2)**2/(2.*c2**2)
      
      dw2 = w2-w2min

      if(dw2.LT.0.0) dw2 = 0.0
      f1mec = Y*exp(-1.0*t1)
      f1mec = f1mec*(dw2)**1.5  !!!  1.254



      f1mec2 = a2*Y*(q2+q20)**1.5*exp(-1.*t2)

      
      f1mec = a1*f1mec+f1mec2
      
      if(nu.LT.numin) f1mec = 0.0
      cof = (ex-numin)**0.5/(0.025-numin)**0.5
      cof = min(1.0,cof)
      cof = max(cof,0.0)
      if(ex.LE.numin) cof = 0.0
       
       F1mec = F1mec*cof

      
      if(dw2.LE.0.0.OR.x.GT.xmax) f1mec = 0.0
      if(f1mec.LE.1.0E-9.OR.x.GE.xmax) f1mec=0.0

      

      

      return
      end


CCC-----------------    

      
