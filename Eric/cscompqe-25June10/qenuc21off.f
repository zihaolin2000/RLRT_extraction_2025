      SUBROUTINE QENUC21OFF(Z, A, Q2, W2, xvalc,F1, F2)

C Calculates quasielastic A(e,e')X structure functions F1 and F2 PER NUCLEUS  
c for A>2. Uses superscaling from Sick, Donnelly, Maieron, nucl-th/0109032
c based on the earlier code F1F2QE09 by P. Bosted.  The superscaling distribution
c shape is determined from the fit to 12C data      
c
c input: Z, A  (real*8) Z and A of nucleus (should be 2.0D0 for deuteron)
c        Q2 (real*8) is 4-vector momentum transfer squared (positive in
c                     chosen metric)
c        W2 (real*8) is invariant mass squared of final state calculated
c                     assuming electron scattered from a free proton
c                 
c outputs: F1, F2 (real*8) are structure functions per nucleus


      IMPLICIT NONE     

      REAL*8 Z, A, avgN, F1, F2, W2, Q2, R
      REAL*8 mp/0.938272/
      REAL*8 PAULI_SUP1, PAULI_SUP2,x2,pb2,pb2L
      REAL*8 GEP, GEN, GMP, GMN
      REAL*8 RMUP/ 2.792782/ ,RMUN/ -1.913148 /                
      REAL*8 nu, nuel, ex, qv, TAU, FY, FYL, FYp, FYn, pauli, nup,num
      REAL*8 kappa, lam, lamp, lampn,taup, xi, ximax, psi, psip, psipn
      REAL*8 psimax, nuL,nut,kf, es, esmin, GM2bar, GE2bar, Delta, GL,GT
      REAL*8 F1ff,F2ff,GMoff,GEoff,moff,xvalc(45),B,xm,m,BL,cof
      integer IA,paulitype

      psimax = 5.0
      moff =  1.0*mp
      paulitype = 1    !!! 1 = Superscaling, 2 = Tsai from Fermi Gas  !!!
      
! return if proton: future change this to allow for
! equivalent W resolution
      F1 = 0.
      F2 = 0.
      IA = int(A)
      avgN = A - Z 
      IF (iA.EQ.1) RETURN
      nuel = q2/2./(0.931494*A)

      
      
! some kinematic factors. Return if nu or q2 is negative
      Nu = (W2 - mp**2 + Q2) / 2. / mp
      if(nu .le. 0.0 .or. Q2 .le. 0.) return
      TAU   = Q2 / 4.0 / mp**2                                        
      qv = sqrt(nu**2 + Q2)

      ex = nu-nuel
      
!     Call New FormFactors  !
      
      call formfacts(Q2,gmp,gep,gmn,gen)

!  Get energy shift and fermi momementum from fit  !      
      
      if(IA.GT.2) then
         Esmin = 0.020           !!! Minimum Es given by removal energy
         if(IA.EQ.12) Esmin = 0.0165
         kf = xvalc(35)
         Es = Esmin + xvalc(36)  !!! Maximum Es
         if(qv.LT.1.5) Es = Es-xvalc(36)*(1.0-qv/1.5)**0.1  !!!  test  !!!        
      endif 
      nup = nu-Es

! Pauli suppression model from Tsai RMP 46,816(74) eq.B54
      IF((QV .GT. 2.* kf).OR.(iA.EQ.1)) THEN
        PAULI_SUP2 =1.0
      ELSE
        PAULI_SUP2 = 0.75 * (QV / kf) * (1.0 - ((QV / kf)**2)/12.)
      ENDIF
      PAULI_SUP1 = PAULI_SUP2

!     structure functions with off shell factors
      
      kappa = qv / 2. / mp
      lam = nu / 2. / mp
      
      lamp = nup / 2. / mp
      lampn = -lamp
      taup = kappa**2 - lamp**2
      xi = sqrt(1. + (kf/moff)**2) -1.
! Very close to treshold, could have a problem
c      if(1.+lamp.le.0.) return
c      if(taup * (1. + taup).le.0.) return

      psi =  (lam  - taup ) / sqrt(xi) /
     >     sqrt((1.+lam )* taup + kappa * sqrt(taup * (1. + taup)))   !!! OK

      psip = (lamp - taup) / sqrt(xi) / 
     >     sqrt((1.+lamp)*taup + kappa * sqrt(taup * (1. + taup)))    !!! OK

      psipn = (lampn - taup) / sqrt(xi) / 
     >     sqrt((1.+lampn)*taup + kappa * sqrt(taup * (1. + taup)))   !!! OK

      
      
      nuL = (Q2/qv/qv)**2
      
c changed definition of nuT from
c      nuT = tau / 2. / kappa**2 + tan(thr/2.)**2
c to this, in order to separate out F1 and F2 (F1 prop. to tan2 term)
      nuT = tau / 2. / kappa**2 


      GM2bar = (Z * GMP**2 + avgN * GMN**2)  
      GE2bar = (Z * GEP**2 + avgN * GEN**2)

      F1ff =  (tau*sqrt(GM2bar)+sqrt(GE2bar))/(1.0+tau) 
      F2ff = (sqrt(GM2bar)-sqrt(GE2bar))/(1.0+tau)
      GEoff = F1ff-tau*moff/mp*F2ff         !!!  Off-shell FF   !!!
      GMoff = F1ff+moff/mp*F2ff             !!!  Off-shell FF   !!!
     
c      write(6,*) "off :",q2,sqrt(GM2bar),GMoff,sqrt(GE2bar),GEoff
      
      Delta = tau/kappa/kappa*xi*(1.-psi**2)*
     & (kappa*sqrt(1.0+1.0/tau)+xi/3.0*(1.-psi**2))

c      write(6,*) w2,q2,Delta
      
      GL = kappa**2 / tau*
     >  (GEoff**2. +(GEoff**2.+tau*GMoff**2.)*Delta/(1.0+tau))
 
      GT = (2.*tau*GMoff**2.+(GEoff**2.+tau*GMoff**2.)*Delta/(1.+tau))
      
      num = 2. *kappa *(1. + xi * (1. + psi**2) / 2.)

      GL = GL/num
      GT = GT/num

c       GL = 2.0*xi*kf*(Moff/kf)**3.0/qv*GL
c       GT = 2.0*xi*kf*(Moff/kf)**3.0/qv*GT
  

! from Maria Barbaro: see Amaro et al., PRC71,015501(2005).
      FY = 1.5576 / (1. + 1.7720**2 * (psip + 0.3014)**2) / 
     >   (1. + exp(-2.4291 * psip))

CCMEC - Fitted Superscaling distribution      CCCC

      FYp = xvalc(7)/ (1. + xvalc(8)**2 * (psip + xvalc(9))**2) / 
     >  (1. + exp(xvalc(10) * psip))*(1.0-abs(psip)/psimax)**2.0*
     >  (1.0+abs(psip)/psimax)**2.0


      FYn = xvalc(7)/ (1. + xvalc(8)**2 * (psipn + xvalc(9))**2) / 
     > (1. + exp(xvalc(10) * psipn))*(1.0-abs(psipn)/psimax)**2.0*
     > (1.0+abs(psipn)/psimax)**2.0


      if(psip.GT.psimax) Fyp = 0.0
      if(psipn.GT.psimax) Fyn = 0.0

      FYp = max(0.0,FYp)     
      FYn = max(0.0,FYn)      
      
      if(paulitype.EQ.1) then        !!! Superscaling
        FY = max(0.0,FYp-FYn)

c        if(qv.LT.0.1) write(6,*) qv,w2,FYp,FYn,FY
        
        x2 = qv/kf
        
        pb2L = 1.0

c        xm = xvalc(40)  !!! 
c        m = 2.0
c        BL = xm**m/(1.0-1.0/xm/xvalc(41))
        
c        if(x2.LT.xm) then
c           pb2L = xvalc(41)*(x2-0.2)*(1.0-(x2**m)/BL)

c           pb2L = pb2L/(1.0-x2/xm)**xvalc(44)
           
c          pb2L = min(1.0,pb2L)
c        endif            




        
CCC   Below is the nominal without the extra (4.0-x2)**1.5 term  CCC
        
        pb2L = 1.0-xvalc(40)*(4.0-x2)**2.5-xvalc(41)*(4.0-x2)**3.5
        
c        pb2L = pb2L-xvalc(6)*exp(-1.0*(x2-2.55)**2.0/xvalc(44)**2)
        pb2L = pb2L-xvalc(44)*(4.0-x2)**1.5
        pb2L = pb2L*(x2-0.2)**2/(x2-0.18)**2.0

CCC        

        if(x2.GT.4.0) pb2L = 1.0
        pb2L = min(pb2L,1.0)
        pb2L = max(pb2L,0.0)      


ccc Next is Mahalia  ccc

c         pb2L = 1.0-0.00051412*(4.0-x2)**2.5-0.0020088*(4.0-x2)**4.5        
c         pb2L = pb2L-0.29654*exp(-1.0*(x2-0.91)**2.0/0.65**2)
c         pb2L  = pb2L*(x2-0.1)**2/(x2-0.07)**2

ccc        
        
        if(x2.GT.4.0) pb2L = 1.0
        pb2L = min(pb2L,1.0)

  
        
        if(psip.GT.psimax) FY = 0.0  !!! effective cutoff 


c        if(x2.LT.4.0) write(6,*) x2,pb2L

        
        FYL = FY
        FYL = pb2L*FYL

      elseif(paulitype.EQ.2) then    !!! Tsai - Fermi Gas
         FY = Pauli_sup2*FYp
      endif
     
      
       F2 = nu/kf * (FYL*nuL*GL + FY*nuT*GT)
       F1 = mp * FY/kf * GT / 2.

       cof = (ex-esmin)**0.5/(0.025-esmin)**0.5
       cof = min(1.0,cof)
       cof = max(cof,0.0)
       if(ex.LE.esmin) cof = 0.0
       
       F1 = F1*cof
       F2 = F2*cof
       
       if(F2.LT.0.0) F2 = 0.0
       if(F1.LT.0.0) F1 = 0.0

c       if(ex.LT.30) write(6,*) q2,nu,ex,cof
       
      return
      end


     
CCC-----------------

      
