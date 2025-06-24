      SUBROUTINE RESMODN(sf,w2,q2,xval,sig)
CCC  Returns proton transverse (sf=1) and longitudinal (sf=2 Cross sections CCC      
CCC  Version from July 15, 2021  -  Author:  M.E. Christy                    CCC
CCC  This routine returns proton photo-absorbtion cross sections            CCC
CCC  for either transverse or longitudinal photons in units of ub/Sr/Gev.   CCC
CCC                                                                         CCC
CCC  Fit form is empirical.  Interpret physics from it at your own risk.    CCC
CCC  replaced 2-pi threshold with eta                                       CCC
      
      IMPLICIT NONE
      REAL*8 W,w2,q2,mp,mp2,mpi2,xb,sig,xval(50),mass(7),width(7)
      REAL*8 height(7),rescoef(6,4)
      REAL*8 nr_coef(3,4),sigr(7),wdif(2),sig_nr,pi2,alpha
      REAL*8 mpi,meta,intwidth(7),k,kcm,kr(7),kcmr(7),ppicm,ppi2cm
      REAL*8 petacm,ppicmr(7),ppi2cmr(7),petacmr(7),epicmr(7),epi2cmr(7)
      REAL*8 eetacmr(7),epicm,epi2cm,eetacm,br(7,3),ang(7)
      REAL*8 pgam(7),pwid(7,3),x0(7),dip,mon,q20,A0
      REAL*8 sig_res,xpr(2),t1,t2
      INTEGER i,j,num,sf
      
      mp = 0.939565
      mpi = 0.134977
      mpi2 = mpi*mpi
      meta = 0.547862
      mp2 = mp*mp
      alpha = 1./137.036
      W = sqrt(w2)
      wdif(1) = w - (mp + mpi)
      wdif(2) = w - (mp + meta)
      
      q20 = 0.05
      q20= xval(50)

       
CCCC   single pion branching ratios  CCCC

      br(1,1) = 1.00       !!!  P33(1232)       
      br(2,1) = 0.45      !!!  S11(1535)   
      br(3,1) = 0.60      !!!  D13(1520)
      br(4,1) = 0.65      !!!  F15(1680)
      br(5,1) = 0.60       !!!  S11(1650)
      br(6,1) = 0.65     !!!  P11(1440) roper 
      br(7,1) = 0.60      !!!  F37(1950)

CCCC  eta branching ratios   CCCC

      br(1,3) = 0.0       !!!  P33(1232)
      br(2,3) = 0.40      !!!  S11(1535) 
      br(3,3) = 0.08      !!!  D13(1520)
      br(4,3) = 0.0       !!!  F15(1680)
      br(5,3) = 0.20      !!!  S11(1650)
      br(6,3) = 0.0       !!!  P11(1440) roper   
      br(7,3) = 0.0       !!!  F37(1950)

CCCC  2-pion branching ratios  CCCC

      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo

CCCC   Meson angular momentum   CCCC

      ang(1) = 1.       !!!  P33(1232)
      ang(2) = 0.       !!!  S11(1535)
      ang(3) = 2.       !!!  D13(1520)
      ang(4) = 3.       !!!  F15(1680)
      ang(5) = 0.       !!!  S15(1650)
      ang(6) = 1.       !!!  P11(1440) roper   
      ang(7) = 3.       !!!  F37(1950)

      do i=1,7     !!!  resonance damping parameter  !!!
         x0(i) = 0.16   !!! 
      enddo

      if(sf.EQ.2) x0(1) = 0.07   !!!  different Delta mass for sigL
  
      
      do i=1,7
        br(i,2) = 1.-br(i,1)-br(i,3)
      enddo
               
      dip = 1./(1.+q2/1.15)**2.             !!!  Dipole parameterization  !!!
      mon = 1./(1.+q2/1.5)**1.

      
      xb = q2/(q2+w2-mp2)
      xpr(1) = 1.00+(w2-(mp+mpi)**2)/(q2+q20)
      xpr(1) = 1./xpr(1)
      xpr(2) = 1.+(w2-(mp+meta)**2)/(q2+q20)
      xpr(2) = 1./xpr(2)
      if(w.LE.(mp+mpi)) xpr(1) = 1.0      
      if(w.LE.(mp+meta)) xpr(2) = 1.0
     

CCC    Calculate kinematics needed for threshold Relativistic B-W  CCC

      k = (w2 - mp2)/2./mp
      kcm = (w2-mp2)/2./w

      epicm = (W2 + mpi**2 -mp2 )/2./w
      ppicm = SQRT(MAX(0.0,(epicm**2 - mpi**2)))
      epi2cm = (W2 + (2.*mpi)**2 -mp2 )/2./w
      ppi2cm = SQRT(MAX(0.0,(epi2cm**2 - (2.*mpi)**2)))
      eetacm = (W2 + meta*meta -mp2 )/2./w
      petacm =  SQRT(MAX(0.0,(eetacm**2 - meta**2)))

      num = 0

      do i=1,6              !!!  Read in resonance masses     !!!
        num = num + 1
        mass(i) = xval(i)
      enddo
      do i=1,6              !!!  Read in resonance widths     !!!
        num = num + 1
        intwidth(i) = xval(num)
        width(i) = intwidth(i)
      enddo

      mass(7) = xval(47)
      intwidth(7) = xval(48)
      width(7) = intwidth(7) 

      do i=1,7
        kr(i) = (mass(i)**2-mp2)/2./mp
        kcmr(i) = (mass(i)**2.-mp2)/2./mass(i)
        epicmr(i) = (mass(i)**2 + mpi**2 -mp2 )/2./mass(i)
        ppicmr(i) = SQRT(MAX(0.0,(epicmr(i)**2 - mpi**2)))
        epi2cmr(i) = (mass(i)**2 + (2.*mpi)**2 -mp2 )/2./mass(i)
        ppi2cmr(i) = SQRT(MAX(0.0,(epi2cmr(i)**2 - (2.*mpi)**2)))
        eetacmr(i) = (mass(i)**2 + meta*meta -mp2 )/2./mass(i)
        petacmr(i) =  SQRT(MAX(0.0,(eetacmr(i)**2 - meta**2)))

CCC   Calculate partial widths   CCC

        pwid(i,1) = intwidth(i)*(ppicm/ppicmr(i))**(2.*ang(i)+1.)
     &           *((ppicmr(i)**2+x0(i)**2)/(ppicm**2+x0(i)**2))**ang(i)
c         !!!  1-pion decay mode


        pwid(i,2) = intwidth(i)*(ppi2cm/ppi2cmr(i))**(2.*ang(i)+4.)
     &         *((ppi2cmr(i)**2+x0(i)**2)/(ppi2cm**2+x0(i)**2))
     &         **(ang(i)+2)   !!!  2-pion decay mode

        pwid(i,2) = W/mass(i)*pwid(i,2)


        pwid(i,3) = 0.          !!!  eta decay mode


        if(i.EQ.2.OR.i.EQ.5) then
          pwid(i,3) =  intwidth(i)*(petacm/petacmr(i))**(2.*ang(i)+1.)
     &          *((petacmr(i)**2+x0(i)**2)/(petacm**2+x0(i)**2))**ang(i)
c         !!!  eta decay only for S11's 
        endif 


        pgam(i) = (kcm/kcmr(i))**2*
     &                   (kcmr(i)**2+x0(i)**2)/(kcm**2+x0(i)**2)

        pgam(i) = intwidth(i)*pgam(i)

        width(i) = br(i,1)*pwid(i,1)+br(i,2)*pwid(i,2)+br(i,3)*pwid(i,3)

      enddo

CCC    End resonance kinematics and Widths calculations   CCC


CCC    Begin resonance Q^2 dependence calculations   CCC
           

      do i=1,6
        do j=1,4
          num = num + 1
          rescoef(i,j)=xval(num)
        enddo

        if(sf.EQ.1) then

          height(i) = rescoef(i,1)*
     &          (1.+rescoef(i,2)*q2/(1.+rescoef(i,3)*q2))*
     &          mon**rescoef(i,4)


        else

           height(i) = (rescoef(i,1)+rescoef(i,2)*q2)
     &           *exp(-1.*rescoef(i,3)*q2)

          
        endif
 
        height(i) = height(i)*height(i)

      enddo


      if(sf.EQ.2) then      !!!  4th resonance region  !!!

        height(7) = (xval(44)+xval(45)*q2)*exp(-1.0*xval(46)*q2)
        
      else
        height(7) = xval(49)*mon
      endif
      
      height(7) = height(7)*height(7)


CCC    End resonance Q^2 dependence calculations   CCC

     
      do i=1,3               !!!  Non-Res coefficients  !!!
        do j=1,4
          num = num + 1
          nr_coef(i,j)=xval(num)
        enddo
      enddo


CCC   Calculate Breit-Wigners for all resonances   CCC

      sig_res = 0.0

      do i=1,7

       sigr(i) = width(i)*pgam(i)/((W2 - mass(i)**2.)**2. 
     &              + (mass(i)*width(i))**2.)
        sigr(i) = height(i)*kr(i)/k*kcmr(i)/kcm*sigr(i)/intwidth(i)
        sig_res = sig_res + sigr(i)   
      enddo

      sig_res = sig_res*w
      if(sf.EQ.2) sig_res = sig_res*q2


CCC    Finish resonances / start non-res background calculation   CCC

 
      sig_nr = 0.

      if(sf.EQ.1.and.xpr(1).LT.1.0) then
c     A0 = xval(37)*(1.0+xval(44)*q2)/(1.+q2/xval(42))**xval(43)  !!! overall amplitude

        A0 = xval(37)/(1.0+q2/xval(42))**xval(43)
         
        t1 = xval(38)*log(1.05+q2)+xval(39)/(1.05+q2)                 !!! exponent of (1-xpr)
c     t2 = xval(40)*log(2.0+q2/xval(41))+xval(45)        !!! exponent of xpr
        t2 = xval(40)*(1.0+q2/xval(41))**xval(44)
 


        if(xpr(1).LE.1.0) then                             !!! 1-pi threshold
          sig_nr = 389.4*A0*(1.-xpr(1))**t1*xpr(1)**t2
        endif
          
        if(xpr(2).LE.1.0) then                             !!! eta threshold 
          sig_nr = sig_nr+xval(46)*389.4*A0*(1.-xpr(2))**t1*xpr(2)**t2
         endif

c        write(6,*) q2,sf,t1,t2
         
      elseif(sf.EQ.2.and.xpr(1).LT.1.0) then
         A0 = xval(37)/(1.0+q2/xval(39))**2.0
         t1 = xval(38)/(1.0+q2/(xval(40)))+xval(32)*log(q2+xval(36))                    !!! exponent of (1-xpr)
         t2 = xval(41)/(1.00+q2/xval(42))**xval(43)       !!! exponent of xpr

c         write(6,*) q2,a0,t1,t2
        
        if(xpr(1).LE.1.0) then
          sig_nr = sig_nr + 389.4*A0*
     &       xb*(1.-xpr(1))**t1*xpr(1)**t2

        endif

      endif
    
      sig = sig_res + sig_nr

      if((w-mp).LT.wdif(1)) sig = 0.0    


 1000  format(8f12.5)
 1001  format(7f12.3)

      RETURN 
      END 



