C=======================================================================

      Subroutine FORMFACTS(q2,gmp,gep,gmn,gen)

CCC   Returns proton and neutron form factors based on new proton fit         CCC
CCC   by M.E. Christy including GMp12 data.  Neutron starts with              CCC
CCC   Kelly fit, but includes correction factors extracted from fit to        CCC      
CCC   inclusive deuteron data.   Version from July 22, 2020                   CCC                      
!--------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 q2,tau,gd,gmp,gep,gmn,gen,mcor,ecor
      REAL*8 mu_p/ 2.792782/ ,mu_n/ -1.913148 /
      REAL*8 mp/ 0.9382727 /, mp2
      
      mp2 = mp*mp
      tau   = q2 / 4.0 / mp2 

CCC   2021 Christy fit to GMp and GEp including GMp12 data  CCC
        
       GMP = mu_p*(1.+0.099481*tau)/
     &  (1.0+11.089*tau+19.374*tau*tau+5.7798*tau**3)
       GEP = (1.+0.24482*tau**2)/
     &  (1.0+11.715*tau+11.964*tau*tau+27.407*tau**3)       
       GD = (1./(1 + q2/0.71))**2

CCC   2021 fit to deuteron inclusive data   CCC       
      
c        GMn = mu_n*(1.0+0.12417E-04*tau)/
c     &   (1.000+11.404*tau+17.014*tau*tau+31.219*tau**3)

c         GEn = (1.5972*tau / (1.0+0.19655*tau)) * GD

        GMn = mu_n           !!!  Kelly Fit
     &       * (1.D0 + 2.330D0*tau)
     &      / (1.D0 + 14.720D0*tau + 24.200D0*tau**2 + 84.100D0*tau**3)
        
        GEn = (1.700D0*tau/(1+ 3.300D0*tau))*GD

        GEn = GEn*((q2+1189.4)/1189.4)**219.73 
        GMn = GMn/((q2+0.35590 )/0.35590 )**0.93020E-01 
        
             
       return
       end
