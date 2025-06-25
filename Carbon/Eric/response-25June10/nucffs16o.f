      SUBROUTINE NUCFFS16O(A,Z,Q2,state,FF)
CCCC  Returns the 12C nuclear scattering form factors for the following states   CCCC
CCCC  State = 1 - nuclear elastic, 2 - 4.4 MeV, 3 - 7.6 MeV, 4 - 9.6 MeV         CCCC
CCCC  5 - 10.84 MeV, 6 -15.1 MeV state, 7 - 16.1 MeV, 8 - 18.6 MeV, 9 - 20 MeV   CCCC
CCCC  10 - GDR, 21-26 MeV (mean at 23?)                                        CCCC
CCCC  11- GDR, 26-37 MeV (mean at 31.5?), 12 - ~15.1 MeV state                    CCCC
      
      IMPLICIT NONE

      Real*8 A,Z,Q2p,Q2,q2f,qf,FF,FF2,x2,alp,char,g1,g2,g3,e1,Radius
      Real*8 a0,a1,alph,g0,h0
      integer state,IA,IZ

      IA = int(A)
      IZ = int(Z)
      q2f = Q2/0.1975/0.1975    !!! Q^2 in inverse Fermi  !!!
      qf = sqrt(q2f)
      Radius = 2.5 
      Radius = 2.45
      
      FF = 0.0
      if(state.EQ.1) then
        x2 = (Q2/0.197328**2)*Radius**2
        alp = (Z-2.)/3.                                               
        char = x2*(2.+3.*alp)/(12.+30.*alp)                             
        if (char.lt.80) FF = exp(-char)*(1.-alp*x2/(6.+15.*alp))
        FF2 = FF*FF
c        write(6,*) "1: ", q2f,FF2

        
CCCC  Arie Bodek's fit using the harmomic well shape from Hofstadter  CCC

        alph = 4.0/3.0
        a0 = 1.65
        g0 = 8.0E-5*exp(-1.0*((q2f-2.9)/0.44)**2.0)
        h0 = (1.0-alph*q2f*a0*a0/2.0/(2.0+3.0*alph))*
     &         exp(-q2f*a0*a0/4.0)
        
c
        FF2 = h0*h0+g0

        
CCCC  The excited states are from a fit by Arie Bodek, 2022  CCCC
        
      elseif(state.EQ.2) then   !!!  6.0494 MeV

        g1 = 0.8E-3*exp(-1.0*(q2f-0.3)**2/0.5/0.5)
        g2 = 0.2E-3*exp(-1.0*(q2f-1.0)**2/0.55/0.55)  
        g3 = 0.0055E-3*exp(-1.0*(q2f-4.3)**2/1.5/1.5)
        e1 = -0.22E-3*exp(-6.0*q2f)     
        
      elseif(state.EQ.3) then     !!!  6.1299 MeV

        g1 = 1.5E-3*exp(-1.0*(q2f-0.90)**2/0.45/0.45)
        g2 = 4.1E-3*exp(-1.0*(q2f-1.1)**2/1.7/1.7)
        g3 = 0.2E-3*exp(-1.0*(q2f-2.55)**2/2.4/2.4)
        e1 = -3.1E-3*exp(-1.1*q2f)
 
      elseif(state.EQ.4) then     !!!  6.9171 MeV

        g1 = 4.4E-3*exp(-1.0*(q2f-0.4)**2/0.6/0.6)
        g2 = 1.0E-3*exp(-1.0*(q2f-1.2)**2/0.88/0.88)
        g3 = 0.003E-3*exp(-1.0*(q2f-5.0)**2/1.5/1.5)
        e1 = -0.5E-3*exp(-1.66*q2f)
        
      elseif(state.EQ.5) then   !!!  7.1169 MeV

        g1 = 1.0E-3*exp(-1.0*(q2f-0.85)**2/0.45/0.45)
        g2 = 0.9E-3*exp(-1.0*(q2f-1.4)**2/1.1/1.1)
        g3 = 0.17E-3*exp(-1.0*(q2f-2.1)**2/1.7/1.7)
        e1 = -0.3E-3*exp(-4.0*q2f)
 
      elseif(state.EQ.6) then   !!!  9.8445 MeV   

        g1 = 0.09E-3*exp(-1.0*(q2f-0.7)**2/0.65/0.65)
        g2 = 0.08E-3*exp(-1.0*(q2f-1.7)**2/1.2/1.2)
        g3 = 0.012E-3*exp(-1.0*(q2f-2.5)**2/2.0/2.0)
        e1 = 0.0007E-3*exp(-2.0*q2f)
  
      elseif(state.EQ.7) then  !!! 10.356 MeV 

        g1 = 0.09E-3*exp(-1.0*(q2f-1.2)**2/0.65/0.65)
        g2 = 0.095E-3*exp(-1.0*(q2f-2.0)**2/1.2/1.2)
        g3 = 0.011E-3*exp(-1.0*(q2f-3.4)**2/1.5/1.5)
        e1 = 0.0007E-3*exp(-2.0*q2f)
        
      elseif(state.EQ.8) then   !!! 11.0967 MeV  - here

        g1 = 0.05E-3*exp(-1.0*(q2f-1.1)**2/0.7/0.7)
        g2 = 0.04E-3*exp(-1.0*(q2f-2.2)**2/1.1/1.1)
        g3 = 0.01E-3*exp(-1.0*(q2f-3.2)**2/1.8/1.8)
        e1 = 0.01E-3*exp(-1.66*q2f)
         
      elseif(state.EQ.9) then   !!! 11.520 MeV

        g1 = 1.9E-3*exp(-1.0*(q2f-0.45)**2/0.6/0.6)
        g2 = 0.6E-3*exp(-1.0*(q2f-1.0)**2/1.0/1.0)
        g3 = 0.003E-3*exp(-1.0*(q2f-5.5)**2/1.8/1.8)
        e1 = 0.08E-3*exp(-1.66*q2f)
 
      elseif(state.EQ.10) then   !!! 12.0490 MeV
        g1 = 1.0E-3*exp(-1.0*(q2f-0.18)**2/0.8/0.8)
        g2 = 0.25E-3*exp(-1.0*(q2f-1.0)**2/0.89/0.89)
        g3 = 0.0015E-3*exp(-1.0*(q2f-5.0)**2/1.5/1.5)
        e1 = 0.2E-3*exp(-4.0*q2f)
 
      elseif(state.EQ.11) then   !!! GDR 12.5-20 MeV - here
        g1 = 2.0E-3*exp(-1.0*(q2f-0.35)**2/0.3/0.3)
        g2 = 8.0E-3*exp(-1.0*(q2f-0.0)**2/1.9/1.9)
        g3 = 0.0
        e1 = 0.007E-3*exp(-1.66*q2f)
        
      elseif(state.EQ.12) then   !!! GDR 20-35 MeV (31.5 MeV ?)
        g1 = 18.0E-3*exp(-1.0*(q2f-0.0)**2/0.4/0.4)
        g2 = 7.0E-3*exp(-1.0*(q2f-0.6)**2/0.8/0.8)
        g3 = 2.5E-3*exp(-1.0*(q2f-1.0)**2/1.5/1.5)
        e1 = 0.0
        
      elseif(state.EQ.13) then   !!! add 35-50 MeV state
        g1 = 0.8E-3*exp(-1.0*(q2f-0.0)**2/0.3/0.3)         
        g2 = 1.3E-3*exp(-1.0*(q2f-0.6)**2/3.0/3.0)        
        g3 = 0.0
        e1 = 0.004E-3*exp(-1.66*q2f)
        
      endif
      
      FF2 = q2f*(g1+g2+g3+e1)
      FF2 = max(0.0,FF2)

      if(qf*qf.GT.12.0) FF2 = 0.0
      FF = sqrt(FF2)

c      if(state.EQ.10) write(6,*) state,qf,FF2
c      if(state.EQ.4) write(6,*) state, q2,q2f,FF2,g1,g2,g3,e1

      return
      end
      
      
