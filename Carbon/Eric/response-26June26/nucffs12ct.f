      SUBROUTINE NUCFFS12CT(A,Z,Q2,state,FF)
CCCC  Returns the 12C nuclear scattering Transverse form factors for the following     CCCC
CCCC  State = 13 - 15.1 MeV, 14 - 16.1 MeV, 15 - 16.6 MeV, 16 - 18.1 MeV               CCCC
CCCC  17 - 19 MeV, 18 - 20.6 MeV, 19  - GDR, 21-26 MeV (mean at 22?)                   CCCC
CCCC  20- GDR, 26-37 MeV (mean at 30?)                                                 CCCC
      
      IMPLICIT NONE

      Real*8 A,Z,Q2p,Q2,q2f,qf,FF,FF2,x2,alp,char,g1,g2,g3,e1,Radius
      Real*8 a0,a1,alph,g0,h0
      integer state,IA,IZ

      IA = int(A)
      IZ = int(Z)
      q2f = Q2/0.1975/0.1975    !!! Q^2 in inverse Fermi  !!!
      qf = sqrt(q2f)

      FF = 0.0
    !!!  Arie Bodek's 2021 fit
      
      if(state.EQ.14) then       !!!  15.1 MeV
        g1 = 2.5E-4*exp(-1.0*(qf-0.63)**2/0.4/0.4)
        g2 = 2.8E-4*exp(-1.0*(qf-0.84)**2/0.2/0.2)
        g3 = 2.4E-5*exp(-1.0*(qf-2.0)**2/0.5/0.5)
        g3 = 0.0
        e1 = -2.5E-5*exp(-1.0*qf)
        FF2 = g1+g2+g3+e1

c        FF2 = FF2*2.0
        
c        FF2 = 0.0
         
      elseif(state.EQ.15) then   !!!   16.1 MeV
c        g1 = 3.5E-4*exp(-1.0*(qf-1.3)**2/0.55/0.55)
        
        g1 = 5.9E-4*exp(-1.0*(qf-1.2)**2/0.55/0.55)
        g2 = 2.4E-4*exp(-1.0*(qf-2.2)**2/0.6/0.6)
        g3 = 0.0
c       g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1

c        FF2 = 0.0

        
      elseif(state.EQ.16) then   !!!   16.6 MeV
        g1 = 2.6E-4*exp(-1.0*(qf-1.6)**2/0.6/0.6)
        g2 = 5.0E-5*exp(-1.0*(qf-2.5)**2/0.35/0.35)
        g3 = 0.0
c       g3 = 0.0
        e1 = 0.0
        
        FF2 = g1+g2+g3+e1
        
      elseif(state.EQ.17) then   !!!   18.1 MeV
        g1 = 2.1E-4*exp(-1.0*(qf-0.8)**2/0.35/0.35)
        g2 = 1.6E-4*exp(-1.0*(qf-1.2)**2/0.425/0.425)        
c        g2 = 8.0E-6*exp(-1.0*(qf-2.1)**2/0.4/0.4)
c        g2  = 0.0
        g3 = 0.0
c       g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
    
      elseif(state.EQ.18) then   !!!   19.3 MeV
        g1 = 9.5E-4*exp(-1.0*(qf-1.27)**2/0.77/0.77)
        g2 = 3.5E-4*exp(-1.0*(qf-1.7)**2/0.6/0.6)
        g3 = 1.0E-4*exp(-1.0*(qf-2.2)**2/0.3/0.3)
c       g3 = 0.0
        e1 = -3.6E-4*exp(-1.0*qf)
        FF2 = g1+g2+g3+e1


      elseif(state.EQ.19) then   !!!   20.6 MeV
        g1 = 2.1E-4*exp(-1.0*(qf-1.45)**2/0.5/0.5)
        g2 = 5.5E-5*exp(-1.0*(qf-2.1)**2/0.4/0.4)
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
c        FF2 = FF2/(1.0+5.08*qf*qf*qf)
        
      elseif(state.EQ.20) then   !!!   23 MeV GDR
        g1 = 1.8E-3*exp(-1.0*(qf-0.8)**2/0.36/0.36)
        g2 = 0.0E-4*exp(-1.0*(qf-1.5)**2/0.6/0.6)

        g2 = 1.0E-4*exp(-1.0*(qf-1.5)**2/0.5/0.5)

c        g1 = 0.0
c        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1

c        if(qf.GT.0.8) FF2 = 0.8*FF2 !!!   TEST!!!

        
      elseif(state.EQ.21) then  !!!   31.5 MeV GDR
        g1 = 9.0E-4*exp(-1.0*(qf-0.35)**2/0.3/0.3)
c        g2 = 3.5E-3*exp(-1.0*(qf-1.5)**2/0.4/0.4)
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1

c        FF2 = 1.3*FF2

        
c        FF2 = FF2*/(1.0+0.5*qf*qf)

c        FF2 = 0.0
        
      endif        
      if(qf*qf.GT.12.0) FF2 = 0.0
      if(FF2.LT.0.0) FF2 = 0.0
      FF = sqrt(FF2)
      
c      if(state.EQ.10) write(6,*) state,qf,FF2
c      if(state.EQ.4) write(6,*) state, q2,q2f,FF2,g1,g2,g3,e1

      return
      end
      
      
