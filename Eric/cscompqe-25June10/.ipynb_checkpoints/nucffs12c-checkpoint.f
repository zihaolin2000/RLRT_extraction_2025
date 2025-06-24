      SUBROUTINE NUCFFS12C(A,Z,q32,state,FF)
CCCC  Returns the 12C nuclear scattering form factors for the following states   CCCC
CCCC  State = 1 - nuclear elastic, 2 - 4.4 MeV, 3 - 7.6 MeV, 4 - 9.6 MeV         CCCC
CCCC  5 - 10.84 MeV, 6 -15.1 MeV state, 7 - 16.1 MeV, 8 - 18.6 MeV, 9 - 20 MeV   CCCC
CCCC  10 - GDR, 21-26 MeV (mean at 23?)                                        CCCC
CCCC  11- GDR, 26-37 MeV (mean at 31.5?), 12 - ~15.1 MeV state                    CCCC
      
      IMPLICIT NONE

      Real*8 A,Z,q3,qf,Q2p,Q2,q3f,q32,q2f,FF,FF2,x2,alp,char,g1,g2,g3,g4
      Real*8 e1,Radius,a0,a1,alph,g0,h0,a2,a3,a4,a5,b
      integer state,IA,IZ

      
      
      IA = int(A)
      IZ = int(Z)
      q2f = q32/0.1975/0.1975 
      qf = sqrt(q2f) !!! q3 in fm^-1
      Radius = 2.5 
      Radius = 2.45
      
      FF = 0.0
      if(state.EQ.1) then
        x2 = (Q2/0.197328**2)*Radius**2
        alp = (Z-2.)/3.                                               
        char = x2*(2.+3.*alp)/(12.+30.*alp)                             
        if (char.lt.80) FF = exp(-char)*(1.-alp*x2/(6.+15.*alp))
        FF2 = FF*FF
cc        write(6,*) "1: ", q2f,FF2

        
CCCC  Arie Bodek's fit using the harmomic well shape from Hofstadter  CCC

        alph = 4.0/3.0
        a0 = 1.65
        g0 = 8.0E-5*exp(-1.0*((q2f-2.9)/0.44)**2.0)
        h0 = (1.0-alph*q2f*a0*a0/2.0/(2.0+3.0*alph))*
     &         exp(-q2f*a0*a0/4.0)
        
c        a1 = a0*sqrt(3.0/2.0*(2.0+5.0*alph)/(2.0+3.0*alph))
c        FF = 1.0-a1*q2f*a0*a0/2.0/(2.0+3.0*alph)          
c        FF = FF*exp(-1.0*q2f*a0*a0/4.0)
c        FF = FF+3.0E-5*exp((q2f-3.5)/2.0)**2.
c        FF2 = FF*FF

        g4 = 1.0E-5*exp(-1.0*((q2f-4.0)/1.2)**2.0)
        
        FF2 = h0*h0+g0+g4

c        FF2 = 2.5*FF2

        
c        write(6,*) "2: ", q2f,FF2
        
CCCC  The excited states are from a fit by Arie Bodek, 2021  CCCC
        
      elseif(state.EQ.2) then   !!!  4.43 MeV
        g1 = 1.41E-2*exp(-1.0*(q2f-1.125)**2/1.71/1.71)
        g2 = 7.0E-4*exp(-1.0*(q2f-3.7)**2/1.6/1.6)  
        g3 = 3.3E-6*exp(-1.0*(q2f-6.5)**2/7.0/7.0)
c     g4 = 1.9E-4*exp(-1.0*(q2f-4.25)**2/1.0/1.0)
        g4 = 0.0
        FF2 = q2f**3./(q2f**3.0+0.1)*(g1+g2+g3+g4)   !!! New fit from Feb 2023



c        write(6,*) q2f,FF2
        
      elseif(state.EQ.3) then     !!!  7.6 MeV

c        g1 = 2.8E-3*exp(-1.0*(qf-0.93)**2/0.42/0.42)
c        g2 = 3.0E-4*exp(-1.0*(qf-1.45)**2/0.24/0.24)
c        g3 = 2.0E-5*exp(-1.0*(qf-2.48)**2/0.53/0.53)
c        e1 = -1.0E-4*exp(-1.0*qf)
c        FF2 = g1+g2+g3+e1        

         b = 1.3457
         a1 = 0.52*(b*qf)**2
         a2 = -0.025*(b*qf)**4
         a3 = -0.7E-2*(b*qf)**6
         a4 = 0.5E-3*(b*qf)**8
         a5 = -0.5E-4*(b*qf)**10

         FF2 = 1/Z*exp(-0.5*b*b*q2f)*(a1+a2+a3+a4+a5)
         FF2 = FF2*FF2

c         if(abs(qf-1.0).LT.2.15) write(6,*) "*** FF2-7.6 = ", qf, FF2
         
        
      elseif(state.EQ.4) then     !!!  9.6 MeV
        g1 = 5.0E-3*exp(-1.0*(q2f-1.46)**2/1.6/1.6)
        g2 = 6.6E-4*exp(-1.0*(q2f-3.46)**2/2.0/2.0)
        g3 = 7.0E-6*exp(-1.0*(q2f-7.0)**2/2.8/2.8)
        FF2 = q2f**3./(q2f**3.0+0.2)*(g1+g2+g3)   !!! New fit from Feb 2023  

      elseif(state.EQ.5) then   !!!  10.84 MeV

        g1 = 5.0E-4*exp(-1.0*(qf-1.0)**2/0.3/0.3)
        g2 = 8.0E-4*exp(-1.0*(qf-1.4)**2/0.4/0.4)
        g3 = 0.0*exp(-1.0*(qf-7.0)**2/2.5/2.5)
        e1 = 0.0
        FF2 = g1+g2+g3+e1

      elseif(state.EQ.6) then   !!!  13.7 MeV

        g1 = 4.0E-4*exp(-1.0*(qf-1.0)**2/0.35/0.35)
        g2 = 8.0E-4*exp(-1.0*(qf-1.75)**2/0.45/0.45)
        g3 = 4.0E-4*exp(-1.0*(qf-0.85)**2/0.65/0.65)
c        g4 = 1.0E-4*exp(-1.0*(qf-1.8)**2/0.6/0.6)
        e1 = 0.0
        FF2 = g1+g2+g3+e1+g4
c        FF2 = FF2/4.0
        

      elseif(state.EQ.7) then  !!! 15.1 MeV
        g1 = 6.0E-4*exp(-1.0*(qf-0.85)**2/0.7/0.7)
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
        
      elseif(state.EQ.8) then   !!! 16.1 MeV
        g1 = 12.0E-4*exp(-1.0*(qf-1.05)**2/0.6/0.6)
c        g1 = 18.0E-4*exp(-1.0*(qf-1.15)**2/0.6/0.6)         
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
        

      elseif(state.EQ.9) then   !!! 18.6 MeV
c     g1 = 4.5E-4*exp(-1.0*(qf-1.5)**2/0.5/0.5)
        g1 = 3.2E-4*exp(-1.0*(qf-1.3)**2/0.5/0.5)
    
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
      elseif(state.EQ.10) then   !!! 20 MeV
        g1 = 1.6E-4*exp(-1.0*(qf-1.2)**2/0.42/0.42)
        g2 = 1.6E-5*exp(-1.0*(qf-1.8)**2/0.4/0.4)
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
      elseif(state.EQ.11) then   !!! GDR 21-26 MeV
        g1 = 2.8E-3*exp(-1.0*(qf-0.60)**2/0.15/0.15)
        g2 = 6.9E-3*exp(-1.0*(qf-0.84)**2/0.55/0.55)
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
        
      elseif(state.EQ.12) then   !!! GDR 26-37 MeV (31.5 MeV ?)
        g1 = 0.0047*exp(-1.0*(qf-1.0)**2/0.48/0.48)
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1
        
        
      elseif(state.EQ.13) then   !!! add 40 MeV state
        g1 = 0.0026*exp(-1.0*(qf-1.49)**2/0.7/0.7)         
        g2 = 0.0
        g3 = 0.0
        e1 = 0.0
        FF2 = g1+g2+g3+e1   
      endif


      if(qf*qf.GT.12.0) FF2 = 0.0
      if(FF2.LT.0.0) FF2 = 0.0
      FF = sqrt(FF2)

c      if(state.EQ.10) write(6,*) state,qf,FF2
c      if(state.EQ.4) write(6,*) state, q2,q2f,FF2,g1,g2,g3,e1

      return
      end
      
      
