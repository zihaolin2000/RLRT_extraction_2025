c--------------------------------- PROGRAM coulomb.F ------------------------------
c
c- P. Solvignon, Dec. 7, 2006
c
c  Coulomb corrections based on D. Gaskell's kumac (correct_emc.kumac)
c  
c  -- update --
c
c   Use method describe in Aste et al. : Eur. Phys. J. A26 (2005) 167
c      --> use the average Coulomb potential:
c          V = B*V0 with 0.75 < B < 0.80
c      --> the cross section is corrected by changing E-->E+V and Ep-->Ep+V
c          in the MOTT and in the spectral function expression
c      --> the corrected cross section is then multiplied by the incoming 
c          focusing factor FF squared: FF = (E+V)/E
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccccccccccccccccccccccccccccccccc
c
      SUBROUTINE VCOUL(A,Z,V)
      IMPLICIT NONE

      REAL*8 A,Z,R0
      REAL*8 C_ASTE,HBARC,V0,V,ALPHA

      HBARC  = 0.197327      ! in GeV.fm
      ALPHA  = 1.0/137.0
      C_ASTE = 0.775
      R0     = 1.1*A**(1./3.) + 0.86*A**(-1./3.)

      ! Coulomb potential at the center of the nucleus
      V0  = (3./2.)*ALPHA*HBARC*(Z-1.)/R0  ! in GeV

      ! Average potential
      V  = C_ASTE*V0      ! from Eur. Phys. J. A26 (2005) 167
      
      ! Ziggy added May 19: use experimentally determined Veff=3.1MeV
      !V = 0.0031

c      write(6,*) V
      
      RETURN
      END
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
