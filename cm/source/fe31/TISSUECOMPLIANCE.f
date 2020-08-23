      SUBROUTINE TISSUECOMPLIANCE(NEELEM,NPNE,NXI,BBM,
     &  CE,CW,dPl,undef,FIRST,ERROR,*)
      IMPLICIT NONE
      INCLUDE 'geom00.cmn' 
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn' 
!     Parameter List
      INTEGER NEELEM(0:NE_R_M),NPNE(NNM,NBFM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 BBM(2,NEM),CE(NMM,NEM),CW,dPl,undef
      LOGICAL FIRST
      CHARACTER ERROR*(*)

      INTEGER ncount,ne,noelem,NTERMS
      REAL*8 a,b,c,cc,CWW,exp_term,lambda,lambda15,min_E,ratio,MAX_PPL,
     &  MIN_PPL,totalC

      CALL ENTERS('TISSUECOMPLIANCE',*9999)

C.....dV/dP=1/[(1/2h^2).c/2.(3a+b)exp().(4h(h^2-1)^2)+(h^2+1)/h^2)]

      a=SEDF_COEFFS(1)  !dimensionless
      b=SEDF_COEFFS(2) !dimensionless
      c=SEDF_COEFFS(3) !Pa
      ncount=0
      NTERMS=0
      totalC=0.d0

      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
C.......calculate a compliance for the tissue unit (if terminal)        
        IF(NXI(1,0,ne).EQ.0)THEN
          ratio=BBM(1,ne)/undef
          
          lambda = ratio**(1.d0/3.d0) !uniform extension ratio
          cc=c
          NTERMS=NTERMS+1
          exp_term=DEXP(0.75d0*(3.d0*a+b)*(lambda**2-1.d0)**2)
          CE(nm_C,ne)=cc*exp_term/6.d0*(3.d0*(3.d0*a+b)**2
     &         *(lambda**2-1.d0)**2/lambda**2+(3.d0*a+b)
     &         *(lambda**2+1.d0)/lambda**4)
          CE(nm_C,ne)=undef/CE(nm_C,ne) !in units of volume/pressure

          CE(nm_Ppl,ne)=cc/2.d0*(3.d0*a+b)*(lambda**2.d0
     &        -1.d0)*exp_term/lambda

          totalC=totalC+CE(nm_C,ne)

        ENDIF
      ENDDO

c      CWW = CW/DBLE(NTERMS)

c      DO noelem=1,NEELEM(0)
c        ne=NEELEM(noelem)
c        IF(NXI(1,0,ne).EQ.0)THEN
C.........Add in the chest wall compliance
c          CE(nm_C,ne)=1.d0/(1.d0/CE(nm_C,ne)+1.d0/CWW)
c        ENDIF
c      ENDDO

      
      CALL EXITS('TISSUECOMPLIANCE')
      RETURN
 9999 CALL ERRORS('TISSUECOMPLIANCE',ERROR)
      CALL EXITS('TISSUECOMPLIANCE')
      RETURN 1
      END
