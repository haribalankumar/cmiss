      SUBROUTINE HMT_SPM_CELL(TIME,Y,DY,CONTROL,MODEL,SIZES,
     '  VARIANT,DERIVED,PARAM,PROTOCOL,AII,AIO,ARI,ARO,ERR_CODE)

C#### Subroutine: HMT_SPM_CELL
C###  Description:
C###    Solve using the HMT cardiac mechanics model coupled to the
C###    Colorado of calcium kinectics.

      IMPLICIT NONE

      INCLUDE 'cell_hmt_spm.inc'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER SIZES(11)
      INTEGER CONTROL(*),MODEL(*),VARIANT,AII(*),AIO(*),ERR_CODE
      REAL*8 TIME(*),Y(*),PARAM(*),PROTOCOL(*),ARI(*),DY(*),
     '  DERIVED(*),ARO(*)
!     Local variables
      REAL*8 ra,ro,rr,koff_s,koff_1,koff_2,test1,test2,test3,
     '  test4,test5,test6,pa1,pa2,pa3,pa4,pa5,pa6,pa7,pb1,
     '  pb2,pb3,pb4,pb5,denom,nom1n,NOM1D,nom2,nom3,denom1,denom2
     '  ,nom1
      REAL*8 stimulus
      LOGICAL DEBUG

      ERR_CODE=0

      DEBUG=.FALSE.

      ! Set the stimulus - kind of ???
      stimulus = 0.0d0
      IF(TIME(TCell).GE.PROTOCOL(Is1start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is1stop))
     '  stimulus = PARAM(kl)
      IF(TIME(TCell).GE.PROTOCOL(Is2start).AND.
     '  TIME(TCell).LT.PROTOCOL(Is2stop))
     '  stimulus = PARAM(kl)
      IF((TIME(TCell).GE.PROTOCOL(IsFreqStart).AND.
     '  TIME(TCell).LT.PROTOCOL(IsFreqStop)).AND.
     '  (PROTOCOL(IsFreqPeriod).GT.1.0d-6)) THEN
        IF(DMOD(TIME(TCell),PROTOCOL(IsFreqPeriod)).LT.
     '    PROTOCOL(IsFreqDuration))
     '    stimulus = PARAM(kl)
      ENDIF

      DY(Vm) = 0.0d0

      ra=(1.0d0-Y(cas1))*(1-Y(cas2))
      ro=Y(cas1)*(1-Y(cas2))
      rr=Y(cas2)*(1-Y(cas1))
      koff_s=PARAM(Kbs)*PARAM(kons)
      koff_1=PARAM(K1)*PARAM(kon1)
      koff_2=PARAM(K2)*PARAM(kon2)

      ! Bidirectional feedback between CQSN and SRRCs ????
      IF (ro.GT.0.15d0) koff_s = koff_s * 20.0d0
      IF (Y(CaCSQ).LT.0.5d0) koff_2 = koff_2 * 0.1d0

      DY(CaCSQ)=PARAM(kons)*Y(Cas)*(PARAM(Bmaxs)-Y(CaCSQ))-
     '  (koff_s*Y(CaCSQ))

      DY(Cas1)=PARAM(kon1)*Y(Caf)*(1.0d0-Y(Cas1))-
     '  ((koff_1**2.0d0)*Y(Cas1)/(PARAM(kon1)*Y(Caf)))

      DY(Cas2)=(PARAM(kon2)*Y(Caf)*(1.0d0-Y(Cas2)))-
     '  (koff_2*Y(Cas2))

      DY(Cai)=(PARAM(kf)*(Y(Caf)-Y(Cai))-
     '  (
     '  (PARAM(Vmaxs)*((Y(Cai)**2.0d0)-
     '  ((Y(Cas)/7000.0d0)**2.0d0)))/
     '  ((PARAM(Kms)**2.0d0)+(Y(Cai)**2.0d0)+
     '  ((Y(Cas)/7000.0d0)**2.0d0))
     '  ))/  !this line and above is the nominator
     '  ((PARAM(Bmaxc)*PARAM(kbc)/
     '  ((Y(Cai)+PARAM(kbc))**2.0d0))+
     '  (PARAM(dyei)*PARAM(kbdye)/
     '  ((Y(Cai)+PARAM(kbdye))**2.0d0))+PARAM(Vc))

      DY(Caf)=((PARAM(ks)*ro*(Y(Cas)-Y(Caf)))
     '  -(PARAM(Rt)*(DY(cas1)+DY(cas2)))
     '  -(PARAM(kf)*(Y(Caf)-Y(Cai)))
     '  +(stimulus*(PARAM(Cao)-Y(Caf)))
     '  -(PARAM(VmaxNaCaX)*Y(Caf)/
     '  (PARAM(KmNaCaX)+Y(Caf))))/   !numerator this line and above
     '  ((PARAM(Bmaxf1)*PARAM(kbf1)/
     '  ((Y(Caf)+PARAM(kbf1))**2.0d0))
     '  +
     '  (PARAM(Bmaxf2)*PARAM(kbf2)/
     '  ((Y(Caf)+PARAM(kbf2))**2.0d0))
     '  + PARAM(Vf))

      DY(Cas)=(((PARAM(Vmaxs)*((Y(Cai)**2.0d0)
     '  -((Y(Cas)/7000.0d0)**2.0d0)))/
     '  ((PARAM(Kms)**2.0d0)+(Y(Cai)**2.0d0)+
     '  ((Y(Cas)/7000.0d0)**2.0d0)))-PARAM(ks)*ro*(Y(Cas)-
     '  Y(Caf)))/(PARAM(Vs)-DY(CaCSQ))


      RETURN
      END


