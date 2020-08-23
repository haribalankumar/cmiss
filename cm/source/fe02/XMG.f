      SUBROUTINE XMG(IBT,IDO,INP,NBJ,nr,GL,GU,XE,XI,ERROR,*)

C#### Subroutine: XMG
C###  Description:
C###    XMG evaluates the covariant (GL) and contravariant (GU) metric
C###    tensors wrt the Xi-coordinate system  and  the deriv.s  of
C###    the Xi-coords wrt the Xj-coords (DXIX) -if NIT=NJT only - at
C###    current Xi point.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NBJ(NJM),nr
      REAL*8 DXXI(3,0:3),G,GL(3,*),GU(3,3),XI(*),XE(NSM,NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER mi,nb,ni,NITB,nj,nu
      REAL*8 AA,PXI,G1,G3,R,RC,RR,RRC,SLX,SMX

      CALL ENTERS('XMG',*9999)
      NITB=NIT(NBJ(1))
      DO ni=0,NITB
        nu=1+ni*(1+ni)/2
        DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
          nb=NBJ(nj)
          DXXI(nj,ni)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,nu,
     '      XI,XE(1,nj))
        ENDDO !nj
      ENDDO !ni
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' DXXI:'',12E10.3)')
     '    ((DXXI(nj,ni),ni=0,NITB),nj=1,NJ_LOC(NJL_GEOM,0,nr))
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      IF(ITYP10(1).EQ.2) THEN
        R=DXXI(1,0)
        RR=R*R
      ELSE IF(ITYP10(1).EQ.3) THEN
        R=DXXI(1,0)
        RR=R*R
        RC=R*DCOS(DXXI(3,0))
        RRC=RC*RC
      ELSE IF(ITYP10(1).EQ.4) THEN
        AA=FOCUS*FOCUS
        SLX=DSINH(DXXI(1,0))
        SMX=DSIN(DXXI(2,0))
        G1=AA*(SLX*SLX+SMX*SMX)
        G3=AA* SLX*SLX*SMX*SMX
      ELSE IF(ITYP10(1).EQ.5) THEN
      ENDIF
      DO mi=1,NITB
        DO ni=1,NITB
          IF(ITYP10(1).NE.4) GL(mi,ni)=DXXI(1,mi)*DXXI(1,ni)
          IF(NJ_LOC(NJL_GEOM,0,nr).GT.1) THEN
            IF(ITYP10(1).EQ.1) THEN
              DO nj=2,NJ_LOC(NJL_GEOM,0,nr)
                GL(mi,ni)=GL(mi,ni)+DXXI(nj,mi)*DXXI(nj,ni)
              ENDDO !nj
            ELSE IF(ITYP10(1).EQ.2) THEN
              GL(mi,ni)=GL(mi,ni)+RR*DXXI(2,mi)*DXXI(2,ni)
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) THEN
                GL(mi,ni)=GL(mi,ni)+DXXI(3,mi)*DXXI(3,ni)
              ENDIF
            ELSE IF(ITYP10(1).EQ.3) THEN
              GL(mi,ni)=GL(mi,ni)+RRC*DXXI(2,mi)*DXXI(2,ni)
     '          +RR *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(1).EQ.4) THEN
              GL(mi,ni)=G1*(DXXI(1,mi)*DXXI(1,ni)
     '          +DXXI(2,mi)*DXXI(2,ni))
              IF(NJ_LOC(NJL_GEOM,0,nr).EQ.3) GL(mi,ni)=GL(mi,ni)+G3
     '          *DXXI(3,mi)*DXXI(3,ni)
            ELSE IF(ITYP10(1).EQ.5) THEN
            ENDIF
          ENDIF
        ENDDO !ni
      ENDDO !mi

      CALL INVERT(NITB,GL,GU,G)


      CALL EXITS('XMG')
      RETURN
 9999 CALL ERRORS('XMG',ERROR)
      CALL EXITS('XMG')
      RETURN 1
      END


