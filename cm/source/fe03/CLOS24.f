      SUBROUTINE CLOS24(IBT,IDO,INP,IT,ITMAX,NBJ,SQ,XE,XI,ZD,ERROR,*)

C#### Subroutine: CLOS24
C###  Description:
C###    CLOS24 finds the XI-coordinates at the closest approach of a
C###    face to a data point with coordinates XD using a modified
C###    Newton algorithm.

C**** ITMAX is the maximum number of iterations.
C**** TOL is the required tolerance of the solution.
C**** VMAX denotes the maximum step length per iteration.
C**** NOTE: Works only for 2-D elements in prolate spheroidal coords.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),IT,ITMAX,NBJ(NJM)
      REAL*8 SQ,XE(NSM,NJM),XI(3),ZD(NJM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER it1,it2,MI,nb,ni,nj
      REAL*8 CCC,CCO,CCS,COO,CSC,CSO,CSS,D2SQXI(2,2),D2XXI(3,2,2),
     '  D2ZXI(3,2,2),DELTA,DET,DSQXI(2),DSQXIV,DXI(2),DXXI(3,2),
     '  DZ(3),DZXI(3,2),OCO,OOC,OOS,OSO,PXI,SCC,SCO,SCS,SOO,SQLIN,
     '  SQOLD,SSC,SSO,SSS,TOL,TOL2,V(2),V2,VMAX,VMAX2,W,X(3),XILIN(2),
     '  Z(3)
      CHARACTER FORMAT*200
      DATA TOL /1.E-3/, VMAX /0.25/

      CALL ENTERS('CLOS24',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(A)') ' >CLOS24 2-D prolate spheroid'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      TOL2=TOL**2
      VMAX2=VMAX**2
      DO nj=1,NJT
        nb=NBJ(nj)
        X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,
     '    XE(1,nj))
      ENDDO
      COO=DCOSH(X(1))*FOCUS
      SOO=DSINH(X(1))*FOCUS
      OCO=DCOS(X(2))
      OSO=DSIN(X(2))
      OOC=DCOS(X(3))
      OOS=DSIN(X(3))
      CCO=COO*OCO
      CSO=COO*OSO
      SCO=SOO*OCO
      SSO=SOO*OSO
      CCC=COO*OCO*OOC
      CCS=COO*OCO*OOS
      CSC=COO*OSO*OOC
      CSS=COO*OSO*OOS
      SCC=SOO*OCO*OOC
      SCS=SOO*OCO*OOS
      SSC=SOO*OSO*OOC
      SSS=SOO*OSO*OOS
      Z(1)=CCO
      Z(2)=SSC
      Z(3)=SSS
      DZ(1)=Z(1)-ZD(1)
      DZ(2)=Z(2)-ZD(2)
      DZ(3)=Z(3)-ZD(3)
      SQ=DZ(1)**2+DZ(2)**2+DZ(3)**2
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        FORMAT='(/'' ZD(nj)='',3(E12.6,4X),/'' XI(ni)='',2(E12.6,4X),'
     '        //'/''  Z(nj)='',3(E12.6,4X),/''     SQ='',E12.6)'
        WRITE(OP_STRING,FORMAT) (ZD(nj),nj=1,NJT),(XI(ni),ni=1,2),
     '    (Z(nj),nj=1,NJT),SQ
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO it1=1,ITMAX
        DO nj=1,NJT
          nb=NBJ(nj)
          DXXI(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,2,XI,XE(1,nj))
          DXXI(nj,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,4,XI,XE(1,nj))
          D2XXI(nj,1,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,3,XI,XE(1,nj))
          D2XXI(nj,1,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,6,XI,XE(1,nj))
          D2XXI(nj,2,2)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '      nb,5,XI,XE(1,nj))
        ENDDO
        DZXI(1,1)=DXXI(1,1)*SCO-DXXI(2,1)*CSO
        DZXI(1,2)=DXXI(1,2)*SCO-DXXI(2,2)*CSO
        DZXI(2,1)=DXXI(1,1)*CSC+DXXI(2,1)*SCC-DXXI(3,1)*SSS
        DZXI(2,2)=DXXI(1,2)*CSC+DXXI(2,2)*SCC-DXXI(3,2)*SSS
        DZXI(3,1)=DXXI(1,1)*CSS+DXXI(2,1)*SCS+DXXI(3,1)*SSC
        DZXI(3,2)=DXXI(1,2)*CSS+DXXI(2,2)*SCS+DXXI(3,2)*SSC
        D2ZXI(1,1,1)=D2XXI(1,1,1)*SCO-D2XXI(2,1,1)*CSO
     '              +DXXI(1,1)*DXXI(1,1)*CCO-2.d0*DXXI(1,1)*DXXI(2,1)*
     '               SSO
     '              -DXXI(2,1)*DXXI(2,1)*CCO
        D2ZXI(1,1,2)=D2XXI(1,1,2)*SCO-D2XXI(2,1,2)*CSO
     '              +DXXI(1,1)*DXXI(1,2)*CCO-2.d0*DXXI(1,1)*DXXI(2,2)*
     '               SSO
     '              -DXXI(2,1)*DXXI(2,2)*CCO
        D2ZXI(1,2,2)=D2XXI(1,2,2)*SCO-D2XXI(2,2,2)*CSO
     '              +DXXI(1,2)*DXXI(1,2)*CCO-2.d0*DXXI(1,2)*DXXI(2,2)*
     '               SSO
     '              -DXXI(2,2)*DXXI(2,2)*CCO
        D2ZXI(2,1,1)=D2XXI(1,1,1)*CSC+D2XXI(2,1,1)*SCC-D2XXI(3,1,1)*SSS
     '              +DXXI(1,1)*DXXI(1,1)*SSC+2.d0*DXXI(1,1)*DXXI(2,1)*
     '               CCC
     '              -DXXI(2,1)*DXXI(2,1)*SSC-2.d0*DXXI(2,1)*DXXI(3,1)*
     '               SCS
     '              -DXXI(3,1)*DXXI(3,1)*SSC-2.d0*DXXI(3,1)*DXXI(1,1)*
     '               CSS
        D2ZXI(2,1,2)=D2XXI(1,1,2)*CSC+D2XXI(2,1,2)*SCC-D2XXI(3,1,2)*SSS
     '              +DXXI(1,1)*DXXI(1,2)*SSC+2.d0*DXXI(1,1)*DXXI(2,2)*
     '               CCC
     '              -DXXI(2,1)*DXXI(2,2)*SSC-2.d0*DXXI(2,1)*DXXI(3,2)*
     '               SCS
     '              -DXXI(3,1)*DXXI(3,2)*SSC-2.d0*DXXI(3,1)*DXXI(1,2)*
     '               CSS
        D2ZXI(2,2,2)=D2XXI(1,2,2)*CSC+D2XXI(2,2,2)*SCC-D2XXI(3,2,2)*SSS
     '              +DXXI(1,2)*DXXI(1,2)*SSC+2.d0*DXXI(1,2)*DXXI(2,2)*
     '               CCC
     '              -DXXI(2,2)*DXXI(2,2)*SSC-2.d0*DXXI(2,2)*DXXI(3,2)*
     '               SCS
     '              -DXXI(3,2)*DXXI(3,2)*SSC-2.d0*DXXI(3,2)*DXXI(1,2)*
     '               CSS
        D2ZXI(3,1,1)=D2XXI(1,1,1)*CSS+D2XXI(2,1,1)*SCS+D2XXI(3,1,1)*SSC
     '              +DXXI(1,1)*DXXI(1,1)*SSS+2.d0*DXXI(1,1)*DXXI(2,1)*
     '               CCS
     '              -DXXI(2,1)*DXXI(2,1)*SSS+2.d0*DXXI(2,1)*DXXI(3,1)*
     '               SCC
     '              -DXXI(3,1)*DXXI(3,1)*SSS+2.d0*DXXI(3,1)*DXXI(1,1)*
     '               CSC
        D2ZXI(3,1,2)=D2XXI(1,1,2)*CSS+D2XXI(2,1,2)*SCS+D2XXI(3,1,2)*SSC
     '              +DXXI(1,1)*DXXI(1,2)*SSS+2.d0*DXXI(1,1)*DXXI(2,2)*
     '               CCS
     '              -DXXI(2,1)*DXXI(2,2)*SSS+2.d0*DXXI(2,1)*DXXI(3,2)*
     '               SCC
     '              -DXXI(3,1)*DXXI(3,2)*SSS+2.d0*DXXI(3,1)*DXXI(1,2)*
     '               CSC
        D2ZXI(3,2,2)=D2XXI(1,2,2)*CSS+D2XXI(2,2,2)*SCS+D2XXI(3,2,2)*SSC
     '              +DXXI(1,2)*DXXI(1,2)*SSS+2.d0*DXXI(1,2)*DXXI(2,2)*
     '               CCS
     '              -DXXI(2,2)*DXXI(2,2)*SSS+2.d0*DXXI(2,2)*DXXI(3,2)*
     '               SCC
     '              -DXXI(3,2)*DXXI(3,2)*SSS+2.d0*DXXI(3,2)*DXXI(1,2)*
     '               CSC
        DSQXI(1)=DZXI(1,1)*DZ(1)+DZXI(2,1)*DZ(2)+DZXI(3,1)*DZ(3)
        DSQXI(2)=DZXI(1,2)*DZ(1)+DZXI(2,2)*DZ(2)+DZXI(3,2)*DZ(3)
        D2SQXI(1,1)=DZXI(1,1)*DZXI(1,1)+D2ZXI(1,1,1)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,1)+D2ZXI(2,1,1)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,1)+D2ZXI(3,1,1)*DZ(3)
        D2SQXI(1,2)=DZXI(1,1)*DZXI(1,2)+D2ZXI(1,1,2)*DZ(1)
     '             +DZXI(2,1)*DZXI(2,2)+D2ZXI(2,1,2)*DZ(2)
     '             +DZXI(3,1)*DZXI(3,2)+D2ZXI(3,1,2)*DZ(3)
        D2SQXI(2,2)=DZXI(1,2)*DZXI(1,2)+D2ZXI(1,2,2)*DZ(1)
     '             +DZXI(2,2)*DZXI(2,2)+D2ZXI(2,2,2)*DZ(2)
     '             +DZXI(3,2)*DZXI(3,2)+D2ZXI(3,2,2)*DZ(3)
        DET=D2SQXI(1,1)*D2SQXI(2,2)-D2SQXI(1,2)**2
        IF(DABS(DET).LE.TOL) THEN
          V(1)=-DSQXI(1)
          V(2)=-DSQXI(2)
        ELSE
          V(1)=-(D2SQXI(2,2)*DSQXI(1)-D2SQXI(1,2)*DSQXI(2))/DET
          V(2)=-(DSQXI(2)+D2SQXI(1,2)*V(1))/D2SQXI(2,2)
          DSQXIV=DSQXI(1)*V(1)+DSQXI(2)*V(2)
          DELTA=DSQRT((DSQXI(1)**2+DSQXI(2)**2)*(V(1)**2+V(2)**2))*TOL
          IF(DABS(DSQXIV).LE.DELTA) THEN
            V(1)=-DSQXI(1)
            V(2)=-DSQXI(2)
          ELSE IF(DSQXIV.GT.DELTA) THEN
            V(1)=-V(1)
            V(2)=-V(2)
          ENDIF
        ENDIF
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          FORMAT='(''    IT1='',I4,10X,''DSQXI(ni)='',2(E12.6,4X),'
     '      //'/18X,''D2SQXI(mi,ni)='',2(E12.6,4X),/32X,2(E12.6,4X),'
     '      //'/26X,''V(ni)='',2(E12.6,4X))'
          WRITE(OP_STRING,FORMAT) it1,(DSQXI(ni),ni=1,2),
     '      ((D2SQXI(MI,ni),MI=1,2),ni=1,2),(V(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        W=1.d0
        V2=V(1)**2+V(2)**2
        IF(V2.GE.VMAX2) W=VMAX/DSQRT(V2)
        DO it2=1,ITMAX
          XILIN(1)=XI(1)+V(1)*W
          XILIN(2)=XI(2)+V(2)*W
          DO nj=1,NJT
            nb=NBJ(nj)
            X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '        XILIN,XE(1,nj))
          ENDDO
          COO=DCOSH(X(1))*FOCUS
          SOO=DSINH(X(1))*FOCUS
          OCO=DCOS(X(2))
          OSO=DSIN(X(2))
          OOC=DCOS(X(3))
          OOS=DSIN(X(3))
          CCO=COO*OCO
          CSO=COO*OSO
          SCO=SOO*OCO
          SSO=SOO*OSO
          CCC=COO*OCO*OOC
          CCS=COO*OCO*OOS
          CSC=COO*OSO*OOC
          CSS=COO*OSO*OOS
          SCC=SOO*OCO*OOC
          SCS=SOO*OCO*OOS
          SSC=SOO*OSO*OOC
          SSS=SOO*OSO*OOS
          Z(1)=CCO
          Z(2)=SSC
          Z(3)=SSS
          DZ(1)=Z(1)-ZD(1)
          DZ(2)=Z(2)-ZD(2)
          DZ(3)=Z(3)-ZD(3)
          SQLIN=DZ(1)**2+DZ(2)**2+DZ(3)**2
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            FORMAT='(8X,''IT2='',I4,14X,''W='',E12.6,/21X,'
     '        //'''XILIN(ni)='','//'2(E12.6,4X)/26X,''Z(nj)='','
     '        //'3(E12.6,4X)/26X,''SQLIN='',E12.6)'
            WRITE(OP_STRING,FORMAT) it2,W,(XILIN(ni),ni=1,2),
     '        (Z(nj),nj=1,NJT),SQLIN
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          IF(SQLIN.LE.SQ) GO TO 5
          W=(DSQXIV*W*W)/(DSQXIV*W+SQ-SQLIN)
        ENDDO
    5   SQOLD=SQ
        SQ=SQLIN
        DXI(1)=V(1)*W
        DXI(2)=V(2)*W
        XI(1)=XI(1)+DXI(1)
        XI(2)=XI(2)+DXI(2)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'(24X,''Xi(ni)='',2(E12.6,4X))')
     '      (XI(ni),ni=1,2)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IT=it1
        IF(XI(1).LT.-0.5D0.OR.XI(1).GT.1.5D0.OR
     '    .XI(2).LT.-0.5D0.OR.XI(2).GT.1.5D0) GO TO 9998
        IF(((DXI(1)**2+DXI(2)**2)/(1.d0+XI(1)**2+XI(2)**2).LE.TOL2)
     '    .AND.((SQOLD-SQ)/(1.d0+SQ).LE.TOL)) GO TO 9998
      ENDDO

 9998 CALL EXITS('CLOS24')
      RETURN
 9999 CALL ERRORS('CLOS24',ERROR)
      CALL EXITS('CLOS24')
      RETURN 1
      END


