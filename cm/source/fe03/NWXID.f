      SUBROUTINE NWXID(ITMAX,l1,LD,LN,NBJ,NBJF,NDDL,NDLT,NKJE,
     '  NKEF,NNF,NPF,NPNE,NPNF,NRE,NVJE,NVJF,NXI,SE,SF,SQ,
     '  XA,XE,XID,XP,ZD,ERROR,*)
C     SMAR009 18/01/99 removed NPL from list
C#### Subroutine: NWXID
C###  Description:
C###    NWXID calculates Xi coordinates corresponding to closest
C###    approach of a data point to a region defined by XP.

C**** If any Xi lie outside [0,1] the region is redefined by XP for a
C**** segment, if it exists, corresponding to the Xi coordinate most
C**** violated in the previous segment.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER ITMAX,LD(NDM),LN(0:NEM),NBJ(NJM,NEM),NBJF(NJM,NFM),
     '  NDDL(NEM,NDEM),NDLT(NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NKEF(0:4,16,6,NBFM),NNF(0:17,6,NBFM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNF(NNM,NBFM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),NVJF(NNM,NBFM,NJM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
C SMAR009 19/01/99 removed NPL(5,0:3,NLM),
      REAL*8 SE(NSM,NBFM,NEM),SF(NSM,NBFM),SQ(NDM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),
     '  XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),ZD(NJM,NDM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IT,l1,L2,l3,LMAX,nd,ndf1,ndl1,ne,nef,NF1,NF2,ni,NIOUT(4),
     '  nj,NJOT,NKJF(NKM,NNM,NJM),NL1,nou,NOUT,nr
      REAL*8 ABSOUT(8)

C LKC 31-OCT_97 unused variables
C     INTEGER njj,njj1,njj2,

      CALL ENTERS('NWXID',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(A)') ' >NWXID'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
CGMH29/8/95 Initialise IT (even though it is not really used)
      IT=0
C     DO 9 l1=1,LN(0)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' L1='',I4,'' LN(L1)='',I4)') l1,LN(l1)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(LN(l1).GT.0) THEN
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' LXI='',5I4)')
     '        ((NXI(ni,1,LMAX),ni=-2,2),LMAX=1,LN(0))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
          NF1=LN(l1)
          NJOT=3
          DO 5 ndf1=1,NDLT(NF1)
            nd=NDDL(NF1,ndf1)
            L2=l1
            DO 4 l3=1,LN(0)
              NF2=LN(L2)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' NDF1='',I4,'' nd='',I4,'
     '            //''' L2='',I4,'' LN(L2)='',I4)') ndf1,nd,L2,NF2
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              ne=NPF(6,nf2)
              nef=NPF(8,nf2)
              nr=NRE(ne)
              CALL CALC_FACE_INFORMATION_IND(NBJ(1,ne),NBJF(1,nf2),
     '          nef,NKJE(1,1,1,ne),NKEF,NKJF,NNF,NPNE(1,1,ne),NPNF,nr,
     '          NVJE(1,1,1,ne),NVJF,SE(1,1,ne),SF,ERROR,*9999)
              CALL XPXE(NBJF(1,nf2),NKJF,NPF(1,NF2),NPNF,nr,NVJF,
     '          SF,XA(1,1,ne),XE,XP,ERROR,*9999)
c Note: PJH 21Nov94 The following commented calls need updating
              IF(ITYP10(nr).EQ.1) THEN
c               CALL CLOS21(IBT,IDO,INP,IT,ITMAX,NJOT,NKE(1,1,1,2),
c    '            NPF(1,NF2),NPNE(1,1,NF2),
c    '            SE(1,1,NF2),SQ(nd),TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.2) THEN
c               CALL CLOS22(IBT,IDO,INP,IT,ITMAX,NPF(1,NF2),SQ(nd),
c    '            TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.3) THEN
c               CALL CLOS23(IBT,IDO,INP,IT,ITMAX,NPF(1,NF2),SQ(nd),
c    '            TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.4) THEN
c               CALL CLOS24(IBT,IDO,INP,IT,ITMAX,NJOT,NKE(1,1,1,NF2),
c    '            NPF(1,NF2),NPNE(1,1,NF2),
c    '            SE(1,1,NF2),SQ(nd),TOL,VMAX,XE,XID(1,nd),ZD(1,nd))
              ENDIF
              IF(IT.GE.ITMAX) THEN
                WRITE(OP_STRING,
     '            '('' WARNING: Convergence not reached in CLOS2X''/'
     '            //''' nd='',I5,'' ZD='',(E12.6,1X)/'
     '            //''' LD='',I4,'' XID'',(E12.6,1X),''SQ='',E12.6/)')
     '            nd,(ZD(nj,nd),nj=1,NJOT),
     '            LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              NOUT=0
              DO 1 ni=1,2
                IF(XID(ni,nd).LT.0.d0) THEN
                  NOUT=NOUT+1
                  NIOUT(NOUT)=-ni
                  ABSOUT(NOUT)=-XID(ni,nd)
                ELSE IF(XID(ni,nd).GT.1.d0) THEN
                  NOUT=NOUT+1
                  NIOUT(NOUT)=ni
                  ABSOUT(NOUT)=XID(ni,nd)-1.d0
                ENDIF
    1         CONTINUE
              IF(NOUT.EQ.0) THEN
                LD(nd)=LN(L2)
                GO TO 5
              ENDIF
              CALL RSORT(NOUT,ABSOUT,NIOUT)
              DO 2 nou=NOUT,1,-1
                IF(ABS(NXI(NIOUT(nou),1,L2)).NE.L2) GO TO 3
                IF(nou.EQ.1) THEN
                  LD(nd)=LN(L2)
                  GO TO 5
                ENDIF
    2         CONTINUE
    3         IF(DOP) THEN
                WRITE(OP_STRING,'('' LXI='',7I4)')
     '            (NXI(ni,1,L2),ni=-2,2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              L2=ABS(NXI(NIOUT(nou),1,L2))
              ni=ABS(NIOUT(nou))
              IF(XID(ni,nd).LT.0.d0) THEN
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' XID(ni,nd)='',E12.4)')
     '              XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                XID(ni,nd)=DMAX1(XID(ni,nd)+1.d0,0.75D0)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Changed to '',E12.4)')
     '              XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ELSE IF(XID(ni,nd).GT.1.d0) THEN
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' XID(ni,nd)='',E12.4)')
     '              XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
                XID(ni,nd)=DMIN1(XID(ni,nd)-1.d0,0.25D0)
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,'('' Changed to '',E12.4)')
     '              XID(ni,nd)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF
              ENDIF
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' L2='',I2,'' NIOUT='',I2,'
     '            //''' XID(ni,nd)='',E12.6)') L2,NIOUT(nou),XID(ni,nd)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
    4       CONTINUE
    5     CONTINUE
        ELSE IF(LN(l1).LT.0) THEN
          NL1=ABS(LN(l1))
          NJOT=NJT
          DO 8 ndl1=1,NDLT(NL1)
            nd=NDDL(NL1,ndl1)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' L1='',I4,'' NL1='',I4,'' NJOT='','
     '          //'I4,'' NDL1='',I4,'' nd='',I4,'' LD(nd)='',I4)')
     '          l1,NL1,NJOT,ndl1,nd,LD(nd)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            DO 7 l3=1,LN(0)
CGMH29/8/95 Unused              NL2=ABS(LN(L2))
              IF(ITYP10(nr).EQ.1) THEN
c               CALL CLOS11(IT,ITMAX,NPL(1,0,NL2),DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.2) THEN
c               CALL CLOS12(IT,ITMAX,NPL(1,0,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.3) THEN
c               CALL CLOS13(IT,ITMAX,NPL(1,0,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ELSE IF(ITYP10(nr).EQ.4) THEN
c               CALL CLOS14(IT,ITMAX,NPL(1,0,NL2),nr,DL(1,NL2),
c    '                      SQ(nd),TOL,VMAX,XP,XID(1,nd),ZD(1,nd))
              ENDIF
              IF(IT.GE.ITMAX) THEN
                WRITE(OP_STRING,
     '            '('' WARNING: Convergence not reached in CLOS1X''/'
     '            //''' nd='',I5,'' ZD='',(E12.6,1X)/'
     '            //''' LD='',I4,'' XID'',(E12.6,1X),''SQ='',E12.6/)')
     '            nd,(ZD(nj,nd),nj=1,NJOT),
     '            LD(nd),(XID(ni,nd),ni=1,1),SQ(nd)
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
              NOUT=0
              IF(XID(1,nd).LT.0.d0) THEN
                NOUT=NOUT+1
                NIOUT(NOUT)=-1
                ABSOUT(NOUT)=-XID(1,nd)
              ELSE IF(XID(1,nd).GT.1.d0) THEN
                NOUT=NOUT+1
                NIOUT(NOUT)=1
                ABSOUT(NOUT)=XID(1,nd)-1.d0
              ELSE
                LD(nd)=LN(L2)
                GO TO 8
              ENDIF
              IF(ABS(NXI(NIOUT(NOUT),1,L2)).NE.L2) GO TO 6
              LD(nd)=LN(L2)
              GO TO 7
    6         L2=ABS(NXI(NIOUT(NOUT),1,L2))
              ni=ABS(NIOUT(NOUT))
              IF(XID(ni,nd).LT.0.d0) XID(ni,nd)=DMAX1(XID(ni,nd)+1.d0,
     '          0.75D0)
              IF(XID(ni,nd).GT.1.d0) XID(ni,nd)=DMIN1(XID(ni,nd)-1.d0,
     '          0.25D0)
    7       CONTINUE
    8     CONTINUE
        ENDIF
C   9 CONTINUE
C     DO 10 nl=1,NLM
C       NDLT(nl)=0
C  10 CONTINUE
C     DO 11 nd=1,NDT
C       nl=ABS(LD(nd))
C       NDLT(nl)=NDLT(nl)+1
C       NDL=NDLT(nl)
C       NDDL(nl,NDL)=nd
C  11 CONTINUE
C     IF(DOP) THEN
C       DO 12 nd=1,NDT
C         IF(LD(nd).GT.0) THEN
C           WRITE(IO4,2001) nd,LD(nd),(XID(ni,nd),ni=1,2),SQ(nd)
C2001       FORMAT('     nd=',I4,4X,'LD(nd)=',I4,9X,
C    '             'XID(ni,nd)=',2(E12.6,4X),4X,'SQ(nd)=',E12.6)
C         ELSE IF(LD(nd).LT.0) THEN
C           WRITE(IO4,2002) nd,LD(nd),XID(1,nd),SQ(nd)
C2002       FORMAT('     nd=',I4,4X,'LD(nd)=',I4,9X,
C    '             'XID(1,nd)=',E12.6,4X,'SQ(nd)=',E12.6)
C         ENDIF
C  12   CONTINUE
C       DO 13 L=1,LN(0)
C         nl=ABS(LN(L))
C         WRITE(IO4,2003) L,nl,NDLT(nl)
C2003     FORMAT(' L=',I4,' nl=',I4,' NDLT(nl)=',I4)
C         IF(NDLT(nl).GT.0) WRITE(IO4,2004)
C    '      (NDDL(nl,NDL),NDL=1,NDLT(nl))
C2004     FORMAT(12X,'NDDL(nl,ndl)=',10I4,16(/24X,10I4))
C  13   CONTINUE
C     ENDIF

      CALL EXITS('NWXID')
      RETURN
 9999 CALL ERRORS('NWXID',ERROR)
      CALL EXITS('NWXID')
      RETURN 1
      END


