      SUBROUTINE PSCOORD2(NBJ,NEELEM,NELEM,NKJE,NPF,
     '  NPNE,NTELEM,NVJE,SE,XA,XD,XE,XI,XP,ERROR,*)

C#### Subroutine: PSCOORD2
C###  Description:
C###    PSCOORD2 finds element Xi coords of a point on Hammer
C###    projection map given by prolate coords XD(nj),nj=1..3.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NELEM(10),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NTELEM,NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XD(3),XE(NSM,NJM),XI(3),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,NITB,nj,noelem,nn,nr,ns1,ns2,ns3,ns4
      REAL*8 XI1,XI2,XXE(4,2)

      CALL ENTERS('PSCOORD2',*9999)
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'('' Prolate coordinates:'',3E11.4)')
     '    (XD(nj),nj=1,3)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF
      DO nr=1,NRT
        NTELEM=0
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '      nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          nb=NBJ(2,ne)
          NITB=NIT(nb)
          IF(NITB.EQ.2) THEN
            ns1=1
            ns2=2
            ns3=3
            ns4=4
          ELSE IF(NITB.EQ.3) THEN
            ns1=5
            ns2=6
            ns3=7
            ns4=8
          ENDIF
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' Element '',I4)') ne
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' XE(1,2&3)='',2E11.3)')
     '        (XE(ns1,nj),nj=2,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' XE(2,2&3)='',2E11.3)')
     '        (XE(ns2,nj),nj=2,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' XE(3,2&3)='',2E11.3)')
     '        (XE(ns3,nj),nj=2,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' XE(4,2&3)='',2E11.3)')
     '        (XE(ns4,nj),nj=2,3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
C         XI1=(XD(3)-XE(ns1,3))/(XE(ns2,3)-XE(ns1,3)) !assumes rect. element
C         XI2=(XD(2)-XE(ns1,2))/(XE(ns3,2)-XE(ns1,2))
          DO nn=1,4
            XXE(nn,1)=XE(nn,3)
            XXE(nn,2)=XE(nn,2)
          ENDDO
          XI1=0.5D0
          XI2=0.5D0
          CALL XYCOORD(XD(3),XD(2),XXE,XI1,XI2,ERROR,*9999)
          IF((XI1.GT.0.0D0.AND.XI1.LE.1.0D0).AND.
     '       (XI2.GT.0.0D0.AND.XI2.LE.1.0D0)) THEN
            NTELEM=NTELEM+1
            NELEM(NTELEM)=ne
            XI(1)=XI1
            XI(2)=XI2
            IF(NITB.EQ.2) GO TO 860
          ENDIF
        ENDDO
 860    CONTINUE
        WRITE(OP_STRING,'('' Xi coordinates:  '',2E11.4)')XI(1),XI(2)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Element numbers '',10I3)')
     '    (NELEM(noelem),noelem=1,NTELEM)
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDDO

      CALL EXITS('PSCOORD2')
      RETURN
 9999 CALL ERRORS('PSCOORD2',ERROR)
      CALL EXITS('PSCOORD2')
      RETURN 1
      END


