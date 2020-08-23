      SUBROUTINE SGCONT(INDEX,IBT,IDO,INP,ISCONO,ISCONT,ISEG,iw,MXI,
     '  NAN,NBH,NBJ,ne,nh,NHE,nr,NTCOVA,nx,
     '  CONTYP,COVA,CSEG,LABEL_CONTOUR,
     '  PG,STATIC,TOL,XE,XG,XICONT,XISTEP,ZE,ZG,ERROR,*)

C#### Subroutine: SGCONT
C###  Description:
C###    SGCONT creates element contour segments ISCONT abd ISCONO.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'time00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INDEX,
     '  INP(NNM,NIM,NBFM),ISCONO,ISCONT,ISEG(*),iw,MXI(2),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM),NBJ(NJM),ne,nh,NHE,
     '  nr,NTCOVA,nx
      REAL*8 COVA(NEM,50),PG(NSM,NUM,NGM,NBM),TOL,
     '  XE(NSM,NJM),XG(NJM,NUM),XICONT,XISTEP,ZE(NSM,NHM),ZG(NHM,NUM)
      CHARACTER CONTYP*(*),CSEG(*)*(*),ERROR*(*)
      LOGICAL LABEL_CONTOUR,STATIC
!     Local Variables
      INTEGER IBEG,ICOORD,IEND,INDEX_OLD,IZVAL,na,nb,NBFF,NBG,
     '  nc,NICONT,nj,nk,nn,nocova,NROOTS,NSF,NSG
      REAL*8 PF1,ROOTS(2,12),TEMP,TEMP2,TWOTOL,XLM(3,50),ZVAL
      CHARACTER CHAR*20,TYPE*8
      LOGICAL BICUBIC,FAST_CONTOUR,LCON(1000)

      CALL ENTERS('SGCONT',*9999)

C GBS 10-NOV-1994
      nc=1   !Temporary

      CALL OPEN_SEGMENT(ISCONT,ISEG,iw,'CONT',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      IF(iw.EQ.4) THEN
        PROJEC=MAP_PROJEC
        IF(PROJEC(1:2).EQ.'XI') THEN
          MXI1=MXI(1) !bottom left coords
          MXI2=MXI(2) !..for Xi map projection
        ENDIF
      ENDIF

      nb=NBH(nh,nc)
      IF(CONTYP(1:5).EQ.'FIELD'.OR.CONTYP(1:9).EQ.'DEPENDENT') THEN
        IF(NIT(nb).EQ.2.AND.NKT(0,nb).EQ.1) THEN
          FAST_CONTOUR=.TRUE.
          TYPE='BILINEAR'
        ELSE IF(NIT(nb).EQ.2.AND.NKT(0,nb).EQ.4) THEN
          BICUBIC=.TRUE.
          DO nn=1,NNT(nb)
            IF(NKT(nn,nb).LT.NKT(0,nb))THEN
              BICUBIC=.FALSE.
            ENDIF
          ENDDO
          IF(BICUBIC)THEN
            FAST_CONTOUR=.TRUE.
            TYPE='BICUBIC'
          ELSE
            FAST_CONTOUR=.FALSE.
          ENDIF
        ELSE
          FAST_CONTOUR=.FALSE.
        ENDIF
      ELSE
        FAST_CONTOUR=.FALSE.
      ENDIF

      IF(.NOT.STATIC) THEN !do time interpolation here
        DO nj=1,NJT
          NBFF=NBH(nj,nc) !time basis
          NBG=NBJ(nj)  !corresponding spatial basis
          NSG=0
          DO nn=1,NNT(NBG)
            DO nk=1,NKT(0,NBG)
              NSG=NSG+1
              TEMP=0.0D0
              TEMP2=0.0D0
              DO na=1,IBT(2,NIT(NBFF),NBFF) !#terms in series
                NSF=(nn-1)*NKT(0,NBFF)+(na-1)*NKT(0,NBG)+nk
                TEMP=TEMP+PF1(na,1,TIME)*ZE(NSF,nj)
                TEMP2=TEMP2+PF1(na,1,TIME_REF)*ZE(NSF,nj)
              ENDDO
              ZE(NSG,nj)=TEMP+XE(NSG,nj)
              XE(NSG,nj)=TEMP2+XE(NSG,nj)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      NICONT=3

      DO nocova=1,NTCOVA
        ZVAL=COVA(ne,nocova)
        IF(DOP) THEN
          WRITE(OP_STRING,'('' ZVAL= '',E11.3)') ZVAL
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        TWOTOL=2.0D0*TOL

        IF(STATIC) THEN
          IF(FAST_CONTOUR) THEN
            IF(TYPE(1:8).EQ.'BILINEAR') THEN
              CALL CROOTS1(100,nh,NROOTS,ROOTS,TWOTOL,ZE,ZVAL,
     '          ERROR,*9999)
            ELSE IF(TYPE(1:7).EQ.'BICUBIC') THEN
              CALL CROOTS2(100,nh,NROOTS,ROOTS,TWOTOL,ZE,ZVAL,
     '          ERROR,*9999)
            ENDIF
          ELSE
            CALL CROOTS(CONTYP,IBT,IDO,INP,100,NAN,NBH,NBJ,
     '        nh,NHE,nj,ICOORD,NICONT,nr,NROOTS,nx,
     '        ROOTS,TWOTOL,XICONT,PG,XE,XG,ZE,ZG,ZVAL,ERROR,*9999)
          ENDIF
        ELSE IF(.NOT.STATIC) THEN !time contour with NBJ in NBH slot
          CALL CROOTS(CONTYP,IBT,IDO,INP,100,NAN,NBJ,NBJ,
     '      nh,NHE,nj,ICOORD,NICONT,nr,NROOTS,nx,
     '      ROOTS,TWOTOL,XICONT,PG,XE,XG,ZE,ZG,ZVAL,ERROR,*9999)
        ENDIF
        IF(NROOTS.GT.0) THEN
          LCON(nocova)=.TRUE.
        ELSE
          LCON(nocova)=.FALSE.
        ENDIF

        IF(STATIC) THEN
          IF(FAST_CONTOUR) THEN
            IF(TYPE(1:8).EQ.'BILINEAR') THEN
              CALL CONTOR1(INDEX,IBT,IDO,INP,
     '          100,iw,NBJ,nh,nj,NROOTS,
     '          ROOTS,TOL,XE,XICONT,XISTEP,
     '          XLM(1,nocova),ZE,ZVAL,ERROR,*9999)
            ELSE IF(TYPE(1:7).EQ.'BICUBIC') THEN
              CALL CONTOR2(INDEX,IBT,IDO,INP,
     '          100,iw,NBJ,nh,nj,NROOTS,
     '          ROOTS,TOL,XE,XICONT,XISTEP,
     '          XLM(1,nocova),ZE,ZVAL,ERROR,*9999)
            ENDIF
          ELSE
            CALL CONTOR(INDEX,CONTYP,IBT,IDO,INP,
     '        ICOORD,100,iw,NICONT,NAN,NBH,NBJ,nh,NHE,nj,
     '        nr,NROOTS,nx,
     '        ROOTS,TOL,PG,XE,XG,XICONT,XISTEP,XLM(1,nocova),
     '        ZE,ZG,ZVAL,ERROR,*9999)
          ENDIF
        ELSE IF(.NOT.STATIC) THEN !time contour called with NBJ in NBH slot
          CALL CONTOR(INDEX,CONTYP,IBT,IDO,INP,
     '      ICOORD,100,iw,NICONT,NAN,NBJ,NBJ,nh,NHE,nj,nr,NROOTS,nx,
     '      ROOTS,TOL,PG,XE,XG,XICONT,XISTEP,XLM(1,nocova),
     '      ZE,ZG,ZVAL,ERROR,*9999)
        ENDIF
      ENDDO

      CALL CLOSE_SEGMENT(ISCONT,iw,ERROR,*9999)

      CALL OPEN_SEGMENT(ISCONO,ISEG,iw,'CONO',INDEX,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      IF(LABEL_CONTOUR) THEN
        DO nocova=1,NTCOVA
          IF(LCON(nocova)) THEN
            ZVAL=COVA(ne,nocova)
            IZVAL=INT(ZVAL)
            WRITE(CHAR,'(I4)') IZVAL
            CALL STRING_TRIM(CHAR,IBEG,IEND)
            CALL TEXT(1,iw,CHAR(IBEG:IEND),XLM(1,nocova),ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL CLOSE_SEGMENT(ISCONO,iw,ERROR,*9999)

      CALL EXITS('SGCONT')
      RETURN
 9999 CALL ERRORS('SGCONT',ERROR)
      CALL EXITS('SGCONT')
      RETURN 1

      END


