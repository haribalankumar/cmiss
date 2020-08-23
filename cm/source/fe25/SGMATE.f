      SUBROUTINE SGMATE(IBT,IDO,INP,ISEG,ISMATE,iw,NAN,NBJ,ne,nm,nr,
     '  AXES_LENGTHS,AXES_XI,CE,CQ,CSEG,TYPE,XE,XG,XQ,ERROR,*)

C#### Subroutine: SGMATE
C###  Description:
C###    SGMATE creates element material segment ISMATE(iw,ne).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'colo00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'graf00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'scal00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  ISEG(*),ISMATE,iw,NAN(NIM,NAM,NBFM),NBJ(NJM),ne,nm,nr
      REAL*8 AXES_LENGTHS(*),AXES_XI(3,*),CE,CQ(NMM,NQM),
     '   XE(NSM,NJM),XG(NJM,NUM),XQ(NJM,NQM)
      CHARACTER CSEG(*)*(*),ERROR*(*),TYPE*(*)
!     Local Variables
      INTEGER INDEX,INDEX1,INDEX2,INDEX3,INDEX_OLD,INDEX_POLYLINE,
     '  nb,nj,nq
      REAL*8 A_VECTOR(3),B_VECTOR(3),C_VECTOR(3),
     '  AXES(3,2),POINT(3),POINTS(3,2),
     '  PXI,SCALE,XI(3),xi_1,xi_2,xi_3,ZVAL

      CALL ENTERS('SGMATE',*9999)
      CALL OPEN_SEGMENT(ISMATE,ISEG,iw,'MATE',1,INDEX_OLD,
     '  ne,1,CSEG,ERROR,*9999)

      INDEX1=INDEX_POLYLINE(0,'SOLID','WIDTH1','RED')
      INDEX2=INDEX_POLYLINE(0,'SOLID','WIDTH1','BLUE')
      INDEX3=INDEX_POLYLINE(0,'SOLID','WIDTH1','GREEN')

      IF(TYPE(1:4).EQ.'AXES') THEN
        SCALE=0.1D0*DBLE(DIAG)
        IF(NJT.EQ.2) THEN
          DO xi_1=AXES_XI(1,1),AXES_XI(2,1),AXES_XI(3,1)
            WRITE(OP_STRING,'('' xi_1='',F5.3)') xi_1
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            XI(1)=xi_1
            DO xi_2=AXES_XI(1,2),AXES_XI(2,2),AXES_XI(3,2)
              WRITE(OP_STRING,'('' xi_2='',F5.3)') xi_2
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              XI(2)=xi_2
              DO nj=1,NJT
                nb=NBJ(nj)
                AXES(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '            nb,1,XI,XE(1,nj))
              ENDDO
              CALL XZ(ITYP10(1),AXES(1,1),AXES(1,1))
              AXES(1,2)=AXES(1,1)+SCALE*AXES_LENGTHS(1)
              AXES(2,2)=AXES(2,1)
              CALL XZ(ITYP10(1),AXES(1,2),AXES(1,2))
              CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
              AXES(1,2)=AXES(1,1)
              AXES(2,2)=AXES(2,1)+SCALE*AXES_LENGTHS(2)
              CALL XZ(ITYP10(1),AXES(1,2),AXES(1,2))
              CALL POLYLINE(INDEX,iw,2,AXES,ERROR,*9999)
            ENDDO !xi_2
          ENDDO !xi_1

        ELSE IF(NJT.EQ.3) THEN
          DO xi_1=AXES_XI(1,1),AXES_XI(2,1),AXES_XI(3,1)
            XI(1)=xi_1
            DO xi_2=AXES_XI(1,2),AXES_XI(2,2),AXES_XI(3,2)
              XI(2)=xi_2
              DO xi_3=AXES_XI(1,3),AXES_XI(2,3),AXES_XI(3,3)
                XI(3)=xi_3
                DO nj=1,NJT
                  nb=NBJ(nj)
                  AXES(nj,1)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '              nb,1,XI,XE(1,nj))
                ENDDO
                CALL XZ(ITYP10(1),AXES(1,1),AXES(1,1))
                CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ,nr,
     '            A_VECTOR,B_VECTOR,C_VECTOR,XE,XG,XI,.TRUE.,
     '            ERROR,*9999)
!  fibre axis
                DO nj=1,NJT
                  AXES(nj,2)=AXES(nj,1)+SCALE*AXES_LENGTHS(nj)
     '              *A_VECTOR(nj)
                ENDDO
                CALL POLYLINE(INDEX1,iw,2,AXES,ERROR,*9999)
!  sheet axis
                DO nj=1,NJT
                  AXES(nj,2)=AXES(nj,1)+SCALE*AXES_LENGTHS(nj)
     '              *B_VECTOR(nj)
                ENDDO
                CALL POLYLINE(INDEX2,iw,2,AXES,ERROR,*9999)
!  trans axis
                DO nj=1,NJT
                  AXES(nj,2)=AXES(nj,1)+SCALE*AXES_LENGTHS(nj)
     '              *C_VECTOR(nj)
                ENDDO
                CALL POLYLINE(INDEX3,iw,2,AXES,ERROR,*9999)
              ENDDO !xi_3
            ENDDO !xi_2
          ENDDO !xi_1
        ENDIF !njt

      ELSE IF(TYPE(1:5).EQ.'FIELD') THEN
        nb=NBJ(1)
        IF(NIT(nb).EQ.1) THEN
          DO nj=1,NJT !to find position of ends of line
            nb=NBJ(nj)
            POINTS(nj,1)=XE(1,nj)
            POINTS(nj,2)=XE(1+NKT(0,nb),nj)
          ENDDO !nj
          ZVAL=CE
          IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
          IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
          IF(COLOUR_WS) THEN
            INDEX=256-NINT((ZVAL-ZMINI)/ZDIFF*215.0D0)
          ELSE
            INDEX=8-NINT((ZVAL-ZMINI)/ZDIFF*3.0D0)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,
     '        '('' Value='',E12.3,'' INDEX='',I3)') ZVAL,INDEX
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          CALL POLYLINE(INDEX,iw,2,POINTS,ERROR,*9999)

        ELSE IF(NIT(nb).EQ.2) THEN
          DO nq=1,NQT
            DO nj=1,NJT
              POINT(nj)=XQ(nj,nq)
            ENDDO !nj
            ZVAL=CQ(nm,nq)
            IF(ZVAL.LT.ZMINI) ZVAL=ZMINI
            IF(ZVAL.GT.ZMAXI) ZVAL=ZMAXI
            IF(COLOUR_WS) THEN
C           Keep index between 45 and 200
              IF(DABS(ZDIFF).GT.1.D-6) THEN
                INDEX=200-NINT((ZVAL-ZMINI)/ZDIFF*155.0d0)
              ELSE
                INDEX=15
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ZVAL='',E12.3,'' Index='',I3)')
     '            ZVAL,INDEX
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF !dop
            ENDIF !COLOUR_WS
            CALL POLYMARKER(INDEX,iw,1,POINT,ERROR,*9999)
          ENDDO !nq
        ENDIF !NIT(nb)
      ENDIF !type

c     X(1)=XWC
c     X(2)=YWC
c     CHAR2=CFROMI(nm,'(I2)')
c     CALL STRING_TRIM(CHAR2,IBEG,IEND)
c     CALL TEXT(1,iw,CHAR2(IBEG:IEND),X,ERROR,*9999)

      CALL CLOSE_SEGMENT(ISMATE,iw,ERROR,*9999)

      CALL EXITS('SGMATE')
      RETURN
 9999 CALL ERRORS('SGMATE',ERROR)
      CALL EXITS('SGMATE')
      RETURN 1
      END


