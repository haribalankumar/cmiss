      SUBROUTINE SHCONT(ISCONO,ISCONT,ISEG,NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: SHCONT
C###  Description:
C###    SHCONT shows element contour segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'map000.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISCONO(NHM,NEM),ISCONT(NHM,NEM,NGRSEGM),ISEG(*),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,IL2(100),iw,N3CO,ne,nh,noil2,nocont,
     '  noelem,nolist,nr,NTIL2
      LOGICAL ABBREV,CBBREV

      CALL ENTERS('SHCONT',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show contours
C###  Description:
C###    Make the specified contour segments in the specified elements
C###    visible.
C###  Parameter:    <of NH_VARIABLE#[1]>
C###    Specify the dependent variable number to draw contours of.
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify which elements to show contours in.
C###  Parameter:    <at CONTOUR_INDICE#s[1]>
C###    Specify which contour segments to show


        OP_STRING(1)=STRING(1:IEND) //' <of NH_VARIABLE#[1]>'
        OP_STRING(2)=BLANK(1:15) //'<in (all/ELEMENT#s)[all]>'
        OP_STRING(3)=BLANK(1:15) //'<at CONTOUR_INDICE#s[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM show contour numbers
C###  Parameter:    <of NH_VARIABLE#[1]>
C###  Description:
C###    Make the specified contour numbers visible.

        OP_STRING(1)=STRING(1:IEND) //' numbers'
        OP_STRING(2)=BLANK(1:15) //'<of NH_VARIABLE#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHCONT',ERROR,*9999)
      ELSE
        iw=2*NJT-3+IMAP
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          IF(ABBREV(CO(noco+1),'NUMBERS',1)) THEN
            IF(ABBREV(CO(noco+2),'OF',1)) THEN
              nh=IFROMC(CO(noco+3))
            ELSE
              nh=1
            ENDIF
            DO nr=1,NRT
              DO noelem=1,NEELEM(0,nr)
                ne=NEELEM(noelem,nr)
                IF(ISEG(ISCONO(nh,ne)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISCONO(nh,ne),'VISIBLE',ERROR,
     '              *9999)
                ENDIF
              ENDDO
            ENDDO
          ELSE
            IF(ABBREV(CO(noco+1),'OF',1)) THEN
              nh=IFROMC(CO(noco+2))
            ELSE
              nh=1
            ENDIF
            IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,
     '          *9999)
            ELSE
              NELIST(0)=0
              DO nr=1,NRT
                DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
                  NELIST(noelem)=NEELEM(noelem,nr)
                ENDDO
                NELIST(0)=NELIST(0)+NEELEM(0,nr)
              ENDDO
            ENDIF
            IF(ABBREV(CO(noco+1),'AT',1)) THEN
              CALL PARSIL(CO(noco+2),100,NTIL2,IL2,ERROR,*9999)
            ELSE IF(ABBREV(CO(noco+3),'AT',1)) THEN
              CALL PARSIL(CO(noco+4),100,NTIL2,IL2,ERROR,*9999)
            ELSE IF(ABBREV(CO(noco+5),'AT',1)) THEN
              CALL PARSIL(CO(noco+6),100,NTIL2,IL2,ERROR,*9999)
            ELSE
              NTIL2=1
              IL2(1)=NTCONT
            ENDIF
C ***       loop over contour indices
            DO noil2=1,NTIL2
              nocont=IL2(noil2)
C ***         loop over elements
              DO nolist=1,NELIST(0)
                ne=NELIST(nolist)
                IF(ISEG(ISCONT(nh,ne,nocont)).EQ.1) THEN
                  CALL VISIB(iw,ISEG,ISCONT(nh,ne,nocont),'VISIBLE',
     '              ERROR,*9999)
                ELSE IF(ISEG(ISCONT(nh,ne,nocont)).EQ.0) THEN
                  WRITE(OP_STRING,
     '              '('' >>Contours not defined on '',I1)') iw
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHCONT')
      RETURN
 9999 CALL ERRORS('SHCONT',ERROR)
      CALL EXITS('SHCONT')
      RETURN 1
      END


