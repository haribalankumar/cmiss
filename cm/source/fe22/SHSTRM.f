      SUBROUTINE SHSTRM(ISEG,ISSTRM,NEELEM,STRING,ERROR,*)

C#### Subroutine: SHSTRM
C###  Description:
C###    SHSTRM shows streamline segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISSTRM(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,iw,ne,noelem,nr

      CALL ENTERS('SHSTRM',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show streamlines
C###  Description:
C###    Make the streamline segments visible.

        OP_STRING(1)=STRING(1:IEND)
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHSTRM',ERROR,*9999)
      ELSE
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          DO nr=1,NRT
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              IF(ISEG(ISSTRM(ne,1)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISSTRM(ne,1),'VISIBLE',ERROR,*9999)
              ELSE IF(ISEG(ISSTRM(ne,1)).EQ.0) THEN
                WRITE(OP_STRING,
     '            '('' >>Streamline is not defined in element '',I1)')
     '            ne
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHSTRM')
      RETURN
 9999 CALL ERRORS('SHSTRM',ERROR)
      CALL EXITS('SHSTRM')
      RETURN 1
      END


