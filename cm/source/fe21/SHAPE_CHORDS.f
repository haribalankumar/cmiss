      SUBROUTINE SHAPE_CHORDS(ISEG,ISLINO,NLCHOR,CSEG,STRING,ERROR,*)

C#### Subroutine: SHAPE
C###  Description:
C###    SHAPE shapes chords.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISLINO(NWM),NLCHOR(0:10,NRM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,INSTAT,IPICK,ISEGM,N,N1CHOR,nl,nochor

      CALL ENTERS('SHAPE_CHORDS',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM shape chords
C###  Parameter:     <CHORD#[pick]>
C###  Description:
C###    Shape sail chords.

        OP_STRING(1)=STRING(1:IEND)//' <CHORD#[pick]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','SHAPE',ERROR,*9999)
      ELSE
        IF(NTCO.GT.noco) THEN
          nochor=IFROMC(CO(noco+1))
          CALL SHAPEC(ERROR,*9999)
        ELSE
          IF(ISEG(ISLINO(1)).EQ.1) THEN
            CALL VISIB(1,ISEG,ISLINO(1),'VISIBLE',ERROR,*9999)
          ENDIF
          CALL DETECT(1,ISEG,ISLINO(1),'DETECTABLE',ERROR,*9999)
          INSTAT=1
          DO WHILE(INSTAT.EQ.1)
            CALL ACWK(1,0,ERROR,*9999)
            CALL PICK(1,'REQUEST',INSTAT,ISEGM,IPICK,ERROR,*9999)
            CALL DAWK(1,0,ERROR,*9999)
            IF(INSTAT.EQ.1) THEN
              nl=IFROMC(CSEG(ISEGM)(53:57))
              DO nochor=1,NTCHOR
                DO N=1,NLCHOR(0,nochor)
                  IF(NLCHOR(N,nochor).eq.nl) THEN
                    N1CHOR=nochor
                    GO TO 201
                  ENDIF
                ENDDO
              ENDDO
 201          WRITE(OP_STRING,'('' Chord line segment '',I4,'
     '          //''' Chord '',I2)') nl,N1CHOR
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              CALL SHAPEC(ERROR,*9999)
            ENDIF
          ENDDO
          CALL DETECT(1,ISEG,ISLINO(1),'UNDETECTABLE',ERROR,*9999)
          IF(ISEG(ISLINO(1)).EQ.2) THEN
            CALL VISIB(1,ISEG,ISLINO(1),'INVISIBLE',ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF

      CALL EXITS('SHAPE_CHORDS')
      RETURN
 9999 CALL ERRORS('SHAPE_CHORDS',ERROR)
      CALL EXITS('SHAPE_CHORDS')
      RETURN 1
      END


