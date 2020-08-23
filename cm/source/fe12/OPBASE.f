      SUBROUTINE OPBASE(IBT,IDO,INP,NAN,NFAM,NGAP,NKEF,NNF,NNL,nu,
     '  PG,XIG,ERROR,*)

C#### Subroutine: OPBASE
C###  Description:
C###    OPBASE outputs basis functions.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NFAM,NGAP(NIM,NBM),NKEF(0:4,16,6,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),nu
      REAL*8 PG(NSM,NUM,NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,nbb,nobase

      CALL ENTERS('OPBASE',*9999)
C CPB 15/7/92 ADDED BASIS FAMILY CODE
      IF(JTYP8.NE.2) THEN        !normal list basis/list basis full
        FORMAT='(''  The number of basis function families = '',I2)'
        WRITE(OP_STRING,FORMAT) NBFT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(JTYP8.EQ.2) THEN   !list basis family
        IF(NFAM.EQ.0) THEN
          FORMAT='(''  The total number of basis functions = '',I2)'
          WRITE(OP_STRING,FORMAT) NBT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          FORMAT='(''  The number of basis functions in family '',I2,'
     '      //''' = '',I2)'
          WRITE(OP_STRING,FORMAT) NFAM,NBASEF(NFAM,0)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(JTYP8.NE.2) THEN !list out family basis functions
        DO nb=1,NBFT  !family basis number
          nbb=NBASEF(nb,1) !global basis number
          CALL OPBASE1(IBT,IDO,INP,NAN,nbb,NFAM,NGAP,NKEF,NNF,NNL,nu,
     '      PG,XIG,ERROR,*9999)
        ENDDO
      ELSE !list basis family
        IF(NFAM.EQ.0) THEN
          DO nb=1,NBFT
            DO nobase=1,NBASEF(nb,0)
              nbb=NBASEF(nb,nobase) !global basis function
              CALL OPBASE1(IBT,IDO,INP,NAN,nbb,NFAM,NGAP,NKEF,NNF,NNL,
     '          nu,PG,XIG,ERROR,*9999)
            ENDDO !nobase
          ENDDO
        ELSE
          nb=NFBASE(1,NFAM) !family basis number
          DO nobase=1,NBASEF(nb,0) !number in family
            nbb=NBASEF(NFAM,nobase) !global basis number
            CALL OPBASE1(IBT,IDO,INP,NAN,nbb,NFAM,NGAP,NKEF,NNF,NNL,nu,
     '        PG,XIG,ERROR,*9999)
          ENDDO !nobase
        ENDIF !NFAM
      ENDIF !jtyp8

      CALL EXITS('OPBASE')
      RETURN
 9999 CALL ERRORS('OPBASE',ERROR)
      CALL EXITS('OPBASE')
      RETURN 1
      END


