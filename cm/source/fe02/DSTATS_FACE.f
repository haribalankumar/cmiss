      SUBROUTINE DSTATS_FACE(FD,LD,NDDL,NDLT,NFF,ERROR,*)

C#### Subroutine: DSTATS_FACE
C###  Description:
C###    DSTATS_FACE determines data point information based on gobal
C###    face no. for face fitting problems of volume meshes.
C**** Written by Kumar Mithraratne, Aug. 2002.

C**** NDDL(nf,ndf)
C**** NDLT(nf)
C**** from LD(nd),FD(nd),NFF(nfe,ne).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
!     Parameter List
      INTEGER FD(NDM),LD(NDM),NDDL(NEM,NDEM),NDLT(NEM),NFF(6,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,IEND,nd,ne,nf,nf_e
      CHARACTER CHAR5*5

      CALL ENTERS('DSTATS_FACE',*9999)

      IF(NEM.LE.NFT) THEN
        WRITE(CHAR5,'(I5)') NFT+1
        CALL STRING_TRIM(CHAR5,IBEG,IEND)
        ERROR=' Increase NEM in .IPPARA file to '//CHAR5(IBEG:IEND)
        GO TO 9999
      ENDIF

      DO nf=1,NFT
        NDLT(nf)=0
      ENDDO

      DO nd=1,NDT
        ne=LD(nd)   !element no.
        nf_e=FD(nd)  !local face no.
        !GDR 14/2/05 if nf_e or ne is zero then access of NFF fails 
        IF((nf_e.NE.0).AND.(ne.NE.0)) THEN
          nf=NFF(nf_e,ne) !global face no
        ELSE
          nf=0
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'('' nd='',I4,'' gl. face ='',I4)') nd,nf
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(nf.EQ.0) THEN
          WRITE(OP_STRING,
     '      '('' No face associated with data point'',I5)') nd
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSEIF (nf.GT.0) THEN
          IF(nf.LE.NFT) THEN
            NDLT(nf)=NDLT(nf)+1
            IF(NDLT(nf).LE.NDEM) THEN
              NDDL(nf,NDLT(nf))=nd
            ELSE
              ERROR=' Too many data points in face: Increase '//
     '          'NDEM in .IPPARA file'
              GO TO 9999
            ENDIF
          ENDIF
        ELSE
          WRITE(OP_STRING,
     '      '('' WARNING : Invalid face associated with nd'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !nd

      CALL EXITS('DSTATS_FACE')
      RETURN
 9999 CALL ERRORS('DSTATS_FACE',ERROR)
      CALL EXITS('DSTATS_FACE')
      RETURN 1
      END



