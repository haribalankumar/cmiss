      SUBROUTINE EXT_FACEFIT_LIST(NBJ,NELIST,NFF,NFLIST,NPF,
     '  EXTERNAL,ERROR,*)

C#### Subroutine: EXT_FACEFIT_LIST
C###  Description:
C###    EXT_FACEFIT_LIST determines the external faces in a volume
C###    mesh in face fitting of volume elements. All external faces
C###    are stored in NFLIST.
C**** Written by Kumar Mithraratne, Aug. 2002.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),NFF(6,NEM),NFLIST(0:NFM),
     '  NPF(9,NFM)
      CHARACTER ERROR*(*)
      LOGICAL EXTERNAL
!     Local Variables
      INTEGER nb,ne,nf_elem(0:1),nf_mis(0:2),no_ne,nf,no_nf


      CALL ENTERS('EXT_FACEFIT_LIST',*9999)

      DO nf=0,NFM
        NFLIST(nf)=0
      ENDDO

!     Creat external face list
      IF (EXTERNAL) THEN
        DO no_ne=1,NELIST(0)
          ne=NELIST(no_ne)
          nb=NBJ(1,ne)

          nf_elem(0)=0
          nf_mis(0)=0
          nf_mis(1)=100
          nf_mis(2)=100

          DO no_nf=1,NFE(nb)
            IF(NFF(no_nf,ne).NE.0) THEN
              nf_elem(0)=nf_elem(0)+1 ! number of faces in the element
            ELSEIF(NFF(no_nf,ne).EQ.0) THEN
              nf_mis(0)=nf_mis(0)+1 ! no. faces missing
              nf_mis(nf_mis(0))=no_nf    ! missing face nos. (local)
            ENDIF
          ENDDO

          DO no_nf=1,nf_elem(0)
            IF(no_nf.GE.nf_mis(1)) THEN
              nf_elem(1)=no_nf+1
            ELSEIF(no_nf.GE.nf_mis(2)) THEN
              nf_elem(1)=no_nf+2
            ELSE
              nf_elem(1)=no_nf
            ENDIF

C MHT 23-11-06: when faces are collapsed NFF=0 (i.e. no global face
C number) so NPF has subscript out of range if the following condition
C is tested for a collapsed face.
            IF(NFF(nf_elem(1),ne).NE.0)THEN
              IF(NPF(5,NFF(nf_elem(1),ne)).EQ.1) THEN !external face
                NFLIST(0)=NFLIST(0)+1 !NFLIST(0) initialised in DEXI
                NFLIST(NFLIST(0))=NFF(nf_elem(1),ne)
              ENDIF
            ENDIF
          ENDDO !no_nf
        ENDDO !no_ne

      ENDIF !EXTERNAL


      CALL EXITS('EXT_FACEFIT_LIST')
      RETURN
 9999 CALL ERRORS('EXT_FACEFIT_LIST',ERROR)
      CALL EXITS('EXT_FACEFIT_LIST')
      RETURN 1
      END



