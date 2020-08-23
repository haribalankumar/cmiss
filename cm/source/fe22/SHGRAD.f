      SUBROUTINE SHGRAD(ISEG,ISGRAD,NEELEM,NELIST,STRING,ERROR,*)

C#### Subroutine: SHGRAD
C###  Description:
C###    SHGRAD shows element gradient segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ntsg00.cmn'
!     Parameter List
      INTEGER ISEG(*),ISGRAD(NEM,NGRSEGM),NEELEM(0:NE_R_M,0:NRM),
     '  NELIST(0:NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IL2(100),iw,N3CO,ne,noelem,nograd,noil2,nolist,
     '  nr,NTIL2
      LOGICAL CBBREV

      CALL ENTERS('SHGRAD',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM show gradient
C###  Parameter:    <in (all/ELEMENT#s)[all]>
C###    Specify which elements to show gradient in.
C###  Parameter:    <at GRADIENT_INDICES[1]>
C###    Specify the workstation (GX window) to show the
C###    gradient on.
C###  Description:
C###    Make the specified gradient segments in the specified elements
C###    visible.

        OP_STRING(1)=STRING(1:IEND) //' <in (all/ELEMENT#s)[all]>'
        OP_STRING(2)=BLANK(1:15) //'<at GRADIENT_INDICES[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','SHGRAD',ERROR,*9999)
      ELSE
        iw=2*NJT-3
        IF(IWKS(iw).GT.0) THEN
          CALL ACWK(iw,1,ERROR,*9999)
          IF(CBBREV(CO,'IN',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),NEM,NELIST(0),NELIST(1),ERROR,*9999)
          ELSE
            NELIST(0)=0
            DO nr=1,NRT
              DO noelem=NELIST(0)+1,NELIST(0)+NEELEM(0,nr)
                NELIST(noelem)=NEELEM(noelem,nr)
              ENDDO
              NELIST(0)=NELIST(0)+NEELEM(0,nr)
            ENDDO
          ENDIF
          IF(CBBREV(CO,'AT',1,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),100,NTIL2,IL2,ERROR,*9999)
          ELSE
            NTIL2=1
            IL2(1)=NTGRAD
          ENDIF
C ***     loop over gradient indices
          DO noil2=1,NTIL2
            nograd=IL2(noil2)
C ***       loop over elements
            DO nolist=1,NELIST(0)
              ne=NELIST(nolist)
              IF(ISEG(ISGRAD(ne,nograd)).EQ.1) THEN
                CALL VISIB(iw,ISEG,ISGRAD(ne,nograd),'VISIBLE',ERROR,
     '            *9999)
              ELSE IF(ISEG(ISGRAD(ne,nograd)).EQ.0) THEN
                WRITE(OP_STRING,'('' >>Gradient in element '',I4,'
     '            //''' is not defined at '',I3)') ne,nograd
                CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO
          ENDDO
          CALL DAWK(iw,1,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('SHGRAD')
      RETURN
 9999 CALL ERRORS('SHGRAD',ERROR)
      CALL EXITS('SHGRAD')
      RETURN 1
      END


