      SUBROUTINE NFNXF(nb,ne,nf,NFF,NNF,NPF,NXF,NXI,XI2,XI3,COLLAPSE,
     '  ERROR,*)

C#### Subroutine: NFNXF
C###  Description:
C###    NFNXF sets up the NXF array from face information.  That is,
C###    the face connectivity is stored in NXF, taking into account
C###    any collapsed faces.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER nb,ne,nf,NFF(6,NEM),NNF(0:17,6,NBFM),NPF(9,NFM),
     '  NXF(-3:3,0:1,NFM),NXI(-NIM:NIM,0:NEIM,0:NEM),XI2,XI3
      LOGICAL COLLAPSE(4)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER i,nef0,nef,ne2,nf2,XID,XI_LIST(4)

      CALL ENTERS('NFNXF',*9999)

      XI_LIST(1)=-XI2
      XI_LIST(2)=XI2
      XI_LIST(3)=-XI3
      XI_LIST(4)=XI3
      DO i=1,4
        IF(.NOT.COLLAPSE(i))THEN !continue, edge not collapsed
          XID=XI_LIST(i)
          IF(NXI(XID,0,ne).NE.0)THEN !neighbouring element
            ne2=NXI(XID,1,ne)
            nef=NPF(8,nf)
            nf2=NFF(nef,ne2) !same local face in ne2
            IF(NPF(5,nf2).EQ.1)THEN !boundary face
              NXF(XID,0,nf)=1
              NXF(XID,1,nf)=nf2
            ENDIF
          ELSE !record neighbouring face in element ne
            DO nef=1,NFE(nb) !for each face
              nf2=NFF(nef,ne)
              IF(nf2.NE.0)THEN
                IF(nf2.NE.nf)THEN !a different face
                  IF(NPF(5,nf2).EQ.1)THEN !boundary face
                    IF(NNF(1,nef,nb).EQ.ABS(XID))THEN !XID normal face
                      IF(XID.LT.0)THEN
                        IF(nef.EQ.1.OR.nef.EQ.3.OR.nef.EQ.5)THEN !record
                          NXF(XID,0,nf)=1
                          NXF(XID,1,nf)=nf2
                        ENDIF
                      ELSE IF(XID.GT.0)THEN
                        IF(nef.EQ.2.OR.nef.EQ.4.OR.nef.EQ.6)THEN !record
                          NXF(XID,0,nf)=1
                          NXF(XID,1,nf)=nf2
                        ENDIF
                      ENDIF !XID
                    ENDIF !NNF
                  ENDIF !NPF
                ENDIF !nf2
              ELSE IF(nf2.EQ.0)THEN
                IF(NNF(1,nef,nb).EQ.ABS(XID))THEN !record opposite face
                  nef0=NPF(8,nf)
                  IF(XID.LT.0)THEN
                    IF(nef.EQ.1.OR.nef.EQ.3.OR.nef.EQ.5)THEN !record
                      NXF(XID,0,nf)=1
                      IF(nef0.EQ.1.OR.nef0.EQ.3.OR.nef0.EQ.5)THEN
                        NXF(XID,1,nf)=NFF(nef0+1,ne)
                      ELSE
                        NXF(XID,1,nf)=NFF(nef0-1,ne)
                      ENDIF
                    ENDIF
                  ELSE IF(XID.GT.0)THEN
                    IF(nef.EQ.2.OR.nef.EQ.4.OR.nef.EQ.6)THEN !record
                      NXF(XID,0,nf)=1
                      IF(nef0.EQ.1.OR.nef0.EQ.3.OR.nef0.EQ.5)THEN
                        NXF(XID,1,nf)=NFF(nef0+1,ne)
                      ELSE
                        NXF(XID,1,nf)=NFF(nef0-1,ne)
                      ENDIF !nef0
                    ENDIF !nef
                  ENDIF !XID
                ENDIF !NNF
              ENDIF !nf2
            ENDDO !nef
          ENDIF !NXI
        ENDIF !.NOT.COLLAPSE
      ENDDO !i

      CALL EXITS('NFNXF')
      RETURN
 9999 CALL ERRORS('NFNXF',ERROR)
      CALL EXITS('NFNXF')
      RETURN 1
      END


