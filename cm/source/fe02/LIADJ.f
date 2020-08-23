      INTEGER FUNCTION LIADJ(LL,nl,NPL,NVJL)

C#### Function: LIADJ
C###  Type: INTEGER
C###  Description:
C###    LIADJ searches for line adjacent to nl at the Xi=0 end (LL=1)
C###    or Xi=1 end (LL=2) of nl.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
!     Parameter List
      INTEGER LL,nl,NPL(5,0:3,NLM),NVJL(4,NJM,NLM)
!     Local Variables
      INTEGER COUNT,lowest_npp2,nll,n,nn,nnn,np,npp,npp2,NUMNN,NUMNNN

      LIADJ=0
      lowest_npp2=IMAX
      IF(NPL(1,1,nl).EQ.4.OR.NPL(1,1,nl).EQ.6.OR.
     '  NPL(1,1,nl).EQ.7) THEN
        NUMNN=2
      ELSE
        NUMNN=1+NPL(1,1,nl)
      ENDIF
      IF(LL.EQ.1) THEN
        nn=1
      ELSE
        nn=NUMNN
      ENDIF
      np=NPL(1+nn,1,nl)
      DO nll=1,NLT
        IF(nll.NE.nl.AND.NPL(1,0,nll).EQ.NPL(1,0,nl)) THEN
          IF(NPL(1,1,nll).EQ.4.OR.NPL(1,1,nll).EQ.6.OR.
     '      NPL(1,1,nll).EQ.7) THEN
            NUMNNN=2
          ELSE
            NUMNNN=1+NPL(1,1,nll)
          ENDIF
          IF(LL.EQ.1) THEN
            nnn=NUMNNN
          ELSE
            nnn=1
          ENDIF
          npp=NPL(1+nnn,1,nll)
          IF(npp.EQ.np) THEN
C           Check that the versions are equivalent (geom var 1 only)
            IF(NVJL(nnn,1,nll).EQ.NVJL(nn,1,nl)) THEN
              IF(NPL(1,1,nll).EQ.NPL(1,1,nl)) THEN
C KAT 2001-02-07: This seems to check for equivalent lines, but I can't
C             see why we would have found an equivalent line with one
C             node at both ends?
C              COUNT=0
C              DO nn=1,NUMNN
C                IF(NPL(1+nn,1,nll).EQ.NPL(1+nn,1,nl)) COUNT=COUNT+1
C              ENDDO
C              IF(COUNT.LT.NUMNN) THEN
C KAT         This seems to check for a coexistant line running the
C             other way.
                COUNT=0
                DO n=1,NUMNN
                  IF(NPL(1+n,1,nll).EQ.NPL(NUMNN-n+2,1,nl))
     '              COUNT=COUNT+1
                ENDDO
                IF(COUNT.LT.NUMNN) THEN
                  IF(JTYP2C.EQ.1) THEN
C CS 5/10/98 If there are multiple adjacent lines use the one with
C            the lowest node number at the other end
                    IF(LL.EQ.1) THEN
                      npp2=NPL(2,1,nll)
                    ELSE
                      npp2=NPL(1+NUMNNN,1,nll)
                    ENDIF
                    IF(npp2.LT.lowest_npp2) THEN
                      LIADJ=nll
                      lowest_npp2=npp2
                    ENDIF
                  ELSE
                    LIADJ=nll
                    GOTO 200
                  ENDIF
                ENDIF
C              ENDIF
              ELSE
                LIADJ=nll
                GOTO 200
              ENDIF
            ENDIF !versions
          ENDIF
        ENDIF
      ENDDO

 200  RETURN
      END


