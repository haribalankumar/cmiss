      SUBROUTINE UPLINE(NBJ,NEL,NLL,NPL,NPNE,NSB,SE,ERROR,*)

C#### Subroutine: UPLINE
C###  Description: Sets the scale_factors the same for common lines

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEL(0:NELM,NLM),NLL(12,NEM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NSB(NKM,NNM,NBFM)
      REAL*8 SE(NSM,NBFM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nae1,nae2,nb,ne1,ne2,ne11,ne22,nl,nl1,nl2,noelem,
     '  ns,VTOE(3,12),VTOE2d(3,4),xi_line_1,xi_line_2
      REAL*8 SE2(18,2)
      LOGICAL PROGRESS

      CALL ENTERS('UPLINE',*9999)

        IF(NJT.EQ.3) THEN ! Only implemented for 3D
C Array of local nodes for local line numbers, and the xi direction
          DATA VTOE /1,2,1,
     '               3,4,1,
     '               5,6,1,
     '               7,8,1,
     '               1,3,2,
     '               5,7,2,
     '               2,4,2,
     '               6,8,2,
     '               1,5,3,
     '               2,6,3,
     '               3,7,3,
     '               4,8,3/
          DATA VTOE2d /1,2,1,
     '                 3,4,1,
     '                 1,3,2,
     '                 2,4,2/


C  Should initialise SE2 to zeroes.
          SE2(16,1)=0
          SE2(17,1)=0
          DO nl=1,NLT
            noelem=NEL(0,nl)
C Only accepting cubic hermite basis in all 3 direction
            IF(noelem.GT.1.AND.NPL(1,1,nl).EQ.4.AND.
     '        NPL(1,2,nl).EQ.4.AND.NPL(1,3,nl).EQ.4) THEN
C! CS 26.4.2005 I don't think this really does only let 
C! tri-Cubic Hermite though, adding support for bicubics
C  Due to above conditions (3D and cubic hermite) nb should
C  now be a 3cubic hermite basis with 8 lines per elem.
              nb=NBJ(1,NEL(1,nl)) ! all nj's should be the same
              DO ne11=1,noelem-1
                ne1=NEL(ne11,nl)
                CALL SE2INPUT(nb,ne1,1,nl,ns,NLL,NPL,NPNE,NSB,SE,SE2,
     '               VTOE,VTOE2d,ERROR,*9999)
                IF(SE2(15,1).NE.0) THEN ! if 0, no line was found
                  DO ne22=ne11+1,noelem
                    ne2=NEL(ne22,nl)
                    CALL SE2INPUT(nb,ne2,2,nl,ns,NLL,NPL,NPNE,NSB,SE,
     '                   SE2,VTOE,VTOE2d,ERROR,*9999)
                    IF(SE2(15,2).NE.0) THEN ! if 0, no line was found
C  Search for other common line between elements (regarding
C  common cross derivatives)
                      PROGRESS=.TRUE.
                      nae1=0
                      SE2(18,1)=0
                      SE2(18,2)=0
                      DO WHILE(nae1.LE.(NLE(nb)-1).AND.PROGRESS)
                        nae1=nae1+1
                        nl1=NLL(nae1,ne1)
                        nae2=0
                        DO WHILE(nae2.LE.(NLE(nb)-1).AND.PROGRESS)
                          nae2=nae2+1
                          nl2=NLL(nae2,ne2)
                          if(NIT(nb).EQ.3) THEN
                            xi_line_1 = VTOE(3,nae1)
                            xi_line_2 = VTOE(3,nae2)
                          ELSE
                            xi_line_1 = VTOE2d(3,nae1)
                            xi_line_2 = VTOE2d(3,nae2)
                          ENDIF
                          IF(SE2(15,1).NE.xi_line_1) THEN
                            IF(nl1.EQ.nl2) THEN
                              SE2(18,1)=xi_line_1
                              SE2(18,2)=xi_line_2
                              PROGRESS=.FALSE.
                            ENDIF
                          ENDIF
                        ENDDO !while nae2
                      ENDDO !while nae1
                      CALL SE2CHANGE(SE2,ERROR,*9999)
                      CALL SE2OUTPUT(nb,ne1,1,nl,ns,NLL,NSB,SE,
     '                     SE2,VTOE,VTOE2d,ERROR,*9999)
                      CALL SE2OUTPUT(nb,ne2,2,nl,ns,NLL,NSB,SE,
     '                     SE2,VTOE,VTOE2d,ERROR,*9999)
                    ENDIF !SE2(15,2)
                  ENDDO !ne22
                ENDIF !SE2(15,1)
              ENDDO !ne11
            ENDIF !noelem
          ENDDO !nl
        ENDIF

      CALL EXITS('UPLINE')
      RETURN
 9999 CALL ERRORS('UPLINE',ERROR)
      CALL EXITS('UPLINE')
      RETURN 1
      END



