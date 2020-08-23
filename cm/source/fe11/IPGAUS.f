      SUBROUTINE IPGAUS(NBJ,nr,NT_YG,YG,ERROR,*)

C#### Subroutine: IPGAUS
C###  Description:
C###    IPGAUS inputs Gauss point parameters YG(niyg,ng,ne) for region NR.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),nr,NT_YG
      REAL*8 YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG1,IBEG2,ICHAR,IEND1,IEND2,INFO,nb,ne,ng,NO_YG,
     '  NOQUES,yg_loop
      CHARACTER CHAR1*5,CHAR2*3
      LOGICAL FILEIP

      CALL ENTERS('IPGAUS',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='($,'' Specify number of variables/Gauss point '//
     '  ' [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=NT_YG
C JWF 15/08/04 Now allows a maximum of 3 variables per Gauss point rather than
C the 2 previously allowed.
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,3,
     &  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)  
      IF(IOTYPE.NE.3) NT_YG=IDATA(1)

C JWF 15/08/04 This code only asks for the variable no. if NT_YG is 1.
C Modified so will loop over this question and ask for the variable no.
C based on the value of NT_YG assigned above. 
      DO yg_loop=1,NT_YG
C      IF(NT_YG.EQ.1) THEN
        FORMAT='($,'' Enter the variable number [1]: '',I1)'
        IF(IOTYPE.EQ.3) THEN
          NO_YG=yg_loop
          IDATA(1)=NO_YG
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &    1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     &    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NO_YG=IDATA(1)
C      ENDIF

        DO ne=1,NET(nr)
          WRITE(CHAR1,'(I5)') ne
          CALL STRING_TRIM(CHAR1,IBEG1,IEND1)
          nb=NBJ(NO_YG,ne)
          DO ng=1,NGT(nb)
            WRITE(CHAR2,'(I3)') ng
            CALL STRING_TRIM(CHAR2,IBEG2,IEND2)
            RDEFLT(1)=0.0D0
            FORMAT='($,'' YG at ne='//CHAR1(1:IEND1)
     '        //' ng='//CHAR2(1:IEND2)//' is [0]: '',E12.5)'
            IF(IOTYPE.EQ.3) THEN
              RDATA(1)=YG(NO_YG,ng,ne)
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) THEN
              YG(NO_YG,ng,ne)=RDATA(1)
            ENDIF
          ENDDO !ng
        ENDDO !ne
      ENDDO !yg_loop
      
      CALL EXITS('IPGAUS')
      RETURN
 9999 CALL ERRORS('IPGAUS',ERROR)
      CALL EXITS('IPGAUS')
      RETURN 1
      END


