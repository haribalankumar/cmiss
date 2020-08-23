      SUBROUTINE IPIMPO(ERROR,*)

C#### Subroutine: IPIMPO
C###  Description:
C###    IPIMPO does input for import parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'emap00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'impo00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,nj,NOQUES
      LOGICAL FILEIP

      CALL ENTERS('IPIMPO',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='(/'' Specify the import type [1]: '''
     '  //'/''   (1) Signal'''
     '  //'/''   (2) Unused'''
     '  //'/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=DEFIMPO_TYPE
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,1,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DEFIMPO_TYPE=IDATA(1)

      IF(DEFIMPO_TYPE.EQ.1) THEN !Signal

        FORMAT='(/'' Specify the signal import type [1]: '''
     '    //'/''   (1) UNEMAP'''
     '    //'/''   (2) Unused'''
     '    //'/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=SIGIMPO_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,1,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) SIGIMPO_TYPE=IDATA(1)

        IF(SIGIMPO_TYPE.EQ.1) THEN !UNEMAP

C cpb 28/10/95 Only doing this for importing into rectangular cartesian
C for the moment.

          FORMAT='($,'' Enter the sock origin [0.0,0.0,0.0]: '',3D11.4)'
          IF(IOTYPE.EQ.3) THEN
            DO nj=1,NJT
              RDATA(nj)=EMAP_SOCKORIGIN(nj)
            ENDDO !nj
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              EMAP_SOCKORIGIN(nj)=RDATA(nj)
            ENDDO !nj
          ENDIF
          FORMAT='($,'' Enter the sock apex [0.0,0.0,-1.0]: '',3D11.4)'
          RDEFLT(1)=0.0d0
          RDEFLT(2)=0.0d0
          RDEFLT(3)=-1.0d0
          IF(IOTYPE.EQ.3) THEN
            DO nj=1,NJT
              RDATA(nj)=EMAP_SOCKAPEX(nj)
            ENDDO !nj
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              EMAP_SOCKAPEX(nj)=RDATA(nj)
            ENDDO !nj
          ENDIF
          FORMAT='($,'' Enter the sock theta origin [0.0,-1.0,0.0]: '','
     '      //'3D11.4)'
          RDEFLT(1)=0.0d0
          RDEFLT(2)=-1.0d0
          RDEFLT(3)=0.0d0
          IF(IOTYPE.EQ.3) THEN
            DO nj=1,NJT
              RDATA(nj)=EMAP_SOCKTHETAZERO(nj)
            ENDDO !nj
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RZERO,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO nj=1,NJT
              EMAP_SOCKTHETAZERO(nj)=RDATA(nj)
            ENDDO !nj
          ENDIF
          FORMAT='($,'' Enter the sock focus [1.0]: '',D11.4)'
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=EMAP_SOCKFOCUS(0)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) EMAP_SOCKFOCUS(0)=RDATA(1)

        ENDIF
      ENDIF

      CALL_IMPO=.TRUE.

      CALL EXITS('IPIMPO')
      RETURN
 9999 CALL ERRORS('IPIMPO',ERROR)
      CALL EXITS('IPIMPO')
      RETURN 1
      END


