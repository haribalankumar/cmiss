      SUBROUTINE IPAERO(NEELEM,NLL,NPL,nr,ERROR,*)

C#### Subroutine: IPAERO
C###  Description:
C###    IPAERO inputs aerofoil and wake parameters.

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NLL(12,NEM),NPL(5,0:3,NLM),nr
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,ne,nl,NOQUES,
     '  no_aero,noelem,no_entry,no_exit,no_wake
      LOGICAL FILEIP

      CALL ENTERS('IPAERO',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

! Aerofoil parameters
      FORMAT='(/$,'' Enter line#s on upper surface of aerofoil: '''
     '  //',20I4)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(0)=NL_AERO(0,1)
        DO no_aero=1,NL_AERO(0,1)
          IDATA(no_aero)=NL_AERO(no_aero,1)
        ENDDO
      ENDIF
      CDATA(1)='LINES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NL_AERO(0,1)=IDATA(0)
        DO no_aero=1,NL_AERO(0,1)
          NL_AERO(no_aero,1)=IDATA(no_aero)
        ENDDO
      ENDIF

      !find elements adjacent to upper aerofoil surface
      NE_AERO(0,1)=NL_AERO(0,1) !#elements=#lines
      DO no_aero=1,NL_AERO(0,1) !loop over upper aerofoil lines
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NLL(1,ne).EQ.NL_AERO(no_aero,1)) NE_AERO(no_aero,1)=ne
        ENDDO !noelem
      ENDDO !no_aero

      FORMAT='($,'' Enter line#s on lower surface of aerofoil: '''
     '  //',20I4)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(0)=NL_AERO(0,2)
        DO no_aero=1,NL_AERO(0,2)
          IDATA(no_aero)=NL_AERO(no_aero,2)
        ENDDO
      ENDIF
      CDATA(1)='LINES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NL_AERO(0,2)=IDATA(0)
        DO no_aero=1,NL_AERO(0,2)
          NL_AERO(no_aero,2)=IDATA(no_aero)
        ENDDO
      ENDIF

      !find elements adjacent to lower aerofoil surface
      NE_AERO(0,2)=NL_AERO(0,2) !#elements=#lines
      DO no_aero=1,NL_AERO(0,2) !loop over lower aerofoil lines
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NLL(2,ne).EQ.NL_AERO(no_aero,2)) NE_AERO(no_aero,2)=ne
        ENDDO !noelem
      ENDDO !no_aero

      NP_aero_LE =NPL(2,1,NL_AERO(1,1))            !leading  edge node
      NP_aero_TE1=NPL(3,1,NL_AERO(NL_AERO(0,1),1)) !trailing edge node 1
      NP_aero_TE2=NPL(3,1,NL_AERO(NL_AERO(0,2),2)) !trailing edge node 2

! Wake parameters
      FORMAT='($,'' Enter line#s on the wake: '',20I4)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(0)=NL_WAKE(0,1)
        DO no_wake=1,NL_WAKE(0,1)
          IDATA(no_wake)=NL_WAKE(no_wake,1)
        ENDDO
      ENDIF
      CDATA(1)='LINES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NL_WAKE(0,1)=IDATA(0)
        NL_WAKE(0,2)=NL_WAKE(0,1)
        DO no_wake=1,NL_WAKE(0,1)
          nl=IDATA(no_wake)
          NL_WAKE(no_wake,1)=nl
          NP_WAKE(no_wake,1)=NPL(2,1,nl) !is node# at Xi=0 end of line
        ENDDO
        NP_WAKE(NL_WAKE(0,1)+1,1)=NPL(3,1,nl) !is last node#
        NP_WAKE(0,1)=NL_WAKE(0,1)+1 !is total# nodes on upper wake
      ENDIF

      !find elements adjacent to upper & lower surfaces of wake
      NE_WAKE(0,1)=NL_WAKE(0,1) !#elements=#lines
      NE_WAKE(0,2)=NL_WAKE(0,2) !#elements=#lines
      DO no_wake=1,NL_WAKE(0,1) !loop over wake lines
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          IF(NLL(1,ne).EQ.NL_WAKE(no_wake,1)) THEN
            NE_WAKE(no_wake,1)=ne
          ELSE IF(NLL(2,ne).EQ.NL_WAKE(no_wake,1)) THEN
            NE_WAKE(no_wake,2)=ne
          ENDIF
        ENDDO !noelem
      ENDDO !no_wake

! Downstream outflow face parameters
      FORMAT='($,'' Enter line#s on the downstream outflow face: '''
     '  //',20I4)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(0)=NL_EXIT(0,1)
        DO no_exit=1,NL_EXIT(0,1)
          IDATA(no_exit)=NL_EXIT(no_exit,1)
        ENDDO
      ENDIF
      CDATA(1)='LINES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NL_EXIT(0,1)=IDATA(0)
        DO no_exit=1,NL_EXIT(0,1)
          NL_EXIT(no_exit,1)=IDATA(no_exit)
        ENDDO
      ENDIF

! Upstream inflow face nodes
      FORMAT='($,'' Enter node#s on the upstream inflow face: '',20I4)'
      IF(IOTYPE.EQ.3) THEN
        IDATA(0)=NP_ENTRY(0)
        DO no_entry=1,NP_ENTRY(0)
          IDATA(no_entry)=NP_ENTRY(no_entry)
        ENDDO
      ENDIF
      CDATA(1)='NODES' !for use with group input
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        NP_ENTRY(0)=IDATA(0)
        DO no_entry=1,NP_ENTRY(0)
          NP_ENTRY(no_entry)=IDATA(no_entry)
        ENDDO
      ENDIF

      CALL_AERO=.TRUE.

      CALL EXITS('IPAERO')
      RETURN
 9999 CALL ERRORS('IPAERO',ERROR)
      CALL EXITS('IPAERO')
      RETURN 1
      END


