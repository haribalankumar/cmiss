      SUBROUTINE IPNONL(nr,nx,SEND_GLOBAL,ERROR,*)

C#### Subroutine: IPNONL
C###  Description:
C###    IPNONL inputs nonlinear analysis parameters.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'fsklib.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'host00.cmn'
      INCLUDE 'host00.inc'
      INCLUDE 'inout00.cmn'
      INCLUDE 'iwrit00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER nr,nx
      CHARACTER ERROR*(*)
      LOGICAL SEND_GLOBAL
!     Local Variables
      INTEGER IBEG,IBEG1,ICHAR,IEND,IEND1,INFO,nhost,NOQUES
      CHARACTER CONNIDSTR*100,PORTSTR*100
      LOGICAL FILEIP,OLD_HOST(N_HOSTS_MX) !N_HOSTS_MX in host00.cmn

      CALL ENTERS('IPNONL',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(.NOT.IS_COUPLED(nx).OR.nr.EQ.COUP_NRLIST(1,nx)) THEN
        FORMAT='('' Specify type of equilibrium iteration [1]:'''//
     '    '/''   (1) Newton-Raphson  (full update at every iter.)'''//
     '    '/''   (2) Modified Newton (with adaptive full update)'''//
     '    '/''   (3) BFGS inverse    (with adaptive rank 2 update)'''//
     '    '/''   (4) Sequential quadratic programming (E04UPF) '''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=ITYP9(nr,nx)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) ITYP9(nr,nx)=IDATA(1)

        IF(ITYP6(nr,nx).EQ.1.AND.ITYP9(nr,nx).EQ.3) THEN !multigrid
          CALL ASSERT(NCM.GE.4,'>>Need NCM=4 ',ERROR,*9999)
        ENDIF

        IF(ITYP9(nr,nx).LE.3) THEN
          FORMAT='('' Specify line search option [1]:'''//
     '      '/''   (1) No search '''//
     '      '/''   (2) Linear search '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP10
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP10=IDATA(1)

C XSL NEWS 18Aug2010 Add option to enter com file name for line search with contact mech
C so that projection codes can be read in from com file
          IF(KTYP10.EQ.2.AND.KTYP5G(nr).EQ.1) THEN !line search and contact
            FORMAT='($,'' Enter the command file name [SOLVE]: '',A)'
            CDEFLT(1)='SOLVE'
            IF(IOTYPE.EQ.3) CDATA(1)=COM_FILE
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) COM_FILE=CDATA(1)(1:255)
          ENDIF
C XSL NEWE


C *** DPN 09 May 2001 - Need to be able to perform a different test for
C *** convergence when using the couple cell/tissue active contraction,
C *** i.e. Want to use the ratio of unconstrained residuals to the
C *** maximum Gauss point active tension.
          FORMAT='('' Specify the convergence criteria [1]:'''//
     &     '/''   (1) Sum of differentiated ratios of unconstrained '//
     '      'to constrained residuals '''//
     '      '/''   (2) Ratio of unconstrained residuals to maximum '//
     '      'Gauss point value '''//
C 21/07/08 XSL Add new options to convergence check
     &     '/''   (3) Ratio of sum of unconstrained residuals to '//
     &     'sum of constrained residuals '''//
C XSL energy norm
     &     '/''   (4) Normalised scalar product of residuals '//
     &     'and solution increments '''//
     '      '/$,''    '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP007
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,4,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP007=IDATA(1)

          IF(KTYP007.EQ.2) THEN
            FORMAT='($,'' Specify the Gauss point field [1]: '',I3)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP007a
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     &        NIYGM,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) KTYP007a=IDATA(1)
          ENDIF

C news MPN 13-Sep-94: Parallel element stiffness matrix computations
          FORMAT='($,'' Do you want to use parallel element stiffness '
     '      //'computations [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(KTYP1A.EQ.2) THEN !KTYP1A=2  -> use parallel computations
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              KTYP1A=2 !use parallel elem stiff computations
            ELSE
              KTYP1A=1 !do not use parallel elem stiff computations
            ENDIF
          ENDIF

          IF(KTYP1A.EQ.2) THEN !use parallel elem stiff computations
            FORMAT='($,'' Enter the number of host machines [2]: '',I2)'
            IDEFLT(1)=2
            IF(IOTYPE.EQ.3) IDATA(1)=NUMHOSTS
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,
     '        MIN(N_HOSTS_MX,IOCMX),
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) NUMHOSTS=IDATA(1)

C           Enter remote host machine names
            CDEFLT(1)='esv25'
            CDEFLT(2)='esv26'
            CDEFLT(3)='esv27'
            CDEFLT(4)='esv28'
            CDEFLT(5)='esv29'
            CDEFLT(6)='esv30'
            CDEFLT(7)='esv31'
            CDEFLT(8)='esv32'
            FORMAT='($,'' Enter host machine names [esv25...]: '',A)'
            IF(IOTYPE.EQ.3) THEN
              DO nhost=1,NUMHOSTS
                CALL STRING_TRIM(HOSTS(nhost),IBEG,IEND)
                CDATA(nhost)=HOSTS(nhost)(IBEG:IEND)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NUMHOSTS,
     '        ADATA,ADEFLT,CDATA,CDEFLT,MXHOSTCH,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nhost=1,NUMHOSTS
                CALL STRING_TRIM(CDATA(nhost),IBEG,IEND)
                CALL STRING_TRIM(HOSTS(nhost),IBEG1,IEND1)
C               Check for previously defined hosts
                IF(IBEG1.LE.IEND1.AND. !HOSTS(nhost) is defined
     '            HOSTS(nhost)(IBEG1:IEND1).EQ.CDATA(nhost)(IBEG:IEND))
     '            THEN
                  OLD_HOST(nhost)=.TRUE.
                ELSE
                  HOSTS(nhost)=CDATA(nhost)(IBEG:IEND)
                  OLD_HOST(nhost)=.FALSE.
C                 If the socket is already  open, close it.
                  IF(SOCKET_OPEN(nhost)) THEN
C                   Signal to slave processes to stop
                    IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,
     '                ICONNID(nhost)).EQ.-1) GOTO 9999
C                    IRET=FSKCLOSE(ICONNID(nhost)) !close socket
                    SOCKET_OPEN(nhost)=.FALSE.
                  ENDIF
                ENDIF
              ENDDO
            ENDIF

C           Enter connection ID's for remote host machines
            DO nhost=1,NUMHOSTS
              IDEFLT(nhost)=nhost
            ENDDO
            IF(NUMHOSTS.GT.1) THEN
              WRITE(CONNIDSTR,'(''1..'',I2)') IDEFLT(NUMHOSTS)
            ELSE
              CONNIDSTR='1'
            ENDIF
            CALL STRING_TRIM(CONNIDSTR,IBEG,IEND)
            FORMAT='($,'' Enter connection id #s for remote hosts ['
     '        //CONNIDSTR(IBEG:IEND)//']: '',20I3)'
            IF(IOTYPE.EQ.3) THEN
              DO nhost=1,NUMHOSTS
                IDATA(nhost)=ICONNID(nhost)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NUMHOSTS,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nhost=1,NUMHOSTS
C               Check that connection IDs match up for existing hosts.
                IF(ICONNID(nhost).NE.IDATA(nhost)) THEN
                  IF(OLD_HOST(nhost).AND.SOCKET_OPEN(nhost)) THEN
C                   If conn ID doesn't match for open socket, close it.
C                   Signal to slave processes to stop
                    IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,
     '                ICONNID(nhost)).EQ.-1) GOTO 9999
C                    IRET=FSKCLOSE(ICONNID(nhost)) !close socket
                    SOCKET_OPEN(nhost)=.FALSE.
                  ENDIF
                  ICONNID(nhost)=IDATA(nhost)
                ENDIF
              ENDDO
            ENDIF

C           Enter port numbers for remote host machines
            DO nhost=1,NUMHOSTS
              IDEFLT(nhost)=9000+ICONNID(nhost)
            ENDDO
            IF(NUMHOSTS.GT.1) THEN
              WRITE(PORTSTR,'(''9001..'',I4)') IDEFLT(NUMHOSTS)
            ELSE
              PORTSTR='9001'
            ENDIF
            CALL STRING_TRIM(PORTSTR,IBEG,IEND)
            FORMAT='($,'' Enter port numbers for remote hosts ['
     '        //PORTSTR(IBEG:IEND)//']: '',20I5)'
            IF(IOTYPE.EQ.3) THEN
              DO nhost=1,NUMHOSTS
                IDATA(nhost)=IPORT(nhost)
              ENDDO
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,NUMHOSTS,
     '        ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '        LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nhost=1,NUMHOSTS
C               Check that port #s match up for existing hosts.
                IF(IPORT(nhost).NE.IDATA(nhost)) THEN
                  IF(OLD_HOST(nhost).AND.SOCKET_OPEN(nhost)) THEN
C                   If port # doesn't match for open socket, close it.
C                   Signal to slave processes to stop
                    IF(FSKWRITE(QUIT_PROCESS,SK_LONG_INT,1,
     '                ICONNID(nhost)).EQ.-1) GOTO 9999
C                    IRET=FSKCLOSE(ICONNID(nhost)) !close socket
                    SOCKET_OPEN(nhost)=.FALSE.
                  ENDIF
                  IPORT(nhost)=IDATA(nhost)
                ENDIF
              ENDDO
            ENDIF

            FORMAT='($,'' Do you want slave processes to spawn '
     '        //'automatically[Y]? '',A)'
            IF(IOTYPE.EQ.3) THEN
              IF(AUTO_SPAWN) THEN
                ADATA(1)='Y'
              ELSE
                ADATA(1)='N'
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,AYES,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              IF(ADATA(1).EQ.'Y') THEN
                AUTO_SPAWN=.TRUE. !spawn slave processes automatically
              ELSE
                AUTO_SPAWN=.FALSE. !dont spawn slave procs automatically
              ENDIF
            ENDIF

            IF(AUTO_SPAWN) THEN
              IDEFLT(1)=8
              FORMAT='($,'' Enter integer delay before opening sockets '
     '          //'(secs) [8]: '',I3)'
              IF(IOTYPE.EQ.3) IDATA(1)=IDELAY_SOCKET
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '          FORMAT,1,
     '          ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '          LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) IDELAY_SOCKET=IDATA(1)
            ENDIF

            IF(SEND_GLOBAL) THEN
C             (re)initialise GLOBAL_SENT to resend global arrays/common
C             blocks as they may have changed
              DO nhost=1,NUMHOSTS
                GLOBAL_SENT(nhost)=.FALSE.
              ENDDO
            ENDIF
          ENDIF
C newe

          FORMAT='('' Specify derivative calculation method [2]:'''//
     '      '/''   (1) Algebraic'''//
     '      '/''   (2) One-sided finite differences'''//
     '      '/''   (3) Central finite differences'''//
     '      '/$,''    '',I1)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP1D
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP1D=IDATA(1)
        ELSE
          IDEFLT(1)=10
          FORMAT='($,'' Enter print level required [10]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP10
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,30,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP10=IDATA(1)

          IDEFLT(1)=0
          FORMAT='($,'' Enter derivative level [0]: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=KTYP1D
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) KTYP1D=IDATA(1)

          IF(KTYP1D.EQ.0) THEN
            RDEFLT(1)=0.0D0
            FORMAT='($,'' Enter the difference interval [computed]: '','
     '        //'E10.3)'
            IF(IOTYPE.EQ.3) RDATA(1)=DIFF_INTERVAL
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,1.0D0,INFO,ERROR,
     &        *9999)
            IF(IOTYPE.NE.3) DIFF_INTERVAL=RDATA(1)
          ELSE
            IDEFLT(1)=0
            FORMAT='($,'' Enter verify level [0]: '',I2)'
            IF(IOTYPE.EQ.3) IDATA(1)=KTYP15
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &        20,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
            IF(IOTYPE.NE.3) KTYP15=IDATA(1)
          ENDIF

        ENDIF
      ELSE
        ITYP9(nr,nx)=ITYP9(COUP_NRLIST(1,nx),nx)
      ENDIF

      FORMAT='('' Specify depth of output required [1]:'''//
     '  '/''   (1) Equilibrium solution information only'''//
     '  '/''   (2) Intermediate solution information also'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=IWRIT2(nr,nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) IWRIT2(nr,nx)=IDATA(1)

      FORMAT='('' Specify type of output required [0]:'''//
     '  '/''   (0) Solver progress only'''//
     '  '/''   (1) Solution vectors'''//
     '  '/''   (2) Solution and residual vectors'''//
     '  '/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=IWRIT3(nr,nx)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) IWRIT3(nr,nx)=IDATA(1)

      CALL EXITS('IPNONL')
      RETURN
 9999 CALL ERRORS('IPNONL',ERROR)
      CALL EXITS('IPNONL')
      RETURN 1
      END


