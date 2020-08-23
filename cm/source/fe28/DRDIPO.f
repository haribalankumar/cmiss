      SUBROUTINE DRDIPO(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,DIPOLE_LIST,
     '  ISEG,ISDIPO,ISDIPA,NDIPOLES,NRLIST,NXLIST,CSEG,DIPOLE_CEN,
     '  DIPOLE_DIR,STRING,ERROR,*)

C#### Subroutine: DRDIPO
C###  Description:
C###    DRDIPO draws dipoles and their paths.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM),DIPOLE_LIST(0:NDIPOLEM),
     '  ISEG(*),ISDIPO(NWM,NDIPOLEM,NGRSEGM),
     '  ISDIPA(NWM,NDIPOLEM,NGRSEGM),NDIPOLES(NRM,NXM),NRLIST(0:NRM),
     '  NXLIST(0:NXM)
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM)
      CHARACTER CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER DIPOLE_NUM,IBEG,IEND,IFROMC,iw,
     '  IWK(6),N3CO,n,noiw,no_n,no_nrlist,nr,NTIW,nx,nxc
      REAL*8 RFROMC,SCALE
      LOGICAL ALL_REGIONS,CBBREV,PATH,VECTOR

      CALL ENTERS('DRDIPO',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM draw dipole
C###  Description:
C###    This command draws dipoles on a 2d window.
C###  Parameter:  <on (all/WS#s)>[all]
C###    Specify the window number on which to draw the dipoles.
C###    The default is to draw dipoles on all windows.
C###  Parameter:  <dipole (all/#)[all]>
C###    Specify specific dipoles to display. The default is to
C###    display all dipoles defined in the current region.
C###  Parameter:  <(vector/path)[vector]>
C###    Specify if you want to display the current dipoles as
C###    vectors for the current time or if you wish to draw the
C###    end points of the dipole vector over time.
C###  Parameter:  <scale #[1.0]>
C###    This parameter scales the size of the drawn dipole
C###    vectors but not the calculated vectors.
C###  Parameter:  <region (#s/all)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:  <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//' <on (all/WS#s)[all]>'
        OP_STRING(2)=BLANK(1:15)//'<dipole (all/#)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(vector/path)[vector]>'
        OP_STRING(4)=BLANK(1:15)//'<scale #[1.0]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(6)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe28','doc','DRDIPO',ERROR,*9999)
      ELSE
        CALL ASSERT(USE_GRAPHICS.EQ.1,
     '    '>>Set USE_GRAPHICS=1',ERROR,*9999)
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.NE.0,'>>No nx defined for this solve class',
     '    ERROR,*9999)

        IF(CBBREV(CO,'DIPOLE',1,noco+1,NTCO,N3CO)) THEN
          DIPOLE_NUM=IFROMC(CO(N3CO+1))
        ELSE
          DIPOLE_NUM=0
        ENDIF

        IF(CBBREV(CO,'VECTOR',1,noco+1,NTCO,N3CO)) THEN
          VECTOR=.TRUE.
          PATH=.FALSE.
        ELSE IF(CBBREV(CO,'PATH',1,noco+1,NTCO,N3CO)) THEN
          PATH=.TRUE.
          VECTOR=.FALSE.
        ELSE
          VECTOR=.TRUE.
          PATH=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
          SCALE=RFROMC(CO(N3CO+1))
        ELSE
          SCALE=1.0d0
        ENDIF

        CALL ASSERT(NRM.LE.NGRSEGM,'>>Increase NGRSEGM to NRM',ERROR,
     '    *9999)

        DO noiw=1,NTIW
          IW=IWK(noiw)
          CALL ACWK(iw,1,ERROR,*9999)
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(NDIPOLES(nr,nx).GT.0,
     '        '>>No dipoles in this region',ERROR,*9999)
            IF(DIPOLE_NUM.EQ.0) THEN !All dipoles in the region
              DIPOLE_LIST(0)=NDIPOLES(nr,nx)
              DO n=1,NDIPOLES(nr,nx)
                DIPOLE_LIST(n)=n
              ENDDO !n
            ELSE
              CALL ASSERT(DIPOLE_NUM.GT.0.AND.DIPOLE_NUM.LE.
     '          NDIPOLES(nr,nx),'>>Invalid dipole number',ERROR,*9999)
              DIPOLE_LIST(0)=1
              DIPOLE_LIST(1)=DIPOLE_NUM
            ENDIF
            DO no_n=1,DIPOLE_LIST(0)
              n=DIPOLE_LIST(no_n)
              IF(PATH) THEN
                CALL ASSERT(DIPOLE_CEN_NTIME(n,nr,nx).GT.0.OR.
     '            DIPOLE_DIR_NTIME(n,nr,nx).GT.0,
     '            '>>No time dependance defined for this dipole',
     '            ERROR,*9999)
              ENDIF
              CALL SGDIPO(DIPOLE_CEN_NTIME(1,1,nx),
     '          DIPOLE_DIR_NTIME(1,1,nx),ISEG,ISDIPO(iw,n,nr),
     '          ISDIPA(iw,n,nr),iw,n,nr,SCALE,CSEG,
     '          DIPOLE_CEN(1,0,1,1,nx),DIPOLE_DIR(1,0,1,1,nx),
     '          PATH,VECTOR,ERROR,*9999)
            ENDDO !no_n (n)
          ENDDO !no_nrlist (nr)
          CALL DAWK(iw,1,ERROR,*9999)
        ENDDO !noiw
      ENDIF

      CALL EXITS('DRDIPO')
      RETURN
 9999 CALL ERRORS('DRDIPO',ERROR)
      CALL EXITS('DRDIPO')
      RETURN 1
      END


