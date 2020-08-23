      SUBROUTINE HIDIPO(DIPOLE_LIST,ISEG,ISDIPO,ISDIPA,NDIPOLES,NRLIST,
     '  NXLIST,STRING,ERROR,*)

C#### Subroutine: HIDIPO
C###  Description:
C###    HIDIPO hides dipole segments.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER DIPOLE_LIST(0:NDIPOLEM),ISEG(*),
     '  ISDIPO(NWM,NDIPOLEM,NGRSEGM),ISDIPA(NWM,NDIPOLEM,NGRSEGM),
     '  NDIPOLES(NRM,NXM),NRLIST(0:NRM),NXLIST(0:NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER DIPOLE_NUM,IBEG,IEND,IFROMC,iw,
     '  IWK(6),N3CO,n,noiw,no_n,no_nrlist,nr,NTIW,nx,nxc
      LOGICAL ALL_REGIONS,CBBREV,PATH,VECTOR

      CALL ENTERS('HIDIPO',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM hide dipole
C###  Description:
C###    Hide dipole segments on specified workstations.
C###  Parameter:    <on (all/WS#s)[all]>
C###    Specify which windows to hide the dipole segments on.
C###    The default is to hide the dipole segments on all windows.
C###  Parameter:    <dipole (all/DIPOLE#)[all]>
C###    This parameter allows you to sepcify to hide either an
C###    individual dipole or all dipoles.
C###  Parameter:    <(vector/path)[vector]>
C###    Specify whether the dipole vectors are to be hidden or
C###    the dipole paths are to be hidden.
C###  Parameter:    <region (all/#s)[1]>
C###    Specify the region numbers to use. The "all" keyword will use
C###    all currently defined regions.
C###  Parameter:    <class #[1]>
C###    Specify the class number (of solve type) to use.

        OP_STRING(1)=STRING(1:IEND)//' <on (all/WS#s)[all>]'
        OP_STRING(2)=BLANK(1:15)//'<dipole (all/#s)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<(vector/path)[vector]>'
        OP_STRING(4)=BLANK(1:15)//'<region (all/#s)[1]>'
        OP_STRING(5)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe22','doc','HIDIPO',ERROR,*9999)
      ELSE
        CALL WS_LIST(IWK,0,NTIW,noco,NTCO,CO,ERROR,*9999)
        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)

        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nxc=NXLIST(1)
        CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
        CALL ASSERT(nx.GT.0,'>>No nx defined for this class',
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

        DO noiw=1,NTIW
          iw=IWK(noiw)
          IF(IWKS(iw).GT.0) THEN
            CALL ACWK(iw,1,ERROR,*9999)
            DO no_nrlist=1,NRLIST(0)
              nr=NRLIST(no_nrlist)
              CALL ASSERT(NDIPOLES(nr,nx).GT.0,
     '          '>>No dipoles in this region',ERROR,*9999)
              IF(DIPOLE_NUM.EQ.0) THEN !All dipoles in the region
                DIPOLE_LIST(0)=NDIPOLES(nr,nx)
                DO n=1,NDIPOLES(nr,nx)
                  DIPOLE_LIST(n)=n
                ENDDO !n
              ELSE
                CALL ASSERT(DIPOLE_NUM.GT.0.AND.DIPOLE_NUM.LE.
     '            NDIPOLES(nr,nx),'>>Invalid dipole number',ERROR,*9999)
                DIPOLE_LIST(0)=1
                DIPOLE_LIST(1)=DIPOLE_NUM
              ENDIF
              DO no_n=1,DIPOLE_LIST(0)
                n=DIPOLE_LIST(no_n)
                IF(VECTOR) THEN
                  IF(ISDIPO(iw,n,nr).GT.0) THEN
                    CALL VISIB(iw,ISEG,ISDIPO(iw,n,nr),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ELSE IF(PATH) THEN
                  IF(ISDIPA(iw,n,nr).GT.0) THEN
                    CALL VISIB(iw,ISEG,ISDIPA(iw,n,nr),'INVISIBLE',
     '                ERROR,*9999)
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            CALL DAWK(iw,1,ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF

      CALL EXITS('HIDIPO')
      RETURN
 9999 CALL ERRORS('HIDIPO',ERROR)
      CALL EXITS('HIDIPO')
      RETURN 1
      END


