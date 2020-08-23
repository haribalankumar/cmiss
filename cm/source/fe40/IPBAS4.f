      SUBROUTINE IPBAS4(NBJ,NBH,NEELEM,NHE,NHP,NPNE,NPNODE,
     '  nr,NW,nx,ERROR,*)

C#### Subroutine: IPBAS4
C###  Description:
C###    IPBAS4 inputs basis functions for dependent variables and
C###    calculates NHE and NHP.
C**** For linear problems the following arrays are used:
C**** NHT40(nje,nve,ie), NHV(nvar,ie), NLV(ie,njt), NVE(ie).

C#### Variable: NHV(nvar,ie)
C###  Type: INTEGER
C###  Set_up: IPBAS4
C###  Description:
C###    NHV(nvar,ie) is global variable number from 1 to NHE(ne),
C###    corresponding to dependent variable nvar in element type ie.

C#### Variable: NVE(ie)
C###  Type: INTEGER
C###  Set_up: IPBAS4
C###  Description:
C###    NVE(ie) is the number of local variables used in
C###    element type ie.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'b40.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'titl40.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),NHP(NPM,0:NRM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,NW(NEM,3),nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER IBEG,ICHAR,ie,IEND,INFO,
     '  nc,ne,nh,nhx,nhh2,NH_MAX,NHT40(3,6,12),nhx_MAX,NLV(12,3),
     '  noelem,NOQUES,nvar,NBH_LOCAL(3,12)
      CHARACTER CHAR1*30
      LOGICAL FILEIP

C#### Variable: NHT40(nje,nve,ie)
C###  Type: INTEGER
C###  Set_up: BLK40
C###  Description:
C###    NHT40(nje,nve,ie) is number of global variables required for
C###    element type ie when number of local variables is nve and
C###    number of dimensions is nje.

      DATA NHT40      !NHT40(nje,nve,ie)=#global variables
     ' /1,2,3, 0,2,3, 0,0,3, 9*0,                 !ie=1  truss,cable
     '  1,2,3, 2,2,3, 0,0,3, 9*0,                 !ie=2  batten
     '  0,0,0, 2,2,3, 3,2,3, 4,4,6, 0,0,6, 0,0,6, !ie=3  beam
     '  18*1,                                     !ie=4  link
     '  18*1,                                     !ie=5  membrane
     '  0,1,3, 0,2,3, 0,3,3,9*0,                  !ie=6  plate
     '  0,0,1, 0,0,2, 0,0,3,9*0,                  !ie=7  shell
     '  0,0,1, 0,0,2, 0,0,3,9*0,                  !ie=8  shell/fluid
     '  18*3,                                     !ie=9  3D elasticity
     '  18*1,                                     !ie=10 tank bottom
     '  18*2,                                     !ie=11 plane stress
     '  18*2/                                     !ie=12 plane strain

C#### Variable: NLV(ie,njt)
C###  Type: INTEGER
C###  Set_up: IPBAS4
C###  Description:
C###    NLV(ie,njt) is the number of local variables possible in
C###    element type ie.

      DATA NLV/1,3,2,1,3,3,3,3,1, 1, 2, 2, !njt=1
     '         2,3,2,1,3,3,3,3,2, 1, 2, 2, !njt=2
     '         3,3,3,1,3,3,3,3,3, 1, 2, 2/ !njt=3

      CALL ENTERS('IPBAS4',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IF(IOTYPE.NE.3) THEN !initialize NVE
        DO ie=1,12
          NVE(ie)=0
        ENDDO
      ENDIF

C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
      IF(IOTYPE.EQ.3) THEN
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          ie=NW(ne,1)
          DO nhx=1,NHE(ne)
            nh=NH_LOC(nhx,nx)
            DO nc=1,2
              NBH_LOCAL(nhx,ie)=NBH(nh,nc,ne)
            ENDDO !nc
          ENDDO !nhx
        ENDDO !noelem (ne)
      ENDIF !IOTYPE==3
CC CPB 27/3/96 Trying to fix this mess
CC! Allocate temporary NH_LOC  PJH 2Feb96
CC      DO nhx=1,NH_LOC_MX
CC        NH_LOC(nhx,nx)=NH_LOC(0,0)+nhx
CC      ENDDO !nhx

! Set NBH, NHV and NVE.
      IF(ITYP6(nr,nx).EQ.1) THEN !linear analysis
        DO ie=1,12
          IF(ETYP(ie)) THEN
            CHAR1=TITL42(ie)
            CALL STRING_TRIM(CHAR1,IBEG,IEND)
            FORMAT='('' '//CHAR1(IBEG:IEND)//' '')'
            CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,1,0,NOQUES,FILEIP,
     '        FORMAT,1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     '        0,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
c           CALL INOUT(IOTYPE,IVDU,IFILE,FORMAT,ERROR,*9999)
            nvar=0
            DO nhx=1,NLV(ie,NJT)
C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
CC cpb 27/3/96
CC              nh=NH_LOC(nhx,nx)
C              DO noelem=1,NEELEM(0,nr)
C                ne=NEELEM(noelem,nr)
C                IF(NW(ne,1).EQ.ie) THEN
C                  IF(IOTYPE.EQ.3) THEN
C                    ne_ie=ne
C                  ELSE IF(IOTYPE.NE.3) THEN
CC cpb 27/3/96
CC                    DO nc=1,2
CC                      NBH(nh,nc,ne)=0
CC                    ENDDO
C                    NBH_LOCAL(nhx)=0
C                  ENDIF
C                ENDIF
C              ENDDO
              CHAR1=TITL43(nhx,ie)
              CALL STRING_TRIM(CHAR1,IBEG,IEND)
              FORMAT='($,'' Enter basis function type number for '
     '          //CHAR1(IBEG:IEND)//' [OMIT]: '',I1)'
C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
              IF(IOTYPE.EQ.3) IDATA(1)=NBH_LOCAL(nhx,ie)
CC CPB 27/3/96
CC              IF(IOTYPE.EQ.3) IDATA(1)=NBH(nh,1,ne_ie)
C              IF(IOTYPE.EQ.3) IDATA(1)=NBH(NH_LOC(nhx,nx),1,ne_ie)
              CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,
     '          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,
     '          ICHAR,IDATA,IONE,1,9,LDATA,LDEFLT,
     '          RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
                NBH_LOCAL(nhx,ie)=IDATA(1)
                IF(IDATA(1).NE.0) THEN
                  nvar=nvar+1
C cpb 28/4/96 I think this line is meant to be like this
C                  NHV(nvar,ie)=nh
                  NHV(nvar,ie)=nhx
                  NVE(ie)=NVE(ie)+1
C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
C                  DO noelem=1,NEELEM(0,nr)
C                    ne=NEELEM(noelem,nr)
C                    IF(NW(ne,1).EQ.ie) THEN
CC cpb 27/3/96
CC                      DO nc=1,2
CC                        NBH(nh,nc,ne)=IDATA(1)
CC                      ENDDO
C                      NBH_LOCAL(nhx)=IDATA(1)
C                    ENDIF
C                  ENDDO !noelem
                ENDIF !idata
              ENDIF !iotype
            ENDDO !nh
          ENDIF !etyp(ie)
        ENDDO !ie

      ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear analysis

C cpb 27/3/96 Below will have to be fixed up.

        DO ie=1,12
          NVE(ie)=NJT
          DO nhh2=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NJ_LOC(NJL_GEOM,nhh2,nr)
            NHV(nh,ie)=nh
          ENDDO
        ENDDO
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nhh2=1,NJ_LOC(NJL_GEOM,0,nr)
            nh=NJ_LOC(NJL_GEOM,nhh2,nr)
            DO nc=1,2
              NBH(nh,nc,ne)=NBJ(nh,ne)
            ENDDO
          ENDDO !nhh2
        ENDDO !noelem
      ENDIF !ityp6

      IF(ETYP(3).AND.NJT.LE.2.AND.NVE(3).EQ.1) THEN !Beam elements
        ERROR='>>Not implemented'
        GOTO 9999
c       DO noelem=1,NEELEM(0,nr)             !..with 1 dof
c         ne=NEELEM(noelem,nr)
c         IF(NW(ne,1).EQ.3) THEN
c           nhx=NHV(1,3)
c           DO nc=1,2
c             NBH(NH_LOC(1,nx),nc,ne)=NBH(NH_LOC(nhx,nx),nc,ne) !reset nh to 1
c           ENDDO
c         ENDIF
c       ENDDO
c       NHV(1,3)=1
c     ENDIF

      ELSE IF(ETYP(4).AND.NJT.EQ.1.AND.NVE(4).EQ.2) THEN !some link elements
C cpb 27/3/96 I have no idea what this is meant to be doing. As far as
C I can see NBH(3/4,,) is not even set up.

        NHV(1,4)=1
        NHV(2,4)=2
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nc=1,2
            NBH(NH_LOC(1,nx),nc,ne)=NBH(NH_LOC(3,nx),nc,ne)
            NBH(NH_LOC(2,nx),nc,ne)=NBH(NH_LOC(4,nx),nc,ne)
          ENDDO
        ENDDO !noelem
      ENDIF

      IF(ETYP(6).AND.NJT.EQ.2.AND.NVE(6).EQ.1) THEN !Plate elements

C cpb 27/3/96 I have no idea what this is meant to be doing. As far as
C I can see NBH(3/4,,) is not even set up.

        NHV(1,6)=1
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nc=1,2
            NBH(NH_LOC(1,nx),nc,ne)=NBH(NH_LOC(3,nx),nc,ne)
          ENDDO
        ENDDO
      ENDIF

! Set NHE(ne) from NVE and NHT40.
      NH_MAX=0
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)

        IF(ITYP5(nr,nx).EQ.1) THEN !static analysis
          IF(ITYP6(nr,nx).EQ.1) THEN      !linear analysis
            ie=NW(ne,1)
            IF(ie.EQ.3) THEN       !beam
              NHE(ne)=3
            ELSE IF(ie.EQ.9) THEN  !3D elasticity
              NHE(ne)=NVE(ie)
            ELSE IF(ie.EQ.11) THEN !plane stress temporary
              NHE(ne)=NVE(ie)
            ELSE
              NHE(ne)=NHT40(NJ_LOC(NJL_GEOM,0,nr),NVE(ie),ie)
            ENDIF
          ELSE IF(ITYP6(nr,nx).EQ.2) THEN !nonlinear analysis
            ie=NW(ne,1)
            NHE(ne)=NVE(ie)
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
          IF(ITYP2(nr,nx).LE.2) THEN !Linear or nonlinear elasticity
            NHE(ne)=NJT
          ENDIF

        ELSE IF(ITYP5(nr,nx).EQ.3) THEN !modal analysis
          IF(ITYP2(nr,nx).LE.2) THEN
            NHE(ne)=NJT
          ELSE
            NHE(ne)=1
          ENDIF
        ENDIF

        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' ne='',I2,'' NW='',I2,'
     '      //''' NVE='',I2,'' NHE='',I2)')
     '      ne,NW(ne,1),NVE(ie),NHE(ne)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF
        IF(NHE(ne).GT.NH_MAX) NH_MAX=NHE(ne)
      ENDDO !noelem


C GMH 21-Apr-95 The following code taken from IPBAS5
! Determine max# nhx's
      nhx_MAX=0
      DO noelem=1,NEELEM(0,nr) !to set #dependent variables
        ne=NEELEM(noelem,nr)
        IF(NHE(ne).GT.nhx_MAX) nhx_MAX=NHE(ne)
      ENDDO !noelem
      CALL CALC_NH_LOC(nhx_MAX,nx,ERROR,*9999)

C cpb 27/3/96 Now set up NBH

C KAT 2001-12-07: Need a different NBH_LOCAL for each ie.
C!!! KAT 2001-12-07: There seem to be problems here with differences
C!!! between NLV, NHT40, and NVE.
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        ie=NW(ne,1)
        DO nhx=1,NHE(ne)
          nh=NH_LOC(nhx,nx)
          DO nc=1,2
            NBH(nh,nc,ne)=NBH_LOCAL(nhx,ie)
          ENDDO !nc
        ENDDO !nhx
      ENDDO !noelem (ne)

C     Set NHP(np) from NHE(ne). (Needs to be done after nh_loc setup)
      IF(IOTYPE.NE.3) THEN
        CALL CALC_NHP(NBH,NEELEM,NHE,NHP,NPNE,NPNODE,nr,nx,
     '    ERROR,*9999)
      ENDIF !iotype.ne.3

      CALL EXITS('IPBAS4')
      RETURN
 9999 CALL ERRORS('IPBAS4',ERROR)
      IF(FILEIP) CLOSE(UNIT=IFILE)
      CALL EXITS('IPBAS4')
      RETURN 1
      END


