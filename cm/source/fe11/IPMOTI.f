      SUBROUTINE IPMOTI(IBT,NBH,NEELEM,NHE,NHP,
     '  NKH,NKJ,nr,YP,FIX,ERROR,*)

C#### Subroutine: IPMOTI
C###  Description:
C###    IPMOTI inputs motion parameters for region nr.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b12.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
C      INCLUDE 'cmiss$reference:mesh00.cmn'
      INCLUDE 'moti00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBH(NHM,NCM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),NKJ(NJM,NPM),nr
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER I,ICHAR,INFO,N,nb,nc,ne,nh,nhx,nk,
     '  no_coeffs,noelem,NOQUES,np,NT_COEFFS,nx,ny,NYTOT
      CHARACTER CHAR1*1,CHAR2*2,CHAR3*1
      LOGICAL FILEIP

      CALL ENTERS('IPMOTI',*9999)
      nx=1 !temporary
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      CALC_INIT_VOL=.TRUE. !set this to true in IPMOTI and FALSE IN MESH_FLOW

      FORMAT='('' Specify type of variable [2]: '''//
     '  '/''   (1) Geometric coordinates'''//
     '  '/''   (2) Displacements'''//
     '  '/''   (3) Flow'''//
     '  '/''   (4) Lung gas flow'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=2
      IF(IOTYPE.EQ.3) IDATA(1)=KTYP58(nr)
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,4,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) KTYP58(nr)=IDATA(1)

      IF(KTYP58(nr).LE.3) THEN !Mot type is geom coords, displs or flow
        FORMAT='('' Enter motion type [1]:'''//
     '    '/''   (1) Fourier coefficients'''//
     '    '/''   (2) Spreadsheet column'''//
     '    '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=MOTION_TYPE
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MOTION_TYPE=IDATA(1)
        IF(MOTION_TYPE.EQ.1.AND.IOTYPE.EQ.1) THEN
          FORMAT='('' Note: Coeffs are entered in order '
     '      //'a0,a1,b1,a2,b2,.. (a=cos,b=sin)'')'
          WRITE(OP_STRING,FORMAT)
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(KTYP58(nr).LE.2) THEN !Mot type is Geom coords or Displacements
        IDEFLT(1)=1
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter index (IY) for motion array ['//CHAR1
     '    //']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MOTION_IY
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,16,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MOTION_IY=IDATA(1)

        IDEFLT(1)=NJ_LOC(NJL_GEOM,0,nr)
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter the number of variables ['//CHAR1
     '    //']: '',I1)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NHE(NEELEM(1,nr))
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NHM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            NHE(ne)=IDATA(1)
          ENDDO
        ENDIF
        CALL ASSERT(NHE(NEELEM(1,nr)).LE.NHM,
     '    '>>Invalid array dimension. Respecify NHM',ERROR,*9999)
        DO nhx=1,NHE(NEELEM(1,nr))
          nh=NH_LOC(nhx,nx)
          WRITE(CHAR1,'(I1)') nh
          FORMAT='($,'' The basis function type number for motion'
     '      //' variable '//CHAR1(1:1)//' is [1]: '',I1)'
          IF(IOTYPE.EQ.3) THEN
            nc=1 !Temporary AJP 17-12-91
            nb=NBH(nh,nc,NEELEM(1,nr))
            IDATA(1)=NB
          ENDIF
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            nb=IDATA(1)
            CALL ASSERT(IBT(1,NIT(nb),nb).EQ.9,
     '        '>>Not a Fourier basis',ERROR,*9999)
            DO noelem=1,NEELEM(0,nr)
              ne=NEELEM(noelem,nr)
              nc=1 !Temporary AJP 17-12-91
              NBH(nh,nc,ne)=NB
            ENDDO
          ENDIF
        ENDDO
      ELSE IF(KTYP58(nr).EQ.3) THEN !Motion type is Flow
        FORMAT='(/$,'' Enter Fourier basis number [1]: '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NB_MOTION
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBFM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NB_MOTION=IDATA(1)
      ENDIF

      IF(KTYP58(nr).LE.2) THEN !Mot type is Geom coords or Displacements
        np=0
 6100   FORMAT='(/$,'' Enter node number [EXIT]: '',I3)'
        IF(IOTYPE.EQ.3) THEN
          np=np+1
          IF(np.LE.NPT(1)) THEN
            IDATA(1)=NP
          ELSE
            IDATA(1)=0
          ENDIF
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,0,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IDATA(1).NE.0) THEN
          IF(IOTYPE.NE.3) np=IDATA(1)
          NYTOT=0
          nc=1 !Temporary AJP 17-12-91
          DO N=1,NP-1
            DO nhx=1,NHP(N)
              nh=NH_LOC(nhx,nx)
              NYTOT=NYTOT+NKH(nh,N,nc)
            ENDDO
          ENDDO
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nhx=1,NHP(np)
            nh=NH_LOC(nhx,nx)
            WRITE(CHAR1,'(I1)') nh
            nc=1 !Temporary AJP 17-12-91
            nb=NBH(nh,nc,NEELEM(1,nr))
            NT_COEFFS=IBT(2,NIT(nb),nb)
            DO nk=1,NKJ(nh,np) !loops over spatial derivatives
              WRITE(CHAR3,'(I1)') nk
              DO no_coeffs=1,NT_COEFFS
                ny=NYTOT+(no_coeffs-1)*NKJ(nh,np)+NK
                RDEFLT(1)=1.1d6
                WRITE(CHAR2,'(I2)') no_coeffs
                FORMAT='($,'' Fourier coefficient '//CHAR2
     '            //' for variable '//CHAR1//' spatial deriv '//CHAR3
     '            //' is [zero]: '',G12.5)'
                IF(IOTYPE.EQ.3) THEN
                  RDATA(1)=YP(ny,MOTION_IY)
                ENDIF
                CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '            FILEIP,FORMAT,1,
     '            ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     '            IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     '            *9999)
                IF(RDATA(1).LT.1.D6) THEN
                  IF(IOTYPE.NE.3) THEN
                    YP(ny,MOTION_IY)=RDATA(1)
                    FIX(ny,5)=.TRUE.
                  ENDIF
                ENDIF
              ENDDO
            ENDDO
            NYTOT=NYTOT+NKH(nh,NP,nc)
          ENDDO
          GO TO 6100
        ENDIF

      ELSE IF(KTYP58(nr).EQ.3) THEN !Motion type is Flow
        CALL_MESH=.TRUE.
        FORMAT='(/$,'' Enter node number [1]: '',I3)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NP_MOTION
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NP_MOTION=IDATA(1)
        NT_COEFFS=IBT(2,NIT(NB_MOTION),NB_MOTION)
        DO no_coeffs=1,NT_COEFFS
          RDEFLT(1)=0.0D0
          WRITE(CHAR2,'(I2)') no_coeffs
          FORMAT='($,'' Enter Fourier coefficient '//CHAR2//' [0]: '','
     '      //'G12.5)'
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=FLOW_COEFFS(no_coeffs)
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '      FILEIP,FORMAT,1,
     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            FLOW_COEFFS(no_coeffs)=RDATA(1)
          ENDIF
        ENDDO

      ELSE IF(KTYP58(nr).EQ.4) THEN !Motion type is Lung gas flow
        CALL_MESH=.TRUE.
        FORMAT='(/$,'' Enter node number [1]: '',I3)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NP_MOTION
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NPM,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NP_MOTION=IDATA(1)

        FORMAT='(/$,'' Enter number of breaths [1]: '',I3)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NT_CYCLES
        ENDIF
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NT_CYCLES=IDATA(1)

        IF(IOTYPE.EQ.3) THEN
          RDATA(1)=FLOW_COEFFS(1)
        ENDIF
        RDEFLT(1)=0.60d0
        FORMAT='($,'' Enter lung tidal volume [0.6 Litres]: '',E12.3)'
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          FLOW_COEFFS(1)=RDATA(1)
        ENDIF

        IF(IOTYPE.EQ.3) THEN
          DO I=1,3
            RDATA(I)=FLOW_COEFFS(I+1)
          ENDDO
        ENDIF
        RDEFLT(1)=2.0d0
        RDEFLT(2)=0.0d0
        RDEFLT(3)=2.0d0
        RMIN=0.0d0
        FORMAT=
     '   '($,'' Enter durations of inspiration,breath-hold & expiration'
     '    //' [2,0,2 sec]: '',3E12.3)'
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,
     '    FILEIP,FORMAT,3,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          DO I=1,3
            FLOW_COEFFS(I+1)=RDATA(I)
          ENDDO
        ENDIF
        T1=DBLE(NT_CYCLES)*(FLOW_COEFFS(2)+FLOW_COEFFS(3)
     '    +FLOW_COEFFS(4))
        FORMAT='('' Specify type of uniform flow [2]: '''//
     '       '/''   (1) Constant (flow=dv/dt)'''//
     '       '/''   (2) Fit to sinusoid'''//
     '       '/$,''    '',I1)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=IDEFLT(1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '       ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '       LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) FLOW_TYPE=IDATA(1)
      ENDIF

      IF(FILEIP) CALL CLOSEF(7,ERROR,*9999)
      CALL EXITS('IPMOTI')
      RETURN
 9999 CALL ERRORS('IPMOTI',ERROR)
      IF(FILEIP) CLOSE(UNIT=7)
      CALL EXITS('IPMOTI')
      RETURN 1
      END


