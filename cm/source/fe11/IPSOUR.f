      SUBROUTINE IPSOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,INP,
     '  NAN,NBH,NDIPOLES,NEELEM,NENQ,nr,nr_grid,NQET,NQLIST,NQNE,NQS,
     '  NQSCNB,NQXI,NWQ,nx,NXQ,AQ,CG,CQ,DIPOLE_CEN,DIPOLE_DIR,DXDXIQ,PG,
     '  PROPQ,SETTIME,WG,XE,XG,XQ,YQ,CALCULATE,ERROR,*)

C#### Subroutine: IPSOUR
C###  Description:
C###    Ipsour Inputs source parameters.  Currently enter dipole
C###    parameters for EEG studies.
C###    IPSOUR now also calculates equivalent dipole sources
C###    from the transmembrane potential in coupled grid problems.

      IMPLICIT NONE

      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'sour00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM),
     '  DIPOLE_DIR_NTIME(NDIPOLEM,NRM),IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),
     '  NBH(NHM,NCM,NEM),NDIPOLES(NRM),NEELEM(0:NE_R_M,0:NRM),
     '  NENQ(0:8,NQM),nr,nr_grid,NQET(NQSCM),NQLIST(0:NQM),NQNE(NEQM,
     '  NQEM),NQS(NEQM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NWQ(8,0:NQM),nx,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CG(NMM,NGM),CQ(NMM,NQM),
     '  DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '  DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM),DXDXIQ(3,3,NQM),
     '  PG(NSM,NUM,NGM,NBM),PROPQ(3,3,4,2,NQM),SETTIME,WG(NGM,NBM),
     '  XE(NSM,NJM),XG(NJM,NUM),XQ(NJM,NQM),YQ(NYQM,NIQM)
      CHARACTER ERROR*(*)
      LOGICAL CALCULATE

!     Local Variables
      INTEGER DEND(3),dip1,dip2,dip3,DIPLIST(0:3),DSIZE(3),DSTART(3),i,
     '  ICHAR,IEND,INFO,j,k,LEN_TRIM,maqrhox,maqrhoy,maqrhoz,n,nb,
     &  NCOUNT,ne,neq,neqx,neqy,neqz,ng,ni,niqV,NITB,nj,njj,njtype,
     &  NLOCFEX,NLOCFEY,NLOCFEZ,nn,NNLIST(0:8),no_neelem,NOQUES,npts,nq,
     &  nql,nqq,nqcount,NQXILOC(3),nq1,nq2,nq3,nu,NU1(0:3),NUMGRID,N3CO,
     &  PT1,SCHEME
      REAL*8 CGEFF(3,3),CPE(8),DDOT,DET,DIPMAG,DNUDX(3,3),DPHIDX(3),
     '  DPHIDXI(3),DPHIDS,DPSIDX(8,3),DS,DSDXI,DVm(3),DXDNU(3,3),
     '  DXDXI(3,3),DXI,DXIDS,DXIDX(3,3),DXIDXI(3),DXNQ(3),
     '  LOCAL_DIPOLE_CENX,LOCAL_DIPOLE_CENY,LOCAL_DIPOLE_CENZ,
     '  LOCAL_DIPOLE_DIRX,LOCAL_DIPOLE_DIRY,LOCAL_DIPOLE_DIRZ,
     '  LOC_CEN(3),LOC_DIR(3),PFACTOR,PFXI,PHI(NQGM),RG,SIGDPHIDXI(3),
     '  SUM,TOTMAG,XI(3)
      CHARACTER CHAR1*1,CHAR2*3,CHAR3*3
      LOGICAL ALLDIPOLES,CBBREV,DIPOLESPERELEM,ERROR_FLAG,FILEIP,
     '  FIXEDPOS,GRID,INTBREAK,OLDCALC,ONEDIPOLE

      DATA NU1/1,2,4,7/

      CALL ENTERS('IPSOUR',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      GRID=.FALSE.
      IF(nr_grid.GT.0) THEN
        GRID=.TRUE.
        IF(IOTYPE.EQ.3) THEN
          DO n=1,NDIPOLES(nr)
            DIPOLE_CEN_NTIME(n,nr)=NUM_DIP_TIMES
            DIPOLE_DIR_NTIME(n,nr)=NUM_DIP_TIMES
          ENDDO
        ENDIF
      ENDIF

      !Any ny value could have a dipole component associated with it
      CALL ASSERT(NZ_GD_M.GE.NYM,'>>Increase NZ_GD_M >= NYM',
     '  ERROR,*9999)

      IF(.NOT.CALCULATE) THEN
        CALL ASSERT(nr.LE.NRM,'>>Increase NRM',ERROR,
     '    *9999)
        WRITE(CHAR1,'(I1)') nr
        FORMAT='($,'' The number of dipole sources in '
     '    //'region '//CHAR1//' is [1]: '',I3)'
        IDEFLT(1)=1
        IF(IOTYPE.EQ.3) IDATA(1)=NDIPOLES(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NDIPOLES(nr)=IDATA(1)

C LKC 4-JAN-2006 Improved error message
C        CALL ASSERT(NDIPOLES(nr).LE.NDIPOLEM,
C     '    '>>Increase NDIPOLEM',ERROR,*9999)
        IF(NDIPOLES(nr).GT.NDIPOLEM) THEN
          IEND=0
          CALL APPENDC(IEND,'>>Increase NDIPOLEM to ',OP_STRING(1))
          CALL APPENDI(IEND,NDIPOLES(nr),OP_STRING(1))
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF
        
        WRITE(CHAR1,'(I1)') NJT
        DO n=1,NDIPOLES(nr)
          WRITE(CHAR2,'(I3)') N
          CDEFLT(1)='N'
          FORMAT='(/$,'' Does the centre of dipole '//CHAR2
     '      //' move [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(DIPOLE_CEN_NTIME(n,nr).GT.0) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              DIPOLE_CEN_NTIME(n,nr)=1 !Will be set later
            ELSE
              DIPOLE_CEN_NTIME(n,nr)=0
            ENDIF
          ENDIF
          FORMAT='($,'' Does the dipole vector of dipole '//CHAR2
     '      //' move [N]? '',A)'
          IF(IOTYPE.EQ.3) THEN
            IF(DIPOLE_DIR_NTIME(n,nr).GT.0) THEN
              ADATA(1)='Y'
            ELSE
              ADATA(1)='N'
            ENDIF
          ENDIF
          CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            IF(ADATA(1).EQ.'Y') THEN
              DIPOLE_DIR_NTIME(n,nr)=1 !Will be set later
            ELSE
              DIPOLE_DIR_NTIME(n,nr)=0
            ENDIF
          ENDIF
          IF(DIPOLE_CEN_NTIME(n,nr).EQ.0) THEN !center does not move
            FORMAT='(/$,'' Enter the centre of dipole '//CHAR2
     '        //': '','//CHAR1//'D12.4)'
            DO nj=1,NJT
              RDEFLT(nj)=0.0d0
            ENDDO
            IF(IOTYPE.EQ.3) THEN
              DO nj=1,NJT
                RDATA(nj)=DIPOLE_CEN(nj,0,N,nr)
              ENDDO !nj
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nj=1,NJT
                DIPOLE_CEN(nj,0,N,nr)=RDATA(nj)
              ENDDO !nj
            ENDIF
          ELSE !centre does move
C cpb 2/5/95 Assume Linear interpolation for now.
            IDEFLT(1)=1
            FORMAT='(/$,'' Enter the number of time points in the '
     '        //'centre path [1]: '',I3)'
            IF(IOTYPE.EQ.3) THEN
              IF(GRID) THEN
                IDATA(1)=DIPOLE_CEN_NTIME(n,nr)-1
              ELSE
                IDATA(1)=DIPOLE_CEN_NTIME(n,nr)
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     '        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) DIPOLE_CEN_NTIME(n,nr)=IDATA(1)

C LKC 4-JAN-2006 Improved error message            
C            CALL ASSERT(DIPOLE_CEN_NTIME(n,nr).LE.NDIPTIMM,
C     '        '>>Increase NDIPTIMM',ERROR,*9999)
            IF(DIPOLE_CEN_NTIME(n,nr).GT.NDIPTIMM) THEN
              IEND=0
              CALL APPENDC(IEND,'>>Increase NDIPTIMM to ',OP_STRING(1))
              CALL APPENDI(IEND,DIPOLE_CEN_NTIME(n,nr),OP_STRING(1))
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              GOTO 9999
            ENDIF
            
            FORMAT='(/$,'' Enter the initial dipole centre for '
     '        //'dipole '//CHAR2//': '','//CHAR1//'D12.4)'
            DO nj=1,NJT
              RDEFLT(nj)=0.0d0
            ENDDO !nj
            IF(GRID) THEN
              npts=1
            ELSE
              npts=0
            ENDIF !GRID
            IF(IOTYPE.EQ.3) THEN
              DO nj=1,NJT
                RDATA(nj)=DIPOLE_CEN(nj,npts,n,nr)
              ENDDO !nj
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
C KAT 12Feb01: using npts instead of 0 for GRID.
              DO nj=1,NJT
                DIPOLE_CEN(nj,npts,n,nr)=RDATA(nj)
              ENDDO !nj
              DIPOLE_CEN(4,npts,n,nr)=0.0d0
            ENDIF
            PT1=npts+1
            DO npts=PT1,DIPOLE_CEN_NTIME(n,nr)
              IF(GRID) THEN
                WRITE(CHAR3,'(I3)') npts-1
              ELSE
                WRITE(CHAR3,'(I3)') npts
              ENDIF
              RDEFLT(1)=0.0d0
              FORMAT='($,'' Enter the time (in s) for time '
     '          //'point '//CHAR3//'for '
     '          //'dipole '//CHAR2//' centre: '',D12.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=DIPOLE_CEN(4,npts,n,nr)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) DIPOLE_CEN(4,npts,n,nr)=RDATA(1)
              FORMAT='($,'' Enter the dipole centre at time point '
     '          //CHAR3//' for dipole '//CHAR2//': '','//CHAR1//'D12.4)'
              DO nj=1,NJT
                RDEFLT(nj)=0.0d0
              ENDDO !nj
              IF(IOTYPE.EQ.3) THEN
                DO nj=1,NJT
                  RDATA(nj)=DIPOLE_CEN(nj,npts,N,nr)
                ENDDO !nj
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO nj=1,NJT
                  DIPOLE_CEN(nj,npts,n,nr)=RDATA(nj)
                ENDDO !nj
              ENDIF
            ENDDO !npts
          ENDIF
          IF(DIPOLE_DIR_NTIME(n,nr).EQ.0) THEN !direction does not move
            FORMAT='(/$,'' Enter the dipole vector p for '
     '        //'dipole '//CHAR2//': '','//CHAR1//'D12.4)'
            DO nj=1,NJT
              RDEFLT(nj)=0.0d0
            ENDDO !nj
            IF(IOTYPE.EQ.3) THEN
              DO nj=1,NJT
                RDATA(nj)=DIPOLE_DIR(nj,0,n,nr)
              ENDDO !nj
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) THEN
              DO nj=1,NJT
                DIPOLE_DIR(nj,0,n,nr)=RDATA(nj)
              ENDDO !nj
            ENDIF
          ELSE
C cpb 2/5/95 Assume Linear interpolation for now.
            IDEFLT(1)=1
            FORMAT='(/$,'' Enter the number of time points in the '
     '        //'vector path [1]: '',I3)'
            IF(IOTYPE.EQ.3) THEN
              IF(GRID) THEN
                IDATA(1)=DIPOLE_DIR_NTIME(n,nr)-1
              ELSE
                IDATA(1)=DIPOLE_DIR_NTIME(n,nr)
              ENDIF
            ENDIF
            CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,
     '        IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,
     '        ERROR,*9999)
            IF(IOTYPE.NE.3) DIPOLE_DIR_NTIME(n,nr)=IDATA(1)

C LKC 4-JAN-2006 Improved error message            
C            CALL ASSERT(DIPOLE_DIR_NTIME(n,nr).LE.NDIPTIMM,
C     '        '>>Increase NDIPTIMM',ERROR,*9999)
            IF(DIPOLE_DIR_NTIME(n,nr).GT.NDIPTIMM) THEN
              IEND=0
              CALL APPENDC(IEND,'>>Increase NDIPTIMM to ',OP_STRING(1))
              CALL APPENDI(IEND,DIPOLE_DIR_NTIME(n,nr),OP_STRING(1))
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              GOTO 9999
            ENDIF
            
            FORMAT='(/$,'' Enter the initial dipole vector p for '
     '        //'dipole '//CHAR2//': '','//CHAR1//'D12.4)'
            DO nj=1,NJT
              RDEFLT(nj)=0.0d0
            ENDDO !nj
            IF(GRID) THEN
              npts=1
            ELSE
              npts=0
            ENDIF !GRID
            IF(IOTYPE.EQ.3) THEN
              DO nj=1,NJT
                RDATA(nj)=DIPOLE_DIR(nj,npts,n,nr)
              ENDDO !nj
            ENDIF
            CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &        FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &        IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,
     '        *9999)
            IF(IOTYPE.NE.3) THEN
              DO nj=1,NJT
C KAT 12Feb01: using npts instead of 0 for GRID.
                DIPOLE_DIR(nj,npts,n,nr)=RDATA(nj)
              ENDDO !nj
C KAT 12Feb01: assuming initial time is zero.
              DIPOLE_DIR(4,npts,n,nr)=0.0d0
            ENDIF
            PT1=npts+1
            DO npts=PT1,DIPOLE_DIR_NTIME(n,nr)
              IF(GRID) THEN
                WRITE(CHAR3,'(I3)') npts-1
              ELSE
                WRITE(CHAR3,'(I3)') npts
              ENDIF
              RDEFLT(1)=0.0d0
              FORMAT='($,'' Enter the time (in s) for time '
     '          //'point '//CHAR3
     '          //' for dipole '//CHAR2//' vector: '',D12.4)'
              IF(IOTYPE.EQ.3) RDATA(1)=DIPOLE_DIR(4,npts,n,nr)
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RZERO(1),RMAX,INFO,
     '          ERROR,*9999)
              IF(IOTYPE.NE.3) DIPOLE_DIR(4,npts,n,nr)=RDATA(1)
              FORMAT='($,'' Enter the dipole vector p at time point '
     '          //CHAR3//' for dipole '//CHAR2//': '','//CHAR1//'D12.4)'
              DO nj=1,NJT
                RDEFLT(nj)=0.0d0
              ENDDO !nj
              IF(IOTYPE.EQ.3) THEN
                DO nj=1,NJT
                  RDATA(nj)=DIPOLE_DIR(nj,npts,n,nr)
                ENDDO !nj
              ENDIF
              CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,
     &          FORMAT,NJT,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,
     &          IMIN,IMAX,LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,
     &          ERROR,*9999)
              IF(IOTYPE.NE.3) THEN
                DO nj=1,NJT
                  DIPOLE_DIR(nj,npts,n,nr)=RDATA(nj)
                ENDDO !nj
              ENDIF
            ENDDO !npts
          ENDIF
        ENDDO !n (# dipoles)


C LKC 29-APR-2003 adding the ability to scale the magnitudes of
C       dipoles which are read in.        
        DO dip1=1,NDIPOLES(nr)          
          DO npts=0,DIPOLE_DIR_NTIME(dip1,nr)
            DO nj=1,NJT
              DIPOLE_DIR(nj,npts,dip1,nr)=DIPOLE_DIR(nj,npts,dip1,nr)*
     '          DIPOLE_SCALE_FACTOR
            ENDDO !nj
          ENDDO !npts
        ENDDO !dipoles
        
      ELSE !calculate dipoles from grid potentials

        IF(CBBREV(CO,'OLD_CALCULATIONS',7,noco+1,NTCO,N3CO)) THEN
          OLDCALC=.TRUE.
        ELSE
          OLDCALC=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ONE_DIPOLE',3,noco+1,NTCO,N3CO)) THEN
          ONEDIPOLE=.TRUE.
          ALLDIPOLES=.FALSE.
        ELSE
          ONEDIPOLE=.FALSE.
        ENDIF

        IF(CBBREV(CO,'ALL_DIPOLES',5,noco+1,NTCO,N3CO)) THEN
          ALLDIPOLES=.TRUE.
          ONEDIPOLE=.FALSE.
        ELSE
          ALLDIPOLES=.FALSE.
        ENDIF

        IF(CBBREV(CO,'FIXED_POSITION',3,noco+1,NTCO,N3CO)) THEN
          FIXEDPOS=.TRUE.
        ELSE
          FIXEDPOS=.FALSE.
        ENDIF

        NITB=NQXI(0,NQS(NEELEM(1,nr_grid)))

        IF(CBBREV(CO,'DIPOLES_PER_ELEMENT',9,noco+1,NTCO,N3CO)) THEN
          DIPOLESPERELEM=.TRUE.
          CALL PARSIL(CO(N3CO+1),3,DIPLIST(0),DIPLIST(1),ERROR,*9999)
          CALL ASSERT(NITB.EQ.DIPLIST(0),
     '      '>>Invalid number of list elements',ERROR,*9999)
        ELSE
          DIPOLESPERELEM=.FALSE.
        ENDIF

        NUM_DIP_TIMES=NUM_DIP_TIMES+1

C LKC 4-JAN-2006 Improved error message        
C        CALL ASSERT(NUM_DIP_TIMES.LE.NDIPTIMM,'>>Increase NDIPTIMM',
C     '    ERROR,*9999)
        IF(NUM_DIP_TIMES.GT.NDIPTIMM) THEN
          IEND=0
          CALL APPENDC(IEND,
     &      '>>Increase NDIPTIMM to AT LEAST ',OP_STRING(1))
          CALL APPENDI(IEND,NUM_DIP_TIMES,OP_STRING(1))
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF

            
        IF(ONEDIPOLE) THEN
          NDIPOLES(nr)=1
        ELSE IF(ALLDIPOLES) THEN
          NDIPOLES(nr)=NQR(2,nr_grid)-NQR(1,nr_grid)+1
        ELSE IF(DIPOLESPERELEM) THEN
          IF(NITB.EQ.1) THEN
            NDIPOLES(nr)=NEELEM(0,nr_grid)*DIPLIST(1)
          ELSE IF(NITB.EQ.2) THEN
            NDIPOLES(nr)=NEELEM(0,nr_grid)*DIPLIST(1)*DIPLIST(2)
          ELSE IF(NITB.EQ.3) THEN
            NDIPOLES(nr)=NEELEM(0,nr_grid)*DIPLIST(1)*DIPLIST(2)*
     '        DIPLIST(3)
          ENDIF
        ELSE
          NDIPOLES(nr)=NEELEM(0,nr_grid)
        ENDIF

C LKC 4-JAN-2006 Improved error message
C        CALL ASSERT(NDIPOLES(nr).LE.NDIPOLEM,'>>Increase NDIPOLEM',
C     '    ERROR,*9999)
        IF(NDIPOLES(nr).GT.NDIPOLEM) THEN
          IEND=0
          CALL APPENDC(IEND,'>>Increase NDIPOLEM to ',OP_STRING(1))
          CALL APPENDI(IEND,NDIPOLES(nr),OP_STRING(1))
          CALL WRITES(IOER,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF
        
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niqV,NIQ_V,ERROR,*9999)
        CALL ASSERT(niqV.GT.0,'>>Zero YQ index',ERROR,*9999)

        IF(.NOT.OLDCALC) THEN !not old calculations

          !Allocate/find AQ indices
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqrhox,MAQ_RHO_X,
     '      ERROR,*9999)
          IF(maqrhox.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqrhox,
     '        MAQ_RHO_X,ERROR,*9999)
          ENDIF
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqrhoy,MAQ_RHO_Y,
     '      ERROR,*9999)
          IF(maqrhoy.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqrhoy,
     '        MAQ_RHO_Y,ERROR,*9999)
          ENDIF
          CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maqrhoz,MAQ_RHO_Z,
     '      ERROR,*9999)
          IF(maqrhoz.EQ.0) THEN
            CALL MAQ_LOC(MAQ_ALLOCATE_AND_LOCK,MAQ_COORD,maqrhoz,
     '        MAQ_RHO_Z,ERROR,*9999)
          ENDIF

          !Initialise grid dipole vectors
C$OMP     PARALLEL DO
C$OMP&    PRIVATE(nq)
C$OMP&    SHARED(AQ,maqrhox,maqrhoy,maqrhoz,NQR,nr_grid)
          DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
            AQ(maqrhox,nq)=0.0d0
            AQ(maqrhoy,nq)=0.0d0
            AQ(maqrhoz,nq)=0.0d0
          ENDDO
C$OMP     END PARALLEL DO

          !Calculate a dipole source at each grid point
          IF(ITYP4(nr_grid,nx).EQ.6) THEN !grid based FE

            nb=NBH(NH_LOC(1,nx),1,NEELEM(1,nr_grid))
            NITB=NIT(nb)

            CALL ASSERT(NITB.GT.1,'>>1D dipoles not defined for GBF',
     '        ERROR,*9999)

            !Loop over global FE's in current region
            ERROR_FLAG=.FALSE.
C$OMP       PARALLEL DO
C$OMP&      PRIVATE(CG,CGEFF,CPE,DET,DNUDX,DPHIDX,DPSIDX,DXDNU,DXDXI,
C$OMP&        DXIDX,i,INTBREAK,j,k,ne,neq,neqx,neqy,neqz,ng,ni,nj,
C$OMP&        njj,njtype,nn,no_neelem,NLOCFEX,NLOCFEY,NLOCFEZ,NNLIST,
C$OMP&        nq,nu,RG,SCHEME,SUM,XE,XG)
C$OMP&      SHARED(AQ,CQ,ERROR,ERROR_FLAG,KTYP32,maqrhox,maqrhoy,
C$OMP&        maqrhoz,nb,NEELEM,NGT,NITB,niqV,NJT,NJ_LOC,NQNE,NQS,
C$OMP&        NQXI,nr_grid,NST,NUT,NU1,NXQ,PG,WG,XQ,YQ)
            DO no_neelem=1,NEELEM(0,nr_grid)
              IF(.NOT.ERROR_FLAG) THEN
                ne=NEELEM(no_neelem,nr_grid)
                SCHEME=NQS(ne)

                !Find the number of local FE's in this element
                IF(NITB.EQ.1) THEN
                  NLOCFEX=NQXI(1,SCHEME)-1
                  NLOCFEY=1
                  NLOCFEZ=1
                ELSE IF(NITB.EQ.2) THEN
                  NLOCFEX=NQXI(1,SCHEME)-1
                  NLOCFEY=NQXI(2,SCHEME)-1
                  NLOCFEZ=1
                ELSE IF(NITB.EQ.3) THEN
                  NLOCFEX=NQXI(1,SCHEME)-1
                  NLOCFEY=NQXI(2,SCHEME)-1
                 NLOCFEZ=NQXI(3,SCHEME)-1
                ENDIF

                !Loop over the local FE's in the global FE
                DO neqz=1,NLOCFEZ
                  DO neqy=1,NLOCFEY
                    DO neqx=1,NLOCFEX
                      neq=neqx+((neqy-1)*NQXI(1,SCHEME))+
     &                  ((neqz-1)*NQXI(1,SCHEME)*NQXI(2,SCHEME))

                      !Create a list of local nodes and
                      !check the current element is valid
                      INTBREAK=.FALSE.
                      IF(NITB.EQ.1) THEN
                        nq=NQNE(ne,neq)
                        IF(NXQ(1,1,nq,1).GT.0) THEN 
                          NNLIST(0)=2
                          NNLIST(1)=nq
                          NNLIST(2)=NXQ(1,1,nq,1)
                        ELSE
                          INTBREAK=.TRUE.
                        ENDIF
                      ELSE IF(NITB.EQ.2) THEN
                        nq=NQNE(ne,neq)
                        IF((NXQ(1,1,nq,1).GT.0).AND.
     &                    (NXQ(2,1,nq,1).GT.0).AND.
     &                    (NXQ(1,1,NXQ(2,1,nq,1),1).GT.0).AND.
     &                    (NXQ(2,1,NXQ(1,1,nq,1),1).GT.0)) THEN
                          NNLIST(0)=4
                          NNLIST(1)=nq
                          NNLIST(2)=NXQ(1,1,nq,1)
                          NNLIST(3)=NXQ(2,1,nq,1)
                          NNLIST(4)=NXQ(1,1,NXQ(2,1,nq,1),1)
                        ELSE
                          INTBREAK=.TRUE.
                        ENDIF
                      ELSE IF(NITB.EQ.3) THEN
                        nq=NQNE(ne,neq)

C LKC 8-MAY-2003 Can get a subscript out of range if
C                       1 of the first 3 terms are -ve due to breakages
C                       Then the subsequent terms will access -ve grid
C                       points.
C MLB Fixed with ABS                    

                       IF((NXQ(1,1,nq,1).GT.0).AND.
     &                    (NXQ(2,1,nq,1).GT.0).AND.
     &                    (NXQ(3,1,nq,1).GT.0).AND.
     &                    (NXQ(1,1,ABS(NXQ(2,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(2,1,ABS(NXQ(1,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(1,1,ABS(NXQ(3,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(3,1,ABS(NXQ(1,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(2,1,ABS(NXQ(3,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(3,1,ABS(NXQ(2,1,nq,1)),1).GT.0).AND.
     &                    (NXQ(1,1,ABS(NXQ(2,1,ABS(NXQ(3,1,nq,1)),1)),1)
     &                    .GT.0).AND.
     &                    (NXQ(1,1,ABS(NXQ(3,1,ABS(NXQ(2,1,nq,1)),1)),1)
     &                    .GT.0).AND.
     &                    (NXQ(3,1,ABS(NXQ(2,1,ABS(NXQ(1,1,nq,1)),1)),1)
     &                    .GT.0))THEN
                      
                          NNLIST(0)=8
                          NNLIST(1)=nq
                          NNLIST(2)=NXQ(1,1,nq,1)
                          NNLIST(3)=NXQ(2,1,nq,1)
                          NNLIST(4)=NXQ(1,1,NXQ(2,1,nq,1),1)
                          NNLIST(5)=NXQ(3,1,nq,1)
                          NNLIST(6)=NXQ(3,1,NXQ(1,1,nq,1),1)
                          NNLIST(7)=NXQ(3,1,NXQ(2,1,nq,1),1)
                          NNLIST(8)=NXQ(3,1,NXQ(1,1,NXQ(2,1,nq,1),1),1)
                        ELSE
                          INTBREAK=.TRUE.
                        ENDIF
                      ENDIF

                      IF(.NOT.INTBREAK) THEN !Element exists
                        !Put the geometry,fibres into XE
                        DO njtype=1,2 !NJL_GEOM,NJL_FIBR
                          DO njj=1,NJ_LOC(njtype,0,nr_grid)
                            nj=NJ_LOC(njtype,njj,nr_grid)
                            DO nn=1,NNLIST(0)
                              XE(nn,nj)=XQ(nj,NNLIST(nn))
                            ENDDO !nn
                          ENDDO !njj
                        ENDDO !njtype

                        !Interpolate the conductivities into CG (CPCG)
                        DO ni=1,NITB
                          IF(KTYP32.EQ.1) THEN !Monodomain
                            DO nn=1,NNLIST(0)
                              nq=NNLIST(nn)
                              CPE(nn)=CQ(ni+2,nq)
                            ENDDO !nn
                          ELSE !Bidomain
                            DO nn=1,NNLIST(0)
                              nq=NNLIST(nn)
                              CPE(nn)=CQ(ni+2,nq)*CQ(ni+5,nq)/
     &                          (CQ(ni+2,nq)+CQ(ni+5,nq))
                            ENDDO !nn
                          ENDIF
                          DO ng=1,NGT(nb)
                            SUM=0.0d0
                            DO nn=1,NNLIST(0)
                              SUM=SUM+PG(nn,1,ng,nb)*CPE(nn)
                            ENDDO !nn
                            CG(ni,ng)=SUM
                          ENDDO !ng
                        ENDDO !ni

                        !Loop over gauss points 
                        DO ng=1,NGT(nb)
                          !Calculate XG from XE (XEXG)
                          DO njtype=1,2 !Loop over geometry and fibres
                            DO njj=1,NJ_LOC(njtype,0,nr_grid)
                              nj=NJ_LOC(njtype,njj,nr_grid)
                              DO nu=1,NUT(nb)
                                XG(nj,nu)=DDOT(NST(nb),PG(1,nu,ng,nb),1,
     &                            XE(1,nj),1)
                              ENDDO !nu
                            ENDDO !njj
                          ENDDO !njtype

                          !Calculate dxi/dx and J (RG) (XGMG)
                          DO ni=1,NITB
                            nu=1+ni*(1+ni)/2
                            DO njj=1,NJ_LOC(NJL_GEOM,0,nr_grid)
                              DXDXI(njj,ni)=XG(njj,nu)
                            ENDDO !njj
                          ENDDO !ni
                          CALL INVERT(NITB,DXDXI,DXIDX,RG)

                          !Calculate dx/dnu and dnu/dx
                          CALL MAT_VEC(NITB,nr_grid,DXDNU(1,1),
     &                      DXDNU(1,2),DXDNU(1,3),DXDXI,XG,ERROR,*100)
                          CALL INVERT(NITB,DXDNU,DNUDX,DET)

                          !Calculate effective conductivity tensor
                          DO j=1,NITB
                            DO i=1,NITB
                              CGEFF(i,j)=0.0d0
                              DO k=1,NITB
                                CGEFF(i,j)=CGEFF(i,j)+CG(k,ng)*
     &                            DXDNU(i,k)*DNUDX(k,j)
                              ENDDO !k
                            ENDDO !j
                          ENDDO !i

                          !Calculate dPsi(n)/dx(j)
                          DO j=1,NITB
                            DO nn=1,NNLIST(0)
                              DPSIDX(nn,j)=0.0d0
                              DO k=1,NITB
                                DPSIDX(nn,j)=DPSIDX(nn,j)+
     &                            PG(nn,NU1(k),ng,nb)*DXIDX(k,j)
                              ENDDO !k
                            ENDDO !nn
                          ENDDO !i

                          !Calculate dPhi/dx(j) = dPsi(n)/dx(j) * Phi(n)
                          DO j=1,NITB
                            DPHIDX(j)=0.0d0
                            DO nn=1,NNLIST(0)
                              DPHIDX(j)=DPHIDX(j)+DPSIDX(nn,j)*
     &                          YQ(NNLIST(nn),niqV)
                            ENDDO !nn
                          ENDDO !j
                      
                          !Calculate the integrand
                          DO nn=1,NNLIST(0)
                            nq=NNLIST(nn)
                            IF(NITB.EQ.2) THEN
                              DO j=1,NJT
                                AQ(maqrhox,nq)=AQ(maqrhox,nq)-
     &                            CGEFF(1,j)*DPHIDX(j)*PG(nn,1,ng,nb)*
     &                            RG*WG(ng,nb)
                                AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-
     &                            CGEFF(2,j)*DPHIDX(j)*PG(nn,1,ng,nb)*
     &                            RG*WG(ng,nb)
                              ENDDO !j
                            ELSE IF(NITB.EQ.3) THEN
                              DO j=1,NJT
                                AQ(maqrhox,nq)=AQ(maqrhox,nq)-
     &                            CGEFF(1,j)*DPHIDX(j)*PG(nn,1,ng,nb)*
     &                            RG*WG(ng,nb)
                                AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-
     &                            CGEFF(2,j)*DPHIDX(j)*PG(nn,1,ng,nb)*
     &                            RG*WG(ng,nb)
                                AQ(maqrhoz,nq)=AQ(maqrhoz,nq)-
     &                            CGEFF(3,j)*DPHIDX(j)*PG(nn,1,ng,nb)*
     &                            RG*WG(ng,nb)
                              ENDDO !j
                            ENDIF !ni
                          ENDDO !nn
                        ENDDO !ng
                      ENDIF !intbreak
                    ENDDO !neqx
                  ENDDO !neqy
                ENDDO !neqz
                GOTO 102
 100            CONTINUE
                ERROR_FLAG=.TRUE.
                IF(ERROR.NE.' ') THEN
                  CALL FLAG_ERROR(0,ERROR(:LEN_TRIM(ERROR)))
                ENDIF
 102            CONTINUE
              ENDIF !not error_flag
            ENDDO !no_neelem
C$OMP       END PARALLEL DO

            IF(ERROR_FLAG) THEN
              ERROR=' '
              GOTO 9999
            ENDIF
            
          ELSE IF(ITYP4(nr_grid,nx).EQ.4) THEN !collocation

            !Find a quadratic basis fxn to use for interpolation
            i=0
            DO nb=1,NBT
              IF(NIT(nb).EQ.NIT(NQSCNB(NQS(NEELEM(1,nr_grid))))) THEN
                DO ni=1,NIT(nb)
                  IF(IBT(1,ni,nb).EQ.1) THEN
                    IF(IBT(2,ni,nb).EQ.2) THEN
                      i=nb
                    ENDIF !quadratic
                  ENDIF !lagrange
                ENDDO !ni
              ENDIF !same # xi's
            ENDDO !nb
            nb=i
            CALL ASSERT(nb.GT.0,'>>Error - no quadratic basis found',
     '        ERROR,*9999)

            !Loop over the global FEs in this region
            DO no_neelem=1,NEELEM(0,nr_grid)
              ne=NEELEM(no_neelem,nr_grid)
              SCHEME=NQS(ne)
              NITB=NQXI(0,SCHEME)

              !Loop over the grid points in this FE
              DO nql=1,NQET(SCHEME)
                nq=NQNE(ne,nql)

                !Need to do this to avoid duplications
                IF(NENQ(1,nq).EQ.ne) THEN
                  
                  !PFACTOR corrects for the overlapping of the 
                  !quadratic patches in the calculations.
                  PFACTOR=2.0d0**DBLE(NITB)

                  !Transferring potential into an array so can 
                  !use PFXI to calculate dphi/dxi.
                  IF(NWQ(1,nq).EQ.0) THEN !internal
                    nqq=nq
                  ELSE !external
                    nqq=NWQ(1,nq)
                  ENDIF

                  !Calculate Vm at each local element node
                  CALL YQZQE(NENQ,NITB,nqq,NXQ(-NIM,0,0,1),NQGM,
     '              YQ(1,niqV),PHI,ERROR,*9999)

                  !Calculate the xi position of the current nq
                  CALL CALC_GRID_XI(NITB,NWQ,nq,NXQ(-NIM,0,0,1),XI,
     '              ERROR,*9999)

                  !Use grid basis fxns to calculate dPhi/dxi
                  DO ni=1,NITB
                    DPHIDXI(ni)=PFXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                INP(1,1,nb),NAN(1,1,nb),nb,NU1(ni),PHI,XI)
                  ENDDO !ni 

                  !Calculate sigma*dPhi/dxi
                  IF(KTYP32.EQ.1) THEN !Monodomain
                    !Use s_i for source scaling
                    DO ni=1,NITB
                      SUM=0.0d0
                      DO nj=1,NITB
                        SUM=SUM+(DPHIDXI(nj)*PROPQ(ni,nj,1,1,nq))
                      ENDDO !nj
                      SIGDPHIDXI(ni)=SUM
                    ENDDO !ni
                  ELSE !Bidomain
                    !Use s_i*s_e/(s_i+s_e) source scaling
                    DO ni=1,NITB
                      SUM=0.0d0
                      DO nj=1,NITB
                        IF(PROPQ(ni,nj,1,1,nq)+PROPQ(ni,nj,1,2,nq).GT.
     '                    ZERO_TOL) THEN
                          SUM=SUM+(DPHIDXI(nj)*PROPQ(ni,nj,1,1,nq)*
     '                      PROPQ(ni,nj,1,2,nq)/(PROPQ(ni,nj,1,1,nq)+
     '                      PROPQ(ni,nj,1,2,nq)))
                        ENDIF
                      ENDDO !nj
                      SIGDPHIDXI(ni)=SUM
                    ENDDO !ni
                  ENDIF

                  !Calculate dXi/dX and the transformation Jacobian
                  DO ni=1,NITB
                    DO nj=1,NITB
                      DXDXI(nj,ni)=DXDXIQ(nj,ni,nq)
                    ENDDO !nj
                  ENDDO !ni
                  CALL INVERT(NITB,DXDXI,DXIDX,DET)
                  DET=DET/PFACTOR

                  !Calculate the dipole vector components
                  IF(NITB.EQ.1) THEN
                    IF(NXQ(1,1,nq,1).NE.0) THEN
                      DO nj=1,NJT
                        DXNQ(nj)=XQ(nj,NXQ(1,1,nq,1))-XQ(nj,nq)
                      ENDDO
                    ELSE
                      DO nj=1,NJT
                        DXNQ(nj)=XQ(nj,nq)-XQ(nj,NXQ(-1,1,nq,1))
                      ENDDO
                    ENDIF
                    CALL NORMALISE(NJT,DXNQ,ERROR,*9999)
                    IF(NJT.EQ.1) THEN
                      AQ(maqrhox,nq)=AQ(maqrhox,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(1)
                    ELSE IF(NJT.EQ.2) THEN
                      AQ(maqrhox,nq)=AQ(maqrhox,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(1)
                      AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(2)
                    ELSE IF(NJT.EQ.3) THEN
                      AQ(maqrhox,nq)=AQ(maqrhox,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(1)
                      AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(2)
                      AQ(maqrhoz,nq)=AQ(maqrhoz,nq)-SIGDPHIDXI(1)*
     '                  DXIDX(1,1)*DET*DXNQ(3)
                    ENDIF
                  ELSE IF(NITB.EQ.2) THEN
                    DO ni=1,NITB
                      AQ(maqrhox,nq)=AQ(maqrhox,nq)-SIGDPHIDXI(ni)*
     '                  DXIDX(ni,1)*DET
                      AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-SIGDPHIDXI(ni)*
     '                  DXIDX(ni,2)*DET
                    ENDDO !ni
                  ELSE IF(NITB.EQ.3) THEN
                    DO ni=1,NITB
                      AQ(maqrhox,nq)=AQ(maqrhox,nq)-SIGDPHIDXI(ni)*
     '                  DXIDX(ni,1)*DET
                      AQ(maqrhoy,nq)=AQ(maqrhoy,nq)-SIGDPHIDXI(ni)*
     '                  DXIDX(ni,2)*DET
                      AQ(maqrhoz,nq)=AQ(maqrhoz,nq)-SIGDPHIDXI(ni)*
     '                  DXIDX(ni,3)*DET
                    ENDDO !ni
                  ENDIF
                ENDIF
              ENDDO !nq
            ENDDO !ne
 
          ENDIF !GBFe/collocation

          !Sum point dipoles as required
          IF(ONEDIPOLE) THEN
            !Create one dipole from all grid point dipoles
            TOTMAG=0.0d0
            DIPOLE_CEN_NTIME(1,nr)=0
            DIPOLE_DIR_NTIME(1,nr)=0
            DO nj=1,NJT
              DIPOLE_CEN(nj,0,1,nr)=0.0d0
              DIPOLE_DIR(nj,0,1,nr)=0.0d0
              LOC_CEN(nj)=0.0d0
              LOC_DIR(nj)=0.0d0
            ENDDO
           
            !Loop over all grid points in this region
            DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
              LOC_DIR(1)=LOC_DIR(1)+AQ(maqrhox,nq)
              IF(NJT.GE.2) LOC_DIR(2)=LOC_DIR(2)+AQ(maqrhoy,nq)
              IF(NJT.GE.3) LOC_DIR(3)=LOC_DIR(3)+AQ(maqrhoz,nq)
 
              DIPMAG=DSQRT((AQ(maqrhox,nq)*AQ(maqrhox,nq))+
     '          (AQ(maqrhoy,nq)*AQ(maqrhoy,nq))+(AQ(maqrhoz,nq)*
     '          AQ(maqrhoz,nq)))
              TOTMAG=TOTMAG+DIPMAG

              IF(FIXEDPOS) THEN
                DO nj=1,NJT
                  LOC_CEN(nj)=LOC_CEN(nj)+XQ(nj,nq)
                ENDDO !nj
              ELSE
                DO nj=1,NJT
                  LOC_CEN(nj)=LOC_CEN(nj)+XQ(nj,nq)*DIPMAG
                ENDDO !nj
              ENDIF
            ENDDO !nq
            
            DO nj=1,NJT
              DIPOLE_DIR(nj,0,1,nr)=LOC_DIR(nj)
              DIPOLE_CEN(nj,0,1,nr)=LOC_CEN(nj)
            ENDDO !nj

            IF(FIXEDPOS) THEN
              DO nj=1,NJT
                DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/
     '            DBLE(NQR(2,nr_grid)-NQR(1,nr_grid)+1)
              ENDDO
            ELSE
              IF(TOTMAG.GT.LOOSE_TOL) THEN
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/TOTMAG
                ENDDO !nj
              ELSE
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=0.0d0
                  DIPOLE_DIR(nj,0,1,nr)=0.0d0
                ENDDO
                DO nq=1,NQT
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)+
     '                XQ(nj,nq)
                  ENDDO !nj
                ENDDO !nq
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/
     '              DBLE(NQR(2,nr_grid)-NQR(1,nr_grid)+1)
                ENDDO
              ENDIF
            ENDIF

C!!! LKC 14-MAY-2003 Note that the new dipole center/orienation
C!!! is always at the last time position specified by NUM_DIP_TIMES
            !Set dipole center for time SETTIME
            DO nj=1,NJT
              DIPOLE_CEN(nj,NUM_DIP_TIMES,1,nr)=DIPOLE_CEN(nj,0,1,nr)
            ENDDO !nj
            DIPOLE_CEN(4,NUM_DIP_TIMES,1,nr)=SETTIME

            !Set dipole vector for time SETTIME
            DO nj=1,NJT
              DIPOLE_DIR(nj,NUM_DIP_TIMES,1,nr)=DIPOLE_DIR(nj,0,1,nr)
            ENDDO !nj
            DIPOLE_DIR(4,NUM_DIP_TIMES,1,nr)=SETTIME

          ELSE IF(ALLDIPOLES) THEN !One dipole per grid point

            NCOUNT=0
            DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
              NCOUNT=NCOUNT+1

              DIPOLE_CEN_NTIME(NCOUNT,nr)=0
              DIPOLE_DIR_NTIME(NCOUNT,nr)=0
           
              LOC_DIR(1)=AQ(maqrhox,nq)
              LOC_DIR(2)=AQ(maqrhoy,nq)
              LOC_DIR(3)=AQ(maqrhoz,nq)
 
              DO nj=1,NJT
                DIPOLE_DIR(nj,0,NCOUNT,nr)=LOC_DIR(nj)
                DIPOLE_CEN(nj,0,NCOUNT,nr)=XQ(nj,nq)
              ENDDO !nj
              
C!!! LKC 14-MAY-2003 Note that the new dipole center/orienation
C!!! is always at the last time position specified by NUM_DIP_TIMES
              !Set dipole center for time SETTIME
              DO nj=1,NJT
                DIPOLE_CEN(nj,NUM_DIP_TIMES,NCOUNT,nr)=
     '            DIPOLE_CEN(nj,0,NCOUNT,nr)
              ENDDO !nj
              DIPOLE_CEN(4,NUM_DIP_TIMES,NCOUNT,nr)=SETTIME

              !Set dipole vector for time SETTIME
              DO nj=1,NJT
                DIPOLE_DIR(nj,NUM_DIP_TIMES,NCOUNT,nr)=
     '            DIPOLE_DIR(nj,0,NCOUNT,nr)
              ENDDO !nj
              DIPOLE_DIR(4,NUM_DIP_TIMES,NCOUNT,nr)=SETTIME
            ENDDO !nq

          ELSE !Multiple dipoles per global FE

            IF(DOP) THEN
              WRITE(OP_STRING,'(''Number of directions '',I2)') 
     '          DIPLIST(0)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              DO nj=1,DIPLIST(0)
                WRITE(OP_STRING,'(''Direction '',I2,'' Dipoles '',I2)')
     '            nj,DIPLIST(nj)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO !nj
            ENDIF

            NCOUNT=0
            IF(.NOT.DIPOLESPERELEM) THEN
              !Default to 1 dipole per elem
              DIPLIST(1)=1
              DIPLIST(2)=1
              DIPLIST(3)=1
            ENDIF

            DO no_neelem=1,NEELEM(0,nr_grid)
              ne=NEELEM(no_neelem,nr_grid)
              SCHEME=NQS(ne)
                
              IF(NITB.EQ.1) THEN
                DSIZE(1)=NQXI(1,SCHEME)/DIPLIST(1)
                DSIZE(2)=1
                DSIZE(3)=1
                DIPLIST(2)=1
                DIPLIST(3)=1
                NQXILOC(1)=NQXI(1,SCHEME)
                NQXILOC(2)=1
                NQXILOC(3)=1
              ELSE IF(NITB.EQ.2) THEN
                DSIZE(1)=NQXI(1,SCHEME)/DIPLIST(1)
                DSIZE(2)=NQXI(2,SCHEME)/DIPLIST(2)
                DSIZE(3)=1
                DIPLIST(3)=1
                NQXILOC(1)=NQXI(1,SCHEME)
                NQXILOC(2)=NQXI(2,SCHEME)
                NQXILOC(3)=1
              ELSE IF(NITB.EQ.3) THEN
                DSIZE(1)=NQXI(1,SCHEME)/DIPLIST(1)
                DSIZE(2)=NQXI(2,SCHEME)/DIPLIST(2)
                DSIZE(3)=NQXI(3,SCHEME)/DIPLIST(3)
                NQXILOC(1)=NQXI(1,SCHEME)
                NQXILOC(2)=NQXI(2,SCHEME)
                NQXILOC(3)=NQXI(3,SCHEME)
              ENDIF

              DO dip3=1,DIPLIST(3)
                DO dip2=1,DIPLIST(2)
                  DO dip1=1,DIPLIST(1)
                    NCOUNT=NCOUNT+1

                    !Make a list of local nq numbers for this dipole
                    DSTART(1)=(dip1-1)*DSIZE(1)+1
                    DSTART(2)=(dip2-1)*DSIZE(2)+1
                    DSTART(3)=(dip3-1)*DSIZE(3)+1
                    DEND(1)=dip1*DSIZE(1)
                    DEND(2)=dip2*DSIZE(2)
                    DEND(3)=dip3*DSIZE(3)
                    IF(dip1.EQ.DIPLIST(1)) DEND(1)=NQXILOC(1)
                    IF(dip2.EQ.DIPLIST(2)) DEND(2)=NQXILOC(2)
                    IF(dip3.EQ.DIPLIST(3)) DEND(3)=NQXILOC(3)

                    NQLIST(0)=0
                    DO nq3=DSTART(3),DEND(3)
                      DO nq2=DSTART(2),DEND(2)
                        DO nq1=DSTART(1),DEND(1)
                          nqq=nq1+(nq2-1)*NQXI(1,SCHEME)+((nq3-1)*
     '                      NQXI(1,SCHEME)*NQXI(2,SCHEME))
                          NQLIST(0)=NQLIST(0)+1
                          NQLIST(NQLIST(0))=nqq
                        ENDDO !nq1
                      ENDDO !nq2
                    ENDDO !nq3

                    !Initialise for this dipole
                    DIPOLE_CEN_NTIME(NCOUNT,nr)=0
                    DIPOLE_DIR_NTIME(NCOUNT,nr)=0
                    TOTMAG=0.0d0
                    nqcount=0
                    DO nj=1,NJT
                      DIPOLE_CEN(nj,0,NCOUNT,nr)=0.0d0
                      DIPOLE_DIR(nj,0,NCOUNT,nr)=0.0d0
                      LOC_CEN(nj)=0.0d0
                      LOC_DIR(nj)=0.0d0
                    ENDDO !nj

                    DO nql=1,NQLIST(0)
                      nq=NQNE(ne,NQLIST(nql))

                      !Need to do this to avoid duplications
                      IF(NENQ(1,nq).EQ.ne) THEN
                        LOC_DIR(1)=LOC_DIR(1)+AQ(maqrhox,nq)
                        IF(NJT.GE.2) LOC_DIR(2)=LOC_DIR(2)+
     '                    AQ(maqrhoy,nq)
                        IF(NJT.GE.3) LOC_DIR(3)=LOC_DIR(3)+
     '                    AQ(maqrhoz,nq)
 
                        DIPMAG=DSQRT((AQ(maqrhox,nq)*AQ(maqrhox,nq))+
     '                    (AQ(maqrhoy,nq)*AQ(maqrhoy,nq))+
     '                    (AQ(maqrhoz,nq)*AQ(maqrhoz,nq)))
                        TOTMAG=TOTMAG+DIPMAG
  
                        IF(FIXEDPOS) THEN
                          DO nj=1,NJT
                            LOC_CEN(nj)=LOC_CEN(nj)+XQ(nj,nq)
                          ENDDO !nj
                        ELSE
                          DO nj=1,NJT
                            LOC_CEN(nj)=LOC_CEN(nj)+XQ(nj,nq)*DIPMAG
                          ENDDO !nj
                        ENDIF
                      ELSE
                        nqcount=nqcount+1
                      ENDIF !first elem for this grid
                    ENDDO !nq

                    DO nj=1,NJT
                      DIPOLE_DIR(nj,0,NCOUNT,nr)=LOC_DIR(nj)
                      DIPOLE_CEN(nj,0,NCOUNT,nr)=LOC_CEN(nj)
                    ENDDO !nj
                    IF(FIXEDPOS) THEN
                      DO nj=1,NJT
                        DIPOLE_CEN(nj,0,NCOUNT,nr)=
     '                    DIPOLE_CEN(nj,0,NCOUNT,nr)/
     '                    DBLE(NQLIST(0)-nqcount)
                      ENDDO !nj
                    ELSE
                      IF(TOTMAG.GT.LOOSE_TOL) THEN
                        DO nj=1,NJT
                          DIPOLE_CEN(nj,0,NCOUNT,nr)=
     '                      DIPOLE_CEN(nj,0,NCOUNT,nr)/TOTMAG
                        ENDDO !nj
                      ELSE
                        DO nj=1,NJT
                          DIPOLE_CEN(nj,0,NCOUNT,nr)=0.0d0
                          DIPOLE_DIR(nj,0,NCOUNT,nr)=0.0d0
                        ENDDO !nj
                        DO nql=1,NQLIST(0)
                          nq=NQNE(ne,NQLIST(nql))
                          IF(NENQ(1,nq).EQ.ne) THEN
                            DO nj=1,NJT
                              DIPOLE_CEN(nj,0,NCOUNT,nr)=
     '                          DIPOLE_CEN(nj,0,NCOUNT,nr)+XQ(nj,nq)
                            ENDDO !nj
                          ENDIF
                        ENDDO !nql
                        DO nj=1,NJT
                          DIPOLE_CEN(nj,0,NCOUNT,nr)=
     '                      DIPOLE_CEN(nj,0,NCOUNT,nr)/
     '                      DBLE(NQLIST(0)-nqcount)
                        ENDDO !nj
                      ENDIF
                    ENDIF

C!!! LKC 14-MAY-2003 Note that the new dipole center/orienation
C!!! is always at the last time position specified by NUM_DIP_TIMES
C                   Set dipole center for time SETTIME
                    DO nj=1,NJT
                      DIPOLE_CEN(nj,NUM_DIP_TIMES,NCOUNT,nr)=
     '                  DIPOLE_CEN(nj,0,NCOUNT,nr)
                    ENDDO !nj
                    DIPOLE_CEN(4,NUM_DIP_TIMES,NCOUNT,nr)=SETTIME

C                   Set dipole vector for time SETTIME
                    DO nj=1,NJT
                      DIPOLE_DIR(nj,NUM_DIP_TIMES,NCOUNT,nr)=
     '                  DIPOLE_DIR(nj,0,NCOUNT,nr)
                    ENDDO !nj
                    DIPOLE_DIR(4,NUM_DIP_TIMES,NCOUNT,nr)=SETTIME
                  ENDDO !dip1
                ENDDO !dip2
              ENDDO !dip3
            ENDDO !ne               

          ENDIF !Number of dipoles

        ELSE !old calculations

          WRITE(OP_STRING,'(''>>Warning - using older calculations'')') 
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          IF(NITB.EQ.1) THEN !1d

            !MLB July 2001
            !Note that for now the 1d dipole calculations create
            !dipoles that are aligned with the Y axis in two dimensions
            !and with the Z axis in three dimensions. Also currently
            !restricted to multiple moving dipoles
            !MLB 26 July 2001 - now generalised

            DO no_neelem=1,NEELEM(0,nr_grid)
              ne=NEELEM(no_neelem,nr_grid)

              DIPOLE_CEN_NTIME(no_neelem,nr)=0
              DIPOLE_DIR_NTIME(no_neelem,nr)=0
              TOTMAG=0.0d0
              DO nj=1,NJT
                DIPOLE_CEN(nj,0,no_neelem,nr)=0.0d0
                DIPOLE_DIR(nj,0,no_neelem,nr)=0.0d0
              ENDDO
              SCHEME=NQS(ne)
              NITB=NQXI(0,SCHEME)
              NUMGRID=NQET(SCHEME)
              LOCAL_DIPOLE_CENX=0.0d0
              LOCAL_DIPOLE_CENY=0.0d0
              LOCAL_DIPOLE_CENZ=0.0d0
              LOCAL_DIPOLE_DIRX=0.0d0
              LOCAL_DIPOLE_DIRY=0.0d0
              LOCAL_DIPOLE_DIRZ=0.0d0

              DO nq=1,NUMGRID
                nq1=NXQ(-1,1,NQNE(ne,nq),1)
                nq2=NXQ(1,1,NQNE(ne,nq),1)
                IF((nq1.GT.0).AND.(nq2.GT.0)) THEN
                  DXI=2.0d0/(DBLE(NQXI(1,NQS(ne)))-1.0d0)
                  DPHIDXI(1)=(YQ(nq2,niqV)-YQ(nq1,niqV))/DXI
                  DS=0.0d0
                  DO nj=1,NJT
                    DXNQ(nj)=XQ(nj,nq2)-XQ(nj,nq1)
                    DS=DS+DXNQ(nj)**2
                  ENDDO
                  CALL NORMALISE(NJT,DXNQ,ERROR,*9999)
                  DS=DSQRT(DS)
                  DSDXI=DS/DXI
                  DXIDS=1.0d0/DSDXI
                  DPHIDS=DPHIDXI(1)*DXIDS
                  DO nj=1,NJT
                    DVm(nj)=DPHIDS*DXNQ(nj)
                  ENDDO
                ELSE
                  DO nj=1,NJT
                    DVm(nj)=0.0d0
                  ENDDO
                ENDIF

                DIPMAG=0.0d0
                IF(NJT.EQ.3) THEN
                  DO nj=1,NJT
                    DIPMAG=DIPMAG+(DVm(nj)**2)
                  ENDDO
                  LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '              ((-DVm(1)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '              ((-DVm(2)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRZ=LOCAL_DIPOLE_DIRZ+
     '              ((-DVm(3)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '            (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  ELSEIF(NJT.EQ.2) THEN
                    DO nj=1,NJT
                      DIPMAG=DIPMAG+(DVm(nj)**2)
                    ENDDO
                  LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '              ((DVm(1)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '              ((DVm(2)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                ENDIF
                DIPMAG=DSQRT(DIPMAG)
                TOTMAG=TOTMAG+DIPMAG

                IF(NJT.EQ.3) THEN
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,NQNE(ne,nq))
     '              *DIPMAG
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,NQNE(ne,nq))
     '              *DIPMAG
                  LOCAL_DIPOLE_CENZ=LOCAL_DIPOLE_CENZ+XQ(3,NQNE(ne,nq))
     '              *DIPMAG
                ELSE
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,NQNE(ne,nq))
     '              *DIPMAG
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,NQNE(ne,nq))
     '              *DIPMAG
                ENDIF
              ENDDO !nq

              IF(NJT.EQ.3) THEN
                DIPOLE_DIR(1,0,no_neelem,nr)=LOCAL_DIPOLE_DIRX
                DIPOLE_DIR(2,0,no_neelem,nr)=LOCAL_DIPOLE_DIRY
                DIPOLE_DIR(3,0,no_neelem,nr)=LOCAL_DIPOLE_DIRZ
                DIPOLE_CEN(1,0,no_neelem,nr)=LOCAL_DIPOLE_CENX
                DIPOLE_CEN(2,0,no_neelem,nr)=LOCAL_DIPOLE_CENY
                DIPOLE_CEN(3,0,no_neelem,nr)=LOCAL_DIPOLE_CENZ
              ELSE
                DIPOLE_DIR(1,0,no_neelem,nr)=LOCAL_DIPOLE_DIRX
                DIPOLE_DIR(2,0,no_neelem,nr)=LOCAL_DIPOLE_DIRY
                DIPOLE_CEN(1,0,no_neelem,nr)=LOCAL_DIPOLE_CENX
                DIPOLE_CEN(2,0,no_neelem,nr)=LOCAL_DIPOLE_CENY
              ENDIF

              IF(TOTMAG.GT.LOOSE_TOL) THEN
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,no_neelem,nr)=
     '              DIPOLE_CEN(nj,0,no_neelem,nr)/TOTMAG
                ENDDO !nj
              ELSE
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,no_neelem,nr)=0.0d0
                  DIPOLE_DIR(nj,0,no_neelem,nr)=0.0d0
                ENDDO
                DO nq=1,NUMGRID
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,no_neelem,nr)=
     '                DIPOLE_CEN(nj,0,no_neelem,nr)+XQ(nj,NQNE(ne,nq))
                  ENDDO !nj
                ENDDO !nq
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,no_neelem,nr)=
     '              DIPOLE_CEN(nj,0,no_neelem,nr)/DBLE(NUMGRID)
                ENDDO
              ENDIF

              !Set dipole center for time SETTIME
              DO nj=1,NJT
                DIPOLE_CEN(nj,NUM_DIP_TIMES,no_neelem,nr)=
     '            DIPOLE_CEN(nj,0,no_neelem,nr)
              ENDDO !nj
              DIPOLE_CEN(4,NUM_DIP_TIMES,no_neelem,nr)=SETTIME

              !Set dipole vector for time SETTIME
              DO nj=1,NJT
                DIPOLE_DIR(nj,NUM_DIP_TIMES,no_neelem,nr)=
     '            DIPOLE_DIR(nj,0,no_neelem,nr)
              ENDDO !nj
              DIPOLE_DIR(4,NUM_DIP_TIMES,no_neelem,nr)=SETTIME

            ENDDO !elem

          ELSE IF(ONEDIPOLE) THEN

            DIPOLE_CEN_NTIME(1,nr)=0
            DIPOLE_DIR_NTIME(1,nr)=0
            TOTMAG=0.0d0
            DO nj=1,NJT
              DIPOLE_CEN(nj,0,1,nr)=0.0d0
              DIPOLE_DIR(nj,0,1,nr)=0.0d0
            ENDDO
            NITB=NJT
            LOCAL_DIPOLE_CENX=0.0d0
            LOCAL_DIPOLE_CENY=0.0d0
            LOCAL_DIPOLE_CENZ=0.0d0
            LOCAL_DIPOLE_DIRX=0.0d0
            LOCAL_DIPOLE_DIRY=0.0d0
            LOCAL_DIPOLE_DIRZ=0.0d0

CC$OMP     PARALLEL DO
CC$&       PRIVATE(DH,DHSQ,DIPMAG,DVm,DX,ni,nj,nq,nq1,nq2,Vm1,Vm2)
CC$&       SHARED(CQ,FIXEDPOS,niqV,NITB,nr_grid,NXQ,XQ,YQ)
CC$&       REDUCTION(+:LOCAL_DIPOLE_CENX,LOCAL_DIPOLE_CENY,
CC$&         LOCAL_DIPOLE_CENZ,LOCAL_DIPOLE_DIRX,
CC$&         LOCAL_DIPOLE_DIRY,LOCAL_DIPOLE_DIRZ,TOTMAG)
            DO nq=NQR(1,nr_grid),NQR(2,nr_grid)
              ne=NENQ(1,nq)
              DO ni=1,NITB
                DXIDXI(ni)=1.0d0/(DBLE(NQXI(ni,NQS(ne)))-1.0d0)
                nq1=NXQ(-ni,1,nq,1)
                nq2=NXQ(ni,1,nq,1)
                IF((nq1.GT.0).AND.(nq2.GT.0)) THEN
                  DPHIDXI(ni)=(YQ(nq2,niqV)-YQ(nq1,niqV))/2.0d0
                ELSE
                  DPHIDXI(ni)=0.0d0
                ENDIF
                DO nj=1,NJT
                  DXDXI(nj,ni)=DXDXIQ(nj,ni,nq)
                ENDDO
              ENDDO
              CALL INVERT(NJT,DXDXI,DXIDX,DET)
              DO nj=1,NJT
                DVm(nj)=0.0d0
                DO ni=1,NITB
                  DVm(nj)=DVm(nj)+(DPHIDXI(ni)*DXIDXI(ni)*DXIDX(ni,nj))
                ENDDO
              ENDDO

              DIPMAG=0.0d0
              IF(NJT.EQ.3) THEN
                DO nj=1,NJT
                  DIPMAG=DIPMAG+(DVm(nj)**2)
                ENDDO
                LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '            ((-DVm(1)*CQ(3,nq)*CQ(6,nq))/(CQ(3,nq)+CQ(6,nq)))
                LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '            ((-DVm(2)*CQ(3,nq)*CQ(6,nq))/(CQ(3,nq)+CQ(6,nq)))
                LOCAL_DIPOLE_DIRZ=LOCAL_DIPOLE_DIRZ+
     '            ((-DVm(3)*CQ(3,nq)*CQ(6,nq))/(CQ(3,nq)+CQ(6,nq)))
              ELSEIF(NJT.EQ.2) THEN
                DO nj=1,NJT
                  DIPMAG=DIPMAG+(DVm(nj)**2)
                ENDDO
                LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '            ((DVm(1)*CQ(3,nq)*CQ(6,nq))/(CQ(3,nq)+CQ(6,nq)))
                LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '            ((DVm(2)*CQ(3,nq)*CQ(6,nq))/(CQ(3,nq)+CQ(6,nq)))
              ENDIF
              DIPMAG=DSQRT(DIPMAG)
              TOTMAG=TOTMAG+DIPMAG

              IF(FIXEDPOS) THEN
                IF(NJT.EQ.3) THEN
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,nq)
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,nq)
                  LOCAL_DIPOLE_CENZ=LOCAL_DIPOLE_CENZ+XQ(3,nq)
                ELSE
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,nq)
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,nq)
                ENDIF
              ELSE
                IF(NJT.EQ.3) THEN
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,nq)*DIPMAG
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,nq)*DIPMAG
                  LOCAL_DIPOLE_CENZ=LOCAL_DIPOLE_CENZ+XQ(3,nq)*DIPMAG
                ELSE
                  LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+XQ(1,nq)*DIPMAG
                  LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+XQ(2,nq)*DIPMAG
                ENDIF
              ENDIF
            ENDDO !nq
CC$OMP     END PARALLEL DO

            IF(NJT.EQ.3) THEN
              DIPOLE_DIR(1,0,1,nr)=LOCAL_DIPOLE_DIRX
              DIPOLE_DIR(2,0,1,nr)=LOCAL_DIPOLE_DIRY
              DIPOLE_DIR(3,0,1,nr)=LOCAL_DIPOLE_DIRZ
              DIPOLE_CEN(1,0,1,nr)=LOCAL_DIPOLE_CENX
              DIPOLE_CEN(2,0,1,nr)=LOCAL_DIPOLE_CENY
              DIPOLE_CEN(3,0,1,nr)=LOCAL_DIPOLE_CENZ
            ELSE
              DIPOLE_DIR(1,0,1,nr)=LOCAL_DIPOLE_DIRX
              DIPOLE_DIR(2,0,1,nr)=LOCAL_DIPOLE_DIRY
              DIPOLE_CEN(1,0,1,nr)=LOCAL_DIPOLE_CENX
              DIPOLE_CEN(2,0,1,nr)=LOCAL_DIPOLE_CENY
            ENDIF

            IF(FIXEDPOS) THEN
              DO nj=1,NJT
                DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/
     '            DBLE(NQR(2,nr_grid)-NQR(1,nr_grid)+1)
              ENDDO
            ELSE
              IF(TOTMAG.GT.LOOSE_TOL) THEN
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/TOTMAG
                ENDDO !nj
              ELSE
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=0.0d0
                  DIPOLE_DIR(nj,0,1,nr)=0.0d0
                ENDDO
                DO nq=1,NQT
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)+
     '                XQ(nj,nq)
                  ENDDO !nj
                ENDDO !nq
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,1,nr)=DIPOLE_CEN(nj,0,1,nr)/
     '              DBLE(NQR(2,nr_grid)-NQR(1,nr_grid)+1)
                ENDDO
              ENDIF
            ENDIF

            !Set dipole center for time SETTIME
            DO nj=1,NJT
              DIPOLE_CEN(nj,NUM_DIP_TIMES,1,nr)=
     '          DIPOLE_CEN(nj,0,1,nr)
            ENDDO !nj
            DIPOLE_CEN(4,NUM_DIP_TIMES,1,nr)=SETTIME

            !Set dipole vector for time SETTIME
            DO nj=1,NJT
              DIPOLE_DIR(nj,NUM_DIP_TIMES,1,nr)=
     '          DIPOLE_DIR(nj,0,1,nr)
            ENDDO !nj
            DIPOLE_DIR(4,NUM_DIP_TIMES,1,nr)=SETTIME

          ELSE
            DO no_neelem=1,NEELEM(0,nr_grid)
              ne=NEELEM(no_neelem,nr_grid)

              DIPOLE_CEN_NTIME(no_neelem,nr)=0
              DIPOLE_DIR_NTIME(no_neelem,nr)=0
              TOTMAG=0.0d0
              DO nj=1,NJT
                DIPOLE_CEN(nj,0,no_neelem,nr)=0.0d0
                DIPOLE_DIR(nj,0,no_neelem,nr)=0.0d0
              ENDDO
              SCHEME=NQS(ne)
              NITB=NQXI(0,SCHEME)
              NUMGRID=NQET(SCHEME)
              LOCAL_DIPOLE_CENX=0.0d0
              LOCAL_DIPOLE_CENY=0.0d0
              LOCAL_DIPOLE_CENZ=0.0d0
              LOCAL_DIPOLE_DIRX=0.0d0
              LOCAL_DIPOLE_DIRY=0.0d0
              LOCAL_DIPOLE_DIRZ=0.0d0

CC$OMP       PARALLEL DO
CC$&         PRIVATE(DH,DHSQ,DIPMAG,DVm,DX,ni,nj,nq,nq1,nq2,Vm1,Vm2)
CC$&         SHARED(CQ,FIXEDPOS,ne,niqV,NITB,no_neelem,NQNE,nr,NUMGRID,
CC$&           XQ,YQ)
CC$&         REDUCTION(+:LOCAL_DIPOLE_CENX,LOCAL_DIPOLE_CENY,
CC$&           LOCAL_DIPOLE_CENZ,LOCAL_DIPOLE_DIRX,
CC$&           LOCAL_DIPOLE_DIRY,LOCAL_DIPOLE_DIRZ,TOTMAG)
              DO nq=1,NUMGRID
                DO ni=1,NITB
                  DXIDXI(ni)=1.0d0/(DBLE(NQXI(ni,NQS(ne)))-1.0d0)
                  nq1=NXQ(-ni,1,NQNE(ne,nq),1)
                  nq2=NXQ(ni,1,NQNE(ne,nq),1)
                  IF((nq1.GT.0).AND.(nq2.GT.0)) THEN
                    DPHIDXI(ni)=(YQ(nq2,niqV)-YQ(nq1,niqV))/2.0d0
                  ELSE
                    DPHIDXI(ni)=0.0d0
                  ENDIF
                  DO nj=1,NJT
                    DXDXI(nj,ni)=DXDXIQ(nj,ni,NQNE(ne,nq))
                  ENDDO
                ENDDO
                CALL INVERT(NJT,DXDXI,DXIDX,DET)
                DO nj=1,NJT
                  DVm(nj)=0.0d0
                  DO ni=1,NITB
                    DVm(nj)=DVm(nj)+
     '                (DPHIDXI(ni)*DXIDXI(ni)*DXIDX(ni,nj))
                  ENDDO
                ENDDO

C MLB 5-Nov-1998
C Note: The freespace potential in 2D from a given point with
C               respect to a dipole has a sign opposite to the same case
C               in 3D. This comes from the sign opposites in the Greens
C               functions.
                DIPMAG=0.0d0
                IF(NJT.EQ.3) THEN
                  DO nj=1,NJT
                    DIPMAG=DIPMAG+(DVm(nj)**2)
                  ENDDO
                  LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '              ((-DVm(1)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '              ((-DVm(2)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRZ=LOCAL_DIPOLE_DIRZ+
     '              ((-DVm(3)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                ELSEIF(NJT.EQ.2) THEN
                  DO nj=1,NJT
                    DIPMAG=DIPMAG+(DVm(nj)**2)
                  ENDDO
                  LOCAL_DIPOLE_DIRX=LOCAL_DIPOLE_DIRX+
     '              ((DVm(1)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                  LOCAL_DIPOLE_DIRY=LOCAL_DIPOLE_DIRY+
     '              ((DVm(2)*CQ(3,NQNE(ne,nq))*CQ(6,NQNE(ne,nq)))/
     '              (CQ(3,NQNE(ne,nq))+CQ(6,NQNE(ne,nq))))
                ENDIF
                DIPMAG=DSQRT(DIPMAG)
                TOTMAG=TOTMAG+DIPMAG

                IF(FIXEDPOS) THEN
                  IF(NJT.EQ.3) THEN
                    LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+
     '                XQ(1,NQNE(ne,nq))
                    LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+
     '                XQ(2,NQNE(ne,nq))
                    LOCAL_DIPOLE_CENZ=LOCAL_DIPOLE_CENZ+
     '                XQ(3,NQNE(ne,nq))
                  ELSE
                    LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+
     '                XQ(1,NQNE(ne,nq))
                    LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+
     '                XQ(2,NQNE(ne,nq))
                  ENDIF
                ELSE
                  IF(NJT.EQ.3) THEN
                    LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+
     '                XQ(1,NQNE(ne,nq))*DIPMAG
                    LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+
     '                XQ(2,NQNE(ne,nq))*DIPMAG
                    LOCAL_DIPOLE_CENZ=LOCAL_DIPOLE_CENZ+
     '                XQ(3,NQNE(ne,nq))*DIPMAG
                  ELSE
                    LOCAL_DIPOLE_CENX=LOCAL_DIPOLE_CENX+
     '                XQ(1,NQNE(ne,nq))*DIPMAG
                    LOCAL_DIPOLE_CENY=LOCAL_DIPOLE_CENY+
     '                XQ(2,NQNE(ne,nq))*DIPMAG
                  ENDIF
                ENDIF
              ENDDO !nq
CC$OMP       END PARALLEL DO

              IF(NJT.EQ.3) THEN
                DIPOLE_DIR(1,0,no_neelem,nr)=LOCAL_DIPOLE_DIRX
                DIPOLE_DIR(2,0,no_neelem,nr)=LOCAL_DIPOLE_DIRY
                DIPOLE_DIR(3,0,no_neelem,nr)=LOCAL_DIPOLE_DIRZ
                DIPOLE_CEN(1,0,no_neelem,nr)=LOCAL_DIPOLE_CENX
                DIPOLE_CEN(2,0,no_neelem,nr)=LOCAL_DIPOLE_CENY
                DIPOLE_CEN(3,0,no_neelem,nr)=LOCAL_DIPOLE_CENZ
              ELSE
                DIPOLE_DIR(1,0,no_neelem,nr)=LOCAL_DIPOLE_DIRX
                DIPOLE_DIR(2,0,no_neelem,nr)=LOCAL_DIPOLE_DIRY
                DIPOLE_CEN(1,0,no_neelem,nr)=LOCAL_DIPOLE_CENX
                DIPOLE_CEN(2,0,no_neelem,nr)=LOCAL_DIPOLE_CENY
              ENDIF

              IF(FIXEDPOS) THEN
                DO nj=1,NJT
                  DIPOLE_CEN(nj,0,no_neelem,nr)=
     '              DIPOLE_CEN(nj,0,no_neelem,nr)/DBLE(NUMGRID)
                ENDDO !nj
              ELSE
                IF(TOTMAG.GT.LOOSE_TOL) THEN
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,no_neelem,nr)=
     '                DIPOLE_CEN(nj,0,no_neelem,nr)/TOTMAG
                  ENDDO !nj
                ELSE
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,no_neelem,nr)=0.0d0
                    DIPOLE_DIR(nj,0,no_neelem,nr)=0.0d0
                  ENDDO
                  DO nq=1,NUMGRID
                    DO nj=1,NJT
                      DIPOLE_CEN(nj,0,no_neelem,nr)=
     '                  DIPOLE_CEN(nj,0,no_neelem,nr)+XQ(nj,NQNE(ne,nq))
                    ENDDO !nj
                  ENDDO !nq
                  DO nj=1,NJT
                    DIPOLE_CEN(nj,0,no_neelem,nr)=
     '                DIPOLE_CEN(nj,0,no_neelem,nr)/DBLE(NUMGRID)
                  ENDDO
                ENDIF
              ENDIF

              !Set dipole center for time SETTIME
              DO nj=1,NJT
                DIPOLE_CEN(nj,NUM_DIP_TIMES,no_neelem,nr)=
     '            DIPOLE_CEN(nj,0,no_neelem,nr)
              ENDDO !nj
              DIPOLE_CEN(4,NUM_DIP_TIMES,no_neelem,nr)=SETTIME

              !Set dipole vector for time SETTIME
              DO nj=1,NJT
                DIPOLE_DIR(nj,NUM_DIP_TIMES,no_neelem,nr)=
     '            DIPOLE_DIR(nj,0,no_neelem,nr)
              ENDDO !nj
              DIPOLE_DIR(4,NUM_DIP_TIMES,no_neelem,nr)=SETTIME

            ENDDO
          ENDIF
        ENDIF !old/new calcs

        DO ne=1,NDIPOLES(nr)
          DO nj=1,NJT
            DIPOLE_DIR(nj,0,ne,nr)=DIPOLE_DIR(nj,0,ne,nr)*
     '        DIPOLE_SCALE_FACTOR
            DIPOLE_DIR(nj,NUM_DIP_TIMES,ne,nr)=DIPOLE_DIR(nj,0,ne,nr)*
     '        DIPOLE_SCALE_FACTOR
          ENDDO !nj
        ENDDO !dipoles

      ENDIF

      CALL EXITS('IPSOUR')
      RETURN
 9999 CALL ERRORS('IPSOUR',ERROR)
      CALL EXITS('IPSOUR')
      RETURN 1
      END


