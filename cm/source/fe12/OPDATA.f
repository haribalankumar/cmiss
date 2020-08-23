      SUBROUTINE OPDATA(GROUPS,IBT,IDO,INP,LD,NBH,NBJ,NDDL,NDLIST,
     '  NDLT,
     '  NDP,NEELEM,NELIST,NFLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '  NPNODE,nr,NRE,NVHE,NVHP,NVJE,nx,NYNE,NYNP,CURVCORRECT,
     '  EDD,EDDMIN,MU,RADIUS,SE,SQ,THETA,TIME,WD,XA,XE,XID,XP,YP,
     '  ZA,ZD,ZE,ZP,BETWEEN,DATA,Dot_product,ELEMENTS,ERR,FULL,OPSTAT,
     '  UNDEFORMED,RANGE,STATIC,ERROR,*)

C#### Subroutine: OPDATA
C###  Description:
C###    OPDATA outputs data.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'data00.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD(NDM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NDDL(NEM,NDEM),NDLIST(2),NDLT(NEM),NDP(NDM),
     '  NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NFLIST(0:NFM),NHE(NEM),NHP(NPM),NKH(NHM,NPM,NCM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),EDD(NDM),EDDMIN,MU(2),RADIUS,
     '  SE(NSM,NBFM,NEM),SQ(NDM),THETA(2),TIME,WD(NJM,NDM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM),ZA(NAM,NHM,NCM,NEM),ZD(NJM,NDM),ZE(NSM,NHM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),RANGE*9
      LOGICAL BETWEEN,DATA,Dot_product,GROUPS,ELEMENTS,ERR,FULL,
     '  OPSTAT,
     '  STATIC,UNDEFORMED
!     Local Variables
      INTEGER elem,nb,nd,NDDT,nde,NDTOT,NDTOTEL,ne,nf,ngr,
     '  NITB,ni,nj,njj,njo,NJTT,nodata,nolist,nwf
      REAL*8 dMU,dTHETA,GAMA,PXI,Rad_Mu_Theta_2,
     '  SAED,SAEDEL,SMED,SMEDEL,SUM,SQEDEL,STAT(3,5),
     '  X(6),XD(3),Z(6),Z2(3)
      CHARACTER CHAR1*1,CHAR2*1
      LOGICAL DATAOP,FOUND,FOUND1

      CALL ENTERS('OPDATA',*9999)

      WRITE(OP_STRING,'('' Total number of data points defined ='','
     '  //'I8)') NDT
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(CALL_DATA_FIELD) THEN
        WRITE(OP_STRING,'('' Field values defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(CALL_DATA_FIBRE) THEN
        WRITE(OP_STRING,'('' Fibre angles defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(CALL_DATA_SHEET) THEN
        WRITE(OP_STRING,'('' Sheet angles defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(CALC_XI) THEN
        WRITE(OP_STRING,'('' Xi positions calculated'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        WRITE(OP_STRING,'('' Xi positions not calculated'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(CALL_DATA_SHEET) THEN
        IF(CALC_SHEET) THEN
          WRITE(OP_STRING,'('' Sheet angles corrected'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE IF(.NOT.CALC_SHEET) THEN
          WRITE(OP_STRING,'('' Sheet angles not corrected'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      IF(DATA) THEN
        IF(CALL_DATA_FIBRE) THEN
          NJTT=NJT+1
        ELSE IF(CALL_DATA_SHEET) THEN
          NJTT=NJT+3
        ELSE IF(NJ_LOC(njl_fibr,0,nr).GT.0) THEN
          NJTT=NJT+NJ_LOC(njl_fibr,0,nr)
        ELSE IF(CALL_DATA_FIELD) THEN
C LKC 6-JUL-1999 Not right
C          NJTT=NJT+1
          NJTT=NJT+NJ_LOC(NJL_FIEL,0,nr)
        ELSE
          NJTT=NJT
        ENDIF

        WRITE(CHAR1,'(I1)') NJT
        WRITE(CHAR2,'(I1)') NJTT

        FORMAT='(1X,I5,'') nd='',I5,'' NDP='',I5,'' LD='',I5,'
     '    //''' ERR='',D12.4,'
     '    //'/16X,'' ZD:  '','//CHAR2(1:1)//'D11.3,'
     '    //'/20X,'' ('','//CHAR1(1:1)//'D11.3,'')'','
     '    //'/16X,'' WD:  '','//CHAR2(1:1)//'D11.3,'
     '    //'/16X,'' XID: '',3D11.3)'

        SMED=0.0d0
        SQED=0.0d0
        nodata=0
        NDDT=0
        DATAOP=.FALSE.

        IF(NBJ(1,1).GT.0) THEN
          IF(NXIDEF.GT.NIT(NBJ(1,1)))THEN
            !extra xi value assumed to be time.
            CALL YPZP(1,NBH,NEELEM,NHE,NHP,NKH,NPNODE,
     '        nr,NVHP,nx,NYNE,NYNP,YP,ZA,ZP,ERROR,*9999)
          ENDIF
        ENDIF

        IF(CALC_XI) THEN !Xi positions calculated

C LKC 17-APR-2002 Need to setup the "book-keeping" arrays here
C  as they are not always setup, eg. IOSIGN.
C  May as well do that here as those arrays are not really used
C  else where anyway.
          CALL DSTATS(LD,NDDL,NDLT,ERROR,*9999)


          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            nr=NRE(ne)
            IF(NDLT(ne).GT.0) THEN !some data lies within element
              DATAOP=.TRUE.
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              IF(NXIDEF.GT.NIT(NBJ(1,ne)))THEN
                !extra xi value assumed to be time.
                nwf=1
                !temp nw TVK 27/03/2000
                CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVHE(1,1,1,ne),nwf,nx,
     '            CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA(1,1,1,ne),
     '            ZE,ZP,ERROR,*9999)
              ENDIF
              NITB=NIT(NBJ(1,ne))

              IF(DOP) THEN
                WRITE(OP_STRING,*)'nitb=',NITB
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

              DO nde=1,NDLT(ne)
                nd=NDDL(ne,nde)

                IF(DOP) THEN
                  WRITE(OP_STRING,*) 'nd=',nd
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF

C MLB July 2001 - breaks in 2d as xid(3,nd) is not set!
C               IF((STATIC.AND.(NITB.LT.3.OR.
C     '            XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX)).OR.
C     '            (.NOT.STATIC.AND.DABS(XID(NITB,nd)-TIME).LT.1.0d-06.
C     '            AND.(NITB.LE.3.OR.XID(3,nd).GE.XI3MIN.AND.XID(3,nd).
C     '            LE.XI3MAX))) THEN

                FOUND1=.FALSE.
                IF(NITB.EQ.3) THEN
                  IF(XID(3,nd).GE.XI3MIN.AND.XID(3,nd).LE.XI3MAX)
     '              FOUND1=.TRUE.
                ENDIF

                IF((STATIC.AND.(NITB.LT.3.OR.FOUND1)).OR.
     '            (.NOT.STATIC.AND.DABS(XID(NITB,nd)-TIME).LT.
     '            1.0d-06.AND.(NITB.LE.3.OR.FOUND1))) THEN
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    nb=NBJ(nj,ne)

                    IF(DOP) THEN
                      WRITE(OP_STRING,*) 'nb=',nb
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF

                    X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,1,XID(1,nd),XE(1,nj))
                  ENDDO
                  CALL XZ(ITYP10(nr),X,Z)

                  IF(DOP) THEN
                    WRITE(OP_STRING,*) 'Call to XZ Ok.'
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF

                  IF(NXIDEF.GT.NIT(NBJ(1,ne)))THEN
                    !extra xi value assumed to be time.
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      nb=NBH(nj,1,ne)

                      IF(DOP) THEN
                        WRITE(OP_STRING,*) '*nb=',nb
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF

                      Z2(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                  nb,1,XID(1,nd),ZE(1,nj))
                      IF(KTYP58(nr).EQ.2)THEN !displacement
                        Z(nj)=Z(nj)+Z2(nj)
                      ELSE
                        Z(nj)=Z2(nj)
                      ENDIF
                    ENDDO
                  ENDIF
                  SQ(nd)=0.0d0
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    SQ(nd)=SQ(nd)+(Z(nj)-ZD(nj,nd))**2
                  ENDDO
                  EDD(nd)=DSQRT(SQ(nd))
                  SMED=SMED+EDD(nd)
                  SQED=SQED+SQ(nd)
                  nb=NBJ(1,ne)

                  IF(DOP) THEN
                    WRITE(OP_STRING,*) '**nb=',nb
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF

                  IF(EDD(nd).GE.EDDMIN)THEN
                    nodata=nodata+1
                    CALL ZX(ITYP10(nr),ZD(1,nd),XD)
                    WRITE(OP_STRING,FORMAT) nodata,nd,NDP(nd),LD(nd),
     '                EDD(nd),(ZD(nj,nd),nj=1,NJTT),(XD(nj),nj=1,NJT),
     '                (WD(nj,nd),nj=1,NJTT),(XID(ni,nd),ni=1,NXIDEF)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF !static etc
              ENDDO !nde
              NDDT=NDDT+NDLT(ne)
            ENDIF !some data lies within element
          ENDDO !nolist
        ENDIF !calc_xi

        IF(.NOT.DATAOP) THEN
          nr=1
          WRITE(CHAR1,'(I1)') NJT
          WRITE(CHAR2,'(I1)') NJTT

          FORMAT='(''  nd='',I5,'' NDP='',I5,'
     '      //'/10X,'' ZD:  '','//CHAR2(1:1)//'D11.3,'
     '      //'/14X,'' ('','//CHAR1(1:1)//'D11.3,'')'','
     '      //'/10X,'' WD:  '','//CHAR2(1:1)//'D11.3)'

          DO nd=1,NDT
            CALL ZX(ITYP10(nr),ZD(1,nd),XD)
            WRITE(OP_STRING,FORMAT) nd,NDP(nd),
     '        (ZD(nj,nd),nj=1,NJTT),
     '        (XD(nj),nj=1,NJT),
     '        (WD(nj,nd),nj=1,NJTT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO !nd
        ENDIF !dataop

        IF(FULL) THEN
          FORMAT='('' ne='',I5,'' NDLT(ne)='',I4,'' NDDL(ne,nde):'','
     '      //'10I4,:,/(34X,10I4))'
          DO nolist=1,NELIST(0)
            ne=NELIST(nolist)
            WRITE(OP_STRING,FORMAT) ne,NDLT(ne),(NDDL(ne,nde),nde=1,
     '        NDLT(ne))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF !full
      ENDIF !data


C LKC 14-NOV-98 Adding statistics output
C*** Output statistics about the data

      IF(OPSTAT) THEN

C*** Stores the info in STAT(nj,ni) array
C*** where ni=1 (max), ni=2 (min),
C***       ni=3 (average), ni=4 (range)
C***       ni=5 (distance)

        IF(.NOT.BETWEEN) THEN !do all the data
          nd=1
          DO ni=1,3 !max,min average
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              STAT(nj,ni)=ZD(nj,nd)
            ENDDO !nj
          ENDDO !ni
          DO nd=2,NDT ! calc max, min, average
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              IF(STAT(nj,1).LT.ZD(nj,nd)) STAT(nj,1)=ZD(nj,nd)
              IF(STAT(nj,2).GT.ZD(nj,nd)) STAT(nj,2)=ZD(nj,nd)
              STAT(nj,3)=STAT(nj,3)+ZD(nj,nd)
            ENDDO !nj
          ENDDO !nd
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            STAT(nj,3)=STAT(nj,3)/DBLE(NDT)
            STAT(nj,4)=DABS(STAT(nj,2)-STAT(nj,1))
          ENDDO

        ELSE !BETWEEN

C!!! LKC 16-NOV-1998
C!!! Performs a search to do reverse mapping of NDP
C!!! To be changed for f90 when data point are tidied up
          nd=0
          FOUND=.FALSE.
          DO WHILE ((nd.LE.NDT).AND..NOT.FOUND)
            nd=nd+1
            IF(NDLIST(1).EQ.NDP(nd)) FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) THEN
            ERROR='>> First Data point invalid'
            GOTO 9999
          ENDIF

          DO ni=1,3 !max,min average
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
              STAT(nj,ni)=ZD(nj,nd)
            ENDDO !nj
          ENDDO !ni

C!!! LKC 16-NOV-1998
C!!! Performs a search to do reverse mapping of NDP
C!!! To be changed for f90 when data point are tidied up
          nd=0
          FOUND=.FALSE.
          DO WHILE ((nd.LE.NDT).AND..NOT.FOUND)
            nd=nd+1
            IF(NDLIST(2).EQ.NDP(nd)) FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) THEN
            ERROR='>> Second Data point invalid'
            GOTO 9999
          ENDIF
          CALL ASSERT(nd.LE.NDT,
     '      '>> Second data point invalid',ERROR,*9999)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            IF(STAT(nj,1).LT.ZD(nj,nd)) STAT(nj,1)=ZD(nj,nd)
            IF(STAT(nj,2).GT.ZD(nj,nd)) STAT(nj,2)=ZD(nj,nd)
            STAT(nj,3)=STAT(nj,3)+ZD(nj,nd)
          ENDDO !nj
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            STAT(nj,3)=STAT(nj,3)/DBLE(2)
            STAT(nj,4)=DABS(STAT(nj,2)-STAT(nj,1))
          ENDDO
        ENDIF !BETWEEN

        WRITE(OP_STRING(1),'('' '')')
        WRITE(OP_STRING(2),'('' Data Statistics : '')')
        WRITE(OP_STRING(3),'(''   Maximum  : '',3F12.3)')
     '    (STAT(nj,1), nj=1,NJ_LOC(NJL_GEOM,0,nr))
        WRITE(OP_STRING(4),'(''   Minimum  : '',3F12.3)')
     '    (STAT(nj,2), nj=1,NJ_LOC(NJL_GEOM,0,nr))
        WRITE(OP_STRING(5),'(''   Average  : '',3F12.3)')
     '    (STAT(nj,3), nj=1,NJ_LOC(NJL_GEOM,0,nr))
        WRITE(OP_STRING(6),'(''   Range    : '',3F12.3)')
     '    (STAT(nj,4), nj=1,NJ_LOC(NJL_GEOM,0,nr))

        IF(BETWEEN) THEN
          STAT(1,5)=0.D0
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            STAT(1,5)=STAT(1,5)+STAT(nj,4)**2
          ENDDO
          STAT(1,5)=DSQRT(STAT(1,5))
          WRITE(OP_STRING(7),'(''   Distance: '',F12.3)') STAT(1,5)
        ENDIF !BETWEEN
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ENDIF !OPSTAT

      IF(GROUPS) THEN
        WRITE(OP_STRING,'(''Number of data groups = '',I2)')
     '    NTGRDA
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO ngr=1,NTGRDA
          WRITE(OP_STRING, 
     '     '('' Group '',I2,'' Label= '',A30,'' Data= '',I3,'
     '     //''' .. '',I3)') ngr,LAGRDA(ngr),LIGRDA(1,ngr),
     '     LIGRDA(2,ngr)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
      ENDIF

      IF(ERR) THEN !list data errors
        SMED=0.0d0
        SAED=0.0d0
        SQED=0.0d0
        NDTOT=0

C PM 14Aug02 : for face fitting
        IF(KTYP11.EQ.2) THEN !face fitting
          DO nolist=0,NFLIST(0)
            NELIST(nolist)=NFLIST(nolist)
          ENDDO
        ENDIF

        DO nolist=1,NELIST(0)
          IF(KTYP11.EQ.1) THEN !element fitting
            ne=NELIST(nolist)
            nr=NRE(ne)
          ELSEIF(KTYP11.EQ.2) THEN !face fitting
            nf=NELIST(nolist)
            ne=NPF(6,nf)
            nr=NRE(ne)
          ENDIF

C KFA 2001-08-01
C  added code to allow error to be calculated from deformed
C  geometery.
          IF(.NOT.UNDEFORMED) THEN
C            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),
C     '        NKHE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),nr,
C     '        NVHE(1,1,1,ne),NW(ne,1,nx),nx,
c     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),
C     '        ZA(1,1,1,ne),XE,ZP,ERROR,*9999)
            CALL ZPZE(NBH(1,1,ne),1,NHE(ne),NKHE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVHE(1,1,1,ne),1,nx,
     '        CURVCORRECT(1,1,1,ne),SE(1,1,ne),ZA,
     '        ZE,ZP,ERROR,*9999)
          ELSE
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF,NPNE(1,1,ne),
     '        nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          ENDIF
          IF(ELEMENTS) THEN
            SMEDEL=0.0d0
            SAEDEL=0.0d0
            SQEDEL=0.0d0
            NDTOTEL=0
            WRITE(OP_STRING,'('' Element# : '',I5)') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          IF(KTYP11.EQ.1) THEN !element fitting
            elem=ne
          ELSEIF(KTYP11.EQ.2) THEN !face fitting
            elem=nf
          ENDIF

          DO nde=1,NDLT(elem)
            nd=NDDL(elem,nde)
            IF(KTYP8.EQ.0.OR.(KTYP8.LE.1.AND.ITYP6(1,1).LE.1)
     '        .OR.KTYP8.EQ.6) THEN
C             linear geometric fit or fit by optimisation (or no fit)
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                nb=NBJ(nj,ne)
              IF(.NOT.UNDEFORMED) THEN
                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XID(1,nd),ZE(1,nj))
              ELSE
                X(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '            XID(1,nd),XE(1,nj))
              ENDIF
              ENDDO
              CALL XZ(ITYP10(nr),X,Z) !transforms coords to r.c.
              SUM=0.0d0
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                IF(Dot_product) THEN
                  SUM=SUM+(Z(nj)-ZD(nj,nd))*WD(nj,nd)
                ELSE
                  SUM=SUM+(Z(nj)-ZD(nj,nd))**2
                ENDIF
              ENDDO !nj
              IF(Dot_product) THEN
                EDD(nd)=DABS(SUM)
              ELSE
                EDD(nd)=DSQRT(SUM)
              ENDIF
            ELSE
c cpb 7/3/95 This needs to be generalised
              njo=NLH_FIT(1,1,1)
              nb=NBJ(njo,ne)
              CALL ASSERT(nb.GT.0,'>>No field or fibre basis defined',
     '          ERROR,*9999)
              Z(njo)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          XID(1,nd),XE(1,njo))
              EDD(nd)=Z(njo)-ZD(njo,nd)
            ENDIF
C ***       Note: this should check that fibre is actually being fitted
            IF(KTYP8.EQ.2.AND.NJ_LOC(njl_fibr,0,nr).LE.2) THEN
              !fibre fit
              njo=NLH_FIT(1,1,1)
              IF(EDD(nd).GT.PI/2.0d0) THEN !
                EDD(nd)=EDD(nd)-PI
                WRITE(OP_STRING,
     '            '('' Warning!!! Fibre error has pi subtracted'','
     '            //'/'' nd='',I5,'' NDP='',I5,'' Angle='',E12.3,'
     '            //''' ne='',I3,'' XID:'',3F6.2)')
     '            nd,NDP(nd),ZD(njo,nd),ne,(XID(ni,nd),ni=1,3)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ELSE IF(EDD(nd).LT.-PI/2.0D0) THEN
                EDD(nd)=EDD(nd)+PI
                WRITE(OP_STRING,
     '            '('' Warning!!! Fibre error has pi added'','
     '            //'/'' nd='',I5,'' NDP='',I5,'' Angle='',E12.3,'
     '            //''' ne='',I3,'' XID:'',3F6.2)')
     '            nd,NDP(nd),ZD(njo,nd),ne,(XID(ni,nd),ni=1,3)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            IF(DABS(EDD(nd)).GE.EDDMIN) THEN
              IF(FULL) THEN
C*** VYW 15/03/2010 Add more digits to copye with large amount of data
          WRITE(OP_STRING,'('' nd='',I8,'' EDD(nd)='','
     '            //'E12.3)') nd,EDD(nd)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
              SMED=SMED+EDD(nd)
              SAED=SAED+DABS(EDD(nd))
              SQED=SQED+EDD(nd)**2
              NDTOT=NDTOT+1
              IF(ELEMENTS) THEN
                SMEDEL=SMEDEL+EDD(nd)
                SAEDEL=SAEDEL+DABS(EDD(nd))
                SQEDEL=SQEDEL+EDD(nd)**2
                NDTOTEL=NDTOTEL+1
              ENDIF
            ENDIF
          ENDDO !nde
          IF(ELEMENTS) THEN
            IF(NDTOTEL.GT.1) THEN
              WRITE(OP_STRING,'('' Number of data points in fit in '
     '          //'element = '',I8)') NDTOTEL
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '           '(/'' Average error           : '',D12.6,'' +/- '','
     '          //'D12.6,'
     '          //'/'' Average absolute error  : '',D12.6,'' +/- '','
     '          //'D12.6,'
     '          //'/'' Root mean squared error : '',D12.6)')
     '          SMEDEL/DBLE(NDTOTEL),
     '          DSQRT((SQEDEL-SMEDEL**2/DBLE(NDTOTEL))/DBLE(NDTOTEL-1)),
     '          SAEDEL/DBLE(NDTOTEL),DSQRT((SQEDEL-SAEDEL**2/
     '          DBLE(NDTOTEL))/DBLE(NDTOTEL-1)),DSQRT(SQEDEL/
     '          DBLE(NDTOTEL))
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              WRITE(OP_STRING,'('' no data points in fit in element'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO !list of elements
        IF(NDTOT.GT.1) THEN
          WRITE(OP_STRING,'(/'' Number of data points in fit ='',I8)')
     '      NDTOT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          IF(KTYP8.EQ.1.AND.ITYP6(nr,nx).EQ.1) THEN !linear geom fit
            WRITE(OP_STRING,'('' Errors based on Euclidean norm'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,
     '       '(/'' Average error           : '',D12.6,'' +/- '',D12.6,'
     '      //'/'' Average absolute error  : '',D12.6,'' +/- '',D12.6,'
     '      //'/'' Root mean squared error : '',D12.6)')
     '      SMED/DBLE(NDTOT),
     '      DSQRT((SQED-SMED**2/DBLE(NDTOT))/DBLE(NDTOT-1)),
     '      SAED/DBLE(NDTOT),DSQRT((SQED-SAED**2/DBLE(NDTOT))/
     '      DBLE(NDTOT-1)),DSQRT(SQED/DBLE(NDTOT))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C LKC 11-APR-1999 added else
        ELSE
          WRITE(OP_STRING,'('' No data points in any elements'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF !ndtot>1
      ENDIF !list data errors

      IF(RANGE(1:9).EQ.'RECTANGLE') THEN !list data within rectangle
        WRITE(OP_STRING,'('' Data within range mu='',F6.3,'','','
     '    //'F6.3,'' and theta='',F6.3,'','',F6.3)')
     '    MU(1),MU(2),THETA(1),THETA(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nd=1,NDT
          CALL ZX(ITYP10(nr),ZD(1,nd),X) !X=pos.n in prolate coords
          IF((X(2).GE.   MU(1).AND.X(2).LE.   MU(2)).AND   !mu in range
     '      .(X(3).GE.THETA(1).AND.X(3).LE.THETA(2))) THEN !& theta
            WRITE(OP_STRING,'(/'' nd='',I6,'' Prolate coords:'','
     '        //'3E12.3)') nd,X(1),X(2),X(3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' (Element'',I5,'' Xi-coords:'','
     '        //'3F5.2,'')'')') LD(nd),(XID(ni,nd),ni=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Angles:'',3E12.3)')
     '        (ZD(NJ_LOC(njl_fibr,njj,nr),nd),njj=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !data within specified range
        ENDDO !nd

      ELSE IF(RANGE(1:6).EQ.'CIRCLE') THEN !list data within circle
        WRITE(OP_STRING,'('' Data within circle radius='',F6.3,'
     '    //''' centre mu='',F6.3,'',theta='',F6.3)')
     '    RADIUS,MU(1),THETA(1)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nd=1,NDT
          CALL ZX(ITYP10(nr),ZD(1,nd),X) !X=pos.n in prolate coords
          dMU   =X(2)-MU(1)    !Mu    diff betw data & specified pt
          dTHETA=X(3)-THETA(1) !Theta diff betw data & specified pt
          Rad_Mu_Theta_2=dMU*dMU+dTHETA*dTHETA !radius squared
          IF(Rad_Mu_Theta_2.LE.RADIUS*RADIUS) THEN !in range
            WRITE(OP_STRING,'(/'' nd='',I6,'' Prolate coords:'','
     '        //'3E12.3)') nd,X(1),X(2),X(3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' (Element'',I5,'' Xi-coords:'','
     '        //'3F5.2,'')'')') LD(nd),(XID(ni,nd),ni=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' Angles:'',3E12.3)')
     '        (ZD(NJ_LOC(NJL_FIBR,njj,nr),nd),njj=1,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            ne=LD(nd)             !element#
            nj=NJ_LOC(NJL_FIBR,3,nr) !nj index for sheet angle
            nb=NBJ(nj,ne)         !basis#  for sheet angle
            CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            GAMA=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '       nb,1,XID(1,nd),XE(1,nj))
            EDD(nd)=ZD(nj,nd)-GAMA
            WRITE(OP_STRING,'('' Gama ='',E12.3,'' ZD-GAMA ='',E12.3)')
     '        GAMA,EDD(nd)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ENDIF !data within specified range
        ENDDO !nd
      ENDIF !range

      CALL EXITS('OPDATA')
      RETURN
 9999 CALL ERRORS('OPDATA',ERROR)
      CALL EXITS('OPDATA')
      RETURN 1
      END


