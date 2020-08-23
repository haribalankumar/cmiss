      SUBROUTINE IPMESH10_SUB(IBT,IDO,INP,LD,LARG_TOL,
     '  MAP_OLD2NEW,NBJ,NEELEM,
     '  NENP,NEP,NKJE,NPE,NPF,NP_INTERFACE,NPNE,NPNODE,
     '  nr,nr_coronary,NRE,NVJE,SE,XA,XE,XID,
     '  XIP,XP,ZD,ERROR,*)

C#### Subroutine: IPMESH10_SUB
C###  Description:
C###    defines mesh parameters for 3D asysmetric coronary meshs
C###     on the surface of the heart.iod host mesh


      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),ICHAR,IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),LD(NDM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NEP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NP_INTERFACE(0:NPM,0:3),
     '  NPE(0:NPM,NEM,NRM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nr_coronary,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XIP(NIM,NPM),
     '  XP(NKM,NVM,NJM,NPM),XE(NSM,NJM),XID(NIM,NDM),
     '  ZD(NJM,NDM)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER INFO,MAP_OLD2NEW(NP_R_M),nb,nb_coronary,
     '  NE,ne_gobal,
     '  ne_max_host,nj,nk,
     '  nn,NOQUES,np_max_host,ns,nodata,noelem,
     '  noelem_fitted,nonode,np,np1,np2,np_gobal,nd
      REAL*8 DIFF,ERR,PXI,TOL,X(3),XI(3),Z(3)
      LOGICAL FILEIP,FITTED,LARG_TOL(NDM)

      CALL ENTERS('IPMESH10_SUB',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      DIFF=0.0015d0
      nb_coronary=9
      ne_max_host=60
      np_max_host=99
      TOL=10.0d0

      FORMAT='('' Enter basis number for elements in the '''//
     '  '/$,''coronary region [9]:    '',I2)'
      IDEFLT(1)=9
      IF(IOTYPE.EQ.3) IDATA(1)=nb_coronary
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) nb_coronary=IDATA(1)

      FORMAT='('' Enter element number to start coronary elements'''//
     '  '/$,'' from [60]:    '',I2)'
      IDEFLT(1)=60
      IF(IOTYPE.EQ.3) IDATA(1)=ne_max_host
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ne_max_host=IDATA(1)

      FORMAT='('' Enter node number to start coronary elements '''//
     '  '/$,''from [99]:    '',I2)'
      IDEFLT(1)=99
      IF(IOTYPE.EQ.3) IDATA(1)=np_max_host
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) np_max_host=IDATA(1)

      FORMAT='('' Enter value for tolerance [10.0]:'''//
     '  '/$,''    '',F4.2)'
      RDEFLT(1)=10.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=TOL
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,20.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TOL=RDATA(1)

      ITYP10(nr_coronary)=ITYP10(nr)
      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      DO ne=1,NET(nr)
        NPE(0,ne,nr_coronary)=0
      ENDDO

      nonode=1
      DO nodata=1,NDT
        LARG_TOL(nodata)=.TRUE.
        ne=LD(nodata)
        np_gobal=NPNODE(nodata,nr_coronary)
        IF(ne.NE.0) THEN
          NPNODE(nonode,nr_coronary)=nonode+np_max_host !NPNODE(0,nr)
          np=NPNODE(nonode,nr_coronary)
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '      NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '      SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
          DO nj=1,NJT
            xi(nj)=xid(nj,nodata)
            XIP(nj,np)=xid(nj,nodata)
          ENDDO
          DO nj=1,NJT
            nb=NBJ(nj,ne)
            XP(1,1,nj,np)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '        nb,1,XI,XE(1,nj))
            IF(nj.EQ.NJ_LOC(NJL_GEOM,3,nr)) THEN
              DO WHILE (XP(1,1,nj,np).GT.(2.0d0*PI))
                XP(1,1,nj,np)=XP(1,1,nj,np)-(2.0d0*PI)
              ENDDO
            ENDIF
            X(nj)=XP(1,1,nj,np)
          ENDDO !nj
          MAP_OLD2NEW(np_gobal)=np
          NEP(np)=ne
          CALL XZ(ITYP10(nr_coronary),X,Z)
          ERR=0.0d0
          DO nj=1,NJT
            ERR=ERR+((ZD(nj,nodata)-Z(nj))**2.0d0)
          ENDDO
          ERR=DSQRT(ERR)
          IF(ERR.LT.TOL) THEN
            nonode=nonode+1
            LARG_TOL(nodata)=.FALSE.
            NPE(0,ne,nr_coronary)=NPE(0,ne,nr_coronary)+1
            NPE(NPE(0,ne,nr_coronary),ne,nr_coronary)=np
          ENDIF
        ENDIF ! (ne.NE.0)
      ENDDO !nodata

      NPNODE(0,nr_coronary)=nonode-1
      NPNODE(0,0)=np_max_host+NPNODE(0,nr_coronary)
      npt(nr_coronary)=NPNODE((nonode-1),nr_coronary)
      npt(0)=NPNODE(0,0)
      noelem_fitted=0

      DO noelem=1,NEELEM(0,nr_coronary)
        ne_gobal=NEELEM(noelem,nr_coronary)
        nb=NBJ(1,ne_gobal)
        NBJ(1,ne_gobal)=nb_coronary
        DO nonode=1,NNT(nb_coronary)
          NPNE(nonode,nb_coronary,ne_gobal)=NPNE(nonode,nb,ne_gobal)
        ENDDO
       ENDDO

      DO noelem=1,NEELEM(0,nr_coronary)
        ne_gobal=NEELEM(noelem,nr_coronary)
        nb=NBJ(1,ne_gobal)
        FITTED=.TRUE.
        DO nonode=1,NNT(nb)
          nd=NPNE(nonode,nb,ne_gobal)-np_max_host !NPNODE(0,nr)
          IF((LD(nd).EQ.0).OR.LARG_TOL(nd)) THEN
            FITTED=.FALSE.
          ENDIF
        ENDDO !nonode
        IF(FITTED) THEN
          noelem_fitted=noelem_fitted+1
          NEELEM(noelem_fitted,nr_coronary)=
     '      noelem_fitted+ne_max_host !NEELEM(0,nr)
          ne=NEELEM(noelem_fitted,nr_coronary)
          np1=MAP_OLD2NEW(NPNE(1,nb,ne_gobal))
          np2=MAP_OLD2NEW(NPNE(2,nb,ne_gobal))
          nj=NJ_LOC(NJL_GEOM,3,nr)
          IF((DABS(XP(1,1,nj,np1)-XP(1,1,nj,np2)))
     '      .LT.DIFF) THEN
C should ne able to take this out when backend becomes not strictly
C decreasing
C GMH 8/1/97 Update cmgui link
          CALL NODE_CHANGE(np1,.FALSE.,ERROR,*9999)
            XP(1,1,nj,np1)=XP(1,1,nj,np1)+DIFF !
C checks theat angles are not equal
          ENDIF
          IF(DABS(XP(1,1,nj,np2)-XP(1,1,nj,np1)).LT.PI) THEN
            IF(XP(1,1,nj,np1).LT.XP(1,1,nj,np2)) THEN
              NPNE(1,nb,ne)=np2
              NPNE(2,nb,ne)=np1
C swaps locals nodes such that theta is decreasing in x1 to be
C consistant with cordinate system
            ELSE
              NPNE(1,nb,ne)=np1
              NPNE(2,nb,ne)=np2
            ENDIF
          ELSE
           IF(XP(1,1,nj,np1).LT.XP(1,1,nj,np2)) THEN
              NPNE(1,nb,ne)=np1
              NPNE(2,nb,ne)=np2
C swaps locals nodes such that theta is decreasing in x1 to be
C consistant with cordinate system
            ELSE
              NPNE(1,nb,ne)=np2
              NPNE(2,nb,ne)=np1
            ENDIF
          ENDIF ! checks if nodes cross theta equals zero
          DO nj=1,NJT
            NBJ(nj,ne)=NBJ(nj,ne_gobal)
          ENDDO
          NRE(ne)=nr_coronary
          DO nn=1,NNT(nb)
C KAT 23Feb01: now handled by NKJE below
C            DO nk=1,NKT(nn,nb)
C              NKE(nk,nn,nb,ne)=nk
C            ENDDO !nk
            DO nj=1,NJ_LOC(NJL_GEOM,0,nr_coronary)+1 !geom+fibre
              NVJE(nn,nb,nj,ne)=1 !version 1 of nn,nj in elem ne
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nj
          ENDDO !nn
          DO ns=1,NST(nb)+NAT(nb)
            SE(ns,nb,ne)=1.0d0 !scaling factor is 1
          ENDDO !ns
        ENDIF !fitted
      ENDDO !noelem
      NET(nr_coronary)=NEELEM(noelem_fitted,nr_coronary)
      NEELEM(0,nr_coronary)=noelem_fitted
      NEELEM(0,0)=ne_max_host+NEELEM(0,nr_coronary)
      NET(0)= NEELEM(0,0)

      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)

      CALL EXITS('IPMESH10_SUB')
      RETURN
 9999 CALL ERRORS('IPMESH10_SUB',ERROR)
      CALL EXITS('IPMESH10_SUB')
      RETURN 1
      END


