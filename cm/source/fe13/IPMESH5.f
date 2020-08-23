      SUBROUTINE IPMESH5(NBJ,NEELEM,NENP,NKJE,NKJ,
     '  NP_INTERFACE,NPNE,NPNODE,NRE,NVJE,NVJP,SE,XP,ERROR,*)

C#### Subroutine: IPMESH5
C###  Description:
C###    IPMESH5 defines mesh parameters for 3d open-ended cylindrical
C###    mesh(s).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh05.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKJ(NJM,NPM),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,m1,m2,nb,ne,nj,nk,nn,
     '  ns,np,NP2,NPTOP,noelem,nonode,NOQUES,nocylinder,
     '  NP_INCR,NESTART,NPSTART,nr
      REAL*8 THETA,S1_DIR(3),S2_DIR(3)
      CHARACTER CHAR1*1
      LOGICAL FILEIP,FIRST

      CALL ENTERS('IPMESH5',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999
      FORMAT='(/$,'' Enter basis function number of the '
     '  //'bicubic hermite BE basis [1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=MESH5_NB
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,9,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MESH5_NB=IDATA(1)
      FORMAT='(/$,'' Enter the number of cylinders [1]: '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=NCYLINDERS
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,10,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NCYLINDERS=IDATA(1)
      FORMAT='(/$,'' Enter the height of the cylinders [1.0]: '',D12.4)'
      RDEFLT(1)=1.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=HEIGHT
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,10,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) HEIGHT=RDATA(1)

      DO nocylinder=1,NCYLINDERS
        WRITE(CHAR1,'(I1)') nocylinder
        FORMAT='(/$,'' Enter the radius of cylinder '//CHAR1
     '    //' ['//CHAR1//']: '',D12.4)'
        RDEFLT(1)=DBLE(nocylinder)
        IF(IOTYPE.EQ.3) RDATA(1)=MESH5_RAD(nocylinder)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH5_RAD(nocylinder)=RDATA(1)

        FORMAT='(/$,'' Enter the number of elements around cylinder '
     '   //CHAR1//' [4]: '',I2)'
        IDEFLT(1)=4
        IF(IOTYPE.EQ.3) IDATA(1)=MESH5_S(nocylinder,1)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,100,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH5_S(nocylinder,1)=IDATA(1)
        FORMAT='(/$,'' Enter the number of elements up cylinder '
     '   //CHAR1//' [2]: '',I2)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=MESH5_S(nocylinder,2)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,100,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH5_S(nocylinder,2)=IDATA(1)
      ENDDO !End of input for all cylinders

C***  Calculate nodal coordinates
      np=1
      DO nocylinder=1,NCYLINDERS
        DO m2=0,MESH5_S(nocylinder,2) !step up
          DO m1=0,MESH5_S(nocylinder,1)-1 !step around
            THETA=DBLE(m1)/DBLE(MESH5_S(nocylinder,1))*2.0d0*PI
            IF(DOP) THEN
              WRITE(OP_STRING,'('' theta ='',D12.4)') THETA
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
            XP(1,1,1,np)=MESH5_RAD(nocylinder)*DCOS(THETA)
            XP(1,1,2,np)=MESH5_RAD(nocylinder)*DSIN(THETA)
            XP(1,1,3,np)=DBLE(m2)/DBLE(MESH5_S(nocylinder,2))*HEIGHT
C           Calculate derivatives
            S1_DIR(1)=-DSIN(THETA)
            S1_DIR(2)=DCOS(THETA)
            S1_DIR(3)=0.0d0
            S2_DIR(1)=0.0d0
            S2_DIR(2)=0.0d0
            S2_DIR(3)=1.0d0
            IF(nocylinder.EQ.1.AND.NCYLINDERS.GT.1)THEN
C             s2 is in the opposite direction
              DO nj=1,NJT
                S2_DIR(nj)=-S2_DIR(nj)
              ENDDO
            ENDIF
            DO nj=1,NJT
              XP(2,1,nj,np)=S1_DIR(nj)
              XP(3,1,nj,np)=S2_DIR(nj)
              XP(4,1,nj,np)=0.0d0 !no cross derivatives
            ENDDO
            np=np+1
            IF(DOP) THEN
              WRITE(OP_STRING,'('' m1='',I3)')m1
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' np='',I5,'' xp(nk,1,1,np):'','
     '          //'2E12.3)') np,(XP(nk,1,1,np),nk=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' np='',I5,'' xp(nk,1,2,np):'','
     '          //'2E12.3)') np,(XP(nk,1,2,np),nk=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' np='',I5,'' xp(nk,1,3,np):'','
     '          //'3E12.3)') np,(XP(nk,1,3,np),nk=1,4)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO
        ENDDO
        np=np-1
        IF(NCYLINDERS.EQ.1..OR.nocylinder.GT.1)THEN
          nr=nocylinder-1
          IF(nr.EQ.0)nr=1
          CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
          IF(NPNODE(0,nr).GT.0) THEN
            CALL REGION_DESTROY(nr,ERROR,*9999)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' nr='',I2,'' nocylinder='',I3)')
     '        nr,nocylinder
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          NPT(nr)=np !Maximum node number in region NOCYLINDER/2
          IF(nr.EQ.1)THEN
            NPNODE(0,nr)=NPT(nr) !number nodes in region NR
            NPSTART=0
          ELSE
            NPSTART=NPT(NR-1)-MESH5_S(nocylinder-1,1)*
     '        (MESH5_S(nocylinder-1,2)+1)
C           Node numbers on interface are unique
            NPNODE(0,nr)=NPT(nr)-NPSTART
          ENDIF
          DO nonode=1,NPNODE(0,nr)
            NPNODE(nonode,nr)=NPSTART+nonode
          ENDDO
          IF(DOP) THEN
            WRITE(OP_STRING,'(''  npnode(0,nr)='',I5)')NPNODE(0,nr)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'(''  npnode(1..,nr)='',10I5)')
     '        (NPNODE(nonode,nr),nonode=1,NPNODE(0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        np=np+1
      ENDDO !End of node construction
      NPT(0)=np-1
      NPNODE(0,0)=np-1
      DO np=1,NPT(0)
        DO nj=1,NJT
          NKJ(nj,np)=4
        ENDDO
      ENDDO
C     Construct elements and npne array
      np=1
      nocylinder=1
      ne=1
      NPTOP=0
      FIRST=.TRUE.
      DO WHILE (nocylinder.LE.NCYLINDERS)
        IF(nocylinder.EQ.1.AND.NCYLINDERS.GT.1)THEN
C         Need to have s2 direction reversed
          np=MESH5_S(nocylinder,1)*MESH5_S(nocylinder,2)+1
          !Number from top down in this case
          NP2=MESH5_S(nocylinder,1)
          DO m2=1,MESH5_S(nocylinder,2) !vertical step
            DO m1=1,MESH5_S(nocylinder,1) !theta steps
              DO nj=1,NJT
                NBJ(nj,ne)=MESH5_NB
              ENDDO
              IF(m1.EQ.MESH5_S(nocylinder,1))THEN
                NP_INCR=1-MESH5_S(nocylinder,1)
              ELSE
                NP_INCR=1
              ENDIF
              NPNE(1,MESH5_NB,ne)=np
              NPNE(2,MESH5_NB,ne)=np+NP_INCR
              NPNE(3,MESH5_NB,ne)=np-NP2
              NPNE(4,MESH5_NB,ne)=np-NP2+NP_INCR
              IF(DOP) THEN
                WRITE(OP_STRING,'('' np='',I5,'' ne='',I5,'' m2='','
     '            //'I3,'' m1='',I3)') np,ne,m2,m1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' npne(nn,'',I2,'',ne):'',4I4)')
     '            MESH5_NB,(NPNE(nn,MESH5_NB,ne),nn=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              ne=ne+1
              IF(m1.LT.MESH5_S(nocylinder,1))THEN
                np=np+1
              ELSEIF(m2.LT.MESH5_S(nocylinder,2))THEN
                np=np-2*NP2+1
              ELSE !at bottom of a cylinder
                NPTOP=MESH5_S(nocylinder,1)*
     '            (MESH5_S(nocylinder,2)+1)+1
                np=NPTOP
              ENDIF
            ENDDO !end m1
          ENDDO !end m2
        ELSE
          DO m2=1,MESH5_S(nocylinder,2) !vertical step
            DO m1=1,MESH5_S(nocylinder,1) !theta steps
              IF(FIRST.AND.m1.EQ.1.AND.m2.EQ.1) THEN
                NPTOP=NPTOP+MESH5_S(nocylinder,1)*
     '            (MESH5_S(nocylinder,2)+1)+2
              ENDIF
              DO nj=1,NJT
                NBJ(nj,ne)=MESH5_NB
              ENDDO
              NP2=MESH5_S(nocylinder,1)
              IF(m1.EQ.MESH5_S(nocylinder,1))THEN
                NP_INCR=1-MESH5_S(nocylinder,1)
              ELSE
                NP_INCR=1
              ENDIF
              NPNE(1,MESH5_NB,ne)=np
              NPNE(2,MESH5_NB,ne)=np+NP_INCR
              NPNE(3,MESH5_NB,ne)=np+NP2
              NPNE(4,MESH5_NB,ne)=np+NP2+NP_INCR
              IF(DOP) THEN
                WRITE(OP_STRING,'('' np='',I5,'' ne='',I5,'' m2='','
     '            //'I3,'' m1='',I3)') np,ne,m2,m1
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' npne(nn,'',I2,'',ne):'',4I4)')
     '            MESH5_NB,(NPNE(nn,MESH5_NB,ne),nn=1,4)
              ENDIF
              ne=ne+1
              IF(m2.LT.MESH5_S(nocylinder,2))THEN
                np=np+1
              ELSEIF(m1.LT.MESH5_S(nocylinder,1))THEN
                np=np+1
              ELSE !at top of a cylinder
                IF(FIRST)THEN
                  np=NPTOP-MESH5_S(nocylinder,1)*
     '             (MESH5_S(nocylinder,2)+1)-1
                ELSE
                  np=NPTOP+1
                ENDIF
              ENDIF
            ENDDO !end m1
          ENDDO !end m2
        ENDIF !end choice of first or other cylinders
        IF(FIRST)THEN
          nr=nocylinder-1
          IF(NCYLINDERS.EQ.1)nr=1
          IF(nr.GE.1)THEN
            NET(nr)=ne-1 !Maximum element number in region nr
            IF(nr.EQ.1)THEN
              NESTART=0
              NEELEM(0,nr)=NET(nr) !number elements in region nr
            ELSE
              NESTART=NET(nr-1)
              NEELEM(0,nr)=NET(nr)-NET(nr-1)
            ENDIF
            DO noelem=1,NEELEM(0,nr)
              NEELEM(noelem,nr)=noelem+NESTART
              NRE(noelem)=nr  !sb NRE(ne)? PJH
            ENDDO
          ENDIF
        ENDIF
        IF(nocylinder.EQ.1)THEN
          nocylinder=nocylinder+1
        ELSE IF(nocylinder.GT.1.AND.nocylinder.LT.NCYLINDERS)THEN
          IF(FIRST)THEN
            FIRST=.FALSE.
C           Need to redo the above since cylinder NOCYLINDER is at the
C           interface of two regions
          ELSE
            nocylinder=nocylinder+1
            FIRST=.TRUE.
          ENDIF
        ELSE
          FIRST=.TRUE.
          nocylinder=nocylinder+1
        ENDIF
      ENDDO !nocylinder

      NET(0)=ne-1
      NEELEM(0,0)=ne-1 !number of elements in all regions

      IF(NCYLINDERS.EQ.1)THEN
        NRT=1
      ELSE
        NRT=NCYLINDERS-1
      ENDIF
      DO nr=1,NRT
        DO nb=1,NBFT
          DO noelem=1,NEELEM(0,nr)
            ne=NEELEM(noelem,nr)
            DO ns=1,NST(nb)+NAT(nb)
              SE(ns,nb,ne)=1.0d0
            ENDDO
            DO nn=1,NNT(nb)
C KAT 23Feb01: now handled by NKJE below
C              DO nk=1,NKT(nn,nb)
C                NKE(nk,nn,nb,ne)=nk
C              ENDDO !nk
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                NVJE(nn,nb,nj,ne)=1 !version one of nn,nj in elem ne
              ENDDO !nj
            ENDDO !nn
          ENDDO !ne
        ENDDO !nb
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            nb=MESH5_NB
            DO nn=1,NNT(nb)
              DO nk=1,NKT(nn,nb)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDDO !nj
        ENDDO !ne
      ENDDO !nr

      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      DO nr=1,NRT
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nj=1,NJT
            NVJP(nj,np)=1
          ENDDO !nj
        ENDDO !nonode
C GMH 14/2/97 create region
        IF(NPNODE(0,nr).GT.0) THEN
          CALL REGION_CREATE(nr,ERROR,*9999)
        ENDIF
      ENDDO !nr

      CALL EXITS('IPMESH5')
      RETURN
 9999 CALL ERRORS('IPMESH5',ERROR)
      CALL EXITS('IPMESH5')
      RETURN 1
      END


