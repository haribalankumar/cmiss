      SUBROUTINE IPMESH6_SUB(AVOID_HOST,AVOID_PT,AVOID_WT,
     '  CART_NODE,ELE_LENGTH,IBT,IDO,INP,NBJ,NEELEM,
     '  NENP,NEP,NKJE,NKJ,NODE_ORDER,NPE,
     '  NPF,NP_INTERFACE,NPNE,NPNODE,nr,nr_coronary,
     '  NRE,NXI,NVJE,NVJP,SE,SITE_NUM,SITE_BRANCH,
     '  SITE_VALUE,SITE_VECTOR,
     '  XA,XE,XIP,XP,INPUT_PATH,ERROR,*)

C#### Subroutine: IPMESH6_SUB
C###  Description:
C###    defines mesh parameters for 3D asysmetric coronary meshs
C###    within the heart.iod host mesh grown the plans
C###    defined by the fibres and sheets.
C###    This code is based around growing from the free sites produced
C###    by the existing branchs such that asymetry of orders can be
C###    acomidated as is used in Stralar ordering. The 3D branchs are
C###    confined within the plane of the fibres and sheets and step
C###    sizes for segments are adjusted to maintain an agnle between
C###    branchs and planes of less than 2 degrees. Branchs grown
C###    outside the host mesh have there free sites canceled.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),
     '  NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NEP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),NPE(0:NPM,NEM,NRM),
     '  NP_INTERFACE(0:NPM,0:3),NODE_ORDER(NP_R_M),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr,nr_coronary,NRE(NEM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM),ORDER,SEGORDER,SEGNUM,SITE_BRANCH(NORM),
     '  SITE_NUM(4,NORM)
      REAL*8 CART_NODE(3,NP_R_M),ELE_LENGTH(0:NE_R_M),
     '  SE(NSM,NBFM,NEM),SEGLENGTH,SITE_VALUE(4,NORM),
     '  SITE_VECTOR(3,NORM),XA(NAM,NJM,NEM),XIP(NIM,NPM),
     '  XP(NKM,NVM,NJM,NPM),XE(NSM,NJM)
      CHARACTER*(*) INPUT_PATH
      CHARACTER ERROR*(*)
!     Local Variables
C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
      CHARACTER*256 FILE_NAME
C DBe
      INTEGER AVOID_HOST(0:NORM,120),AVOID_PT_COUNT,
     '  COUNT,ELE_COUNT,
     '  ELE_LIST(27),HOST_COORD,I,IBEG,ICHAR,IEND,
     '  INFO,ISEED,J,K,mb,nb,nb_coronary,ne,ne1,
     '  ne_curent,ne_daughter,ne_loop,ne_nr,ni,ni1,ni2,ni3,NITB,nj,
     '  nk,nn,NOQUES,np,np1,np2,np_daughter,np_list,
     '  np_nr,ns,NUM_SITES,PARENT_ORDER,SITE,SITE_COUNT,SWEEP,
     '  SITE_IN_FILE,TOT,TOT_ORDERS,choice,XI_DIR(3)
      REAL*8 a,A_VECTOR(3),ALPHA1,ALPHA_NET,ALPHA2,
     '  AVOID_VECTOR(3),AVOID_PT(3,NORM*7),AVOID_WT(NORM*7),
     '  b,BETA1,BETA2,
     '  BOUND_COORD(3),B_VECTOR(3),c,
     '  C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,
     '  CART_BOUND_COORD(3),CONECT(0:13,11,3),DAUGHTER_VECTOR(3),DENOM,
     '  DERIV(9),DIR_COR(9),DIST,DOT,EVAL(3),
     '  EVEC(3,3),HESSIAN(3,3),MIN,MOMENTUM_W1,MOMENTUM_W2,
     '  PXI,NODE(3),NODE_STEP,POINT(3),OLDVALUE,
     '  OPTMIN,Q,R,
     '  RAD,RADIUS,STEP,STEP_PROP,T,U,V,
     '  VALUE,SUM,VECTOR(3),x,X_OPT,XI(3),XI_POSITION(3),
     '  y,Y_OPT,z,Z_OPT
      LOGICAL FILEIP,FINISH,OPTFOUND

      CALL ENTERS('IPMESH6_SUB',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ISEED=1
      ICHAR=999
C this is where sites are initialized from surface vessels
C at the moment start with one to test eventually will have to
C pass sites from ipmesh10 NPS 4/3/97


      ITYP10(nr_coronary)=ITYP10(nr)
      TOT_ORDERS=4
      BETA1=0.8d0 !boundary exponent
      BETA2=0.5d0 !network exponent
      ALPHA1=10.0d0 !bounadry multiplier
      ALPHA_NET=1.0d0 !network multiplier


      FORMAT='($,'' Enter the filename1 [current]: '',A30)'
      CDEFLT(1)=FILE00
      IF(IOTYPE.EQ.3) CDATA(1)=FILE03
      CALL GINOUT(IOTYPE,IPCHAR,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,30,IDATA,IDEFLT,0,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) THEN
        CALL STRING_TRIM(CDATA(1),IBEG,IEND)
        FILE03=CDATA(1)(IBEG:IEND)
      ENDIF
      CALL STRING_TRIM(FILE03,IBEG,IEND)
C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
C     in argument
C???DB.  F77 will truncate if necessary
      FILE_NAME=INPUT_PATH//FILE03(IBEG:IEND)//'.coronary'
      CALL OPENF(IOFILE1,
     '  'DISK',FILE_NAME,
     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C      CALL OPENF(IOFILE1,
C     '  'DISK',INPUT_PATH//FILE03(IBEG:IEND)//'.coronary',
C     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C DBe
      READ(IOFILE1,*) NUM_SITES
      DO SITE=1,NUM_SITES
        DO SITE_COUNT=1,4
          READ(IOFILE1,*) SITE_NUM(SITE_COUNT,SITE)
        ENDDO
        DO SITE_COUNT=1,3
          READ(IOFILE1,*) SITE_VALUE(SITE_COUNT,SITE)
        ENDDO
        SITE_VALUE(4,SITE)=10.0d0
        DO SITE_COUNT=1,3
          READ(IOFILE1,*) SITE_VECTOR(SITE_COUNT,SITE)
        ENDDO
        READ(IOFILE1,*) SITE_BRANCH(SITE) ! 1=RCA,2=LAD,3=CX
      ENDDO
      CALL CLOSEF(IOFILE1,ERROR,*9999)

C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
C     in argument
C???DB.  F77 will truncate if necessary
      FILE_NAME=INPUT_PATH//'rca.data'
      CALL OPENF(IOFILE1,'DISK',FILE_NAME,
     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C      CALL OPENF(IOFILE1,'DISK',INPUT_PATH//'rca.data',
C     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C DBe
      DO i=0,13 !parnet order for RCA
        READ(IOFILE1,*) (CONECT(i,j,1),j=1,11)
      ENDDO

      CALL CLOSEF(IOFILE1,ERROR,*9999)

C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
C     in argument
C???DB.  F77 will truncate if necessary
      FILE_NAME=INPUT_PATH//'lad.data'
      CALL OPENF(IOFILE1,'DISK',FILE_NAME,
     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C      CALL OPENF(IOFILE1,'DISK',INPUT_PATH//'lad.data',
C     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C DBe
      DO i=0,13 !parnet order for LAD
        READ(IOFILE1,*) (CONECT(i,j,2),j=1,11)
      ENDDO
      CALL CLOSEF(IOFILE1,ERROR,*9999)

C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
C     in argument
C???DB.  F77 will truncate if necessary
      FILE_NAME=INPUT_PATH//'cx.data'
      CALL OPENF(IOFILE1,'DISK',FILE_NAME,
     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C      CALL OPENF(IOFILE1,'DISK',INPUT_PATH//'cx.data',
C     '  'OLD','SEQUEN','FORMATTED',160,ERROR,*9999)
C DBe

      DO i=0,13 !parnet order for CX
        READ(IOFILE1,*) (CONECT(i,j,3),j=1,11)
      ENDDO
      CALL CLOSEF(IOFILE1,ERROR,*9999)

      DO I=1,3
        DO J=1,11
          SUM=0.0D0
          DO K=0,11
            SUM=SUM+CONECT(K,J,I)
          ENDDO
           DO K=0,11
           CONECT(K,J,I)=CONECT(K,J,I)/SUM
          ENDDO
          SUM=0.0D0
          DO K=12,13
            SUM=SUM+CONECT(K,J,I)
          ENDDO
          DO K=12,13
            CONECT(K,J,I)=CONECT(K,J,I)/SUM
          ENDDO
        ENDDO
      ENDDO

      FORMAT='($,'' Enter the number of orders [3]: '',I2)'
      IDEFLT(1)=3
      IF(IOTYPE.EQ.3) IDATA(1)=TOT_ORDERS
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TOT_ORDERS=IDATA(1)

      FORMAT='($,'' Enter stepsize [0.4]: '',F4.2)'
      RDEFLT(1)=0.40d0
      STEP=0.4d0
      IF(IOTYPE.EQ.3) RDATA(1)=STEP
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,20.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) STEP=RDATA(1)

      FORMAT='($,'' Enter step proportion [4.0]: '',F4.2)'
      RDEFLT(1)=4.0d0
      STEP_PROP=4.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=STEP_PROP
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,20.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) STEP_PROP=RDATA(1)

      FORMAT='($,'' Enter value for node step [3.0]: '',F4.2)'
      RDEFLT(1)=3.0d0
      NODE_STEP=3.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=NODE_STEP
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,20.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NODE_STEP=RDATA(1)

      FORMAT='($,'' Enter linear momentum wieght [2.0]: '',F4.2)'
      RDEFLT(1)=2.0d0
      MOMENTUM_W1=2.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=MOMENTUM_W1
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,500.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MOMENTUM_W1=RDATA(1)

      FORMAT='($,'' Enter quadractic momentum wieght [3.0]: '',F4.2)'
      RDEFLT(1)=3.0d0
      MOMENTUM_W2=3.0d0
      IF(IOTYPE.EQ.3) RDATA(1)=MOMENTUM_W2
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,0.0D0,500.0D0,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MOMENTUM_W2=RDATA(1)

      MOMENTUM_W2=MOMENTUM_W2/(STEP**2.0d0) !scalled by square of step

c          SITE_NUM(1,NUM_SITES)=np1 !hard coded temporarialy
c          SITE_NUM(2,NUM_SITES)=8 !order number
c          SITE_NUM(3,NUM_SITES)=0 !segment number
c          SITE_NUM(4,NUM_SITES)=3 !number of segments
c          SITE_VALUE(1,NUM_SITES)= SEGLENGTH(SITE_NUM(2,NUM_SITES))
c          SITE_VALUE(2,NUM_SITES)=STEP !current step length for site
c          SITE_VALUE(3,NUM_SITES)=0.0d0 !current length for segment
c          SITE_VALUE(4,NUM_SITES)=0.0d0 !current length for step


      NITB=NIT(NBJ(1,NEELEM(1,1)))
      nb_coronary=NBJ(1,NEELEM(1,nr_coronary))

      DO SITE=1,NUM_SITES !moves surface nodes inwards
        DO nj=1,NITB
          IF(DABS(XIP(nj,SITE_NUM(1,SITE))-1.0d0).LT.ZERO_TOL) THEN
            XIP(nj,SITE_NUM(1,SITE))=XIP(nj,SITE_NUM(1,SITE))-0.01d0
          ENDIF
          XI(nj)=XIP(nj,SITE_NUM(1,SITE))
        ENDDO
        ne=NEP(SITE_NUM(1,SITE))
        CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '    nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
        DO nj=1,NITB
          nb=NBJ(nj,ne)
          XP(1,1,nj,SITE_NUM(1,SITE))=PXI(IBT(1,1,nb),
     '      IDO(1,1,0,nb),INP(1,1,nb),nb,1,XI,XE(1,nj))
        ENDDO
        nj=NJ_LOC(NJL_GEOM,3,nr)
        np=SITE_NUM(1,SITE)
        IF(ITYP10(nr_coronary).EQ.4) THEN
          DO WHILE (XP(1,1,nj,np).GT.(2.0d0*PI))
            XP(1,1,nj,np)=XP(1,1,nj,np)-(2.0d0*PI)
          ENDDO
        ENDIF
      ENDDO




      HOST_COORD=ITYP10(nr_coronary) ! NCO(NEP(NPNODE(1,nr_coronary)))
      DO np=1,NPNODE(0,nr_coronary)  ! above replaced by stuff on left
!calculates initial shadow cartesian coords
        np1=NPNODE(np,nr_coronary)
        NODE_ORDER(np1)=0
        DO nj=1,NITB
         NODE(nj)=XP(1,1,nj,np1)
        ENDDO
        CALL COORD(HOST_COORD,1,NODE,
     '    CART_NODE(1,np1),ERROR,*9999)
        DO nj=1,NITB
          AVOID_PT(nj,np1)=CART_NODE(nj,np1)
        ENDDO
      ENDDO

      AVOID_PT_COUNT=NPNODE(NPNODE(0,nr_coronary),nr_coronary)

      DO ne=1,NEELEM(0,nr_coronary)
!calculates initial shadow cartesian coords
        ne1=NEELEM(ne,nr_coronary)
        nb=NBJ(1,ne1)
        np1=NPNE(1,nb,ne1)
        np2=NPNE(2,nb,ne1)
        ELE_LENGTH(ne1)=0.0d0
        DO nj=1,NITB
          ELE_LENGTH(ne1)= ELE_LENGTH(ne1)
     '      +(CART_NODE(nj,np1)-CART_NODE(nj,np2))**2.0d0
        ENDDO
        ELE_LENGTH(ne1)=ELE_LENGTH(ne1)**0.5d0
      ENDDO

      DO ne=1,NEELEM(0,nr)
       ne1=NEELEM(ne,nr)
       DO np1=0,NPE(0,ne1,nr_coronary)
         AVOID_HOST(np1,ne1)=NPE(np1,ne1,nr_coronary)
       ENDDO
      ENDDO

      DO np=1,NPNODE(0,nr_coronary)
        np1=NPNODE(np,nr_coronary)
C KAT 2002-02-01: avoiding uninitialized components of NENP
        SUM=0.0d0
        DO i=1,NENP(np1,0,nr_coronary)
          SUM=SUM+ELE_LENGTH(NENP(np1,i,nr_coronary))
        ENDDO !i
        AVOID_WT(np1)=SUM/2.0d0
C        AVOID_WT(np1)=(ELE_LENGTH(NENP(np1,1,nr_coronary))
C     '    +ELE_LENGTH(NENP(np1,2,nr_coronary)))/2.0d0
      ENDDO !np

      CHOICE=1
      FORMAT='($,'' Enter 1 to read avoidance file '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=CHOICE
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CHOICE=IDATA(1)

      IF (CHOICE.EQ.1)THEN
C DBs 22/3/01: g77 doesn't like unknown length string being concatenated
C     in argument
C???DB.  F77 will truncate if necessary
        FILE_NAME=INPUT_PATH//'avoid.data'
        CALL OPENF(IOFILE1,'DISK',FILE_NAME,'OLD',
     '    'SEQUEN','FORMATTED',160,ERROR,*9999)
C        CALL OPENF(IOFILE1,'DISK',INPUT_PATH//'avoid.data','OLD',
C     '    'SEQUEN','FORMATTED',160,ERROR,*9999)
C DBe
        READ(IOFILE1,*) AVOID_PT_COUNT
        DO I=1,AVOID_PT_COUNT
          READ(IOFILE1,*) (AVOID_PT(j,i),j=1,3)
          READ(IOFILE1,*) AVOID_WT(i)
        ENDDO
        DO ne=1,NEELEM(0,nr)
          ne1=NEELEM(ne,nr)
          READ(IOFILE1,*) AVOID_HOST(0,ne1)
          READ(IOFILE1,*) (AVOID_HOST(np1,ne1),np1=1,AVOID_HOST(0,ne1))
        ENDDO
        CALL CLOSEF(IOFILE1,ERROR,*9999)
      ENDIF
      SWEEP=0



      DO ORDER=11,(11-TOT_ORDERS),-1
        FINISH=.FALSE.
        IF (ORDER.EQ.6)THEN
          STEP=2.5d0
c          AVOID_PT_COUNT=AVOID_PT_COUNT+1
c          AVOID_WT(AVOID_PT_COUNT)=0.0d0
        ENDIF
        DO WHILE (.NOT.FINISH)
C checks no new sites have been created of order equal
C to the generation
          FINISH=.TRUE.
          SITE=0
          SWEEP=0
          DO WHILE (SITE.LT.NUM_SITES)
            SITE=SITE+1
            IF(SITE_NUM(2,SITE).EQ.ORDER) THEN
C checkes no more sites of order equal to the generation
              FINISH=.FALSE.
              SWEEP=SWEEP+1
              ne=NEP(SITE_NUM(1,SITE))
              ELE_COUNT=0

              IF (ORDER.GT.8)THEN
                DO ni3=-3,3,3
                  DO ni2=-2,2,2
                    DO ni1=-1,1,1
                      ELE_LIST(ELE_COUNT+1)=
     '                  NXI(ni3,1,NXI(ni2,1,NXI(ni1,1,ne)))
                      IF (ELE_LIST(ELE_COUNT+1).NE.0) THEN
                        ELE_COUNT=ELE_COUNT+1
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ELSE
                DO ni=1,NITB
                  IF (XIP(ni,SITE_NUM(1,SITE)).LT.0.5d0)THEN
                    XI_DIR(ni)=-1
                  ELSE
                    XI_DIR(ni)=1
                  ENDIF
                ENDDO
                DO ni3=0,3,3
                  DO ni2=0,2,2
                    DO ni1=0,1,1
                      ELE_LIST(ELE_COUNT+1)=
     '                  NXI((XI_DIR(3)*ni3),1,
     '                  NXI((XI_DIR(2)*ni2),1,
     '                  NXI((XI_DIR(1)*ni1),1,ne)))
                      IF (ELE_LIST(ELE_COUNT+1).NE.0) THEN
                        ELE_COUNT=ELE_COUNT+1
                      ENDIF
                    ENDDO
                  ENDDO
                ENDDO
              ENDIF
C creat element list with which to search for avoidance
C creates smaller list for lower orders


              DERIV(1)=0.0d0 ! wrt x
              DERIV(2)=0.0d0 ! wrt y
              DERIV(3)=0.0d0 ! wrt z
              DERIV(4)=0.0d0 ! wrt x^2
              DERIV(5)=0.0d0 ! wrt y^2
              DERIV(6)=0.0d0 ! wrt z^2
              DERIV(7)=0.0d0 ! wrt xy
              DERIV(8)=0.0d0 ! wrt xz
              DERIV(9)=0.0d0 ! wrt yz

              DO ni1=1,NITB
                DO ni2=-1,1,2
                  COUNT=0
                  ne_curent=ne
                  DO WHILE ((NXI((ni2*ni1),1,ne_curent).NE.0)
     '              .AND.(COUNT.LE.2))
                    COUNT=COUNT+1
                    ne_curent=NXI((ni2*ni1),1,ne_curent)
                  ENDDO
C                  BOUND(ni1*ni2)=ne_curent
                  IF (COUNT.NE.3)THEN
                    np=SITE_NUM(1,SITE)
                    DIST=0.0d0
                    DO nj=1,NITB
                      IF (nj.NE.ni1)THEN
                        XI(nj)=XIP(nj,np)
C cordinates of point on boundary
                      ELSE
                        XI(nj)=MAX(ni2,0)
                      ENDIF
                    ENDDO
                    CALL XPXE(NBJ(1,ne_curent),
     '                NKJE(1,1,1,ne_curent),NPF(1,1),
     '                NPNE(1,1,ne_curent),
     '                nr,NVJE(1,1,1,ne_curent),
     '                SE(1,1,ne_curent),XA(1,1,ne_curent),XE,XP,ERROR,
     '                *9999)
                      nb=NBJ(nj,ne_curent)
                    DO nj=1,NITB
                      nb=NBJ(nj,ne_curent)
                      BOUND_COORD(nj)=PXI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                  INP(1,1,nb),nb,1,XI,XE(1,nj))
                    ENDDO
                    CALL COORD(ITYP10(NRE(ne_curent)),1,BOUND_COORD,
     '                CART_BOUND_COORD,ERROR,*9999)
                    DO nj=1,NITB
                      DIST=DIST+((CART_NODE(nj,np)-
     '                  CART_BOUND_COORD(nj))**2)
                      VECTOR(nj)=CART_NODE(nj,np)-CART_BOUND_COORD(nj)
                    ENDDO
                    DIST=DIST**0.5d0
                    CALL NORMALISE(NITB,VECTOR,ERROR,*9999)

c                    IF (DIST.LT.SITE_VALUE(2,np))THEN
c                      DIST=DIST**2.0d0/SITE_VALUE(2,np)
c                    ENDIF
C the above increasing the boundary avoidance when close tp the boundary
C put in after testing


                    IF (DIST.GT.ZERO_TOL)THEN


                      a=CART_BOUND_COORD(1)
                      b=CART_BOUND_COORD(2)
                      c=CART_BOUND_COORD(3)
                      x=CART_NODE(1,np)
                      y=CART_NODE(2,np)
                      z=CART_NODE(3,np)
                      DENOM=(VECTOR(1)*(x-a))+(VECTOR(2)*(y-b))+
     '                  (VECTOR(3)*(z-c))

                      DERIV(1)=DERIV(1)+
     '                  (ALPHA1*(-BETA1)*VECTOR(1))/(DENOM**(BETA1+1))
                      DERIV(2)=DERIV(2)+
     '                  (ALPHA1*(-BETA1)*VECTOR(2))/(DENOM**(BETA1+1))
                      DERIV(3)=DERIV(3)+
     '                  (ALPHA1*(-BETA1)*VECTOR(3))/(DENOM**(BETA1+1))

                   ENDIF
                  ENDIF ! COUNT.NE.2
                ENDDO !ni2
              ENDDO ! ni1


C calculate boundary avoidance vector ABOVE DONE

              np=SITE_NUM(1,SITE)
              DO ne=1,ELE_COUNT
                ne_curent=ELE_LIST(ne)
                DO np_list=1,AVOID_HOST(0,ne_curent)
                  np1=AVOID_HOST(np_list,ne_curent)
c                  IF (NODE_ORDER(np1).LE.SITE_NUM(2,SITE))THEN
                    DIST=0.0d0
                    DO nj=1,NITB
                      DIST=DIST+((CART_NODE(nj,np)-
     '                  AVOID_PT(nj,np1))**2.0d0)
                    ENDDO !nj
                    DIST=DIST**0.5d0
                    IF (DIST.GT.ZERO_TOL)THEN
                      ALPHA2=ALPHA_NET*AVOID_WT(np1)
                      a=AVOID_PT(1,np1)
                      b=AVOID_PT(2,np1)
                      c=AVOID_PT(3,np1)
                      x=CART_NODE(1,np)
                      y=CART_NODE(2,np)
                      z=CART_NODE(3,np)

                      DERIV(1)=DERIV(1)+
     '                  (ALPHA2*(-BETA2)*(DIST**(-BETA2-2.0d0))*
     '                  (x-a))
                      DERIV(2)=DERIV(2)+
     '                  ((ALPHA2*(-BETA2*DIST**(-BETA2-2.0d0)))*
     '                  (y-b))
                      DERIV(3)=DERIV(3)+
     '                  ((ALPHA2*(-BETA2*DIST**(-BETA2-2.0d0)))*
     '                  (z-c))
                      IF (DIST.GT.(0.3d0*SITE_VALUE(2,SITE)))THEN
                        DERIV(4)=DERIV(4)+
     '                    (ALPHA2*(-BETA2*DIST**(-BETA2-2.0d0)))+
     '                    (ALPHA2*((x-a)**2.0d0)*BETA2*(BETA2+2.0d0)*
     '                    (DIST**(-BETA2-4.0d0)))
                        DERIV(5)=DERIV(5)+
     '                    (ALPHA2*(-BETA2*DIST**(-BETA2-2.0d0)))+
     '                    (ALPHA2*((y-b)**2.0d0)*BETA2*(BETA2+2.0d0)*
     '                    (DIST**(-BETA2-4.0d0)))
                        DERIV(6)=DERIV(6)+
     '                    (ALPHA2*(-BETA2*DIST**(-BETA2-2.0d0)))+
     '                    (ALPHA2*((z-c)**2.0d0)*BETA2*(BETA2+2.0d0)*
     '                    (DIST**(-BETA2-4.0d0)))
                        DERIV(7)=DERIV(7)+
     '                    (ALPHA2*BETA2*(BETA2+2.0d0)*(x-a)*(y-b)
     '                    *(DIST**(-BETA2-4.0d0)))
                        DERIV(8)=DERIV(8)+
     '                    (ALPHA2*BETA2*(BETA2+2.0d0)*(x-a)*(z-c)
     '                    *(DIST**(-BETA2-4.0d0)))
                        DERIV(9)=DERIV(9)+
     '                    (ALPHA2*BETA2*(BETA2+2.0d0)*(y-b)*(z-c)
     '                    *(DIST**(-BETA2-4.0d0)))
                      ENDIF
                    ENDIF
c                  ENDIF
                ENDDO ! np_list
              ENDDO ! ne

C now need to calculate netword avoidance vector ABOVE

              RAD=SITE_VALUE(2,SITE)/STEP_PROP

              DO nj=1,NJT
                AVOID_VECTOR(nj)=-DERIV(nj)
              ENDDO
              CALL NORMALISE(NITB,AVOID_VECTOR,ERROR,*9999)
              DO nj=1,NJT
                AVOID_VECTOR(nj)=AVOID_VECTOR(nj)*RAD
              ENDDO


              CALL NORMALISE(NITB,SITE_VECTOR(1,SITE),ERROR,*9999)


C This sectrion of code adjusts the coefficents of the talyor series
C expansion around the point to basis it in the flow direction

              A_VECTOR(1)=-SITE_VECTOR(2,SITE)
              A_VECTOR(2)=SITE_VECTOR(1,SITE)
              A_VECTOR(3)=0

              B_VECTOR(1)=-(SITE_VECTOR(1,SITE)*SITE_VECTOR(3,SITE))
              B_VECTOR(2)=-(SITE_VECTOR(2,SITE)*SITE_VECTOR(3,SITE))
              B_VECTOR(3)=(SITE_VECTOR(1,SITE)**2.0D0)+
     '          (SITE_VECTOR(2,SITE)**2.0D0)

              CALL NORMALISE(NITB,A_VECTOR(1),ERROR,*9999)
              CALL NORMALISE(NITB,B_VECTOR(1),ERROR,*9999)

              Q=A_VECTOR(1)
              R=A_VECTOR(2)
              T=B_VECTOR(1)
              U=B_VECTOR(2)
              V=B_VECTOR(3)

              IF ((SITE_VALUE(3,SITE).LE.ZERO_TOL).AND.
     '          (SITE_NUM(3,SITE).EQ.0))THEN !branching

                DIR_COR(1)=MOMENTUM_W1*(-SITE_VECTOR(1,SITE))
                DIR_COR(2)=MOMENTUM_W1*(-SITE_VECTOR(2,SITE))
                DIR_COR(3)=MOMENTUM_W1*(-SITE_VECTOR(3,SITE))
                DIR_COR(4)=-MOMENTUM_W2*((Q**2)+(T**2)) !X^2
                DIR_COR(5)=-MOMENTUM_W2*((R**2)+(U**2)) !Y^2
                DIR_COR(6)=-MOMENTUM_W2*(V**2) !Z^2
                DIR_COR(7)=-MOMENTUM_W2*(2.0d0*((Q*R)+(T*U))) !xy
                DIR_COR(8)=-MOMENTUM_W2*(2.0d0*T*V) !xz
                DIR_COR(9)=-MOMENTUM_W2*(2.0d0*U*V) !yz

c look for next branch and calculates penality if too close to the
C parent branch
                IF (NENP(SITE_NUM(1,SITE),0,nr_coronary).GT.1)THEN
                  ne_daughter=NENP(SITE_NUM(1,SITE),1,nr_coronary)
                  IF (NENP(SITE_NUM(1,SITE),2,nr_coronary)
     '              .GT.ne_daughter)THEN
                    ne_daughter=NENP(SITE_NUM(1,SITE),2,nr_coronary)
                  ENDIF
                  mb=NBJ(1,ne_daughter)
                  IF (NPNE(1,mb,ne_daughter).EQ.SITE_NUM(1,SITE))THEN
                    np_daughter=NPNE(2,mb,ne_daughter)
                  ELSE
                    np_daughter=NPNE(1,mb,ne_daughter)
                  ENDIF
                  DO nj=1,NITB
                    DAUGHTER_VECTOR(nj)=CART_NODE(nj,np_daughter)-
     '                CART_NODE(nj,SITE_NUM(1,SITE))
                  ENDDO
                  A_VECTOR(1)=-DAUGHTER_VECTOR(2)
                  A_VECTOR(2)=DAUGHTER_VECTOR(1)
                  A_VECTOR(3)=0

                  B_VECTOR(1)=-(DAUGHTER_VECTOR(1)*
     '              DAUGHTER_VECTOR(3))
                  B_VECTOR(2)=-(DAUGHTER_VECTOR(2)*
     '              DAUGHTER_VECTOR(3))
                  B_VECTOR(3)=(DAUGHTER_VECTOR(1)**2.0D0)+
     '              (DAUGHTER_VECTOR(2)**2.0D0)

                  CALL NORMALISE(NITB,A_VECTOR,ERROR,*9999)
                  CALL NORMALISE(NITB,B_VECTOR,ERROR,*9999)

                  Q=A_VECTOR(1)
                  R=A_VECTOR(2)
                  T=B_VECTOR(1)
                  U=B_VECTOR(2)
                  V=B_VECTOR(3)
                  DIR_COR(4)=DIR_COR(4)-MOMENTUM_W2*((Q**2)+(T**2)) !X^2
                  DIR_COR(5)=DIR_COR(5)-MOMENTUM_W2*((R**2)+(U**2)) !Y^2
                  DIR_COR(6)=DIR_COR(6)-MOMENTUM_W2*(V**2) !Z^2
                  DIR_COR(7)=DIR_COR(7)-MOMENTUM_W2*
     '              (2.0d0*((Q*R)+(T*U)))
                  DIR_COR(8)=DIR_COR(8)-MOMENTUM_W2*(2.0d0*T*V) !xz
                  DIR_COR(9)=DIR_COR(9)-MOMENTUM_W2*(2.0d0*U*V) !yz

                ENDIF
              ELSE
                DIR_COR(1)=MOMENTUM_W1*(-SITE_VECTOR(1,SITE))
                DIR_COR(2)=MOMENTUM_W1*(-SITE_VECTOR(2,SITE))
                DIR_COR(3)=MOMENTUM_W1*(-SITE_VECTOR(3,SITE))
                DIR_COR(4)=MOMENTUM_W2*((Q**2)+(T**2)) !X^2
                DIR_COR(5)=MOMENTUM_W2*((R**2)+(U**2)) !Y^2
                DIR_COR(6)=MOMENTUM_W2*(V**2) !Z^2
                DIR_COR(7)=MOMENTUM_W2*(2.0d0*((Q*R)+(T*U))) !xy
                DIR_COR(8)=MOMENTUM_W2*(2.0d0*T*V) !xz
                DIR_COR(9)=MOMENTUM_W2*(2.0d0*U*V) !yz
              ENDIF
              DO I=1,9
                DERIV(I)=DERIV(I)+DIR_COR(I)
              ENDDO
              HESSIAN(1,1)=DERIV(4)
              HESSIAN(1,2)=DERIV(7)
              HESSIAN(1,3)=DERIV(8)
              HESSIAN(2,1)=DERIV(7)
              HESSIAN(2,2)=DERIV(5)
              HESSIAN(2,3)=DERIV(9)
              HESSIAN(3,1)=DERIV(8)
              HESSIAN(3,2)=DERIV(9)
              HESSIAN(3,3)=DERIV(6)

              CALL EIGEN1(3,3,HESSIAN,EVAL,EVEC,ERROR,*9999)

C clac egen value and vectors here FROM HESSIAN


              C1=0.5d0*EVAL(1)
              C2=0.5d0*EVAL(2)
              C3=0.5d0*EVAL(3)
              C4=DERIV(1)*EVEC(1,1)+DERIV(2)*EVEC(2,1)+
     '          DERIV(3)*EVEC(3,1)
              C5=DERIV(1)*EVEC(1,2)+DERIV(2)*EVEC(2,2)+
     '          DERIV(3)*EVEC(3,2)
              C6=DERIV(1)*EVEC(1,3)+DERIV(2)*EVEC(2,3)+
     '          DERIV(3)*EVEC(3,3)
              C7=2*C1-2*C2
              C8=2*C1-2*C3
              TOT=1000
              COUNT=0

              X_OPT=EVEC(1,1)*AVOID_VECTOR(1)+
     '          EVEC(2,1)*AVOID_VECTOR(2)+
     '          EVEC(3,1)*AVOID_VECTOR(3)
              Y_OPT=EVEC(1,2)*AVOID_VECTOR(1)+
     '          EVEC(2,2)*AVOID_VECTOR(2)+
     '          EVEC(3,2)*AVOID_VECTOR(3)
              Z_OPT=EVEC(1,3)*AVOID_VECTOR(1)+
     '          EVEC(2,3)*AVOID_VECTOR(2)+
     '          EVEC(3,3)*AVOID_VECTOR(3)

              OPTMIN=(C1*(X_OPT**2))+(C2*(Y_OPT**2))+
     '          (C3*(Z_OPT**2))+(C4*X_OPT)+(C5*Y_OPT)+(C6*Z_OPT)
              OPTFOUND=.FALSE.
              X=-(DBLE(TOT+1))/(DBLE(TOT))
              C9=((C7*X)+C4)**2
              C10=((C8*X)+C4)**2
              OLDVALUE=((X**2)*C9*C10)+(((C5*X)**2)*C10)+
     '          (((C6*X)**2)*C9)-
     '          ((C9*C10)*(RAD**2))
              DO I=-TOT,TOT
                X=(DBLE(I))/(DBLE(TOT))
                C9=((C7*X)+C4)**2
                C10=((C8*X)+C4)**2
                VALUE=((X**2)*C9*C10)+(((C5*X)**2)*C10)+(((C6*X)**2)
     '            *C9)-
     '            ((C9*C10)*(RAD**2))
                IF ((VALUE*OLDVALUE).LE.0.0d0)THEN
                  Y=(C5*X)/((C7*X)+C4)
                  Z=(C6*X)/((C8*X)+C4)
                  MIN=(C1*(X**2))+(C2*(Y**2))+(C3*(Z**2))+(C4*X)+
     '              (C5*Y)+(C6*Z)
                  IF (MIN.LT.OPTMIN)THEN
                    X_OPT=X
                    Y_OPT=Y
                    Z_OPT=Z
                    OPTMIN=MIN
                    OPTFOUND=.TRUE.
                  ENDIF
                ENDIF
                OLDVALUE=VALUE
              ENDDO
C NOW NEED TO TRANSFORM CORDINATES BACK

              IF (OPTFOUND)THEN
                DO nj=1,NJT
                  AVOID_VECTOR(nj)=EVEC(nj,1)*X_OPT+
     '              EVEC(nj,2)*Y_OPT+EVEC(nj,3)*Z_OPT
                ENDDO
              ENDIF

              CALL NORMALISE(NITB,AVOID_VECTOR,ERROR,*9999)
              CALL NORMALISE(NITB,SITE_VECTOR(1,SITE),ERROR,*9999)

              IF (ORDER.NE.6)THEN
                DOT=AVOID_VECTOR(1)*SITE_VECTOR(1,SITE)+
     '            AVOID_VECTOR(2)*SITE_VECTOR(2,SITE)+
     '            AVOID_VECTOR(3)*SITE_VECTOR(3,SITE)
                IF ((SITE_VALUE(3,SITE).LE.ZERO_TOL).AND.
     '            (SITE_NUM(3,SITE).EQ.0))THEN !branching
                  DO nj=1,NITB
                    AVOID_VECTOR(nj)=AVOID_VECTOR(nj)+
     '                (((((1.0d0-DOT)/2.0D0))**1.0d0)*
     '                SITE_VECTOR(nj,SITE))
                  ENDDO
                ENDIF
              ENDIF

              CALL NORMALISE(NITB,AVOID_VECTOR,ERROR,*9999)

c              IF (ORDER.GT.6)THEN
                AVOID_PT_COUNT=AVOID_PT_COUNT+1
c              ENDIF
                DO nj=1,NJT
                  AVOID_PT(nj,AVOID_PT_COUNT)=
     '            CART_NODE(nj,SITE_NUM(1,SITE))+
     '              (SITE_VALUE(2,SITE)*AVOID_VECTOR(nj))
                ENDDO
              ne=0
              ne_loop=0
              DO WHILE((ne.EQ.0).AND.(ne_loop.LT.NEELEM(0,nr)))
                ne_loop=ne_loop+1
                ne_curent=NEELEM(ne_loop,nr)
                DO nj=1,NJT
                  POINT(nj)=AVOID_PT(nj,AVOID_PT_COUNT)
                ENDDO
                CALL XPXE(NBJ(1,ne_loop),
     '            NKJE(1,1,1,ne_loop),NPF(1,1),
     '            NPNE(1,1,ne_loop),
     '            nr,NVJE(1,1,1,ne_loop),
     '            SE(1,1,ne_loop),XA(1,1,ne_loop),XE,XP,ERROR,*9999)
                DO nj=1,NITB !initialising XI
                  XI_POSITION(nj)=0.5d0
                ENDDO
               CALL DEXI_POINT(IBT,IDO,INP,ne,NBJ,
     '            ne_curent,NITB,NRE(ne_curent),10.d0,XE,XI_POSITION,
     '            XI_POSITION,POINT,.TRUE.,ERROR,*9999)
              ENDDO


              IF (ne.NE.0)THEN !if inside host mesh
                AVOID_HOST(0,ne)=AVOID_HOST(0,ne)+1
                AVOID_HOST(AVOID_HOST(0,ne),ne)=AVOID_PT_COUNT
                AVOID_WT(AVOID_PT_COUNT)=(SITE_VALUE(2,SITE)+STEP)/2.0d0
                SITE_VALUE(3,SITE)=SITE_VALUE(3,SITE)+
     '             SITE_VALUE(2,SITE)
                SITE_VALUE(2,SITE)=STEP
                SITE_VALUE(4,SITE)=SITE_VALUE(4,SITE)+SITE_VALUE(2,SITE)
                DO nj=1,NITB
                  SITE_VECTOR(nj,SITE)=AVOID_VECTOR(nj)
                ENDDO
                IF (SITE_VALUE(3,SITE).GE.SITE_VALUE(1,SITE))THEN
C segment end
                  SITE_NUM(3,SITE)= SITE_NUM(3,SITE)+1
                  SITE_VALUE(3,SITE)=0.0d0
                  SITE_VALUE(4,SITE)=NODE_STEP-ZERO_TOL
                  NUM_SITES=NUM_SITES+1
                  SITE_NUM(1,NUM_SITES)=np
                  SITE_NUM(3,NUM_SITES)=0 !segment number
                  SITE_VALUE(2,NUM_SITES)=STEP
                  SITE_VALUE(3,NUM_SITES)=0.0d0 !current length for seg
                  SITE_VALUE(4,NUM_SITES)=10.0d0+NODE_STEP!node length for seg

                  DO nj=1,NITB
                    SITE_VECTOR(nj,NUM_SITES)=AVOID_VECTOR(nj)
                  ENDDO
                  SITE_BRANCH(NUM_SITES)=SITE_BRANCH(SITE)
                  IF (SITE_NUM(3,SITE).LT.SITE_NUM(4,SITE))THEN
C not vessel end
C NEW SITE
                    SITE_NUM(2,NUM_SITES)=SEGORDER(SITE_NUM(2,SITE),
     '                .FALSE.,
     '                CONECT(0,SITE_NUM(2,SITE),SITE_BRANCH(SITE))
     '                ,ISEED)
                    SITE_NUM(4,NUM_SITES)
     '                =SEGNUM(SITE_NUM(2,NUM_SITES),SITE_BRANCH(SITE))
C number of segments
                    SITE_VALUE(1,NUM_SITES)=
     '                SEGLENGTH(SITE_NUM(2,NUM_SITES),SITE_BRANCH(SITE))
C total length for segment
                  ELSE
C VESSEL END
                    PARENT_ORDER=SITE_NUM(2,SITE)
                    SITE_NUM(2,SITE)=
     '                SEGORDER(PARENT_ORDER,.TRUE.,
     '                CONECT(0,PARENT_ORDER,SITE_BRANCH(SITE)),ISEED)
                    SITE_NUM(4,SITE)=
     '                SEGNUM(SITE_NUM(2,SITE),SITE_BRANCH(SITE))
!number of segments
                    SITE_VALUE(1,SITE)=
     '                SEGLENGTH(SITE_NUM(2,SITE),SITE_BRANCH(SITE))
!total length for segment
                    SITE_NUM(3,SITE)=0 !segment number
                    SITE_NUM(2,NUM_SITES)=
     '                SEGORDER(PARENT_ORDER,.TRUE.,
     '                CONECT(0,PARENT_ORDER,SITE_BRANCH(SITE)),ISEED)
                    SITE_NUM(4,NUM_SITES)=
     '                SEGNUM(SITE_NUM(2,NUM_SITES),SITE_BRANCH(SITE))
!number of segments
                    SITE_VALUE(1,NUM_SITES)=
     '                SEGLENGTH(SITE_NUM(2,NUM_SITES),SITE_BRANCH(SITE))
!total length for segment
                  ENDIF
                ELSE
                  IF (SITE_VALUE(3,SITE)+
     '              (1.50d0*SITE_VALUE(2,SITE)).GE.
     '              SITE_VALUE(1,SITE))THEN
C adjusts step size with comming close to the end of the segment
C makes sure its not too close
                    SITE_VALUE(2,SITE)=SITE_VALUE(1,SITE)-
     '                SITE_VALUE(3,SITE)+0.001d0
                  ENDIF
                ENDIF

C if site_value(4,site) gt step or at bifurcation then new node and ele
C ment
                IF (SITE_VALUE(4,SITE).GT.NODE_STEP)THEN
                  SITE_VALUE(4,SITE)=0.0d0
c new node and element
                  np_nr=NPNODE(0,nr_coronary)+1
C np_nr is the node number in region nr_coronary
                  NPNODE(np_nr,nr_coronary)=NPT(nr_coronary)+1
C                   node#s in region nr
                  np=NPNODE(np_nr,nr_coronary)
                  DO nj=1,NJT
                    CART_NODE(nj,np)=AVOID_PT(nj,AVOID_PT_COUNT)
                  ENDDO
                  CALL COORD(1,HOST_COORD,CART_NODE(1,np),
     '              NODE,ERROR,*9999)
                  DO nj=1,NITB
                    XP(1,1,nj,np)=NODE(nj)
                  ENDDO
                  DO nj=1,NITB
                    XIP(nj,np)=XI_POSITION(nj)
                  ENDDO
                  NPNODE(np,0)=NPNODE(0,0)+1
                  XP(1,1,nj_loc(njl_fiel,1,nr_coronary),np)=
     '              RADIUS(SITE_NUM(2,SITE),SITE_BRANCH(SITE))
C node radius
                  DO nj=1,NJT+NJ_LOC(NJL_FIEL,0,nr_coronary) !geom & radius
                    NKJ(nj,np)=1 !no derivatives
                  ENDDO !nj
                  DO nj=1,NJT+NJ_LOC(NJL_FIEL,0,nr_coronary) !geom+radius
                    NVJP(nj,np)=1
                  ENDDO !nj
                  CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)

                  NEP(np)=ne
                  NP_INTERFACE(np,0)=1
                  NP_INTERFACE(np,1)=nr_coronary
                  IF (NE.NE.0)THEN
                    NPE(0,ne,nr_coronary)=NPE(0,ne,nr_coronary)+1
                    NPE(NPE(0,ne,nr_coronary),ne,nr_coronary)=np
                  ENDIF
                  NPNODE(0,0)=NPNODE(0,0)+1
                  NPNODE(0,nr_coronary)= NPNODE(0,nr_coronary)+1
                  NPT(nr_coronary)= NPT(nr_coronary)+1
                  NPT(0)=NPT(0)+1
                  NODE_ORDER(np)=SITE_NUM(2,SITE)
C need to creat new node with vector
C   highest element# in region nr
                  ITYP10(nr_coronary)=ITYP10(nr)
                  ne_nr=NEELEM(0,nr_coronary)+1
C ne_nr is the element number in region nr_coronary
                  NEELEM(ne_nr,nr_coronary)=
     '              NEELEM(NEELEM(0,nr_coronary),nr_coronary)+1
C   element#s in region nr
                  ne=NEELEM(ne_nr,nr_coronary)
                  NRE(ne)=nr_coronary !region#
                  DO nj=1,
     '              NJ_LOC(NJL_GEOM,0,nr_coronary)+!RGB replaced NJE(ne)
     '              NJ_LOC(NJL_FIEL,0,nr_coronary)
C geom & radius
                    NBJ(nj,ne)=nb_coronary
                  ENDDO !nj
                  nb=NBJ(1,ne)
                  DO ns=1,NST(nb)+NAT(nb)
                    SE(ns,nb,ne)=1.0d0 !scaling factor is 1
                  ENDDO !ns
                  DO nn=1,NNT(nb)
C KAT 23Feb01: now handled by NKJE below
C                    DO nk=1,NKT(nn,nb)
C                      NKE(nk,nn,nb,ne)=nk
C                    ENDDO !nk
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr_coronary)+ !RGB again
     '                NJ_LOC(NJL_FIEL,0,nr_coronary)
C geom & radius
                      NVJE(nn,nb,nj,ne)=1 !version 1 of nn,nj in elem ne
                      DO nk=1,NKT(nn,nb)
                        NKJE(nk,nn,nj,ne)=nk
                      ENDDO !nk
                    ENDDO !nj
                  ENDDO !nn
                  nb=NBJ(1,ne)
                  np1=SITE_NUM(1,SITE)
                  np2=np
                  nj=NJ_LOC(NJL_GEOM,3,nr)
                  NPNE(1,nb,ne)=np1
                  NPNE(2,nb,ne)=np2
c-------------
                  IF(ITYP10(nr_coronary).EQ.4) THEN
                    DO WHILE (XP(1,1,nj,np).GT.(2.0d0*PI))
                      XP(1,1,nj,np)=XP(1,1,nj,np)-(2.0d0*PI)
                    ENDDO
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
                  ENDIF
c------

                  NEELEM(0,nr_coronary)=NEELEM(0,nr_coronary)+1
C #elements in region nr
                  NET(nr_coronary)= NET(nr_coronary)+1
                  NEELEM(0,0)=NEELEM(0,0)+1
                  NET(0)=NET(0)+1
                  CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)
C need to creat new element with node
                  SITE_NUM(1,SITE)=np
                ELSE
                  np=SITE_NUM(1,SITE)
                  IF ((SITE_VALUE(4,SITE)+SITE_VALUE(2,SITE)).GT.
     '              NODE_STEP)THEN
C this is a node at a bifurcation or end of node step and is there for
C perminent
                    NPE(0,ne,nr_coronary)=NPE(0,ne,nr_coronary)+1
                    NPE(NPE(0,ne,nr_coronary),ne,nr_coronary)=np
                  ENDIF
C move old node
                  np=SITE_NUM(1,SITE)
                  DO nj=1,NJT
                    CART_NODE(nj,np)=AVOID_PT(nj,AVOID_PT_COUNT)
                  ENDDO
                  CALL COORD(1,HOST_COORD,CART_NODE(1,np),
     '              NODE,ERROR,*9999)
                  DO nj=1,NITB
                    XP(1,1,nj,np)=NODE(nj)
                  ENDDO
                  DO nj=1,NITB
                    XIP(nj,np)=XI_POSITION(nj)
                  ENDDO
                  NEP(np)=ne
                  ne=NENP(np,1,nr_coronary)
                  nb=NBJ(1,ne)
                  np1=NPNE(1,nb,ne)
                  np2=NPNE(2,nb,ne)
                  nj=NJ_LOC(NJL_GEOM,3,nr)
                  NPNE(1,nb,ne)=np1
                  NPNE(2,nb,ne)=np2
c---------------
                  IF(ITYP10(nr_coronary).EQ.4) THEN
                    DO WHILE (XP(1,1,nj,np).GT.(2.0d0*PI))
                      XP(1,1,nj,np)=XP(1,1,nj,np)-(2.0d0*PI)
                    ENDDO
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
                  ENDIF
c--------------------------
                ENDIF
              ELSE
                NUM_SITES=NUM_SITES-1
                DO np1=SITE,NUM_SITES
                  SITE_NUM(1,np1)=SITE_NUM(1,np1+1)
                  SITE_NUM(2,np1)=SITE_NUM(2,np1+1)
                  SITE_NUM(3,np1)=SITE_NUM(3,np1+1)
                  SITE_NUM(4,np1)=SITE_NUM(4,np1+1)
                  SITE_VALUE(1,np1)=SITE_VALUE(1,np1+1)
                  SITE_VALUE(2,np1)=SITE_VALUE(2,np1+1)
                  SITE_VALUE(3,np1)=SITE_VALUE(3,np1+1)
                  SITE_VECTOR(1,np1)=SITE_VECTOR(1,np1+1)
                  SITE_VECTOR(2,np1)=SITE_VECTOR(2,np1+1)
                  SITE_VECTOR(3,np1)=SITE_VECTOR(3,np1+1)
                  SITE_BRANCH(np1)= SITE_BRANCH(np1+1)
                ENDDO
              ENDIF
              IF (SITE_NUM(2,NUM_SITES).LT.6)THEN
                NUM_SITES=NUM_SITES-1
              ENDIF
              IF (SITE_NUM(2,SITE).LT.6)THEN
                DO np1=SITE,NUM_SITES
                  SITE_NUM(1,np1)=SITE_NUM(1,np1+1)
                  SITE_NUM(2,np1)=SITE_NUM(2,np1+1)
                  SITE_NUM(3,np1)=SITE_NUM(3,np1+1)
                  SITE_NUM(4,np1)=SITE_NUM(4,np1+1)
                  SITE_VALUE(1,np1)=SITE_VALUE(1,np1+1)
                  SITE_VALUE(2,np1)=SITE_VALUE(2,np1+1)
                  SITE_VALUE(3,np1)=SITE_VALUE(3,np1+1)
                  SITE_VECTOR(1,np1)=SITE_VECTOR(1,np1+1)
                  SITE_VECTOR(2,np1)=SITE_VECTOR(2,np1+1)
                  SITE_VECTOR(3,np1)=SITE_VECTOR(3,np1+1)
                  SITE_BRANCH(np1)= SITE_BRANCH(np1+1)
                ENDDO
                SITE=SITE-1
                NUM_SITES=NUM_SITES-1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDDO


      choice=3
      if (choice.eq.1)then
C----------------------
        CALL OPENF(IOFILE1,'DISK',FILE03(IBEG:IEND)//'.coronary',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

        WRITE(IFILE,*) num_sites
        DO SITE=1,NUM_SITES
          DO SITE_COUNT=1,4
            WRITE(IOFILE1,*) SITE_NUM(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VALUE(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VECTOR(SITE_COUNT,SITE)
          ENDDO
          WRITE(IOFILE1,*) SITE_BRANCH(SITE) ! 1=RCA,2=LAD,3=CX

        ENDDO
        CALL CLOSEF(IOFILE1,ERROR,*9999)
      endif


      if (choice.eq.2)then
C----------------------
        SITE_IN_FILE=INT(NUM_SITES/3)

        CALL OPENF(IOFILE1,'DISK',FILE03(IBEG:IEND)//'1.coronary',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

        WRITE(IOFILE1,*) SITE_IN_FILE
        DO SITE=1,NUM_SITES
          DO SITE_COUNT=1,4
            WRITE(IOFILE1,*) SITE_NUM(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VALUE(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VECTOR(SITE_COUNT,SITE)
          ENDDO
          WRITE(IOFILE1,*) SITE_BRANCH(SITE) ! 1=RCA,2=LAD,3=CX

        ENDDO
        CALL CLOSEF(IOFILE1,ERROR,*9999)


        CALL OPENF(IOFILE1,'DISK',FILE03(IBEG:IEND)//'2.coronary',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

        WRITE(IOFILE1,*) SITE_IN_FILE
        DO SITE=(SITE_IN_FILE+1),(2*SITE_IN_FILE)
          DO SITE_COUNT=1,4
            WRITE(IOFILE1,*) SITE_NUM(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VALUE(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VECTOR(SITE_COUNT,SITE)
          ENDDO
          WRITE(IOFILE1,*) SITE_BRANCH(SITE) ! 1=RCA,2=LAD,3=CX

        ENDDO
        CALL CLOSEF(IOFILE1,ERROR,*9999)

       CALL OPENF(IOFILE1,'DISK',FILE03(IBEG:IEND)//'3.coronary',
     '    'UNKNOWN','SEQUEN','FORMATTED',160,ERROR,*9999)

       WRITE(IOFILE1,*) NUM_SITES-(2*SITE_IN_FILE)
        DO SITE=(2*SITE_IN_FILE+1),NUM_SITES
          DO SITE_COUNT=1,4
            WRITE(IOFILE1,*) SITE_NUM(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VALUE(SITE_COUNT,SITE)
          ENDDO
          DO SITE_COUNT=1,3
            WRITE(IOFILE1,*) SITE_VECTOR(SITE_COUNT,SITE)
          ENDDO
          WRITE(IOFILE1,*) SITE_BRANCH(SITE) ! 1=RCA,2=LAD,3=CX
        ENDDO
        CALL CLOSEF(IOFILE1,ERROR,*9999)
C-----------------------
      endif

c      write(*,*) 'enter 1 write avoid form file 2 for not'
c      read (*,*) choice
      choice=1
      IF (choice.eq.1)THEN
       CALL OPENF(IOFILE1,'DISK','avoid_next.data','UNKNOWN',
     '  'SEQUEN','FORMATTED',160,ERROR,*9999)
      write(IOFILE1,*) AVOID_PT_COUNT
      DO I=1,AVOID_PT_COUNT
        write(IOFILE1,*) (AVOID_PT(j,i),j=1,3)
        write(IOFILE1,*) AVOID_WT(i)
      ENDDO
      DO ne=1,NEELEM(0,nr)
       ne1=NEELEM(ne,nr)
       write(IOFILE1,*) AVOID_HOST(0,ne1)
       write(IOFILE1,*) (AVOID_HOST(np1,ne1),np1=1,AVOID_HOST(0,ne1))
      ENDDO
      CALL CLOSEF(IFILE,ERROR,*9999)
      ENDIF

      CALL EXITS('IPMESH6_SUB')
      RETURN
 9999 CALL ERRORS('IPMESH6_SUB',ERROR)
      CALL EXITS('IPMESH6_SUB')
      RETURN 1
      END


