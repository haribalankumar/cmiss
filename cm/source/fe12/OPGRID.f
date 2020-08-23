      SUBROUTINE OPGRID(TYPE,FULL,GRID_LATTICE,HISTORY,ARRAY_INDEX,na,
     &  NQLIST,nr,NAQ,NEELEM,NENQ,NLATNE,NLATNQ,NLATPNQ,NLQ,NQGP,
     &  NQGP_PIVOT,NQNLAT,NQNP,NQS,NQSCNB,NQXI,NWQ,nx,nxc,NXQ,
     &  AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,PROPQ,RADIUS,XQ,YQ,YQS,
     &  ARC_LENGTH,FLUX,NIDS,POTENTIAL,RMS,ERROR,*)

C#### Subroutine: OPGRID
C###  Description:
C###    OPGRID outputs finite difference grid point arrays

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      !INCLUDE 'cmiss$reference:cell02.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
C      INCLUDE 'cmiss$reference:inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'maqloc00.inc'
      INCLUDE 'nqloc00.inc'
      INCLUDE 'solv00.cmn'
      INCLUDE 'tol00.cmn'

!     Parameter List
      INTEGER ARRAY_INDEX,na,NQLIST(0:NQM),nr,NAQ(NQM,NAM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NLATNE(NEQM+1),
     '  NLATNQ(NEQM*NQEM),NLATPNQ(NQM),NLQ(NQM),NQGP(0:NQGM,NQM),
     '  NQGP_PIVOT(NQGM,NQM),NQNLAT(NEQM*NQEM),NQNP(NPM),NQS(NEQM),
     '  NQSCNB(NQSCM),NQXI(0:NIM,NQSCM),NWQ(8,0:NQM,NAM),nx,nxc,
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 AQ(NMAQM,NQM),CQ(NMM,NQM),DNUDXQ(3,3,NQM),DXDXIQ(3,3,NQM),
     '  DXDXIQ2(3,3,NQM),GCHQ(3,NQM),GUQ(3,3,NQM),PROPQ(3,3,4,2,NQM),
     '  RADIUS(2),XQ(NJM,NQM),YQ(NYQM,NIQM),YQS(NIQSM,NQM)
      CHARACTER ERROR*(*),TYPE*(*)
      LOGICAL ARC_LENGTH,FLUX,FULL,GRID_LATTICE,HISTORY,NIDS,POTENTIAL,
     &  RMS 
!     Local Variables
      INTEGER atime,DX_NUMBER(3),DX_NUMBER_TOT,elem_latpt,end_latpt,
     '  I,IBEG,IEND,IFROMC,IJ,IK,latpt,maq,maqlist(20),maqtot,
     '  na1,na2,ne,ni,niq,nii,nij,nik,
     '  NITB,nj,njj,NJJT,NJMAX,nk,nm,nm1,nm2,nnq,np,nq,nq1,
     '  SUB_LAT(0:NQEM),entry,count,NQBOUND,nqh1(2,3),nqh1_tot(2),
     '  nqh2(2,3),nqh2_tot(2),nzero
      REAL*8 ANALYTIC,AV_DX(3),AV_DX_TOT,CALCULATED,COEFF(22),
     '  DX(3),DX_TOT,FIRSTACTIV,GDPHIDN,LASTACTIV,MAX_DX(3),
     '  MAX_DX_TOT,MIN_DX(3),MIN_DX_TOT,NIDSERROR,RMSERROR,
     '  SIZE,SUM,SUMERRORSQ,SUMEXACTSQ,TEMP(3),TNLOCAL(2),VLENGTH,
     '  XNLOCAL(3)
      CHARACTER CFROMI*2,CHAR2*2,maqheader(20)*10,NIQ_STRG*9,NM_STRG*9,
     &  NQ_STRG*5,STRG*4
C cpb 10/12/99 Compiler bug with this statement
C      VOLATILE nqh1_tot,nqh2_tot

      CALL ENTERS('OPGRID',*9999)

      IF(TYPE(1:12).NE.'CONNECTIVITY') THEN
        WRITE(OP_STRING,'(/'' Class '',I1,'' (nx='',I1,'
     '    //''') Region '',I1,'' Grid level '',I2,'':'')') nxc,nx,nr,na
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE !IF(TYPE(1:12).EQ.'CONNECTIVITY') THEN
        WRITE(OP_STRING,'('' Grid level'',I3,'' #points listed ='',I8)')
     '    na,NQLIST(0)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(TYPE(1:8).EQ.'GEOMETRY') THEN
        WRITE(OP_STRING,'('' Grid point positions:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NJJT=NJ_LOC(NJL_GEOM,0,nr)
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(OP_STRING,'('' grid point '',I8,'': '',3D12.4)')
     '      nq,(XQ(NJ_LOC(NJL_GEOM,njj,nr),nq),njj=1,NJJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nqq

      ELSE IF(TYPE(1:8).EQ.'MATERIAL') THEN
        IF(ARRAY_INDEX.NE.0) THEN !only list specified index
          nm1=ARRAY_INDEX
          nm2=ARRAY_INDEX
          WRITE(NM_STRG,'(I3)') nm1
        ELSE !list all indices
          nm1=ARRAY_INDEX
          nm2=ILT(1,nr,nx)
          WRITE(NM_STRG,'(I1,''..'',I3)') nm1,nm2
        ENDIF
        CALL STRING_TRIM(NM_STRG,IBEG,IEND)
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(OP_STRING,'('' CQ(nm='//NM_STRG(IBEG:IEND)
     '      //',nq='',I5,''): '',20D12.4)')
     '      nq,(CQ(nm,nq),nm=nm1,nm2)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nnq

      ELSE IF(TYPE(1:10).EQ.'ACTIVTIMES') THEN
        CALL ASSERT(CALC_ACTIV_TIMES,
     '    '>>You must save activation times during solve',ERROR,*9999)
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_TIME,atime,MAQ_ACTIV_TIME,
     '    ERROR,*9999)
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(OP_STRING,'('' Grid number: '',I8,'' Activation time'',
     '      F12.4)') nq,AQ(atime,nq)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        FIRSTACTIV=1.0d6
        LASTACTIV=1.0d-6
        DO nq=1,NQT
          IF(AQ(atime,nq).LT.FIRSTACTIV) FIRSTACTIV=AQ(atime,nq)
          IF(AQ(atime,nq).GT.LASTACTIV) LASTACTIV=AQ(atime,nq)
        ENDDO
        WRITE(OP_STRING,'('' Earliest activation time '',F12.4)')
     '    FIRSTACTIV
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Latest   activation time '',F12.4)')
     '    LASTACTIV
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(TYPE(1:9).EQ.'AUXILIARY') THEN
        i=0
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_X,
     '    ERROR,*9999)
        IF (maq.GT.0) THEN
          i=i+1
          maqlist(i)=maq
          maqheader(i)='  Normal_x'
        ENDIF
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Y,
     '    ERROR,*9999)
        IF (maq.GT.0) THEN
          i=i+1
          maqlist(i)=maq
          maqheader(i)='  Normal_y'
        ENDIF
        CALL MAQ_LOC(MAQ_INQUIRE,MAQ_COORD,maq,MAQ_NORMAL_Z,
     '    ERROR,*9999)
        IF (maq.GT.0) THEN
          i=i+1
          maqlist(i)=maq
          maqheader(i)='  Normal_z'
        ENDIF
        maqtot=i
        WRITE(OP_STRING,'(''    nq    '',20A10)')
     &    (maqheader(i),i=1,maqtot)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(OP_STRING,'(I10,20F10.2)') 
     &      nq,(AQ(maqlist(i),nq),i=1,maqtot)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

      ELSE IF(TYPE(1:12).EQ.'CONNECTIVITY') THEN
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(OP_STRING,'('' NXQ(-3..3,1,nq,na):'',7I7,'
     '      //'''   NWQ(1..2,nq,na):'',2I7)')
     '      (NXQ(i,1,nq,na),i=-3,3),NWQ(1,nq,na),NWQ(2,nq,na)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nnq

      ELSE IF(TYPE(1:4).EQ.'FLUX') THEN
        CALL NIQ_LOC(NIQ_INQUIRE,NIQ_ION,niq,NIQ_V,ERROR,*9999)
        IF(niq.LE.0) niq=1
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          IF(NWQ(1,nq,na).NE.0) THEN
            IF(ITYP2(nr,nx).EQ.3) THEN
              TEMP(1)=CQ(1,nq)
              TEMP(2)=CQ(2,nq)
              TEMP(3)=CQ(3,nq)
C              CQ(1,nq)=1.0d-3
C              CQ(2,nq)=1.0d-3
C              CQ(3,nq)=1.0d-3
              CQ(1,nq)=1.0d0
              CQ(2,nq)=1.0d0
              CQ(3,nq)=1.0d0
              CALL GGRADPHIQDN(NENQ,niq,nq,NQS,NQXI,NXQ(-NIM,0,0,1),AQ,
     &          CQ(1,nq),DNUDXQ,DXDXIQ,DXDXIQ2,GDPHIDN,YQ,ERROR,*9999)
              CQ(1,nq)=TEMP(1)
              CQ(2,nq)=TEMP(2)
              CQ(3,nq)=TEMP(3)
              CALL NORM31(NJT,nq,NXQ(-NIM,0,0,1),DXDXIQ,DXDXIQ2,XNLOCAL,
     '          ERROR,*9999)
              IF(NIT(NBT).EQ.1) THEN
                ANALYTIC=XNLOCAL(1)
              ELSE IF(NIT(NBT).EQ.2) THEN
                ANALYTIC=(XNLOCAL(1)*XQ(1,nq)*2.0d0)-
     '            (XNLOCAL(2)*XQ(2,nq)*2.0d0)
              ELSE IF(NIT(NBT).EQ.3) THEN
                ANALYTIC=(XNLOCAL(1)*XQ(1,nq)*2.0d0)+
     '            (XNLOCAL(2)*XQ(2,nq)*2.0d0)-
     '            (XNLOCAL(3)*XQ(3,nq)*4.0d0)
              ENDIF

              !MLB new, check CALC_GRID_BOUND_COEF get the correct
              !values, should end up with the same as GGRADPHIQDN
              TEMP(1)=CQ(1,nq)
              TEMP(2)=CQ(2,nq)
              TEMP(3)=CQ(3,nq)
              CQ(1,nq)=1.0d0
              CQ(2,nq)=1.0d0
              CQ(3,nq)=1.0d0
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,
     '          NXQ(-NIM,0,0,1),COEFF,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          ERROR,*9999)
              CQ(1,nq)=TEMP(1)
              CQ(2,nq)=TEMP(2)
              CQ(3,nq)=TEMP(3)
              SUM=0.0d0
              DO nzero=1,NQGP(0,nq)
                SUM=SUM+COEFF(NQGP_PIVOT(nzero,nq))*
     '            YQ(NQGP(nzero,nq),1)
              ENDDO

              WRITE(OP_STRING,
     '          '('' Grid: '',I8,'' g*dphi/dn: '',3F12.4)')
     '          nq,GDPHIDN,ANALYTIC,SUM
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              CALL GGRADPHIQDN(NENQ,niq,nq,NQS,NQXI,NXQ(-NIM,0,0,1),AQ,
     &          CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,GDPHIDN,YQ,ERROR,*9999)
              CALL CALC_GRID_BOUND_COEF(NENQ,nq,NQS,NQXI,
     '          NXQ(-NIM,0,0,1),COEFF,CQ(6,nq),DNUDXQ,DXDXIQ,DXDXIQ2,
     '          ERROR,*9999)
              SUM=0.0d0
              DO nzero=1,NQGP(0,nq)
                SUM=SUM+COEFF(NQGP_PIVOT(nzero,nq))*
     '            YQ(NQGP(nzero,nq),1)
              ENDDO
              WRITE(OP_STRING,
     '          '('' Grid: '',I8,'' g*dphi/dn: '',2F12.4)')
     '          nq,GDPHIDN,SUM
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

C            CALL NORM31(nq,NXQ(-NIM,0,0,1),DXDXIQ,XNLOCAL,ERROR,*9999)
C            ANALYTIC=(XNLOCAL(1)*XQ(1,nq)*2.0d0)-
C     '        (XNLOCAL(2)*XQ(2,nq)*2.0d0)
C            WRITE(OP_STRING,'('' Grid: '',I8,'' g*dphi/dn: '',2F12.4)')
C     '        nq,GDPHIDN,ANALYTIC
C            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          ENDIF !NWQ
        ENDDO !nnq

      ELSE IF(TYPE(1:8).EQ.'ADAPTIVE') THEN
        WRITE(OP_STRING,'('' Grid levels for adaptive residual'
     '    //' calculations:'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        NITB=2 !temporary
        IK=MAX(0,NITB-2) !zero for 1,2D, one for 3D
        IJ=MIN(NITB-1,1) !zero for 1D, one for 2,3D
        DO nq=1,NQT
          IF(NLQ(nq).GT.0.AND.NLQ(nq).LT.10) THEN !residual pt
            na1=NLQ(nq) !is level for interpolation of residual
            IF(NWQ(1,nq,na1).EQ.0) THEN !nq is in interior
              STRG='    '
            ELSE                        !nq is on boundary
              STRG='Bdry'
            ENDIF
            WRITE(OP_STRING,'('' nq= '',I6,1X,A,'
     '        //'''  Compute residual at level na= NLQ(nq)= '',I2,'
     '        //'16X,'' NXQ(..,1,nq,'',I1,''): '',9I4)')
     '        nq,STRG,NLQ(nq),na1,
     '        (((NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,na1),na1),
     '        na1),nii=-1,1),nij=-IJ,IJ),nik=-IK,IK)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(NLQ(nq).GT.10) THEN             !interpolated pt
            CHAR2=CFROMI(NLQ(nq),'(I2)')
            na1=IFROMC(CHAR2(1:1))
            na2=IFROMC(CHAR2(2:2))
            IF(NWQ(1,nq,na1).EQ.0) THEN !nq is in interior
              STRG='    '
            ELSE                        !nq is on boundary
              STRG='Bdry'
            ENDIF
            IF(NAQ(nq,na2).LE.3) THEN !nq is betw 2 pts in Xi(NAQ(nq,na2)) dir.n
              WRITE(OP_STRING,'('' nq= '',I6,1X,A,'
     '          //'''  Interpolate from level above: NLQ(nq)= '',I2,'
     '          //''', NAQ(nq,'',I1,'')= '',I2,'
     '          //''', NXQ(..,1,nq,'',I1,''): '',2I4)')
     '          nq,STRG,NLQ(nq),na2,NAQ(nq,na2),na1,(NXQ(ni,1,nq,na1),
     '          ni=-NAQ(nq,na2),NAQ(nq,na2),2*NAQ(nq,na2))
            ELSE IF(NAQ(nq,na2).EQ.4) THEN !nq is betw 4 pts in Xi1-Xi2 plane
              WRITE(OP_STRING,'('' nq= '',I6,1X,A,'
     '          //'''  Interpolate from level above: NLQ(nq)= '',I2,'
     '          //''', NAQ(nq,'',I1,'')= '',I2,'
     '          //''', NXQ(..,1,nq,'',I1,''): '',4I4)')
     '          nq,STRG,NLQ(nq),na2,NAQ(nq,na2),na1,
     '          ((NXQ(nii,1,NXQ(nij*2,1,nq,na1),na1),
     '          nii=-1,1,2),nij=-1,1,2)
            ENDIF !NAQ
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nq

      ELSE IF(TYPE(1:10).EQ.'STATISTICS') THEN
        DX_NUMBER_TOT=0
        DX_TOT=0.0d0
        MIN_DX_TOT=RMAX
        MAX_DX_TOT=-RMAX
        nqh1_tot(1)=0
        nqh1_tot(2)=0
        nqh2_tot(1)=0
        nqh2_tot(2)=0
        DO ni=1,NIM
          DX_NUMBER(ni)=0
          DX(ni)=0.0d0
          MIN_DX(ni)=RMAX
          MAX_DX(ni)=-RMAX
          nqh1(1,ni)=0
          nqh1(2,ni)=0
          nqh2(1,ni)=0
          nqh2(2,ni)=0
        ENDDO
        NITB=NIT(NQSCNB(NQS(NEELEM(1,nr))))

        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          DO ni=1,NITB
            IF(NXQ(ni,1,nq,1).GT.0) THEN
              nq1=NXQ(ni,1,nq,1)
              SIZE=0.0d0
              DO nj=1,NJT
                SIZE=SIZE+(XQ(nj,nq)-XQ(nj,nq1))**2
              ENDDO
              SIZE=DSQRT(SIZE)
              DX_NUMBER_TOT=DX_NUMBER_TOT+1
              DX_NUMBER(ni)=DX_NUMBER(ni)+1
              DX_TOT=DX_TOT+SIZE
              DX(ni)=DX(ni)+SIZE
              IF(SIZE.GT.MAX_DX_TOT) THEN
                MAX_DX_TOT=SIZE
                nqh1_tot(1)=nq
                nqh1_tot(2)=nq1
              ENDIF
              IF(SIZE.GT.MAX_DX(ni)) THEN
                MAX_DX(ni)=SIZE
                nqh1(1,ni)=nq
                nqh1(2,ni)=nq1
              ENDIF
              IF(SIZE.LT.MIN_DX_TOT) THEN
                MIN_DX_TOT=SIZE
                nqh2_tot(1)=nq
                nqh2_tot(2)=nq1
              ENDIF
              IF(SIZE.LT.MIN_DX(ni)) THEN
                MIN_DX(ni)=SIZE
                nqh2(1,ni)=nq
                nqh2(2,ni)=nq1
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        AV_DX_TOT=DX_TOT/DBLE(DX_NUMBER_TOT)
        DO ni=1,NITB
          AV_DX(ni)=DX(ni)/DBLE(DX_NUMBER(ni))
        ENDDO

C *** DPN 02 December 1999 - Without these stupid lines, all optimised
C ***   versions seg fault on the fem list grid statistics command,
C ***   at least for large, 64bit problems! I'm sure there is a
C ***   better way to do this, but I don't know it!! Maybe Munce will be
C ***   able to work it out ?!?
        nqh1_tot(1)=nqh1_tot(1)
        nqh1_tot(2)=nqh1_tot(2)
        nqh2_tot(1)=nqh2_tot(1)
        nqh2_tot(2)=nqh2_tot(2)
C *** DPN 02 December 1999 - end of stupid lines!

        WRITE(OP_STRING,'('' Grid point statistics '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO ni=1,NITB
          WRITE(OP_STRING,'(/'' Xi DIRECTION '',I1)') ni
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' Minimum dx : '',F12.8)') MIN_DX(ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Occurs between : '',2I8)') nqh2(1,ni),
     '      nqh2(2,ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' In elements : '',2I8)')
     '      NENQ(1,nqh2(1,ni)),NENQ(1,nqh2(2,ni))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING,'(/'' Maximum dx : '',F12.8)') MAX_DX(ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Occurs between : '',2I8)') nqh1(1,ni),
     '      nqh1(2,ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' In elements : '',2I8)')
     '      NENQ(1,nqh1(1,ni)),NENQ(1,nqh1(2,ni))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          WRITE(OP_STRING,'(/'' Average dx : '',F12.8)') AV_DX(ni)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        WRITE(OP_STRING,'(/'' TOTAL'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Minimum dx : '',F12.8)') MIN_DX_TOT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Occurs between : '',2I8)') nqh2_tot(1),
     '    nqh2_tot(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' In elements : '',2I8)')
     '    NENQ(1,nqh2_tot(1)),NENQ(1,nqh2_tot(2))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'(/'' Maximum dx : '',F12.8)') MAX_DX_TOT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Occurs between : '',2I8)') nqh1_tot(1),
     '    nqh1_tot(2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' In elements : '',2I8)')
     '    NENQ(1,nqh1_tot(1)),NENQ(1,nqh1_tot(2))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        WRITE(OP_STRING,'(/'' Average dx : '',F12.8)') AV_DX_TOT
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      ELSE IF(TYPE(1:10).EQ.'LATTICE') THEN !List the lattice grid points
        IF(GRID_LATTICE) THEN
          WRITE(OP_STRING,'( /''Grid to lattice mapping'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( /''Grid point     '','
     '      //'''Principal & Sub lattice points'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            SUB_LAT(0)=0
            entry=NLATPNQ(nq)
            DO WHILE(NLATNQ(entry).NE.0)
              SUB_LAT(0)=SUB_LAT(0)+1
              SUB_LAT(SUB_LAT(0))=NLATNQ(entry)
              entry=NLATNQ(entry)
            ENDDO
            WRITE(OP_STRING,'(I8,7X,(100I8))')
     '        nq,NLATPNQ(nq),(SUB_LAT(entry),entry=1,SUB_LAT(0))
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)                
          ENDDO
        ELSE
          WRITE(OP_STRING,'( /''Lattice to grid mapping'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO ne=1,NEQM
            WRITE(OP_STRING,'( /''Element: '',I8'
     '        //'/''================='')') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            elem_latpt=NLATNE(ne)
            count=1
            DO WHILE(NLATNE(ne+count).EQ.0)
              count=count+1
            ENDDO
            end_latpt=(NLATNE(ne+count)-1)
            WRITE(OP_STRING,
     '        '(''Number of lattice points in this element: '',I5)')
     '        end_latpt-elem_latpt+1
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'( /''Lattice   =>  Grid point'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO latpt=elem_latpt,end_latpt
              WRITE(OP_STRING,'(I8,''  =>  '''
     '          //',I8)') latpt,NQNLAT(latpt)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF

      ELSE IF(TYPE(1:10).EQ.'NODES') THEN
        WRITE(OP_STRING,'(/'' Node number : Grid number'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        DO np=1,NPM
          IF((NQNP(np).GE.1).AND.(NQNP(np).LE.NQT)) THEN
            WRITE(OP_STRING,'('' Node number '',I6,'
     '        //''' : Grid number '',I8)') np,NQNP(np)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !np

      ELSE IF(TYPE(1:3).EQ.'YQS'.OR.TYPE(1:2).EQ.'YQ') THEN
        DO nnq=1,NQLIST(0)
          nq=NQLIST(nnq)
          WRITE(NQ_STRG,'(I5)') nq
          IF(TYPE(1:3).EQ.'YQS') THEN
            IF(ARRAY_INDEX.NE.0) THEN !only list specified index
              WRITE(OP_STRING,'('' YQS('',I5,'',nq='',I5,''): '','
     '          //'D12.4)') ARRAY_INDEX,nq,YQS(ARRAY_INDEX,nq)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE !list all indices
 !NIQST=CELL_NUM_STATE+CELL_NUM_DERIVED - DPN this should already be set
              WRITE(NIQ_STRG,'(''1..'',I3)') NIQST
              CALL STRING_TRIM(NIQ_STRG,IBEG,IEND)
              CALL WRITE_LONG(DPTYPE,(nnq-1)*NIQSM+1,1,IOFI,
     '          (nnq-1)*NIQSM+NIQST,5,5,%VAL(0),YQS,
     '          '('' YQS('//NIQ_STRG(IBEG:IEND)//',nq='//NQ_STRG(1:5)
     '          //'): '',5D12.4)',
     '          '(23X,5D12.4)',ERROR,*9999)
            ENDIF
          ELSE IF(TYPE(1:2).EQ.'YQ') THEN
            IF(ARRAY_INDEX.NE.0) THEN !only list specified index
              WRITE(OP_STRING,'('' grid point '',I8,'': '','
     '          //'D12.4)') nq,YQ(nq,ARRAY_INDEX)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE !list all indices
C ??? Could use sth like this (as above) but it needs checking!!!
C              WRITE(NIQ_STRG,'(''1..'',I3)') NIQM
C              CALL STRING_TRIM(NIQ_STRG,IBEG,IEND)
C              CALL WRITE_LONG(DPTYPE,1,NQM,IOFI,
C     '          NQM*NIQM,5,5,%VAL(0),YQ,
C     '          '('' YQ(nq='//NQ_STRG(1:5)//','//NIQ_STRG(IBEG:IEND)
C     '          //',na): '',5D12.4)',
C     '          '(25X,5D12.4)',ERROR,*9999)
C TEMP for now use:
              WRITE(OP_STRING,'('' YQ(nq='',I5,'',niq=1..): '','
     '          //'20D12.4)') nq,(YQ(nq,niq),niq=1,NIQM)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
        ENDDO !nnq

      ELSE IF(TYPE(1:7).EQ.'CURRENT') THEN
        IF(HISTORY) THEN  !Write in history file format
          WRITE(IOFI,'(I8)') NQLIST(0)
          WRITE(IOFI,'(''  1'')')
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            WRITE(IOFI,'(3F12.5)') (XQ(nj,nq), nj=1,3)
          ENDDO
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            WRITE(IOFI,'(8F10.4)') YQ(nq,1)
          ENDDO

        ELSE !Write all information in output format
          IF(ITYP5(nr,nx).EQ.2) THEN !Time integration
            IF(ITYP2(nr,nx).EQ.8.OR.ITYP2(nr,nx).EQ.9) THEN !activation
              DO nnq=1,NQLIST(0)
                nq=NQLIST(nnq)
                NJMAX=NJ_LOC(0,0,0)
                WRITE(OP_STRING,'('' nq= '',I7)') nq
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'(''    XQ(1..,nq): '',8D12.4)')
     '            (XQ(nj,nq),nj=1,NJMAX)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' YQ(nq,1..,na): '',8D12.4)')
     '            (YQ(nq,niq),niq=1,NIQM)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO
            ENDIF !ityp2
          ELSE
            WRITE(OP_STRING,'('' Grid point positions:'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            DO nnq=1,NQLIST(0)
              nq=NQLIST(nnq)
              WRITE(OP_STRING,'('' XQ(nj,'',I8,''): '',3D12.4)')
     '          nq,(XQ(nj,nq),nj=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nqq
          ENDIF !ityp5

          IF(FULL) THEN !List additional grid point arrays
            DO nnq=1,NQLIST(0)
              nq=NQLIST(nnq)
              WRITE(OP_STRING,'(  /'' **** Grid point nq = '',I8)') nq
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'(/''            NLQ(nq): '',I3)') NLQ(nq)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

              DO na=1,MAX(NMGT,1)
                WRITE(OP_STRING,'(/'' Multigrid level na = '',I2)') na
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,  '(''         NAQ(nq,na): '', I8)')
     '        NAQ(nq,na)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                IF(NAQ(nq,na).EQ.0) THEN !nq is in grid na
              WRITE(OP_STRING,'(''   NWQ( 1..6,nq,na): '',6I8)')
     '          (NWQ(i,nq,na),i=1,6)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NXQ(-3..3,1,nq,na): '',7I8)')
     '          (NXQ(i,1,nq,na),i=-3,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  IF(NXQ(-3,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X,    I8)') NXQ(-3,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXQ(-2,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X, 8X,I8)') NXQ(-2,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXQ(-1,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X,16X,I8)') NXQ(-1,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXQ( 1,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X,32X,I8)') NXQ( 1,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXQ( 2,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X,40X,I8)') NXQ( 2,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                  IF(NXQ( 3,0,nq,na).GT.1) THEN
                    WRITE(OP_STRING,'(21X,48X,I8)') NXQ( 3,2,nq,na)
                    CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF !NAQ
              ENDDO !na

              WRITE(OP_STRING,'(/''      GUQ(i,j,nq): '',9D12.4)')
     '          GUQ(1,1,nq),GUQ(1,2,nq),GUQ(1,3,nq),
     '          GUQ(2,1,nq),GUQ(2,2,nq),GUQ(2,3,nq),
     '          GUQ(3,1,nq),GUQ(3,2,nq),GUQ(3,3,nq)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''       GCHQ(k,nq): '',3D12.4)')
     '          GCHQ(1,nq),GCHQ(2,nq),GCHQ(3,nq)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''   dNUdXQ(i,j,nq): '',9D12.4)')
     '          ((dNUdXQ(ni,nj,nq),nj=1,3),ni=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''   dXdXiQ(j,i,nq): '',9D12.4)')
     '          ((dXdXiQ(nj,ni,nq),ni=1,3),nj=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''     C(I)(i,j,nq): '',9D12.4)')
     '          ((PROPQ(ni,nj,1,1,nq),nj=1,3),ni=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''     C(E)(i,j,nq): '',9D12.4)')
     '          ((PROPQ(ni,nj,1,2,nq),nj=1,3),ni=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''   C(I),k(i,j,nq): '',9D12.4)')
     '          (((PROPQ(ni,nj,nk+1,1,nq),nk=1,3),nj=1,3),ni=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'( ''   C(E),k(i,j,nq): '',9D12.4)')
     '          (((PROPQ(ni,nj,nk+1,2,nq),nk=1,3),nj=1,3),ni=1,3)
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDDO !nnq
          ENDIF !full
        ENDIF !history
      ELSE IF(TYPE(1:7).EQ.'ERRORS') THEN
        CALL ASSERT(ITYP6(nr,nx).EQ.1,'>>Invalid equation type',
     '    ERROR,*9999)
C SGM 26Oct2000 grid-based Finite element also
C MLT 29Nov02 Grid Finite Volume also
        CALL ASSERT(ITYP4(nr,nx).EQ.4.OR.ITYP4(nr,nx).EQ.6.OR.
     '              ITYP4(nr,nx).EQ.7,
     '    '>>Invalid equation type',ERROR,*9999)
        CALL ASSERT(ITYP2(nr,nx).EQ.3,'>>Invalid equation type',
     '    ERROR,*9999)

        NITB=NIT(NQSCNB(NQS(NEELEM(1,nr))))

        IF(POTENTIAL) THEN
          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' Number of grid points selected: '',I8)') NQLIST(0)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          SUMERRORSQ=0.0d0
          SUMEXACTSQ=0.0d0
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            IF(NITB.EQ.1) THEN
              ANALYTIC=0.0d0
              DO nj=1,NJT
                ANALYTIC=ANALYTIC+XQ(nj,nq)**2.0d0
              ENDDO
              ANALYTIC=DSQRT(ANALYTIC)
            ELSE IF(NITB.EQ.2) THEN
              ANALYTIC=XQ(1,nq)**2.0d0-XQ(2,nq)**2.0d0
            ELSE IF(NITB.EQ.3) THEN
              ANALYTIC=XQ(1,nq)**2.0d0+XQ(2,nq)**2.0d0-
     '          (2.0d0*(XQ(3,nq)**2.0d0))
            ELSE
              ERROR='>>NIT must be >=1 and <=3'
              GOTO 9999
            ENDIF
            SUMEXACTSQ=SUMEXACTSQ+ANALYTIC**2.0d0
            SUMERRORSQ=SUMERRORSQ+(ANALYTIC-YQ(nq,1))**2.0d0

            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' Grid,calculated,analytic '',I8,2F12.8)')
     '          nq,YQ(nq,1),ANALYTIC
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nq

          IF(RMS) THEN
            IF(NQLIST(0).GT.0) THEN
              RMSERROR=DSQRT(SUMERRORSQ/DBLE(NQLIST(0)))
            ELSE
              ERROR='>>No grid points selected'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The RMS error in potential is: '',F12.8)') RMSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !rms
          IF(NIDS) THEN
            IF(SUMEXACTSQ.GT.ZERO_TOL) THEN
              NIDSERROR=DSQRT(SUMERRORSQ/SUMEXACTSQ)
            ELSE
              ERROR='>>Divide by zero in OPGRID'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The NIDS error in potential is: '',F12.8)') NIDSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF !nids
        ENDIF

        IF(FLUX) THEN
          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          NQBOUND=0
          SUMERRORSQ=0.0d0
          SUMEXACTSQ=0.0d0
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            IF(NWQ(1,nq,na).GE.1) THEN !boundary point
              NQBOUND=NQBOUND+1
              CALL NORM31(NJT,nq,NXQ(-NIM,0,0,1),DXDXIQ,DXDXIQ2,XNLOCAL,
     '          ERROR,*9999)
              IF(NITB.EQ.1) THEN
                IF(NWQ(1,nq,na).GT.nq) THEN
                  ANALYTIC=-1.0d0
                ELSE
                  ANALYTIC=1.0d0
                ENDIF
              ELSE IF(NITB.EQ.2) THEN
                ANALYTIC=(XNLOCAL(1)*XQ(1,nq)*2.0d0)-
     '            (XNLOCAL(2)*XQ(2,nq)*2.0d0)
              ELSE IF(NITB.EQ.3) THEN
                ANALYTIC=(XNLOCAL(1)*XQ(1,nq)*2.0d0)+
     '            (XNLOCAL(2)*XQ(2,nq)*2.0d0)-
     '            (XNLOCAL(3)*XQ(3,nq)*4.0d0)
              ELSE
                ERROR='>>NIT must be >=1 and <=3'
                GOTO 9999
              ENDIF

              IF(NITB.GT.1) THEN
                CALL GGRADPHIQDN(NENQ,1,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),
     '            DNUDXQ,DXDXIQ,DXDXIQ2,CALCULATED,YQ,ERROR,*9999)
              ELSE
                CALCULATED=(3.0d0*YQ(nq,1))-
     '            (4.0d0*YQ(NWQ(1,nq,na),1))+
     '            (1.0d0*YQ(NWQ(2,nq,na),1))
                CALCULATED=CALCULATED/DXDXIQ(1,1,NWQ(1,nq,na))
              ENDIF

              SUMEXACTSQ=SUMEXACTSQ+ANALYTIC**2.0d0
              SUMERRORSQ=SUMERRORSQ+(ANALYTIC-CALCULATED)**2.0d0

              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' Grid,calculated,analytic '',I8,2F12.8)')
     '            nq,CALCULATED,ANALYTIC
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF !nwq
          ENDDO !nq

          WRITE(OP_STRING,
     '      '('' Number of boundary grid points selected: '',I8)')
     '      NQBOUND
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(RMS) THEN
            IF(NQLIST(0).GT.0) THEN
              RMSERROR=DSQRT(SUMERRORSQ/DBLE(NQBOUND))
            ELSE
              ERROR='>>No grid points selected'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The RMS error in boundary flux is: '',F12.8)')
     '        RMSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NIDS) THEN
            IF(SUMEXACTSQ.GT.ZERO_TOL) THEN
              NIDSERROR=DSQRT(SUMERRORSQ/SUMEXACTSQ)
            ELSE
              ERROR='>>Divide by zero in OPGRID'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The NIDS error in boundary flux is: '',F12.8)')
     '        NIDSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF

        IF(ARC_LENGTH) THEN
          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,
     '      '('' WARNING: valid in circular domains only'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)

          NQBOUND=0
          SUMERRORSQ=0.0d0
          SUMEXACTSQ=0.0d0
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            IF(NWQ(1,nq,na).GE.1) THEN !boundary point
              NQBOUND=NQBOUND+1

              !Create a tangent vector
              IF(NXQ(-1,1,nq,na).EQ.0) THEN
                TNLOCAL(1)=DXDXIQ(1,2,nq)
                TNLOCAL(2)=DXDXIQ(2,2,nq)
              ELSEIF(NXQ(1,1,nq,na).EQ.0) THEN
                TNLOCAL(1)=DXDXIQ(1,2,nq)
                TNLOCAL(2)=DXDXIQ(2,2,nq)
              ELSEIF(NXQ(-2,1,nq,na).EQ.0) THEN
                TNLOCAL(1)=DXDXIQ(1,1,nq)
                TNLOCAL(2)=DXDXIQ(2,1,nq)
              ELSEIF(NXQ(2,1,nq,na).EQ.0) THEN
                TNLOCAL(1)=DXDXIQ(1,1,nq)
                TNLOCAL(2)=DXDXIQ(2,1,nq)
              ENDIF
              VLENGTH=DSQRT((TNLOCAL(1)**2)+(TNLOCAL(2)**2))
              IF(VLENGTH.GE.ZERO_TOL) TNLOCAL(1)=TNLOCAL(1)/VLENGTH
              IF(VLENGTH.GE.ZERO_TOL) TNLOCAL(2)=TNLOCAL(2)/VLENGTH

              ANALYTIC=(2.0d0*XQ(1,nq)*TNLOCAL(1))-
     '          (2.0d0*XQ(2,nq)*TNLOCAL(2))

              CALL NQDS(2,1,NENQ,1,1,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),DNUDXQ,
     &          CALCULATED,DXDXIQ,DXDXIQ2,XQ,YQ,.FALSE.,ERROR,*9999)

              SUMERRORSQ=SUMERRORSQ+(ANALYTIC-CALCULATED)**2.0d0
              SUMEXACTSQ=SUMEXACTSQ+ANALYTIC**2.0d0

              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' Grid,calculated,analytic '',I8,2F12.8)')
     '            nq,CALCULATED,ANALYTIC
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF

            ENDIF !nwq
          ENDDO !nq

          WRITE(OP_STRING,
     '      '('' Number of boundary grid points selected: '',I8)')
     '      NQBOUND
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          IF(RMS) THEN
            IF(NQLIST(0).GT.0) THEN
              RMSERROR=DSQRT(SUMERRORSQ/DBLE(NQBOUND))
            ELSE
              ERROR='>>No grid points selected'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The RMS error in boundary arc length '
     '        //'derivative is: '',F12.8)') RMSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NIDS) THEN
            IF(SUMEXACTSQ.GT.ZERO_TOL) THEN
              NIDSERROR=DSQRT(SUMERRORSQ/SUMEXACTSQ)
            ELSE
              ERROR='>>Divide by zero in OPGRID'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The NIDS error in boundary arc length '
     '        //'derivative is: '',F12.8)') NIDSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF

          WRITE(OP_STRING,'('' '')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

          SUMERRORSQ=0.0d0
          SUMEXACTSQ=0.0d0
          DO nnq=1,NQLIST(0)
            nq=NQLIST(nnq)
            IF(NWQ(1,nq,na).GE.1) THEN !boundary point

              IF(nq.LT.NWQ(1,nq,na)) THEN
                IF(RADIUS(1).GT.ZERO_TOL) THEN
                  ANALYTIC=(8.0d0*XQ(1,nq)*XQ(2,nq))/
     '              (RADIUS(1)**2.0d0)
                ELSE
                  ERROR='>>Divide by zero in opgrid'
                  GOTO 9999
                ENDIF
              ELSE
                IF(RADIUS(2).GT.ZERO_TOL) THEN
                  ANALYTIC=-(8.0d0*XQ(1,nq)*XQ(2,nq))/
     '              (RADIUS(2)**2.0d0)
                ELSE
                  ERROR='>>Divide by zero in opgrid'
                  GOTO 9999
                ENDIF
              ENDIF

              CALL NQDS(2,1,NENQ,1,1,nq,NQS,NQXI,NXQ,AQ,CQ(6,nq),DNUDXQ,
     &          CALCULATED,DXDXIQ,DXDXIQ2,XQ,YQ,.TRUE.,ERROR,*9999)

              SUMERRORSQ=SUMERRORSQ+(ANALYTIC-CALCULATED)**2.0d0
              SUMEXACTSQ=SUMEXACTSQ+ANALYTIC**2.0d0

              IF(DOP) THEN
                WRITE(OP_STRING,
     '            '('' Grid,calculated,analytic '',I8,2F12.8)')
     '            nq,CALCULATED,ANALYTIC
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF

            ENDIF !nwq
          ENDDO !nq

          IF(RMS) THEN
            IF(NQLIST(0).GT.0) THEN
              RMSERROR=DSQRT(SUMERRORSQ/DBLE(NQBOUND))
            ELSE
              ERROR='>>No grid points selected'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The RMS error in the cross derivative'
     '        //' is: '',F12.8)') RMSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(NIDS) THEN
            IF(SUMEXACTSQ.GT.ZERO_TOL) THEN
              NIDSERROR=DSQRT(SUMERRORSQ/SUMEXACTSQ)
            ELSE
              ERROR='>>Divide by zero in OPGRID'
              GOTO 9999
            ENDIF
            WRITE(OP_STRING,
     '        '('' The NIDS error in the cross derivative'
     '        //' is: '',F12.8)') NIDSERROR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
      ENDIF !type

      CALL EXITS('OPGRID')
      RETURN
 9999 CALL ERRORS('OPGRID',ERROR)
      CALL EXITS('OPGRID')
      RETURN 1
      END
