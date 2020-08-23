      SUBROUTINE IPFLOW(Gdirn,ITER_USER,MECHANICS_FILETYPE,NBREATHS,
     &  NEELEM,NORD,NPLIST,NPNE,NVJE,nr,COV,CW,dt,dt_init,dt_max,
     &  ERR_USER,FlowIN,FlowOUT,I_TO_E_RATIO,MeanCompliance,Pmax,Pmin,
     &  Ppl_step,
     &  PressureIN,
     &  T_interval,refvol,RMaxMean,RMinMean,volume_target,undef,XP,
     &  DIAG_OP,NORMALISE,PATHLENGTHS,PRINT,READ_VOLUMES,UNIFORM,
     &  EXTEND_FRC,filename,P_TYPE,ERROR,*)
      
C#### Subroutine: IPFLOW
C###  Description:
C###    IPFLOW inputs parameters for evaluating a flow field

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grou00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'
      INCLUDE 'mach00.inc'

!     Parameter List
      INTEGER BELOW,C_ORDER,GDIRN,ITER_USER,MECHANICS_FILETYPE,
     &  NBREATHS,NEELEM(0:NE_R_M,0:NRM),NJJ_VRATIO,
     &  NORD(5,NE_R_M),NPLIST(0:NPM),
     &  NPNE(NNM,NBFM,NEM),nr,NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 constrict_by,COV,CW,dt,dt_init,dt_max,ERR_USER,
     &  FlowIN,FlowOUT,I_TO_E_RATIO,Pmax,Pmin,
     &  MeanCOmpliance,PPL_STEP,PressureIN,REFVOL,RMAXMEAN,RMINMEAN,
     &  T_INTERVAL,UNDEF,
     &  VOLUME_TARGET,XP(NKM,NVM,NJM,NPM)
      LOGICAL DIAG_OP,NORMALISE,PATHLENGTHS,PRINT,
     &  READ_VOLUMES,UNIFORM
      CHARACTER EXTEND_FRC*1,filename*255,ERROR*(*),P_TYPE*50
!     Local Variables
      INTEGER IBEG,ICHAR,IEND,IFROMC,INFO,N3CO,ne,njj,
     &  nn,noelem,NOQUES,np,nv,order
      REAL*8 RFROMC
      LOGICAL CBBREV,FILEIP

      CALL ENTERS('IPFLOW',*9999)

C.....SET DEFAULT VALUES
      UNIFORM=.FALSE.
      READ_VOLUMES=.FALSE.
      Ppl_step=-10.d0 * 98.0665d0 !-10 cmH2O converted to Pa
      T_interval = 2.d0      !s
      dt=0.2d0
      MeanCompliance=0.d0
      Gdirn=3
      FlowIN=0.d0
      PressureIN=0.d0
      FlowOUT=0.d0
      COV=0.d0
      BELOW=0
      C_ORDER=0
      constrict_by=1.d0
      RMaxMean=1.d0 !default ratio max to mean volume
      RMinMean=1.d0 !default ratio min to mean volume
C.....Convert default chest wall compliance of 0.2 L/cmH2O
C.....to mm^3/Pa (L/cmH2O * 1.d6mm^3/L * cmH2O/(98.0665 Pa))      
c      CW=0.2d0*1.d6/98.0665d0
      NBREATHS=5
      iter_user = 20
      err_user=1.d-6
      P_TYPE='SINUSOID'
      I_TO_E_RATIO = 0.5d0
      Pmax=1.d0
      Pmin=1.d0

      FORMAT='('' The tissue volumes are [1]:'''//
     '  '/''   (1) Calculated from a linear gradient'''//
     '  '/''   (2) Read in'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=1
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        IF(IDATA(1).EQ.1)THEN
          READ_VOLUMES=.FALSE.
        ELSE IF(IDATA(1).EQ.2)THEN
          READ_VOLUMES=.TRUE.
        ENDIF
      ENDIF
 
      IF(READ_VOLUMES)THEN
        IF(CBBREV(CO,'UNDEFORMED',3,noco+1,NTCO,N3CO)) THEN
          undef=RFROMC(CO(N3CO+1))
        ELSE
          undef=1.d0
        ENDIF
        MECHANICS_FILETYPE=1 !default
        IF(CBBREV(CO,'DATAPOINT',3,noco+1,NTCO,N3CO))THEN
          nj_Vratio=IFROMC(CO(N3CO+1))
          MECHANICS_FILETYPE=1
        ELSEIF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO))THEN
          NJJ_VRATIO=IFROMC(CO(N3CO+1))
          nj_Vratio=NJ_LOC(NJL_FIEL,NJJ_VRATIO,nr) !deformed to undeformed acinar volume ratio
          MECHANICS_FILETYPE=2
        ENDIF
      ELSEIF(CBBREV(CO,'READ_VOLUMES',4,noco+1,NTCO,N3CO)) THEN
        READ_VOLUMES=.TRUE.
        IF(CBBREV(CO,'UNDEFORMED',3,noco+1,NTCO,N3CO)) THEN
          undef=RFROMC(CO(N3CO+1))
        ELSE
          undef=1.d0
        ENDIF
        MECHANICS_FILETYPE=1 !default
        IF(CBBREV(CO,'DATAPOINT',3,noco+1,NTCO,N3CO))THEN
          nj_Vratio=IFROMC(CO(N3CO+1))
          MECHANICS_FILETYPE=1
        ELSEIF(CBBREV(CO,'FIELD',3,noco+1,NTCO,N3CO))THEN
          NJJ_VRATIO=IFROMC(CO(N3CO+1))
          nj_Vratio=NJ_LOC(NJL_FIEL,NJJ_VRATIO,nr) !deformed to undeformed acinar volume ratio
          MECHANICS_FILETYPE=2
        ENDIF

      ELSE

        FORMAT='($,'' Enter the COV for compliance [0.1]: '',D12.4)'
        RDEFLT(1)=0.1d0
        IF(IOTYPE.EQ.3) RDATA(1)=COV
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) COV=RDATA(1)

        FORMAT='($,'' Enter the ratio of maximum to mean '//
     &   'expansion [1.0]: '',D12.4)'
        RDEFLT(1)=1.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=RMaxMean
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) RMaxMean=RDATA(1)

        FORMAT='($,'' Enter the ratio of minimum to mean '//
     &   'expansion [1.0]: '',D12.4)'
        RDEFLT(1)=1.0d0
        IF(IOTYPE.EQ.3) RDATA(1)=RMinMean
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &    FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &    IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) RMinMean=RDATA(1)

      ENDIF

      FORMAT='($,'' Enter the ratio of reference to initial '//
     &  'volume [0.5]: '',D12.4)'
      RDEFLT(1)=0.5d0
      IF(IOTYPE.EQ.3) RDATA(1)=refvol
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) refvol=RDATA(1)

      FORMAT='($,'' Enter the chest wall compliance '//
     &  'in L/cmH2O [0.2]: '',D12.4)'
      RDEFLT(1)=0.2d0
      IF(IOTYPE.EQ.3) RDATA(1)=CW/1.d6*98.0665d0
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CW=RDATA(1)*1.d6/98.0665d0

      FORMAT='($,'' Enter the initial time-step size [0.01]: '',D12.4)'
      RDEFLT(1)=0.01d0
      IF(IOTYPE.EQ.3) RDATA(1)=DT
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DT=RDATA(1)
      DT_INIT = DT

      FORMAT='($,'' Enter the maximum time-step size [0.01]: '',D12.4)'
      RDEFLT(1)=0.01d0
      IF(IOTYPE.EQ.3) RDATA(1)=DT
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) DT_MAX=RDATA(1)

      FORMAT='($,'' Enter the tidal volume in L [1]: '',D12.4)'
      RDEFLT(1)=1.d0
      IF(IOTYPE.EQ.3) RDATA(1)=volume_target/1.d6
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) volume_target=RDATA(1)*1.d6

      FORMAT='($,'' Enter the change in driving pressure '//
     &  'in cmH2O [-5]: '',D12.4)'
      RDEFLT(1)=-5.d0
      IF(IOTYPE.EQ.3) RDATA(1)=Ppl_step/98.0665d0
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) Ppl_step=RDATA(1)* 98.0665d0

      FORMAT='($,'' Enter the weighting for non-dependent '//
     &  'driving pressure in cmH2O [1.0]: '',D12.4)'
      RDEFLT(1)=1.d0
      IF(IOTYPE.EQ.3) RDATA(1)=Pmin/98.0665d0 !convert to cmH2O
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) Pmin=RDATA(1)*Ppl_step

      FORMAT='($,'' Enter the weighting for dependent '//
     &  'driving pressure in cmH2O [1.0]: '',D12.4)'
      RDEFLT(1)=1.d0
      IF(IOTYPE.EQ.3) RDATA(1)=Pmax/98.0665d0 !convert to cmH2O
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) Pmax=RDATA(1)*Ppl_step

      FORMAT='($,'' The model orientation is [1]:'''//
     '  '/''   (1) Upright'''//
     '  '/''   (2) Supine'''//
     '  '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IF(GDirn.EQ.2) IDATA(1)=GDirn
         IF(GDirn.EQ.3) IDATA(1)=1
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        IF(IDATA(1).EQ.2) GDirn=2
        IF(IDATA(1).EQ.1) GDirn=3
      ENDIF         

      FORMAT='($,'' Enter the breath duration in s [4.0] '',D12.4)'
      RDEFLT(1)=4.d0
      IF(IOTYPE.EQ.3) RDATA(1)=T_interval
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) T_interval=RDATA(1)

      FORMAT='($,'' Enter the insp:expn ratio      [0.5] '',D12.4)'
      RDEFLT(1)=0.5d0
      IF(IOTYPE.EQ.3) RDATA(1)=I_TO_E_RATIO
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) I_TO_E_RATIO=RDATA(1)

      FORMAT='($,'' Return to initial FRC by extending expn? '',D12.4)'
      ADEFLT(1)='N'
      IF(IOTYPE.EQ.3) ADATA(1)=EXTEND_FRC
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,
     &  FORMAT,1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) EXTEND_FRC=ADATA(1)

ccc nj_pressure
      FORMAT='($,'' Enter the field number for pressure:      '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_pressure !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_pressure=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_flow
      FORMAT='($,'' Enter the field number for flow:          '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_flow !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_flow=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_sV
      FORMAT='($,'' Enter the field number for specific vent: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_sV !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_sV=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_Ppos
      FORMAT='($,'' Enter the field number for max -ve press: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_Ppos !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_Ppos=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_Pneg
      FORMAT='($,'' Enter the field number for max +ve press: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_Pneg !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_Pneg=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_Vratio
      FORMAT='($,'' Enter the field number for ratio def:und: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_Vratio !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_Vratio=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_Pe
      FORMAT='($,'' Enter the field number for Peclet number: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_Pe !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_Pe=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
ccc nj_mu
      FORMAT='($,'' Enter the field number for shear stress : '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_mu !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_mu=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF
      write(*,*) 'SHEAR', nj_mu,njj
ccc   nj_vent
      FORMAT='($,'' Enter the field number for time-ave vent: '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nj_vent !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        njj=IDATA(1)
        nj_vent=NJ_LOC(NJL_FIEL,njj,nr)
      ENDIF

ccc nm_R
      FORMAT='($,'' Enter the material index for R_airway  :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_R !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_R=IDATA(1)
      ENDIF
ccc nm_Rt
      FORMAT='($,'' Enter the material index for R_total   :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_Rt !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_Rt=IDATA(1)
      ENDIF
ccc nm_Ppl
      FORMAT='($,'' Enter the material index for P_pleural :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_Ppl !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_Ppl=IDATA(1)
      ENDIF
ccc nm_C
      FORMAT='($,'' Enter the material index for compliance:  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_C !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_C=IDATA(1)
      ENDIF
ccc nm_vol_bel
      FORMAT='($,'' Enter the material index for vol below :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_vol_bel !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_vol_bel=IDATA(1)
      ENDIF
ccc nm_volumes
      FORMAT='($,'' Enter the material index for airway vol:  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_volumes !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_volumes=IDATA(1)
      ENDIF
ccc nm_dpl
      FORMAT='($,'' Enter the material index for D_Ppl     :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_dpl !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_dpl=IDATA(1)
      ENDIF
ccc nm_vinit
      FORMAT='($,'' Enter the material index for vol_init  :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_vinit !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_vinit=IDATA(1)
      ENDIF
ccc nm_lambda
      FORMAT='($,'' Enter the material index for lambda    :  '',I2)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3)THEN
         IDATA(1)=nm_lambda !wrong, needs to be NJJ
      ENDIF
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMM,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3)THEN
        nm_lambda=IDATA(1)
      ENDIF
      write(*,*) 'material index', nm_dpl
C get the following from command line only

      IF(CBBREV(CO,'FILENAME',4,noco+1,NTCO,N3CO)) THEN
        filename=CO(N3CO+1)
        CALL STRING_TRIM(filename,IBEG,IEND)
        PRINT=.TRUE.
      ELSE
        PRINT=.FALSE.
      ENDIF

      IF(CBBREV(CO,'NORMALISE',3,noco+1,NTCO,N3CO)) THEN
        NORMALISE=.TRUE.
      ELSE
        NORMALISE=.FALSE.
      ENDIF
      IF(CBBREV(CO,'NBREATHS',3,noco+1,NTCO,N3CO)) THEN
        NBREATHS=IFROMC(CO(N3CO+1))
      ENDIF
      write(*,*) 'NBREATHS', NBREATHS      


      IF(CBBREV(CO,'PATHLENGTHS',3,noco+1,NTCO,N3CO)) THEN
        PATHLENGTHS=.TRUE.
      ELSE
        PATHLENGTHS=.FALSE.
      ENDIF

      IF(CBBREV(CO,'PEEP',3,noco+1,NTCO,N3CO)) THEN
        PEEP=RFROMC(CO(N3CO+1)) !*98.0665d0
      ELSE
	PEEP=0.0d0
      ENDIF

C parameters that are defined in the file, but can be changed on
C the command line

      DIAG_OP=.FALSE.
      IF(CBBREV(CO,'DIAGNOSTICS',3,noco+1,NTCO,N3CO)) 
     &  DIAG_OP=.TRUE.
      IF(CBBREV(CO,'COV',3,noco+1,NTCO,N3CO)) 
     &  COV=RFROMC(CO(N3CO+1)) !cov for tissue compliance
          write(*,*) 'COV', COV
      IF(CBBREV(CO,'REFERENCE',3,noco+1,NTCO,N3CO))
     &  refvol=RFROMC(CO(N3CO+1)) !ratio of reference volume to initial volume
      IF(CBBREV(CO,'CHESTWALL',3,noco+1,NTCO,N3CO))
     &  CW=RFROMC(CO(N3CO+1))*1.d6/98.0665d0 !chest wall compliance
      IF(CBBREV(CO,'DT_INIT',4,noco+1,NTCO,N3CO))
     &  dt=RFROMC(CO(N3CO+1))
      dt_init=dt
      IF(CBBREV(CO,'DT_MAX',4,noco+1,NTCO,N3CO))
     &  dt_max=RFROMC(CO(N3CO+1))
      IF(CBBREV(CO,'ERROR',3,noco+1,NTCO,N3CO))
     &  err_user=RFROMC(CO(N3CO+1))
      IF(CBBREV(CO,'ITERATIONS',3,noco+1,NTCO,N3CO))
     &  iter_user=IFROMC(CO(N3CO+1))
      IF(CBBREV(CO,'RMAX',3,noco+1,NTCO,N3CO))
     &  RMaxMean=RFROMC(CO(N3CO+1)) !ratio of maximum to mean volume
        write(*,*) 'RMAX', RMaxMean
      IF(CBBREV(CO,'RMIN',3,noco+1,NTCO,N3CO))
     &  RMinMean=RFROMC(CO(N3CO+1)) !ratio of maximum to minimum volume
      IF(CBBREV(CO,'TARGET',3,noco+1,NTCO,N3CO))
     &  volume_target=RFROMC(CO(N3CO+1))*1.d6
        write(*,*) 'target',volume_target
      IF(CBBREV(CO,'PPL_STEP',3,noco+1,NTCO,N3CO))
     &  Ppl_step=RFROMC(CO(N3CO+1))*98.0665d0
        write(*,*) 'ppl_step', Ppl_step
      IF(CBBREV(CO,'DIRECTION',3,noco+1,NTCO,N3CO))
     &  GDirn=IFROMC(CO(N3CO+1)) !direction of gravity: 1,2,3
       write(*,*) 'direction',GDirn
      IF(CBBREV(CO,'INTERVAL',3,noco+1,NTCO,N3CO))
     &   T_interval=RFROMC(CO(N3CO+1))
      IF(CBBREV(CO,'SINUSOID',3,noco+1,NTCO,N3CO))
     &   P_TYPE='SINUSOID'
      IF(CBBREV(CO,'LINEAR',3,noco+1,NTCO,N3CO))
     &   P_TYPE='LINEAR'
      IF(CBBREV(CO,'CURVE',3,noco+1,NTCO,N3CO))
     &   P_TYPE='CURVE'
     

c      IF(CBBREV(CO,'INFLOW',3,noco+1,NTCO,N3CO)) THEN
c        FlowIN=RFROMC(CO(N3CO+1))
c        FIX_Q_IN=.TRUE.
c      ENDIF
c      IF(CBBREV(CO,'INPRESSURE',3,noco+1,NTCO,N3CO)) THEN
c        PressureIN=RFROMC(CO(N3CO+1)) !Pa
c        FIX_P_IN=.TRUE.
c      ENDIF
c      IF(CBBREV(CO,'ZERO_FLOW',3,noco+1,NTCO,N3CO)) THEN
c        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
c        CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
c     '    ERROR,*9999)
c      ELSE
c        NELIST(0)=0
c        NPLIST(0) = 0
c      ENDIF

      IF(CBBREV(CO,'DEPOSITION',3,noco+1,NTCO,N3CO)) THEN
        UNIFORM=.TRUE. ! Super Hack, export particles
        T_interval=-1.0d0
      ENDIF

      IF(CBBREV(CO,'UNIFORM',3,noco+1,NTCO,N3CO)) THEN
        UNIFORM=.TRUE. 
      ENDIF

      IF(CBBREV(CO,'CONSTRICT',3,noco+1,NTCO,N3CO)) THEN
        constrict_by=RFROMC(CO(N3CO+1))
        IF(CBBREV(CO,'BELOW',3,noco+1,NTCO,N3CO)) THEN
          BELOW=IFROMC(CO(N3CO+1))
        ENDIF
        IF(CBBREV(CO,'ORDER',3,noco+1,NTCO,N3CO)) THEN
          C_ORDER=IFROMC(CO(N3CO+1))
        ENDIF



        DO noelem=1,NEELEM(0,1)
          ne=NEELEM(noelem,1)
          order=NORD(2,ne) !Horsfield order
          IF(BELOW.GT.0)THEN
            IF(order.LE.BELOW)THEN !constrict
              DO nn=1,2
                np=NPNE(nn,1,ne) !which node at start and end of element
                nv=NVJE(nn,1,nj_radius,ne) !which version of node in the element
                XP(1,nv,nj_radius,np)=XP(1,nv,nj_radius,np)*constrict_by
              ENDDO !nn
            ENDIF !order
          ELSE IF(C_ORDER.GT.0)THEN
            IF(order.EQ.C_ORDER)THEN !constrict
              DO nn=1,2
                np=NPNE(nn,1,ne) !which node at start and end of element
                nv=NVJE(nn,1,nj_radius,ne) !which version of node in the element
                XP(1,nv,nj_radius,np)=XP(1,nv,nj_radius,np)*constrict_by
              ENDDO !nn
            ENDIF !order
          ENDIF
        ENDDO !noelem
      ENDIF !constrict

      CALL EXITS('IPFLOW')
      RETURN
 9999 CALL ERRORS('IPFLOW',ERROR)
      CALL EXITS('IPFLOW')
      RETURN 1
      END


