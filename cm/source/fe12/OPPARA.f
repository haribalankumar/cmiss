      SUBROUTINE OPPARA(TYPE,PARAM,NEELEM,NENP,NFFACE,NLLINE,NPNODE,
     '  nx,NYNR,NX_DEFINED,ERROR,*)

C#### Subroutine: OPPARA
C###  Description:
C###    OPPARA outputs parameters.

      IMPLICIT NONE
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'ktyp40.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp60.cmn'
      INCLUDE 'ktyp70.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
      INCLUDE 'mxbc00.cmn'
      INCLUDE 'ntsg00.cmn'
      INCLUDE 'ofst00.cmn'
      INCLUDE 'opti00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NFFACE(0:NF_R_M,NRM),NLLINE(0:NL_R_M,0:NRM),
     '  NPNODE(0:NP_R_M,0:NRM),nx,NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      CHARACTER TYPE*(*),PARAM*(*),ERROR*(*)
      LOGICAL NX_DEFINED
!     Local Variables
      INTEGER i,maxNENP,nb,nc,nh,nhx,njj1,njj2,nonode,np,nr,nrc,nx1,
     '  nxx,NY_MAX

      CALL ENTERS('OPPARA',*9999)

      WRITE(OP_STRING,'(/'' Parameters for nx = '',I1,'':'')') nx
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(TYPE(1:10).EQ.'DIMENSIONS') THEN
!na
        IF(PARAM(1:2).EQ.'na'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NAM    ='',I10,'' NAT(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NAM,(NAT(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nb
        IF(PARAM(1:2).EQ.'nb'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NBM    ='',I10,'' NBT             ='','
     '      //'I10)') NBM,NBT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NBFM   ='',I10,'' NBFT            ='','
     '      //'I10)') NBFM,NBFT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nc
        IF(PARAM(1:2).EQ.'nc'.OR.PARAM(1:3).EQ.'all') THEN
          IF(NX_DEFINED) THEN
            WRITE(OP_STRING,'('' NCM    ='',I10,'' NCT(0..nr)      ='','
     '        //'4I10,/:(37X,4I10))') NCM,(NCT(nr,nx),nr=0,NRT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' NCM    ='',I10)') NCM
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' NCOM   ='',I10,'' NTCNTR          ='','
     '      //'I10)') NCOM,NTCNTR
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nd
        IF(PARAM(1:2).EQ.'nd'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NDM    ='',I10,'' NDT             ='','
     '      //'I10)') NDM,NDT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NDEM   ='',I10)') NDEM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!ne
        IF(PARAM(1:2).EQ.'ne'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NEM    ='',I10,'' NET(0..nr)      ='','
     '      //'4I10,/:(37X,4I10))') NEM,(NET(nr),nr=0,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NE_R_M ='',I10,'' NEELEM(0,0..nr) ='','
     '      //'4I10,/:(37X,4I10))') NE_R_M,(NEELEM(0,nr),nr=0,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NELM   ='',I10)') NELM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          maxNENP=0
          DO nr=1,NRT
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
              IF(NENP(np,0,0).GT.maxNENP) maxNENP=NENP(np,0,0)
            ENDDO !nonode (np)
          ENDDO !nr
          WRITE(OP_STRING,'('' NEPM   ='',I10,'' Max NENP        ='','
     '      //'4I10,/:(37X,4I10))') NEPM,maxNENP
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nf
        IF(PARAM(1:2).EQ.'nf'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NFM    ='',I10,'' NFT             ='','
     '      //'I10)') NFM ,NFT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NF_R_M ='',I10,'' NFFACE(0,nr)    ='','
     '      //'4I10,/:(37X,4I10))') NF_R_M,(NFFACE(0,nr),nr=1,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!ng
        IF(PARAM(1:2).EQ.'ng'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NGM    ='',I10,'' NGT(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NGM,(NGT(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NGRSEGM='',I10)') NGRSEGM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nh
        IF(PARAM(1:2).EQ.'nh'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NHM    ='',I10)') NHM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!ni
        IF(PARAM(1:2).EQ.'ni'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NIM    ='',I10,'' NIT(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NIM,(NIT(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIQM   ='',I10)') NIQM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIFEXTM='',I10)') NIFEXTM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIYM   ='',I10)') NIYM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIYFIXM='',I10)') NIYFIXM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIYGM  ='',I10)') NIYGM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIYGFM ='',I10)') NIYGFM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nj
        IF(PARAM(1:2).EQ.'nj'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NJM    ='',I10,'' NJT             ='','
     '      //'I10)') NJM,NJT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nk
        IF(PARAM(1:2).EQ.'nk'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NKM    ='',I10,'' NKT(0,nb)       ='','
     '      //'4I10,/:(37X,4I10))') NKM,(NKT(0,nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nl
        IF(PARAM(1:2).EQ.'nl'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NLM    ='',I10,'' NLT             ='','
     '      //'I10)') NLM,NLT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NLCM   ='',I10)') NLCM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NL_R_M ='',I10,'' NLLINE(0,0..nr) ='','
     '      //'4I10,/:(37X,4I10))') NL_R_M,(NLLINE(0,nr),nr=0,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nm
        IF(PARAM(1:2).EQ.'nm'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NMM    ='',I10)') NMM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NMGT   ='',I10)') NMGT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nn
        IF(PARAM(1:2).EQ.'nn'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NNM    ='',I10,'' NNT(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NNM,(NNT(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!no
        IF(PARAM(1:2).EQ.'no'.OR.PARAM(1:3).EQ.'all') THEN
          IF(NX_DEFINED) THEN
            WRITE(OP_STRING,'('' NOM    ='',I10,'
     '        //''' NOT(1,1,0..nr,nx) ='','
     '        //'4I10,/:(37X,4I10))') NOM,(NOT(1,1,nr,nx),nr=0,NRT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NOM    ='',I10,'
     '        //''' NOT(2,1,0..nr,nx) ='','
     '        //'4I10,/:(37X,4I10))') NOM,(NOT(2,1,nr,nx),nr=0,NRT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' NOM    ='',I10)') NOM
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' NOPM   ='',I10,'' NTOPTI          ='','
     '      //'I10)') NOPM,NTOPTI
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NOYM   ='',I10)') NOYM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!np
        IF(PARAM(1:2).EQ.'np'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NPM    ='',I10,'' NPT(0..nr)      ='','
     '      //'4I10,/:(37X,4I10))') NPM,(NPT(nr),nr=0,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NP_R_M ='',I10,'' NPNODE(0,0..nr) ='','
     '      //'4I10,/:(37X,4I10))') NP_R_M,(NPNODE(0,nr),nr=0,NRT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nq
        IF(PARAM(1:2).EQ.'nq'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NQM    ='',I10,'' NQT             ='','
     '      //'I10)') NQM,NQT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NQEM   ='',I10)') NQEM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NQGM   ='',I10)') NQGM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NQSCM  ='',I10)') NQSCM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NMAQM  ='',I10)') NMAQM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nr
        IF(PARAM(1:2).EQ.'nr'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NRM    ='',I10,'' NRT             ='','
     '      //'I10)') NRM,NRT
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NRCM   ='',I10)') NRCM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NREM   ='',I10,'' NT_RES          ='','
     '      //'I10)') NREM,NT_RES
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!ns
        IF(PARAM(1:2).EQ.'ns'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NSM    ='',I10,'' NST(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NSM,(NST(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NSFM   ='',I10)') NSFM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nt
        IF(PARAM(1:2).EQ.'nt'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NTM    ='',I10)') NTM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

CC AJPs 13-11-97
!nts
        IF(PARAM(1:3).EQ.'nts'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NTSM   ='',I10,'' NTST            ='','
     '      //'I10)') NTSM,NTST
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
CC AJPe

!nu
        IF(PARAM(1:2).EQ.'nu'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NUM    ='',I10,'' NUT(nb)         ='','
     '      //'4I10,/:(37X,4I10))') NUM,(NUT(nb),nb=1,NBFT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nv
        IF(PARAM(1:2).EQ.'nv'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NVM    ='',I10)') NVM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nw
        IF(PARAM(1:2).EQ.'nw'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NWM    ='',I10)') NWM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!nx
        IF(PARAM(1:2).EQ.'nx'.OR.PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NXM    ='',I10,'' NX_LIST(0)      ='','
     '      //'I10)') NXM,NX_LIST(0)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
!ny
        IF(PARAM(1:2).EQ.'ny'.OR.PARAM(1:3).EQ.'all') THEN
          IF(NX_DEFINED) THEN
            NY_MAX=0
            DO nr=0,NRT
              DO nc=1,NCT(nr,nx)
                DO nrc=0,NRCM
                  IF(NYNR(NYNR(0,nrc,nc,nr,nx),nrc,nc,nr,nx).GT.NY_MAX)
     '              NY_MAX=NYNR(NYNR(0,nrc,nc,nr,nx),nrc,nc,nr,nx)
                ENDDO !nrc
              ENDDO !nc
            ENDDO !nr
            WRITE(OP_STRING,'('' NYM    ='',I10,'' NY_MAX          ='','
     '        //'I10)') NYM,NY_MAX
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYROWM ='',I10,'' NYT(1,nc,nx)    ='','
     '        //'4I10,/:(37X,4I10))') NYROWM,(NYT(1,nc,nx),nc=1,NCM)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYROWM ='',I10,'' NYT(2,nc,nx)    ='','
     '        //'4I10,/:(37X,4I10))') NYROWM,(NYT(2,nc,nx),nc=1,NCM)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NY_R_M ='',I10,'
     '        //''' NYNR(0,0,1,0..nr,nx)='',3I10,/:(37X,4I10))')
     '        NY_R_M,(NYNR(0,0,1,nr,nx),nr=0,NRT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE
            WRITE(OP_STRING,'('' NYM    ='',I10)') NYM
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NYROWM ='',I10)') NYROWM
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NY_R_M ='',I10)') NY_R_M
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'('' NYOM   ='',I10)') NYOM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NYQM   ='',I10)') NYQM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
CC AJPs 13-11-97
          WRITE(OP_STRING,'('' NY_TRANSFER_M   ='',I10)') NY_TRANSFER_M
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
CC AJPe
        ENDIF
!nz
        IF(PARAM(1:2).EQ.'nz'.OR.PARAM(1:3).EQ.'all') THEN
          IF(NX_DEFINED) THEN
            WRITE(OP_STRING,'('' NZT(nc,nx) ='',4I10)')
     '        (NZT(nc,nx),nc=1,NCM)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NZZT(1,0..nr,nx)='','
     '        //'4I10,/:(37X,4I10))') (NZZT(1,nr,nx),nr=0,NRT)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
          WRITE(OP_STRING,'(/'' NZ_MINOS='',I10)') NZ_MINOSM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' NZ_FNY_M='',I10)')
     '      NZ_FNY_M
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NZ_GK_M='',I10,''  NZ_GKK_M='',I10)')
     '      NZ_GK_M,NZ_GKK_M
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NZ_GD_M='',I10,''  NZ_GQ_M='',I10)')
     '      NZ_GD_M,NZ_GQ_M
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NZ_GM_M='',I10,''  NZ_GMM_M='',I10)')
     '      NZ_GM_M,NZ_GMM_M
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NISC_GDM='',I10,''   NISR_GDM='',I10)')
     '      NISC_GDM,NISR_GDM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NISC_GKM='',I10,''   NISR_GKM='',I10)')
     '      NISC_GKM,NISR_GKM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' NISC_GKKM='',I10,''  NISR_GKKM='',I10)')
     '      NISC_GKKM,NISR_GKKM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NISC_GMM='',I10,''   NISR_GMM='',I10)')
     '      NISC_GMM,NISR_GMM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( '' NISC_GMMM='',I10,''  NISR_GMMM='',I10)')
     '      NISC_GMMM,NISR_GMMM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'( ''  NISC_GQM='',I10,''   NISR_GQM='',I10)')
     '      NISC_GQM,NISR_GQM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(PARAM(1:3).EQ.'all') THEN
          WRITE(OP_STRING,'('' NFVCM  ='',I10)') NFVCM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NVCBM  ='',I10)') NVCBM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NVCM   ='',I10)') NVCM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NORM   ='',I10)') NORM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NPDM   ='',I10)') NPDM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' NIMAGEM='',I10)') NIMAGEM
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF


!itype
      ELSE IF(TYPE(1:5).EQ.'ITYPE') THEN
        DO nr=1,NRT
          FORMAT='(/'' Region '',I1,/'
     '      //''' ITYP1 ='',I2,'' ITYP2 ='',I2,'' ITYP3 ='',I2,'
     '      //''' ITYP4 ='',I2,'' ITYP5 ='',I2,/'
     '      //''' ITYP6 ='',I2,'' ITYP7 ='',I2,'' ITYP8 ='',I2,'
     '      //''' ITYP9 ='',I2,'' ITYP10='',I2,/'
     '      //''' ITYP11='',I2,'' ITYP12='',I2,'' ITYP13='',I2,'
     '      //''' ITYP14='',I2,'' ITYP15='',I2)'
          WRITE(OP_STRING,FORMAT) nr,ITYP1(nr,nx),ITYP2(nr,nx),
     '      ITYP3(nr,nx),ITYP4(nr,nx),ITYP5(nr,nx),ITYP6(nr,nx),
     '      ITYP7(nr,nx),ITYP8(nr,nx),ITYP9(nr,nx),
     '      ITYP10(nr),ITYP11(nr),ITYP12(nr,nx),ITYP13(nr),
     '      ITYP14(nr),ITYP15(nr,nx)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nr
!jtype
      ELSE IF(TYPE(1:5).EQ.'JTYPE') THEN
        FORMAT='(/'' JTYP1 ='',I2,'' JTYP2A ='',I2,'' JTYP2B ='',I2,'
     '      //'   '' JTYP2C ='',I2,'' JTYP3 ='',I2,'
     '      //'   '' JTYP4 ='',I2,'' JTYP5 ='',I2,'' JTYP6 ='',I2,'
     '      //'/, '' JTYP7 ='',I2,'' JTYP8 ='',I2,'
     '      //'   '' JTYP10='',I2,'' JTYP12='',I2,'
     '      //'/, '' JTYP13='',I2,'' JTYP14='',I2,'' JTYP15='',I2)'
        WRITE(OP_STRING,FORMAT) JTYP1,JTYP2A,JTYP2B,JTYP2C,
     '    JTYP3,JTYP4,JTYP5,
     '    JTYP6,JTYP7,JTYP8,JTYP10,JTYP12,JTYP13,
     '    JTYP14,JTYP15
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!ktype
      ELSE IF(TYPE(1:5).EQ.'KTYPE') THEN
        FORMAT='(/'' KTYP1 ='',I2,'' KTYP2 ='',I2,'' KTYP3 ='',I2,'
     '      //'   '' KTYP4 ='',I2,'' KTYP5 ='',I2,'' KTYP6 ='',I2,'
     '      //'/, '' KTYP7 ='',I2,'' KTYP8 ='',I2,'' KTYP9 ='',I2,'
     '      //'   '' KTYP10='',I2,'' KTYP11='',I2,'' KTYP12='',I2,'
     '      //'/, '' KTYP13='',I2,'' KTYP14='',I2,'' KTYP15='',I2,'
     '      //'   '' KTYP16='',I2,'' KTYP17='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP1,KTYP2,KTYP3 ,KTYP4 ,KTYP5 ,
     '    KTYP6,KTYP7 ,KTYP8 ,KTYP9 ,KTYP10,KTYP11,KTYP12,KTYP13,
     '    KTYP14,KTYP15,KTYP16,KTYP17
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' KTYP19='',I2,'' KTYP1A='',I2,'' KTYP1B='',I2,'
     '      //'   '' KTYP1C='',I2,'' KTYP1D='',I2,'' KTYP1E='',I2,'
     '      //'/, '' KTYP20='',I2,'' KTYP21='',I2,'' KTYP22='',I2,'
     '      //'   '' KTYP23='',I2,'' KTYP24='',I2,'' KTYP25='',I2,'
     '      //'/, '' KTYP26='',I2,'' KTYP27='',I2,'' KTYP28='',I2,'
     '      //'   '' KTYP29='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP19,KTYP1A,KTYP1B,KTYP1C,KTYP1D,
     '    KTYP1E,KTYP20,KTYP21,KTYP22,KTYP23,KTYP24,KTYP25,KTYP26,
     '    KTYP27,KTYP28,KTYP29
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' KTYP30='',I2,'' KTYP31='',I2,'' KTYP32='',I2,'
     '      //'   '' KTYP33='',I2,'' KTYP34='',I2,'' KTYP35='',I2,'
     '      //'/, '' KTYP36='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP30,KTYP31,KTYP32,
     '    KTYP33,KTYP34,KTYP35,KTYP36
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' KTYP40='',I2,'' KTYP41='',I2,'' KTYP42='',I2,'
     '      //'   '' KTYP43='',I2,'' KTYP44='',I2,'' KTYP45='',I2)'
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,FORMAT) KTYP40,KTYP41,KTYP42,KTYP43,KTYP44,
     '    KTYP45
        DO nr=1,NRT
          FORMAT='(/'' Region '',I1)'
          WRITE(OP_STRING,FORMAT) nr
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='( '' KTYP50='',I2,'' KTYP51='',I2,'' KTYP52='',I2,'
     '      //'   '' KTYP53='',I2,'' KTYP54='',I2,'' KTYP55='',I2,'
     '      //'   '' KTYP55a='',I2,'
     '      //'/  '' KTYP56='',I2,'' KTYP57='',I2,'' KTYP58='',I2,'
     '      //'   '' KTYP59='',I2,'' KTYP5A='',I2,'' KTYP5B='',I2,'
     '      //'/  '' KTYP5C='',I2,'' KTYP5D='',I2,'' KTYP5E='',I2,'
     '      //'   '' KTYP5F='',I2,'' KTYP5G='',I2,'' KTYP5H='',I2)'
          WRITE(OP_STRING,FORMAT) KTYP50(nr),KTYP51(nr),KTYP52(nr),
     '      KTYP53(nr),KTYP54(nr),KTYP55(nr),KTYP55a(nr),KTYP56(nr),
     '      KTYP57(nr),KTYP58(nr),KTYP59(nr),
     '      KTYP5A(nr),KTYP5B(nr),KTYP5C(nr),KTYP5D(nr),KTYP5E(nr),
     '      KTYP5F(nr),KTYP5G(nr),KTYP5H(nr)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        FORMAT='(  '' KTYP60='',I2,'' KTYP61='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP60,KTYP61
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(  '' KTYP70='',I2,'' KTYP71='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP70,KTYP71
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(  '' KTYP90='',I2,'' KTYP91='',I2,'' KTYP92='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP90,KTYP91,KTYP92
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(  '' KTYP93(1,nr)='',9I2,'' KTYP93(2,nr)='',9I2)'
        WRITE(OP_STRING,FORMAT) (KTYP93(1,nr),nr=1,NRT),
     '                          (KTYP93(2,nr),nr=1,NRT)
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(  '' KTYP94='',I2,'' KTYPMBC='',I2)'
        WRITE(OP_STRING,FORMAT) KTYP94,KTYPMBC
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!nh_loc
      ELSE IF(TYPE(1:6).EQ.'NH_LOC') THEN
        FORMAT='(/'' NH_LOC(0,0..): '',10I2)'
        WRITE(OP_STRING,FORMAT) (NH_LOC(0,nxx),nxx=0,NXM)
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nxx=1,NX_LIST(0)
          nx1=NX_LIST(nxx)
          FORMAT='('' nh=NH_LOC(0..,'',I1,''): '',10I2)'
          WRITE(OP_STRING,FORMAT)
     '      nx1,(NH_LOC(nhx,nx1),nhx=0,NH_LOC(0,nx1))
           CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nxx
        FORMAT='('' nhx=NH_TYPE(nh,1): '',10I2)'
        WRITE(OP_STRING,FORMAT) (NH_TYPE(nh,1),nh=1,NHM)
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(''  nx=NH_TYPE(nh,2): '',10I2)'
        WRITE(OP_STRING,FORMAT) (NH_TYPE(nh,2),nh=1,NHM)
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!nj_loc
      ELSE IF(TYPE(1:6).EQ.'NJ_LOC') THEN
        FORMAT='(/'' NJL_geom='',I1,'' NJL_fibr='',I1,'
     '    //''' NJL_fiel='',I1)'
        WRITE(OP_STRING,FORMAT) NJL_geom,NJL_fibr,NJL_fiel
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' NJ_LOC(0..3,0,0): '',4I2)'
        WRITE(OP_STRING,FORMAT) (NJ_LOC(njj1,0,0),njj1=0,3)
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO nr=1,NRT
          FORMAT='( '' Region number: '',I2)'
          WRITE(OP_STRING,FORMAT) nr
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='( '' NJ_LOC(1,1..,nr):  '',2X,12I2)'
          WRITE(OP_STRING,FORMAT)
     '      (NJ_LOC(1,njj2,nr),njj2=1,NJ_LOC(1,0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='( '' NJ_LOC(2,1..,nr):  '',2X,12I2)'
          WRITE(OP_STRING,FORMAT)
     '      (NJ_LOC(2,njj2,nr),njj2=1,NJ_LOC(2,0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          FORMAT='( '' NJ_LOC(3,1..,nr):  '',2X,12I2)'
          WRITE(OP_STRING,FORMAT)
     '      (NJ_LOC(3,njj2,nr),njj2=1,NJ_LOC(3,0,nr))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nr
!nx_list
      ELSE IF(TYPE(1:7).EQ.'NX_LIST') THEN
        FORMAT='(/'' Types: NX_FIT='',I1,'' NX_OPTI='',I1,'
     '    //''' NX_SOLVE='',I1)'
        WRITE(OP_STRING,FORMAT) NX_FIT,NX_OPTI,NX_SOLVE
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(/'' Locks: NX_LOCK='',I1,'' NX_NOLOCK='',I1)'
        WRITE(OP_STRING,FORMAT) NX_LOCK,NX_NOLOCK
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' NX_LIST(0..): '',10I2)'
        WRITE(OP_STRING,FORMAT) NX_LIST(0),(NX_LIST(i),i=1,NX_LIST(0))
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' NX_CLASS(1..):'',2X,9I2)'
        WRITE(OP_STRING,FORMAT) (NX_CLASS(NX_LIST(i)),i=1,NX_LIST(0))
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' NX_TYPE(1..):'',3X,9I2)'
        WRITE(OP_STRING,FORMAT) (NX_TYPE(NX_LIST(i)),i=1,NX_LIST(0))
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' NX_LOCKS(1..):'',2X,9I2)'
        WRITE(OP_STRING,FORMAT) (NX_LOCKS(NX_LIST(i)),i=1,NX_LIST(0))
         CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!calls
      ELSE IF(TYPE(1:5).EQ.'CALLS') THEN
        FORMAT='('
     '    //'/'' CALL_ACTI='',L1,'
     '    //' '' CALL_AERO='',L1,'
     '    //' '' CALL_BASE='',L1,'
     '    //' '' CALL_COUP='',L1,'
     '    //' '' CALL_DATA='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_ACTI,CALL_AERO,CALL_BASE,CALL_COUP,CALL_DATA
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_ELEM='',L1,'
     '    //' '' CALL_ELFB='',L1,'
     '    //' '' CALL_ELFD='',L1,'
     '    //' '' CALL_EQUA='',L1,'
     '    //' '' CALL_EXPO='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_ELEM,CALL_ELFB,CALL_ELFD,CALL_EQUA,CALL_EXPO
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_FACE='',L1,'
     '    //' '' CALL_FIBR='',L1,'
     '    //' '' CALL_FIEL='',L1,'
     '    //' '' CALL_FIT ='',L1,'
     '    //' '' CALL_GRID='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_FACE,CALL_FIBR,CALL_FIEL,CALL_FIT,CALL_GRID
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_GROW='',L1,'
     '    //' '' CALL_INIT='',L1,'
     '    //' '' CALL_LINE='',L1,'
     '    //' '' CALL_MATE='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_GROW,CALL_INIT,CALL_LINE,CALL_MATE
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_MESH='',L1,'
     '    //' '' CALL_MOTI='',L1,'
     '    //' '' CALL_NODE='',L1,'
     '    //' '' CALL_OBJE='',L1,'
     '    //' '' CALL_OPTI='',L1,'
     '    //' '' CALL_SOLV='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_MESH,CALL_MOTI,CALL_NODE,CALL_OBJE,CALL_OPTI
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_SOLV='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_SOLV
       CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('
     '    //' '' CALL_DATA_FIBRE='',L1,'
     '    //' '' CALL_DATA_FIELD='',L1,'
     '    //' '' CALL_DATA_SHEET='',L1)'
        WRITE(OP_STRING,FORMAT)
     '    CALL_DATA_FIBRE,CALL_DATA_FIELD,CALL_DATA_SHEET
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!offsets
      ELSE IF(TYPE(1:7).EQ.'OFFSETS') THEN
        FORMAT='(/''  OS_ISTATE='',I7,''  OS_IWORK='',I7,'
     '        //'/'' OS_NEEDCON='',I7,''    OS_KAM='',I7,'
     '        //'/''     OS_HAM='',I7,''    OS_HSM='',I7)'
        WRITE(OP_STRING,FORMAT) OS_ISTATE,OS_IWORK,OS_NEEDCON,OS_KAM,
     '    OS_HAM,OS_HSM
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(/''  OS_CONJAC='',I7,'
     '        //'/''   OS_CONTR='',I7,'' OS_PAOPTI='',I7,'
     '        //'/''  OS_PBOPTI='',I7,''   OS_PMIN='',I7,'
     '        //'/''    OS_PMAX='',I7,''      OS_R='',I7,'
     '        //'/''   OS_RESID='',I7,'' OS_RESJAC='',I7,'
     '        //'/''    OS_WORK='',I7,''      OS_C='',I7,'
     '        //'/''       OS_F='',I7,''   OS_CJAC='',I7,'
     '        //'/''    OS_FJAC='',I7,''     OS_XC='',I7,'
     '        //'/''    L_IWORK='',I7,''    L_WORK='',I7)'
        WRITE(OP_STRING,FORMAT) OS_CONJAC,OS_CONTR,OS_PAOPTI,OS_PBOPTI,
     '    OS_PMIN,OS_PMAX,OS_R,OS_RESID,OS_RESJAC,OS_WORK,OS_C,OS_F,
     '    OS_CJAC,OS_FJAC,OS_XC,L_IWORK,L_WORK
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='(/''      OS_AM='',I7,''    OS_BLM='',I7,'
     '        //'/''     OS_BUM='',I7,''    OS_XNM='',I7,'
     '        //'/''     OS_PIM='',I7,''    OS_RCM='',I7,'
     '        //'/''  OS_RESIDM='',I7,'' OS_FGRADM='',I7,'
     '        //'/''      OS_CM='',I7,''  OS_CJACM='',I7,'
     '        //'/''      OS_ZM='',I7,''      L_ZM='',I7,'
     '        //'/''       NIWM='',I7,''      NRWM='',I7)'
        WRITE(OP_STRING,FORMAT) OS_AM,OS_BLM,OS_BUM,OS_XNM,OS_PIM,
     '    OS_RCM,OS_RESIDM,OS_FGRADM,OS_CM,OS_CJACM,OS_ZM,L_ZM,
     '    NIWM,NRWM
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
!use
      ELSE IF(TYPE(1:3).EQ.'USE') THEN
        WRITE(OP_STRING,'(/'' USE_BEM     ='',I1,T30,'
     '    //'''USE_DATA      ='',I1)') USE_BEM,USE_DATA
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_GRAPHICS   ='',I1)') USE_GRAPHICS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_GRID    ='',I1,T30,'
     '    //'''USE_LUNG      ='',I1)') USE_GRID,USE_LUNG
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_MINOS   ='',I1,T30,'
     '    //'''USE_NONLIN    ='',I1)') USE_MINOS,USE_NONLIN
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_NPSOL   ='',I1,T30,'
     '    //'''USE_OPTI      ='',I1)') USE_NPSOL,USE_OPTI
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_SPARSE  ='',I1,T30,'
     '    //'''USE_TRANSFER  ='',I1)') USE_SPARSE,USE_TRANSFER
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_NLSTRIPE  ='',I1)') USE_NLSTRIPE
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_TIME  ='',I1,T30,'
     '    //'''USE_VORONOI  ='',I1)') USE_TIME,USE_VORONOI
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_DIPOLE  ='',I1,T30,'
     '    //'''USE_CELL     ='',I1)') USE_DIPOLE,USE_CELL
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'( '' USE_GAUSS_PT_MATERIALS ='',I1)')
     '    USE_GAUSS_PT_MATERIALS
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)


CC AJPs 13-11-97
C        WRITE(OP_STRING,'( '' USE_TIME  ='',I1)') USE_TIME
C        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
CC AJPe

!versions
      ELSE IF(TYPE(1:8).EQ.'VERSIONS') THEN
        FORMAT='(/'' NTCHOR='',I2,'' NTCONT='',I2,'' NTCROS='',I2,'
     '    //''' NTDATA='',I2,'' NTFIBR='',I2,'' NTHIST='',I2,'
     '    //'/'' NTISOC='',I2,'' NTLINE='',I2,'' NTMAP ='',I2,'
     '    //''' NTPLIN='',I2,'' NTSECT='',I2,'' NTSTRA='',I2,'
     '    //'/'' NTSTRE='',I2,'' NTVELO='',I2)'
        WRITE(OP_STRING,FORMAT) NTCHOR,NTCONT,NTCROS,NTDATA,NTFIBR,
     '    NTHIST,NTISOC,NTLINE,NTMAP,NTPLIN,NTSECT,NTSTRA,NTSTRE,NTVELO
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!dipoles
      ELSE IF(TYPE(1:7).EQ.'DIPOLES') THEN
        FORMAT='(/'' NDIPOLEM ='',I4,'' NDIPTIMM ='',I4,'
     '    //''' USE_DIPOLE ='',I2)'
        WRITE(OP_STRING,FORMAT) NDIPOLEM,NDIPTIMM,USE_DIPOLE
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

!time variables
      ELSE IF(TYPE(1:8).EQ.'TIMEVARS') THEN
        FORMAT='(/'' NTIMEPOINTSM ='',I4,'' NTIMEPOINTST ='',I4)'
        WRITE(OP_STRING,FORMAT) NTIMEPOINTSM,NTIMEPOINTST
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='('' NTIMEVARSM   ='',I4,'' NTIMEVARST   ='',I4)'
        WRITE(OP_STRING,FORMAT) NTIMEVARSM,NTIMEVARST
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      CALL EXITS('OPPARA')
      RETURN
 9999 CALL ERRORS('OPPARA',ERROR)
      CALL EXITS('OPPARA')
      RETURN 1
      END

      
