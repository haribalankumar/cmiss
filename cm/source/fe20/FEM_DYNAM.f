      SUBROUTINE FEM_DYNAM(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '  CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,DIPOLE_CEN_NTIME,
     '  DIPOLE_DIR_NTIME,DIPOLE_LIST,FD,GRNGLIST,IBT,IICQS_SPATIAL,
     '  IRCQS_SPATIAL,ICQS,ICQS_SPATIAL,IDO,ILPIN,ILTIN,INP,ISC_GD,
     '  ISC_GKK,ISC_GM,ISC_GMM,ISR_GD,ISR_GKK,ISR_GM,ISR_GMM,
     '  ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,ITHRES,LD,LD_NP,LDR,
     &  LIST,LGE,LN,MAP_ART_VEIN,MXI,NAN,NBH,NBHF,NBJ,NBJF,NDADJ,NDDATA,
     &  NDDL,NDET,NDIPOLES,NDLT,NDP,NEELEM,NEL,NELIST,NELIST2,NELIST3,
     &  NENFVC,NENP,NENQ,NEP,NFF,NFFACE,NFLIST,NFLIST1,NFVC,NGAP,NGLIST,
     &  NHE,NHP,NHQ,NKB,NKEF,NKH,NKHE,NKJ,NKJE,NLCHOR,NLF,NLL,NLLINE,
     &  NLLIST,NLNO,NLQ,NLS_NDATA_CONT,NLS_NDATA_IMAG,NMBIN,NMNO,NNB,
     &  NNF,NNL,NODENVC,NODENVCB,NONL,NONM,NONY,NORD,NP1OPT,NP2OPT,
     &  NP3OPT,NPB,NPF,NPINTER,NP_INTERFACE,NPL,NPLIST1,NPLIST2,NPLIST3,
     &  NPLIST4,NPLIST5,NPNE,NPNF,NPNODE,NPNY,NPQ,NQET,NQGP,NQGP_PIVOT,
     &  NQLIST,NQNP,NQNY,NQSCNB,NQXI,NRE,NRLIST,NRLIST2,NSB,NTCOVA,
     '  NTIME_INTERP,NTIME_POINTS,NTIME_NR,NUNK,NVCB,NVCNODE,NVCNP,
     '  NVHE,NVHF,NVHP,NVJE,NVJF,NVJL,NVJP,NAQ,NW,NWP,NWQ,NXI,NXLIST,
     '  NXQ,NYNE,NYNO,NYNP,NYNQ,NYNR,NYNY,NYQNR,PAOPTY,IPIVOT,
     '  TV_BC_SET,VOLTC,Z_CONT_LIST,ACINUS,AQ,BBM,CE,CELL_CP,
     &  CELL_RCQS_VALUE,
     '  CELL_YQS_VALUE,CG,CGE,CIN,CONY,COVA,CP,CQ,CURVCORRECT,CYNO,CYNY,
     '  DET,DF,DIPOLE_CEN,DIPOLE_DIR,DL,DLL,DRDN,DRDNO,D_RE,D_RI3,D_RP,
     '  D_TG,D_ZG,DXDXIQ,DXDXIQ2,DNUDXQ,ED,EDD,EIGVAL,EIGVEC,EM,ER,ES,
     '  FEXT,GCHQ,GD,GKK,GM,GMM,GRR,GR,GUQ,LAPL,LAPLSQR,MFI,
     '  NEERR,NLS_SURF_XI,NLS_SURF_PSI,NQGW,PE,PF,PG,PGNQE,
     '  PHI,PHI_H,PHI_H_EXACT,PROPQ,RAD,RCQS,RCQS_SPATIAL,RD,RE1,RE2,
     '  REG_PARAMETER,RG,RHS,SE,SF,SP,SIGMA_PHI,SIGMA_T_BH,SQ,T_BH,
     '  T_BH_INV,THRES,TIME_VALUES,U_PHI,U_T_BH,VC,VC_INIT,VOL,VOLT,
     '  VT_PHI,VT_T_BH,WD,WDL,WG,WK1_INV,WK2_INV,WK3_INV,WK4_INV,
     '  WK5_INV,WU,XA,XAB,XB,XE,XG,XG1,XID,XIDL,XIG,XIP,
     &  XIQ,XN,XN_GRAD,XNFV,XO,XP,XR,XR_GRAD,YD,XQ,YG,YGF,YP,YQ,
     &  YQS,ZA,ZA1,ZC,Z_CONT,ZCROSSING,ZD,ZD2,ZDD,ZDL,ZE,ZF,ZG,ZG1,ZNFV,
     '  ZP,ZP1,FIX,FIXP,FIXQ,ISEG,CELL_ICQS_NAMES,
     &  CELL_RCQS_NAMES,CELL_YQS_NAMES,TIME_VARIABLE_NAMES,CSEG,END,
     &  STRING,INTWORK,REALWORK,ERROR,*)

C#### Subroutine: FEM_DYNAM
C###  Description:
C##     FEM_DYNAM Checks command string in finite element environment.

C**** LIST_1 is .true. if the 1st command list is displayed
C**** LIST_2 is .true. if the 2nd command list is displayed
C**** LIST_3 is .true. if the 3rd command list is displayed

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'cbfe01.cmn'
      INCLUDE 'cbwk01.cmn'
      INCLUDE 'cmgui00.cmn'
      INCLUDE 'comm00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lead00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'mach00.inc'
      INCLUDE 'parameters.inc'
      INCLUDE 'ptr00.cmn'
      INCLUDE 'ofst00.cmn'

!     Parameter List
      INTEGER CELL_ICQS_SPATIAL(NQIM,NQVM) !CELL_ICQS_SPATIAL(nmqi,nqv)
      INTEGER CELL_ICQS_VALUE(NQIM,NQVM)  !CELL_ICQS_VALUE(nmqi,nqv)
      INTEGER CELL_RCQS_SPATIAL(NQRM,NQVM) !CELL_RCQS_SPATIAL(nmqi,nqv)
      INTEGER CELL_YQS_SPATIAL(NIQSM,NQVM) !CELL_YQS_SPATIAL(niqs,nqv)
      INTEGER DIPOLE_CEN_NTIME(NDIPOLEM,NRM,NXM)
                                          !DIPOLE_CEN_NTIME(i,nr,nx)
      INTEGER DIPOLE_DIR_NTIME(NDIPOLEM,NRM,NXM)
                                          !DIPOLE_DIR_NTIME(i,nr,nx)
      INTEGER DIPOLE_LIST(0:NDIPOLEM)     !DIPOLE_LIST(0:ndipole)
      INTEGER FD(NDM)                     !FD(nd)
      INTEGER GRNGLIST(0:NEGM)            !GRNGLIST(0:negm)
      INTEGER IBT(3,NIM,NBFM)             !IBT(i,ni,nb)
      INTEGER IICQS_SPATIAL(0:NQISVM,NQVM) !IICQS_SPATIAL(nmqisv,nqv)
      INTEGER IRCQS_SPATIAL(0:NQRSVM,NQVM) !IRCQS_SPATIAL(nmqrsv,nqv)
      INTEGER ICQS(NQIM)                  !ICQS(nmqi)
      INTEGER ICQS_SPATIAL(NQISVM,NQM)    !ICQS_SPATIAL(nmqisv,nq)
      INTEGER IDO(NKM,NNM,0:NIM,NBFM)     !IDO(nk,nn,ni,nb)
      INTEGER ILPIN(NMM,NRM,NXM)          !ILPIN(nm,nr,nx)
      INTEGER ILTIN(NRM,NXM)              !ILTIN(nr,nx)
      INTEGER INP(NNM,NIM,NBFM)           !INP(nn,ni,nb)
      INTEGER ISC_GD(NISC_GDM)            !ISC_GD(ny)
C     INTEGER ISC_GK(NISC_GKM)            !ISC_GK(ny)
      INTEGER ISC_GKK(NISC_GKKM,NXM)      !ISC_GKK(ny,nx)
      INTEGER ISC_GM(NISC_GMM)            !ISC_GM(ny)
      INTEGER ISC_GMM(NISC_GMMM)          !ISC_GMM(no)
C     INTEGER ISC_GQ(NISC_GQM)            !ISC_GQ(ny)
      INTEGER ISR_GD(NISR_GDM)            !ISR_GD(ny)
C     INTEGER ISR_GK(NISR_GKM)            !ISR_GK(ny)
      INTEGER ISR_GKK(NISR_GKKM,NXM)      !ISR_GKK(no,nx)
      INTEGER ISR_GM(NISR_GMM)            !ISR_GM(ny)
      INTEGER ISR_GMM(NISR_GMMM)          !ISR_GMM(no)
C     INTEGER ISR_GQ(NISR_GQM)            !ISR_GQ(ny)
      INTEGER ISIZE_MFI(3,NSSM)           !ISIZE_MFI(3,nss)
      INTEGER ISIZE_PHI(2)                !ISIZE_PHI(2)
      INTEGER ISIZE_PHIH(2)               !ISIZE_PHIH(2)
      INTEGER ISIZE_TBH(2)                !ISIZE_TBH(2)
      INTEGER ITHRES(3,NGM,NEM)           !ITHRES(i,ng,ne)
      INTEGER LD(NDM)                     !LD(nd)
      INTEGER LD_NP(NDM)                  !LD(nd)
      INTEGER LDR(0:NDM)                  !LDR(nd)
      INTEGER LIST(0:NLISTM)              !LIST(nl)
      INTEGER LGE(NHM*NSM,NRCM)           !LGE(nh*ns,nrc)
      INTEGER LN(0:NEM)                   !LN(ne)
      INTEGER MAP_ART_VEIN(0:NDM,NRM)  !MAP_ART_VEIN(nd,nr)
      INTEGER MXI(2,NEM)                  !MXI(i,ne)
      INTEGER NAN(NIM,NAM,NBFM)           !NAN(ni,na,nb)
      INTEGER NBH(NHM,NCM,NEM)            !NBH(nh,nc,ne)
      INTEGER NBHF(NHM,NCM,NFM)           !NBHF(nh,nc,nf)
      INTEGER NBJ(NJM,NEM)                !NBJ(nj,ne)
      INTEGER NBJF(NJM,NFM)               !NBJF(nj,nf)
      INTEGER NDADJ(6,NDM)                !NDADJ(i,nd)
      INTEGER NDDATA(0:NDM,0:NRM)         !NDDATA(nd,0:nr)
      INTEGER NDDL(NEM,NDEM)              !NDDL(ne,nde)
      INTEGER NDET(NBFM,0:NNM)            !NDET(nbf,nn)
      INTEGER NDIPOLES(NRM,NXM)           !NDIPOLES(nr,nx)
      INTEGER NDLT(NEM)                   !NDLT(ne)
      INTEGER NDP(NDM)                    !NDP(nd)
      INTEGER NEELEM(0:NE_R_M,0:NRM)      !NEELEM(ne,nr)
      INTEGER NEL(0:NELM,NLM)             !NEL(i,nl)
      INTEGER NELIST(0:NEM)               !NELIST(ne)
      INTEGER NELIST2(0:NEM)              !NELIST2(ne)
      INTEGER NELIST3(0:NEM)              !NELIST3(ne)
      INTEGER NENFVC(0:NFVCM,NFVM)        !NENFVC(noelem,nfl)
      INTEGER NENP(NPM,0:NEPM,0:NRM)      !NENP(np,nep,nr)
      INTEGER NENQ(0:8,NQM)               !NENQ(ne,nq)
      INTEGER NEP(NPM)                    !NEP(np)
      INTEGER NFF(6,NEM)                  !NFF(i,ne)
      INTEGER NFFACE(0:NF_R_M,NRM)        !NFFACE(noface,nr)
      INTEGER NFLIST(0:NFM)               !NFLIST(nf)
      INTEGER NFLIST1(0:NFM)              !NFLIST1(nf)
      INTEGER NFVC(2,0:NFVCM,NVCM)        !NFVC(2,nfvl,nvc)
      INTEGER NGAP(NIM,NBM)               !NGAP(ni,nb)
      INTEGER NGLIST(0:NGM)               !NGLIST(ng)
C      INTEGER NGNQE(NGM,0:NQEM,NQSCM)     !NGNQE(ng,nqq,nqsc)
      INTEGER NHE(NEM,NXM)                !NHE(ne,nx)
      INTEGER NHP(NPM,0:NRM,NXM)          !NHP(np,nr,nx)
      INTEGER NHQ(NRM,NXM)                !NHQ(nr,nx)
      INTEGER NKB(2,2,2,NNM,NBFM)         !NKB(ido1,ido2,ido3,nn,nb)
      INTEGER NKEF(0:4,16,6,NBFM)         !NKEF(i,nn,nfe,nb)
      INTEGER NKH(NHM,NPM,NCM,0:NRM)      !NKH(nh,np,nc,nr)
      INTEGER NKHE(NKM,NNM,NHM,NEM)       !NKHE(nk,nn,nh,ne)
      INTEGER NKJ(NJM,NPM)                !NKJ(nj,np)
      INTEGER NKJE(NKM,NNM,NJM,NEM)       !NKJE(nk,nn,nj,ne)
      INTEGER NLCHOR(0:10,NRM)            !NLCHOR(i,nr)
      INTEGER NLF(4,NFM)                  !NLF(i,nf)
      INTEGER NLL(12,NEM)                 !NLL(i,ne)
      INTEGER NLLINE(0:NL_R_M,0:NRM)      !NLLINE(nl,nr)
      INTEGER NLLIST(0:NLM)               !NLLIST(nl)
      INTEGER NLNO(NOPM,NXM)              !NLNO(nop,nx)
      INTEGER NLQ(NQM)                    !NLQ(nq)
      INTEGER NLS_NDATA_CONT(NDM)         !NLS_NDATA_CONT(nd)
      INTEGER NLS_NDATA_IMAG(NDM)         !NLS_NDATA_IMAG(nd)
      INTEGER NMBIN(NMM,NRM,NXM)          !NMBIN(nm,nr,nx)
      INTEGER NMNO(1:2,0:NOPM,NXM)        !NMNO(1:2,nop,nx)
      INTEGER NNB(4,4,4,NBFM)             !NNB(inp1,inp2,inp3,nb)
      INTEGER NNF(0:17,6,NBFM)            !NNF(i,nfe,nb)
      INTEGER NNL(0:4,12,NBFM)            !NNL(i,j,nb)
      INTEGER NODENVC(NVCM)               !NODENVC(nvc)
      INTEGER NODENVCB(NVCBM)             !NODENVCB(bnvc)
      INTEGER NONL(NLM,NXM)               !NONL(nl,nx)
      INTEGER NONM(NMM,NPM,NXM)           !NONM(nm,np,nx)
      INTEGER NONY(0:NOYM,NYM,NRCM,0:NRM,NXM) !NONY(noy,ny,nrc,nr,nx)
      INTEGER NORD(5,NE_R_M)              !NORD(i,ne)
      INTEGER NP1OPT(NOPM)                !NP1OPT(nop)
      INTEGER NP2OPT(NOPM)                !NP2OPT(nop)
      INTEGER NP3OPT(NOPM)                !NP3OPT(nop)
      INTEGER NPB(0:NP_R_M,5)             !NPB(np,5)
      INTEGER NPF(9,NFM)                  !NPF(i,nf)
      INTEGER NPINTER(3,0:20)             !NPINTER(i,j)
      INTEGER NP_INTERFACE(0:NPM,0:3)     !NP_INTERFACE(np,i)
      INTEGER NPL(5,0:3,NLM)              !NPL(i,0:3,nl)
      INTEGER NPLIST1(0:NPM)              !NPLIST1(np)
      INTEGER NPLIST2(0:NPM)              !NPLIST2(np)
      INTEGER NPLIST3(0:NPM)              !NPLIST3(np)
      INTEGER NPLIST4(0:NPM)              !NPLIST4(np)
      INTEGER NPLIST5(0:NPM)              !NPLIST5(np)
      INTEGER NPNE(NNM,NBFM,NEM)          !NPNE(nn,nbf,ne)
      INTEGER NPNF(NNM,NBFM)              !NPNF(nn,nbf)
      INTEGER NPNODE(0:NP_R_M,0:NRM)      !NPNODE(np,nr)
      INTEGER NPNY(0:6,NYM,0:NRCM,NXM)    !NPNY(i,ny,nrc,nx)
      INTEGER NPQ(NQM)                    !NPQ(nq)
      INTEGER NQET(NQSCM)                 !NQET(nqsc)
      INTEGER NQGP(0:NQGM,NQM)            !NQGP(i,nq)
      INTEGER NQGP_PIVOT(NQGM,NQM)        !NQGP_PIVOT(i,nq)
      INTEGER NQLIST(0:NQM)               !NQLIST(nq)
C      INTEGER NQNE(NEQM,NQEM)              !NQNE(ne,nqe)
      INTEGER NQNP(NPM)                   !NQNP(np)
      INTEGER NQNY(2,NYM,0:NRCM,NXM)      !NQNY(2,ny,nx)
C      INTEGER NQS(NEQM)                    !NQS(ne)
      INTEGER NQSCNB(NQSCM)               !NQSCNB(nqsc)
      INTEGER NQXI(0:NIM,NQSCM)           !NQXI(ni,nqs)
      INTEGER NRE(NEM)                    !NRE(ne)
      INTEGER NRLIST(0:NRM)               !NRLIST(nr)
      INTEGER NRLIST2(0:NRM)              !NRLIST2(nr)
      INTEGER NSB(NKM,NNM,NBFM)           !NSB(nk,nn,nb)
      INTEGER NTCOVA(NEM)                 !NTCOVA(ne)
      INTEGER NTIME_INTERP(NTIMEVARSM)    !NTIME_INTERP(ntv)
      INTEGER NTIME_POINTS(NTIMEVARSM)    !NTIME_POINTS(ntv)
      INTEGER NTIME_NR(0:NTIMEVARSM,NRM)  !NTIME_NR(ntv,nr)
      INTEGER NUNK(NKM,NJM,NPM)           !NUNK(nk,nj,np)
      INTEGER NVCB(-1:3,NVCBM)            !NVCB(-1:3,nvc)
      INTEGER NVCNODE(2,NP_R_M)           !NVCNODE(2,nonode)
      INTEGER NVCNP(NP_R_M)               !NVCNP(nonode)
      INTEGER NVHE(NNM,NBFM,NHM,NEM)      !NVHE(nn,nb,nh,ne)
      INTEGER NVHF(NNM,NBFM,NHM)          !NVHF(nn,nb,nh)
      INTEGER NVHP(NHM,NPM,NCM,0:NRM)     !NVHP(nh,np,nc,nr)
      INTEGER NVJE(NNM,NBFM,NJM,NEM)      !NVJE(nn,nb,nj,ne)
      INTEGER NVJF(NNM,NBFM,NJM)          !NVJF(nn,nb,nj)
      INTEGER NVJL(4,NJM,NLM)             !NVJL(nn,nj,nl)
      INTEGER NVJP(NJM,NPM)               !NVJP(nj,np)
      INTEGER NAQ(NQM,NAM)                !NAQ(nq,na)
      INTEGER NW(NEM,3,NXM)               !NW(ne,i,nx)
      INTEGER NWP(NPM,2)                  !NWP(np,i)
      INTEGER NWQ(8,0:NQM,NAM)            !NWQ(8,nq,na)
      INTEGER NXI(-NIM:NIM,0:NEIM,0:NEM)  !NXI(-ni:ni,1,ne)
      INTEGER NXLIST(0:NXM)               !NXLIST(nx)
      INTEGER NXQ(-NIM:NIM,0:4,0:NQM,NAM) !NXQ(-ni:ni,i,nq,na)
      INTEGER NYNE(NAM,NHM,0:NRCM,NCM,NEM) !NYNE(na,nh,nrc,nc,ne)
      INTEGER NYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM)
                                           !NYNO(nyo,noop,nrc,nr,nx)
      INTEGER NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
                                           !NYNP(nk,nv,nh,np,nrc,nc,nr)
      INTEGER NYNQ(NHM,NQM,0:NRCM,NXM)      !NYNQ(nh,nq,nx)
      INTEGER NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
                                          !NYNR(ny,nrc,nc,nrnx)
      INTEGER NYNY(0:NYYM,NYM,NRM,NXM)    !NYNY(0:nyy,ny,nr,nx)
      INTEGER NYQNR(0:NYQM,0:NRCM,NCM,0:NRM,NXM)
                                          !NYQNR(ny,nrc,nc,nrnx)
      INTEGER PAOPTY(NOPM)                !PAOPTY(no?)
c      INTEGER PN_TEMP(NE_R_M)              !PN_TEMP(ne)
      INTEGER IPIVOT(NOM)
      INTEGER TV_BC_SET(0:NIQSM,0:NQM)     !TV_BC_SET(niqs,nq)
      INTEGER VOLTC(NBFM)                  !VOLTC(nb)
C*** 19/02/08 JHC increased array size
C      INTEGER Z_CONT_LIST(NDM,2,4)           !Z_CONT_LIST(nd,i,j)
      INTEGER Z_CONT_LIST(NDM,2,7)           !Z_CONT_LIST(nd,i,j)


      REAL*8 AQ(NMAQM,NQM)                     !AQ(maq,nq)
      REAL*8 ACINUS(4,NEM)
      REAL*8 BBM(2,NEM)
      REAL*8 CE(NMM,NEM,NXM)                   !CE(nm,ne,nx)
      REAL*8 CELL_CP(NMQM,NPM)                 !CELL_CP(nmqm,np)
      REAL*8 CELL_RCQS_VALUE(NQRM,NQVM)       !CELL_RCQS_VALUE(nmqr,nqv)
      REAL*8 CELL_YQS_VALUE(NIQSM,NQVM)        !CELL_YQS_VALUE(niqs,nqv)
      REAL*8 CG(NMM,NGM)                       !CG(nm,ng)
      REAL*8 CGE(NMM,NGM,NEM,NXM)              !CGE(nm,ng,ne,nx)
      REAL*8 CIN(NMM,0:NGM,NNEPM)              !CIN(nm,0:ng,nnep)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM,NXM)   !CONY(noy,ny,nrc,nr,nx)
      REAL*8 COVA(NEM,50)                      !COVA(ne,i)
      REAL*8 CP(NMM,NPM,NXM)                   !CP(nm,np,nx)
      REAL*8 CQ(NMM,NQM,NXM)                   !CQ(nm,nq,nx)
      REAL*8 CURVCORRECT(2,2,NNM,NEM)          !CURVCORRECT(2,2,nn,ne)
      REAL*8 CYNO(0:NYOM,NOOPM,NRCM,0:NRM,NXM) !CYNO(nyo,noop,nrc,nr,nx)
      REAL*8 CYNY(0:NYYM,NYM,NRM,NXM)          !CYNY(0:nyy,ny,nr,nx)
      REAL*8 DET(NBFM,0:NNM,NGM,6)             !DET(nbf,nn,ng,i)
      REAL*8 DF(NFM)                           !DF(nf)
      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM)
                                               !DIPOLE_CEN(i,j,k,nr,nx)
      REAL*8 DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM,NXM)
                                               !DIPOLE_DIR(i,j,k,nr,nx)
      REAL*8 DL(3,NLM)                         !DL(i,nl)
      REAL*8 DLL(3,NLM)                        !DLL(i,nl)
      REAL*8 DRDN(NGM)                         !DRDN(ng)
      REAL*8 DRDNO(NGM,NKM)                    !DRDNO(ng,nk)
      REAL*8 D_RE(NSM,NHM,NOPM)                !D_RE(ns,nh,nop)
      REAL*8 D_RI3(NHM*NSM)                    !D_RI3(nhs)
      REAL*8 D_RP(NYM,NYM)                     !D_RP(ny,ny)
      REAL*8 D_TG(3,3,NHM*NSM)                 !D_TG(i,j,nhs)
      REAL*8 D_ZG(NHM,NUM,NHM*NSM)             !D_ZG(nh,nu,nhs)
      REAL*8 DXDXIQ(3,3,NQM)                   !DXDXIQ(i,j,nq)
      REAL*8 DXDXIQ2(3,3,NQM)                  !DXDXIQ2(i,j,nq)
      REAL*8 DNUDXQ(3,3,NQM)                   !DNUDXQ(i,j,nq)
      REAL*8 ED(NHM*NSM,NHM*NSM)               !ED(nhs,nhs)
      REAL*8 EDD(NDM)                          !EDD(nd)
      REAL*8 EIGVAL(NTM,2)                     !EIGVAL(nt,i)
      REAL*8 EIGVEC(NOM,NTM,2)                 !EIGVEC(no,nt,i)
      REAL*8 EM(NHM*NSM,NHM*NSM)               !EM(nhs,nhs)
      REAL*8 ER(NHM*NSM)                       !ER(nhs)
      REAL*8 ES(NHM*NSM,NHM*NSM)               !ES(nhs,nhs)
      REAL*8 FEXT(NIFEXTM,NGM,NEM)             !FEXT(i,ng,ne)
      REAL*8 GCHQ(3,NQM)                       !GCHQ(i,nq)
      REAL*8 GD(NZ_GD_M)                       !GD(nz)
C      REAL*8 GK(NZ_GK_M)                       !GK(nz)
      REAL*8 GKK(NZ_GKK_M,NXM)                 !GKK(nz,nx)
      REAL*8 GM(NZ_GM_M)                       !GM(nz)
      REAL*8 GMM(NZ_GMM_M)                     !GMM(nz)
C      REAL*8 GQ(NZ_GQ_M)                       !GQ(nz)
      REAL*8 GR(NYROWM)                        !GR(nyrow)
      REAL*8 GRR(NOM)                          !GRR(no)
      REAL*8 GUQ(3,3,NQM)                      !GUQ(i,j,nq)
      REAL*8 LAPL(NY_TRANSFER_M,NY_TRANSFER_M) !LAPL(np,np2)
      REAL*8 LAPLSQR(NY_TRANSFER_M,NY_TRANSFER_M) !LAPLSQR(np,np2)
      REAL*8 MFI(NDM,NTSM,3,NSSM)              !MFI(nd,nts,3,nss)
      REAL*8 NEERR(NEM,3)                      !NEERR(ne,1..3)
      REAL*8 NLS_SURF_XI(3,26,26,4)            !NLS_SURF_XI(,,,)
      REAL*8 NLS_SURF_PSI(16,26,26,3,4)        !NLS_SURF_PSI(,,,,)
      REAL*8 NQGW(NQGM,NQM)                    !NQGW(i,nq)
      REAL*8 PE(2,NEM)                         !PE(i,ne)
      REAL*8 PF(2,NEM)                         !PF(i,ne)
      REAL*8 PG(NSM,NUM,NGM,NBM)               !PG(ns,nu,ng,nb)
      REAL*8 PGNQE(NGM,NQEM,NQSCM)             !PGNQE(ng,nqe,nqsc)
      REAL*8 PHI(NY_TRANSFER_M,NTSM)           !PHI(nytr,nts)
      REAL*8 PHI_H(NY_TRANSFER_M,NTSM)         !PHI_H(nytr,nts)
      REAL*8 PHI_H_EXACT(NY_TRANSFER_M,NTSM)   !PHI_H_EXACT(nytr,nts)
      REAL*8 PROPQ(3,3,4,2,NQM,NXM)            !PROPQ(,,,,nq,nx)
      REAL*8 RAD(NGM)                          !RAD(ng)
      REAL*8 RCQS(NQRM)                        !RCQS(nmqr)
      REAL*8 RCQS_SPATIAL(NQRSVM,NQM)          !RCQS_SPATIAL(nmqrsv,nq)
      REAL*8 RD(NGM)                           !RD(ng)
      REAL*8 RE1(NSM,NHM)                      !RE1(ns,nh)
      REAL*8 RE2(NSM,NHM)                      !RE2(ns,nh)
      REAL*8 REG_PARAMETER(0:NTSM)             !REG_PARAMETER(nts)
      REAL*8 RG(NGM)                           !RG(ng)
      REAL*8 RHS(NQM)                          !RHS(nq)
      REAL*8 SE(NSM,NBFM,NEM)                  !SE(ns,nbf,ne)
      REAL*8 SF(NSM,NBFM)                      !SF(ns,nbf)
      REAL*8 SP(NKM,NBFM,NPM)                  !SP(nk,nbf,np)
      REAL*8 SIGMA_PHI(NY_TRANSFER_M)          !SIGMA_PHI(nytr)
      REAL*8 SIGMA_T_BH(NY_TRANSFER_M)         !SIGMA_PHI(nytr)
      REAL*8 SQ(NDM)                           !SQ(nd)
      REAL*8 T_BH(NY_TRANSFER_M,NY_TRANSFER_M)     !T_BH(nytr,nytr)
      REAL*8 T_BH_INV(NY_TRANSFER_M,NY_TRANSFER_M) !T_BH_INV(nytr,nytr)
      REAL*8 THRES(3,NGM,NEM)                  !THRES(i,ng,ne)
      REAL*8 TIME_VALUES(2,0:NTIMEPOINTSM+1,NTIMEVARSM)
      REAL*8 U_PHI(NY_TRANSFER_M,NY_TRANSFER_M) !U_PHI(nytr,nytr)
      REAL*8 U_T_BH(NY_TRANSFER_M,NY_TRANSFER_M) !U_T_BH(nytr,nytr)
      REAL*8 VC(0:NVCM)                        !VC(nvc)
      REAL*8 VC_INIT(2,NVCM)                   !VC_INIT(i,nvc)
      REAL*8 VOL(NBFM)                         !VOL(nb)
      REAL*8 VOLT(NBFM)                        !VOLT(nb)
      REAL*8 VT_PHI(NTSM,NTSM)                 !VT_PHI(nts,nts)
      REAL*8 VT_T_BH(NY_TRANSFER_M,NY_TRANSFER_M) !VT_T_BH(nytr,nytr)
      REAL*8 WD(NJM,NDM)                       !WD(nj,nd)
      REAL*8 WDL(NHM,NDEM)                     !WDL(nh,nde)
      REAL*8 WG(NGM,NBM)                       !WG(ng,nb)
      REAL*8 WK1_INV(NY_TRANSFER_M,NY_TRANSFER_M)
      REAL*8 WK2_INV(NY_TRANSFER_M,NY_TRANSFER_M)
      REAL*8 WK3_INV(NY_TRANSFER_M,NY_TRANSFER_M)
      REAL*8 WK4_INV(NY_TRANSFER_M)
      REAL*8 WK5_INV(NTSM)
      REAL*8 WU(0:NUM+1,NEM)                   !WU(0:10,ne)
      REAL*8 XA(NAM,NJM,NEM)                   !XA(na,nj,ne)
      REAL*8 XAB(NORM,NEM)                     !XAB(i,ne)
      REAL*8 XB(2,NJM,NLM)                     !XB(i,nj,nl)
      REAL*8 XE(NSM,NJM)                       !XE(ns,nj)
      REAL*8 XG(NJM,NUM)                       !XG(nj,nu)
      REAL*8 XG1(NJM,NUM,NGM)                  !XG1(nj,nu,ng)
      REAL*8 XID(NIM,NDM)                      !XID(ni,nd)
      REAL*8 XIDL(NIM,NDEM)                    !XIDL(ni,nde)
      REAL*8 XIG(NIM,NGM,NBM)                  !XIG(ni,ng,nb)
      REAL*8 XIP(NIM,NPM)                      !XIP(ni,np)
      REAL*8 XIQ(NIM,NQM)                      !XIQ(ni,nq)
      REAL*8 XN(NJM,NGM)                       !XN(nj,ng)
      REAL*8 XN_GRAD(NJM,NGM)                  !XN_GRAD(nj,ng)
      REAL*8 XNFV(-(NJM+1):NJM,NFVM)           !XNFV(:,nfv)
      REAL*8 XO(NOM,NXM)                       !XO(no,nx)
      REAL*8 XP(NKM,NVM,NJM,NPM)               !XP(nk,nv,nj,np)
      REAL*8 XR(NJM,NGM)                       !XR(nj,ng)
      REAL*8 XR_GRAD(NJM,NGM)                  !XR_GRAD(nj,ng)
      REAL*8 YD(NHM)                           !YD(nh)
      REAL*8 XQ(NJM,NQM)                       !XQ(nj,nq)
      REAL*8 YG(NIYGM,NGM,NEM)                 !YG(niyg,ng,ne)
      REAL*8 YGF(NIYGFM,NGFM,NFM)              !YGF(niygf,ng,nf)
      REAL*8 YP(NYM,NIYM,NXM)                  !YP(ny,niy,nx)
      REAL*8 YQ(NYQM,NIQM,NAM,NXM)             !YQ(nyq,niq,na,nx)
      REAL*8 YQS(NIQSM,NQM)                    !YQS(niqs,nq)
      REAL*8 ZA(NAM,NHM,NCM,NEM)               !ZA(na,nh,nc,ne)
      REAL*8 ZA1(NAM,NHM,NCM,NEM)              !ZA1(na,nh,nc,ne)
      REAL*8 ZC(NJM,NEM)                       !ZC(nj,ne)
C*** 19/02/08 JHC increased array size
C      REAL*8 Z_CONT(NDM,2,25)                  !Z_CONT(nd,i,j)
      REAL*8 Z_CONT(NDM,2,67)                  !Z_CONT(nd,i,j)
      REAL*8 ZCROSSING(NY_TRANSFER_M,NTSM)     !ZCROSSING(nytr,nts)
      REAL*8 ZD(NJM,NDM)                       !ZD(nj,nd)
      REAL*8 ZD2(NJM,NDM)                      !ZD2(nj,nd)
      REAL*8 ZDD(NJM,NDM)                      !ZDD(nj,nd)
      REAL*8 ZDL(NHM,NDEM)                     !ZDL(nh,nde)
      REAL*8 ZE(NSM,NHM)                       !ZE(ns,nh)
      REAL*8 ZF(NSM,NHM)                       !ZF(ns,nh)
      REAL*8 ZG(NHM,NUM)                       !ZG(nh,nu)
      REAL*8 ZG1(NHM,NUM)                      !ZG1(nh,nu)
      REAL*8 ZNFV(NFVM)
      REAL*8 ZP(NKM,NVM,NHM,NPM,NCM)           !ZP(nk,nv,nh,np,nc)
      REAL*8 ZP1(NKM,NVM,NHM,NPM,NCM)          !ZP1(nk,nv,nh,np,nc)

      LOGICAL FIX(NYM,NIYFIXM,NXM)             !FIX(ny,niy_fix,nx)
      LOGICAL FIXP(2,NEM)                      !FIXP(i,ne)
      LOGICAL FIXQ(NYQM,NIYFIXM,NXM)           !FIX(nyq,niy_fix,nx)

      INTEGER ISEG(*),INTWORK(*),VORO_SIZE
      INTEGER*4 BC_POINTS_PTR,BRANCH_PTR,CALCULATED_PTR,
     '  CONECT_PTR,IPIV_PTR,S_PTR,VORO_PTR,YP1_PTR
      REAL*8 REALWORK(*)
      CHARACTER
     '  CELL_ICQS_NAMES(NQIM,NQVM)*(*),        !CELL_ICQS_NAMES(nqi,nqv)
     '  CELL_RCQS_NAMES(NQRM,NQVM)*(*),        !CELL_RCQS_NAMES(nqr,nqv)
     '  CELL_YQS_NAMES(NIQSM,NQVM)*(*),        !CELL_YQS_NAMES(niqs,nqv)
     '  TIME_VARIABLE_NAMES(NTIMEVARSM)*(*),
     '  CSEG(*)*(*),ERROR*(*),STRING*(MXCH)
      LOGICAL END

!     Local Variables
      INTEGER ISFIPR,iw,MAXNITW,
     '  noch2,noch3,NTCH2,NTCH3,NVARINDEX
      REAL*8 EVALTIME,VALUE
      LOGICAL ABBREV,CONTINUE2,CONTINUE3,DOCUMENT
      LOGICAL LIST_2,LIST_3
      CHARACTER DATYPE(20)*11,OPTION_3(99)*80,OPTION_2(99)*80

      noch2=0
      noch3=0

      CALL ENTERS('FEM_DYNAM',*9999)

      IF(CMGUI_LINK) THEN
C cpb 4/3/97 temporary disable for class version
C db 7/10/97 temporary disable disabled
C cs 16/11/99 temporary disable disable disabled
C       Check for any information updates from CMGUI
C        CALL CMGUI_LINK_GET_DATA(IBT,IDO,INP,LD,NBJ,NBJF,NEELEM,
C     '    NEL,NENP,NFF,NFFACE,NKJE,NKEF,NKH,NKJ,NLF,NLL,NLLINE,
C     '    NNF,NNL,NPF,NPL,NPNE,NPNF,NPNODE,NRE,NVHP,NVJE,NVJF,NVJL,NVJP,
C     '    NYNP,DF,DL,PG,RG,SE,SF,WD,WG,XA,XE,XG,XID,XP,YP,ZD,
C     '    FIX,ERROR,*9999)
      ENDIF

      IF(ABBREV(CO(2),'DOCUMENT',3)) THEN
C#### Command: FEM document
C###  Description:
C###    Writes out CMISS command help for all commands.
C###    The help produced is the same as if the '?' qualifier was used.
C       LIST_1=.FALSE.
        LIST_2=.TRUE.
        LIST_3=.TRUE.
        COMMAND=.FALSE.
        DOCUMENT=.TRUE.
        noch2=0
        noch3=0
      ELSE IF(ABBREV(CO(2),'?',1)) THEN
C       LIST_1=.FALSE.
        LIST_2=.TRUE.
        LIST_3=.FALSE.
        COMMAND=.FALSE.
        DOCUMENT=.FALSE.
      ELSE IF(ABBREV(CO(3),'?',1)) THEN
C       LIST_1=.FALSE.
        LIST_2=.FALSE.
        LIST_3=.TRUE.
        COMMAND=.FALSE.
        DOCUMENT=.FALSE.
      ELSE
C       LIST_1=.FALSE.
        LIST_2=.FALSE.
        LIST_3=.FALSE.
        COMMAND=.TRUE.
        DOCUMENT=.FALSE.
      ENDIF
      noco=2

      CONTINUE2=.TRUE.
      DO WHILE (CONTINUE2)
        IF(LIST_2) THEN
          OPTION_2( 1)='Add'
          OPTION_2( 2)='Apply'
          OPTION_2( 3)='Cancel'
          OPTION_2( 4)='Change'
          OPTION_2( 5)='Check'
          OPTION_2( 6)='Close'
          OPTION_2( 7)='Combine'
          OPTION_2( 8)='Compare'
          OPTION_2( 9)='Convert'
          OPTION_2(10)='Copy'
C Temporary ? AJP
          OPTION_2(11)='Crash'
          OPTION_2(12)='Define'
          OPTION_2(13)='Display'
          OPTION_2(14)='Document'
          OPTION_2(15)='Draw'
          OPTION_2(16)='Duplicate'
          OPTION_2(17)='Evaluate'
          OPTION_2(18)='Export'
          OPTION_2(19)='Fit'
          OPTION_2(20)='Group'
C MPN 26Mar2002: changed >print to >fem gxprint
          OPTION_2(21)='Gxprint'
          OPTION_2(22)='Help'
C MPN 26Mar2002: fixed bug ('Hide' originally omitted here)
          OPTION_2(23)='Hide'
          OPTION_2(24)='Import'
          OPTION_2(25)='Inquire'
          OPTION_2(26)='List'
          OPTION_2(27)='Open'
          OPTION_2(28)='Quit'
          OPTION_2(29)='Read'
          OPTION_2(30)='Reallocate'
          OPTION_2(31)='Refine'
          OPTION_2(32)='Renumber'
C LKC 6-OCT-1999 Archived
C          OPTION_2(31)='Send'
          OPTION_2(33)='Shape'
          OPTION_2(34)='Show'
          OPTION_2(35)='Solve'
          OPTION_2(36)='Step'
          OPTION_2(37)='Track'
          OPTION_2(38)='Truncate'
          OPTION_2(39)='Update'
          OPTION_2(40)='Write'
          OPTION_2(41)='Exit'
          NTCH2=41
          IF(DOCUMENT) THEN
            noch2=noch2+1
            CO(2)=OPTION_2(noch2)
            CO(3)='?'
          ELSE
            CALL LIST_COMMANDS(2,NTCH2,OPTION_2,ERROR,*9999)
          ENDIF
        ENDIF !list_2

        IF(ABBREV(CO(2),'ADD',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Field'
              OPTION_3( 2)='Signal'
              NTCH3=2
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'FIELD',1)) THEN
              CALL ADDFIEL(NKH,NKJ,NPNODE,NRLIST,NVHP,NVJP,XP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SIGNAL',1)) THEN
              CALL ADDSIG(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3
        ELSE IF(ABBREV(CO(2),'APPLY',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Noise'
              OPTION_3( 2)='Reference'
              OPTION_3( 3)='Transfer'
              OPTION_3( 4)='Inverse'
              OPTION_3( 5)='Exit'
              NTCH3=5
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'NOISE',2)) THEN
              CALL APNOIS(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REFERENCE',2)) THEN
              CALL APREFE(LD,NBJ,NDDATA,NHQ,NP_INTERFACE,NPNY,NQNY,
     '          NRLIST,NRLIST2,
     '          NXLIST,NYNP,NYNR,NYQNR,WD,XID,YP,YQ,YQS,ZD,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TRANSFER',3)) THEN
              CALL APTRSF(ISIZE_PHI,ISIZE_PHIH,
     '          ISIZE_TBH,NHP,NHQ,NKH,
     '          NPLIST3,NPLIST4,NPLIST5,NPNY,NQNY,NVHP,NXLIST,NYNP,NYNR,
     '          NYQNR,PHI_H,T_BH,YP,YQ,YQS,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INVERSE',2)) THEN
              S_PTR=0 
              CALL ALLOCATE_MEMORY(NY_TRANSFER_M,0,3,S_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL APINVE(ISIZE_TBH,NHP,NHQ,NKH,NPLIST3,NPLIST4,NPNY,
     '          NQNY,NVHP,NXLIST,NYNP,NYNR,NYQNR,%VAL(S_PTR),
     '          STRING,T_BH,T_BH_INV,WK1_INV,WK4_INV,YP,YQ,YQS,
     '          ERROR,*9999)
              CALL FREE_MEMORY(S_PTR,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'CANCEL',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Alignment'
              OPTION_3( 2)='Axes'
              OPTION_3( 3)='Bases'
              OPTION_3( 4)='Boundary'
              OPTION_3( 5)='Clock'
              OPTION_3( 6)='Contours'
              OPTION_3( 7)='Cross-section'
              OPTION_3( 8)='Data'
              OPTION_3( 9)='Elements'
              OPTION_3(10)='Faces'
              OPTION_3(11)='Fibres'
              OPTION_3(12)='Field'
              OPTION_3(13)='Gauss'
              OPTION_3(14)='Gradient'
              OPTION_3(15)='Grid'
              OPTION_3(16)='History'
              OPTION_3(17)='Increments'
              OPTION_3(18)='Isochrones'
              OPTION_3(19)='Lines'
              OPTION_3(20)='Map'
              OPTION_3(21)='Materials'
              OPTION_3(22)='Nodes'
              OPTION_3(24)='Object'
              OPTION_3(25)='Plot'
              OPTION_3(26)='Polyline'
              OPTION_3(27)='Profile'
              OPTION_3(28)='Reactions'
              OPTION_3(29)='Residuals'
              OPTION_3(30)='Rule'
              OPTION_3(31)='Section'
              OPTION_3(32)='Sheets'
              OPTION_3(33)='Strain'
              OPTION_3(34)='Stress'
              OPTION_3(35)='Exit'
              NTCH3=35
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF  
            IF(ABBREV(CO(3),'ALIGNMENT',2)) THEN
              CALL CAALIG(%VAL(ISALIG_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'AXES',2)) THEN
              CALL CAAXES(%VAL(ISAXES_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BASES',2)) THEN
              CALL CABASE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BOUNDARY',2)) THEN
C             CALL CABOUN(%VAL(ISBOUN_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CLOCK',2)) THEN
              CALL CACLOC(%VAL(ISCLOC_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONTOURS',2)) THEN
              CALL CACONT(%VAL(ISCONO_PTR),%VAL(ISCONT_PTR),ISEG,NEELEM,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CROSS-SECTION',2)) THEN
              CALL CACROS(%VAL(ISCROS_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',1)) THEN
              CALL CADATA(%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,LD,NDDL,NDLT,NDP,
     &          NEELEM,NELIST,CSEG,EDD,SQ,STRING,WD,XID,ZD,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL CAELEM(ISEG,%VAL(ISELNO_PTR),NEELEM,NELIST,NRLIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FACES',2)) THEN
              CALL CAFACE(ISEG,%VAL(ISFACE_PTR),%VAL(ISFANO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRES',3)) THEN
              CALL CAFIBR(ISEG,%VAL(ISFIBR_PTR),NEELEM,NRLIST,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL CAFIEL(ISEG,%VAL(ISFIEL_PTR),NEELEM,NRLIST,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL CAGAUS(ISEG,%VAL(ISGAUS_PTR),NBJ,NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'GRADIENTS',3)) THEN
              CALL CAGRAD(ISEG,%VAL(ISGRAD_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
              CALL CAGRID(ISEG,%VAL(ISGRID_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HISTORY',2)) THEN
              CALL CAHIST(ISEG,%VAL(ISHIST_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',2)) THEN
              CALL CAINCR(ISEG,%VAL(ISINCR_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ISOCHRONES',2)) THEN
              CALL CAISOC(ISEG,%VAL(ISISOC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',1)) THEN
              CALL CALINE(ISEG,%VAL(ISLINE_PTR),%VAL(ISLINO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAP',3)) THEN
              CALL CAMAP(ISEG,%VAL(ISMAP_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL CAMATE(ISEG,%VAL(ISMATE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL CANODE(ISEG,%VAL(ISNONO_PTR),NPNODE,
     '          NPLIST1,NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OBJECT',1)) THEN
              CALL CAOBJE(ISEG,%VAL(ISOBJE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PLOT',2)) THEN
              CALL CAPLOT(ISEG,%VAL(ISPLOT_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'POLYLINE',2)) THEN
              CALL CAPLIN(ISEG,%VAL(ISPLIN_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PROFILE',2)) THEN
              CALL CAPROF(ISEG,%VAL(ISPROF_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REACTIONS',3)) THEN
              CALL CAREAC(ISEG,%VAL(ISREAC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',3)) THEN
              CALL CARESI(ISEG,%VAL(ISRESI_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RULE',2)) THEN
              CALL CARULE(ISEG,%VAL(ISRULE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SECTION',2)) THEN
              CALL CASECT(ISEG,%VAL(ISSECT_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SHEETS',2)) THEN
              CALL CASHEE(ISEG,%VAL(ISSHEE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL CASTRA(ISEG,%VAL(ISSTRA_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STRESS',4)) THEN
              CALL CASTRE(ISEG,%VAL(ISSTRE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'SURFACE',2)) THEN
              CALL CASURF(ISEG,%VAL(ISSURF_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3
        ELSE IF(ABBREV(CO(2),'CHANGE',3)) THEN
          noco=3
          ADD=.FALSE.
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Bases'
c GBS 2/11/99
c              OPTION_3( 2)='Clock'
c              OPTION_3( 3)='Colours'
              OPTION_3( 2)='Coordinates'
              OPTION_3( 3)='Data'
C alin 22/7/05	      
              OPTION_3(4)='Source'

C OR 22-06-05
C   Adding new option to change displacements 
C   Moved other one one up to keep alphabetical order
C
              OPTION_3( 5)='Displacement_BC' 
              OPTION_3( 6)='Focus'
              OPTION_3( 7)='Lines'
c              OPTION_3( 8)='Lut'

              OPTION_3( 8)='Material'
              OPTION_3( 9)='Mesh'
              OPTION_3( 10)='Nodes'
              OPTION_3(11)='Polylines'
              OPTION_3(12)='Reg-Parameters'
              OPTION_3(13)='Exit'
              NTCH3=13

              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'BASES',1)) THEN
              CALL CHBASE(IBT,IDO,INP,NAN,NBJ,NBJF,NDET,NEELEM,NEL,NENP,
     '          NFF,NFFACE,NGAP,NKB,NKJE,NKEF,NLF,NLL,NLLINE,
     '          NNB,NNF,NNL,NPF,NPL,NPNE,NPNF,NPNODE,NRE,NSB,NVJE,
     '          NVJF,NVJL,DET,DF,DL,PG,RG,SE,SF,WG,XA,XE,XG,XIG,XP,
     '          STRING,ERROR,*9999)
c GBS 2/11/99
c            ELSE IF(ABBREV(CO(3),'CLOCK',2)) THEN
c             CALL CHCLOC(ISCLOC,ISEG,CSEG,STRING,ERROR,*9999)
c GBS 2/11/99
c            ELSE IF(ABBREV(CO(3),'COLOURS',3)) THEN
c             CALL CHCOLO(ISEG,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COORDINATES',2)) THEN
              CALL CHCOOR(NBJ,NENP,NKJ,NPNE,NPNODE,NUNK,
     '          NVJE,NVJP,XG,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',1)) THEN
              CALL CHDATA(IBT,IDO,INP,LD,LN,NBJ,NDDL,NDLT,
     '          NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,NVJE,NXI,SE,
     '          WD,WDL,XA,XE,XID,XP,ZD,STRING,ERROR,*9999)
C ALIN 22/7/05
            ELSE IF(ABBREV(CO(3),'SOURCE',2)) THEN
              CALL CHSOUR(DIPOLE_DIR_NTIME,DIPOLE_DIR, DIPOLE_CEN,
     '              NDIPOLES,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FOCUS',2)) THEN
              CALL CHFOCU(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',2)) THEN
              CALL CHLINE(IBT,IDO,INP,ISEG,%VAL(ISELNO_PTR),
     &          %VAL(ISFIBR_PTR),%VAL(ISFIEL_PTR),%VAL(ISLINE_PTR),
     &          %VAL(ISLINO_PTR),%VAL(ISL2BE_PTR),%VAL(ISL3BE_PTR),
     &          %VAL(ISNONO_PTR),%VAL(ISN2BE_PTR),%VAL(ISN3BE_PTR),MXI,
     &          NAN,NBH,NBJ,NEELEM,NEL,NGAP,NHE,NHP,NKH,NKHE,NKJ,NKJE,
     &          %VAL(NLATNE_PTR),NLL,NLLIST,NPF,NPL,NPNE,NPNODE,
     &          %VAL(NQNE_PTR),%VAL(NQNLAT_PTR),%VAL(NQS_PTR),NQXI,NRE,
     &          NVHE,NVHP,NVJE,NVJL,NW,NYNE,NYNP,CURVCORRECT,DL,SE,XA,
     &          XB,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZE,ZP,CSEG,STRING,FIX,ERROR,
     &          *9999)
c GBS 2/11/99
c            ELSE IF(ABBREV(CO(3),'LUT',2)) THEN
c             CALL CHLUT(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIAL',1)) THEN
              CALL CHMATE(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,GRNGLIST,IBT,ICQS,
     '          ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,ILPIN,
     '          ILTIN,INP,NBJ,NEELEM,NELIST,NENQ,NGLIST,NMBIN,NPLIST1,
     '          NPNE,NPNODE,NQLIST,%VAL(NQNE_PTR),%VAL(NQS_PTR),NQXI,
     '          NRLIST,NW,NXI,NXLIST,CE,CELL_CP,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CGE,CIN,CP,CQ,RCQS,RCQS_SPATIAL,XE,XIG,
     '          YG,YP,YQS,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MESH',2)) THEN
              CALL CHMESH(IBT,IDO,INP,NBJF,NBJ,NEELEM,NEL,NELIST,NENP,
     '          NFF,NFFACE,NGAP,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNF,NNL,
     '          NP_INTERFACE,NPF,NPL,NPLIST1,NPNE,NPNF,NPNODE,NRE,
     '          NRLIST,NUNK,NVJE,NVJF,NVJL,NVJP,DF,DL,
     '          REALWORK(OS_PAOPTI),PG,RG,SE,SF,WG,XA,XE,XG,XIG,XP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL CHNODS(IBT,IDO,INP,ISEG,%VAL(ISELNO_PTR),
     &          %VAL(ISFIBR_PTR),%VAL(ISFIEL_PTR),%VAL(ISLINE_PTR),
     &          %VAL(ISLINO_PTR),%VAL(ISNONO_PTR),MXI,NAN,NBH,NBJ,
     &          NEELEM,NEL,NELIST,NENP,NGAP,NHE,NHP,NKH,NKHE,NKJ,NKJE,
     &          %VAL(NLATNE_PTR),NLL,NLLINE,NLLIST,NNL,NPF,NPL,NPLIST1,
     &          NPNE,NPNODE,%VAL(NQNE_PTR),%VAL(NQNLAT_PTR),
     &          %VAL(NQS_PTR),NQXI,NRE,NRLIST,NVHE,NVHP,NVJE,NVJL,NVJP,
     &          NW,NWP,NXLIST,NYNE,NYNP,CURVCORRECT,DL,
     &          REALWORK(OS_PAOPTI),SE,XA,XE,XG,XP,XQ,YG,YP,YQ,ZA,ZC,ZD,
     &          ZE,ZP,CSEG,STRING,FIX,ERROR,*9999)
C OR 22-06-05
C     Included new command to change displacements according to rotations
C     and translations from its current positioning.
C 
            ELSE IF(ABBREV(CO(3),'DISPLACEMENTS',1)) THEN
              CALL CHDISP(NHP,NKJ,NPLIST1,NPNODE,NRLIST,NVJP,NYNP,
     &             REALWORK(OS_PAOPTI),XP,YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINES',1)) THEN
              CALL CHPLIN(ISEG,%VAL(ISPLIN_PTR),CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REG-PARAMETERS',1)) THEN
              CALL CHREGP(REG_PARAMETER,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'CHECK',3)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Solution'
              OPTION_3( 2)='Convergence'
              OPTION_3( 2)='Voronoi'
              NTCH3=2
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'SOLUTION',1)) THEN
              CALL CHKSOL(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,NBH,
     '          NBHF,NBJ,NBJF,NDIPOLES,NEELEM,NFF,NHE,
     '          NHP,NKEF,NKH,NKHE,NKJE,NNF,NPF,NP_INTERFACE,NPNE,NPNF,
     '          NPNODE,NRLIST,NVHE,NVHF,NVJE,NVJF,NVHP,NW,NXLIST,
     '          NYNE,NYNP,CE,CURVCORRECT,DIPOLE_CEN,DIPOLE_DIR,
     '          PG,RG,SE,SF,WG,XA,XE,XG,XP,
     '          YP,ZA,ZE,ZP,STRING,ERROR,*9999)
            ELSEIF(ABBREV(CO(3),'CONVERGENCE',1)) THEN
              CALL CHECK_CONV(NENQ,NPNODE,NP_INTERFACE,NQNP,
     &          %VAL(NQS_PTR),NQXI,NXQ,NYNP,AQ,CQ,DNUDXQ,DXDXIQ,DXDXIQ2,
     &          YP,YQ,STRING,ERROR,*9999)
            ELSEIF(ABBREV(CO(3),'VORONOI',1)) THEN
              CALL CHECKMSH(NBJ,NEELEM,NFVC,NPLIST1,NPNE,NPNODE,1,
     '          NVCNODE,VC,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

C MPN 26Mar2002: alphabetical fix
        ELSE IF(ABBREV(CO(2),'CLOSE',2)) THEN
          noco=2
          CALL CLOSE_FILES(NHQ,NPNY,NQNY,NRLIST,NRLIST2,NXLIST,
     '      NYNR,NYQNR,YP,YQ,YQS,STRING,ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'COMBINE',4)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Signal'
              NTCH3=1
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'SIGNAL',3)) THEN
              CALL COMBSIG(LD,NBJ,NDDATA,NRLIST,WD,XID,ZD,STRING,
     '          ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'COMPARE',4)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Data'
              OPTION_3( 2)='Gaussvars'
              OPTION_3( 3)='Gridvars'
              OPTION_3( 4)='Signal'
              NTCH3=4
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'DATA',3)) THEN
C*** 10/10/08 JHC Added Z_CONT in order to pass it into IODATA.f
              CALL COMPDAT(NDDATA,NDP,NRLIST,WD,Z_CONT,ZD,ZD2,STRING,
     '          ERROR,*9999)
            ELSEIF(ABBREV(CO(3),'GAUSSVARS',3)) THEN
              CALL COMPGAUSSVAR(IBT,IDO,INP,NEELEM,NELIST,NGLIST,NQET,
     '          NQLIST,%VAL(NQNE_PTR),%VAL(NQS_PTR),NQSCNB,NQXI,NRLIST,
     '          PGNQE,XIG,YG,YQS,ERROR,*9999)
            ELSEIF(ABBREV(CO(3),'GRIDVARS',3)) THEN
              CALL COMPGRIDVAR(NQLIST,NRLIST,YQS,STRING,ERROR,*9999)
            ELSEIF(ABBREV(CO(3),'SIGNAL',3)) THEN
              CALL COMPSIG(LD,LIST,NBJ,NDDATA,NRLIST,WD,XID,ZD,ZD2,
     '          STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

C LKC 23-NOV-97 added convert group
        ELSE IF(ABBREV(CO(2),'CONVERT',3)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Signal'
              NTCH3=1
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'SIGNAL',1)) THEN
              CALL CONVSIG(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'COPY',3)) THEN
          noco=3
          ADD=.FALSE.
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Elements'
              OPTION_3( 2)='Nodes'
              OPTION_3( 3)='Exit'
              NTCH3=3
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ELEMENTS',1)) THEN
              CALL COELEM(IBT,IDO,INP,NBH,NBJ,NEELEM,NEL,NELIST,
     '          NENP,NHE,NKHE,NKJ,NKJE,NLL,NLLINE,
     '          NNL,NPL,NPLIST1,NPNE,NPNODE,NVJE,NVJL,NW,
     '          DL,SE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL CONODE(NHP,NKH,NKJ,NPLIST1,NPNODE,NVJP,
     '          XP,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

!Temporary AJP 6/6/95
        ELSE IF(ABBREV(CO(2),'CRASH',5)) THEN
C#### Command: FEM crash
C###  Description:
C###    Causes crash by assigning 3 parameters in call to two
C###    parameters in subroutine. Used for deliberately crashing
C###    from CMISS.
          IF(.NOT.DOCUMENT) THEN
            CALL CRASH(ERROR,*9999)
          ENDIF

        ELSE IF(ABBREV(CO(2),'DEFINE',2)) THEN
          noco=3
          IF(ABBREV(COQU(2,1),'ADD',1)) THEN
            ADD=.TRUE.
          ELSE
            ADD=.FALSE.
          ENDIF
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Active'
              OPTION_3( 2)='Analytic'
              OPTION_3( 3)='Bases'
C LKC 6-OCT-1999
C              OPTION_3( 4)='Bib'
              OPTION_3( 4)='Boundary'
              OPTION_3( 5)='Cell'
              OPTION_3( 6)='Clock'
              OPTION_3( 7)='Constant'
C*** 22/04/08 JHC Added option 'contact' to define contact parameter values
              OPTION_3( 8)='Contact'
              OPTION_3( 9)='Coordinates'
              OPTION_3(10)='Corners'
              OPTION_3(11)='Coupling'
              OPTION_3(12)='Customisation'
              OPTION_3(13)='Data'
              OPTION_3(14)='Elements'
              OPTION_3(15)='Equation'
              OPTION_3(16)='Export'
              OPTION_3(17)='Faces'
              OPTION_3(18)='Fibres'
C LKC 6-OCT-1999
C              OPTION_3(17)='Fiducial'
              OPTION_3(19)='Field'
              OPTION_3(20)='File..'
              OPTION_3(21)='Fit'
              OPTION_3(22)='Gauss'
              OPTION_3(23)='Grid'
              OPTION_3(24)='Growth'
              OPTION_3(25)='Heading'
              OPTION_3(26)='Import'
              OPTION_3(27)='Increments'
              OPTION_3(28)='Initial'
              OPTION_3(29)='Inverse'
              OPTION_3(30)='Label'
              OPTION_3(31)='Leads'
              OPTION_3(32)='Lines'
              OPTION_3(33)='Map'
              OPTION_3(34)='Mapping'
              OPTION_3(35)='Materials'
              OPTION_3(36)='Maths'
              OPTION_3(37)='Mesh'
              OPTION_3(38)='Motion'
              OPTION_3(39)='Nodes'
              OPTION_3(40)='Noise'
              OPTION_3(41)='Normals'
              OPTION_3(42)='Optimise'
              OPTION_3(43)='Parameters'
              OPTION_3(44)='Regions'
              OPTION_3(45)='Reference'
              OPTION_3(46)='Refinement'
              OPTION_3(47)='Sail'
              OPTION_3(48)='Singularity'
              OPTION_3(49)='Solve'
              OPTION_3(50)='Source'
              OPTION_3(51)='Strain'
              OPTION_3(52)='Surface'
              OPTION_3(53)='Time_variable'
              OPTION_3(54)='Transfer'
              OPTION_3(55)='Window'
              OPTION_3(56)='Xi'
              OPTION_3(57)='Exit'
              NTCH3=57
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ACTIVE',2)) THEN
              CALL DEACTI(NBH,NBJ,NEELEM,NELIST,NGLIST,NRLIST,NXLIST,
     '          FEXT,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ANALYTIC',2)) THEN
CC AJPs 23/10/97
              CALL DEANAL(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     '          INP,NBJ,NDIPOLES,NEELEM,NENP,NKJE,NKH,
     '          NPF,NP_INTERFACE,NPNE,NPNODE,NRE,NRLIST,NVHP,NVJE,
     '          NXLIST,NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,
     '          SE,XA,XE,XP,YP,YQ,STRING,ERROR,*9999)
CC AJPe

            ELSE IF(ABBREV(CO(3),'BASES',2)) THEN
              CALL DEBASE(IBT,IDO,INP,NAN,NDET,NGAP,NKB,NNB,NSB,
     '          DET,PG,WG,XIG,STRING,ERROR,*9999)

C LKC 14-SEP-1999 Routine archived
C            ELSE IF(ABBREV(CO(3),'BIB',3)) THEN
C              CALL DEBIB(IBT,IDO,INP,NBJ,NENP,NKE,NPF,
C     '          NP_INTERFACE,NPNE,NPNODE,NRE,NVJE,NVJP,SE,XA,XE,
C     '          XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BOUNDARY',2)) THEN
              CALL DEBOUN(NBJ,NEELEM,NENP,NNB,NPNE,NRLIST,NXI,NYNP,
     '          YP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CELL',2)) THEN
              CALL DECELL(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,IBT,IDO,INP,NEELEM,
     '          NELIST,NENQ,NPNE,NPNODE,NQET,%VAL(NQNE_PTR),
     '          %VAL(NQS_PTR),NQXI,NRLIST,NXLIST,CE,CELL_RCQS_VALUE,
     '          CELL_YQS_VALUE,CELL_ICQS_NAMES,CELL_RCQS_NAMES,
     '          CELL_YQS_NAMES,CP,CQ,XE,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CHORDS',2)) THEN
              CALL DECHOR(NKJ,NLCHOR,NPL,DL,XP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CLOCK',2)) THEN
              CALL DECLOC(ISEG,%VAL(ISCLOC_PTR),CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONSTANT',5)) THEN
              CALL DEFINE_CONSTANT(STRING,ERROR,*9999)
C*** 22/04/08 JHC Added option contact to define contact parameter values
            ELSE IF(ABBREV(CO(3),'CONTACT',4)) THEN
              CALL DECONT(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COORDINATES',3)) THEN
              CALL DECOOR(NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CORNERS',3)) THEN
              CALL DECORN(INP,NBJ,NKJE,NPF,NPNE,NRLIST,NVHP,NVJE,NW,
     '          PG,SE,XA,XE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COUPLING',3)) THEN
              CALL DECOUP(NEELEM,NP_INTERFACE,NXLIST,NXQ,STRING,
     '          XQ,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CUSTOMISATION',3)) THEN
              CALL DECUST(NXLIST,STRING,ERROR,*9999)
C  SMAR009 18/01/99 removed NRLIST, from list
            ELSE IF(ABBREV(CO(3),'DATA',1)) THEN
C*** 10/10/08 JHC added Z_CONT to DEDATA.f
              CALL DEDATA(IBT,IDO,INP,%VAL(ISDATA_PTR),ISEG,LD,LDR,LN,
     &          MXI,NAN,NBJ,NBH,NDADJ,NDDATA,NDDL,NDLT,NDP,NEELEM,
     &          NELIST,NENP,NHE,NHP,NKH,NKHE,NKJ,NKJE,NLL,
     &          NLS_NDATA_CONT,NLS_NDATA_IMAG,NP_INTERFACE,NPF,NPLIST1,
     &          NPNE,NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNP,
     &          CE,CURVCORRECT,DET,DL,DRDN,
     &          PG,RAD,RD,RG,SE,WD,WG,XA,XE,XG,XG1,XID,XIG,XN,XP,
     &          XR,XQ,YD,YP,ZA,Z_CONT,ZD,ZE,ZF,ZP,CSEG,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL DEELEM(IBT,IDO,INP,ISEG,%VAL(ISELNO_PTR),
     &          %VAL(ISLINE_PTR),%VAL(ISLINO_PTR),%VAL(ISNONO_PTR),MXI,
     &          NBJ,NBJF,NEELEM,NEL,NELIST,NENP,NFF,NFFACE,NKJE,NKEF,
     &          NKJ,NLF,NLL,NLLINE,NLLIST,NNB,NNF,NNL,NPF,NPL,NPLIST1,
     &          NPLIST2,NPNE,NPNF,NPNODE,NRE,NRLIST,NUNK,NVJE,NVJF,NVJL,
     &          NVJP,NXI,DF,DL,PG,RG,
     '          SE,SF,WG,XA,XE,XG,XP,ZC,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EQUATIONS',2)) THEN
              CALL DEEQUA(IBT,IDO,INP,LGE,NBH,NBHF,NBJ,NEELEM,NELIST,
     '          NENP,NFF,NHE,NHP,NHQ,NKH,NKHE,NKJE,NLL,NNF,
     '          NPF,NP_INTERFACE,NPL,NPNE,NPNODE,NPNY,NQNY,NRLIST,NVHE,
     '          NVHP,NVJE,NVJP,NW,NXI,NXLIST,NYNE,NYNP,NYNQ,NYNR,NYQNR,
     '          CE,CURVCORRECT,SE,XA,XAB,XE,XP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EXPORT',2)) THEN
              CALL DEEXPO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FACES',2)) THEN
              CALL DEFACE(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '          NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,
     '          NVJE,NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRES',3)) THEN
              CALL DEFIBR(NKJ,NPNODE,NRLIST,NVJP,XP,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL DEFIEL(NKJ,NEELEM,NENP,NPNODE,NRLIST,NVJP,XAB,XP,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FILE',3)) THEN
              CALL DEFILE(NBH,NBJ,NEELEM,NENP,NHE,NHP,
     '          NKH,NNB,NPL,NPNE,NPNODE,NVHP,NW,NXI,NYNE,
     '          NYNP,CE,DL,XP,YP,ZA,ZP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIT',3)) THEN
              CALL DEFIT(IBT,IDO,INP,LD,LGE,LN,NBH,
     '          NBJ,NBJF,NEELEM,NELIST,NFLIST,NFFACE,
     '          NHE,NHP,NKEF,NKH,NKHE,NKJE,NMNO,NNF,
     '          NPF,NP_INTERFACE,NPLIST1,NPNE,NPNF,NPNODE,NPNY,
     '          NRLIST,NVHE,NVJE,NVJF,NVHP,
     '          NVJP,NXLIST,NYNE,NYNP,NYNR,SE,SF,WU,XA,XE,XID,XP,ZD,ZP,
     '          STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL DEGAUS(NBJ,NRLIST,YG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
C PM 26-JUL-01 : more arguments added
              CALL DEGRID(NAQ,NBJ,NEELEM,NELIST,NENP,NENQ,NLL,NLQ,NNB,
     &          NPNE,NQET,NQLIST,NQSCNB,NQXI,NRLIST,NWQ,NXI,NXQ,
     &          NLATNE_PTR,NLATNQ_PTR,NLATPNQ_PTR,NQNE_PTR,NQNLAT_PTR,
     &          NQS_PTR,STRING,DL,XIQ,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GROWTH',3)) THEN
              CALL DEGROW(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HEADING',1)) THEN
              CALL DEHEAD(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'IMPORT',3)) THEN
              CALL DEIMPO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',3)) THEN
              CALL DEINCR(NHP,NKH,YP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INITIAL',3)) THEN
              CALL DEINIT(IBT,IDO,INP,ITHRES,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NEL,NELIST,NENP,NENQ,NFF,NGAP,NHE,NHP,NHQ,
     '          NKEF,NKH,NKHE,NKJE,NLL,NNF,NNL,NODENVCB,NPL,NPLIST1,
     '          NP_INTERFACE,NPF,NPNE,NPNF,NPNODE,NPNY,NQET,NQLIST,
     '          %VAL(NQNE_PTR),%VAL(NQS_PTR),NRLIST,NTIME_POINTS,
     &          NTIME_NR,NVCB,NVCNODE,NVHE,NVHF,NVHP,NVJE,NVJF,NW,NWQ,
     &          NXI,NXLIST,NYNE,NYNP,NYNQ,NYNR,TV_BC_SET,AQ,BBM,CE,CG,
     &          CGE,CP,CQ,CURVCORRECT,DF,DL,PG,RCQS,RG,SE,SF,THRES,WG,
     &          XA,XE,XG,XIG,XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZG,ZP,FIX,FIXQ,
     '          STRING,TIME_VARIABLE_NAMES,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INVERSE',3)) THEN
              CALL DEINVE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LABEL',5)) THEN
              CALL DEFINE_LABEL(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LEADS',2)) THEN
              CALL DELEAD(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',2)) THEN
              CALL DELINE(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
     '          NKJE,NLL,NLLINE,NNL,NPINTER,NPL,NPNE,NPNODE,
     '          NRE,NRLIST,NVJE,NVJL,DL,SE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAP',3)) THEN
              CALL DEMAP(ISEG,%VAL(ISMAP_PTR),MXI,NBJ,NEELEM,NELIST,
     '          NENP,NNB,NPNE,NRLIST,NXI,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAPPING',4)) THEN
              CALL DEMAPPING(NBJ,NEELEM,NKJ,NPNODE,NPNY,
     '          NRLIST,NVJP,NXLIST,NYNP,NYNY,CYNY,XP,
     '          ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL DEMATE(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,GRNGLIST,IBT,
     '          ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,ILPIN,
     '          ILTIN,INP,NBJ,NEELEM,NELIST,NENP,NENQ,NGLIST,NMBIN,NNB,
     '          NPLIST1,NPNE,NPNODE,NQET,NQLIST,%VAL(NQNE_PTR),
     '          %VAL(NQS_PTR),NQXI,NRLIST,NW,NXI,NXLIST,
     '          CE,CELL_CP,CELL_RCQS_VALUE,CELL_YQS_VALUE,CGE,
     '          CIN,CP,CQ,RCQS_SPATIAL,XE,XIG,YG,YP,YQS,FIX,
     '          CELL_ICQS_NAMES,CELL_RCQS_NAMES,CELL_YQS_NAMES,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATHS',4)) THEN
              CALL DEFINE_MATHS(CELL_RCQS_NAMES,CELL_RCQS_VALUE,
     &          CELL_YQS_NAMES,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MESH',2)) THEN
              CALL DEMESH(IBT,IDO,INP,LD,LD_NP,NBJ,NBJF,NEELEM,
     '          NEL,NELIST,NELIST2,NENFVC,NENP,NEP,NFF,
     &          NFFACE,NFVC,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NLLIST,NNB,NNF,
     &          NNL,NODENVC,NODENVCB,NORD,NPF,NPL,NPLIST1,NPLIST2,NPNE,
     &          NPNF,NPNODE,NPQ,NQLIST,NRE,NRLIST,NUNK,NVCNODE,NVJE,
     &          NVJF,NVJL,NVJP,NXI,NXQ,NP_INTERFACE,BBM,CE,DF,DL,PG,RG,
     &          SE,SF,VC,VC_INIT,WD,WG,XA,XAB,XE,XG,XID,XIP,XNFV,XP,XQ,
     &          ZA,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MOTION',2)) THEN
              CALL DEMOTI(IBT,NBH,NBJ,NEELEM,NELIST,NELIST2,NENP,NHE,
     &          NHP,NKH,NKJ,NKJE,NORD,NPNE,NPNODE,NRLIST,NVJE,NVJP,NXI,
     &          NXLIST,BBM,CE,XAB,XP,YP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',2)) THEN
              CALL DENODS(IBT,IDO,INP,ISEG,%VAL(ISNONO_PTR),NBH,NBJ,
     &          NBJF,NEELEM,NEL,NELIST,NENP,NFF,NFFACE,NFLIST,NHE,NHP,
     &          NKJE,NKEF,NKH,NKJ,NLF,NLL,NLLINE,NNF,NNL,NPF,NPINTER,
     &          NP_INTERFACE,NPL,NPLIST1,NPNE,NPNF,NPNODE,NPNY,NRE,
     &          NRLIST,NVCNP,NVHE,NVHP,NVJE,NVJF,NVJL,NVJP,NW,NXI,
     &          NXLIST,NYNE,NYNP,NYNR,DF,DL,PG,RG,SE,SF,WG,XA,XE,XG,XP,
     &          YP,ZD,ZD2,ZP,FIX,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NOISE',3)) THEN
              CALL DENOIS(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NORMALS',3)) THEN
              CALL DENORM(NELIST,NW,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OPTIMISE',2)) THEN
              CALL DEOPTI(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '          IDO,INP,ISIZE_MFI,ISIZE_PHI,LD,LDR,NBH,NBJ,
     '          NDIPOLES,NEELEM,NENP,NHP,NKHE,NKB,NKH,NKJ,
     '          NLNO,NMNO,NNB,NONL,NONM,NONY,NP1OPT,NP2OPT,NP3OPT,
     '          NPL,NPLIST4,NPNE,NPNODE,NPNY,NRLIST,NVHE,NVHP,
     '          NVJP,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,NYNY,PAOPTY,
     '          CONY,CYNO,CYNY,DIPOLE_CEN,DIPOLE_DIR,
     '          DL,REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PBOPTI),REALWORK(OS_PMAX),REALWORK(OS_PMIN),
     '          XP,YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PARAMETERS',2)) THEN
              CALL DEPARA(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '          ISIZE_MFI,ISIZE_TBH,NDIPOLES,NDLT,NEELEM,NEL,NENP,
     '          NFFACE,NLLINE,NONY,NPNODE,NQET,NQXI,NVHP,NVJP,NYNO,
     '          STRING,YP,ERROR,*9999)
C KAT 2001-04-09 Archived
C            ELSE IF(ABBREV(CO(3),'PLANES',2)) THEN
C              CALL DEPLANE(STRING,ERROR,*9999)
C LKC 4-NOV-98 Archived
C            ELSE IF(ABBREV(CO(3),'POTENTIAL',2)) THEN
C              CALL DEPOTE(NDP,WD,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REFERENCE',4)) THEN
              CALL DEREFE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REFINEMENT',4)) THEN
              CALL DEREFI(IBT,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,
     '          NVJE,NEERR,PG,RG,SE,WG,XA,XE,XG,XP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REGIONS',3)) THEN
              CALL DEREGI(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SAIL',2)) THEN
c             CALL DESAIL(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
c    '          NKE,NKJ,NLL,NLLINE,NNL,NPL,NPLIST1,
c    '          NPNE,NPNODE,NVJE,NVJP,NXI,
c    '          DL,SE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SHEETS',2)) THEN
C             Use define fibre now cs 20-March-1997
C              CALL DESHEE(IDO,NBJ,NEELEM,NKJ,NPNODE,
C     '          NVJE,NVJP,XA,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SINGULARITY',2)) THEN
              CALL DESING(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOLVE',3)) THEN
              CALL DESOLV(IBT,IDO,INP,NAN,NAQ,NBH,NBJ,NEELEM,
     '          NELIST,NENP,NHE,NKB,NKHE,NKJE,NLL,NNB,NNF,NNL,
     '          NONY,NP_INTERFACE,NPF,NPL,NPNE,NPNODE,NPNY,NRE,
     '          NRLIST,NVHE,NVHP,NVJE,NW,NWP,NWQ,NXI,NXLIST,NXQ,
     '          NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,AQ,CONY,CYNO,
     '          CYNY,SE,SP,XA,XE,XP,YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOURCE',3)) THEN
              CALL DESOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     '          INP,NAN,NBH,NDIPOLES,NEELEM,NENQ,NRLIST,NQET,NQLIST,
     '          %VAL(NQNE_PTR),%VAL(NQS_PTR),
     '          NQSCNB,NQXI,NWQ,NXLIST,NXQ,AQ,CG,CQ,DIPOLE_CEN,
     '          DIPOLE_DIR,DXDXIQ,GD,PG,PROPQ,WG,XE,XG,XQ,YQ,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SPAMM',2)) THEN
!             CALL DESPAM(STRING,WD,XP,ZD,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL DESTRA(IBT,IDO,INP,LD,NAN,NBJ,NKJE,
     '          NPF,NPNE,NRE,NVJE,SE,XA,XE,XID,XP,ZP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SURFACE',2)) THEN
              CALL DESURF(IBT,IDO,INP,NBJ,NELIST,NEELEM,NLS_SURF_PSI,
     '          NLS_SURF_XI,NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TIME_VARIABLE',4)) THEN
C              CALL DETIME(IBT,IDO,INP,NAN,NGAP,XE,STRING,ERROR,*9999)
              CALL DETIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,
     '          STRING,TIME_VARIABLE_NAMES,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TRANSFER',3)) THEN
              CALL DETRSF(NELIST3,NHP,NKH,NPLIST3,NPLIST4,NPLIST5,NPNY,
     '          NVHP,NXLIST,NYNP,NYNR,YP,FIX,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VALUES',3)) THEN
              CALL DEVALU(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'WINDOWS',1)) THEN
              CALL DEWIND(NPNODE,XP,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'XI',1)) THEN
              CALL DEXI(FD,IBT,IDO,INP,LD,LD_NP,LN,NBH,NBHF,NBJ,NBJF,
     '          NDDL,NDLT,NDP,NEELEM,NELIST,NENP,NEP,NFF,
     '          NFFACE,NFLIST,NHE,NKEF,NKHE,NKJE,NNB,
     '          NNF,NPF,NPNE,NPNF,NPNODE,NRE,NRLIST,NRLIST2,NVHE,
     '          NVJE,NVJF,NW,NXI,Z_CONT_LIST,CURVCORRECT,SE,SF,SQ,
     '          XA,XAB,XE,XID,XIP,XP,ZA,ZD,ZDD,ZE,ZP,
     &          STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'DISPLAY',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Fibre'
              OPTION_3( 2)='History'
              OPTION_3( 3)='Leads'
              OPTION_3( 4)='Profile'
              OPTION_3( 5)='Section'
              OPTION_3( 6)='Exit'
              NTCH3=6
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'FIBRE',1)) THEN
              CALL DIFIBR(IBT,IDO,INP,ISEG,ISFIPR,NAN,NBJ,
     '          NEELEM,NELIST,NKJE,NPF,NPNE,NRE,NVJE,
     '          SE,XA,XE,XG,XP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HISTORY',1)) THEN
              CALL DIHIST(IBT,IDO,INP,ISEG,%VAL(ISHIST_PTR),NBH,NBJ,
     &          NEELEM,NHE,NHP,NHQ,NKH,NKHE,NKJE,NPF,NPLIST1,NPNE,
     &          NPNODE,NPNY,NQNY,NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NW,
     &          NXLIST,NYNE,NYNP,NYNR,NYQNR,CURVCORRECT,SE,XA,XE,XID,XP,
     &          YP,YQ,YQS,ZA,ZE,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LEADS',1)) THEN
              CALL DILEAD(ISEG,%VAL(ISLEAD_PTR),LD,LIST,NBJ,NDDATA,
     &          NRLIST,WD,XID,ZD,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PROFILE',1)) THEN
              CALL DIPROF(IBT,IDO,INP,ISEG,%VAL(ISPROF_PTR),LD,NAN,NBH,
     &          NBJ,NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,
     '          NPF,NPNE,NPNODE,NRE,NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '          CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RG,SE,XA,XE,XG,XID,
     '          XP,YG,YP,ZA,ZD,ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SECTION',1)) THEN
              CALL DISECT(IBT,IDO,INP,ISEG,%VAL(ISSECT_PTR),NBH,NBJ,
     '          NEELEM,NHE,NHP,NKJE,NKH,NPF,NPLIST1,
     '          NPNE,NPNODE,NRE,NVHP,NVJE,NYNE,NYNP,
     '          CSEG,SE,STRING,XA,XE,XP,YP,ZA,ZD,ZP,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'DRAW',2)) THEN
          noco=3
          IF(ABBREV(COQU(2,1),'ADD',1)) THEN
            ADD=.TRUE.
          ELSE
            ADD=.FALSE.
          ENDIF
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Alignment'
              OPTION_3( 2)='Axes'
              OPTION_3( 3)='Clock'
              OPTION_3( 4)='Contours'
              OPTION_3( 5)='Cross-section'
              OPTION_3( 6)='Data'
              OPTION_3( 7)='Dipole'
              OPTION_3( 8)='Elements'
              OPTION_3( 9)='Fibres'
              OPTION_3(10)='Field'
              OPTION_3(11)='Gauss'
              OPTION_3(12)='Gradients'
              OPTION_3(13)='Grid'
              OPTION_3(14)='Increments'
              OPTION_3(15)='Isochrones'
              OPTION_3(16)='Lines'
              OPTION_3(17)='L-curve'
              OPTION_3(18)='Materials'
              OPTION_3(19)='Nodes'
              OPTION_3(20)='Object'
              OPTION_3(21)='Plot'
              OPTION_3(22)='Polyline'
              OPTION_3(23)='Reactions'
              OPTION_3(24)='Reg-Parameters'
              OPTION_3(25)='Residuals'
              OPTION_3(26)='Rule'
              OPTION_3(27)='Scale'
              OPTION_3(28)='Sheets'
              OPTION_3(29)='Streamlines'
              OPTION_3(30)='Strain'
              OPTION_3(31)='Stress'
              OPTION_3(32)='Surface'
              OPTION_3(33)='Velocity'
              OPTION_3(34)='Xi'
              OPTION_3(35)='Exit'
              NTCH3=35
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ALIGNMENT',2)) THEN
              CALL DRALIG(%VAL(ISALIG_PTR),ISEG,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'AXES',2)) THEN
              CALL DRAXES(%VAL(ISAXES_PTR),ISEG,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CLOCK',2)) THEN
              CALL DRCLOC(%VAL(ISCLOC_PTR),ISEG,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONTOURS',3)) THEN
              CALL DRCONT(IBT,IDO,INP,%VAL(ISCONO_PTR),%VAL(ISCONT_PTR),
     &          ISEG,MXI,NAN,NBJ,NBH,NEELEM,NELIST,NHE,NHP,NKH,NKHE,
     &          NKJE,NPF,NPNE,NPNODE,NRE,NRLIST,NTCOVA,NVHE,NVHP,NVJE,
     &          NW,NXLIST,NYNE,NYNP,COVA,CURVCORRECT,PG,SE,XA,XE,XG,XP,
     '          YG,YP,ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CROSS-SECTION',2)) THEN
              CALL DRCROS(IBT,IDO,INP,%VAL(ISCROS_PTR),ISEG,NAN,NBJ,
     '          NEELEM,NENP,NHE,NKHE,NNB,NPF,NPNE,NRE,NVJE,
     '          NXI,NXLIST,PG,SE,XA,XE,XG,XP,ZD,ZE,ZG,CSEG,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',1)) THEN
              CALL DRDATA(IBT,IDO,INP,%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,LD,LN,MXI,NBJ,
     &          NBH,NDDL,NDLT,NDP,NEELEM,NELIST,NHP,NKH,NKHE,NKJE,NPF,
     &          NPNE,NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NYNE,NYNP,CE,
     &          CG,CGE,CP,CURVCORRECT,PG,SE,WD,WDL,XA,XE,XG,XID,XIDL,XP,
     '          YP,ZA,ZD,ZDD,ZDL,ZE,ZP,CSEG,DATYPE,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DIPOLE',1)) THEN
              CALL DRDIPO(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,DIPOLE_LIST,
     '          ISEG,%VAL(ISDIPO_PTR),%VAL(ISDIPA_PTR),NDIPOLES,NRLIST,
     &          NXLIST,CSEG,DIPOLE_CEN,DIPOLE_DIR,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL DRELEM(IBT,IDO,INP,ISEG,%VAL(ISELNO_PTR),
     &          %VAL(ISERR_PTR),MXI,NAN,NBJ,NEELEM,NELIST,NKJE,NLL,
     '          NPF,NPL,NPNE,NRE,NRLIST,NVJE,
     '          DL,NEERR,SE,XA,XE,XP,XG,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRES',3)) THEN
              CALL DRFIBR(IBT,IDO,INP,ISEG,%VAL(ISFIBR_PTR),MXI,NAN,NBH,
     &          NBJ,NEELEM,NELIST,NENP,NHE,NHP,NKH,NKHE,NKJE,NNB,NPF,
     &          NPNE,NPNODE,NRLIST,NVHE,NVHP,NVJE,NW,NXI,NXLIST,NYNE,
     &          NYNP,CURVCORRECT,SE,XA,XE,XG,XP,YP,ZA,ZE,ZP,CSEG,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL DRFIEL(IBT,IDO,INP,ISEG,%VAL(ISFIEL_PTR),
     &          %VAL(ISSCAL_PTR),MXI,NBH,NBJ,NEELEM,NELIST,NGAP,NHE,NHP,
     &          NHQ,NKH,NKHE,NKJ,NKJE,%VAL(NLATNE_PTR),NPF,NPNE,NPNODE,
     &          NPNY,%VAL(NQNLAT_PTR),%VAL(NQNE_PTR),NQNY,%VAL(NQS_PTR),
     '          NQXI,NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NW,NXLIST,
     '          NYNE,NYNP,NYNR,NYQNR,CURVCORRECT,PG,SE,XA,XE,XG,XP,XQ,
     '          YG,YQ,YQS,YP,ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL DRGAUS(ISEG,%VAL(ISGAUS_PTR),MXI,NBH,NBJ,
     '          NEELEM,NELIST,NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,
     '          NYNE,NYNP,CURVCORRECT,SE,PG,XA,XE,XG,XIG,XP,
     '          YP,ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRADIENTS',3)) THEN
              CALL DRGRAD(ISEG,%VAL(ISGRAD_PTR),
     '          NBH,NBJ,NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '          CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,XA,XE,XG,XP,YP,ZA,
     '          ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
              CALL DRGRID(ISEG,%VAL(ISGRID_PTR),NAQ,NELIST,NLQ,NQET,
     &          NQLIST,%VAL(NQNE_PTR),%VAL(NQS_PTR),NWQ,NXLIST,NXQ,
     &          DXDXIQ,DXDXIQ2,XQ,YQ,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',3)) THEN
              CALL DRINCR(ISEG,%VAL(ISINCR_PTR),NHP,NKH,NPNODE,NRLIST,
     &          NXLIST,XP,ZP,CSEG,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ISOCHRONES',2)) THEN
              CALL DRISOC(ISEG,%VAL(ISISOC_PTR),ITHRES,NBJ,NEELEM,NGAP,
     &          NKJE,NPF,NPNE,NRE,NVJE,
     '          SE,THRES,XA,XE,XP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'L-CURVE',2)) THEN
              CALL DRLCURVE(ISEG,ISIZE_TBH,%VAL(ISPLOTXY_PTR),PHI,
     &          REG_PARAMETER,SIGMA_T_BH,U_PHI,U_T_BH,WK1_INV,CSEG,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',2)) THEN
              CALL DRLINE(IBT,ISEG,%VAL(ISLINE_PTR),%VAL(ISLINO_PTR),
     &          NBH,NBJ,NEELEM,NELIST,NHE,NHP,NKH,NKJ,NLL,NLLINE,NLLIST,
     '          NONY,NPL,NPNODE,NPNY,NRLIST,NVHP,NW,NXLIST,
     '          NYNE,NYNP,NYNR,CONY,
     '          DL,EIGVEC,XP,YP,ZA,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL DRMATE(IBT,IDO,INP,ISEG,%VAL(ISMATE_PTR),
     '          NAN,NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRE,NVJE,NXLIST,
     '          CE,CQ,SE,XA,XE,XG,XP,XQ,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL DRNODS(ISEG,%VAL(ISNONO_PTR),NBH,NEELEM,NHE,NHP,NKH,
     '          NPNODE,NRLIST,NVHP,NXLIST,NYNE,NYNP,
     '          XP,YP,ZA,ZP,CSEG,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OBJECT',2)) THEN
              CALL DROBJE(ISEG,%VAL(ISOBJE_PTR),LD,MXI,NDDL,NDLT,NEELEM,
     '          XID,ZD,ZDD,WD,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PLOT',2)) THEN
              CALL DRPLOT(ISEG,%VAL(ISPLOT_PTR),NBJ,NBH,NEELEM,NELIST,
     &          NHE,NHP,NKH,NKHE,NKJE,NPF,NPL,NPNE,NPNODE,NRE,
     '          NVHE,NVHP,NVJE,NW,NYNE,NYNP,CURVCORRECT,
     '          DL,SE,XA,XE,XP,YP,ZA,ZE,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINE',5)) THEN
              CALL DRPLIN(ISEG,%VAL(ISPLIN_PTR),REALWORK(OS_PAOPTI),WD,
     &          ZD,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYMARKER',5)) THEN
              CALL DRPMAR(ISEG,%VAL(ISPMAR_PTR),CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REACTIONS',3)) THEN
              CALL DRREAC(ISEG,%VAL(ISREAC_PTR),NBH,NEELEM,NHE,NHP,NKH,
     '          NPL,NPNODE,NRLIST,NVHP,NXLIST,NYNE,NYNP,NYNR,
     '          XP,YP,ZA,ZP,CSEG,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REG-PARAMETERS',3)) THEN
              CALL DRREGP(ISEG,%VAL(ISPLOTXY_PTR),REG_PARAMETER,WK5_INV,
     &          CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',3)) THEN
              CALL DRRESI(ISEG,%VAL(ISRESI_PTR),NBH,NEELEM,NHE,NHP,NKH,
     '          NPNODE,NPNY,NVHP,NXLIST,NYNE,NYNO,NYNP,
     '          REALWORK(OS_RESID),YP,ZA,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RULE',2)) THEN
              CALL DRRULE(ISEG,%VAL(ISRULE_PTR),CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SCALE',2)) THEN
              CALL DRSCAL(ISEG,%VAL(ISSCAL_PTR),CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SHEETS',2)) THEN
              CALL DRSHEE(IBT,IDO,INP,ISEG,%VAL(ISSHEE_PTR),NAN,NBH,NBJ,
     '          NEELEM,NELIST,NENP,NHE,NHP,NKJE,NKH,NNB,NPF,NPNE,
     '          NPNODE,NRE,NVHP,NVJE,NXI,NYNE,NYNP,
     '          SE,XA,XE,XG,XP,YP,ZA,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL DRSTRA(IBT,IDO,INP,ISEG,%VAL(ISSTRA_PTR),LD,NAN,NBH,
     &          NBJ,NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,
     '          NYNE,NYNP,CURVCORRECT,PG,RG,SE,XA,XE,XG,XID,XIG,XP,YP,
     '          ZA,ZE,ZG,ZP,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRESS',5)) THEN
              CALL DRSTRE(IBT,IDO,INP,ISEG,%VAL(ISSTRE_PTR),NAN,NBH,NBJ,
     '          NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,
     '          NYNE,NYNP,CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RG,SE,
     '          XA,XE,XG,XIG,XP,YG,YP,ZA,ZE,ZG,ZP,CSEG,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STREAMLINES',5)) THEN
              CALL DRSTRM(ISEG,%VAL(ISSTRM_PTR),NEELEM,NELIST,
     '          CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SURFACE',2)) THEN
              CALL DRSURF(IBT,IDO,INP,ISEG,%VAL(ISSURF_PTR),NAN,NBH,NBJ,
     '          NEELEM,NELIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NVHE,NVHP,NVJE,NW,NYNE,NYNP,
     '          CURVCORRECT,PG,RG,SE,XA,XE,XG,XP,YP,ZA,ZE,ZG,ZP,
     '          CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VELOCITY',2)) THEN
              CALL DRVELO(NEELEM,NELIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'XI',1)) THEN
              CALL DRXI(STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'DUPLICATE',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3(1)='Mesh'
              NTCH3=1
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'MESH',1)) THEN
              CALL DUMESH(IBT,IDO,INP,NBJ,NBJF,NEELEM,NEL,
     '          NELIST,NENP,NFF,NFFACE,NKJE,NKEF,NKJ,NLF,NLL,
     '          NLLINE,NNB,NNF,NNL,NPF,NP_INTERFACE,NPL,NPNE,NPNF,
     '          NPNODE,NRE,NRLIST,NVJE,NVJF,NVJL,NVJP,NXI,
     '          DF,DL,PG,RG,SE,SF,STRING,WG,XA,XE,XG,XP,
     '          ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'EVALUATE',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3(1)='Aerofoil'
              OPTION_3(2)='Compressible'
              OPTION_3(3)='Constraints'
              OPTION_3(4)='Coronary'
              OPTION_3(5)='Coupling'
              OPTION_3(6)='Electrodes'
              OPTION_3(7)='Error'
              OPTION_3(8)='Event'
              OPTION_3(9)='Fibre'
              OPTION_3(10)='Field'
              OPTION_3(11)='Flow'
              OPTION_3(12)='Grid'
              OPTION_3(13)='Integral'
              OPTION_3(14)='Inverse'
              OPTION_3(15)='Noise'
              OPTION_3(16)='Laplacian'
              OPTION_3(17)='Objective'
              OPTION_3(18)='Ordering'
              OPTION_3(19)='MFI'
              OPTION_3(20)='Phi'
              OPTION_3(21)='Pressure'
              OPTION_3(22)='Pulmonary'
              OPTION_3(23)='Reactions'
              OPTION_3(24)='Residuals'
              OPTION_3(25)='Sensitivity'
              OPTION_3(26)='Signal'
              OPTION_3(27)='Solution'
              OPTION_3(28)='Strain'
              OPTION_3(29)='Time_variable'
              OPTION_3(30)='Transfer'
              OPTION_3(31)='Zeroxing'
              OPTION_3(32)='Exit'
              NTCH3=32
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'AEROFOIL',1)) THEN
              CALL EVAERO(NBH,NEELEM,NHE,NHP,
     '          NKH,NKJ,NPL,NPNODE,NVHP,NYNE,NYNP,
     '          CE,DL,DLL,WG,XIG,XP,YP,ZA,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COMPRESSIBLE',3)) THEN
              CALL EVCOMP(IBT,IDO,INP,LD,NAN,NBH,NBJ,NBJF,NDDATA,
     &          NEELEM,NEP,NFF,NFFACE,NFLIST,NHE,NHP,NKEF,NKH,NKHE,NKJE,
     &          NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,NP1OPT,NPF,
     &          NP_INTERFACE,NPL,NPLIST1,NPNE,NPNF,NPNODE,NPNY,NRE,
     &          NRLIST,NSB,NVHE,NVHP,NVJE,NVJF,NVJL,NVJP,NW,NXI,NXLIST,
     &          NYNE,NYNO,NYNP,NYNR,PAOPTY,AQ,CE,CELL_RCQS_VALUE,CG,CGE,
     &          CONY,CP,CURVCORRECT,DF,FEXT,PG,RG,SE,SF,WG,XA,XE,XG,XID,
     &          XIG,XIP,XN,XP,YG,YGF,YP,ZA,ZA1,ZD,ZE,ZG,ZG1,ZP,ZP1,
     &          STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONSTRAINTS',3)) THEN
              CALL EVCONT(IBT,IDO,NBJ,INTWORK(OS_NEEDCON),
     '          NEELEM,NEL,NKH,NLL,NLNO,NNL,NONL,
     '          NONY,NPL,NPNE,NPNODE,NPNY,NRLIST,NVHP,NVJL,
     '          NVJP,NXLIST,NYNO,NYNP,PAOPTY,
     '          REALWORK(OS_CM),REALWORK(OS_CJACM),
     '          REALWORK(OS_CONTR),REALWORK(OS_CONJAC),DL,
     '          REALWORK(OS_PAOPTI),SE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CORONARY',3)) THEN
              CONECT_PTR=0
              CALCULATED_PTR=0
              BC_POINTS_PTR=0
              BRANCH_PTR=0
              CALL ALLOCATE_MEMORY(3*3*NQM,0,INTTYPE,CONECT_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL ALLOCATE_MEMORY(NQM,0,INTTYPE,BC_POINTS_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL ALLOCATE_MEMORY(NQM,0,LOGTYPE,BRANCH_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL ALLOCATE_MEMORY(NQM,0,LOGTYPE,CALCULATED_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL EVCORO(%VAL(BC_POINTS_PTR),%VAL(BRANCH_PTR),
     '          %VAL(CALCULATED_PTR),%VAL(CONECT_PTR),NBJ,NENQ,
     '          NEP,NPNE,NQET,%VAL(NQNE_PTR),%VAL(NQS_PTR),
     '          NXQ,NYNQ,XIP,XQ,YQ,STRING,ERROR,*9999)
              CALL FREE_MEMORY(CALCULATED_PTR,ERROR,*9999)
              CALL FREE_MEMORY(CONECT_PTR,ERROR,*9999)
              CALL FREE_MEMORY(BRANCH_PTR,ERROR,*9999)
              CALL FREE_MEMORY(BC_POINTS_PTR,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COUPLING',3)) THEN
              CALL EVCOUP(NPNODE,NRLIST,NYNP,STRING,XP,YP,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ELECTRODES',3)) THEN
              CALL EVELEC(IBT,IDO,INP,ISIZE_PHIH,ISIZE_MFI,
     '          LD,LIST,NBH,NBJ,NDDATA,NEELEM,NELIST,NENP,
     '          NENQ,NHE,NHP,NHQ,NKH,NKHE,NKJE,NPF,NP_INTERFACE,NPLIST1,
     '          NPLIST3,NPLIST4,NPNE,NPNODE,NPNY,
     '          NQLIST,%VAL(NQNE_PTR),NQNY,%VAL(NQS_PTR),NQSCNB,NQXI,
     '          NRE,NRLIST,NRLIST2,NVHE,NVJE,NVHP,NW,NXLIST,NYNE,NYNP,
     '          NYNR,NYQNR,CURVCORRECT,MFI,PHI_H,SE,WD,XA,XE,XID,XP,XQ,
     '          YP,YQ,YQS,ZA,ZD,ZE,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ERROR',3)) THEN
              CALL EVERROR(NEELEM,NELIST,NKJE,NPF,NPNE,
     '          NRLIST,NVJE,NW,CG,NEERR,PG,SE,
     '          WG,XA,XE,XG,XP,YG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EVENT',3)) THEN
              CALL EVEVENT(NBJ,NDDATA,ISIZE_PHIH(1),ISIZE_PHIH(2),
     '          LD,PHI_H,WD,XID,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRE',3)) THEN
              CALL EVFIBR(IBT,IDO,INP,LD,NAN,NBJ,NEELEM,NENP,NKJE,NKJ,
     '          NNB,NPF,NPNE,NPNODE,NVJE,NXI,SE,STRING,WD,XA,XE,XG,XG1,
     '          XID,XP,ZD,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL EVFIEL(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NELIST,
     '          NGLIST,NHE,NHP,NKH,NKHE,NKJE,NLL,NPF,NPL,NPNE,NPNODE,
     '          NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,CURVCORRECT,
     '          DL,SE,XA,XE,XG,XIG,XP,YP,ZA,ZE,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FLOW',2)) THEN
              CALL EVFLOW(NBJ,NDP,NEELEM,NELIST,NENP,NORD,NPLIST1,
     &          NPNE,NPNODE,NRLIST,NVJE,NVJP,NXI,
     &          BBM,CE,XAB,XP,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',1)) THEN
              CALL EVGRID(NAQ,NLQ,NWQ,NXQ,
     '          GCHQ,GUQ,PROPQ,YQ,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INTEGRAL',3)) THEN
              CALL EVINTE(NLLINE,NLLIST,NPL,NRLIST,
     &          XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INVERSE',3)) THEN
              IPIV_PTR=0
              CALL ALLOCATE_MEMORY(NY_TRANSFER_M,0,1,IPIV_PTR,.TRUE.,
     &          ERROR,*9999)
              CALL EVINVE(%VAL(IPIV_PTR),ISEG,ISIZE_PHI,ISIZE_PHIH,
     &          ISIZE_TBH,%VAL(ISPLOTXY_PTR),NXLIST,PHI,PHI_H,
     &          PHI_H_EXACT,REG_PARAMETER,SIGMA_PHI,SIGMA_T_BH,T_BH,
     &          T_BH_INV,U_PHI,U_T_BH,VT_PHI,VT_T_BH,WK1_INV,WK2_INV,
     &          WK3_INV,WK4_INV,CSEG,STRING,ERROR,*9999)
              CALL FREE_MEMORY(IPIV_PTR,ERROR,*9999)
C GBS 28-NOV-2001 Adding evlapl
            ELSE IF(ABBREV(CO(3),'LAPLACIAN',3)) THEN
              CALL EVLAPL(IBT,NBH,NENP,NPLIST3,NPNE,NXI,NXLIST,
     '          LAPL,LAPLSQR,XP,STRING,ERROR,*9999)
C LKC 10-FEB-2003 adding evmfi
            ELSE IF(ABBREV(CO(3),'MFI',3)) THEN
              CALL EVMFI(ISIZE_MFI,LD,NBJ,NDDATA,MFI,WD,XID,ZD,
     '          STRING,ERROR,*9999)
C LKC 7-AUG-2000 Adding evnois              
            ELSE IF(ABBREV(CO(3),'NOISE',3)) THEN
              CALL EVNOIS(LD,LIST,NBJ,NDDATA,WD,XID,ZD,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OBJECTIVE',1)) THEN
              CALL EVOBJE(IBT,IDO,INP,LD,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,
     '          NKHE,NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,
     '          NP_INTERFACE,NP1OPT,NPF,NPL,NPNE,NPNODE,NPNY,NRE,
     '          NSB,NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNO,NYNP,NYNR,PAOPTY,
     '          Z_CONT_LIST,
     '          CE,CG,CGE,CP,CURVCORRECT,DL,FEXT,
     '          REALWORK(OS_FGRADM),REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PBOPTI),PG,REALWORK(OS_PMIN),
     '          REALWORK(OS_PMAX),RE1,REALWORK(OS_RESIDM),
     '          REALWORK(OS_RESJAC),RG,SE,WG,WU,XA,XE,XG,
     '          XID,XIG,XP,
     '          YG,YGF,YP,ZA,ZA1,Z_CONT,ZD,ZE,ZG,ZP,ZP1,STRING,
     '          FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ORDERING',3)) THEN
              CALL ASSERT(USE_LUNG.EQ.1,
     &          '>>Must set USE_LUNG to 1 in .ippara file',ERROR,*9999)
              CALL EVORDR(NBJ,NEELEM,NELIST,NENP,NORD,NPNE,NPNODE,
     &          NRLIST,NXI,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PHI',3)) THEN
              CALL EVPHI(ISIZE_PHI,ISIZE_PHIH,LD,NBJ,
     '          NDDATA,NHP,NHQ,NKH,NPLIST3,
     '          NPLIST4,NPLIST5,NPNY,NQNY,NRLIST,NVHP,NXLIST,NYNP,NYNR,
     '          NYQNR,PHI,PHI_H,PHI_H_EXACT,SIGMA_PHI,U_PHI,VT_PHI,
     '          WD,WK1_INV,
     '          WK3_INV,XID,XP,YP,YQ,YQS,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PRESSURE',3)) THEN
              CALL EVPRESS(IBT,IDO,INP,LD,NAN,NBH,NBJ,NBJF,NDDATA,
     &          NEELEM,NEP,NFF,NFFACE,NFLIST,NHE,NHP,NKEF,NKH,NKHE,
     &          NKJE,NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,NP1OPT,
     &          NPF,NP_INTERFACE,NPL,NPLIST1,NPNE,NPNF,NPNODE,NPNY,NRE,
     &          NRLIST,NSB,NVHE,NVHP,NVJE,NVJF,NVJL,NVJP,NW,NXI,NXLIST,
     &          NYNE,NYNO,NYNP,NYNR,PAOPTY,AQ,CE,CELL_RCQS_VALUE,CG,CGE,
     &          CONY,CP,CURVCORRECT,DF,FEXT,PG,RG,SE,SF,WG,XA,XE,XG,XID,
     &          XIG,XIP,XN,XP,YG,YGF,YP,ZA,ZA1,ZD,ZE,ZG,ZG1,ZP,ZP1,
     &          STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PULMONARY',3)) THEN
              CALL EVPULM(NBJ,NEELEM,NELIST,NPNE,NRLIST,NVJE,NXI,NXLIST,
     &          NYNP,CE,CP,XP,YP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REACTIONS',3)) THEN
              CALL EVREAC(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '          NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '          NNF,NONY,NPLIST1,NPF,NPNE,
     '          NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,
     '          NVJE,NW,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,Z_CONT_LIST,CE,
     '          CG,CGE,CONY,CP,CURVCORRECT,FEXT,FIX,GRR,PG,
     '          RE1,RG,SE,
     '          WG,XA,XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZP,ZP1,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',3)) THEN
              CALL EVRESI(IBT,IDO,INP,ISIZE_MFI,
     '          ISIZE_PHI,ISIZE_TBH,LD,LDR,LGE,
     '          NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NELIST,NENP,NFF,
     '          NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,
     '          NLL,NLNO,NMNO,NNB,NNF,NNL,NONL,NONM,NONY,NP1OPT,NPF,
     '          NP_INTERFACE,NPL,NPLIST3,NPNE,NPNODE,NPNY,NRE,
     '          NRLIST,NSB,NVHE,NVHP,NVJE,NVJL,NW,NXI,NXLIST,
     '          NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,
     '          AQ,CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,
     '          CURVCORRECT,DL,
     '          D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,REALWORK(OS_FGRADM),
     '          LAPL,LAPLSQR,
     '          MFI,REALWORK(OS_PAOPTI),REALWORK(OS_PBOPTI),
     '          PG,PHI,PHI_H,REALWORK(OS_PMIN),REALWORK(OS_PMAX),
     '          RE1,RE2,REALWORK(OS_RESID),REALWORK(OS_RESJAC),RG,
     '          SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,XN,
     '          XP,YG,YGF,YP,
     '          ZA,ZA1,Z_CONT,ZD,ZE,ZG,ZG1,ZP,ZP1,STRING,FIX,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SENSITIVITY',3)) THEN
              CALL EVSENS(INTWORK(OS_ISTATE),INTWORK,
     '          NEELEM,NMNO,NPNODE,NVJP,NYNP,
     '          REALWORK(OS_PAOPTI),PE,PF,
     '          REALWORK(OS_PMIN),REALWORK(OS_PMAX),REALWORK(OS_R),
     '          REALWORK(OS_F),REALWORK(OS_FJAC),REALWORK,
     '          REALWORK(OS_XC),XP,YP,STRING,FIX,FIXP,ERROR,*9999)

            ELSE IF(ABBREV(CO(3),'SIGNAL',3)) THEN
              CALL EVSIGN(LD,NBJ,NDDATA,WD,XID,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOLUTION',3)) THEN
              CALL EVSOLU(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,IBT,IDO,
     &          INP,ISIZE_MFI,NBH,NBJ,NDDATA,NDET,NDIPOLES,NEELEM,NENP,
     &          NGAP,NHE,NHP,NKH,NKHE,NKJE,NLL,NPF,NP_INTERFACE,
     &          NPNE,NPNODE,NQET,%VAL(NQNE_PTR),%VAL(NQS_PTR),NRE,
     &          NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,CE,
     &          CURVCORRECT,DET,DIPOLE_CEN,DIPOLE_DIR,DL,DRDN,PG,MFI,
     &          RAD,RD,RG,SE,WD,WG,XA,XE,
     &          XG1,XIG,XN,XP,XQ,XR,YD,YP,YQ,ZA,ZD,ZE,
     &          ZF,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',3)) THEN
              CALL EVSTRA(NBJ,NEELEM,NELIST,NGLIST,%VAL(NQNE_PTR),
     &          %VAL(NQS_PTR),NQSCNB,NQXI,NRLIST,NXLIST,YQS,XIG,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TIME_VARIABLE',3)) THEN
              CALL EVTIME(NTIME_INTERP,NTIME_POINTS,NVARINDEX,EVALTIME,
     '          TIME_VALUES,VALUE,STRING,TIME_VARIABLE_NAMES,.TRUE.,
     '          YP,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TRANSFER',3)) THEN
              IPIV_PTR=0
              CALL ALLOCATE_MEMORY(NY_TRANSFER_M,0,1,IPIV_PTR,.TRUE.,
     '          ERROR,*9999)
              CALL EVTRSF(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '          IBT,IDO,INP,%VAL(IPIV_PTR),ISC_GKK,
     '          ISIZE_TBH,ISR_GKK,LD,LGE,NAN,NBH,NBJ,
     '          NDET,NDIPOLES,NEELEM,NELIST,NENP,NGAP,NHE,NHP,
     '          NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NORD,NPF,
     '          NP_INTERFACE,NPL,NPLIST3,NPLIST4,NPNE,
     '          NPNODE,NPNY,NRE,NVHE,NVHP,NVJE,NW,NWP,NXI,NXLIST,
     '          NYNE,NYNO,
     '          NYNP,NYNR,NYNY,NYQNR,CE,CGE,CONY,CP,CURVCORRECT,CYNO,
     '          CYNY,DET,DIPOLE_CEN,DIPOLE_DIR,DL,
     '          GD,GKK,GR,GRR,PG,SE,SP,T_BH,WG,
     '          WK1_INV,WK2_INV,WK3_INV,WK4_INV,XA,XE,XG,XID,XIG,
     '          XO,XP,YG,YP,ZA,ZE,ZP,STRING,FIX,
     '          ERROR,*9999)
              CALL FREE_MEMORY(IPIV_PTR,ERROR,*9999)

            ELSE IF(ABBREV(CO(3),'ZEROXING',3)) THEN
              CALL EVZEROXING(ISIZE_PHI,ISIZE_TBH,NPLIST4,
     '          NXLIST,PHI,SIGMA_PHI,
     '          T_BH,U_PHI,VT_PHI,WK1_INV,WK2_INV,WK4_INV,ZCROSSING,
     '          STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'EXPORT',2)) THEN
          noco=3
          IF(ABBREV(COQU(2,1),'ADD',1)) THEN
            ADD=.TRUE.
          ELSE
            ADD=.FALSE.
          ENDIF
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Data'
              OPTION_3( 2)='Elements'
              OPTION_3( 3)='FieldML'
              OPTION_3( 4)='GAUSS'
              OPTION_3( 5)='Grid'
              OPTION_3( 6)='Nodes'
              OPTION_3( 7)='Points'
              OPTION_3( 8)='Properties'
              OPTION_3( 9)='Signal'
              OPTION_3(10)='Source'
              OPTION_3(11)='Textures'
              OPTION_3(12)='Voronoi'
              OPTION_3(13)='Exit'
              NTCH3=13
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'DATA',1)) THEN
C*** 19/02/08 JHC Added Z_CONT to exdata for outputing contact force data
C              CALL EXDATA(IBT,IDO,INP,LD,LIST,NBJ,NDDATA,NDP,NENQ,NKJE,
C     &          NPF,NPNE,NRE,NRLIST,NVJE,NWQ,NXQ,AQ,DXDXIQ,DXDXIQ2,SE,
C     &          XA,XE,XID,XIQ,XP,ZD,STRING,ERROR,*9999)
              CALL EXDATA(IBT,IDO,INP,LD,LIST,NBJ,NDDATA,NDP,NENQ,NKJE,
     &          NPF,NPNE,NRE,NRLIST,NVJE,NWQ,NXQ,AQ,DXDXIQ,DXDXIQ2,SE,
     &          XA,XE,XID,XIQ,XP,Z_CONT,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',1)) THEN
              CALL EXELEM(IBT,IDO,INP,NBH,NBJ,NBJF,NEELEM,NELIST,
     '          NFF,NFLIST,NHE,NHP,NKB,NKHE,NKJ,NKJE,%VAL(NLATNE_PTR),
     '          NLF,NLL,NLLIST,NNB,NNF,NPNE,NQET,%VAL(NQNE_PTR),
     '          %VAL(NQNLAT_PTR),%VAL(NQS_PTR),NQXI,NRLIST,NSB,
     '          NVHE,NVJE,NVJP,NW,NXLIST,NYNQ,SE,YQ,YQS,
     '          CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '          CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '          CELL_YQS_VALUE ,CELL_YQS_SPATIAL ,CELL_YQS_NAMES ,
     '          ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '          STRING,ERROR,*9999)
C DPN 07-03-2003 - Adding export to FieldML
            ELSE IF(ABBREV(CO(3),'FIELDML',1)) THEN
              CALL EXFIELDML(LIST,NQLIST,NXLIST,XQ,YQS,
     '          CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '          CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '          CELL_YQS_VALUE,CELL_YQS_SPATIAL ,CELL_YQS_NAMES,
     '          ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL EXGAUS(NBJ,NEELEM,NRLIST,NXLIST,CGE,XIG,YG,
     '          STRING,ERROR,*9999)
C LKC 8-JUN-1999 archived routine
C            ELSE IF(ABBREV(CO(3),'GEOMETRY',2)) THEN
C              CALL EXGEOM(IBT,IDO,INP,NBJ,NEELEM,NELIST,NKE,
C     '          NPF,NPLIST1,NPNE,NRE,NRLIST,NVJE,
C     '          SE,STRING,XA,XE,XP,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',2)) THEN
              CALL EXGRID(NEELEM,%VAL(NLATNE_PTR),NQLIST,
     &          %VAL(NQNLAT_PTR),%VAL(NQS_PTR),NQSCNB,NQXI,
     &          NRLIST,NWQ,NXQ,XQ,YQ,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL EXNODE(IBT,NBH,NEELEM,NHE,NHP,NHQ,
     '          NKH,NKJ,NP_INTERFACE,NPLIST1,NPNODE,NPNY,NQNY,
     '          NRLIST,NRLIST2,NVHP,NVJP,NXLIST,NYNE,NYNP,NYNR,
     '          NYQNR,NW,XP,YP,YQ,YQS,ZA,ZA1,ZP,ZP1,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POINTS',1)) THEN
              CALL EXPOIN(LIST,NBH,NBJ,NHE,NKHE,NKJE,
     '          NPF,NPNE,NPNODE,NRE,NVHE,NVJE,NW,NXLIST,
     '          NYNQ,AQ,CE,CG,CGE,CP,CURVCORRECT,PG,SE,XA,XE,XG,XP,
     '          XQ,YQ,ZA,ZE,ZG,ZP,
     '          YQS,
     '          CELL_ICQS_VALUE,CELL_ICQS_SPATIAL,CELL_ICQS_NAMES,
     '          CELL_RCQS_VALUE,CELL_RCQS_SPATIAL,CELL_RCQS_NAMES,
     '          CELL_YQS_VALUE ,CELL_YQS_SPATIAL ,CELL_YQS_NAMES ,
     '          ICQS_SPATIAL,IICQS_SPATIAL,IRCQS_SPATIAL,RCQS_SPATIAL,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PROPERTIES',1)) THEN
              CALL EXPROP(NBJ,NEELEM,NELIST,NRLIST,NXLIST,CE,STRING,
     '          XAB,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SIGNAL',3)) THEN
C *** DPN 19 February 2000 - adding cell names
              CALL EXSIGN(ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,LD,
     &          MFI,NBJ,NDDATA,NDP,NPLIST3,PHI,PHI_H,WD,XID,ZCROSSING,
     '          ZD,CELL_YQS_NAMES,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOURCE',3)) THEN
              CALL EXSOUR(NDIPOLES,DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '          NRLIST,NXLIST,
     '          DIPOLE_DIR,DIPOLE_CEN,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TEXTURES',1)) THEN
              CALL EXTEXT(NEELEM,NELIST,%VAL(NQNE_PTR),%VAL(NQS_PTR),
     &          NQXI,NRLIST,YQ,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VORONOI',2)) THEN
              VORO_PTR=0
              VORO_SIZE=8*NEM
              CALL ALLOCATE_MEMORY(VORO_SIZE,0,INTTYPE,VORO_PTR,
     '          .TRUE.,ERROR,*9999)
              CALL EXVORO(NBJ,NEELEM,NXI,%VAL(VORO_PTR),VORO_SIZE,ZA,
     '          STRING,ERROR,*9999)
              CALL FREE_MEMORY(VORO_PTR,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'FIT',1)) THEN
          noco=2
          CALL FIT(IBT,IDO,INP,IPIVOT,ISC_GKK,
     '      ISR_GKK,
     '      LD,LGE,LN,NAN,NBH,NBHF,NBJ,NBJF,NDDATA,NDDL,
     '      NDLT,NEELEM,NELIST,NENP,NENQ,
     '      NEP,NFF,NFFACE,NGAP,NHE,
     '      NHP,NHQ,NKB,NKEF,NKH,NKHE,NKJE,NLL,NMNO,NNB,NNF,NNL,
     '      NONY,NPF,NP_INTERFACE,NPL,NPLIST1,NPLIST2,NPNE,NPNF,
     '      NPNODE,NPNY,%VAL(NQNE_PTR),NQNY,%VAL(NQS_PTR),NQXI,NRE,
     '      NRLIST,NRLIST2,NSB,NVHE,NVHP,NVJE,NVJF,NW,NWP,NXI,NXLIST,
     '      NYNE,NYNO,NYNP,NYNR,NYNY,NYQNR,Z_CONT_LIST,CE,CG,CGE,
     '      CONY,CYNY,CP,CURVCORRECT,CYNO,EDD,ER,ES,FEXT,GKK,
     '      GR,GRR,PG,RE1,RG,SE,SF,SP,WD,WDL,WG,WU,
     '      XA,XE,XG,XID,XIDL,XIG,XIP,XO,XP,YG,YGF,YP,YQ,YQS,
     '      ZA,ZA1,Z_CONT,ZD,ZDL,ZE,ZP,ZP1,
     '      STRING,FIX,ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'GROUP',2)) THEN
          noco=3
          IF(ABBREV(COQU(2,1),'ADD',1)) THEN
            ADD=.TRUE.
          ELSE
            ADD=.FALSE.
          ENDIF
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Elements'
              OPTION_3( 2)='Faces'
              OPTION_3( 3)='Gauss'
              OPTION_3( 4)='Grids'
              OPTION_3( 5)='Lines'
              OPTION_3( 6)='Nodes'
              OPTION_3( 7)='Polylines'
              OPTION_3( 8)='Exit'
              NTCH3=8
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ELEMENTS',1)) THEN
              CALL GRELEM(NEELEM,NEL,NELIST,NELIST2,NLLINE,NLLIST,NORD,
     '          NRLIST,NXI,NFF,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FACES',1)) THEN
              CALL GRFACE(NEELEM,NELIST,NFFACE,NFLIST,NFLIST1,NPF,NFF,
     '          NXI,NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL GRGAUS(GRNGLIST,NBJ,NEELEM,NELIST,
     '          NGLIST,NKJE,NPF,NPNE,NRLIST,NVJE,PG,SE,XA,
     '          XE,XG,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRIDS',1)) THEN
              CALL GRGRID(IBT,IDO,INP,LD,NBJ,NBH,NEELEM,NELIST,
     '          NENQ,NHE,NKHE,NKJE,%VAL(NLATNE_PTR),%VAL(NLATPNQ_PTR),
     '          NPF,NPNE,NQET,NQLIST,%VAL(NQNE_PTR),%VAL(NQNLAT_PTR),
     '          %VAL(NQS_PTR),NQXI,NRE,NRLIST,NVHE,NVJE,NW,NWQ,
     '          NXQ,CURVCORRECT,SE,SQ,XA,XE,XID,XIQ,XP,XQ,ZA,ZD,ZP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',1)) THEN
              CALL GRLINE(NLLINE,NLLIST,NPL,NPLIST1,NRLIST,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
C TVK 10/03/2000 Adding NVCNODE to CALL
              CALL GRNODE(NBJ,NEELEM,NEL,NELIST,NPL,NPNE,NPNODE,NPLIST1,
     &          NPLIST2,NRLIST,NVCNODE,NXI,CE,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINES',1)) THEN
              CALL GRPLIN(STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

C MPN 26Mar2002: changed >print to >fem gxprint
        ELSE IF(ABBREV(CO(2),'GXPRINT',2)) THEN
          noco=2
          CALL PRSCRN(STRING,ERROR,*9999)

C rgb 23/09/99 not used
C        ELSE IF(ABBREV(CO(2),'HELP',2)) THEN
C          noco=3
C          CONTINUE3=.TRUE.
C          DO WHILE(CONTINUE3)
C            IF(LIST_3) THEN
C              OPTION_3( 1)='Character_arrays'
C              OPTION_3( 2)='Commands'
C              OPTION_3( 3)='Common_blocks'
C              OPTION_3( 4)='Examples'
C              OPTION_3( 5)='Files'
C              OPTION_3( 6)='Image_commands'
C              OPTION_3( 7)='Image_trans'
C              OPTION_3( 8)='Integer_arrays'
C              OPTION_3( 9)='Integers'
C              OPTION_3(10)='Modules'
C              OPTION_3(11)='New_features'
C              OPTION_3(12)='Quit'
C              OPTION_3(13)='Real_arrays'
C              OPTION_3(14)='Exit'
C              NTCH3=14
C              IF(DOCUMENT) THEN
C                noch3=noch3+1
C                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
C              ELSE
C                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
C              ENDIF
C            ENDIF
Cc           IF(ABBREV(CO(3),'CHARACTER_ARRAYS',2)) THEN
Cc             CALL DOCUM('help','doc','CHARACTER_ARRAYS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'COMMANDS',5)) THEN
Cc             CALL DOCUM('help','doc','COMMANDS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'COMMON_BLOCKS',5)) THEN
Cc             CALL DOCUM('help','doc','COMMON_BLOCKS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'EXAMPLES',1)) THEN
Cc             CALL DOCUM('help','doc','EXAMPLES',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'FILES',1)) THEN
Cc             CALL DOCUM('help','doc','FILES',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'IMAGE_COMMANDS',7)) THEN
Cc             CALL DOCUM('help','doc','IMAGE_COMMANDS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'IMAGE_TRANS',5)) THEN
Cc             CALL DOCUM('help','doc','IMAGE_TRANS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'INTEGER_ARRAYS',1)) THEN
Cc             CALL DOCUM('help','doc','INTEGER_ARRAYS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'INTEGERS',8)) THEN
Cc             CALL DOCUM('help','doc','INTEGERS',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'MODULES',1)) THEN
Cc             CALL DOCUM('help','doc','MODULES',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'NEW_FEATURES',1)) THEN
Cc             CALL DOCUM('help','doc','NEW_FEATURES',ERROR,*9999)
Cc           ELSE IF(ABBREV(CO(3),'REAL_ARRAYS',1)) THEN
Cc             CALL DOCUM('help','doc','REAL_ARRAYS',ERROR,*9999)
Cc           ELSE
Cc             CALL STRING_TRIM(CO(3),IBEG,IEND)
Cc             WRITE(OP_STRING,'('' >>no help available on '',A)')
Cc    '          CO(3)(IBEG:IEND)
Cc             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
Cc           ENDIF
C            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
C              CONTINUE3=.TRUE.
C            ELSE
C              noch3=0
C              CONTINUE3=.FALSE.
C            ENDIF !document
C          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'HIDE',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Alignment'
              OPTION_3( 2)='Axes'
              OPTION_3( 3)='Boundary'
              OPTION_3( 4)='Clock'
              OPTION_3( 5)='Contours'
              OPTION_3( 6)='Cross-section'
              OPTION_3( 7)='Data'
              OPTION_3( 8)='Dipole'
              OPTION_3( 9)='Elements'
              OPTION_3(10)='Faces'
              OPTION_3(11)='Fibres'
              OPTION_3(12)='Field'
              OPTION_3(13)='Gauss'
              OPTION_3(14)='Gradient'
              OPTION_3(15)='Grid'
              OPTION_3(16)='History'
              OPTION_3(17)='Increments'
              OPTION_3(18)='Isochrones'
              OPTION_3(19)='Lines'
              OPTION_3(20)='Map'
              OPTION_3(21)='Materials'
              OPTION_3(22)='Nodes'
              OPTION_3(23)='Objects'
              OPTION_3(24)='Polyline'
              OPTION_3(25)='Profile'
              OPTION_3(26)='Reactions'
              OPTION_3(27)='Residuals'
              OPTION_3(28)='Rule'
              OPTION_3(29)='Scale'
              OPTION_3(30)='Section'
              OPTION_3(31)='Sheets'
              OPTION_3(32)='Strain'
              OPTION_3(33)='Streamlines'
              OPTION_3(34)='Stress'
              OPTION_3(35)=' '
              OPTION_3(36)='Velocity_field'
              OPTION_3(37)='Exit'
              NTCH3=37
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ALIGNMENT',2)) THEN
              CALL HIALIG(%VAL(ISALIG_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'AXES',2)) THEN
              CALL HIAXES(%VAL(ISAXES_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BOUNDARY',1)) THEN
C             CALL HIBOUN(%VAL(ISBOUN_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CLOCK',2)) THEN
              CALL HICLOC(%VAL(ISCLOC_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONTOURS',2)) THEN
              CALL HICONT(%VAL(ISCONO_PTR),%VAL(ISCONT_PTR),ISEG,NEELEM,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CROSS-SECTION',2)) THEN
              CALL HICROS(%VAL(ISCROS_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',2)) THEN
              CALL HIDATA(%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,NEELEM,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DIPOLE',2)) THEN
              CALL HIDIPO(DIPOLE_LIST,ISEG,%VAL(ISDIPO_PTR),
     &          %VAL(ISDIPA_PTR),NDIPOLES,NRLIST,NXLIST,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL HIELEM(ISEG,%VAL(ISELNO_PTR),%VAL(ISERR_PTR),NEELEM,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FACES',2)) THEN
              CALL HIFACE(ISEG,%VAL(ISFACE_PTR),%VAL(ISFANO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRES',3)) THEN
              CALL HIFIBR(ISEG,%VAL(ISFIBR_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL HIFIEL(ISEG,%VAL(ISFIEL_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL HIGAUS(ISEG,%VAL(ISGAUS_PTR),NBJ,NEELEM,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRADIENT',3)) THEN
              CALL HIGRAD(ISEG,%VAL(ISGRAD_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
              CALL HIGRID(ISEG,%VAL(ISGRID_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HISTORY',2)) THEN
              CALL HIHIST(ISEG,%VAL(ISHIST_PTR),NPNODE,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',2)) THEN
              CALL HIINCR(ISEG,%VAL(ISINCR_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ISOCHRONES',2)) THEN
              CALL HIISOC(ISEG,%VAL(ISISOC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',1)) THEN
              CALL HILINE(ISEG,%VAL(ISLINE_PTR),%VAL(ISLINO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAP',3)) THEN
              CALL HIMAP(ISEG,%VAL(ISMAP_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL HIMATE(ISEG,%VAL(ISMATE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL HINODE(ISEG,%VAL(ISNONO_PTR),NPNODE,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'OBJECTS',1)) THEN
              CALL HIOBJE(ISEG,%VAL(ISOBJE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINE',2)) THEN
              CALL HIPLIN(ISEG,%VAL(ISPLIN_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PROFILE',2)) THEN
              CALL HIPROF(ISEG,%VAL(ISPROF_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REACTIONS',3)) THEN
              CALL HIREAC(ISEG,%VAL(ISREAC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',3)) THEN
              CALL HIRESI(ISEG,%VAL(ISRESI_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RULE',2)) THEN
              CALL HIRULE(ISEG,%VAL(ISRULE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SCALE',2)) THEN
              CALL HISCAL(ISEG,%VAL(ISSCAL_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SECTION',2)) THEN
              CALL HISECT(ISEG,%VAL(ISSECT_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SHEETS',2)) THEN
              CALL HISHEE(ISEG,%VAL(ISSHEE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL HISTRA(ISEG,%VAL(ISSTRA_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STREAMLINES',5)) THEN
              CALL HISTRM(ISEG,%VAL(ISSTRM_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STRESS',5)) THEN
              CALL HISTRE(ISEG,%VAL(ISSTRE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'SURFACE',2)) THEN
              CALL HISURF(ISEG,%VAL(ISSURF_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'VELOCITY_FIELD',1)) THEN
              CALL HIVELO(ISEG,%VAL(ISVELO_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'IMPORT',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Grid'
              OPTION_3( 2)='Signal'
              OPTION_3( 3)='Exit'
              NTCH3=3
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'GRID',1)) THEN
              CALL IMGRID(NBJ,NEELEM,NENQ,NLQ,NQET,NQNE_PTR,NQS_PTR,
     '          NQSCNB,NQXI,NRLIST,NWQ,NXQ,NLATNE_PTR,NLATNQ_PTR,
     '          NLATPNQ_PTR,NQNLAT_PTR,AQ,DNUDXQ,XQ,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SIGNAL',1)) THEN
              CALL IMSIGN(LD,NBJ,NDDATA,NRLIST,NRLIST2,WD,XID,ZD,ZD2,
     '          STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'INQUIRE',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Cell_variable'
              NTCH3=1
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'CELL_VARIABLE',1)) THEN
              CALL INQUIRE_CELL_VARIABLE(CELL_ICQS_NAMES,
     '          CELL_RCQS_NAMES,CELL_YQS_NAMES,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'LIST',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Active'
              OPTION_3( 2)='Analytic'
              OPTION_3( 3)='Bases'
              OPTION_3( 4)='Bib'
              OPTION_3( 5)='Boundary'
              OPTION_3( 6)='Cell'
              OPTION_3( 7)='Colours'
              OPTION_3( 8)='Corners'
              OPTION_3( 9)='Coordinates'
              OPTION_3(10)='Coupling'
              OPTION_3(11)='Customisation'
              OPTION_3(12)='Data'
              OPTION_3(13)='Delaunay'
              OPTION_3(14)='Eigenvalues'
              OPTION_3(15)='Elements'
              OPTION_3(16)='Equations'
              OPTION_3(17)='Export'
              OPTION_3(18)='Faces'
              OPTION_3(19)='Fibres'
              OPTION_3(20)='Fit'
              OPTION_3(21)='Function'
              OPTION_3(22)='Gauss'
              OPTION_3(23)='Grid'
              OPTION_3(24)='Growth'
              OPTION_3(25)='Heading'
              OPTION_3(26)='History'
              OPTION_3(27)='Import'
              OPTION_3(28)='Increments'
              OPTION_3(29)='Initial'
              OPTION_3(30)='Inverse'
              OPTION_3(31)='Labels'
              OPTION_3(32)='Leads'
              OPTION_3(33)='Lines'
              OPTION_3(34)='Materials'
              OPTION_3(35)='Matrix'
              OPTION_3(36)='Mesh'
              OPTION_3(37)='Modal'
              OPTION_3(38)='Motion'
              OPTION_3(39)='Nodes'
              OPTION_3(40)='Noise'
              OPTION_3(41)='Normals'
              OPTION_3(42)='Objects'
              OPTION_3(43)='Optimise'
              OPTION_3(44)='Output'
              OPTION_3(45)='Parameters'
              OPTION_3(46)='Polyline'
              OPTION_3(47)='Record'
              OPTION_3(48)='Reference'
              OPTION_3(49)='Regions'
              OPTION_3(50)='Reg-parameter'
              OPTION_3(51)='Sail'
              OPTION_3(52)='Scale'
              OPTION_3(53)='Singularity'
              OPTION_3(54)='Signal'
              OPTION_3(55)='Segments'
              OPTION_3(56)='Solve'
              OPTION_3(57)='Source'
              OPTION_3(58)='Strain'
              OPTION_3(59)='Stress'
              OPTION_3(60)='Structure'
              OPTION_3(61)='Time_variable'
              OPTION_3(62)='Transfer'
              OPTION_3(63)='User'
              OPTION_3(64)='Variables'
              OPTION_3(65)='Volume'
              OPTION_3(66)='Voronoi'
              OPTION_3(67)='Window'
              OPTION_3(68)='Xi'
              OPTION_3(69)='Exit'
              NTCH3=69
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ACTIVE',2)) THEN
              CALL LIACTI(NBJ,NEELEM,NRLIST,FEXT,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ANALYTIC',2)) THEN
              CALL LIANAL(NPNODE,NRLIST,NXLIST,NYNP,STRING,YP,ERROR,
     '          *9999)
            ELSE IF(ABBREV(CO(3),'BASES',2)) THEN
              CALL LIBASE(IBT,IDO,INP,NAN,NGAP,NKEF,NNF,NNL,
     '          PG,XIG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BOUNDARY',2)) THEN
C!!! LKC there is no call here ?!?!  25-MAR-1998
              CALL ASSERT(1.LT.0,'>>Not implemented',ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BIB',3)) THEN
C              CALL LIBIB(NEELEM,NPNE,NPNODE,NVJE,NVJP,NXI,NYNP,XP,
C     '          YP(1,1,1),ZA,FIX(1,1,1),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CELL',2)) THEN
              CALL LICELL(CELL_ICQS_SPATIAL,CELL_ICQS_VALUE,
     '          CELL_RCQS_SPATIAL,CELL_YQS_SPATIAL,NEELEM,NPNODE,
     '          NRLIST,NXLIST,CE,
     '          CELL_RCQS_VALUE,CELL_YQS_VALUE,CELL_ICQS_NAMES,
     '          CELL_RCQS_NAMES,CELL_YQS_NAMES,CP,CQ,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COLOURS',3)) THEN
              CALL LICOLO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COORDINATES',3)) THEN
              CALL LICOOR(NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CORNERS',3)) THEN
              CALL LICORN(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COUPLING',3)) THEN
              CALL LICOUP(NEELEM,NP_INTERFACE,NPNODE,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CUSTOMISATION',3)) THEN
              CALL LICUST(NRLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',2)) THEN
              CALL LIDATA(IBT,IDO,INP,LD,NBH,NBJ,NDDL,NDLT,NDP,
     '          NEELEM,NELIST,NFLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRE,NRLIST,NVHE,NVHP,NVJE,NYNE,
     '          NYNP,CURVCORRECT,EDD,SE,SQ,WD,XA,XE,XID,XP,YP,
     '          ZA,ZD,ZE,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DELAUNAY',2)) THEN
              CALL LIDELA(NBJ,NEELEM,NPNE,NVJE,NXI,XP,ZA,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EIGENVALUES',2)) THEN
              CALL LIEIGE(EIGVEC,YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',2)) THEN
              CALL LIELEM(IBT,NBH,NBJ,NEELEM,NELIST,
     '          NFF,NHE,NHP,NKH,NKHE,NKJE,
     '          NLL,NP_INTERFACE,NPF,NPLIST1,NPNE,NPNODE,NRE,NRLIST,
     '          NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,VOLTC,
     '          CURVCORRECT,PG,RG,SE,VOL,VOLT,WG,XA,XAB,XE,XG,XP,YP,
     '          ZA,ZE,ZG,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EQUATIONS',2)) THEN
              CALL LIEQUA(NBH,NEELEM,NELIST,NHE,NHP,NKH,
     '          NPL,NPNODE,NRLIST,NVHE,NVHP,NW,NXLIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'EXPORT',2)) THEN
              CALL LIEXPO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FACES',2)) THEN
              CALL LIFACE(NBJ,NBJF,NLF,NNF,NPF,NPNE,NPNF,NRLIST,
     '          DF,SF,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRE',3)) THEN
C!!! LKC there is no call here ?!?!  8-MAY-1998
              CALL ASSERT(1.LT.0,'>>Not implemented',ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIT',3)) THEN
              CALL LIFIT(NBJ,NEELEM,NKH,NMNO,NPNODE,NRLIST,
     '          NVHP,NXLIST,NYNP,WU,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FUNCTION',2)) THEN
              CALL LIFUNC(IBT,IDO,INP,LD,NBH,NBJ,NEELEM,NHE,NHP,
     '          NKH,NKHE,NKJE,NPF,NPNE,NPNODE,NVHE,NVHP,NVJE,NW,NYNE,
     '          NYNP,CURVCORRECT,SE,XA,XE,XID,XP,YP,ZA,ZE,ZP,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL LIGAUS(NBH,NBJ,NEELEM,NELIST,NGLIST,NHE,NHP,NKH,
     '          NKHE,NKJE,NPF,NPNE,NPNODE,NRLIST,
     '          NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,CURVCORRECT,
     '          PG,RG,SE,XA,XE,XG,XIG,XP,YG,YP,ZA,ZE,ZG,ZP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
              CALL LIGRID(NQLIST,NRLIST,NAQ,NEELEM,NENQ,
     '          %VAL(NLATNE_PTR),%VAL(NLATNQ_PTR),%VAL(NLATPNQ_PTR),NLQ,
     '          NQET,NQGP,NQGP_PIVOT,%VAL(NQNE_PTR),%VAL(NQNLAT_PTR),
     '          NQNP,%VAL(NQS_PTR),NQSCNB,NQXI,NWQ,NXQ,NXLIST,AQ,CQ,
     '          DNUDXQ,DXDXIQ,DXDXIQ2,GCHQ,GUQ,PROPQ,XQ,YQ,YQS,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GROWTH',3)) THEN
              CALL LIGROW(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HEADING',2)) THEN
              CALL LIHEAD(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HISTORY',2)) THEN
              CALL LIHIST(NBH,NEELEM,NHE,NHP,NKH,
     '          NPLIST1,NPNODE,NVHP,NYNE,NYNP,
     '          YP,ZA,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'IMPORT',2)) THEN
              CALL LIIMPO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',3)) THEN
              CALL LIINCR(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INITIAL',3)) THEN
              CALL LIINIT(ITHRES,NBH,NBJ,NEELEM,NHE,NHP,NKH,NODENVCB,
     '          NPNODE,NRLIST,NVCB,NVHP,NW,
     '          NWQ,NXLIST,NYNE,NYNP,AQ,THRES,
     '          YG,YP,YQ,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INVERSE',3)) THEN
              CALL LIINVE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LABELS',3)) THEN
              CALL LIST_LABELS(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LEADS',2)) THEN
              CALL LILEAD(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',2)) THEN
              CALL LILINE(NEL,NLLINE,NPL,DL,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAPPING',3)) THEN
              CALL LIMAP(IBT,IDO,INP,NBH,NBJ,
     '          NEELEM,NENP,NHE,NHP,NKB,NKHE,NKH,NKJ,NKJE,
     '          NLL,NNB,NNF,NNL,NONY,NPF,NPL,
     '          NPNE,NPNODE,NPNY,NRLIST,NVHE,NVHP,NVJE,NVJP,
     '          NWP,NXI,NXLIST,NYNE,NYNO,
     '          NYNP,NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,STRING,FIX,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATRIX',4)) THEN
              CALL LIMATR(ISC_GD,ISC_GKK,ISC_GM,ISC_GMM,ISIZE_MFI,
     '          ISIZE_PHI,ISIZE_TBH,ISR_GD,ISR_GKK,ISR_GM,
     '          ISR_GMM,LD_NP,NRLIST,NXLIST,GD,GKK,GM,GMM,GR,GRR,
     '          MFI,PHI,PHI_H,
     '          PHI_H_EXACT,T_BH,T_BH_INV,YP,ZCROSSING,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL LIMATE(ILPIN,LIST,NBJ,NEELEM,NMBIN,
     '          NPNODE,NRLIST,NW,NXLIST,
     '          CE,CGE,CIN,CP,CQ,YG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MESH',2)) THEN
              CALL LIMESH(NBJ,NEELEM,NELIST,NENP,NORD,NPLIST1,NPNE,
     &          NPNODE,NRLIST,NVJE,NXI,NXLIST,NYNP,CE,CP,XP,YP,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MODAL',2)) THEN
              CALL LIMODA(LIST,NRLIST,NXLIST,
     '          EIGVAL,EIGVEC,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MOTION',2)) THEN
              CALL LIMOTI(IBT,NBH,NEELEM,NHP,NKH,NPNODE,
     '          YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',2)) THEN
              YP1_PTR=0
              CALL ALLOCATE_MEMORY(NYM*NIYM*NXM,1,DPTYPE,YP1_PTR,
     '          .TRUE.,ERROR,*9999)
              CALL LINODE(NBH,NEELEM,NHE,NHP,NKH,NKJ,
     '          NP_INTERFACE,NPNODE,NPLIST1,NRLIST,NVHP,
     '          NVJP,NWP,NXLIST,NYNE,NYNP,XP,YP,%VAL(YP1_PTR),ZA,ZP,
     '          STRING,ERROR,*9999)
              CALL FREE_MEMORY(YP1_PTR,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NOISE',3)) THEN
              CALL LINOIS(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NORMALS',3)) THEN
              CALL LINORM(NEELEM,NRLIST,NW,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OBJECTS',2)) THEN
              CALL LIOBJE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OPTIMISE',2)) THEN
              CALL LIOPTI(LDR,NLNO,NMNO,NP1OPT,NP2OPT,NP3OPT,
     '          NRLIST,NXLIST,NYNO,PAOPTY,REALWORK(OS_CM),
     '          REALWORK(OS_CONTR),REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PMIN),REALWORK(OS_PMAX),
     '          REALWORK(OS_RESID),REALWORK(OS_RESIDM),
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OUTPUT',2)) THEN
              CALL LIOUTP(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,
     '          NEELEM,NELIST,NENP,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,
     '          NKH,NKHE,NKJE,NLL,NNB,NNF,NP_INTERFACE,NPF,NPNE,
     '          NPNODE,NPNY,NRE,NRLIST,NSB,NVHE,NVHP,NVJE,NVJP,NW,
     '          NXI,NXLIST,NYNE,NYNP,NYNR,Z_CONT_LIST,
     '          CE,CG,CGE,CP,CURVCORRECT,DET,DL,DRDN,FEXT,
     '          PG,RAD,RD,RE1,
     '          RG,SE,WG,XA,XE,XG,XG1,XIG,XN,XP,XR,YD,YG,YGF,
     '          YP,ZA,ZA1,Z_CONT,ZE,ZF,ZG,
     '          ZP,ZP1,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PARAMETERS',2)) THEN
              CALL LIPARA(NEELEM,NENP,NFFACE,NLLINE,NPNODE,NXLIST,
     '          NYNR,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINE',2)) THEN
              CALL LIPLIN(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RECORD',3)) THEN
C old MPN unused?
c              CALL LIRECO(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REFERENCE',3)) THEN
              CALL LIREFE(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REGIONS',4)) THEN
              CALL LIREGI(NEELEM,NPNODE,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REG-PARAMETER',4)) THEN
              CALL LIREGP(REG_PARAMETER,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SAIL',2)) THEN
              CALL LISAIL(XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SCALE',2)) THEN
            ELSE IF(ABBREV(CO(3),'SEGMENTS',2)) THEN
              CALL LISEGM(%VAL(ISAXES_PTR),%VAL(ISBASE_PTR),
     &          %VAL(ISCONO_PTR),%VAL(ISCONT_PTR),%VAL(ISDANO_PTR),
     &          %VAL(ISDAPR_PTR),%VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,
     &          %VAL(ISELNO_PTR),%VAL(ISFACE_PTR),%VAL(ISFANO_PTR),
     &          %VAL(ISFIBR_PTR),%VAL(ISGAUS_PTR),%VAL(ISGRAD_PTR),
     &          %VAL(ISGRID_PTR),%VAL(ISHIST_PTR),%VAL(ISINCR_PTR),
     &          %VAL(ISLINE_PTR),%VAL(ISLINO_PTR),%VAL(ISMAP_PTR),
     &          %VAL(ISMATE_PTR),%VAL(ISNONO_PTR),%VAL(ISREAC_PTR),
     &          %VAL(ISRULE_PTR),%VAL(ISSECT_PTR),%VAL(ISSTRA_PTR),
     &          %VAL(ISSTRE_PTR),%VAL(ISSTRM_PTR),%VAL(ISSURF_PTR),
     &          %VAL(ISVELO_PTR),NEELEM,NPNODE,CSEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SINGULARITY',3)) THEN
              CALL LISING(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SIGNAL',3)) THEN
              CALL LISIGN(LD,NBJ,NDDATA,WD,XID,ZD,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOLVE',3)) THEN
              CALL LISOLV(NPNY,NRLIST,NXLIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOURCE',3)) THEN
              CALL LISOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,NDIPOLES,
     '          NRLIST,NXLIST,DIPOLE_CEN,DIPOLE_DIR,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL LISTRA(IBT,IDO,INP,NAN,NBH,NBJ,NDDL,NDLT,
     '          NEELEM,NELIST,NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,
     '          NPNODE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '          CE,CG,CGE,CP,CURVCORRECT,PG,RG,SE,XA,XE,XG,
     '          XID,XIG,XP,YP,ZA,ZE,ZG,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRESS',4)) THEN
              CALL LISTRE(IBT,IDO,INP,NAN,NBH,NBJ,
     '          NDDL,NDLT,NEELEM,NELIST,
     '          NGLIST,NHE,NHP,NKH,NKHE,NKJE,NPF,NPNE,NPNODE,
     '          NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,
     '          CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,
     '          RG,SE,WG,XA,XE,
     '          XG,XID,XIG,XP,YG,YP,ZA,ZE,ZG,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRUCTURE',4)) THEN
c             CALL LISTRU(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TIME_VARIABLE',4)) THEN
C              CALL LITIME(STRING,ERROR,*9999)
              CALL LITIME(NTIME_INTERP,NTIME_POINTS,TIME_VALUES,STRING,
     '          TIME_VARIABLE_NAMES,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'TRANSFER',3)) THEN
              CALL LITRSF(NPLIST3,NPLIST4,NPLIST5,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VARIABLES',2)) THEN
C              CALL LIVARI(NBH,NEELEM,NHE,NHP,NKH,NPNODE,NP_INTERFACE,
C     '          NQNP,NRLIST,NWQ,NVHP,NXLIST,NXQ,NYNE,NYNP,NYNR,XQ,YP,
C     '          YQ,ZA,ZP,STRING,FIX,FIXQ,ERROR,*9999)
              CALL LIVARI(NBH,NEELEM,NENQ,NHE,NHP,
     '          NKH,NPNODE,NP_INTERFACE,NQNP,%VAL(NQS_PTR),NQXI,
     '          NRLIST,NVHP,NXLIST,NXQ,NYNE,NYNP,NYNR,AQ,CQ,
     '          DNUDXQ,DXDXIQ,DXDXIQ2,XQ,YP,YQ,ZA,ZP,STRING,
     '          FIX,FIXQ,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VOLUME',3)) THEN
               CALL LIVOLU(NBJ,NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,
     '          NVJE,NW,PG,RG,SE,WG,XA,XE,XG,XN,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VORONOI',3)) THEN
               CALL LIVORO(NFVC,NODENVC,NODENVCB,NPNODE,NVCB,NVCNODE,
     '          VC,XNFV,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'WINDOW',1)) THEN
              CALL LIWIND(STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'XI',1)) THEN
C LKC 25-MAY-1998 rewrite
C              CALL LIXI(NPNODE,NRLIST,XP,STRING,ERROR,*9999)
C GBS 26-Oct-2000 Redone again
C              CALL LIXI(NEP,NPNODE,LD,XID,XIP,STRING,ERROR,*9999)
              CALL LIXI(NEP,NPNODE,NRLIST,LD,XID,XIP,STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'OPEN',1)) THEN
C cpb 17/7/96 Old
C          noco=3
C          CONTINUE3=.TRUE.
C          DO WHILE(CONTINUE3)
C            IF(LIST_3) THEN
C              OPTION_3( 1)='Metafile'
C              OPTION_3( 2)='Postcript'
C              OPTION_3( 3)='Exit'
C              NTCH3=3
C              IF(DOCUMENT) THEN
C                noch3=noch3+1
C                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
C              ELSE
C                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
C              ENDIF
C            ENDIF
C            IF(ABBREV(CO(4),'?',1)) THEN
C              CALL STRING_TRIM(STRING,IBEG,IEND)
C              WRITE(OP_STRING,'(1X,A)') STRING(1:IEND)//' <on WS_ID>[1]'
C             CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
C            ELSE
C              IF(CBBREV(CO,'ON',1,noco+1,NTCO,N3CO)) THEN
C                iw=IFROMC(CO(N3CO+1))
C              ELSE
C                iw=1
C              ENDIF
C              IF(ABBREV(CO(3),'METAFILE',1)) THEN
Cc               CALL OPEN_PRINT_FILE(iw,'METAFILE',ERROR,*9999)
C              ELSE IF(ABBREV(CO(3),'POSTSCRIPT',1)) THEN
CC                IF(CBBREV(CO,'LANDSCAPE',1,noco+1,NTCO,N3CO)) THEN
CC                  PORTRT=.FALSE.
CC                ELSE
CC                  PORTRT=.TRUE.
CC                ENDIF
Cc               CALL OPEN_PRINT_FILE(iw,'POSTSCRIPT',ERROR,*9999)
C              ELSE
C                IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
C              ENDIF
C              IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
C                CONTINUE3=.TRUE.
C              ELSE
C                noch3=0
C                CONTINUE3=.FALSE.
C              ENDIF !document
C            ENDIF
C          ENDDO !continue3
          noco=2
          CALL OPEN_FILES(LD,NBJ,NDDATA,NHQ,NPNY,NQNY,NRLIST,NRLIST2,
     '      NXLIST,NYNR,NYQNR,WD,XID,YP,YQ,YQS,ZD,STRING,ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'READ',3)) THEN
          noco=2
          CALL READF(IBT,IDO,INP,ISC_GKK,ISIZE_MFI,
     '      ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,ISR_GKK,
     '      ITHRES,LD_NP,MAP_ART_VEIN,NAN,NBH,NBJ,NBJF,NDET,NEELEM,
     &      NELIST,NENP,NEL,NFF,NGAP,NHE,NHP,NHQ,NKJE,NKEF,NKH,NKJ,NLF,
     &      NLL,NNF,NNL,NONY,NPF,NPL,NPNE,NPNODE,NPNY,NP_INTERFACE,
     &      %VAL(NQNE_PTR),NQNY,%VAL(NQS_PTR),NQXI,NRE,NRLIST,NRLIST2,
     &      NVHE,NVHP,NVJE,NVJP,NW,NWQ,NXI,NXLIST,NYNE,NYNO,NYNP,NYNR,
     &      NYQNR,CE,CONY,CP,DET,DL,FEXT,GKK,GR,GRR,
     '      MFI,PE,PF,PG,PHI,PHI_H,PHI_H_EXACT,SE,T_BH,T_BH_INV,
     '      THRES,WG,XA,XIG,XP,YP,YQ,YQS,ZCROSSING,STRING,FIX,FIXP,
     '      ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'REALLOCATE',3)) THEN
C#### Command: FEM reallocate
C###  Description:
C###    Reallocates arrays dimensioned by include parameters.
          REALLOCATE_FEM=.TRUE.
          REALLOCATE_SYNTAX=.TRUE.
          MAXNITW=6 ! max # windows that can be parsed
          DO iw=1,MAXNITW
            IF(IWKG(iw).EQ.1) THEN
              IWKG(iw)=0
              CALL CLOSE_WS(iw,ERROR,*9999)
            ENDIF
          ENDDO !iw

        ELSE IF(ABBREV(CO(2),'REFINE',3)) THEN
          noco=2
C LKC 5-NOV-97 not ref ,NPLIST1
          CALL REFINE(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,
     '      NELIST,NELIST2,NENP,NFF,NFFACE,NHE,NHP,NKEF,NKH,NKHE,NKJ,
     '      NKJE,NLF,NLL,NLLINE,NNB,NNF,NNL,NONY,NPF,NP_INTERFACE,NPL,
     '      NPLIST1,NPNE,NPNF,NPNODE,NPNY,NRE,NRLIST,NUNK,NVHE,NVHP,
     '      NVJE,NVJF,NVJL,NVJP,NW,NWP,NXI,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '      CE,CONY,CYNO,DF,
     '      DL,PG,RG,SE,SF,SP,WG,XA,XE,XG,XP,YP,ZA,FIX,
     '      STRING,ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'RENUMBER',3)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3(1)='Mesh'
              NTCH3=1
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'MESH',1)) THEN
              CALL RNMESH(IBT,IDO,INP,NBJ,NBJF,NEELEM,NEL,NENP,
     '          NFF,NFFACE,NKJE,NKEF,NKJ,NLF,NLL,NLLINE,NNB,NNF,
     '          NNL,NPF,NPL,NPNE,NPNF,NPNODE,NRE,NRLIST,NVJE,
     '          NVJF,NVJL,NVJP,NXI,DF,DL,PG,RG,SE,SF,STRING,WG,XA,
     '          XE,XG,XP,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

C        ELSE IF(ABBREV(CO(2),'RESET',5)) THEN
C MHT archived 18-Feb-99


C LKC 6-OCT-1999 Archived fiducial
C
C        ELSE IF(ABBREV(CO(2),'SEND',2)) THEN
C          noco=3
C          CONTINUE3=.TRUE.
C          DO WHILE(CONTINUE3)
C            IF(LIST_3) THEN
C              OPTION_3( 1)='Fiducial'
C              OPTION_3( 2)='Exit'
C              NTCH3=2
C              IF(DOCUMENT) THEN
C                noch3=noch3+1
C                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
C              ELSE
C                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
C              ENDIF
C            ENDIF
C            IF(ABBREV(CO(3),'FIDUCIAL',1)) THEN
CC#### Command: FEM send fiducial
CC###  Description:
CC###    Send information about the fiducial markers from a Motif
CC###    front end, through the socket interface.
C              CALL SENDFID(STRING,ERROR,*9999)
C            ELSE
C              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
C            ENDIF
C            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
C              CONTINUE3=.TRUE.
C            ELSE
C              noch3=0
C              CONTINUE3=.FALSE.
C            ENDIF !document
C          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'SHAPE',3)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Chords'
              OPTION_3( 2)='Exit'
              NTCH3=2
              CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'CHORDS',1)) THEN
              CALL SHAPE_CHORDS(ISEG,%VAL(ISLINO_PTR),NLCHOR,CSEG,
     &          STRING,ERROR,*9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'SHOW',2)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Alignment'
              OPTION_3( 2)='Axes'
              OPTION_3( 3)='Boundary'
              OPTION_3( 4)='Clock'
              OPTION_3( 5)='Contours'
              OPTION_3( 6)='Cross-section'
              OPTION_3( 7)='Data'
              OPTION_3( 8)='Dipole'
              OPTION_3( 9)='Elements'
              OPTION_3(10)='Faces'
              OPTION_3(11)='Fibres'
              OPTION_3(12)='Field'
              OPTION_3(13)='Gauss'
              OPTION_3(14)='Gradient'
              OPTION_3(15)='Grid'
              OPTION_3(16)='History'
              OPTION_3(17)='Increments'
              OPTION_3(18)='Isochrones'
              OPTION_3(19)='Lines'
              OPTION_3(20)='Map'
              OPTION_3(21)='Materials'
              OPTION_3(22)='Nodes'
              OPTION_3(23)='Objects'
              OPTION_3(24)='Polyline'
              OPTION_3(25)='Profile'
              OPTION_3(26)='Reactions'
              OPTION_3(27)='Residuals'
              OPTION_3(28)='Rule'
              OPTION_3(29)='Scale'
              OPTION_3(30)='Section'
              OPTION_3(31)='Sheets'
              OPTION_3(32)='Strain'
              OPTION_3(33)='Streamlines'
              OPTION_3(34)='Stress'
              OPTION_3(35)='Surface'
              OPTION_3(36)='Velocity_field'
              OPTION_3(37)='Exit'
              NTCH3=38
              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'ALIGNMENT',2)) THEN
              CALL SHALIG(%VAL(ISALIG_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'AXES',2)) THEN
              CALL SHAXES(%VAL(ISAXES_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'BOUNDARY',1)) THEN
C             CALL SHBOUN(%VAL(ISBOUN_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CLOCK',1)) THEN
              CALL SHCLOC(%VAL(ISCLOC_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CONTOURS',2)) THEN
              CALL SHCONT(%VAL(ISCONO_PTR),%VAL(ISCONT_PTR),ISEG,NEELEM,
     &          NELIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'CROSS-SECTION',2)) THEN
              CALL SHCROS(%VAL(ISCROS_PTR),ISEG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',2)) THEN
              CALL SHDATA(%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,NEELEM,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DIPOLE',2)) THEN
              CALL SHDIPO(DIPOLE_LIST,ISEG,%VAL(ISDIPO_PTR),
     &          %VAL(ISDIPA_PTR),NDIPOLES,NRLIST,NXLIST,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL SHELEM(ISEG,%VAL(ISELNO_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'FACES',2)) THEN
              CALL SHFACE(ISEG,%VAL(ISFACE_PTR),%VAL(ISFANO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRES',3)) THEN
              CALL SHFIBR(ISEG,%VAL(ISFIBR_PTR),NEELEM,NELIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL SHFIEL(ISEG,%VAL(ISFIEL_PTR),NEELEM,NELIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL SHGAUS(ISEG,%VAL(ISGAUS_PTR),NBJ,NEELEM,NELIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRADIENT',3)) THEN
              CALL SHGRAD(ISEG,%VAL(ISGRAD_PTR),NEELEM,NELIST,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
              CALL SHGRID(ISEG,%VAL(ISGRID_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'HISTORY',2)) THEN
              CALL SHHIST(ISEG,%VAL(ISHIST_PTR),NPLIST1,NPNODE,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INCREMENTS',2)) THEN
              CALL SHINCR(ISEG,%VAL(ISINCR_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'%VAL(ISOCHRONES_PTR)',2)) THEN
              CALL SHISOC(ISEG,%VAL(ISISOC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',1)) THEN
              CALL SHLINE(ISEG,%VAL(ISLINE_PTR),%VAL(ISLINO_PTR),STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MAP',3)) THEN
              CALL SHMAP(ISEG,%VAL(ISMAP_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIALS',3)) THEN
              CALL SHMATE(ISEG,%VAL(ISMATE_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'NODES',1)) THEN
              CALL SHNODE(ISEG,%VAL(ISNONO_PTR),NPNODE,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'OBJECTS',1)) THEN
              CALL SHOBJE(ISEG,%VAL(ISOBJE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'POLYLINE',2)) THEN
              CALL SHPLIN(ISEG,%VAL(ISPLIN_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PROFILE',2)) THEN
              CALL SHPROF(ISEG,%VAL(ISPROF_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'REACTIONS',3)) THEN
              CALL SHREAC(ISEG,%VAL(ISREAC_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',3)) THEN
              CALL SHRESI(ISEG,%VAL(ISRESI_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RULE',2)) THEN
              CALL SHRULE(ISEG,%VAL(ISRULE_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SCALE',2)) THEN
              CALL SHSCAL(ISEG,%VAL(ISSCAL_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SECTION',2)) THEN
              CALL SHSECT(ISEG,%VAL(ISSECT_PTR),STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SHEETS',2)) THEN
              CALL SHSHEE(ISEG,%VAL(ISSHEE_PTR),NEELEM,NELIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STRAIN',4)) THEN
              CALL SHSTRA(ISEG,%VAL(ISSTRA_PTR),NEELEM,NELIST,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'STREAMLINES',5)) THEN
              CALL SHSTRM(ISEG,%VAL(ISSTRM_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'STRESS',5)) THEN
              CALL SHSTRE(ISEG,%VAL(ISSTRE_PTR),NEELEM,NELIST,STRING,
     &          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SURFACE',1)) THEN
              CALL SHSURF(ISEG,%VAL(ISSURF_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'VELOCITY_FIELD',1)) THEN
              CALL SHVELO(ISEG,%VAL(ISVELO_PTR),NEELEM,STRING,ERROR,
     &          *9999)
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'SOLVE',2)) THEN
          noco=2
          CALL SOLVE(CELL_ICQS_VALUE,DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,
     '      IBT,ICQS,ICQS_SPATIAL,IDO,IICQS_SPATIAL,IRCQS_SPATIAL,INP,
     '      ISC_GD,ISC_GKK,ISC_GM,ISC_GMM,ISR_GD,ISR_GKK,ISR_GM,ISR_GMM,
     '      ISEG,%VAL(ISELNO_PTR),%VAL(ISFIBR_PTR),%VAL(ISFIEL_PTR),
     &      %VAL(ISLINE_PTR),%VAL(ISLINO_PTR),%VAL(ISNONO_PTR),
     &      INTWORK(OS_ISTATE),INTWORK,LGE,MXI,NAN,NAQ,NBH,NBHF,NBJ,
     &      NBJF,NDET,NDIPOLES,NEELEM,NELIST,NENFVC,NENP,NENQ,
     &      NEP,NFF,NFFACE,NFVC,NGAP,NHE,NHP,NHQ,NKB,NKEF,NKH,NKHE,NKJE,
     &      %VAL(NLATNE_PTR),%VAL(NLATNQ_PTR),%VAL(NLATPNQ_PTR),NLL,
     &      NLLIST,NLQ,NNB,NNF,NODENVC,NODENVCB,NONY,NORD,NMNO,
     &      NP_INTERFACE,NPB,NPF,NPL,NPLIST1,NPLIST2,NPNE,NPNODE,NPNY,
     &      NQET,NQGP,NQGP_PIVOT,NQGW,NQLIST,%VAL(NQNE_PTR),
     '      %VAL(NQNLAT_PTR),NQNP,NQNY,%VAL(NQS_PTR),NQSCNB,NQXI,NRE,
     '      NRLIST,NRLIST2,NSB,NTIME_INTERP,NTIME_POINTS,NTIME_NR,NVCB,
     '      NVCNODE,NVHE,NVHP,NVJE,NVJP,NW,NWQ,NXI,NXLIST,NXQ,NYNE,
     '      NYNO,NYNP,NYNQ,NYNR,NYQNR,TV_BC_SET,Z_CONT_LIST,ACINUS,
     '      AQ,BBM,CE,CELL_RCQS_VALUE,CG,CGE,
     '      CONY,CP,CQ,CURVCORRECT,CYNO,DET,
     '      DIPOLE_CEN,DIPOLE_DIR,DL,DNUDXQ,DRDN,DRDNO,
     '      DXDXIQ,DXDXIQ2,ED,EIGVAL,EIGVEC,EM,ER,ES,FEXT,GCHQ,GD,
     '      GKK,GM,GMM,GR,GRR,GUQ,REALWORK(OS_PAOPTI),PG,
     '      REALWORK(OS_PMIN),REALWORK(OS_PMAX),PROPQ,REALWORK(OS_R),
     '      RAD,RCQS,RCQS_SPATIAL,RD,RE1,REALWORK(OS_F),
     '      REALWORK(OS_RESJAC),RG,RHS,SE,TIME_VALUES,REALWORK,VC,
     '      VC_INIT,WG,XA,XAB,REALWORK(OS_XC),XE,XG,XG1,XIG,XIP,XIQ,XN,
     &      XNFV,XN_GRAD,XO,XP,XQ,XR,XR_GRAD,YG,YGF,YP,YQ,YQS,ZA,ZA1,
     &      Z_CONT,ZE,ZG,ZNFV,ZP,ZP1,CSEG,STRING,FIX,FIXQ,ERROR,
     &      *9999)

        ELSE IF(ABBREV(CO(2),'STEP',2)) THEN
          noco=2
          CALL STEP(IBT,%VAL(ISCLOC_PTR),ISEG,%VAL(ISLINE_PTR),
     &      %VAL(ISLINO_PTR),NBH,NEELEM,NHE,NHP,NKH,NKJ,NLLIST,NPL,
     &      NPNODE,NVHP,NXLIST,NYNE,NYNP,DL,XP,YP,ZA,ZP,CSEG,STRING,
     &      ERROR,*9999)

        ELSE IF(ABBREV(CO(2),'TRACK',4)) THEN
C#### Command: FEM track
C###  Description:
C###    Track current lines through solution domain.
          noco=2
          CALL TRACK(NEELEM,NLL,NRLIST,NXLIST,DL,STRING,ERROR,*9999)
c##### Command fem truncate 
C##### Description Truncate a 1D airway mesh to simplify solution            
        ELSEIF(ABBREV(CO(2),'TRUNCATE',3)) THEN
          CALL TRUNCATE(NBJ,NDP,NEELEM,NELIST,NENP,NORD,NPNE,
     &      NPNODE,NRLIST,NVJE,NVJP,NXI,BBM,CE,XAB,XP,
     &      ZD,ERROR,*9999)  

        ELSE IF(ABBREV(CO(2),'UPDATE',1)) THEN
          noco=3
          CONTINUE3=.TRUE.
          DO WHILE(CONTINUE3)
            IF(LIST_3) THEN
              OPTION_3( 1)='Aerofoil'
              OPTION_3( 2)='Coupling'
              OPTION_3( 3)='Data'
              OPTION_3( 4)='Delaunay'
              OPTION_3( 5)='Elements'
              OPTION_3( 6)='Fibre'
              OPTION_3( 7)='Field'
              OPTION_3( 8)='Flux'
              OPTION_3( 9)='Gauss'
              OPTION_3(10)='Geometry'
              OPTION_3(11)='Grid'
              OPTION_3(12)='Growth'
              OPTION_3(13)='Initial'
              OPTION_3(14)='Initial_Voronoi_volume'
              OPTION_3(15)='Mass'
              OPTION_3(16)='Material'
              OPTION_3(17)='Mesh'
              OPTION_3(18)='MFI'
              OPTION_3(19)='Motion' !AJS 11/2010
              OPTION_3(20)='Nodes'
              OPTION_3(21)='Optimisation'
              OPTION_3(22)='Ordering'
              OPTION_3(23)='PHI'
              OPTION_3(24)='Pressure'
              OPTION_3(25)='Residuals'
              OPTION_3(26)='Scale_factors'
              OPTION_3(27)='Signal'
              OPTION_3(28)='Sobolev'
              OPTION_3(29)='Solution'
              OPTION_3(30)='Source'
              OPTION_3(31)='Variable'
              OPTION_3(32)='Vertices'
              OPTION_3(33)='View'
              OPTION_3(34)='Voronoi'
              OPTION_3(35)='Xi'
              OPTION_3(36)='Exit'
              NTCH3=36

              IF(DOCUMENT) THEN
                noch3=noch3+1
                CALL COM_DOC3(CO,OPTION_3(noch3),STRING,ERROR,*9999)
              ELSE
                CALL LIST_COMMANDS(3,NTCH3,OPTION_3,ERROR,*9999)
              ENDIF
            ENDIF
            IF(ABBREV(CO(3),'AEROFOIL',1)) THEN
              CALL UPAERO(NPL,NPNODE,NYNP,NYNR,
     '          DL,XP,YP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'COUPLING',3)) THEN
              CALL UPCOUP(NBH,NEELEM,NHE,NHP,NKH,NKJ,NPNODE,
     &          NP_INTERFACE,NVHP,NVJP,NYNE,NYNP,PF,XP,YP,ZA,ZP,STRING,
     &          ERROR,*9999)
c MHT new call to upcoup in progress 13-03-06              
c              CALL UPCOUP(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NEP,NHE,NHP,
c     &          NKH,NKHE,NKJ,NKJE,NPF,NPNE,NPNODE,NP_INTERFACE,NVHE,
c     &          NVHP,NVJE,NVJP,NW,NXLIST,NYNE,NYNP,CE,CG,CGE,CP,
c     &          CURVCORRECT,FEXT,PF,PG,RG,SE,XA,XE,XG,XIP,XP,YG,YP,ZA,
c     &          ZE,ZG,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'DATA',1)) THEN
              CALL UPDATA(FD,IBT,IDO,INP,ISIZE_MFI,LD,LD_NP,NAN,NBH,NBJ,
     &          NBJF,NDDATA,NEELEM,NELIST,NFF,NHE,NKEF,NKHE,NKJE,NLL,
     &          NNF,NPF,NPNE,NPNF,NVHE,NVJE,NVJF,NRE,NRLIST,NW,
     &          Z_CONT_LIST,CURVCORRECT,MFI,PG,RG,SE,SF,SQ,WD,WG,XA,XAB,
     &          XE,XG,XID,XP,YQS,ZA,Z_CONT,ZD,ZE,ZG,ZP,STRING,ERROR,
     &          *9999)
            ELSE IF(ABBREV(CO(3),'DELAUNAY',1)) THEN
              CALL UPDELA(IBT,NBJ,NEELEM,NELIST,NENP,NKJE,NPLIST1,NPNE,
     '          NPNODE,NRE,NVJE,NVJP,NXI,SE,XP,ZA,STRING,ERROR,*9999)
C               .. Update lines and faces if not using Voronoi cells
              IF(USE_VORONOI.NE.1) THEN
                CALL NENXI(NBJ,NEELEM,NENP,NNB,NPNE,NXI,ERROR,*9999)
                CALL LINCAL(IBT,IDO,INP,0,NBJ,NEELEM,NEL,NENP,NKJE,NLL,
     '            NLLINE,NNL,NPL,NPNE,NPNODE,NRE,NVJE,NVJL,DL,SE,XP,
     '            ERROR,*9999)
                CALL FACCAL(IBT,IDO,INP,NBJ,NBJF,NEELEM,NFF,NFFACE,
     '            NKJE,NKEF,NLF,NNF,NNL,NPF,NPL,NPNE,NPNF,NRE,NVJE,
     '            NVJF,DF,PG,RG,SE,SF,WG,XA,XE,XG,XP,ERROR,*9999)
              ENDIF
            ELSE IF(ABBREV(CO(3),'ELEMENTS',3)) THEN
              CALL UPELEM(IBT,IDO,INP,NBJ,NEELEM,NENQ,NKJE,NPF,NPNE,
     '          NQLIST,NRLIST,NVJE,NWQ,NXLIST,NYNP,SE,XA,XE,XIQ,XP,XQ,
     '          STRING,FIX,FIXQ,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIBRE',3)) THEN
              CALL UPFIBR(CP,NKH,NKJ,NPLIST1,NPNODE,NRLIST,NVHP,
     '          NVJP,NXLIST,NYNP,XP,YP,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FIELD',3)) THEN
              CALL UPFIEL(IBT,IDO,INP,NBH,NBHF,
     '          NBJ,NBJF,NEELEM,NEL,NENP,NFF,NFFACE,
     '          NFLIST,NKB,NKEF,NKH,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NLLIST,
     '          NNF,NNL,NPF,NPL,NPLIST1,NPNE,NPNF,NPNODE,NPNY,NRE,
     '          NRLIST,NVHE,NVHP,NVJE,NVJF,NVJP,NWP,NXI,NXLIST,
     '          NYNO,NYNE,NYNP,CP,DF,REALWORK(OS_PAOPTI),PG,RG,
     '          SE,SF,WG,XA,XAB,XE,XG,XP,YP,ZD,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'FLUX',2)) THEN
              CALL UPFLUX(FIX,NPNODE,NRLIST,NYNP,
     '          STRING,XP,YP,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GAUSS',2)) THEN
              CALL UPGAUS(IBT,IDO,ILPIN,INP,NAN,
     '          NBH,NBHF,NBJ,NBJF,NEELEM,NELIST,NELIST2,NFFACE,NGLIST,
     '          NHE,NHP,NKB,NKH,NKHE,NKJE,NNF,NPF,NPNE,NPNODE,NQET,
     '          NQLIST,%VAL(NQNE_PTR),%VAL(NQS_PTR),NQSCNB,NQXI,NRE,
     '          NRLIST,NSB,NVHE,NVHP,NVJE,NW,NXLIST,NYNE,NYNP,AQ,CE,
     '          CG,CGE,CP,CURVCORRECT,FEXT,PG,PGNQE,RG,SE,WG,XA,XE,
     '          XG,XIG,XP,YG,YGF,YP,YQ,YQS,ZA,ZE,ZG,ZP,STRING,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GEOMETRY',2)) THEN
              CALL UPGEOM(CP,NKH,NKJ,NPLIST1,
     '          NPNODE,NPNY,NRLIST,NVHP,NVJP,
     '          NXLIST,NYNP,NYNR,XP,YP,ZD,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GRID',3)) THEN
C news VJ 21Jan2004 Added NELIST to the UPGRID parameter list
C news VJ 29Jan2004 Added ICQS_SPATIAL to UPGRID param list
C  VJ 9Feb2004 Added NHP,NKH,NPNODE,NVHP,NYNE,YP to UPGRID param list
C                  to be able to call YPZP inside UPGRID
              CALL UPGRID(IBT,IDO,INP,ICQS_SPATIAL,IRCQS_SPATIAL,NAN,
     '          NAQ,NBH,NBJ,NEELEM,NELIST,NENQ,NGAP,NHE,NHP,NKH,
     '          NKHE,NKJE,NLL,NLQ,NPF,NPL,NPNE,NPNODE,NQGP,
     '          NQLIST,%VAL(NQNE_PTR),NQXI,%VAL(NQS_PTR),
     '          NRLIST,NVHE,NVHP,NVJE,NW,NWQ,NXLIST,NXQ,NYNE,NYNP,AQ,
     '          CE,CP,CQ,CURVCORRECT,DL,DNUDXQ,DXDXIQ,DXDXIQ2,
     '          FEXT,GCHQ,GUQ,PG,PROPQ,RCQS_SPATIAL,SE,XA,XE,XG,XIQ,
     '          XP,XQ,YG,YP,YQ,YQS,ZA,ZE,ZG,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'GROWTH',3)) THEN
              CALL UPGROW(NBJ,NEELEM,NRLIST,NXLIST,
     '          XP,YG,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INITIAL',3)) THEN
              CALL UPINIT(IBT,IDO,INP,LD_NP,NBH,NBJ,NEELEM,NELIST,NENP,
     &          NENQ,NFFACE,NFLIST,NHE,NHP,NKH,NKHE,NKJ,NKJE,NNF,NONY,
     &          NPF,NPLIST1,NPLIST3,NPNE,NPNODE,NP_INTERFACE,NQNP,
     &          %VAL(NQS_PTR),NQXI,NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXLIST,
     &          NXQ,NYNE,NYNP,AQ,CQ,CURVCORRECT,DF,DNUDXQ,DXDXIQ,
     &          DXDXIQ2,SE,XA,XAB,XE,XP,XQ,YP,YQ,ZA,ZCROSSING,ZD,ZE,ZP,
     &          STRING,FIX,FIXQ,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INITIAL_VORONOI_VOLUME',1)) THEN
              CALL INITVC(NPNODE,NRLIST,NVCNODE,VC,VC_init,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'INVERSE',3)) THEN
              CALL UPINVE(ISIZE_PHI,ISIZE_TBH(1),ISIZE_TBH(2),NHP,NKH,
     '          NPLIST3,NVHP,NXLIST,NYNP,NYNR,LAPL,LAPLSQR,PHI,T_BH,YP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'LINES',1)) THEN
              CALL UPLINE(NBJ,NEL,NLL,NPL,NPNE,NSB,SE,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MASS',1)) THEN
              CALL UPMASS(NEELEM,NRLIST,NXLIST,BBM,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MATERIAL',2)) THEN
C PM 26-JUL-01 : more arguments added
              CALL UPMATE(IBT,IDO,INP,IRCQS_SPATIAL,NBJ,NEELEM,NELIST,
     &          NENQ,NKJ,NPLIST1,NPNODE,NRLIST,NKH,NKJE,NPF,NPNE,
     '          NQET,%VAL(NQNE_PTR),%VAL(NQS_PTR),NVHP,NVJE,NVJP,
     '          NXI,NXLIST,NYNP,CE,CELL_RCQS_VALUE,CP,CQ,RCQS_SPATIAL,
     '          REALWORK(OS_PAOPTI),SE,XA,XE,
     '          XP,XQ,YP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MESH',2)) THEN
              CALL UPMESH(IBT,IDO,INP,ISEG,%VAL(ISELNO_PTR),
     &          %VAL(ISFIBR_PTR),%VAL(ISFIEL_PTR),%VAL(ISLINE_PTR),
     &          %VAL(ISLINO_PTR),%VAL(ISNONO_PTR),MXI,NAN,NBH,NBJ,
     &          NEELEM,NELIST,NEP,NENP,NGAP,NHE,NHP,NKH,NKHE,NKJE,
     &          %VAL(NLATNE_PTR),NLL,NLLIST,NORD,NPF,NPL,NPNE,NPNODE,
     &          %VAL(NQNE_PTR),%VAL(NQNLAT_PTR),%VAL(NQS_PTR),NQXI,NRE,
     &          NVHE,NVHP,NVJE,NVJP,NW,NXI,NYNE,NYNP,CE,CURVCORRECT,DL,
     &          SE,XA,XAB,XE,XG,XIP,XP,XQ,YG,YP,YQ,ZA,ZE,ZP,CSEG,STRING,
     &          FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'MFI',3)) THEN
C LKC 14-JAN-2011 - Remove for now.              
C              CALL UPMFI(ERROR,*9999)
C AJS 11/2010 Add update motion option
            ELSE IF(ABBREV(CO(3),'MOTION',3)) THEN
              CALL UPMOTI(NXLIST,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'NODES',3)) THEN
              CALL UPNODE(IBT,IDO,INP,NBH,NBJ,NBJF,NEELEM,NEL,NELIST,
     '          NENP,NFF,NFFACE,NHE,NHP,NKJE,NKEF,NKH,NKJ,NLF,NLL,
     '          NLLINE,NNF,NNL,NPF,NP_INTERFACE,NPL,NPLIST1,NPNE,
     '          NPNF,NPNODE,NQNP,%VAL(NQS_PTR),NQXI,NRE,NRLIST,
     '          NVHP,NVJE,NVJF,NVJL,NVJP,NWP,NWQ,NXLIST,NXQ,NYNE,
     '          NYNP,CE,CP,CURVCORRECT,DF,DL,DLL,REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PBOPTI),PG,RG,SE,SF,SP,WG,XA,XE,XG,
     '          XP,XQ,YP,YQ,ZA,ZP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'OPTIMISATION',3)) THEN
              CALL UPOPTI(CE,CELL_RCQS_VALUE,CP,DIPOLE_CEN_NTIME,
     &          DIPOLE_DIR_NTIME,
     '          ILPIN,ISIZE_PHI,NDIPOLES,NEELEM,NMNO,NONY,
     '          NP1OPT,NP2OPT,NP3OPT,NP_INTERFACE,NPNY,NRLIST,
     '          NXLIST,NYNO,NYNP,NYNR,PAOPTY,
     '          CQ,DIPOLE_CEN,DIPOLE_DIR,REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PBOPTI),XP,YP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'ORDERING',3)) THEN
              CALL UPORDR(NBJ,NEELEM,NELIST,NENP,NKJ,NKJE,NORD,NPLIST1,
     &          NPLIST2,NPNE,NPNODE,NRE,NRLIST,NVJE,NVJP,NXI,SE,XP,
     &          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PHI',3)) THEN
              CALL UPPHI(IBT,IDO,INP,ISIZE_PHI,ISIZE_PHIH,
     '          LIST,LD,NBH,NBJ,NDDATA,
     '          NEELEM,NENP,NPF,NHE,NHP,NKH,NKHE,NPNE,NPNODE,NRE,
     '          NRLIST,NVHE,NVHP,NW,NYNE,NYNP,
     '          CURVCORRECT,PHI,PHI_H,SE,XID,YP,ZA,ZE,ZP,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'PRESSURE',3)) THEN
              CALL UPPRES(IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,NELIST,
     '          NGAP,NHE,NHP,NKH,NKHE,NKJE,NMNO,NPF,NPNE,NPNODE,NPNY,
     '          NRE,NRLIST,NVHE,NVHP,NVJE,NW,NXI,NYNE,NYNP,NYNR,
     '          CE,CP,CURVCORRECT,FEXT,REALWORK(OS_PAOPTI),
     '          PF,PG,SE,WG,XA,XE,XG,XP,
     '          YG,YP,ZA,ZE,ZG,ZP,STRING,FIX,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'RESIDUALS',1)) THEN
              CALL UPRESI(IBT,IDO,INP,LGE,NAN,NBH,NBHF,NBJ,NBJF,NEELEM,
     '          NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,NKH,NKHE,NKJE,NNB,
     '          NNF,NPF,NPNE,NPNODE,NPNY,NRE,NSB,NVHE,NVHP,NVJE,NW,
     '          NXI,NXLIST,NYNE,NYNP,NYNR,Z_CONT_LIST,
     '          CE,CG,CGE,CP,CURVCORRECT,FEXT,PG,RE1,RG,SE,
     '          WG,XA,
     '          XG,XP,YG,YGF,YP,ZA,ZA1,Z_CONT,ZE,ZD,ZP,ZP1,STRING,FIX,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SCALE_FACTORS',2)) THEN
              CALL UPSCAL(IBT,IDO,INP,NBJ,NEELEM,NEL,NENP,
     '          NKJE,NKJ,NLL,NLLINE,NNL,NPL,NPNE,NPNODE,
     '          NRE,NRLIST,NSB,NVJE,NVJL,NVJP,DL,DLL,SE,XP,STRING,ERROR,
     '          *9999)
            ELSE IF(ABBREV(CO(3),'SIGNAL',2)) THEN
              CALL UPSIGN(LD,NBJ,NDDATA,NEELEM,WD,XID,ZD,
     '          STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOBOLEV',2)) THEN
              CALL UPSOBE(IBT,IDO,INP,ISIZE_MFI,ISIZE_PHI,ISIZE_TBH,
     '          LD,LDR,LGE,NAN,NBH,NBHF,NBJ,NBJF,NDDL,NDLT,NEELEM,
     '          NEL,NELIST,NFF,NFFACE,NGAP,NHE,NHP,NKB,NKEF,
     '          NKH,NKHE,NKJE,NLL,NLNO,NMNO,
     '          NNB,NNF,NNL,NONL,NONM,NONY,NP_INTERFACE,NP1OPT,
     '          NPF,NPL,NPLIST3,NPNE,NPNODE,NPNY,NRE,NRLIST,NSB,
     '          NVHE,NVHP,NVJE,NVJL,NW,NXI,NXLIST,
     '          NYNE,NYNO,NYNP,NYNR,PAOPTY,Z_CONT_LIST,AQ,
     '          CE,CELL_RCQS_VALUE,CG,CGE,CONY,CP,
     '          CURVCORRECT,DL,
     &          D_RE,D_RI3,D_RP,D_TG,D_ZG,ES,FEXT,
     '          REALWORK(OS_FGRADM),LAPL,LAPLSQR,
     '          MFI,REALWORK(OS_PAOPTI),
     '          REALWORK(OS_PBOPTI),PG,PHI,PHI_H,REALWORK(OS_PMIN),
     '          REALWORK(OS_PMAX),RE1,RE2,REALWORK(OS_RESID),
     '          REALWORK(OS_RESIDM),REALWORK(OS_RESJAC),RG,
     '          SE,T_BH,WG,WK1_INV,WU,XA,XE,XG,XID,XIG,
     '          XN,XP,YG,YGF,YP,
     '          ZA,ZA1,Z_CONT,ZD,ZE,ZG,ZG1,ZP,ZP1,STRING,FIX,
     '          ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'SOLUTION',3)) THEN
              CALL UPSOLU(CP,IBT,IDO,INP,NAN,NBH,NBJ,NEELEM,
     '          NEP,NHE,NHP,NHQ,NKHE,NKH,NKJ,NONY,NPF,NPLIST1,NPLIST3,
     '          NPLIST4,NPNE,NPNODE,NPNY,NQNY,NRLIST,NRLIST2,NVHE,NVHP,
     '          NVJE,NVJP,NW,NXLIST,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '          CURVCORRECT,PHI,PHI_H,SE,XAB,XIP,XP,YG,YP,YQ,YQS,
     '          ZA,ZD,ZE,ZG,ZP,STRING,ERROR,*9999)
C SMAR009 18/01/99 removed NP_INTERFACE,
            ELSE IF(ABBREV(CO(3),'SOURCE',3)) THEN
              CALL UPSOUR(DIPOLE_CEN_NTIME,DIPOLE_DIR_NTIME,NBH,
     '          NDIPOLES,NENP,NKH,NP_INTERFACE,NPNODE,
     '          NRLIST,NW,NXLIST,NYNP,CE,DIPOLE_CEN,DIPOLE_DIR,
     '          GD,XG,XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VARIABLE',3)) THEN
              CALL UPVARI(XP,STRING,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'VORONOI',7)) THEN
              CALL UPVORO(NENFVC,NFVC,NODENVC,NPNODE,NRLIST,
     '          NVCNODE,VC,VC_INIT,XNFV,XP,ZA,STRING,ERROR,*9999)
c            ELSE IF(ABBREV(CO(3),'VORO_DEFORM',7)) THEN
c              CALL CALC_VORO_LUNG(IBT,IDO,INP,NBH,NBJ,NEELEM,NENFVC,
c     '          NENP,NEP,NEPZ,NFVC,NHE,NKHE,NKJE,NODENVC,NPF,
c     '          NPLIST1,NPNE,NPNODE,NVCBBM,NVCNODE,NVHE,NVJE,NW,NXI,BBM,
c     '          CURVCORRECT,SE,VC,VC_init,XA,XE,XIP,XNFV,XP,ZA,ZAV,ZD,
c     '          ZE,ZP,ERROR,*9999)
C              CALL UPVOROFROMDEFORMED(NBH,NHE,NKE,NPF,NPNE,1,NXLIST(0),
C     '          NVHE,NW,CURVCORRECT,SE,ZA,ZP,IBT,IDO,INP,NPNODE,
C     '          XIP,XP,NEP,ERROR,*9999)
c            ELSE IF(ABBREV(CO(3),'VORO_PV',7)) THEN
c              CALL UPPRESSVOROPVCURVE(NPNODE,NRLIST,NVCNODE,VC,VC_init,
c     '          XP,ERROR,*9999)
c           ELSE IF(ABBREV(CO(3),'VORO_IDEAL',7)) THEN
c              CALL UPPRESSVOROIDEAL(IBT,IDO,INP,NBH,NBJ,NEELEM,NEPZ,
c     '         NFVC,NHE,NKHE,NODENVC,NPF,NPNE,NPNODE,NRLIST,
c     '         NVCNODE,NVHE,NW,CURVCORRECT,SE,VC,VC_INIT,XP,ZA,ZD,
c     '         ZE,ZP,ERROR,*9999)
            ELSE IF(ABBREV(CO(3),'XI',1)) THEN
C MHT 09/11/10 
              CALL UPXI(IBT,IDO,INP,%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,FD,LD,LN,MXI,
     &          NBJ,NBJF,NBH,NBHF,NDDL,NDLT,NDP,NEELEM,NELIST,NFF,
     &          NFFACE,NFLIST,NHE,NKEF,NKHE,NKJE,NNF,NPF,NPLIST1,NPNE,
     &          NPNF,NRE,NRLIST,NVHE,NVJE,NVJF,NVJP,NW,NXI,CE,CG,CGE,CP,
     &          CURVCORRECT,PG,SE,SF,SQ,WD,WDL,XA,XE,XG,XID,XIDL,XP,ZA,
     &          ZD,ZDD,ZDL,ZE,ZP,CSEG,STRING,ERROR,*9999)
c              CALL UPXI(IBT,IDO,INP,%VAL(ISDANO_PTR),%VAL(ISDAPR_PTR),
c     &          %VAL(ISDATA_PTR),%VAL(ISDATR_PTR),ISEG,LD,LN,MXI,NBJ,
c     &          NBJF,NBH,NDDL,NDLT,NDP,NEELEM,NELIST,NKEF,NKHE,NKJE,NNF,
c     &          NPF,NPNE,NPNF,NRE,NVHE,NVJE,NVJF,NW,NXI,CE,CG,CGE,CP,
c     &          CURVCORRECT,PG,SE,SF,SQ,WD,WDL,XA,XE,XG,XID,XIDL,XP,ZA,
c     &          ZD,ZDD,ZDL,ZE,ZP,CSEG,STRING,ERROR,*9999)
C SMAR009 19/01/99 removed NPL,
            ELSE
              IF(COMMAND) CALL STAND(END,STRING,ERROR,*9999)
            ENDIF
            IF(DOCUMENT.AND.noch3.LT.NTCH3) THEN
              CONTINUE3=.TRUE.
            ELSE
              noch3=0
              CONTINUE3=.FALSE.
            ENDIF !document
          ENDDO !continue3

        ELSE IF(ABBREV(CO(2),'WRITE',1)) THEN
          noco=2
C*** 10/10/08 JHC added Z_CONT in order to pass it into IODATA.f
          CALL WRITEF(IBT,IDO,INP,ISC_GKK,
     '      ISIZE_MFI,ISIZE_PHI,ISIZE_PHIH,ISIZE_TBH,
     '      ISR_GKK,ITHRES,LD_NP,MAP_ART_VEIN,NAN,NBH,NBJ,NBJF,
     '      NDDL,NDET,NDLT,NDP,NEELEM,NEL,NELIST,NFF,NGAP,NHE,NHP,
     '      NHQ,NKEF,NKH,NKHE,NKJ,NKJE,NLF,NLL,NNF,NNL,
     '      NONY,NPF,NPL,NPNE,NPNODE,NPNY,NP_INTERFACE,
     '      NQNY,NRE,NRLIST,NRLIST2,NVHE,NVHP,NVJE,NVJP,NW,NXI,
     '      NXLIST,NYNE,NYNO,NYNP,NYNR,NYQNR,
     '      CE,CG,CONY,CP,CURVCORRECT,DET,DL,FEXT,GKK,GR,GRR,MFI,
     '      PE,PF,PG,PHI,PHI_H,PHI_H_EXACT,SE,T_BH,T_BH_INV,THRES,
     '      WD,WDL,WG,XA,XE,XG,XID,XIDL,XIG,XP,YP,YQ,YQS,ZA,Z_CONT,
     '      ZCROSSING,ZD,ZD2,ZDL,ZE,
     '      ZG,ZP,STRING,FIX,FIXP,ERROR,*9999)

C rgb 14/09/1999 All archived into archive_fe01.f
C        ELSE IF(ABBREV(CO(2),'ZOOM',1)) THEN
C          IF(COMMAND) THEN
C            noco=2
Cc MS Archived wouldn't comile on linux
Cc            CALL ZOOM(STRING,ERROR,*9999)
C          ELSE
Cc           CALL DOCUM('fe20','doc','ZOOM',ERROR,*9999)
C          ENDIF

        ELSE
          IF(COMMAND) THEN
            CALL STAND(END,STRING,ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'ASSIGN',1)) THEN
c            CALL DOCUM('fe20','doc','ASSIGN',ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'DEASSIGN',2)) THEN
c            CALL DOCUM('fe20','doc','DEASSIGN',ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'DISPLAY',3)) THEN
c            CALL DOCUM('fe20','doc','DISPLAY',ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'LEARN',3)) THEN
c            CALL DOCUM('fe20','doc','LEARN',ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'QUIT',1)) THEN
c            CALL DOCUM('fe20','doc','QUIT',ERROR,*9999)
c          ELSE IF(ABBREV(CO(2),'TRACE',1)) THEN
c            CALL DOCUM('fe20','doc','TRACE',ERROR,*9999)
          ENDIF
        ENDIF

        IF(DOCUMENT.AND.noch2.LT.NTCH2) THEN
          CONTINUE2=.TRUE.
        ELSE
          CONTINUE2=.FALSE.
        ENDIF !document

      ENDDO !continue2
      IF(CMGUI_LINK) THEN
C cpb 4/3/97 temporary disable for class version
C db 7/10/97 temporary disable disabled
C cs 16/11/99 temporary disable disable disabled
C       Check for any information updates to send to CMGUI
C        CALL CMGUI_LINK_UPDATE(LD,NKH,NKJ,NPNODE,NVHP,NVJP,
C     '    NP_INTERFACE,NYNP,WD,XID,XP,YP,ZD,FIX,ERROR,*9999)
      ENDIF

      CALL EXITS('FEM_DYNAM')
      RETURN
 9999 CALL ERRORS('FEM_DYNAM',ERROR)
      CALL EXITS('FEM_DYNAM')
      RETURN 1
      END
