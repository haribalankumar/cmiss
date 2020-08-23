      SUBROUTINE XPEQ30(ISC_GQ,ISR_GQ,LGE,NBH,NBJ,ne,NEELEM,NHP,NHST,
     '  NKJE,NPF,NP_INTERFACE,NPNE,nr,NRE,NVHE,NVJE,nx,
     '  NYNP,CURVCORRECT,ES,EM,GQ,PG,RG,SE,WG,XA,XE,XG1,XP,ERROR,*)

C#### Subroutine: XPEQ30
C###  Description:
C###    XPEQ30 calculates element flux matrices (ES(nhs1,nhs2)) for
C###    linear static equations.
C****  AJP 21/2/95 Currently only implemented for Laplaces equation.

C cpb 23/1/97 adding curvature corrections

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER ISC_GQ(NISC_GQM),ISR_GQ(NISR_GQM),LGE(NHM*NSM,NRCM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),ne,NEELEM(0:NE_R_M,0:NRM),
     '  NHP(NPM),NHST(*),NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  nr,NRE(NEM),NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CURVCORRECT(2,2,NNM,NEM),ES(NHM*NSM,NHM*NSM),
     '  EM(NHM*NSM,NHM*NSM),GQ(NZ_GQ_M),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),WG(NGM,NBM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG1(NJM,NUM,NGM),
     '  XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER MATCH,nb1j,nbbe,nbbe2,nbfe,nebe,ng,nh,nhx,ni,nj,nk1,
     '  nk2,nn,nn2,nnbe,nnbe2,nnbe3,nnfe,noelem,
     '  nonr,np,npbe,npfe,nr1,nrbe,ns1,nsbe1,ns2,nsbe2,nu,NU1(0:3),
     '  nv,ny1,ny2,nz
      REAL*8 DDOT,DXXI(3,3),RWG,SUM
      LOGICAL FOUND_NE,FOUND_NR,ONINTERFACE

      EXTERNAL DDOT

      DATA NU1/1,2,4,7/

      CALL ENTERS('XPEQ30',*9999)

C***  Construct equivalent BE matrix in place of FE load vector

      nbfe=NBH(NH_LOC(1,nx),2,ne)
      DO nn=1,NNT(nbfe)
        np=NPNE(nn,nbfe,ne)
C       Check if np is shared by a BE region.
        FOUND_NR=.FALSE.
        IF((NP_INTERFACE(np,0).GT.1))THEN
C         Find a BE region sharing node np
          DO nonr=1,NP_INTERFACE(np,0)
            nr1=NP_INTERFACE(np,nonr)
            IF(ITYP4(nr1,nx).EQ.2)THEN
              FOUND_NR=.TRUE.
              nrbe=nr1
            ENDIF
          ENDDO !nr
          IF(FOUND_NR) THEN
C           Node np is shared by BE region nrbe.
C           Find appropriate BE element and local node number.
            FOUND_NE=.FALSE.
            noelem=1
            DO WHILE(.NOT.FOUND_NE.AND.noelem.LE.NEELEM(0,nrbe))
              nebe=NEELEM(noelem,nrbe)
              nbbe=NBH(NH_LOC(1,nx),1,nebe) !family number
              MATCH=0
              DO nn2=1,NNT(nbbe)
                IF(np.EQ.NPNE(nn2,nbbe,nebe))THEN
C                 Node np is in element nebe. Check that the other BE
C                 nodes are also in element ne.
                  nnbe=nn2
                  MATCH=1
                  DO nnbe2=1,NNT(nbbe)
                    IF(nnbe2.NE.nnbe)THEN
                      DO nnfe=1,NNT(nbfe)
                        IF(NPNE(nnbe2,nbbe,nebe).EQ.NPNE(nnfe,nbfe,ne))
     '                    MATCH=MATCH+1
                      ENDDO !nnfe
                    ENDIF
                  ENDDO !nnbe2
                ENDIF !np
              ENDDO !nn2
              IF(MATCH.EQ.NNT(nbbe)) THEN
                FOUND_NE=.TRUE.
              ELSE
                noelem=noelem+1
              ENDIF
            ENDDO !Found BE element nebe
            CALL ASSERT(FOUND_NE,'>>Error in finding BE element',
     '        ERROR,*9999)
            CALL XPXE(NBJ(1,nebe),NKJE(1,1,1,nebe),NPF(1,1),
     '        NPNE(1,1,nebe),NRE(nebe),
     '        NVJE(1,1,1,nebe),SE(1,1,nebe),
     '        XA(1,1,nebe),XE,XP,ERROR,*9999)
            nb1j=NBJ(1,nebe) !Assumes same basis used for each nj
C           Calculate Jacobian etc
            DO ng=1,NGT(nb1j)
              DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                DO nu=1,NUT(nb1j)
C cpb 7/9/95 Using BLAS
C                SUMXG=0.0d0
C                DO ns1=1,NST(nb1j)
C                  SUMXG=SUMXG+PG(ns1,nu,ng,nb1j)*XE(ns1,nj)
C                ENDDO
C                XG1(nj,nu,ng)=SUMXG
                  XG1(nj,nu,ng)=DDOT(NST(nb1j),PG(1,nu,ng,nb1j),1,
     '              XE(1,nj),1)
                ENDDO !nu
              ENDDO !nj
              DO ni=1,NIT(nb1j)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  DXXI(nj,ni)=XG1(nj,NU1(ni),ng)
                ENDDO !nj
              ENDDO !ni
              IF(NIT(nbbe).EQ.1)THEN
                RG(ng)=DSQRT(DXXI(1,1)*DXXI(1,1)+DXXI(2,1)*DXXI(2,1))
              ELSE
                RG(ng)=
     '            DSQRT((DXXI(2,1)*DXXI(3,2)-DXXI(2,2)*DXXI(3,1))**2+
     '            (DXXI(1,2)*DXXI(3,1)-DXXI(1,1)*DXXI(3,2))**2+
     '            (DXXI(1,1)*DXXI(2,2)-DXXI(1,2)*DXXI(2,1))**2)
              ENDIF
            ENDDO !End ng loop
            DO nhx=1,NHP(np)
              nh=nh_loc(nhx,nx)
              ns1=0
              nsbe1=0
              DO nnfe=1,nn-1
                ns1=ns1+NKT(nnfe,nbfe)
              ENDDO
              DO nnbe2=1,nnbe-1
                nsbe1=nsbe1+NKT(nnbe2,nbbe)
              ENDDO
              nbbe2=NBH(nh,2,nebe)
              DO nk1=1,NKT(nnbe,nbbe)
                ns1=ns1+1
                nsbe1=nsbe1+1
                ns2=0
                DO nn2=1,NNT(nbfe)
                  ONINTERFACE=.FALSE.
                  nnbe2=1
                  npfe=NPNE(nn2,nbfe,ne)
                  DO WHILE(.NOT.ONINTERFACE.AND.nnbe2.LE.NNT(nbbe))
                    npbe=NPNE(nnbe2,nbbe,nebe)
                    IF(npbe.EQ.npfe) THEN
                      ONINTERFACE=.TRUE.
                    ELSE
                      nnbe2=nnbe2+1
                    ENDIF
                  ENDDO
                  IF(ONINTERFACE) THEN
                    nsbe2=0
                    DO nnbe3=1,nnbe2-1
                      nsbe2=nsbe2+NKT(nnbe3,nbbe)
                    ENDDO
                    DO nk2=1,NKT(nnbe2,nbbe)
                      ns2=ns2+1
                      nsbe2=nsbe2+1
                      SUM=0.0d0
C                     Evaluate integral over element nebe with basis
C                     functions nbbe and nbbe2 as the integrand.
                      DO ng=1,NGT(nbbe)
                        RWG=RG(ng)*WG(ng,nbbe)
                        SUM=SUM+PG(nsbe1,1,ng,nbbe)*
     '                    PG(nsbe2,1,ng,nbbe2)*RWG
                      ENDDO !ng
                      ES(ns1,ns2)=SUM*SE(nsbe1,nbbe,nebe)*
     '                  SE(nsbe2,nbbe2,nebe)
                      IF(BEMCURVATURECORRECTION) THEN
                        IF((nk1.EQ.2.OR.nk1.EQ.3).AND.(nk2.EQ.2.OR.
     '                    nk2.EQ.3)) THEN
                          EM(ns1,ns2)=-SUM*SE(nsbe1,nbbe,nebe)*
     '                      SE(nsbe2,nbbe2,nebe)*
     '                      CURVCORRECT(nk1-1,nk2-1,nnbe2,nebe)
                        ENDIF
                      ENDIF
                    ENDDO !End of nk2 loop
                  ELSE
                    ns2=ns2+NKT(nn2,nbfe) !Correct ns2 ? AJP 26/1/96
cc                    ns2=ns2+NKT(nnbe2,nbbe)
cC!!!! WARNING       assumes NKT the same for each node in the BE element.
cC!!!!               What about simplex-type elements??
                  ENDIF
                ENDDO !End of nn2 loop
              ENDDO !End of nk1 loop
            ENDDO !End of nh loop
          ENDIF !Found_nr
        ELSE !No be region sharing element face
          DO nhx=1,NHP(np)
            nh=nh_loc(nhx,nx)
            ns1=0
            DO nnfe=1,nn-1
              ns1=ns1+NKT(nnfe,nbfe)
            ENDDO !nnfe
            nv=NVHE(nn,nbfe,nh,ne)
            DO nk1=1,NKT(nn,nbfe)
              ns1=ns1+1
              ny1=LGE(ns1,1) !row #
              ny2=NYNP(nk1,nv,nh,np,2,2,nr) !local column #
              CALL SPARSE(ny1,ny2,NYT(1,2,nx),nz,NZ_GQ_M,
     '          NZT(2,nx),ISC_GQ,ISR_GQ,KTYP24,ERROR,*9999)
              IF(nz.NE.0.AND.ABS(GQ(nz)).LT.RDELTA) THEN
                ns2=1 !is there a better way ?
                DO WHILE(ny2.NE.IABS(LGE(ns2,2)))
                  ns2=ns2+1
                  IF(ns2.GT.NHST(2))THEN
                    ERROR='>>Error in value of ns2'
                    GOTO 9999
                  ENDIF
                ENDDO
C news AJP 22/1/96
                ES(ns1,ns2)=1.0d0
c                IF(DABS(ER(ns1)).LT.RDELTA) THEN
c                  ES(ns1,ns2)=1.0d0
cC                 To avoid a zero column in GQ and to ensure that GQ
cC                 only gets a 1 on its diagonal.
c                ELSE
c                  ES(ns1,ns2)=ER(ns1)
c                ENDIF
C newe AJP
              ENDIF
            ENDDO !nk1
          ENDDO !nh
        ENDIF !End of interface np loop
      ENDDO !End of nn loop

      CALL EXITS('XPEQ30')
      RETURN
 9999 CALL ERRORS('XPEQ30',ERROR)
      CALL EXITS('XPEQ30')
      RETURN 1
      END



