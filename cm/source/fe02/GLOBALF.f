      SUBROUTINE GLOBALF(IBT,IDO,INP,NBH,NBJ,
     '  NENP,njj,NKB,NKH,NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NPF,NPL,NPNE,
     '  NPNODE,NPNY,nr,NVHE,NVHP,NVJE,NWP,
     '  nx,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,
     '  FIX,ERROR,*)

C#### Subroutine: GLOBALF
C###  Description:
C###    GLOBALF calculates the mapping arrays NYNO/NONY/CYNO/CONY for
C###    fit variable njj.
C###  See-Also: GLOBALH

      IMPLICIT NONE
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),
     '  IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NBH(NHM,NCM,NEM),NBJ(NJM,NEM),NENP(NPM,0:NEPM,0:NRM),njj,
     '  NKB(2,2,2,NNM,NBFM),NKH(NHM,NPM,NCM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),nr,
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM),
     '  NVJE(NNM,NBFM,NJM,NEM),NWP(NPM,2),nx,
     '  NXI(-NIM:NIM,0:NEIM,0:NEM),NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYNY(0:NYYM,NYM,NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),
     '  CYNO(0:NYOM,NOOPM,NRCM,0:NRM),CYNY(0:NYYM,NYM,NRM),
     '  SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER connected_ne_nn,INCREASE_NOM,
     '  INCREASE_NYOM,LD,MAPPED_TOT,
     '  nae,nb,ne,nffe,nhj,nhx,nh,ni,NITB,nl,nv,nk,nk1,nn,nne,nne1,
     '  no,nonode,noy,NOYT,no_nynr,no_tot(2),np,np1,nrr,nrc,ns,nu,
     '  NUNK(8),ny,ny2,nyr,nyr2,nyy(2),nyy2(2),nyo,nynp1
      REAL*8 COY,RATIO,PSI,weight,XD(3),XI(3)
      CHARACTER CHAR*1
      LOGICAL DONE
!     Functions
      INTEGER GETNYR,IDIGITS

      DATA NUNK/1,2,4,6,7,9,10,11/

      CALL ENTERS('GLOBALF',*9999)

C***  Initialise mapping arrays above current region

      INCREASE_NOM=0
      INCREASE_NYOM=0

      DO nrr=nr,NRT
        DO nrc=1,2
          DO no_nynr=1,NYNR(0,0,1,nrr)
            ny=NYNR(no_nynr,2-nrc,1,nrr)
            NONY(0,ny,nrc,nr)=0
            NONY(1,ny,nrc,nr)=0
            CONY(0,ny,nrc,nr)=0.0d0
            CONY(1,ny,nrc,nr)=0.0d0
          ENDDO
          DO no=1,NOM
            DO nyo=1,NYOM
              NYNO(nyo,no,nrc,nrr)=0
              CYNO(nyo,no,nrc,nrr)=0.0d0
            ENDDO
            NYNO(0,no,nrc,nrr)=0
            CYNO(0,no,nrc,nrr)=0.0d0
          ENDDO
          NOT(nrc,1,nr,nx)=0
        ENDDO !nrc
      ENDDO

C*** Calculate mapping arrays

      DO nrc=1,2
        no_tot(nrc)=0
      ENDDO !nrc
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO nhj=1,NUM_FIT(njj)
          nhx=NLH_FIT(nhj,3,njj)
          nh=NH_LOC(nhx,nx)
          DO nv=1,NVHP(nh,np,1)
            DO nk=1,NKH(nh,np,1)
              ny=NYNP(nk,nv,nh,np,0,1,nr)
              IF(FIX(ny,1)) THEN !delete row and column from the system
C               Do Nothing
              ELSEIF(NWP(np,1).GT.0) THEN
C               Do Nothing, set up hanging nodes below
              ELSE !variable needs to be solved for
                DONE=.FALSE.
                nyy(1)=NYNP(nk,nv,nh,np,1,1,nr) !global row #
                nyy(2)=ny !global variable #
                IF(JTYP2A.EQ.1) THEN !treat versions as coincident
C MPN 3Jul2003: dropping down of indices was all wrong here!!!
                  CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),
     '              NKB,NKHE,NNB,
     '              NONY(0,1,1,nr),NPNE,NPNY,nr,NVHE,
     '              NVHP,NXI,ny,ny2,NYNP(1,1,1,1,0,1,nr),
     '              NYNY(0,1,nr),CYNY(0,1,nr),RATIO,FIX,*9999)
C MPN 3Jul2003 OLD
C                  CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),
C     '              NKB,NKHE,NNB,
C     '              NONY(0,1,1,nr),NPNE,NPNY(0,1,nr),nr,NVHE,
C     '              NVHP(1,1,1),NXI,ny,ny2,NYNP,NYNY(0,1,nr),
C     '              CYNY(0,1,nr),RATIO,FIX,*9999)
                  IF(ny2.NE.ny) THEN
C                   There is an equivalent mesh degree of freedom for
C                   which an no as already been assigned.
                    DONE=.TRUE.
                    IF(ny2.EQ.0) THEN !dof not used
                      FIX(ny,1)=.TRUE.
                    ELSE IF(FIX(ny2,1)) THEN
                      FIX(ny,1)=.TRUE.
                    ELSE
C                     equivalent mesh degree of freedom
                      nyy2(1)=GETNYR(1,NPNY,nr,1,0,ny2,NYNE,NYNP) !row#
                      nyy2(2)=ny2 !global col#
                      DO nrc=1,2 !nrc=1,2 local row and local column
                        nyr=nyy(nrc)
                        nyr2=nyy2(nrc)
                        NOYT=NONY(0,nyr2,nrc,nr)
                        NONY(0,nyr,nrc,nr)=NOYT !CYNO(0) is already 0
                        DO noy=1,NOYT
                          no=NONY(noy,nyr2,nrc,nr)
                          NONY(noy,nyr,nrc,nr)=no
                          COY=RATIO*CONY(noy,nyr2,nrc,nr)
                          CONY(noy,nyr,nrc,nr)=COY
                          IF(no.LT.NOM) THEN
                            nyo=NYNO(0,no,nrc,nr)+1
                            NYNO(0,no,nrc,nr)=nyo
                            IF(nyo.LE.NYOM) THEN
                              NYNO(nyo,no,nrc,nr)=nyr
C                             assuming COY is 1 or -1
                              CYNO(nyo,no,nrc,nr)=COY
                            ELSE IF(nyo.GT.INCREASE_NYOM) THEN
                              INCREASE_NYOM=nyo
                            ENDIF
                          ENDIF
                        ENDDO ! noy
                      ENDDO !nrc
                    ENDIF !ny2=0/FIX
                  ENDIF !ny.NE.ny1
                ENDIF !JTYP2
                IF(.NOT.DONE) THEN
                  DO nrc=1,2 !rows and columns
                    no_tot(nrc)=no_tot(nrc)+1
                    NONY(0,nyy(nrc),nrc,nr)=1
                    NONY(1,nyy(nrc),nrc,nr)=no_tot(nrc)
                    CONY(0,nyy(nrc),nrc,nr)=0.0d0
                    CONY(1,nyy(nrc),nrc,nr)=1.0d0
                    IF(no_tot(nrc).LE.NOM) THEN
                      NYNO(0,no_tot(nrc),nrc,nr)=1
                      NYNO(1,no_tot(nrc),nrc,nr)=nyy(nrc)
                      CYNO(0,no_tot(nrc),nrc,nr)=0.0d0
                      CYNO(1,no_tot(nrc),nrc,nr)=1.0d0
                    ENDIF
                  ENDDO !nrc
                ENDIF !not done
              ENDIF !fix
            ENDDO !nk
          ENDDO !nv
        ENDDO !nhj
      ENDDO !nonode

C new CS 1/8/2001
      IF(JTYP2C.EQ.1) THEN !hanging nodes
C       Set up hanging nodes
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhj=1,NUM_FIT(njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1)
              DO nk=1,NKH(nh,np,1)
                ny=NYNP(nk,nv,nh,np,0,1,nr)
                IF(NWP(np,1).GT.0) THEN ! Hanging node
                  nyy(1)=NYNP(nk,nv,nh,np,1,1,nr) !global row #
                  nyy(2)=ny !global variable #

                  ne=NWP(np,1) ! element ny hangs in
                  nb=NBH(nh,1,ne)
                  NITB=NIT(nb)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '              XA(1,1,ne),XE,XP,ERROR,*9999)
                  LD=0
                  XD(1)=XP(1,1,1,np)
                  XD(2)=XP(1,1,2,np)
                  XD(3)=XP(1,1,3,np)
                  DO ni=1,NITB
                    XI(ni)=0.5d0
                  ENDDO
                  CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,
     '              NITB,nr,0.d0,XE,XI,XI,XD,.FALSE.,ERROR,*9999)
                  IF(NITB.EQ.2) THEN
                    nae=NWP(np,2)
                    nl=NLL(nae,ne) ! global line np hangs on
                    MAPPED_TOT=NNL(0,nae,nb)
                  ELSE
                    nffe=NWP(np,2)
                    MAPPED_TOT=NNF(0,nffe,nb)
                  ENDIF

C                 Set row (nrc=1) and column (nrc=2) mappings
                  DO nrc=1,2 !nrc=1,2 local row and local column
                    NONY(0,nyy(nrc),nrc,nr)=0
                    DO nn=1,MAPPED_TOT !number of nodes ny is mapped to
                      IF(NITB.EQ.2) THEN
                        np1=NPL(1+nn,1,nl)
                      ELSE
                        connected_ne_nn=NNF(nn+1,nffe,nb)
                        np1=NPNE(connected_ne_nn,nb,ne)
                      ENDIF
                      nne=1
                      DO WHILE ((np1.NE.NPNE(nne,nb,ne))
     '                  .AND.(nne.LE.NNT(nb)))
                        nne=nne+1
                      ENDDO
                      DO nk1=1,NKT(nne,nb)
                        nynp1=NYNP(nk1,nv,nh,np1,nrc,1,nr)
                        DO noy=1,NONY(0,nynp1,nrc,nr)
                            no=NONY(noy,nynp1,nrc,nr)
                            nu=NUNK(nk)
                            weight=PSI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                        INP(1,1,nb),nb,nne,nk1,nu,XI)
                            ns=0
                            DO nne1=1,nne-1
                              ns=ns+NKT(nne1,nb)
                            ENDDO
                            ns=ns+nk1
                            weight=weight*SE(ns,nb,ne)

                            ! Devide through by SE since we want
                            ! d/ds not d/dxi
                            ! The trick is to get the correct
                            ! scale factor. The scale factors required
                            ! at the hanging nodes are calculated in
                            ! HANGING_NODE_DETECT and stored in SP
                            ns=0
                            DO nne1=1,nne-1
                              ns=ns+NKT(nne1,nb)
                            ENDDO
                            ns=ns+nk
C!!!
C                         This is a cheap hack because I haven't bothered
C                         to set up SP correctly. Should fill SP for all nj,
C                         currently only doing for nj=1 therefore not correct
C                         for the pressure basis. This corrects it.
                            IF(DABS(SP(nk,nb,np)).LE.LOOSE_TOL) THEN
                              SP(nk,nb,np)=1.0d0
                            ENDIF
                            weight=weight/SP(nk,nb,np)

                            IF(DABS(weight).GT.LOOSE_TOL) THEN
                              NONY(0,nyy(nrc),nrc,nr)=
     '                        NONY(0,nyy(nrc),nrc,nr)+1
                              CONY(0,nyy(nrc),nrc,nr)= 0.0d0

                           NONY(NONY(0,nyy(nrc),nrc,nr),nyy(nrc),nrc,nr)
     '                         =no
                           CONY(NONY(0,nyy(nrc),nrc,nr),nyy(nrc),nrc,nr)
     '                         =weight

                              NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                              CYNO(0,no,nrc,nr)=0.0d0

                              NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                              CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=weight
                            ENDIF
                        ENDDO ! noy
                      ENDDO ! nk
                    ENDDO ! nn
                  ENDDO ! nrc
                ENDIF ! hanging
              ENDDO ! nk
            ENDDO ! nv
          ENDDO ! nhj
        ENDDO ! nonode
      ENDIF ! JTYP2

      DO nrc=1,2
        IF(no_tot(nrc).LE.NOM) THEN
          NOT(nrc,1,nr,nx)=no_tot(nrc)
        ELSE IF(no_tot(nrc).GT.INCREASE_NOM) THEN
          INCREASE_NOM=no_tot(nrc)
        ENDIF
C          CALL ASSERT(NOT(nrc,1,nr,nx).LE.NOM,'>> Increase NOM',
C     '      ERROR,*9999)
      ENDDO !nrc
      IF(INCREASE_NOM.NE.0) THEN
        WRITE(CHAR,'(I1)') IDIGITS(INCREASE_NOM)
        WRITE(ERROR,'(''>>Increase NOM to '',I'//CHAR//')')
     '    INCREASE_NOM
        GO TO 9999
      ELSE IF(INCREASE_NYOM.NE.0) THEN
        WRITE(CHAR,'(I1)') IDIGITS(INCREASE_NYOM)
        WRITE(ERROR,'(''>>Increase NYOM to '',I'//CHAR//')')
     '    INCREASE_NYOM
        GO TO 9999
      ENDIF

      CALL EXITS('GLOBALF')
      RETURN
 9999 CALL ERRORS('GLOBALF',ERROR)
      CALL EXITS('GLOBALF')
      RETURN 1
      END


