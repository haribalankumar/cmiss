      SUBROUTINE GLOBALH(IBT,IDO,INP,NAN,NBH,NBJ,NELIST,NENP,NHE,NKB,
     '  NKHE,NKJE,NLL,NNB,NNF,NNL,NONY,NP_INTERFACE,NPF,NPL,NPNE,NPNY,
     '  nr,NRE,NVHE,NVHP,NVJE,NWP,nx,NXI,NYNE,NYNO,NYNP,NYNR,NYNY,
     '  NYQNR,CONY,CYNO,CYNY,SE,SP,XA,XE,XP,YP,FIX,ERROR,*)

C#### Subroutine: GLOBALH
C###  Description:
C###    GLOBALH calculates the mapping between the mesh dofs, ny, and
C###    the solution dofs, no. It sets up the arrays NONY, CONY and
C###    their inverse arrays NYNO and CYNO for region nr.

C#### Variable: NONY(0:noy,ny,nrc,0:nr,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALH,IPOPTI
C###  Description:
C###    NONY(0,ny,nrc,0:nr,nx) is the number of solution dofs (noy) at
C###    a mesh dof, ny, for a region nr and problem type nx. For
C###    region 0 the complete multiregion mapping is defined. nrc=1
C###    gives the row (equation number) mappings, and nrc=2 gives the
C###    global column (variable number) mappings. For one to one
C###    mappings noy=1, for one ny to many no mappings noy>1, for many
C###    nys to one no mappings noy=1.  NONY(1..noy,ny,nrc,0:nr,nx) are
C###    the solution dof numbers (no) at a mesh dof ny for a region nr
C###    and problem type nx.

C#### Variable: NOT(nrc,nc,nr,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALH
C###  Description:
C###    NOT(nrc,nc,nr,nx) is the number of solution dofs for matrix
C###    nc in a region nr and problem type nx. nrc=1 gives the number
C###    of rows (equations) and nrc=2 gives the number of global
C###    columns (variables).

C#### Variable: NOQT(nrc,nc,nr,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALH
C###  Description:
C###    NOQT(nrc,nc,nr,nx) is the number of grid solution dofs
C###    for matrix nc in a region nr and problem type nx. nrc=1
C###    gives the number of rows (equations) and nrc=2 gives the
C###    number of global columns (variables).

C#### Variable: NYNO(0:nyo,no,nrc,0:nr,nx)
C###  Type: INTEGER
C###  Set_up: GLOBALH
C###  Description:
C###    NYNO(0,no,nrc,0:nr,nx) is the number of mesh dofs (nyo) at a
C###    solution dof, no, for a region nr and problem type nx. For one
C###    to one mappings nyo=1, for one no to many ny mappings nyo>1.
C###    NYNO(1..nyo,no,nrc,0:nr,nx) are the mesh dof numbers (ny) at
C###    a solution dof, no, for a region nr and problem type nx.

C#### Variable: CONY(0:noy,ny,nrc,0:nr,nx)
C###  Type: REAL*8
C###  Set_up: GLOBALH
C###  Description:
C###    CONY(0,ny,nrc,0:nr,nx) is the constant coupling coefficent
C###    applied to the transformation between ny and no for region nr
C###    and problem type nx. For region 0 the complete multiregion
C###    mapping is defined. nrc=1 gives the row (equation number)
C###    mappings, nrc=2 gives the global column (variable number)
C###    mappings. CONY(1..noy,ny,nrc,0:nr,nx) are the coupling
C###    coefficents for the ny to no mappings for region nr and problem
C###    type nx. ie XO(noy,nx)=YP(ny,nx)*CONY(noy,ny,nrc,nr,nx)
C###    +CONY(0,ny,nrc,nr,nx).

C#### Variable: CYNO(0:nyo,no,nrc,0:nr,nx)
C###  Type: REAL*8
C###  Set_up: GLOBALH
C###  Description:
C###    CYNO(0,no,nrc,0:nr,nx) is the constant coupling coefficent
C###    applied to the transformation between no and ny for region nr
C###    and problem type nx. For region 0 the complete multiregion
C###    mapping is defined. nrc=1 gives the row (equation number)
C###    mappings, nrc=2 gives the column (variable number) mappings.
C###    CYNO(1..nyo,no,nrc,0:nr,nx) are the coupling coefficents for the
C###    no to ny mappings for region nr and problem type nx. ie
C###    YP(nyo,nx)=XO(no,nx)*CYNO(nyo,no,nrc,nr,nx)+CYNO(0,no,nrc,nr,nx)

      IMPLICIT NONE
      INCLUDE 'aero00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'jtyp00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'solv00.cmn'
      INCLUDE 'tol00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),
     '  NBJ(NJM,NEM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),
     '  NHE(NEM),NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NLL(12,NEM),NNB(4,4,4,NBFM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NONY(0:NOYM,NYM,NRCM,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),NPF(9,NFM),NPL(5,0:3,NLM),
     '  NPNE(NNM,NBFM,NEM),NPNY(0:6,NYM,0:NRCM),nr,NRE(NEM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NWP(NPM,2),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM),NYNY(0:NYYM,NYM,NRM),
     '  NYQNR(0:NYQM,0:NRCM,NCM,0:NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  CYNY(0:NYYM,NYM,NRM),SE(NSM,NBFM,NEM),SP(NKM,NBFM,NPM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),
     '  YP(NYM,NIYM,NXM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER adj_dirn,adj_orderA,adj_orderB,
     '  connected_ne_nn,INCREASE_NOM,INCREASE_NYOM,LD,
     '  MAPPED_TOT,na,naadj,nae,nb,NBP,nc,ne,neadj,neadj2,ne_inlist,
     '  NELIST_max,nelowest,nepos_inlist, !SMAR009 23/12/98 nf,
     '  nffe,nh,ni,ni_sign,NITB,nk,nk1,
     '  nl,nn,nne,nne1,no,noadj,no_nelist,nonr,
     '  no_nynr,no_tot(2),no_wake,noy,
     '  NOYT,np,np1,nrr,nrbe,nrc,nrc_tot,ns,
     '  nu,NUNK(8),nv,ny,ny1r,ny1v,ny2,ny2v,ny3,nyadj,
     '  nyadj2,nybe,nynp1,nyo,nyr,nyr1,nyy1(2),nyy(2),ny_TE1_pot(2),
     '  ny_TE2_pot(2),ny_wake2_pot(2),ny_wake2_vel(2),nyy_wake(2)
      !SMAR009 22/12/98 connected_ne,nk2,nse,
      REAL*8 COY,RATIO,PSI,weight,XI(3),XD(3)
      CHARACTER CHAR*1
      LOGICAL BEM_REGION,INLIST,FOUND_NEW_ELEM,
     '  SPECIAL
!     External Functions
      INTEGER GETNYR,IDIGITS

      DATA NUNK/1,2,4,6,7,9,10,11/

      CALL ENTERS('GLOBALH',*9999)

      INCREASE_NOM=0
      INCREASE_NYOM=0

C***  Initialise mapping arrays above current region

      DO nrr=nr,NRT
        DO nrc=1,2
          DO no=1,NOM
            DO nyo=0,NYOM
              NYNO(nyo,no,nrc,nrr)=0
              CYNO(nyo,no,nrc,nrr)=0.0d0
            ENDDO !nyo
          ENDDO !no
          NOT(nrc,1,nr,nx)=0
          NOQT(nrc,1,nr,nx)=0
        ENDDO !nrc
        DO nc=1,NCT(nrr,nx) !GK, GQ (and GD) variables
          DO no_nynr=1,NYNR(0,0,nc,nrr) !Loop over variables in nrr
            ny=NYNR(no_nynr,0,nc,nrr) !variable #
            IF(ITYP4(nr,nx).NE.3) THEN !except finite differences
              ny2=GETNYR(1,NPNY,nrr,1,0,ny,NYNE,NYNP) !equiv rhs row #
            ELSE
              ny2=ny
            ENDIF
            nyy(1)=ny2
            nyy(2)=ny
            DO nrc=1,2
              DO noy=0,NOYM
                NONY(noy,nyy(nrc),nrc,nrr)=0
                CONY(noy,nyy(nrc),nrc,nrr)=0.0d0
              ENDDO !noy
            ENDDO !nrc
          ENDDO !no_nynr (ny)
        ENDDO !nc
      ENDDO !nrr

      DO nrc=1,2
        no_tot(nrc)=0
      ENDDO

C-Finite elements-------------------//----------------------------------
      IF(ITYP4(nr,nx).EQ.1) THEN !FEM
C       Standard FEM

C       Look at the FEM LHS variables to determine the number of
C       rows and columns in the reduced system. Flux variables do not
C       affect the mappings.

        nc=1 !Just want GK variables for FEM
        DO no_nynr=1,NYNR(0,0,1,nr) !Loop over the global vars for nr
          ny1v=NYNR(no_nynr,0,1,nr) !global var#
          ny2v=NYNR(no_nynr,0,2,nr) !global flux#
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' no_nynr='',I4,'' ny1v='',I8)')
     '        no_nynr,ny1v
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF

          ny1r=GETNYR(1,NPNY,NR,1,0,ny1v,NYNE,NYNP) !equivalent row number

          nyy(1)=ny1r !row#
          nyy(2)=ny1v !col#

          SPECIAL=.FALSE. !special ny's are mapped to other variables
          IF(CALL_AERO) THEN !aerofoil problem needs special mapping
            CALL ASSERT(NPNY(0,ny1v,0).EQ.1,
     '        '>>Element based ny''s not implemented',ERROR,*9999)

            np=NPNY(4,ny1v,0)
            IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' np='',I4)') np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(INLIST(np,NP_WAKE(1,1),NP_WAKE(0,1),no_wake)) THEN !wake
              nv=NPNY(2,ny1v,0)
              nk=NPNY(1,ny1v,0)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'('' nv='',I1,'' nk='',I1)') nv,nk
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
              IF(nv.EQ.1.AND.nk.EQ.1) THEN !lower wake potential
                IF(np.EQ.NP_WAKE(1,1)) THEN !first node
                  DO nrc=1,2
                    ny_TE2_pot(nrc)=nyy(nrc) !lower trailing edge pot
                  ENDDO !nrc
                  IF(DOP) WRITE(*,'('' ny_TE2_pot(nrc)='',2I4)')
     '              (ny_TE2_pot(nrc),nrc=1,2)
                ELSE !other nodes
                  DO nrc=1,2
                    ny_wake2_pot(nrc)=nyy(nrc) !lower wake pot ny
                  ENDDO !nrc
                  IF(DOP) WRITE(*,'('' ny_wake2_pot(nrc)='',2I4)')
     '              (ny_wake2_pot(nrc),nrc=1,2)
                ENDIF

              ELSE IF(nv.EQ.1.AND.nk.EQ.2) THEN !lower wake velocity
                DO nrc=1,2
                  ny_wake2_vel(nrc)=nyy(nrc) !lower wake velocity ny
                ENDDO !nrc
                IF(DOP) WRITE(*,'('' ny_wake2_vel(nrc)='',2I4)')
     '            (ny_wake2_vel(nrc),nrc=1,2)

              ELSE IF(nv.EQ.2.AND.nk.EQ.1) THEN !upper wake potential
                IF(np.EQ.NP_WAKE(1,1)) THEN !first node
                  DO nrc=1,2
                    ny_TE1_pot(nrc)=nyy(nrc) !upper trailing edge pot
                  ENDDO !nrc
                  IF(DOP) WRITE(*,'('' ny_TE1_pot(nrc)='',2I4)')
     '              (ny_TE1_pot(nrc),nrc=1,2)
                ELSE !other nodes
                  SPECIAL=.TRUE.
C                 Couple to lower wake potential & TE jump
                  DO nrc=1,2
                    nyy_wake(nrc)=ny_wake2_pot(nrc)
                    IF(DOP) WRITE(*,'('' nyy_wake('',I1,'')='',I5)')
     '                nrc,nyy_wake(nrc)
                    no=NONY(1,nyy_wake(nrc),nrc,nr) !is lower wake pot
                    IF(DOP) WRITE(*,'('' nrc='',I1,'' no='',I4)')nrc,no
                    NONY(0,nyy(nrc),nrc,nr)= 3
                    CONY(0,nyy(nrc),nrc,nr)= 0.0d0
                    NONY(1,nyy(nrc),nrc,nr)= no !lower wake potential
                    CONY(1,nyy(nrc),nrc,nr)= 1.0d0
                    NONY(2,nyy(nrc),nrc,nr)
     '                =NONY(1,ny_TE2_pot(nrc),nrc,nr) !lower TE pot
                    CONY(2,nyy(nrc),nrc,nr)=-1.0d0
                    NONY(3,nyy(nrc),nrc,nr)
     '                =NONY(1,ny_TE1_pot(nrc),nrc,nr) !upper TE pot
                    CONY(3,nyy(nrc),nrc,nr)= 1.0d0
                    NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                    CALL ASSERT(NYNO(0,no,nrc,nr).LE.NYOM,
     '                '>>Increase NYOM',ERROR,*9999)
                    NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                    CYNO(0,no,nrc,nr)=0.0d0
                    CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=1.0d0
                  ENDDO !nrc
C                 Couple to lower TE potential
                  DO nrc=1,2
                    nyy_wake(nrc)=ny_TE2_pot(nrc)
                    no=NONY(1,nyy_wake(nrc),nrc,nr)
                    IF(DOP) WRITE(*,'('' nrc='',I1,'' no='',I4)') nrc,no
                    NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                    CALL ASSERT(NYNO(0,no,nrc,nr).LE.NYOM,
     '                '>>Increase NYOM',ERROR,*9999)
                    NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                    CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=-1.0d0
                  ENDDO !nrc
C                 Couple to upper TE potential
                  DO nrc=1,2
                    nyy_wake(nrc)=ny_TE1_pot(nrc)
                    no=NONY(1,nyy_wake(nrc),nrc,nr)
                    IF(DOP) WRITE(*,'('' nrc='',I1,'' no='',I4)')nrc,no
                    NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                    CALL ASSERT(NYNO(0,no,nrc,nr).LE.NYOM,
     '                '>>Increase NYOM',ERROR,*9999)
                    NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                    CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=1.0d0
                  ENDDO !nrc
                ENDIF

              ELSE IF(nv.EQ.2.AND.nk.EQ.2) THEN !upper wake velocity
                SPECIAL=.TRUE. !all wake nodes
                DO nrc=1,2
                  nyy_wake(nrc)=ny_wake2_vel(nrc) !lower wake veloc
                  no=NONY(1,nyy_wake(nrc),nrc,nr)
                  NONY(0,nyy(nrc),nrc,nr)= 1
                  CONY(0,nyy(nrc),nrc,nr)= 0.0d0
                  NONY(1,nyy(nrc),nrc,nr)= no !lower wake veloc
                  CONY(1,nyy(nrc),nrc,nr)= 1.0d0
                  NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                  CALL ASSERT(NYNO(0,no,nrc,nr).LE.NYOM,
     '              '>>Increase NYOM',ERROR,*9999)
                  NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                  CYNO(0,no,nrc,nr)=0.0d0
                  CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=1.0d0
                ENDDO !nrc
              ENDIF !nv,nk
            ENDIF !in wake node list

          ELSE IF(ITYP1(nr,nx).EQ.5) THEN !FE50 problems
            IF(NPNY(0,ny1v,0).EQ.2) THEN !ny is element based
              na=NPNY(1,ny1v,0)
              nh=NPNY(2,ny1v,0)
              nc=NPNY(3,ny1v,0)
              ne=NPNY(4,ny1v,0)
              IF(ITYP2(nr,nx).EQ.8.AND.ITYP3(nr,nx).GE.3)THEN !cnst vol
                IF(nh.EQ.NHE(ne).AND.KTYP51(nr).EQ.3) THEN
C                 ny1v is pressure param; 3D case.
C                 Traverse through all elements adjacent to ne in
C                 current region to find elem with lowest number.
C                 NELIST stores all elements that have been scanned.
                  NELIST(0)=1
                  NELIST(1)=ne
                  nelowest=ne
                  FOUND_NEW_ELEM=.TRUE. !to get into the loop
                  DO WHILE(FOUND_NEW_ELEM)
                    FOUND_NEW_ELEM=.FALSE.
                    NELIST_max=NELIST(0)
                    DO no_nelist=1,NELIST_max
                      ne_inlist=NELIST(no_nelist)
                      DO ni=1,NIT(NBH(nh,nc,ne_inlist))
                        DO ni_sign=-1,1,2
                          neadj=NXI(ni*ni_sign,1,ne_inlist)
                          IF(neadj.NE.0) THEN !adjacent elem exists
                            IF(NRE(neadj).EQ.nr.AND.
     '                        .NOT.INLIST(neadj,NELIST(1),NELIST(0),
     '                        nepos_inlist)) THEN
C                             neadj is in reg nr and has not yet been
C                             traversed so add it to list
                              NELIST(0)=NELIST(0)+1
                              NELIST(NELIST(0))=neadj
                              IF(neadj.LT.nelowest) nelowest=neadj
                              FOUND_NEW_ELEM=.TRUE.
                            ENDIF !NRE(neadj).EQ.nr and .NOT.INLIST
                          ENDIF !neadj.NE.0
                        ENDDO !ni_sign
                      ENDDO !ni
                    ENDDO !no_nelist
                  ENDDO !FOUND_NEW_ELEM
                  IF(nelowest.LT.ne) THEN
C                   If element for current ny1v does not have lowest
C                   element number for adj Xi2 elems then couple ny1v to
C                   already created no (otherwise couple in
C                   normally - see .NOT.SPECIAL below)
                    SPECIAL=.TRUE.
                    DO nrc=1,2
C!!!                  Use 2-nrc in nrc position of NYNE as
C!!!                  want global quantities (nrc=0 global var;
C!!!                  nrc=1 global row)
                      nyadj=NYNE(na,nh,2-nrc,nc,nelowest)
                      noadj=NONY(1,nyadj,nrc,nr)
                      CALL ASSERT(noadj.NE.0,
     '                  '>>noadj is zero',ERROR,*9999)
                      NONY(0,nyy(nrc),nrc,nr)=1 !1 no cup to nyy(nrc)
                      NONY(1,nyy(nrc),nrc,nr)=noadj !nyy(nrc) cup to no
                      CONY(0,nyy(nrc),nrc,nr)=0.0d0 !no+coup 4 nyy(nrc)
                      CONY(1,nyy(nrc),nrc,nr)=1.0d0 !nyy(nrc) no mult
                      NYNO(0,noadj,nrc,nr)=NYNO(0,noadj,nrc,nr)+1
                      CALL ASSERT(NYNO(0,noadj,nrc,nr).LE.NYOM,
     '                  '>>Increase NYOM',ERROR,*9999)
                      NYNO(NYNO(0,noadj,nrc,nr),noadj,nrc,nr)=nyy(nrc)
                      CYNO(0,noadj,nrc,nr)=0.0d0
                      CYNO(NYNO(0,noadj,nrc,nr),noadj,nrc,nr)=1.0d0
                    ENDDO !nrc
                  ENDIF !nelowest.LT.ne
                ENDIF !KTYP51(nr).EQ.3.AND ...
              ELSE !other FE50 problem types
                IF(nh.EQ.NHE(ne).AND.KTYP51(nr).EQ.3.AND.
     '            KTYP52(nr).GE.2.AND.KTYP57(nr).GT.1) THEN
C                 ny1v is pressure param; 3D;
C                 Incomp (+fluid) with pressure bc increm.s entered
                  NBP=NBH(nh,nc,ne) !basis fn for pressure vars
                  IF(NAN(3,na,NBP).EQ.-1.OR.NAN(3,na,NBP).EQ.0) THEN
C                   ny is hyd press var assoc with an Xi3=0 face
                    adj_dirn=-3
                  ELSE IF(NAN(3,na,NBP).EQ.-2.OR.NAN(3,na,NBP).EQ.1
     '                .OR.NAN(3,na,NBP).EQ.3) THEN
C                   ny is hyd press var assoc with an Xi3=1 face
                    adj_dirn=3
                  ELSE
                    adj_dirn=0
                  ENDIF
                  neadj=NXI(adj_dirn,1,ne)
                  IF(neadj.EQ.0) THEN
C                   no element adjacent to associated face
                    IF(NAN(3,na,NBP).EQ.-1.OR. !Xi3=0 face bc
     '                NAN(3,na,NBP).EQ.-2) THEN !Xi3=1 face bc
                      SPECIAL=.TRUE.
C                     Don't couple current ny1v into solution
                      DO nrc=1,2
                        NONY(0,nyy(nrc),nrc,nr)=0
                        CONY(0,nyy(nrc),nrc,nr)=0.0d0
                      ENDDO
                    ENDIF !NAN=-1,-2
                  ELSE IF(neadj.NE.0.AND.neadj.NE.ne.AND.
     '                NRE(neadj).EQ.nr) THEN
C                   there is an adjacent element from current region
C                   in +/- Xi3 dirn
                    IF(NAN(3,na,NBP).EQ.-1.OR. !Xi3=0 face bc
     '                NAN(3,na,NBP).EQ.-2) THEN !Xi3=1 face bc
                      SPECIAL=.TRUE.
C                     Don't couple current ny1v into solution
                      DO nrc=1,2
                        NONY(0,nyy(nrc),nrc,nr)=0
                        CONY(0,nyy(nrc),nrc,nr)=0.0d0
                      ENDDO
                    ELSE IF(NAN(3,na,NBP).EQ.0.OR.
     '                  NAN(3,na,NBP).EQ.1.OR.
     '                  NAN(3,na,NBP).EQ.3) THEN
                      IF(neadj.LT.ne) THEN
C                       If element # for current ny1v is greater than
C                       that for adjacent element then couple ny1v to
C                       already created no (otherwise couple in
C                       normally - see .NOT.SPECIAL below)
                        IF(NAN(3,na,NBP).EQ.0) THEN
C                         hyd press var assoc with Xi3=0 face
                          adj_orderA=1
                          adj_orderB=3
                        ELSE IF(NAN(3,na,NBP).EQ.1.OR.
     '                      NAN(3,na,NBP).EQ.3) THEN
C                         hyd press var assoc with Xi3=1 face
                          adj_orderA=0
                          adj_orderB=0
                        ENDIF
                        naadj=1
                        DO WHILE(naadj.LE.NAT(NBP).AND..NOT.SPECIAL)
                          IF(NAN(3,naadj,NBP).EQ.adj_orderA.OR.
     '                      NAN(3,naadj,NBP).EQ.adj_orderB) THEN
                            SPECIAL=.TRUE.
                          ELSE
                            naadj=naadj+1
                          ENDIF
                        ENDDO
                        IF(SPECIAL) THEN
                          DO nrc=1,2 !global rows/variables
C!!!                        Use 2-nrc in nrc position of NYNE as
C!!!                        want global quantities (nrc=0 global var;
C!!!                        nrc=1 global row)
                            nyadj=NYNE(naadj,nh,2-nrc,nc,neadj)
                            noadj=NONY(1,nyadj,nrc,nr)
                            CALL ASSERT(noadj.NE.0,
     '                        '>>noadj is zero',ERROR,*9999)
                            NONY(0,nyy(nrc),nrc,nr)=1
                            NONY(1,nyy(nrc),nrc,nr)=noadj
                            CONY(0,nyy(nrc),nrc,nr)=0.0d0
                            CONY(1,nyy(nrc),nrc,nr)=1.0d0
                            IF(nrc.EQ.1) THEN
C                             Only set up no->ny map for rows here
                              NYNO(0,noadj,nrc,nr)=
     '                          NYNO(0,noadj,nrc,nr)+1
                              CALL ASSERT(NYNO(0,noadj,nrc,nr).LE.NYOM,
     '                          '>>Increase NYOM',ERROR,*9999)
                              NYNO(NYNO(0,noadj,nrc,nr),noadj,nrc,nr)=
     '                          nyy(nrc)
                              CYNO(0,noadj,nrc,nr)=0.0d0
                              CYNO(NYNO(0,noadj,nrc,nr),noadj,nrc,nr)=
     '                          1.0d0
                            ENDIF !nrc.EQ.1
                          ENDDO !nrc

C                         Set up no->ny map for variables
                          nrc=2 !global vars
                          IF(NAN(3,na,NBP).EQ.0) THEN
C                           ny1v is hyd press var assoc with Xi3=0
C                           face so couple it to all soln vars assoc
C                           with the Xi3=0,1 face hyd press vars in
C                           all elements in the -Xi3 dirn.
                            neadj2=NXI(-3,1,ne)
C CS 3/9/99 fixed bug
C                            DO WHILE(neadj2.NE.0.AND.NRE(neadj2).EQ.nr)
                            DO WHILE(neadj2.NE.0)
                              IF(NRE(neadj2).EQ.nr) THEN
C                             Check that current mapping implementation
C                             can handle element configuration.
C                             Elements must have increasing numbers
C                             in the dirn of +Xi3.
                                CALL ASSERT(NXI(3,1,neadj2).GT.neadj2,
     '                          '>>Soln maps cannot handle elem config.'
     '                            //' Renumber so that elem #s incr'
     '                            //' with +Xi3',ERROR,*9999)
                                nyadj=NYNE(na,nh,2-nrc,nc,neadj2)
                                nyadj2=NYNE(naadj,nh,2-nrc,nc,neadj2)
C!!!                          Use 2-nrc in nrc position of NYNE as
C!!!                          want global quantities (nrc=0 global var;
C!!!                          nrc=1 global row)
                                DO noadj=1,no_tot(nrc)
                                  DO nyo=1,NYNO(0,noadj,nrc,nr)
                                    IF(NYNO(NYNO(0,noadj,nrc,nr),
     '                                noadj,nrc,nr).EQ.nyadj.OR.
     '                                NYNO(NYNO(0,noadj,nrc,nr)
     '                                ,noadj,nrc,nr).EQ.nyadj2) THEN
                                      NYNO(0,noadj,nrc,nr)=
     '                                  NYNO(0,noadj,nrc,nr)+1
                                      CALL ASSERT(NYNO(0,noadj,nrc,nr)
     '                                  .LE.NYOM,'>>Increase NYOM',
     '                                  ERROR,*9999)
                                      NYNO(NYNO(0,noadj,nrc,nr),
     '                                  noadj,nrc,nr)=nyy(nrc)
                                      CYNO(0,noadj,nrc,nr)=0.0d0
                                      CYNO(NYNO(0,noadj,nrc,nr),
     '                                  noadj,nrc,nr)=1.0d0
                                    ENDIF !nyno=nyadj1 or nyadj2
                                  ENDDO !nyo
                                ENDDO !noadj
C                             Move to next -Xi3 element.
                                neadj2=NXI(-3,1,neadj2)
                              ENDIF !NRE(neadj2).EQ.nr
                            ENDDO !while
                          ELSE IF(NAN(3,na,NBP).EQ.1.OR.
     '                        NAN(3,na,NBP).EQ.3) THEN
C                           ny1v is hyd press var assoc with Xi3=1
C                           face so couple it to soln var normally
                            no_tot(nrc)=no_tot(nrc)+1
                            IF(no_tot(nrc).LE.NOM) THEN
                              NYNO(0,no_tot(nrc),nrc,nr)=1
                              NYNO(1,no_tot(nrc),nrc,nr)=nyy(nrc)
                              CYNO(0,no_tot(nrc),nrc,nr)=0.0d0
                              CYNO(1,no_tot(nrc),nrc,nr)=1.0d0
                            ENDIF
                          ENDIF !NAN=0/1,3
                        ENDIF !SPECIAL
                      ENDIF !neadj.LT.ne
                    ENDIF !NAN=-1,-2/0,1,3
                  ENDIF !neadj.NE.0,ne neadj in same region
                ENDIF !KTYP51(nr).EQ.3.AND ...
              ENDIF !const vol/other types
            ENDIF !NPNY(0,ny1v,0).EQ.2 : ny1v element based
          ENDIF !problem type

C***      new CS 15/4/98
          IF(JTYP2C.EQ.1) THEN
            IF(NPNY(0,ny1v,0).NE.2) THEN ! not element based
C             Check if node is hanging. If it is, it is special.
              IF(NWP(NPNY(4,ny1v,0),1).GT.0) SPECIAL=.TRUE.
            ENDIF
          ENDIF

          IF(FIX(ny1v,1)) THEN !ny set as essential bc
C           Delete row and column from the reduced system of eqns.
            SPECIAL=.TRUE.
            DO nrc=1,2 !rows and columns
              NONY(0,nyy(nrc),nrc,nr)=0
              CONY(0,nyy(nrc),nrc,nr)=0.0d0
            ENDDO !nrc
          ENDIF !FIX(ny1v,1)
          
          IF(ITYP5(nr,nx).EQ.2.AND. !lung capillary problem
     &      ITYP2(nr,nx).EQ.11.AND.ITYP3(nr,nx).GE.3) SPECIAL=.FALSE.

          IF(.NOT.SPECIAL.AND.JTYP2A.EQ.1) THEN !treat versions as coincident
C           JTYP2A=1 means dof nonstandard
C MPN 3Jul2003: dropping down of indices was all wrong here!!!
            CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),NKB,NKHE,NNB,
     '        NONY(0,1,1,nr),NPNE,NPNY,nr,NVHE,NVHP(1,1,1,nr),NXI,
     '        ny1v,ny,NYNP(1,1,1,1,0,1,nr),NYNY(0,1,nr),CYNY(0,1,nr),
     '        RATIO,FIX,*9999)
C MPN 3Jul2003 OLD
C            CALL GETEQVNONY(IDO,INP,NBH,NBJ,NENP(1,0,nr),NKB,NKHE,NNB,
C     '        NONY(0,1,1,nr),NPNE,NPNY(0,1,nr),nr,
C     '        NVHE,NVHP(1,1,1,nr),NXI,
C     '        ny1v,ny,NYNP,NYNY(0,1,nr),CYNY(0,1,nr),RATIO,FIX,*9999)
            IF(ny.NE.ny1v) THEN
C             There is an equivalent mesh degree of freedom for which an
C             no as already been assigned.
              SPECIAL=.TRUE.
              IF(ny.EQ.0) THEN !dof not used
                YP(ny1v,1,nx)=0.0d0
                FIX(ny1v,1)=.TRUE.
              ELSE
                YP(ny1v,1,nx)=RATIO*YP(ny,1,nx)
                IF(FIX(ny,1)) THEN
                  FIX(ny1v,1)=.TRUE.
                ELSE
C                 current mesh degree of freedom
                  nyy1(1)=GETNYR(nc,NPNY,nr,1,0,ny1v,NYNE,NYNP) !global row#
                  nyy1(2)=ny1v !global col#
C                 equivalent mesh degree of freedom
                  nyy(1)=GETNYR(nc,NPNY,nr,1,0,ny,NYNE,NYNP) !global row#
                  nyy(2)=ny !global col#
                  DO nrc=1,2 !nrc=1,2 local row and local column
                    nyr1=nyy1(nrc)
                    nyr=nyy(nrc)
                    NOYT=NONY(0,nyr,nrc,nr)
                    NONY(0,nyr1,nrc,nr)=NOYT !CYNO(0) is already 0
                    DO noy=1,NOYT
                      no=NONY(noy,nyr,nrc,nr)
                      NONY(noy,nyr1,nrc,nr)=no
                      COY=RATIO*CONY(noy,nyr,nrc,nr)
                      CONY(noy,nyr1,nrc,nr)=COY

                      IF(no.LT.NOM) THEN
                        nyo=NYNO(0,no,nrc,nr)+1
                        NYNO(0,no,nrc,nr)=nyo
                        IF(nyo.LE.NYOM) THEN
                          NYNO(nyo,no,nrc,nr)=nyr1
C                         assuming COY is 1 or -1
                          CYNO(nyo,no,nrc,nr)=COY
                        ELSE IF(nyo.GT.INCREASE_NYOM) THEN
                          INCREASE_NYOM=nyo
                        ENDIF
                      ENDIF
                    ENDDO ! noy
                  ENDDO !nrc
                ENDIF !FIX
              ENDIF !ny=0
            ENDIF !ny.NE.ny1v

          ELSE IF(IS_COUPLED(nx).AND..NOT.SPECIAL) THEN
            IF(NPNY(0,ny1v,0).EQ.1) THEN !ny is node based
              np=NPNY(4,ny1v,0)
              IF(NP_INTERFACE(np,0).GT.1) THEN !interface node
                BEM_REGION=.FALSE.
                DO nonr=1,NP_INTERFACE(np,0)
                  nrr=NP_INTERFACE(np,nonr)
                  IF(nrr.NE.nr.AND.ITYP4(nrr,nx).EQ.2) THEN
                    BEM_REGION=.TRUE.
                    nrbe=nrr
                  ENDIF
                ENDDO
                IF(BEM_REGION) THEN !FEM coupled to a BEM region
                  SPECIAL=.TRUE.
                  nk=NPNY(1,ny1v,0)
                  nv=NPNY(2,ny1v,0)
                  nh=NPNY(3,ny1v,0)
                  nc=NPNY(5,ny1v,0)
                  nybe=NYNP(nk,nv,nh,np,0,nc,nrbe)
                  IF(nybe.EQ.0) THEN
C                   Var in the other region does not exist so remove
C                   this variable from the system.
                    DO nrc=1,2 !rows and columns
                      NONY(0,nyy(nrc),nrc,nr)=0
                      CONY(0,nyy(nrc),nrc,nr)=0.0d0
                    ENDDO !nrc
                  ELSE
C                   Variable in the other region exists so add in flux
C                   variable to be solved for as well.
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
                    ny3=GETNYR(2,NPNY,nr,0,0,ny1v,NYNE,NYNP)
C                   ny3 is equivalent global rhs variable#
                    IF(ny3.NE.0) THEN
                      IF(.NOT.FIX(ny3,1)) THEN
                        no_tot(2)=no_tot(2)+1
                        NONY(0,ny3,2,nr)=1
                        NONY(1,ny3,2,nr)=no_tot(2)
                        CONY(0,ny3,2,nr)=0.0d0
                        CONY(1,ny3,2,nr)=1.0d0
                        IF(no_tot(2).LE.NOM) THEN
                          NYNO(0,no_tot(2),2,nr)=1
                          NYNO(1,no_tot(2),2,nr)=ny3
                          CYNO(0,no_tot(2),2,nr)=0.0d0
                          CYNO(1,no_tot(2),2,nr)=1.0d0
                        ENDIF
                      ENDIF
                    ENDIF
                  ENDIF
                ENDIF !bem_region
              ENDIF !np_interface(np,0)>1
            ENDIF !ny is element based
          ENDIF !IS_COUPLED(nx).AND..NOT.SPECIAL

          IF(.NOT.SPECIAL) THEN !add ny to list of solve variables
            nrc_tot=2
            IF(FIX(ny1v,1).AND.ITYP5(nr,nx).EQ.2.AND. !lung capillary problem
     &        ITYP2(nr,nx).EQ.11.AND.ITYP3(nr,nx).GE.3) nrc_tot=1 !(only row)
            
            DO nrc=1,nrc_tot !rows and columns
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

C KAT 5Aug99
            IF(ITYP1(nr,nx).EQ.3.AND.NIYFIXM.GE.2) THEN
              IF(FIX(ny2v,1).AND.FIX(ny2v,2)) THEN !flux=au bdry cond
CC PJH 23May99 Flux=au bdry cond
C            IF(ITYP1(nr,nx).EQ.3.AND.FIX(ny2v,1).AND.FIX(ny2v,2)) THEN !flux=au
!              IF(.NOT.SALU_CONSISTENCY(nx)) THEN
!C MLB 24May99 Check for Salu
                ny2v=GETNYR(2,NPNY,nr,0,0,ny1v,NYNE,NYNP) !is equiv global flux var#
                NONY(0,ny2v,2,nr)=1 !#ny's coupled into no for ny2v
                CONY(0,ny2v,2,nr)=0.d0
                NONY(1,ny2v,2,nr)=NONY(1,ny1v,2,nr) !u no# for flux term
                CONY(1,ny2v,2,nr)=1.d0
                IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$                call mp_setlock()
                  WRITE(OP_STRING,
     '              '(/'' Flux=au bc at ny1v,ny2v ='',2I8)') ny1v,ny2v
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' NONY(1,ny2v,2,nr)='',I8)')
     '              NONY(1,ny2v,2,nr)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$                call mp_unsetlock()
                ENDIF !DOP
!              ENDIF !not Salu
              ENDIF !flux=au bdry cond.
            ENDIF !ITYP1(nr,nx).EQ.3.AND.NIYFIXM.GE.2

          ENDIF !.NOT.SPECIAL

          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'(/'' no_nynr='',I4,'' no_tot(nrc)='',2I8)')
     '        no_nynr,no_tot(1),no_tot(2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDDO !no_nynr   --------- end of no_nynr loop for FEM case ----------

        IF(JTYP2C.EQ.1) THEN !hanging nodes
C         Set up hanging nodes
C          DO nc=1,NCT(nr,nx)
          DO nc=1,1 !Needs extended for problems not though SOLVE1
            DO no_nynr=1,NYNR(0,0,nc,nr) !Loop over the global vars
              ny1v=NYNR(no_nynr,0,nc,nr) !Global Variable # (global col)

            IF(NPNY(0,ny1v,0).EQ.2) THEN !ny is element based
C             write(*,*) "element based"
            ELSE

              np=NPNY(4,ny1v,0)
C              FIX_TEMP=.FALSE.
              IF(NWP(np,1).GT.0) THEN !hanging node
                ny1r=GETNYR(nc,NPNY,NR,1,0,ny1v,NYNE,NYNP)
                nyy(1)=ny1r !global row#
                nyy(2)=ny1v !global col#

                ne=NWP(np,1) ! element ny hangs in
                nk=NPNY(1,ny1v,0)
                nv=NPNY(2,ny1v,0)
                nh=NPNY(3,ny1v,0)
                nb=NBH(nh,nc,ne)
                NITB=NIT(nb)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),nr,NVJE(1,1,1,ne),SE(1,1,ne),
     '            XA(1,1,ne),XE,XP,ERROR,*9999)
                LD=0
                XD(1)=XP(1,1,1,np)
                XD(2)=XP(1,1,2,np)
                XD(3)=XP(1,1,3,np)
                DO ni=1,NITB
                  XI(ni)=0.5d0
                ENDDO
                CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr,0.d0,XE,
     '            XI,XI,XD,.FALSE.,ERROR,*9999)
                IF(NITB.EQ.2) THEN
                  nae=NWP(np,2)
                  nl=NLL(nae,ne) ! global line np hangs on
                  MAPPED_TOT=NNL(0,nae,nb)
                ELSE
                  nffe=NWP(np,2)
C                  nf=NFF(nffe,ne) ! global face np hangs on
                  MAPPED_TOT=NNF(0,nffe,nb)
                ENDIF

C               Set row (nrc=1) and column (nrc=2) mappings
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
     '                .AND.(nne.LE.NNT(nb)))
                      nne=nne+1
                    ENDDO
                    DO nk1=1,NKT(nne,nb)
                      nynp1=NYNP(nk1,nv,nh,np1,nrc,nc,nr)
                      IF(NONY(0,nynp1,nrc,nr).NE.0) THEN
                        DO noy=1,NONY(0,nynp1,nrc,nr)
                          no=NONY(noy,nynp1,nrc,nr)
                          nu=NUNK(nk)
                          weight=PSI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,nne,nk1,nu,XI)
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
     '                       NONY(0,nyy(nrc),nrc,nr)+1
                            CONY(0,nyy(nrc),nrc,nr)= 0.0d0

                           NONY(NONY(0,nyy(nrc),nrc,nr),nyy(nrc),nrc,nr)
     '                       =no
                           CONY(NONY(0,nyy(nrc),nrc,nr),nyy(nrc),nrc,nr)
     '                       =weight

                            NYNO(0,no,nrc,nr)=NYNO(0,no,nrc,nr)+1
                            CYNO(0,no,nrc,nr)=0.0d0

                            NYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=nyy(nrc)
                            CYNO(NYNO(0,no,nrc,nr),no,nrc,nr)=weight
                          ENDIF
                        ENDDO ! noy
                      ELSE
                        IF(nrc.EQ.2) THEN ! copied from above
                          nu=NUNK(nk)
                          weight=PSI(IBT(1,1,nb),IDO(1,1,0,nb),
     '                      INP(1,1,nb),nb,nne,nk1,nu,XI)
                          ns=0
                          DO nne1=1,nne-1
                            ns=ns+NKT(nne1,nb)
                          ENDDO
                          ns=ns+nk1
                          weight=weight*SE(ns,nb,ne)
                          ns=0
                          DO nne1=1,nne-1
                            ns=ns+NKT(nne1,nb)
                          ENDDO
                          ns=ns+nk
                          IF(DABS(SP(nk,nb,np)).LE.LOOSE_TOL) THEN
                            SP(nk,nb,np)=1.0d0
                          ENDIF
                          weight=weight/SP(nk,nb,np)
C??? might want this when all nodes mapped to have no no's
C                          YP(ny1v,1,nx)=YP(ny1v,1,nx)+YP(nynp1,1,nx)*
C     '                      weight

C                          FIX_TEMP=.TRUE.
                        ENDIF ! nrc=2
                      ENDIF ! NONY.NE.O
                    ENDDO ! nk
                  ENDDO ! nn
                ENDDO ! nrc
              ENDIF ! hanging
C              IF(FIX_TEMP) FIX(ny1v,1)=.TRUE.

            ENDIF

            ENDDO ! no_nynr
          ENDDO ! nc
        ENDIF ! JTYP2

C-Boundary elements----------------//-----------------------------------
      ELSE IF(ITYP4(nr,nx).EQ.2) THEN !BEM

        BEMOVERDETERMINED(nr,nx)=.FALSE.
        nc=1 !just look at the variables of GK (assumes GQ is no
        !bigger than GK)
        DO no_nynr=1,NYNR(0,0,1,nr) !loop over global variables of GK
          ny=NYNR(no_nynr,0,1,nr) !global variable#

          SPECIAL=.FALSE. !special ny's are mapped to other variables
C         Check if node is hanging. If it is, it is special.
          IF(NWP(NPNY(4,ny,0),1).GT.0) SPECIAL=.TRUE.
          IF(.NOT.SPECIAL) THEN

            ny1v=GETNYR(1,NPNY,nr,1,0,ny,NYNE,NYNP) !is local row# for global var ny
            ny2v=GETNYR(2,NPNY,nr,0,0,ny,NYNE,NYNP) !is equiv global flux variable #
            IF(NPNY(0,ny,0).EQ.1) THEN !nodal based ny
              nk=NPNY(1,ny,0)
            ELSE
              ERROR='>> Element dofs not implemented'
              GOTO 9999
            ENDIF
            IF(ny2v.GT.0) THEN
              IF(nk.EQ.4.AND.KTYP93(nc,nr).GT.0) THEN
                !A cross derivative for 3d bem.  A linearly independent
                !derivative equation cannot be formed.  Skip this case.
                !AJP/CPB 24/1/95
                NONY(0,ny1v,1,nr)=0 !nk=4 variables and equations are not
                NONY(0,ny,2,nr)=0
                NONY(0,ny2v,2,nr)=0
                CONY(0,ny1v,1,nr)=0.0d0
                CONY(0,ny,2,nr)=0.0d0
                CONY(0,ny2v,2,nr)=0.0d0
              ELSE
                IF(FIX(ny,1)) THEN !variable is set as a b.c.
                  IF(FIX(ny2v,1)) THEN !flux is also set as a b.c. e.g.
                                      !inverse problems
                    IF(SALU_CONSISTENCY(nx)) THEN
                      NONY(0,ny,1,nr)=0 !row is not in the system
                      CONY(0,ny,1,nr)=0.0d0 !no additive coupling
                      NONY(0,ny2v,1,nr)=0 !row is not in the system
                      CONY(0,ny2v,1,nr)=0.0d0 !no additive coupling
                      NONY(0,ny,2,nr)=0 !variable is not in the system
                      CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                      NONY(0,ny2v,2,nr)=0 !flux is not in the system
                      CONY(0,ny2v,2,nr)=0.0d0 !no additive coupling
                    ELSE
                      BEMOVERDETERMINED(nr,nx)=.TRUE.
                      no_tot(1)=no_tot(1)+1
                      NONY(0,ny1v,1,nr)=1 !one equation in the system
                      NONY(1,ny1v,1,nr)=no_tot(1)
                      CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                      CONY(1,ny1v,1,nr)=1.0d0 !one to one row coupling
                      IF(no_tot(1).LE.NOM) THEN
                        NYNO(0,no_tot(1),1,nr)=1 !one ny to sol row
                        NYNO(1,no_tot(1),1,nr)=ny1v !row no coupled to ny1v
                        CYNO(0,no_tot(1),1,nr)=0.0d0 !no additive coupling
                        CYNO(1,no_tot(1),1,nr)=1.0d0 !one to one coupling
                      ENDIF
                      NONY(0,ny,2,nr)=0 !variable is not in the system
                      CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                      NONY(0,ny2v,2,nr)=0 !flux is not in the system
                      CONY(0,ny2v,2,nr)=0.0d0 !no additive coupling
                    ENDIF
                  ELSE !standard BEM case
                    no_tot(1)=no_tot(1)+1
                    NONY(0,ny1v,1,nr)=1 !one equation in the system
                    NONY(1,ny1v,1,nr)=no_tot(1)
                    CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                    CONY(1,ny1v,1,nr)=1.0d0 !one to one row coupling
                    IF(no_tot(1).LE.NOM) THEN
                      NYNO(0,no_tot(1),1,nr)=1 !one ny to sol row
                      NYNO(1,no_tot(1),1,nr)=ny1v !row no coupled to ny1v
                      CYNO(0,no_tot(1),1,nr)=0.0d0 !no additive coupling
                      CYNO(1,no_tot(1),1,nr)=1.0d0 !one to one coupling
                    ENDIF
                    no_tot(2)=no_tot(2)+1
                    NONY(0,ny,2,nr)=0 !variable is not in the system
                    CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                    NONY(0,ny2v,2,nr)=1 !flux is in the system
                    NONY(1,ny2v,2,nr)=no_tot(2) !no coupled to the flux var
                    CONY(0,ny2v,2,nr)=0.0d0 !no additive coupling
                    CONY(1,ny2v,2,nr)=1.0d0 !one to one var coupling
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no_tot(2),2,nr)=1 !one ny to sol var
                      NYNO(1,no_tot(2),2,nr)=ny2v !no coupled to flux var.
                      CYNO(0,no_tot(2),2,nr)=0.0d0 !no additive coupling
                      CYNO(1,no_tot(2),2,nr)=1.0d0 !one to one coupling
                    ENDIF
                  ENDIF !FIX(ny2v
                ELSE IF(FIX(ny2v,1)) THEN !flux is set as a b.c.
                  no_tot(1)=no_tot(1)+1
                  NONY(0,ny1v,1,nr)=1 !one equation in the system
                  NONY(1,ny1v,1,nr)=no_tot(1)
                  CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                  CONY(1,ny1v,1,nr)=1.0d0 !one to one row coupling
                  IF(no_tot(1).LE.NOM) THEN
                    NYNO(0,no_tot(1),1,nr)=1 !one ny to sol row
                    NYNO(1,no_tot(1),1,nr)=ny1v !row no coupled to ny1v
                    CYNO(0,no_tot(1),1,nr)=0.0d0 !no additive coupling
                    CYNO(1,no_tot(1),1,nr)=1.0d0 !one to one coupling
                  ENDIF
                  no_tot(2)=no_tot(2)+1
                  NONY(0,ny2v,2,nr)=0 !flux is not in the system
                  CONY(0,ny2v,2,nr)=0.0d0 !no additive coupling
                  NONY(0,ny,2,nr)=1 !variable is in the system
                  NONY(1,ny,2,nr)=no_tot(2) !no coupled to the var
                  CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                  CONY(1,ny,2,nr)=1.0d0 !one to one var coupling
                  IF(no_tot(2).LE.NOM) THEN
                    NYNO(0,no_tot(2),2,nr)=1 !one ny to sol var
                    NYNO(1,no_tot(2),2,nr)=ny !no coupled to variable.
                    CYNO(0,no_tot(2),2,nr)=0.0d0 !no additive coupling
                    CYNO(1,no_tot(2),2,nr)=1.0d0 !one to one coupling
                  ENDIF
                ELSE !both variable and flux are free
                  IF(IS_COUPLED(nx)) THEN !coupled problem
                    no_tot(1)=no_tot(1)+1
                    NONY(0,ny1v,1,nr)=1 !the eqn is in the system
                    NONY(1,ny1v,1,nr)=no_tot(1) !no is coupled to the eqn
                    CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                    CONY(1,ny1v,1,nr)=1.0d0 !one to one row coupling
                    IF(no_tot(1).LE.NOM) THEN
                      NYNO(0,no_tot(1),1,nr)=1 !one ny coupled to the no
                      NYNO(1,no_tot(1),1,nr)=ny1v !no is coupled to the eqn
                      CYNO(0,no_tot(1),1,nr)=0.0d0 !no add soln var coupling
                      CYNO(1,no_tot(1),1,nr)=1.0d0 !one to one coupling
                    ENDIF
                    no_tot(2)=no_tot(2)+1
                    NONY(0,ny,2,nr)=1 !the var needs to be solved for
                    NONY(1,ny,2,nr)=no_tot(2) !no is coupled to the variable
                    CONY(0,ny,2,nr)=0.0d0 !no additive variable coupling
                    CONY(1,ny,2,nr)=1.0d0 !GK goes into GKK
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no_tot(2),2,nr)=1 !one ny coupled to the no
                      NYNO(1,no_tot(2),2,nr)=ny !no is coupled to the
                                                 !variable
                      CYNO(0,no_tot(2),2,nr)=0.0d0 !no add soln var coupling
                      CYNO(1,no_tot(2),2,nr)=1.0d0 !1-1 var coupling
                    ENDIF
                    no_tot(2)=no_tot(2)+1
                    NONY(0,ny2v,2,nr)=1 !the flux is in the system
                    NONY(1,ny2v,2,nr)=no_tot(2) !no is coupled to the flux
                    CONY(0,ny2v,2,nr)=0.0d0 !no add flux coupling
                    CONY(1,ny2v,2,nr)=1.0d0 !one to one var coupling
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no_tot(2),2,nr)=1 !one ny coupled to the no
                      NYNO(1,no_tot(2),2,nr)=ny2v !no is coupled to the
                                                 !flux
                      CYNO(0,no_tot(2),2,nr)=0.0d0 !no add soln var coupling
                      CYNO(1,no_tot(2),2,nr)=1.0d0 !1-1 var coupling
                    ENDIF
                  ELSE !not coupled
                    np=NPNY(4,ny,0)
                    nh=NPNY(3,ny,0)
                    nv=NPNY(2,ny,0)
                    nk=NPNY(1,ny,0)
                    WRITE(ERROR,'('' >> No boundary conditions set '
     '                //'for np='',I4,'', nh='',I1,'', nv='',I2,'
     '                //''', nk='',I1,'', nr='',I1)') np,nh,nv,nk,nr
                    GOTO 9999
                  ENDIF !coupled problem
                ENDIF !IF(FIX(ny) ...
              ENDIF !nk=4
            ELSE !ny2v=0 !(i.e. lower order interpolation used for
                      !normal derivative and hence no corresponding
                      !flux variable.  The only well-posed square
                      !problem is for this variable to be free).
              IF(FIX(ny,1)) THEN
                ERROR='>> Problem overspecified - not implemented yet'
                GOTO 9999
              ELSE
                no_tot(1)=no_tot(1)+1
                NONY(0,ny1v,1,nr)=1 !one equation in the system
                NONY(1,ny1v,1,nr)=no_tot(1)
                CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                CONY(1,ny1v,1,nr)=1.0d0 !one to one row coupling
                IF(no_tot(1).LE.NOM) THEN
                  NYNO(0,no_tot(1),1,nr)=1 !one ny to sol row
                  NYNO(1,no_tot(1),1,nr)=ny1v !row no coupled to ny1v
                  CYNO(0,no_tot(1),1,nr)=0.0d0 !no additive coupling
                  CYNO(1,no_tot(1),1,nr)=1.0d0 !one to one coupling
                ENDIF
                no_tot(2)=no_tot(2)+1
                NONY(0,ny,2,nr)=1 !variable is in the system
                NONY(1,ny,2,nr)=no_tot(2) !no coupled to the var
                CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                CONY(1,ny,2,nr)=1.0d0 !one to one var coupling
                IF(no_tot(2).LE.NOM) THEN
                  NYNO(0,no_tot(2),2,nr)=1 !one ny to sol var
                  NYNO(1,no_tot(2),2,nr)=ny !no coupled to variable.
                  CYNO(0,no_tot(2),2,nr)=0.0d0 !no additive coupling
                  CYNO(1,no_tot(2),2,nr)=1.0d0 !one to one coupling
                ENDIF
              ENDIF !FIX(ny) ...
            ENDIF !ny2v=0 loop
          ENDIF !.NOT.SPECIAL
        ENDDO !no_nynr (ny)


C LKC 8-APR-1999 Don't want to do this if we don't have hanging
C  nodes (SPECIAL mapping)
        IF(SPECIAL) THEN

C!!!! 12/7/2001 CS Beware I doubt this is correct anymore
C new CS 21/12/98 set up mappings for hanging nodes
C!!! Only implemented for standard bem with flux b.c.'s
          DO no_nynr=1,NYNR(0,0,1,nr) !Loop over the global vars
            ny=NYNR(no_nynr,0,1,nr) !Global Variable # (global col)
            np=NPNY(4,ny,0)

            IF(NWP(np,1).GT.0) THEN ! hanging
              ny1v=GETNYR(1,NPNY,nr,1,0,ny,NYNE,NYNP) !is local row# for global var ny
              ny2v=GETNYR(2,NPNY,nr,0,0,ny,NYNE,NYNP) !is equiv global flux variable #
              IF(FIX(ny1v,1)) THEN !variable is set as a b.c.
                CALL ASSERT(.FALSE.,'>> B.C. type for'
     '            //' hanging nodes not implemented',ERROR,
     '            *9999)
              ELSE IF(.NOT.FIX(ny2v,1)) THEN !flux is set as a b.c.
                CALL ASSERT(.FALSE.,'>> B.C. type for'
     '            //' hanging nodes not implemented',ERROR,
     '            *9999)
              ENDIF

              ne=NWP(np,1) ! element ny hangs in
              nk=NPNY(1,ny,0)
              nv=NPNY(2,ny,0)
              nh=NPNY(3,ny,0)
              nb=NBH(nh,nc,ne)
              NITB=NIT(nb)
              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),NPNE(1,1,ne),
     '          nr,NVJE(1,1,1,ne),SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,
     '          *9999)
              LD=0
              XD(1)=XP(1,1,1,np)
              XD(2)=XP(1,1,2,np)
              XD(3)=XP(1,1,3,np)
              DO ni=1,NITB
                XI(ni)=0.5d0
              ENDDO
              CALL DEXI_POINT(IBT,IDO,INP,LD,NBJ,ne,NITB,nr,0.d0,XE,XI,
     '          XI,XD,.FALSE.,ERROR,*9999)
              IF(NITB.EQ.2) THEN
                nae=NWP(np,2)
                nl=NLL(nae,ne) ! global line np hangs on
                MAPPED_TOT=NNL(0,nae,nb)
C!!!! NNL APPEARS TO TO BE WRONG FOR SECTORS SO JUST SET THIS FOR NOW
                MAPPED_TOT=2
              ELSE
                nffe=NWP(np,2)
c        SMAR009 23/12/98  nf=NFF(nffe,ne) ! global face np hangs on
                MAPPED_TOT=NNF(0,nffe,nb)
              ENDIF

C             Set row (nrc=1) and column (nrc=2) mappings
              NONY(0,ny,1,nr)=0
              NONY(0,ny,2,nr)=0
              DO nn=1,MAPPED_TOT !number of nodes ny is mapped to
                IF(NITB.EQ.2) THEN
                  np1=NPL(1+nn,1,nl)
                ELSE
                  connected_ne_nn=NNF(nn+1,nffe,nb)
                  np1=NPNE(connected_ne_nn,nb,ne)
                ENDIF
                nne=1
                DO WHILE ((np1.NE.NPNE(nne,nb,ne))
     '            .AND.(nne.LE.NNT(nb)))
                  nne=nne+1
                ENDDO

                nrc=1
                DO nk1=1,NKT(nne,nb)-KTYP93(1,nr)
                  nynp1=NYNP(nk1,nv,nh,np1,0,1,nr)
                  DO noy=1,NONY(0,nynp1,1,nr)
                    no=NONY(noy,nynp1,nrc,nr)
                    nu=NUNK(nk)
                    weight=PSI(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '                nb,nne,nk1,nu,XI)
                    ns=0
                    DO nne1=1,nne-1
                      ns=ns+NKT(nne1,nb)
                    ENDDO
                    ns=ns+nk1
                    IF(nu.EQ.1) THEN
                      weight=weight*SE(ns,nb,ne)
                    ELSEIF(nk1.eq.1) THEN
C!!!! NOT CORRECT FOR 3D
                      weight=weight/SE(ns+NPL(1,0,nl),nb,ne)
                    ENDIF

                    NONY(0,ny1v,1,nr)=1+NONY(0,ny1v,1,nr) !variable is mapped out
                    NONY(NONY(0,ny1v,1,nr),ny1v,1,nr)=no
                    CONY(0,ny1v,1,nr)=0.0d0 !no additive coupling
                    CONY(NONY(0,ny1v,1,nr),ny1v,1,nr)=weight !row coupling weight
                    IF(no_tot(1).LE.NOM) THEN
                      NYNO(0,no,1,nr)=1+NYNO(0,no,1,nr) !ny's to sol row
                      NYNO(NYNO(0,no,1,nr),no,1,nr)=ny1v !row no coupled to ny1v
                      CYNO(0,no,1,nr)=0.0d0 !no additive coupling
                      CYNO(NYNO(0,no,1,nr),no,1,nr)=weight !row coupling weight
                    ENDIF

                    NONY(0,ny2v,2,nr)=0 !flux is not in the system
                    CONY(0,ny2v,2,nr)=0.0d0 !no additive coupling
                    NONY(0,ny,2,nr)=1+NONY(0,ny,2,nr) !variable is mapped out
                    NONY(NONY(0,ny,2,nr),ny,2,nr)=no !no coupled to the var
                    CONY(0,ny,2,nr)=0.0d0 !no additive coupling
                    CONY(NONY(0,ny,2,nr),ny,2,nr)=weight !row coupling weight
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no,2,nr)=1+NYNO(0,no,2,nr) !ny's to sol var
                      NYNO(NYNO(0,no,2,nr),no,2,nr)=ny !no coupled to variable.
                      CYNO(0,no,2,nr)=0.0d0 !no additive coupling
                      CYNO(NYNO(0,no,2,nr),no,2,nr)=weight !row coupling weight
                    ENDIF
                  ENDDO ! noy
                ENDDO ! nk
              ENDDO ! nn
            ENDIF ! hanging
          ENDDO ! no_nynr
        ENDIF !SPECIAL

C-Finite differences--------------//-----------------------------------
      ELSE IF(ITYP4(nr,nx).EQ.3) THEN !Finite Difference
        DO no_nynr=1,NYQNR(0,0,1,nr)
          DO nrc=1,2 !rows and columns
            no_tot(nrc)=no_tot(nrc)+1
          ENDDO
        ENDDO

C-Collocation---------------------//-----------------------------------
      ELSE IF(ITYP4(nr,nx).EQ.4) THEN !Orthogonal Collocation
      ENDIF !ityp4

      DO nrc=1,2
        IF(ITYP4(nr,nx).EQ.3) THEN !Finite Difference Grid
          NOQT(nrc,1,nr,nx)=no_tot(nrc)
          CALL ASSERT(NOQT(nrc,1,nr,nx).LE.NYQM,'>> Increase NYQM',
     '      ERROR,*9999)
        ELSE
          IF(no_tot(nrc).LE.NOM) THEN
            NOT(nrc,1,nr,nx)=no_tot(nrc)
          ELSE IF(no_tot(nrc).GT.INCREASE_NOM) THEN
            INCREASE_NOM=no_tot(nrc)
          ENDIF
C          CALL ASSERT(NOT(nrc,1,nr,nx).LE.NOM,'>> Increase NOM',
C     '      ERROR,*9999)
        ENDIF
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

      CALL EXITS('GLOBALH')
      RETURN
 9999 CALL ERRORS('GLOBALH',ERROR)
      CALL EXITS('GLOBALH')
      RETURN 1
      END


