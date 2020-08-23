      SUBROUTINE GLOBALC(NAN,NBH,NHE,NONY,NP_INTERFACE,NPNODE,NPNY,
     '  NRE,NW,nx,NXI,NYNE,NYNO,NYNP,NYNR,CONY,CYNO,FIX,ERROR,*)

C#### Subroutine: GLOBALC
C###  Description:
C###    GLOBALC calculates the mapping between the mesh dofs, ny, and
C###    the solution dofs, no, for the entire coupled problem.  It
C###    sets up the arrays NONY, CONY and their inverse arrays NYNO
C###    and CYNO for region 0 (the coupled region).

C**** See GLOBALH for mapping array definitions.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp90.cmn'

!     Parameter List
      INTEGER NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NHE(NEM),
     '  NONY(0:NOYM,NYM,NRCM,0:NRM),NP_INTERFACE(0:NPM,0:3),
     '  NPNODE(0:NP_R_M,0:NRM),NPNY(0:6,NYM,0:NRCM),NRE(NEM),
     '  NW(NEM,3),nx,NXI(-NIM:NIM,0:NEIM,0:NEM),
     '  NYNO(0:NYOM,NOOPM,NRCM,0:NRM),
     '  NYNE(NAM,NHM,0:NRCM,NCM,NEM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM,0:NRM),CYNO(0:NYOM,NOOPM,NRCM,0:NRM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER bcindex(2),COUP_ROW_NRLIST(0:9),GETNYR,i,iface,na,naadj,
     '  nbaux,ncc,nc,ne,neadj(2),nh,nhaux,nk,no,no1,nonode,nonr,nonr2,
     '  nonr3,no_nynr,no_tot(2),noy,np,nr,nrr,nrc,nv,ny,ny1,ny2,nyo,
     '  nyy(2)
      INTEGER NR_MX
      PARAMETER (NR_MX=9)

      LOGICAL ADJACENT,BEMFEMCOUP,FEMFEMCOUP,FOUND,FOUNDny1

      CALL ENTERS('GLOBALC',*9999)

C***  Initialise mapping arrays for nr=0
      DO nrc=1,2
        DO no=1,NOM
          DO nyo=1,NYOM
            NYNO(nyo,no,nrc,0)=0
            CYNO(nyo,no,nrc,0)=0.0d0
          ENDDO
          NYNO(0,no,nrc,0)=0
          CYNO(0,no,nrc,0)=0.0d0
        ENDDO
        NOT(nrc,1,0,nx)=0
        no_tot(nrc)=0
      ENDDO !nrc
      DO nc=1,NCT(0,nx) !GK, GQ (and GD) variables
        DO no_nynr=1,NYNR(0,0,nc,0) !Loop over variables in region 0
          ny=NYNR(no_nynr,0,nc,0) !variable #
          ny2=GETNYR(1,NPNY,0,1,0,ny,NYNE,NYNP) !equiv row #
          nyy(1)=ny2
          nyy(2)=ny
          DO nrc=1,2
            DO noy=0,NOYM
              NONY(noy,nyy(nrc),nrc,0)=0
              CONY(noy,nyy(nrc),nrc,0)=0.0d0
            ENDDO !noy
          ENDDO !nrc
        ENDDO !no_nynr (ny)
      ENDDO !nc

      IF(DOP) THEN
        IF(KTYP90.EQ.2) THEN
          WRITE(OP_STRING,'(/'' Coupling is Coupled Laplace''''s '
     '      //'equation'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF(KTYP90.EQ.7) THEN
          WRITE(OP_STRING,'(/'' Coupling is Coupled pressure on '
     '      //'surface'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

C**** Rows (equations)
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Calculating row couplings'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF

C cpb 5/10/97 For coupled Laplace equations you want have to reorder
C the rows to avoid zero's on the diagonal. In general you will want
C the rows variables to be related to the column variables. For now
C just reorder the region list. This may have to be generalised later.
      CALL ASSERT(COUP_NRLIST(0,nx).LE.NR_MX,
     '  '>>Increase COUP_ROW_NRLIST size',ERROR,*9999)
      IF(KTYP90.EQ.2) THEN
        COUP_ROW_NRLIST(0)=0
        DO nonr=1,COUP_NRLIST(0,nx)
          nr=COUP_NRLIST(nonr,nx)
          DO nonode=1,NPNODE(0,nr)
            np=NPNODE(nonode,nr)
            DO nonr2=1,NP_INTERFACE(np,0)
              nrr=NP_INTERFACE(np,nonr2)
              FOUND=.FALSE.
              DO nonr3=1,COUP_NRLIST(0,nx)
                IF(COUP_NRLIST(nonr3,nx).EQ.nrr) FOUND=.TRUE.
              ENDDO
              IF(FOUND) THEN
                FOUND=.FALSE.
                DO nonr3=1,COUP_ROW_NRLIST(0)
                  IF(COUP_ROW_NRLIST(nonr3).EQ.nrr) FOUND=.TRUE.
                ENDDO !nonr3
                IF(.NOT.FOUND) THEN
                  COUP_ROW_NRLIST(0)=COUP_ROW_NRLIST(0)+1
                  COUP_ROW_NRLIST(COUP_ROW_NRLIST(0))=nrr
                ENDIF
              ENDIF
            ENDDO !nonr2
          ENDDO !nonode
        ENDDO !nonr
      ELSE
        COUP_ROW_NRLIST(0)=COUP_NRLIST(0,nx)
        DO nonr=1,COUP_NRLIST(0,nx)
          COUP_ROW_NRLIST(nonr)=COUP_NRLIST(nonr,nx)
        ENDDO !nonr
      ENDIF

      DO nonr=1,COUP_ROW_NRLIST(0) !loop over coupled regions
        nr=COUP_ROW_NRLIST(nonr) !is coupled nr
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Coupled region: '',I2)') nr
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
        DO no=1,NOT(1,1,nr,nx) !loop over eqtn of region nr
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Solution row: '',I5)') no
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO nyo=1,NYNO(0,no,1,nr)
            ny=NYNO(nyo,no,1,nr) !ny coupled to the solution row
            IF(DOP) THEN
              WRITE(OP_STRING,'('' nyo='',I2,'', ny='',I5)') nyo,ny
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(NPNY(0,ny,1).EQ.1) THEN !ny is node based
              nk=NPNY(1,ny,1)
              nv=NPNY(2,ny,1)
              nh=NPNY(3,ny,1)
              np=NPNY(4,ny,1)
              nc=NPNY(5,ny,1)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ny is a nodal equation'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' nk='',I1,'', nv='',I2,'
     '            //''', nh='',I1,'', np='',I5,'', nc='',I1)')
     '            nk,nv,nh,np,nc
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' NP_INTERFACE(np,0..):'',3(1X,I2))')
     '            (NP_INTERFACE(np,i),i=0,NP_INTERFACE(np,0))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(NP_INTERFACE(np,0).GT.1) THEN !interface node
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' np is an interface node'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(KTYP90.EQ.2.OR.  !Coupled Laplace
     '            KTYP90.EQ.7) THEN !Coupled pressure on surface
                  IF(nr.EQ.NP_INTERFACE(np,1)) THEN !master region
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nr is the master region'')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
C MPN 3Jul2003: 
                    IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
C MPN 3Jul2003 OLD                    no_tot(1)=no_tot(1)+1
                    IF(no_tot(1).LE.NOM) THEN
                      NYNO(0,no_tot(1),1,0)=NYNO(0,no,1,nr)
                      CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                      NYNO(nyo,no_tot(1),1,0)=NYNO(nyo,no,1,nr)
                      CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***                  coupled eqtn mapping is 1-1 with region eqtn mapping
                    ENDIF
                    DO noy=1,NONY(0,ny,1,nr)
                      NONY(noy,ny,1,0)=no_tot(1)
                      CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***                  coupled eqtn mapping is 1-1 with region eqtn mapping
                    ENDDO !noy
                    NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                    CONY(0,ny,1,0)=CONY(0,ny,1,nr)
                  ELSE !region coupled to master region
                    nrr=NP_INTERFACE(np,1)
                    FEMFEMCOUP=ITYP4(nr,nx).EQ.1.AND.ITYP4(nrr,nx).EQ.1
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' nr is not the master '
     '                  //'region, master region is '',I2)') nrr
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(FEMFEMCOUP) THEN
                      ny1=NYNP(nk,nv,nh,np,1,nc,nrr) !master region ny
                      no1=NONY(1,ny1,1,0) !master region no
                      IF(DOP) THEN
                        WRITE(OP_STRING,'('' Interface is FEM-FEM, '
     '                    //'master region ny='',I5,'', master region '
     '                    //'no='',I5)') ny1,no1
                        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      ENDIF
                      IF(KTYP90.EQ.2) THEN !Coupled Laplace
                        NONY(0,ny,1,0)=1 !one no coupled to slave ny
                        NONY(1,ny,1,0)=no1 !slave ny coup to master no
                        NYNO(0,no1,1,0)=NYNO(0,no1,1,0)+1 !add 1 extra
C                                                   ny to the coupled no
                        CALL ASSERT(NYNO(0,no1,1,0).LE.NYOM,
     '                    '>>Increase NYOM',ERROR,*9999)
                        NYNO(NYNO(0,no1,1,0),no1,1,0)=ny !slave ny
                        CONY(0,ny,1,0)=0.0d0 !no add coup for slave ny
                        CONY(1,ny,1,0)=1.0d0 !ny->no mult. const.
                        CYNO(0,no1,1,0)=0.0d0 !no add coup for mast. no
                        CYNO(NYNO(0,no1,1,0),no1,1,0)=1.0d0 !mult const.
                      ELSE IF(KTYP90.EQ.7) THEN !Coup press on surf.
C                       ignore equations from slave region rows
                        NONY(0,ny,1,0)=0 !zero no's coupled to slave ny
                        CONY(0,ny,1,0)=0.0d0 !no add coup for slave ny
C!!! CS 3/11/00 Adding check that no exists for master ny
                        IF (no1.NE.0) THEN
                          CYNO(0,no1,1,0)=0.0d0 !no add coup for mast. no
                        ENDIF
                      ENDIF
                    ELSE
C MPN 3Jul2003: 
                      IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
C MPN 3Jul2003 OLD                    no_tot(1)=no_tot(1)+1
                      IF(no_tot(1).LE.NOM) THEN
                        NYNO(0,no_tot(1),1,0)=NYNO(0,no,1,nr)
                        CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                        NYNO(nyo,no_tot(1),1,0)=NYNO(nyo,no,1,nr)
                        CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***                    coup eqtn mapping is 1-1 with region eqtn mapping
                      ENDIF
                      DO noy=1,NONY(0,ny,1,nr)
                        NONY(noy,ny,1,0)=no_tot(1)
                        CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***                    coup eqtn mapping is 1-1 with region eqtn mapping
                      ENDDO !noy
                      NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                      CONY(0,ny,1,0)=CONY(0,ny,1,nr)
                    ENDIF
                  ENDIF
                ELSE !all other coupling types (KTYP90.NE.2 and 7)
                  ERROR='>>Equation coupling not implemented'
                  GOTO 9999
                ENDIF
              ELSE !not an interface node
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' np is not an interface node'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
C MPN 3Jul2003: 
                IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
C MPN 3Jul2003 OLD                    no_tot(1)=no_tot(1)+1
                IF(no_tot(1).LE.NOM) THEN
                  NYNO(0,no_tot(1),1,0)=NYNO(0,no,1,nr)
                  CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                  NYNO(nyo,no_tot(1),1,0)=NYNO(nyo,no,1,nr)
                  CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***              coupled eqtn mapping is 1-1 with region eqtn mapping
                ENDIF
                DO noy=1,NONY(0,ny,1,nr)
                  NONY(noy,ny,1,0)=no_tot(1)
                  CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***              coupled eqtn mapping is 1-1 with region eqtn mapping
                ENDDO !noy
                NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                CONY(0,ny,1,0)=CONY(0,ny,1,nr)
              ENDIF !Interface
            ELSE !ny is element based
              na=NPNY(1,ny,1)
              nh=NPNY(2,ny,1)
              nc=NPNY(3,ny,1)
              ne=NPNY(4,ny,1)
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ny is an auxiliary element eqn'')')
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' na='',I1,'
     '            //''', nh='',I1,'', nc='',I1,'', ne='',I5)')
     '            na,nh,nc,ne
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(KTYP90.EQ.7) THEN !Coupled pressure on surface
C***            Couple auxiliary element equns from slave region
C***            to pressure bc in master region.
                IF(nr.EQ.COUP_NRLIST(1,nx)) THEN !master region
                  neadj(1)=NXI(-3,1,ne)
                  neadj(2)=NXI( 3,1,ne)
                  ADJACENT=.FALSE.
                  IF(neadj(1).NE.0) THEN
                    IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-1.AND.
     '                NRE(neadj(1)).NE.nr) ADJACENT=.TRUE.
C                   ny is aux param for Xi3=0 face pressure bc and
C                   there is an element from a slave region sharing
C                   the Xi3=0 face.
                  ENDIF
                  IF(neadj(2).NE.0) THEN
                    IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-2.AND.
     '                NRE(neadj(2)).NE.nr) ADJACENT=.TRUE.
C                   ny is aux param for Xi3=1 face pressure bc and
C                   there is an element from a slave region sharing
C                   the Xi3=1 face.
                  ENDIF
                  IF(ADJACENT) THEN
C                   Pressure bdy values comes from slave regions so
C                   do not couple this no into global list
                    NONY(0,ny,1,0)=0 !zero no's coupled to slave ny
                    CONY(0,ny,1,0)=0.0d0 !no add coup for slave ny
                    CYNO(0,no,1,0)=0.0d0 !no add coup for mast. no
C                   Set NW to indicate that a pressure bc is applied to
C                   element in master region
                    IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-1) THEN
C                     pressure bc on Xi3=0 face
                      IF(NW(ne,1).EQ.1) NW(ne,1)=2
                      IF(NW(ne,1).EQ.3) NW(ne,1)=4
                    ELSE IF(NAN(3,na,NBH(nh,nc,ne)).EQ.-2) THEN
C                     pressure bc on Xi3=1 face
                      IF(NW(ne,1).EQ.1) NW(ne,1)=3
                      IF(NW(ne,1).EQ.2) NW(ne,1)=4
                    ENDIF
                  ELSE !slave region is not adjacent to current elem
C                   put aux param equns straight into global mapping
C MPN 3Jul2003: 
                    IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
C MPN 3Jul2003 OLD                    no_tot(1)=no_tot(1)+1
                    IF(no_tot(1).LE.NOM) THEN
                      NYNO(0,no_tot(1),1,0)=NYNO(0,no,1,nr)
                      CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                      NYNO(nyo,no_tot(1),1,0)=NYNO(nyo,no,1,nr)
                      CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***                  coupled eqtn mapping is 1-1 with region eqtn mapping
                    ENDIF
                    DO noy=1,NONY(0,ny,1,nr)
                      NONY(noy,ny,1,0)=no_tot(1)
                      CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***                  coupled eqtn mapping is 1-1 with region eqtn mapping
                    ENDDO !noy
                    NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                    CONY(0,ny,1,0)=CONY(0,ny,1,nr)
                  ENDIF
                ELSE !slave region
C                 put aux param equns straight into global mapping.
C                 only add new no if first nyo for current region
                  IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
                  IF(no_tot(1).LE.NOM) THEN
                    NYNO(0,no_tot(1),1,0)=NYNO(0,no_tot(1),1,0)+1
                    CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                    NYNO(nyo,no_tot(1),1,0)=ny
                    CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***                coupled eqtn mapping is 1-1 with region eqtn mapping
                  ENDIF
                  DO noy=1,NONY(0,ny,1,nr)
                    NONY(noy,ny,1,0)=no_tot(1)
                    CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***                coupled eqtn mapping is 1-1 with region eqtn mapping
                  ENDDO !noy
                  NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                  CONY(0,ny,1,0)=CONY(0,ny,1,nr)
                ENDIF
              ELSE !all other coupling types (KTYP90.NE.7)
C!!!            assumes no coupling between element equns
C MPN 3Jul2003: 
                IF(nyo.EQ.1) no_tot(1)=no_tot(1)+1
C MPN 3Jul2003 OLD                    no_tot(1)=no_tot(1)+1
                IF(no_tot(1).LE.NOM) THEN
                  NYNO(0,no_tot(1),1,0)=NYNO(0,no,1,nr)
                  CYNO(0,no_tot(1),1,0)=CYNO(0,no,1,nr)
                  NYNO(nyo,no_tot(1),1,0)=NYNO(nyo,no,1,nr)
                  CYNO(nyo,no_tot(1),1,0)=CYNO(nyo,no,1,nr)
C***              coupled eqtn mapping is 1-1 with region eqtn mapping
                ENDIF
                DO noy=1,NONY(0,ny,1,nr)
                  NONY(noy,ny,1,0)=no_tot(1)
                  CONY(noy,ny,1,0)=CONY(noy,ny,1,nr)
C***              coupled eqtn mapping is 1-1 with region eqtn mapping
                ENDDO !noy
                NONY(0,ny,1,0)=NONY(0,ny,1,nr)
                CONY(0,ny,1,0)=CONY(0,ny,1,nr)
              ENDIF
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' no_tot(1)='',I5)') no_tot(1)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NYNO(0..,'',I5,'',1,0):'',5(1X,I5))')
     '          no_tot(1),(NYNO(i,no_tot(1),1,0),i=0,
     '          NYNO(0,no_tot(1),1,0))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NONY(0..,'',I5,'',1,0):'',5(1X,I5))')
     '          ny,(NONY(i,ny,1,0),i=0,NONY(0,ny,1,0))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' CYNO(0..,'',I5,'',1,0):'','
     '          //'5(1X,F6.3))') no_tot(1),(CYNO(i,no_tot(1),1,0),
     '          i=0,NYNO(0,no_tot(1),1,0))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' CONY(0..,'',I5,'',1,0):'','
     '          //'5(1X,F6.3))') ny,(CONY(i,ny,1,0),i=0,NONY(0,ny,1,0))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nyo
        ENDDO !no
      ENDDO !nonr

C**** Columns (variables)
      IF(DOP) THEN
        WRITE(OP_STRING,'(/'' Calculating variable couplings'')')
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      ENDIF
      DO nonr=1,COUP_NRLIST(0,nx) !loop over coupled regions
        nr=COUP_NRLIST(nonr,nx) !is coupled nr
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Coupled region: '',I2)') nr
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' Calculating variable couplings'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
C cpb 1/10/97 Reorder the coupled variables into groups of the same
C nc variable for the region i.e. have all the nc=1 variables of the
C region, then all the nc=2 variables, then the nc=1 variables of the
C next region. This generates a block matrix type pattern which can
C then be exploited later on in the solver.
        DO ncc=1,NCT(nr,nx)
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' ncc= '',I1)') ncc
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO no=1,NOT(2,1,nr,nx) !loop over variables of region nr
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Solution variable: '',I5)') no
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nyo=1,NYNO(0,no,2,nr) !loop over ny's coupled to the no
              ny=NYNO(nyo,no,2,nr) !ny coupled to the solution var.
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nyo='',I2,'', ny='',I5)') nyo,ny
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(NPNY(0,ny,0).EQ.1) THEN !ny is node based
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ny is a node dof'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                nk=NPNY(1,ny,0)
                nv=NPNY(2,ny,0)
                nh=NPNY(3,ny,0)
                np=NPNY(4,ny,0)
                nc=NPNY(5,ny,0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' nk='',I1,'', nv='',I2,'
     '              //''', nh='',I1,'', np='',I5,'', nc='',I1)')
     '              nk,nv,nh,np,nc
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' NP_INTERFACE(np,0..):'','
     '              //'3(1X,I2))') (NP_INTERFACE(np,i),i=0,
     '              NP_INTERFACE(np,0))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(nc.EQ.ncc) THEN
                  IF(NP_INTERFACE(np,0).GT.1) THEN !interface node
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' np is an interface node'')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    IF(KTYP90.EQ.2.OR. !Coupled Laplace
     '                KTYP90.EQ.7) THEN !Coupled pressure on surface
                      IF(nr.EQ.NP_INTERFACE(np,1)) THEN !master region
                        IF(DOP) THEN
                          WRITE(OP_STRING,
     '                      '('' nr is the master region'')')
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
C MPN 3Jul2003: 
                        IF(nyo.EQ.1) no_tot(2)=no_tot(2)+1
C MPN 3Jul2003 OLD                    no_tot(2)=no_tot(2)+1
                        IF(no_tot(2).LE.NOM) THEN
                          NYNO(0,no_tot(2),2,0)=NYNO(0,no,2,nr)
                          CYNO(0,no_tot(2),2,0)=CYNO(0,no,2,nr)
                          NYNO(nyo,no_tot(2),2,0)=NYNO(nyo,no,2,nr)
                          CYNO(nyo,no_tot(2),2,0)=CYNO(nyo,no,2,nr)
C                         coup var mapping is 1-1 with reg var mapping
                        ENDIF
                        DO noy=1,NONY(0,ny,2,nr)
                          NONY(noy,ny,2,0)=no_tot(2)
                          CONY(noy,ny,2,0)=CONY(noy,ny,2,nr)
C                         coup var mapping is 1-1 with reg var mapping
                        ENDDO !noy
                        NONY(0,ny,2,0)=NONY(0,ny,2,nr)
                        CONY(0,ny,2,0)=CONY(0,ny,2,nr)
                      ELSE !region coupled to master region
                        nrr=NP_INTERFACE(np,1) !master region
                        ny1=NYNP(nk,nv,nh,np,0,nc,nrr) !master region ny
                        IF(DOP) THEN
                          WRITE(OP_STRING,'('' nr is not the master '
     '                      //'region, master region is '',I2)') nrr
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          WRITE(OP_STRING,'('' Master region ny='','
     '                      //'I5)') ny1
                          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                        ENDIF
                        IF(ny1.EQ.0) THEN
C                         No equivalent master region ny. This can
C                         happen for example if you have FEM-BEM
C                         coupling in 3D with the BEM master region
C                         not having any cross derivatives but the
C                         slave FEM region does.
C cpb 13/5/95 Don't include these variables.
                        ELSE
                          FEMFEMCOUP=ITYP4(nr,nx).EQ.1.AND.
     '                      ITYP4(nrr,nx).EQ.1
                          BEMFEMCOUP=ITYP4(nr,nx).EQ.1.AND.
     '                      ITYP4(nrr,nx).EQ.2 !bem master, fem slave
                          no1=NONY(1,ny1,2,0) !master region no
                          IF(DOP) THEN
                            WRITE(OP_STRING,'('' Master region no='','
     '                        //'I5)') no1
                            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                          ENDIF
                          IF(FEMFEMCOUP) THEN
                            IF(DOP) THEN
                              WRITE(OP_STRING,
     '                          '('' Interface is FEM-FEM'')')
                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                            ENDIF
C!!! CS 3/11/00 Adding check that no exists for master ny
                            IF (no1.NE.0) THEN
                              NONY(0,ny,2,0)=1 !one no coupled to slave ny
                              NONY(1,ny,2,0)=no1 !slave ny is coup to mast no
                              NYNO(0,no1,2,0)=NYNO(0,no1,2,0)+1 !add 1 extra
C                                                     ny to the coupled no
                              CALL ASSERT(NYNO(0,no1,2,0).LE.NYOM,
     '                        '>>Increase NYOM',ERROR,*9999)
                              NYNO(NYNO(0,no1,2,0),no1,2,0)=ny !slave ny
                              CONY(0,ny,2,0)=0.0d0 !no add. coup for slave ny
                              CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                              CYNO(0,no1,2,0)=0.0d0 !no add. coup for mast no
                              CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0 !mult const.
                            ELSE
                              NONY(0,ny,2,0)=0 !ny not coupled to slave ny
                            ENDIF
                          ELSE !.NOT.FEMFEMCOUP
                            IF(DOP.AND.BEMFEMCOUP) THEN
                              WRITE(OP_STRING,
     '                          '('' Interface is BEM-FEM'')')
                              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                            ENDIF
                            NONY(0,ny,2,0)=1 !one no coupled to slave ny
                            NONY(1,ny,2,0)=no1 !slave ny is coup to mast no
                            NYNO(0,no1,2,0)=NYNO(0,no1,2,0)+1 !add 1 extra
C                                                   ny to the coupled no
                            CALL ASSERT(NYNO(0,no1,2,0).LE.NYOM,
     '                        '>>Increase NYOM',ERROR,*9999)
                            NYNO(NYNO(0,no1,2,0),no1,2,0)=ny !slave ny
                            CONY(0,ny,2,0)=0.0d0 !no add. coup for slave ny
                            CYNO(0,no1,2,0)=0.0d0 !no add. coup for mast no

C!!! Assumes that at an interface in 3D it is the s2 derivative that
C!!! has been reversed but in 2D it is the s1 derivative. Therefore when
C!!! mapping solution variables (XO) back to mesh variables (YP) we
C!!! need to reverse the appropriate signs of the solution variables
C!!! for storage in slave region YP locations.

C cpb 27/1/97 Adding bem-fem coupling
                            IF(nc.EQ.1) THEN !dependent variable
                              IF(NJT.EQ.2) THEN !2D
                                IF(nk.EQ.1) THEN
                                  CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                ELSE
                                  CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                                  IF(BEMFEMCOUP) THEN
                                    CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                  ELSE
                                    CYNO(NYNO(0,no1,2,0),no1,2,0)=-1.0d0
                                  ENDIF
                                ENDIF
                              ELSE IF(NJT.EQ.3) THEN
                                IF(nk.GE.3) THEN !s2 and cross
                                  CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                                  IF(BEMFEMCOUP) THEN
                                    CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                  ELSE
                                    CYNO(NYNO(0,no1,2,0),no1,2,0)=-1.0d0
                                  ENDIF
                                ELSE
                                  CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                ENDIF
                              ENDIF
                            ELSE IF(nc.EQ.2) THEN !flux variable
                              IF(NJT.EQ.2) THEN !2D
                                IF(nk.EQ.1) THEN
                                  CONY(1,ny,2,0)=-1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=-1.0d0
                                ELSE
                                  CONY(1,ny,2,0)=-1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                ENDIF
                              ELSE IF(NJT.EQ.3) THEN !3D
                                IF(nk.GE.3) THEN !s2 and cross
                                  CONY(1,ny,2,0)=-1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=1.0d0
                                ELSE
                                  CONY(1,ny,2,0)=-1.0d0 !ny->no mult. const.
                                  CYNO(NYNO(0,no1,2,0),no1,2,0)=-1.0d0
                                ENDIF
                              ENDIF
                            ENDIF
                          ENDIF
                        ENDIF
                      ENDIF
                    ELSE !all other coupling types (KTYP90.NE.2 and 7)
                      ERROR='>>Variable coupling not implemented'
                      GOTO 9999
                    ENDIF
                  ELSE !not an interface node
                    IF(DOP) THEN
                      WRITE(OP_STRING,
     '                  '('' np is not an interface node'')')
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
C MPN 3Jul2003: 
                    IF(nyo.EQ.1) no_tot(2)=no_tot(2)+1
C MPN 3Jul2003 OLD                    no_tot(2)=no_tot(2)+1
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no_tot(2),2,0)=NYNO(0,no,2,nr)
                      CYNO(0,no_tot(2),2,0)=CYNO(0,no,2,nr)
                      NYNO(nyo,no_tot(2),2,0)=NYNO(nyo,no,2,nr)
                      CYNO(nyo,no_tot(2),2,0)=CYNO(nyo,no,2,nr)
C***                  coupled var mapping is 1-1 with region var mapping
                    ENDIF
                    DO noy=1,NONY(0,ny,2,nr)
                      NONY(noy,ny,2,0)=no_tot(2)
                      CONY(noy,ny,2,0)=CONY(noy,ny,2,nr)
C***                  coupled var mapping is 1-1 with region var mapping
                    ENDDO !noy
                    NONY(0,ny,2,0)=NONY(0,ny,2,nr)
                    CONY(0,ny,2,0)=CONY(0,ny,2,nr)
                  ENDIF !interface
                ENDIF
              ELSE !ny is element based
                na=NPNY(1,ny,0)
                nh=NPNY(2,ny,0)
                nc=NPNY(3,ny,0)
                ne=NPNY(4,ny,0)
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ny is an element dof'')')
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' na='',I1,'
     '              //''', nh='',I1,'', nc='',I1,'', ne='',I5)')
     '              na,nh,nc,ne
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(nc.EQ.ncc) THEN
                  IF(KTYP90.EQ.7) THEN !Coupled pressure on surface
                    IF(nr.EQ.COUP_NRLIST(1,nx)) THEN !master region
C new MPN 18Aug2004
C                     Check whether current ny is associated with a 
C                     pressure bc auxiliary param
C                     Get basis # for aux param
                      nhaux=NHE(ne)
                      nbaux=NBH(nhaux,nc,ne)
                      IF(NAN(3,na,nbaux).LT.0) THEN
C                       na is press bc param on Xi3 face of
C                       master reg elem so don't couple it through
C                       here (it is handled when slave region press
C                       param is coupled in below)
                        NONY(0,ny,2,0)=0 !zero no's coupled to slave ny
                        CONY(0,ny,2,0)=0.0d0 !no add coup for slave ny
                        CYNO(0,no,2,0)=0.0d0 !no add coup for mast. no
C end MPN 18Aug2004
                      ELSE
C                       ny is not associated with a pressure bc so
C                       put aux params straight into global mapping
C MPN 3Jul2003: 
                        IF(nyo.EQ.1) no_tot(2)=no_tot(2)+1
C MPN 3Jul2003 OLD                      no_tot(2)=no_tot(2)+1
                        IF(no_tot(2).LE.NOM) THEN
                          NYNO(0,no_tot(2),2,0)=NYNO(0,no,2,nr)
                          CYNO(0,no_tot(2),2,0)=CYNO(0,no,2,nr)
                          NYNO(nyo,no_tot(2),2,0)=NYNO(nyo,no,2,nr)
                          CYNO(nyo,no_tot(2),2,0)=CYNO(nyo,no,2,nr)
C***                      coup var mapping is 1-1 with reg var mapping
                        ENDIF
                        DO noy=1,NONY(0,ny,2,nr)
                          NONY(noy,ny,2,0)=no_tot(2)
                          CONY(noy,ny,2,0)=CONY(noy,ny,2,nr)
C***                      coup var mapping is 1-1 with reg var mapping
                        ENDDO !noy
                        NONY(0,ny,2,0)=NONY(0,ny,2,nr)
                        CONY(0,ny,2,0)=CONY(0,ny,2,nr)
                      ENDIF
                    ELSE IF(nr.NE.COUP_NRLIST(1,nx)) THEN !slave region
C***                  Couple auxiliary element vars from slave region
C***                  to pressure bc vars in master region.
                      nrr=COUP_NRLIST(1,nx) !master region #
C                     only add new no if first nyo for current region
                      IF(nyo.EQ.1) no_tot(2)=no_tot(2)+1
                      IF(no_tot(2).LE.NOM) THEN
                        NYNO(0,no_tot(2),2,0)=NYNO(0,no_tot(2),2,0)+1
                        CYNO(0,no_tot(2),2,0)=CYNO(0,no_tot(2),2,0)+
     '                    CYNO(0,no,2,nr)
                        NYNO(NYNO(0,no_tot(2),2,0),no_tot(2),2,0)=ny
                        CYNO(NYNO(0,no_tot(2),2,0),no_tot(2),2,0)=
     '                    CYNO(nyo,no,2,nr)
C***                    coupled var mapping is 1-1 with reg var mapping
                      ENDIF
                      NONY(0,ny,2,0)=1 !one no coupled to slave ny
                      NONY(1,ny,2,0)=no_tot(2) !slave ny->master no
                      CONY(0,ny,2,0)=0.0d0 !no add coup for slave ny
                      CONY(1,ny,2,0)=1.0d0 !ny->no mult. const.
                      neadj(1)=NXI(-3,1,ne)
                      bcindex(1)=-2
                      neadj(2)=NXI( 3,1,ne)
                      bcindex(2)=-1
                      ADJACENT=.FALSE.
                      IF(neadj(1).NE.0) THEN
                        IF(NRE(neadj(1)).EQ.nrr) ADJACENT=.TRUE.
                      ENDIF
                      IF(neadj(2).NE.0) THEN
                        IF(NRE(neadj(2)).EQ.nrr) ADJACENT=.TRUE.
                      ENDIF
                      IF(ADJACENT) THEN
C                       there is an adjacent elem in the master reg
                        DO iface=1,2
                          IF(neadj(iface).NE.0) THEN
                            IF(NRE(neadj(iface)).EQ.nrr) THEN
C                             Element in master region is coupled through
C                             Xi3=1,0 face of current element so couple
C                             this ny with press bc param on Xi3=0,1 face
C                             of master reg elem.
C                             Get basis # for aux params in adj elem
                              nhaux=NHE(neadj(iface))
                              nbaux=NBH(nhaux,nc,neadj(iface))
C                             Loop through na's for adj element to find
C                             parameter for current face bc
                              FOUNDny1=.FALSE.
                              naadj=1
                              DO WHILE(naadj.LE.NAT(nbaux).AND.
     '                          .NOT.FOUNDny1)
                                IF(NAN(3,naadj,nbaux).EQ.
     '                            bcindex(iface)) THEN
C                                 naadj is press bc par on Xi3 face of
C                                 master reg elem
                                  FOUNDny1=.TRUE.
C                                 get master region ny1 number
                                  ny1=
     '                              NYNE(naadj,nhaux,0,nc,neadj(iface))
C new MPN 20-Jun
                                  IF(NONY(0,ny1,2,nrr).EQ.0) THEN
C old                              IF(FIX(ny1,1)) THEN
C                                   ny1 is not coupled to a soln no for
C                                   master region so couple it in.
C                                   get equivalent row number ny2
                                    ny2=GETNYR(1,NPNY,nrr,1,0,ny1,NYNE,
     '                                NYNP)
                                    nyy(1)=ny2
                                    nyy(2)=ny1
                                    DO nrc=1,2 !rows and columns
                                      NOT(nrc,1,nrr,nx)=
     '                                  NOT(nrc,1,nrr,nx)+1
                                      NONY(0,nyy(nrc),nrc,nrr)=1
                                      NONY(1,nyy(nrc),nrc,nrr)=
     '                                  NOT(nrc,1,nrr,nx)
                                      CONY(0,nyy(nrc),nrc,nrr)=0.0d0
                                      CONY(1,nyy(nrc),nrc,nrr)=1.0d0
                                      IF(NOT(nrc,1,nrr,nx).LE.NOM) THEN
                                        NYNO(0,NOT(nrc,1,nrr,nx),
     '                                    nrc,nrr)=1
                                        NYNO(1,NOT(nrc,1,nrr,nx),
     '                                    nrc,nrr)=nyy(nrc)
                                        CYNO(0,NOT(nrc,1,nrr,nx),
     '                                    nrc,nrr)=0.d0
                                        CYNO(1,NOT(nrc,1,nrr,nx),
     '                                    nrc,nrr)=1.d0
                                      ENDIF
                                    ENDDO !nrc
                                    DO nrc=1,2 !rows and columns
                                      CALL ASSERT(
     '                                  NOT(nrc,1,nrr,nx).LE.NOM,
     '                                  '>>Increase NOM',ERROR,*9999)
                                    ENDDO !nrc
C                                   Set FIX to false for master region
C                                   ny1 since no EXTERNAL pressure
C                                   increm's are now applied on current
C                                   naadj var (pressure bc)
                                    FIX(ny1,1)=.FALSE. !curr. value not fixed
                                    FIX(ny1,2)=.FALSE. !no solution increments
                                  ENDIF !FIX(ny1,1)
                                ENDIF
                                naadj=naadj+1
                              ENDDO !while loop through aux params in adj elem
                              IF(FOUNDny1) THEN
                                NONY(0,ny1,2,0)=1 !one no coupled to mast ny1
                                NONY(1,ny1,2,0)=no_tot(2) !master ny1->mast no
                                CONY(0,ny1,2,0)=0.0d0 !0 add coup for mast ny1
                                CONY(1,ny1,2,0)=1.0d0
                                IF(no_tot(2).LE.NOM) THEN
                                  NYNO(0,no_tot(2),2,0)=
     '                              NYNO(0,no_tot(2),2,0)+1
                                  CALL ASSERT(NYNO(0,no_tot(2),2,0).LE.
     '                              NYOM,'>>Increase NYOM',ERROR,*9999)
                                  no1=NONY(1,ny1,2,nrr) !master region no
                                  CYNO(0,no_tot(2),2,0)=
     '                              CYNO(0,no_tot(2),2,0)+
     '                              CYNO(0,no1,2,nrr)
                                  NYNO(NYNO(0,no_tot(2),2,0),
     '                              no_tot(2),2,0)=ny1
                                  CYNO(NYNO(0,no_tot(2),2,0),
     '                              no_tot(2),2,0)=CYNO(1,no1,2,nrr)
                                ENDIF
                              ENDIF !FOUNDny1
                            ENDIF !adjacent elem in same region
                          ENDIF !neadj(iface).ne.0
                        ENDDO !iface
                      ELSE !master reg is not adjacent to current elem
C                           put aux params straight into global mapping
                            no_tot(2)=no_tot(2)+1
                        IF(no_tot(2).LE.NOM) THEN
                          NYNO(0,no_tot(2),2,0)=NYNO(0,no,2,nr)
                          CYNO(0,no_tot(2),2,0)=CYNO(0,no,2,nr)
                          NYNO(nyo,no_tot(2),2,0)=NYNO(nyo,no,2,nr)
                          CYNO(nyo,no_tot(2),2,0)=CYNO(nyo,no,2,nr)
C***                      coup var mapping is 1-1 with reg var mapping
                        ENDIF
                        DO noy=1,NONY(0,ny,2,nr)
                          NONY(noy,ny,2,0)=no_tot(2)
                          CONY(noy,ny,2,0)=CONY(noy,ny,2,nr)
C***                      coup var mapping is 1-1 with reg var mapping
                        ENDDO !noy
                        NONY(0,ny,2,0)=NONY(0,ny,2,nr)
                        CONY(0,ny,2,0)=CONY(0,ny,2,nr)
                      ENDIF
                    ENDIF
                  ELSE !all other coupling types (KTYP90.NE.7)
C!!!                assumes no coupling between element params
C MPN 3Jul2003: 
                    IF(nyo.EQ.1) no_tot(2)=no_tot(2)+1
C MPN 3Jul2003 OLD                    no_tot(2)=no_tot(2)+1
                    IF(no_tot(2).LE.NOM) THEN
                      NYNO(0,no_tot(2),2,0)=NYNO(0,no,2,nr)
                      CYNO(0,no_tot(2),2,0)=CYNO(0,no,2,nr)
                      NYNO(nyo,no_tot(2),2,0)=NYNO(nyo,no,2,nr)
                      CYNO(nyo,no_tot(2),2,0)=CYNO(nyo,no,2,nr)
C***                  coupled var mapping is 1-1 with region var mapping
                    ENDIF
                    DO noy=1,NONY(0,ny,2,nr)
                      NONY(noy,ny,2,0)=no_tot(2)
                      CONY(noy,ny,2,0)=CONY(noy,ny,2,nr)
C***                  coupled var mapping is 1-1 with region var mapping
                    ENDDO !noy
                    NONY(0,ny,2,0)=NONY(0,ny,2,nr)
                    CONY(0,ny,2,0)=CONY(0,ny,2,nr)
                  ENDIF
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' no_tot(2)='',I5)') no_tot(2)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' NYNO(0..,'',I5,'',2,0):'','
     '            //'5(1X,I5))') no_tot(2),(NYNO(i,no_tot(2),2,0),i=0,
     '            NYNO(0,no_tot(2),2,0))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' NONY(0..,'',I5,'',2,0):'','
     '            //'5(1X,I5))') ny,(NONY(i,ny,2,0),i=0,NONY(0,ny,2,0))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' CYNO(0..,'',I5,'',2,0):'','
     '            //'5(1X,F6.3))') no_tot(2),(CYNO(i,no_tot(2),2,0),
     '            i=0,NYNO(0,no_tot(2),2,0))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' CONY(0..,'',I5,'',2,0):'',5(1X,F6.3))')
     '            ny,(CONY(i,ny,2,0),i=0,NONY(0,ny,2,0))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !nyo
          ENDDO !no
        ENDDO !ncc
      ENDDO !nonr

      DO nrc=1,2
        NOT(nrc,1,0,nx)=no_tot(nrc)
        CALL ASSERT(NOT(nrc,1,0,nx).LE.NOM,'>>Increase NOM',
     '    ERROR,*9999)
      ENDDO !nrc

      CALL EXITS('GLOBALC')
      RETURN
 9999 CALL ERRORS('GLOBALC',ERROR)
      CALL EXITS('GLOBALC')
      RETURN 1
      END

