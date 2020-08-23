      SUBROUTINE UPGEOM(CP,NKH,NKJ,NPLIST, NPNODE,NPNY,NRLIST,NVHP,NVJP
     &     ,NXLIST, NYNP,NYNR,XP,YP,ZD,STRING,ERROR,*)

C#### Subroutine: UPGEOM
C###  Description:
C###    Updates geometry variables from field or solution variables.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     '  NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM),
     '  NPNY(0:6,NYM,0:NRCM,NXM),NRLIST(0:NRM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJP(NJM,NPM),NXLIST(0:NXM),
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM),
     '  NYNR(0:NY_R_M,0:NRCM,NCM,0:NRM,NXM)
      REAL*8 CP(NMM,NPM,NXM),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM),
     &  ZD(NJM,NDM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER FIELLIST(0:3),GEOMLIST(0:3),IBEG,IEND,IFROMC,iy,
     '  N3CO,nj,nj1,nj2,njj_FIEL,njj_GEOM,nk,NJ_OFFSET,
     '  no_njj,no_nrlist,no_nynr,nonode,np,nr,nv,nx,nxc,ny,PART2
      LOGICAL ABBREV,ALL_REGIONS,CBBREV
      CHARACTER OPERATION*16,TYPE*16,UPDATE*16

      CALL ENTERS('UPGEOM',*9999)

 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

        UPDATE='HELPGEOMETRY'
        CALL UPFG(UPDATE,%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     '    %VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),%VAL(0),
     &    STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update geometry <#s/all[all]> from field
C###  Parameter:      <field_numbers #s/all[all]>
C###    Specify the field avriable numbers
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Description: Updates geometry variables from field variables
C###

        OP_STRING(1)=STRING(1:IEND)//' <#s/all[all]> from field'
        OP_STRING(2)=BLANK(1:15)//'<field_numbers (#s/all[all])>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM update geometry from solution
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the region numbers to update.
C###  Parameter:      <class #[1]>
C###    Specify the class number (of solve type) to use.
C###  Parameter:      <to (nj#s)[1..3]>
C###    Specify the nj #'s of XP to update.
C###  Description: Updates geometry variables from solution variables
C###

        OP_STRING(1)=STRING(1:IEND)//' from solution'
        OP_STRING(2)=STRING(1:IEND)//'<YP_index #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<class #[1]>'
        OP_STRING(3)=BLANK(1:15)//'<to (nj#s)[1..3]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe23','doc','UPGEOM',ERROR,*9999)
      ELSE

        IF(CBBREV(CO,'SUBSTITUTE',4,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBSTITUTE'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'ADD',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='ADD'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'SUBTRACT',4,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='SUBTRACT'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'MULTIPLY',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='MULTIPLY'
          PART2=N3CO
        ELSEIF(CBBREV(CO,'DIVIDE',2,noco+1,NTCO,N3CO)) THEN
          TYPE='OPERATION'
          OPERATION='DIVIDE'
          PART2=N3CO

        ELSEIF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN

          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,
     '      ERROR,*9999)
          nr=NRLIST(0)
          CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
          nxc=NXLIST(1)

          IF(CBBREV(CO,'FROM',2,noco,noco+1,N3CO)) THEN
            GEOMLIST(0)=NJ_LOC(NJL_GEOM,0,nr)
            DO nj=1,GEOMLIST(0)
              GEOMLIST(nj)=nj
            ENDDO !nj
          ELSE
            CALL PARSIL(CO(noco+1),3,GEOMLIST(0),GEOMLIST(1),
     '        ERROR,*9999)
          ENDIF

          IF(ABBREV(CO(n3co+1),'FIELD',2)) THEN
            TYPE='FIELD'

C LKC 4-NOV-98 Problem here as the field numbers do not have an
C              identifier. So add a new tag "FIELD_NUMBERS"
C
C            CALL PARSIL(CO(N3CO+2),3,FIELLIST(0),FIELLIST(1),
C     '        ERROR,*9999)
C            IF(FIELLIST(0).EQ.0) THEN
C              FIELLIST(0)=NJ_LOC(NJL_FIEL,0,nr)
C              DO nj=1,FIELLIST(0)
C                FIELLIST(nj)=nj
C              ENDDO !nj

            IF(CBBREV(CO,'FIELD_NUMBERS',7,noco+1,NTCO,N3CO)) THEN
              IF(ABBREV(CO(n3co+1),'ALL',2)) THEN
                FIELLIST(0)=NJ_LOC(NJL_FIEL,0,nr)
                DO nj=1,FIELLIST(0)
                  FIELLIST(nj)=nj
                ENDDO !nj
              ELSE
                CALL PARSIL(CO(N3CO+2),3,FIELLIST(0),FIELLIST(1),
     '            ERROR,*9999)
              ENDIF ! ALL
            ELSE ! default to ALL
              FIELLIST(0)=NJ_LOC(NJL_FIEL,0,nr)
              DO nj=1,FIELLIST(0)
                FIELLIST(nj)=nj
              ENDDO !nj
            ENDIF !FIELD_NUMBERS

          ELSE IF(ABBREV(CO(n3co+1),'SOLUTION',2)) THEN
            TYPE='SOLUTION'
          ELSE
            CO(noco+1)='?'
            GO TO 1
          ENDIF
        ELSE
          CO(noco+1)='?'
          GO TO 1
        ENDIF

        IF(TYPE(1:9).EQ.'OPERATION') THEN

          UPDATE='GEOMETRY'
          CALL UPFG(UPDATE,CP,%VAL(0),%VAL(0),NKH,NKJ,NPLIST,
     '      NPNODE,NRLIST,NVHP,NVJP,NYNP,OPERATION,PART2,%VAL(0),XP,
     '      %VAL(0),YP,ZD,STRING,ERROR,*9999)

        ELSEIF(TYPE(1:5).EQ.'FIELD') THEN
          DO no_nrlist=1,NRLIST(0)
            nr=NRLIST(no_nrlist)
            CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GT.0,
     '        '>>Define field first',ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C CS 13/7/98 modifying to specify which geom and field
C              DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
              DO no_njj=1,GEOMLIST(0)
                njj_GEOM=GEOMLIST(no_njj)
                njj_FIEL=FIELLIST(no_njj)
                nj1=NJ_LOC(NJL_FIEL,njj_FIEL,nr)
                nj2=NJ_LOC(NJL_GEOM,njj_GEOM,nr)
                DO nv=1,NVJP(nj1,np)
                  DO nk=1,NKJ(nj1,np)
                    XP(nk,nv,nj2,np)=XP(nk,nv,nj1,np)
                  ENDDO !nk
                ENDDO !nv
              ENDDO !njj
            ENDDO !nonode (np)
          ENDDO !nr

        ELSE IF(TYPE(1:8).EQ.'SOLUTION') THEN
          IF(CBBREV(CO,'YP_INDEX',2,noco+1,NTCO,N3CO)) THEN
            iy=IFROMC(CO(N3CO+1))
          ELSE
            iy=1
          ENDIF

          CALL ASSERT(iy.GT.0.AND.iy.LE.6,
     '      '>>IY out of range',ERROR,*9999)
          CALL NX_LOC(NX_INQUIRE,nxc,nx,NX_SOLVE,ERROR,*9999)
          CALL ASSERT(nx.GT.0,'>>No nx defined for this solve class',
     '      ERROR,*9999)

C JWF 29.4.02 Allows user to select nj's of XP to update.
          NJ_OFFSET=0
          IF(CBBREV(CO,'TO',2,noco+1,NTCO,N3CO)) THEN
            CALL PARSIL(CO(N3CO+1),9,FIELLIST(0),FIELLIST(1), ERROR,
     &           *9999) 
C
C  OR 02/09/06
C
C     REMOVED THE MORE STRINGENT REQUIREMENT THAT field VARIABLE nj
C     NEEDS TO BE SMALLER THAN 9.
C     
            CALL ASSERT(FIELLIST(1).GT.0.AND.FIELLIST(1).LT.NJ_LOC(0,0
     &           ,nr),'>>ERROR: Define more fields.',ERROR,*9999) 
            IF (NPNY(3,NYNR(1,0,1,nr,nx),0,nx).NE.FIELLIST(1)) THEN
              NJ_OFFSET=FIELLIST(1)-NPNY(3,NYNR(1,0,1,nr,nx),0,nx)
            ENDIF
          ENDIF

          DO no_nrlist=1,NRLIST(0) !loop over regions
            nr=NRLIST(no_nrlist)
            DO no_nynr=1,NYNR(0,0,1,nr,nx) !loop over global vars
              ny=NYNR(no_nynr,0,1,nr,nx) !is global var number
              nk=NPNY(1,ny,0,nx)
              nv=NPNY(2,ny,0,nx)
              nj=NPNY(3,ny,0,nx)
              np=NPNY(4,ny,0,nx)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C
C OR 04/09/06 
C
C     Using trilinear basis functions for representing pressure terms
C     for finite elasticity problems leads to the fact that nj takes a
C     value greater than 3. In this case, a gloabal variable is
C     associated with the pressure variable and not with a geometrical
C     Xj variable. Hence the additional check for
C     nj.LE.NJ_LOC(NJL_GEOM,1..,nr)
C              
C             IF (NPNY(0,ny,0,nx).EQ.1) THEN
              IF (NPNY(0,ny,0,nx).EQ.1.AND.nj.LE.NJ_LOC(NJL_GEOM,0,nr)
     &             ) THEN
                nj=nj+NJ_OFFSET
                XP(nk,nv,nj,np)=YP(ny,iy,nx)
              ENDIF
            ENDDO !no_nynr
          ENDDO !no_nrlist

        ENDIF !type

      ENDIF

      CALL EXITS('UPGEOM')
      RETURN
 9999 CALL ERRORS('UPGEOM',ERROR)
      CALL EXITS('UPGEOM')
      RETURN 1
      END


