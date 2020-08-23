      SUBROUTINE UPFIEL_EIK(INP,NBH,NBHF,NBJ,NBJF,NEELEM,NEL,NFF,NFFACE,
     '  NFLIST,NKB,NKHE,NKJ,NKJE,NLF,NLL,NLLINE,NLLIST,NNF,NNL,NPF,NPNE,
     '  NPNODE,NRLIST,NVHE,NVHP,NVJE,NVJP,nx,NYNP,XP,XP_TEMP,FIX,
     '  ERROR,*)

C#### Subroutine: UPFIEL_EIK
C###  Description:
C###    Sets up field variable 1 from dependent variable versions and
C###    boundary conditions for solution of an eikonal equation.  The
C###    field is C0 continuous, is 1 over most of the domain but
C###    approaches zero at essential boundary conditions and on faces
C###    where out-of-face derivatives are not continous.  It is assumed
C###    that cubic Hermite elements are used for the dependent variable.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INP(NNM,NIM,NBFM),NBH(NHM,NCM,NEM),NBHF(NHM,NCM,NFM),
     '  NBJ(NJM,NEM),NBJF(NJM,NFM),NEELEM(0:NE_R_M,0:NRM),
     '  NEL(0:NELM,NLM),NFF(6,NEM),NFFACE(0:NF_R_M,NRM),NFLIST(0:NFM),
     '  NKB(2,2,2,NNM,NBFM),NKHE(NKM,NNM,NHM,NEM),
     '  NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),NLF(4,NFM),
     '  NLL(12,NEM),NLLINE(0:NL_R_M,0:NRM),NLLIST(0:NLM),
     '  NNF(0:17,6,NBFM),NNL(0:4,12,NBFM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRLIST(0:NRM),
     '  NVHE(NNM,NBFM,NHM,NEM),NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 XP(NKM,NVM,NJM,NPM),XP_TEMP(NKM)
      LOGICAL FIX(NYM,NIYFIXM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FACTOR,IDOI(3),nae,naf,NB_EA(2),nb_f,nb_he,nb_je,nc,
     '  ne,NEA(2),neax,NELT,nf,nf_e,NF_EA(2),nh,nhx,ni,ni2,ni3,NITB,nj,
     '  njj,nk,nk_e,NKTOT,nl,NLET,nn_f,nn_he,nn_je,nn_l,NNLT,NNTB,no_ne,
     '  no_nf,no_nl,no_np,no_nr,np,nr,nv,NVA(2),NVL(4),NVTOT,ny,
     '  NZ_LINE,NZ_FACE,SIGN(3)
      LOGICAL CONTINUOUS,FACE_ZERO(3),FIRST,FIXED,FOUND,INCREASE_NVM,
     '  LINE_ZERO(3),MULTNV

CC     Local adjacent face numbers at each local node.
C      DATA NFEANN/
C     '  1,3,5
C     '  2,3,5
C     '  1,4,5
C     '  2,4,5
C     '  1,3,6
C     '  2,3,6
C     '  1,4,6
C     '  2,4,6/

      CALL ENTERS('UPFIEL_EIK',*9999)

      nc=1 !essential dependent variable boundary conditions
      nhx=NH_LOC(0,nx)
      CALL ASSERT(nhx.EQ.1,'>>Incorrect number of dependent variables',
     '  ERROR,*9999)
      nh=NH_LOC(nhx,nx)

      DO no_nr=1,NRLIST(0)
        nr=NRLIST(no_nr)
        CALL ASSERT(NJ_LOC(NJL_GEOM,0,nr).EQ.3,
     '    '>>Only implemented for 3D',ERROR,*9999)
        njj=1
        CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).GE.njj,
     '    '>>Continuity field not defined',ERROR,*9999)
        nj=NJ_LOC(NJL_FIEL,njj,nr)
C       Initialize NFLIST,NLLIST to 1 to indicate that the field is
C       non-zero on each face,line.
        DO no_nf=1,NFFACE(0,nr)
          NFLIST(NFFACE(no_nf,nr))=1
        ENDDO !nf
        DO no_nl=1,NLLINE(0,nr)
          NLLIST(NLLINE(no_nl,nr))=1
        ENDDO !nl
C       Initialize NVJP to -1 to indicate that the default field value
C       is to be used.
        DO no_np=1,NPNODE(0,nr)
          NVJP(nj,NPNODE(no_np,nr))=-1
        ENDDO
      ENDDO !nr

C     Find out which nodes are not set to the default
      DO no_nr=1,NRLIST(0)
        nr=NRLIST(no_nr)
        njj=1
        nj=NJ_LOC(NJL_FIEL,njj,nr)
C       Check each face to see if derivatives are continuous across it.
        DO no_nf=1,NFFACE(0,nr)
          nf=NFFACE(no_nf,nr)
          IF(NPF(5,nf).EQ.2) THEN !2 adjacent elements
            NEA(1)=NPF(6,nf) !adjacent element numbers
            NEA(2)=NPF(7,nf)
            NF_EA(1)=NPF(8,nf) !local face number in element nea(1)
            NF_EA(2)=NPF(9,nf) !local face number in element nea(2)
            CONTINUOUS=.TRUE.
            nb_f=NBHF(nh,nc,nf)
            DO neax=1,2 !loop over adjacent elements
              NB_EA(neax)=NBH(nh,nc,NEA(neax))
            ENDDO !neax
C           Check that the same versions are used at each node.
            DO nn_f=1,NNT(nb_f)
              DO neax=1,2 !loop over adjacent elements
                nb_he=NB_EA(neax)
                nn_he=NNF(1+nn_f,NF_EA(neax),nb_he)
                NVA(neax)=NVHE(nn_he,nb_he,nh,NEA(neax))
              ENDDO !neax
              IF(NVA(1).NE.NVA(2)) CONTINUOUS=.FALSE.
            ENDDO !nn_f
            IF(.NOT.CONTINUOUS) THEN
C             Field is zero on this face
              NFLIST(nf)=0
              nb_f=NBJF(nj,nf)
C             Field is zero on each adjacent line
              DO naf=1,NLE(nb_f)
                NLLIST(NLF(naf,nf))=0
              ENDDO !naf
C             At each node field is not default
              ne=NEA(1)
              nb_je=NBJ(nj,ne)
              nf_e=NF_EA(1)
              DO nn_f=1,NNT(nb_f)
                nn_je=NNF(1+nn_f,nf_e,nb_je)
                np=NPNE(nn_je,nb_je,ne)
                NVJP(nj,np)=0
              ENDDO !nn_f
            ENDIF !continuous
          ENDIF !2 adj. elements
        ENDDO !nf
C       Check each line to see if derivatives are continuous across it.
        DO no_nl=1,NLLINE(0,nr)
          nl=NLLINE(no_nl,nr)
          IF(NLLIST(nl).NE.0) THEN !field is still default
            CONTINUOUS=.TRUE.
            FIRST=.TRUE.
            MULTNV=.TRUE. !assume more than one version until checked
            NELT=NEL(0,nl)
            neax=0
C           Loop over elements around the line
            DO WHILE(neax.LT.NELT.AND.CONTINUOUS.AND.MULTNV)
              neax=neax+1
              ne=NEL(neax,nl)
              nb_he=NBH(nh,nc,ne)
C             Loop over local lines in the element
              NLET=NLE(nb_he)
              nae=0
              DO WHILE(nae.LT.NLET.AND.CONTINUOUS)
                nae=nae+1
                IF(NLL(nae,ne).EQ.nl) THEN !same line
C                 Check that the same versions are used in each element.
                  IF(FIRST) THEN
                    IF(NELT.GT.4) THEN !at apex
                      CONTINUOUS=.FALSE.
                    ELSE
                      FIRST=.FALSE.
                      MULTNV=.FALSE. !assume 1 version until more than 1 found
                      NNLT=NNL(0,nae,nb_he)
                      DO nn_l=1,NNLT
                        nn_he=NNL(nn_l,nae,nb_he)
                        NVL(nn_l)=NVHE(nn_he,nb_he,nh,ne)
                        np=NPNE(nn_he,nb_he,ne)
                        IF(NVHP(nh,np,nc,nr).GT.1) MULTNV=.TRUE.
                      ENDDO
                      DO nn_l=NNLT+1,4
                        NVL(nn_l)=0
                      ENDDO
                    ENDIF !NELT
                  ELSE
                    DO nn_l=1,NNL(0,nae,nb_he)
                      nn_he=NNL(nn_l,nae,nb_he)
                      IF(NVHE(nn_he,nb_he,nh,ne).NE.NVL(nn_l))
     '                  CONTINUOUS=.FALSE.
                    ENDDO !nn_l
                  ENDIF
                ENDIF !same line
              ENDDO !nae
            ENDDO !neax
            IF(.NOT.CONTINUOUS) THEN
C             Field is zero on this line
              NLLIST(nl)=0
C             At each node field is not default
              nb_je=NBJ(nj,ne)
              DO nn_l=1,NNL(0,nae,nb_je) !assume consistent arc numbering
                nn_je=NNL(nn_l,nae,nb_je)
                np=NPNE(nn_je,nb_je,ne)
                NVJP(nj,np)=0
              ENDDO !nn_f
            ENDIF !continuous
          ENDIF !NLLIST(nl).NE.0
        ENDDO !nl
C       Check each node to see if there are multiple versions or
C       an essential b.c.
        DO no_np=1,NPNODE(0,nr)
          np=NPNODE(no_np,nr)
          IF(NVJP(nj,np).NE.0) THEN !field is still default
            IF(NVHP(nh,np,nc,nr).GT.1) THEN !>1 version
C             Field is zero at this node
              NVJP(nj,np)=0 !not default
            ELSE
C             Check 0th deriv for an essential b.c.
              ny=NYNP(1,1,nh,np,0,nc,nr)
              IF(ny.NE.0) THEN
                IF(FIX(ny,1)) NVJP(nj,np)=0 !field is zero at this node
              ENDIF
            ENDIF !versions
          ENDIF !default
        ENDDO !np
      ENDDO !nr
C     If there is no face or line adjacent to a non-default node then
C     assume associated derivatives are zero.
      NFLIST(0)=0
      NLLIST(0)=0

C     Loop over elements to set NVJE and the non-default field parameters.
      INCREASE_NVM=.FALSE.
      DO no_nr=1,NRLIST(0)
        nr=NRLIST(no_nr)
        njj=1
        nj=NJ_LOC(NJL_FIEL,njj,nr)
        DO no_ne=1,NEELEM(0,nr)
          ne=NEELEM(no_ne,nr)
          nb_je=NBJ(nj,ne)
          NITB=NIT(nb_je)
          DO nn_je=1,NNT(nb_je)
            np=NPNE(nn_je,nb_je,ne)
            IF(NVJP(nj,np).EQ.-1) THEN !default
              NVJE(nn_je,nb_je,nj,ne)=1
            ELSE IF(NITB.EQ.3) THEN
C             Find local node number in dep var basis
              nb_he=NBH(nh,nc,ne)
              NNTB=NNT(nb_he)
              FOUND=.FALSE.
              nn_he=0
              DO WHILE(nn_he.LT.NNTB.AND..NOT.FOUND)
                nn_he=nn_he+1
                FOUND=NPNE(nn_he,nb_he,ne).EQ.np
              ENDDO
C             Find out if the node is at an essential b.c.
              IF(FOUND) THEN
                nv=NVHE(nn_he,nb_he,nh,ne)
                ny=NYNP(1,nv,nh,np,0,nc,nr)
                FIXED=FIX(ny,1)
              ELSE
                FIXED=.FALSE.
              ENDIF
C             Initialize temporary XP.
              DO nk=1,NKJ(nj,np)
                XP_TEMP(nk)=0.0d0
              ENDDO !nk
C             Find out which adjacent lines and faces are zero
C             and the signs for the non-zero values.
C             Count the zero lines and set the first derivs.
              NZ_LINE=0
              DO ni=1,NITB
                ni2=1+MOD(ni,3)
                ni3=1+MOD(ni+1,3)
                IF(INP(nn_je,ni,nb_je).EQ.1) THEN
                  nf_e=2*ni-1 !adjacent face local number
                  SIGN(ni)=1 !indicates xi(ni) goes into element
                ELSE
                  nf_e=2*ni !adjacent face local number
                  SIGN(ni)=-1 !indicates xi(ni) goes out of element
                ENDIF
                FACE_ZERO(ni)=NFLIST(NFF(nf_e,ne)).EQ.0
C               Adjacent line local number (only works for standard bases)
                nae=4*ni
                IF(INP(nn_je,ni2,nb_je).EQ.1) nae=nae-1
                IF(INP(nn_je,ni3,nb_je).EQ.1) nae=nae-2
                LINE_ZERO(ni)=NLLIST(NLL(nae,ne)).EQ.0
                IF(.NOT.LINE_ZERO(ni)) THEN
C                 Find the derivative orders.
                  IDOI(ni)=2
                  IDOI(ni2)=1
                  IDOI(ni3)=1
                  IF(FIXED) THEN
C                   Check boundary conditions on derivative
                    nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_he,nb_he)
                    nk=NKHE(nk_e,nn_he,nh,ne)
                    ny=NYNP(nk,nv,nh,np,0,nc,nr)
                    LINE_ZERO(ni)=FIX(ny,1)
                  ENDIF !FIXED
                ENDIF !not LINE_ZERO
                IF(LINE_ZERO(ni)) THEN !field is zero on line
                  NZ_LINE=NZ_LINE+1
                ELSE
C                 Set the first deriv.
                  nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_je,nb_je)
                  nk=NKJE(nk_e,nn_je,nj,ne)
                  XP_TEMP(nk)=DBLE(SIGN(ni)*3)
                ENDIF !LINE_ZERO
              ENDDO !ni
C             Count the zero faces and set the face cross derivs.
              NZ_FACE=0
              DO ni=1,NITB
                IF(.NOT.FACE_ZERO(ni)) THEN
C                 Find the derivative orders.
                  ni2=1+MOD(ni,3)
                  ni3=1+MOD(ni+1,3)
                  IDOI(ni)=1
                  IDOI(ni2)=2
                  IDOI(ni3)=2
                  IF(FIXED.AND.LINE_ZERO(ni2).AND.LINE_ZERO(ni3)) THEN
C                   Check boundary conditions on cross derivative
                    nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_he,nb_he)
                    nk=NKHE(nk_e,nn_he,nh,ne)
                    ny=NYNP(nk,nv,nh,np,0,nc,nr)
                    FACE_ZERO(ni)=FIX(ny,1)
                  ENDIF !FIXED
                ENDIF !not FACE_ZERO
C               Set the face cross derivs.
                IF(FACE_ZERO(ni)) THEN !field is zero on face
                  NZ_FACE=NZ_FACE+1
                ELSE IF(LINE_ZERO(ni2).EQV.LINE_ZERO(ni3)) THEN
                  nk_e=NKB(IDOI(1),IDOI(2),IDOI(3),nn_je,nb_je)
                  nk=NKJE(nk_e,nn_je,nj,ne)
                  IF(LINE_ZERO(ni2)) THEN
                    FACTOR=9
                  ELSE
                    FACTOR=-9
                  ENDIF !LINE_ZERO
                  XP_TEMP(nk)=DBLE(SIGN(ni2)*SIGN(ni3)*FACTOR)
                ENDIF !FACE_ZERO
              ENDDO !ni
C             Set the volume cross deriv.
C             I don't know why this formula works but it is much simpler
C             than considering each case in turn.
              FACTOR=NZ_FACE+1-NZ_LINE
              IF(FACTOR.NE.0) THEN
                nk_e=NKB(2,2,2,nn_je,nb_je)
                nk=NKJE(nk_e,nn_je,nj,ne)
                XP_TEMP(nk)=DBLE(SIGN(1)*SIGN(2)*SIGN(3)*27*FACTOR)
              ENDIF !non-zero volume cross deriv
C             Search for a version at np with the appropriate values.
              FOUND=.FALSE.
              NVTOT=NVJP(nj,np)
              NKTOT=NKJ(nj,np)
              nv=0
              DO WHILE(nv.LT.NVTOT.AND..NOT.FOUND)
                nv=nv+1
                FOUND=.TRUE.
                nk=0
                DO WHILE(nk.LT.NKTOT.AND.FOUND)
                  nk=nk+1
                  IF(XP_TEMP(nk).NE.XP(nk,nv,nj,np)) FOUND=.FALSE.
                ENDDO !nk
              ENDDO !nv
              IF(.NOT.FOUND) THEN
C               Add the version
                IF(nv.LT.NVM) THEN
                  nv=nv+1
                  NVJP(nj,np)=nv
                  DO nk=1,NKJ(nj,np)
                    XP(nk,nv,nj,np)=XP_TEMP(nk)
                  ENDDO !nk
                ELSE
                  INCREASE_NVM=.TRUE.
                ENDIF
              ENDIF !not found
              NVJE(nn_je,nb_je,nj,ne)=nv
            ENDIF !default
          ENDDO !nn_je
        ENDDO !ne
      ENDDO !nr

C     Set the default nodal values to 1.
      DO no_nr=1,NRLIST(0)
        DO no_np=1,NPNODE(0,nr)
          np=NPNODE(no_np,nr)
          IF(NVJP(nj,np).EQ.-1) THEN !default
            NVJP(nj,np)=1
            XP(1,1,nj,np)=1.0d0
            DO nk=2,NKJ(nj,np)
              XP(nk,1,nj,np)=0.0d0
            ENDDO
          ENDIF !default
        ENDDO !np
      ENDDO !nr

      CALL ASSERT(.NOT.INCREASE_NVM,'>>Increase NVM',ERROR,*9999)

      CALL EXITS('UPFIEL_EIK')
      RETURN
 9999 CALL ERRORS('UPFIEL_EIK',ERROR)
      CALL EXITS('UPFIEL_EIK')
      RETURN 1
      END


