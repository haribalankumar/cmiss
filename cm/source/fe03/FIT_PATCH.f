      SUBROUTINE FIT_PATCH(IPIVOT,NELIST,NENP,NKJE,NKH,NPF,NPLIST,
     '  NPLIST2,NPNE,nr,NVHP,NVJE,nxFIT,GK,GR,PG,SE,XA,XE,XG,
     '  XP,YG,ZP,ERROR,*)

C#### Subroutine: FIT_PATCH
C###  Description:
C###    Implementation of O.C. Zienkiewicz and J.Z. Zhu's Super
C###    Convergent Patch Recovery technique from
C###    Int. J. Nurmer. Methods Eng.,33,1331-1364(1992).

C****   Only bilinear quadrilateral, biquadratic quadrilateral
C****   linear triangular and trilinear elements, are implemented.
C****   CS 11/6/98 I am not sure linear triangular elements is correct

C**** Created by Carey Stevens May 1997

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mach00.inc'
      INCLUDE 'mxch.inc'
!     Parameter List
      INTEGER IPIVOT(NOM),
     '  NELIST(0:NEM),NENP(NPM,0:NEPM,0:NRM),NKJE(NKM,NNM,NJM,NEM),
     '  NKH(NHM,NPM,NCM,0:NRM),NPF(9,NFM),NPLIST(0:NPM),NPLIST2(0:NPM),
     '  NPNE(NNM,NBFM,NEM),nr,NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJE(NNM,NBFM,NJM,NEM),nxFIT
      REAL*8 GK(NZ_GK_M),GR(NYROWM),PG(NSM,NUM,NGM,NBM),
     '  SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM),
     '  YG(NIYGM,NGM,NEM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ASSEMBLY_NODE,nb,i,INFO,j,NBJ_temp(12),ne,
     '  ng,nh,nhj,nhx,nj,njj,njj2,nk,nn,no_nelist,no_ng,
     '  no_nplist,np,NUM_NG_NE,NUM_TERMS,nv,nz,PATTERN2(4,4),
     '  PATTERN3(8,8)
      REAL*8 P(9),PP1

C      CS 22/6/98 vars for 2d stuff temporarily removed
C      INTEGER MAX_ELE,ne1,ne2,ne_main
C      REAL*8 Z
C      LOGICAL MID

C     Patterns pick out which Gauss points are to be fitted in the
C     patch. ng=PATTERNX(nn,no_ng) where nn is the local node number
C     of the assembly node in the considered element, and no_ng
C     is a count of the fitted Gauss points
      DATA PATTERN2/1,2,4,5,
     '  2,3,5,6,
     '  4,5,7,8,
     '  5,6,8,9/
      DATA PATTERN3/1,2,4,5,10,11,13,14,
     '  2,3,5,6,11,12,14,15,
     '  4,5,7,8,13,14,16,17,
     '  5,6,8,9,14,15,17,18,
     '  19,20,22,23,10,11,13,14,
     '  20,21,23,24,11,12,14,15,
     '  22,23,25,26,13,14,16,17,
     '  23,24,26,27,14,15,17,18/

      CALL ENTERS('FIT_PATCH',*9999)

C     Initialise NPLIST2
      DO i=1,NPM
        NPLIST2(i)=0
      ENDDO

      nb=PATCH_BASIS ! specified in IPFIT()

      IF(NIT(nb).EQ.2) THEN
        IF(PATCH_POLYNOMIAL.EQ.1) THEN
          NUM_TERMS=3
        ELSE
          NUM_TERMS=NNT(nb)
        ENDIF
      ELSE
        NUM_TERMS=NNT(nb)
      ENDIF

      IF(DOP) THEN
        WRITE(OP_STRING,
     '    '('' Polynomial terms ='',I3,'' Gauss point basis ='',I3)')
     '    NUM_TERMS,nb
        CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
      ENDIF

C**** Should initialise XP and ZP for iterative calling

      DO njj=1,3 !geom/fibres/field
        DO njj2=1,NJ_LOC(njj,0,nr)
          nj=NJ_LOC(njj,njj2,nr)
          CALL ASSERT(nj.LE.12,'>>ERROR: increase size of '
     '      //'NBJ_temp to NJM',ERROR,*9999)
          NBJ_temp(nj)=nb
        ENDDO !njj2
      ENDDO !njj
      IF(PATCH_PATTERN.EQ.1) THEN ! determine # of gauss points per elem
        NUM_NG_NE=NGT(nb)
      ELSE IF(PATCH_PATTERN.EQ.2) THEN
        NUM_NG_NE=4
      ELSE IF(PATCH_PATTERN.EQ.3) THEN
        NUM_NG_NE=8
      ENDIF

C     Loop over each recovery node
      DO no_nplist=1,NPLIST(0)
        np=NPLIST(no_nplist)

C       Determine the patch assembly node
        NELIST(0)=NENP(np,0,nr)
        ne=NENP(np,1,nr)
C CS not required while 2d stuff removed
C        ne_main=ne ! first connected element

        IF(PATCH_INTERNAL) THEN
C         Assemble about internal nodes if np is on the boundary
C         This is a cumbersome and confusing way to do this but ...
          !****** removed for the moment CS 11/6/98
          !****** iplementation being rethought
          CALL ASSERT(.FALSE.,'>>Internal patches for boundary'
     '      //' nodes is unavailable at the moment',ERROR,*9999)
        ELSE
          ASSEMBLY_NODE=np
        ENDIF

        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Recovery node ='',I3,'' Assembly node ='',I3)')
     '      np,ASSEMBLY_NODE
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C       Determine which elements are in the patch
        NELIST(0)=NENP(ASSEMBLY_NODE,0,nr)
        DO no_nelist=1,NELIST(0)
          NELIST(no_nelist)=NENP(ASSEMBLY_NODE,no_nelist,nr)
        ENDDO

        CALL ASSERT((NELIST(0)*NUM_NG_NE).GE.NUM_TERMS,'>>Insufficient'
     '    //' data points to fit polynomial',ERROR,*9999)

        IF(DOP) THEN
          WRITE(OP_STRING,
     '      '('' Elements in patch : '',4I6)')
     '      (NELIST(no_nelist),no_nelist=1,NELIST(0))
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

C       Fit
        DO njj=1,NJ_LOC(NJL_FIEL,0,nr)

          IF(njj.EQ.1) THEN !Initialise GK & GR
            DO nz=1,NUM_TERMS*NUM_TERMS
              GK(nz)=0.0d0
            ENDDO !nz
            DO nz=1,NUM_TERMS
              GR(nz)=0.0d0
            ENDDO !nz
          ELSE
            DO nz=1,NUM_TERMS
              GR(nz)=0.0d0
            ENDDO !nz
          ENDIF

C         Calculate Matrix and RHS
          DO no_nelist=1,NELIST(0)
            ne=NELIST(no_nelist)
            CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)

            DO i=1,NNT(nb)
              IF(NPNE(i,nb,ne).EQ.ASSEMBLY_NODE) nn=i
            ENDDO

            DO no_ng=1,NUM_NG_NE
              IF(PATCH_PATTERN.EQ.1) THEN
                ng=no_ng
              ELSE IF(PATCH_PATTERN.EQ.2) THEN
                ng=PATTERN2(nn,no_ng)
              ELSE IF(PATCH_PATTERN.EQ.3) THEN
                ng=PATTERN3(nn,no_ng)
              ENDIF

              IF(njj.EQ.1) THEN
                IF(DOP) THEN
                  WRITE(OP_STRING,
     '              '('' Gauss point number '',I2,'' in element '',I6)')
     '              ng,ne
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF

              CALL XEXG(NBJ_temp,ng,nr,PG,XE,XG,ERROR,*9999)
              IF(NNT(nb).EQ.3) THEN
                P(1)=1
                P(2)=XG(1,1)
                P(3)=XG(2,1)
              ELSE IF(NNT(nb).EQ.4) THEN
                P(1)=1
                P(2)=XG(1,1)
                P(3)=XG(2,1)
                IF(NUM_TERMS.EQ.4) THEN
                  P(4)=XG(1,1)*XG(2,1)
                ENDIF
              ELSE IF(NNT(nb).EQ.6) THEN
                P(1)=1
                P(2)=XG(1,1)
                P(3)=XG(2,1)
                P(4)=XG(1,1)*XG(2,1)
                P(5)=XG(1,1)*XG(1,1)
                P(6)=XG(2,1)*XG(2,1)
              ELSE IF(NNT(nb).EQ.8) THEN
                P(1)=1
                P(2)=XG(1,1)
                P(3)=XG(2,1)
                P(4)=XG(3,1)
                P(5)=XG(1,1)*XG(2,1)
                P(6)=XG(1,1)*XG(3,1)
                P(7)=XG(2,1)*XG(3,1)
                P(8)=XG(1,1)*XG(2,1)*XG(3,1)
              ELSE IF(NNT(nb).EQ.9) THEN
                P(1)=1
                P(2)=XG(1,1)
                P(3)=XG(2,1)
                P(4)=XG(1,1)*XG(2,1)
                P(5)=XG(1,1)*XG(1,1)
                P(6)=XG(2,1)*XG(2,1)
                P(7)=XG(1,1)*XG(1,1)*XG(2,1)
                P(8)=XG(1,1)*XG(2,1)*XG(2,1)
                P(9)=XG(1,1)*XG(1,1)*XG(2,1)*XG(2,1)
              ENDIF
              IF(njj.EQ.1) THEN
                DO i=1,NUM_TERMS
                  GR(i)=GR(i)+P(i)*YG(NG_FIT(1,njj),ng,ne)
                  DO j=1,NUM_TERMS
                    GK(j+(i-1)*NUM_TERMS)=GK(j+(i-1)*NUM_TERMS)+P(i)
     '                *P(j)
                  ENDDO
                ENDDO
              ELSE
                DO i=1,NUM_TERMS
                  GR(i)=GR(i)+P(i)*YG(NG_FIT(1,njj),ng,ne)
                ENDDO
              ENDIF
            ENDDO !ng
          ENDDO !no_nelist

C         Factorise the system of linear equations
          IF(njj.EQ.1) THEN !factorise GK
C           Find the LU factorisation of GK.
C           GK is overwritten and the
C           pivots are stored in pivot
            CALL DGETRF(NUM_TERMS,NUM_TERMS,GK,NUM_TERMS,IPIVOT,INFO)
            IF(INFO.NE.0) THEN
              WRITE(ERROR,'('' >>INFO='',I6,'' in  DGETRF'')') INFO
              GOTO 9999
            ENDIF
          ENDIF

C         Solve the factorised system
          CALL DGETRS('N',NUM_TERMS,1,GK,NUM_TERMS,IPIVOT,GR,
     '      NUM_TERMS,INFO)
          IF(INFO.NE.0) THEN
            WRITE(ERROR,'('' >>INFO='',I6,'' in  DGETRS'')') INFO
            GOTO 9999
          ENDIF

C         Calculate nodal value
          CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
          DO nhj=1,NUM_FIT(njj)
            nj=NLH_FIT(nhj,1,njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nxFIT)
            DO nv=1,NVHP(nh,np,1,nr)
              DO nk=1,NKH(nh,np,1,nr)
                ZP(nk,nv,nh,np,1)=PP1(NUM_TERMS,GR,XP(1,1,1,np),
     '            XP(1,1,2,np),XP(1,1,3,np))
                XP(nk,nv,nj,np)=ZP(nk,nv,nh,np,1)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nhj
        ENDDO !njj
      ENDDO !no_nplist

      CALL EXITS('FIT_PATCH')
      RETURN
 9999 CALL ERRORS('FIT_PATCH',ERROR)
      CALL EXITS('FIT_PATCH')
      RETURN 1
      END


