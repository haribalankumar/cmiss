      SUBROUTINE MECH_FLOW_DYNAM(IBT,IDO,INP,NAN,NBH,NBJ,
     '  NEELEM,NENQ,NEP,ne_flow,NHE,NKHE,NKJE,NITOT,
     '  NPF,NPNE,NQNE,NQS,NQXI,nr_mech,
     '  NVHE,NVJE,
     '  NW,nx_flow,nx_mech,NYNQ,
     '  CE,CG,CGE,CP,CQ,CURVCORRECT,FEXT,PG,REF,RG,SE,XA,XE,
     '  XG,XIP,XP,XQ,YG,YQ,ZA,ZE,ZG,ZP,ERROR,*)

C#### Subroutine: MECH_FLOW_DYNAM
C###  Description:
C###    MECH_FLOW_DYNAM is the dynamic subroutine associted with
C###    MECH_FLOW
      IMPLICIT NONE

      INCLUDE 'coup00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'time02.cmn'
C      INCLUDE 'cmiss$reference:tol00.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  NAN(NIM,NAM,NBFM),NBH(NHM,NCM,NEM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENQ(0:8,NQM),NEP(NPM),NHE(NEM),
     '  NKHE(NKM,NNM,NHM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NQS(NEQM),NQXI(0:NIM,NQSCM),NQNE(NEQM,NQEM),
     '  nr_mech,NVHE(NNM,NBFM,NHM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NW(NEM,3),nx_flow,nx_mech,NYNQ(NHM,NQM,0:NRCM,NXM)
      REAL*8 CE(NMM,NEM),CG(NMM,NGM),CGE(NNM,NGM,NEM),
     '  CP(NMM,NPM),CQ(NMM,NQM),
     '  CURVCORRECT(2,2,NNM,NEM),FEXT(NIFEXTM,NGM,NEM),
     '  PG(NSM,NUM,NGM,NBM),RG(NGM),SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XIP(NIM,NPM),
     '  XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM),
     '  YG(NIYGM,NGM,NEM),YQ(NYQM,NIQM,NAM,NXM),
     '  ZA(NAM,NHM,NCM,NEM),ZE(NSM,NHM),ZG(NHM,NUM),
     '  ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,j,m,n,nb,nb1,nb_stress,ne,ne_curent,
     '  ne_flow,ne_loop,ne_mech,ne_np1,ne_np2,ni,NITOT,
     '  nj,nonq,np1,np2,nq,NUM_NQ,ny_art,ny_lambda_art,
     '  ny_lambda_vien,ny_vien,SCHEME
      REAL*8 ART_FRAC,dydNu(3,3),dydx(3,3),MAT_VECTOR(3,3),
     '  CART_NP1(3),CART_NP2(3),COR(3,3),DEF_LEN,DEF_UNLEN,
     '  DIR_DEFORM(3),DIR_UN_DEFORM(3),DXIZN(3,3),DZNXI(3,3),
     '   LAMBDA,NODE1(3),LENGTH,NODE2(3),NODE1_UNDEF(3),
     '  NODE2_UNDEF(3),PHI(3),POINT(3),
     '  PST(3),PXI,REF(3),RG2D,RGZ,RGZ2D,RM(3,3),T1,T2,TC(3,3),SUM,
     '  SUM_XI_NE1,SUM_XI_NE2,
     '  TG(3,3),TG_COR(3,3),TN(3,3),TNA,UND_CART_NP1(3),UND_CART_NP2(3),
     '  XI(3),XI_NE_NP1(3),XI_NE_NP2(3),XI_POS
      CHARACTER STRESSTYPE*17
      LOGICAL FOUND
      EXTERNAL PXI

      DATA STRESSTYPE/'Total'/

      CALL ENTERS('MECH_FLOW_DYNAM',*9999)

        nb=NBJ(1,ne_flow)

        np1=NPNE(1,nb,ne_flow)
        np2=NPNE(NNT(nb),nb,ne_flow)

        ne_np1=NEP(np1)
        ne_np2=NEP(np2)

        CALL ZPZE(NBH(1,1,ne_np1),1,NHE(ne_np1),NKHE(1,1,1,ne_np1),
     '    NPF(1,1),NPNE(1,1,ne_np1),nr_mech,NVHE(1,1,1,ne_np1),
     '    NW(ne_np1,1),nx_mech,CURVCORRECT(1,1,1,ne_np1),
     '    SE(1,1,ne_np1),ZA(1,1,1,ne_np1),ZE,ZP,ERROR,*9999)


        CALL XPXE(NBJ(1,ne_np1),
     '    NKJE(1,1,1,ne_np1),NPF(1,1),
     '    NPNE(1,1,ne_np1),
     '    nr_mech,NVJE(1,1,1,ne_np1),
     '    SE(1,1,ne_np1),XA(1,1,ne_np1),XE,XP,ERROR,*9999)

C calcultes ZE for host element of node one in host element

        DO nj=1,NITOT
          nb1=NBJ(nj,ne_np1)
          NODE1(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '      INP(1,1,nb1),nb1,1,XIP(1,np1),ZE(1,nj))
          NODE1_UNDEF(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '      INP(1,1,nb1),nb1,1,XIP(1,np1),XE(1,nj))
        ENDDO !nj

C interpolates deformed and undeformed coordinates from
C xi positions in the host mesh

        CALL ZPZE(NBH(1,1,ne_np2),1,NHE(ne_np2),NKHE(1,1,1,ne_np2),
     '    NPF(1,1),NPNE(1,1,ne_np2),nr_mech,NVHE(1,1,1,ne_np2),
     '    NW(ne_np2,1),nx_mech,CURVCORRECT(1,1,1,ne_np2),
     '    SE(1,1,ne_np2),ZA(1,1,1,ne_np2),ZE,ZP,ERROR,*9999)


        CALL XPXE(NBJ(1,ne_np2),
     '    NKJE(1,1,1,ne_np2),NPF(1,1),
     '    NPNE(1,1,ne_np2),
     '    nr_mech,NVJE(1,1,1,ne_np2),
     '    SE(1,1,ne_np2),XA(1,1,ne_np2),XE,XP,ERROR,*9999)

C calcultes ZE for host element of node two in host element

        DO nj=1,NITOT
          nb1=NBJ(nj,ne_np2)
          NODE2(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '      INP(1,1,nb1),nb1,1,XIP(1,np2),ZE(1,nj))
          NODE2_UNDEF(nj)=PXI(IBT(1,1,nb1),IDO(1,1,0,nb1),
     '      INP(1,1,nb1),nb1,1,XIP(1,np2),XE(1,nj))
        ENDDO !nj

C interpolates deformed coordinates from xi positions in deform
C host mesh

        CALL COORD(ITYP10(nr_mech),1,NODE1,CART_NP1,ERROR,*9999)
        CALL COORD(ITYP10(nr_mech),1,NODE2,CART_NP2,ERROR,*9999)
        CALL COORD(ITYP10(nr_mech),1,NODE1_UNDEF,
     '    UND_CART_NP1,ERROR,*9999)
        CALL COORD(ITYP10(nr_mech),1,NODE2_UNDEF,
     '    UND_CART_NP2,ERROR,*9999)

C calculates the cartesian coordinates each node contained within the
C element

        DO nj=1,NITOT
          DIR_DEFORM(nj)=CART_NP1(nj)-CART_NP2(nj)
          DIR_UN_DEFORM(nj)=UND_CART_NP1(nj)-UND_CART_NP2(nj)
        ENDDO !nj

C calculates the direction of the coronary element from its
C cartesian coords in the deformed host mesh


C calculates the extention ratio
        DEF_LEN=0.0d0
        DEF_UNLEN=0.0d0

        DO nj=1,NITOT
          DEF_LEN=DEF_LEN+(DIR_DEFORM(nj)**2.0d0)
          DEF_UNLEN=DEF_UNLEN+(DIR_UN_DEFORM(nj)**2.0d0)
        ENDDO !nj

        DEF_LEN=DEF_LEN**0.5d0
        DEF_UNLEN=DEF_UNLEN**0.5d0
        LAMBDA=DEF_LEN/DEF_UNLEN

        DO nj=1,NITOT
          COR(nj,1)=DIR_DEFORM(nj)
        ENDDO !nj

C calculate coronary cordinate system from deformed position of
C coronary element

        SCHEME=NQS(ne_flow)
        NUM_NQ=NQXI(1,SCHEME) !number of grid points in element

        DO nonq=1,NUM_NQ!number of grid points in each coronary element
          nq=NQNE(ne_flow,nonq)
          POINT(1)=XQ(1,nq)
          POINT(2)=XQ(2,nq)
          POINT(3)=XQ(3,nq)

          IF(CALC_GRID_XI) THEN
C calculates the xi value for the gridpoint within the host mesh
            FOUND=.FALSE.
            IF(ne_np1.eq.ne_np2) THEN
C searches for a valid xi postion within the  host element of the
C nodes of the coronary element (where this element is the same
C for both nodes)

              ne_curent=NEP(np1)
              CALL XPXE(NBJ(1,ne_curent),
     '          NKJE(1,1,1,ne_curent),NPF(1,1),
     '          NPNE(1,1,ne_curent),
     '          nr_mech,NVJE(1,1,1,ne_curent),
     '          SE(1,1,ne_curent),XA(1,1,ne_curent),XE,XP,ERROR,*9999)
              XI_POS=(nonq-1.0d0)/(NUM_NQ-1.0d0)
              DO ni=1,3
                XI(ni)=((1.0d0-XI_POS)*XIP(ni,np1))
     '            +(XI_POS*XIP(ni,np2))
              ENDDO

C interplates the inital xi values from nodal values

              ne=0
              CALL DEXI_POINT(IBT,IDO,INP,ne,NBJ,
     '          ne_curent,3,nr_mech,0.d0,XE,XI,XI,POINT,.FALSE.,ERROR,
     '          *9999)

C calls dexi point to iterativly calculate the xi point of the grid
C point position within element ne_curent

              ne=ne_np1
              IF(XI(2).GT.1.0D0) THEN
                XI(2)=1.0D0
              ELSE IF(XI(2).LT.0d0) THEN
                ne=0
              ENDIF
              IF(XI(3).GT.1.0D0) THEN
                XI(3)=1.0D0
              ELSE IF(XI(2).LT.0.0d0) THEN
                ne=0
              ENDIF
              IF(ne.NE.0) THEN
                FOUND=.TRUE.
              ENDIF
            ELSE ! ne_np1.ne.ne_np2
C if both element nodes are not in the same host elemnet try the first
C host element starting at the xi position of that node np1
              ne_curent=NEP(np1)
              CALL XPXE(NBJ(1,ne_curent),
     '          NKJE(1,1,1,ne_curent),NPF(1,1),
     '          NPNE(1,1,ne_curent),
     '          nr_mech,NVJE(1,1,1,ne_curent),
     '          SE(1,1,ne_curent),XA(1,1,ne_curent),XE,XP,ERROR,*9999)
              DO ni=1,3
                XI(ni)=XIP(ni,np1)
              ENDDO
              ne=0
              CALL DEXI_POINT(IBT,IDO,INP,ne,NBJ,
     '          ne_curent,3,nr_mech,0.d0,XE,XI,XI,POINT,.FALSE.,ERROR,
     '          *9999)
C calls dexi point to iterativly calculate the xi point of the grid
C point position within element ne_curent
              IF(ne.NE.0) THEN
                FOUND=.TRUE.
              ELSE
                SUM_XI_NE1=0.0D0
                DO ni=1,3
                  IF(XI(ni).GT.1.0D0) THEN
                    SUM_XI_NE1=SUM_XI_NE1+(XI(ni)-1.0D0)
                  ELSE IF(XI(ni).LT.0.0D0) THEN
                    SUM_XI_NE1=SUM_XI_NE1-XI(ni)
                  ENDIF
                  XI_NE_NP1(ni)=XI(ni)
                ENDDO
              ENDIF
              IF(.NOT.FOUND) THEN
C try the second host element starting at np2
                ne_curent=NEP(np2)
                CALL XPXE(NBJ(1,ne_curent),
     '            NKJE(1,1,1,ne_curent),NPF(1,1),
     '            NPNE(1,1,ne_curent),
     '            nr_mech,NVJE(1,1,1,ne_curent),
     '            SE(1,1,ne_curent),XA(1,1,ne_curent),XE,XP,ERROR,*9999)
                DO ni=1,3
                  XI(ni)=XIP(ni,np2)
                ENDDO !ni
C starts at the xi postion of the element node
                ne=0
                CALL DEXI_POINT(IBT,IDO,INP,ne,NBJ,
     '            ne_curent,3,nr_mech,0.d0,XE,XI,XI,POINT,.FALSE.,
     '            ERROR,*9999)
C calls dexi point to iterativly calculate the xi point of the grid
C point position within element ne_curent
                IF(ne.NE.0) THEN
                  FOUND=.TRUE.
                ELSE
                  SUM_XI_NE2=0.0D0
                  DO ni=1,3
                    IF(XI(ni).GT.1.0D0) THEN
                      SUM_XI_NE2=SUM_XI_NE2+(XI(ni)-1.0D0)
                    ELSE IF(XI(ni).LT.0.0D0) THEN
                     SUM_XI_NE2=SUM_XI_NE2-XI(ni)
                    ENDIF
                    XI_NE_NP2(ni)=XI(ni)
                  ENDDO  !ni
                ENDIF !ne.NE.0
              ENDIF !.NOT.FOUND
            ENDIF ! ne_np1.eq.ne_np2
            IF(.NOT.FOUND) THEN
C if the xi position of the grid point cannot be found in either of
C the host elements of the coronary elements nodes then search though
C all elements
              ne=0
              ne_loop=0
              DO WHILE((ne.EQ.0).AND.(ne_loop.LT.NEELEM(0,nr_mech)))
                ne_loop=ne_loop+1
                ne_curent=NEELEM(ne_loop,nr_mech)
                DO ni=1,3 !initialising XI
                  XI(ni)=0.5d0
                ENDDO
                CALL XPXE(NBJ(1,ne_curent),
     '            NKJE(1,1,1,ne_curent),NPF(1,1),
     '            NPNE(1,1,ne_curent),
     '            nr_mech,NVJE(1,1,1,ne_curent),
     '            SE(1,1,ne_curent),XA(1,1,ne_curent),XE,XP,ERROR,*9999)
                CALL DEXI_POINT(IBT,IDO,INP,ne,NBJ,
     '            ne_curent,3,nr_mech,0.d0,XE,XI,XI,POINT,.FALSE.,
     '            ERROR,*9999)
              ENDDO !((ne.EQ.0).AND.(ne_loop.LT.NEELEM(0,nr_mech)))
            ENDIF !.NOT.FOUND

            IF(NE.EQ.0) THEN
              IF(SUM_XI_NE1.GT.SUM_XI_NE2) THEN
                NE=ne_np2
                XI(1)=XI_NE_NP2(1)
                XI(2)=XI_NE_NP2(2)
                XI(3)=XI_NE_NP2(3)
              ELSE
                NE=ne_np1
                XI(1)=XI_NE_NP1(1)
                XI(2)=XI_NE_NP1(2)
                XI(3)=XI_NE_NP1(3)
              ENDIF
            ENDIF

            CALL ASSERT(ne.NE.0,'>>nq xi value not found',ERROR,*9999)

            YQ(NYNQ(1,nq,0,nx_flow),10,1,nx_flow)=XI(1)
            YQ(NYNQ(2,nq,0,nx_flow),10,1,nx_flow)=XI(2)
            YQ(NYNQ(3,nq,0,nx_flow),10,1,nx_flow)=XI(3)
            CQ(8,nq)=ne
! storese the xi values of a grid point and the host element

          ELSE
            XI(1)=YQ(NYNQ(1,nq,0,nx_flow),10,1,nx_flow)
            XI(2)=YQ(NYNQ(2,nq,0,nx_flow),10,1,nx_flow)
            XI(3)=YQ(NYNQ(3,nq,0,nx_flow),10,1,nx_flow)
            ne=NINT(CQ(8,nq))
C retrieves previously calculated values of the xi values of
C the grid point

          ENDIF

          ne_mech=ne ! the mechanics element in the host mesh
          CALL XPXE(NBJ(1,ne_mech),NKJE(1,1,1,ne_mech),NPF(1,1),
     '      NPNE(1,1,ne_mech),nr_mech,NVJE(1,1,1,ne_mech),
     '      SE(1,1,ne_mech),XA(1,1,ne_mech),XE,XP,ERROR,*9999)

C interpolates the local XE array for the element ne_mech

          IF(STRESS_CALC.EQ.0) THEN
C calculates the stress solution directly from the deformation of the
C host mesh
            CALL ZPZE(NBH(1,1,ne_mech),1,NHE(ne_mech),
     '        NKHE(1,1,1,ne_mech),
     '        NPF(1,1),NPNE(1,1,ne_mech),nr_mech,NVHE(1,1,1,ne_mech),
     '        NW(ne_mech,1),nx_mech,CURVCORRECT(1,1,1,ne_mech),
     '        SE(1,1,ne_mech),ZA(1,1,1,ne_mech),ZE,ZP,ERROR,*9999)
            CALL CPCG(1,NBH(1,1,ne_mech),NPNE(1,1,ne_mech),
     '        nr_mech,nx_mech,
     '        CE(1,ne_mech),CG,CGE(1,1,ne_mech),CP,PG,ERROR,*9999)
            CALL ZETX50('Fibre','Cauchy',STRESSTYPE,IBT,IDO,INP,NAN,
     '        NBH(1,1,ne_mech),NBJ(1,ne_mech),0,NHE(ne_mech),
     '        NPNE(1,1,ne_mech),nr_mech,ne_mech,nx_mech,
     '        CE(1,ne_mech),CG,CP,FEXT(1,1,ne_mech),PG,PHI,PST,
     '        RG(1),RG2D,RGZ,RGZ2D,RM,TC,TG,TN,TNA,XE,XG,XI,
     '        YG(1,1,ne_mech),ZE,ZG,ERROR,*9999)
          ELSE IF(STRESS_CALC.EQ.1) THEN

C interplolates the vessel stresses from nodal fields

            CALL ZPZE(NBH(1,1,ne_mech),1,NHE(ne_mech),
     '        NKHE(1,1,1,ne_mech),
     '        NPF(1,1),NPNE(1,1,ne_mech),nr_mech,NVHE(1,1,1,ne_mech),
     '        NW(ne_mech,1),nx_mech,CURVCORRECT(1,1,1,ne_mech),
     '        SE(1,1,ne_mech),ZA(1,1,1,ne_mech),ZE,ZP,ERROR,*9999)

            nj=NJ_LOC(NJL_FIEL,1,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(1,1)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 11 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the first
C compinent of the nodal field

            nj=NJ_LOC(NJL_FIEL,2,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(2,2)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 22 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the second
C compinent of the nodal field

            nj=NJ_LOC(NJL_FIEL,3,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(3,3)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 33 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the third
C compinent of the nodal field

            nj=NJ_LOC(NJL_FIEL,4,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(1,2)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 12 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the fourth
C compinent of the nodal field

            nj=NJ_LOC(NJL_FIEL,5,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(1,3)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 13 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the fifth
C compinent of the nodal field

            nj=NJ_LOC(NJL_FIEL,6,nr_mech)
            nb_stress=NBJ(nj,ne_mech)
            TG(2,3)=PXI(IBT(1,1,nb_stress),IDO(1,1,0,nb_stress),
     '      INP(1,1,nb_stress),nb_stress,1,XI,XE(1,nj))

C interplolates the 23 component of Cauchy stress tensor
C at the at the xi postion of the grid point using the sixth
C compinent of the nodal field

            TG(2,1)=TG(1,2) ! the stress tensor is symetric
            TG(3,1)=TG(1,3)
            TG(3,2)=TG(2,3)

          ENDIF

          CALL MAT_VEC_XI(IBT,IDO,INP,NAN,NBJ(1,ne_mech),
     '      nr_mech,MAT_VECTOR(1,1),MAT_VECTOR(1,2),MAT_VECTOR(1,3),
     '      XE,XG,XI,.TRUE.,ERROR,*9999)

C finds fibre,sheet etc material axes from undeformed geometry

c          CALL NORMALISE(NITOT,COR(1,1),ERROR,*9999)


          LENGTH=0.0d0
          DO nj=1,NITOT
            LENGTH=LENGTH+(COR(nj,1)**2.0d0)
          ENDDO
          LENGTH=(LENGTH**0.5d0)

          DO nj=1,NITOT
            COR(nj,1)=COR(nj,1)/LENGTH
          ENDDO

          CALL DXIDZM(IBT,IDO,INP,NAN,NBH,NBJ,0,NHE(ne_mech),nr_mech,
     '      nx_mech,DXIZN,DZNXI,PG,XE,XG,XI,ZE,ZG,'Fibre',ERROR,*9999)

          DO nj=1,NITOT
            REF(nj)=0.0d0
          ENDDO
          IF(ABS(COR(1,1)).LT.0.7d0) THEN
C making sure reference direction is not
C coincident with a the vessel dirction
            REF(1)=1.0d0
          ELSE
            REF(2)=1.0d0
          ENDIF

          CALL CROSS(COR(1,1),REF,COR(1,2))
          CALL CROSS(COR(1,1),COR(1,2),COR(1,3))


c          CALL NORMALISE(NITOT,COR(1,2),ERROR,*9999)
c          CALL NORMALISE(NITOT,COR(1,3),ERROR,*9999)

          LENGTH=0.0d0
          DO nj=1,NITOT
            LENGTH=LENGTH+(COR(nj,2)**2.0d0)
          ENDDO
          LENGTH=(LENGTH**0.5d0)
          DO nj=1,NITOT
            COR(nj,2)=COR(nj,2)/LENGTH
          ENDDO

          LENGTH=0.0d0
          DO nj=1,NITOT
            LENGTH=LENGTH+(COR(nj,3)**2.0d0)
          ENDDO
          LENGTH=(LENGTH**0.5d0)
          DO nj=1,NITOT
            COR(nj,3)=COR(nj,3)/LENGTH
          ENDDO

C calculting axes of coronary element

          DO i=1,NITOT
            DO j=1,NITOT
              dydx(i,j)=COR(j,i)
            ENDDO
          ENDDO
C calculating the rotation from reference coords into vessel coords

          CALL CALC_dydNu(IBT,IDO,INP,NAN,NBJ,nr_mech,
     '      dydNu,dydx,XE,XG,XI,ZE,ERROR,*9999)

          DO i=1,NITOT
            DO j=1,NITOT
              SUM=0.0d0
              DO n=1,NITOT
                DO m=1,NITOT
                  SUM=SUM+dydNu(i,n)*dydNu(j,m)*TG(m,n)
                ENDDO
              ENDDO
              TG_COR(i,j)=SUM
            ENDDO
          ENDDO

          T1=0.0d0
          T2=0.0d0
          DO I=1,NITOT
            DO J=1,NITOT
              T1=TG_COR(2,2)
              T2=TG_COR(3,3)
            ENDDO
          ENDDO


          ny_art=NYNQ(1,nq,0,nx_flow)
          ny_vien=NYNQ(4,nq,0,nx_flow)

          ART_FRAC=YQ(NYNQ(2,nq,0,nx_flow),2,1,nx_flow)

C calculating the strech ratio of each grid point
          ny_lambda_art=NYNQ(3,nq,0,nx_flow)
          ny_lambda_vien=NYNQ(6,nq,0,nx_flow)

C stores old value of trace

          IF(NENQ(0,nq).EQ.1) THEN
C the grid point is only in one element thus the full trace is added to
C that grid point

            YQ(ny_art,2,1,nx_flow)=(0.50d0*(T1+T2)*ART_FRAC)
            YQ(ny_lambda_art,2,1,nx_flow)=LAMBDA
            YQ(ny_lambda_vien,2,1,nx_flow)=LAMBDA
            YQ(ny_vien,2,1,nx_flow)=YQ(ny_art,2,1,nx_flow)
          ELSE IF(NENQ(0,nq).EQ.2) THEN
cC$OMP CRITICAL(MECH_FLOW_DYNAM)
C locks the variable to avoid simultaneous writing when
C multi processing
            YQ(ny_art,2,1,nx_flow)=YQ(ny_art,2,1,nx_flow)+
     '        (0.250d0*(T1+T2)*ART_FRAC)
            YQ(ny_lambda_art,2,1,nx_flow)=
     '        YQ(ny_lambda_art,2,1,nx_flow)+
     '        LAMBDA/2.0d0
            YQ(ny_lambda_vien,2,1,nx_flow)=
     '        YQ(ny_lambda_vien,2,1,nx_flow)+
     '        LAMBDA/2.0d0
            YQ(ny_vien,2,1,nx_flow)=YQ(ny_art,2,1,nx_flow)
cC$OMP END CRITICAL(MECH_FLOW_DYNAM)

C takes average of trace from two element and removes half the old
C values of the trace each of the two times
          ELSE
            CALL ASSERT(NENQ(0,nq).LT.2,'>>trifurcation',ERROR,*9999)
          ENDIF

        ENDDO !nonq

      CALL EXITS('MECH_FLOW_DYNAM')
      RETURN
 9999 CALL ERRORS('MECH_FLOW_DYNAM',ERROR)
      CALL EXITS('MECH_FLOW_DYNAM')
      RETURN 1
      END



