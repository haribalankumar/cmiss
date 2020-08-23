      SUBROUTINE UPMESH_GEOMETRY(NBJ,NEELEM,NELIST,ne_start,NORD,NPNE,
     &  NVJE,NXI,AA_R,ACINUS_LENGTH,ACINUS_VOLUME,CE,ld_ratio,
     &  RATIO_DIAMETER,scale_f,value,XAB,XP,CUBE_ROOT,UNSTRAINED,
     &  RADIUS_SCHEME,SCALED,ERROR,*)

C####  Subroutine: UPMESH_GEOMETRY
C###   Description:
C###     UPMESH_GEOMETRY updates the fields associated with 1D mesh
C###     geometry. i.e. the radius and ratio of alveolar cross-sectional
C###     area. The update field command is not used at this stage
C###     (though this will hopefully be merged in the future) because 
C###     1. UPFG can't currently update the correct nodal versions, and
C###     2. we don't have a standard format for order based data.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      INCLUDE 'lung_nej00.cmn'

!     Parameter values
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),NELIST(0:NEM),ne_start,
     &  NORD(5,NE_R_M),NPNE(NNM,NBFM,NEM),NVJE(NNM,NBFM,NJM,NEM),
     &  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 AA_R,ACINUS_LENGTH,ACINUS_VOLUME,CE(NMM,NEM),ld_ratio,
     &  RATIO_DIAMETER,scale_f,value,XAB(NORM,NEM),XP(NKM,NVM,NJM,NPM)
      LOGICAL CUBE_ROOT,UNSTRAINED,SCALED
      CHARACTER RADIUS_SCHEME*255,ERROR*(*)
!     Local variables
      INTEGER N,nb,ne,ne0,ne1,ne2,NGEN_OFFSET,n_generation,
     &  n_generation_ne0,n_horsfield,n_horsfield_p,n_horsfield_s,nindex,
     &  nj,NMAX_ORD,nn,noelem,np,np0,np1,np2,np3,n_order,n_strahler,
     &  NTALLY(GENM),nv,nv1,nv2
      REAL*8 AA_RATIO,angle1,angle2,BRANCH_D,BRANCH_L,compliance_vessel,
     &  diameter,length,length2,lengths(NE_R_M),Ptm,
     &  radius,R_SQ,RatioDmajorToDparent,
     &  RatioLengthToDminor,RatioLengthToDmajor,RatioDminorToDmajor,
     &  RatioDminorToDparent,slope,unstrained_factor,vector0(3),
     &  vector1(3),vector2(3),weight1,weight2,X(GENM),YREGRESS(GENM)
      REAL*8 DOT_PROD
      LOGICAL CHECK,MAJOR_CHILD,TERMINAL
C      INTEGER N,nb,ne,ne0,NGEN_OFFSET,n_generation,n_generation_ne0,
C     &  n_horsfield,nindex,nj,NMAX_ORD,nn,noelem,np,np1,np2,n_order,
C     &  n_strahler,NTALLY(GENM),nv,nv1,nv2
C      REAL*8 AA_RATIO,BRANCH_D,compliance_vessel,diameter,length,
C     &  lengths(NE_R_M),Ptm,radius,R_SQ,slope,unstrained_factor,
C     &  X(GENM),YREGRESS(GENM)

      CALL ENTERS('UPMESH_GEOMETRY',*9999)
      
      IF(nj_radius.EQ.0.AND.nj_alveoli.EQ.0) THEN
        CALL ASSERT(.FALSE.,
     &    '>> Neither radius or alveolar fields set up, use: '
     &    // 'FEM DEFINE MESH;C LUNG FIELD',ERROR,*9999)
      ENDIF
      
      IF(RADIUS_SCHEME.EQ.'CONSTANT')THEN
C     The radius of all nodes in the element list is a constant.
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(nj_radius,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            CALL ASSERT(nv.LE.NVM,'>>Increase NVM',ERROR,*9999)
            XP(1,nv,nj_radius,np)=value*scale_f
          ENDDO !nn
        ENDDO !noelem
        
      ELSE IF(RADIUS_SCHEME(1:8).EQ.'LD_RATIO')THEN
C     The radius is equal to the TOTAL branch length * L:D ratio

C       Calculate the length of each element in the list        
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          n_generation=NORD(1,ne)
          nb=NBJ(nj_radius,ne)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          nv1=NVJE(1,nb,1,ne)
          nv2=NVJE(2,nb,1,ne)
          length=0.d0
          DO nj=1,3
            length=length+(XP(1,nv1,nj,np1)-XP(1,nv2,nj,np2))**2
          ENDDO !nj
          lengths(ne)=DSQRT(length)
        ENDDO !noelem

C       Calculate the total branch length for each element 
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          n_generation=NORD(1,ne) !generation
C         Add preceeding element lengths in same generation
          length=lengths(ne) !element length
          ne0=NXI(-1,1,ne) !parent element
          IF(ne0.NE.0)THEN !not stem
            n_generation_ne0=NORD(1,ne0) !parent generation
            DO WHILE(n_generation_ne0.EQ.n_generation) !same generation
              length=length+lengths(ne0) !add parent length
              ne0=NXI(-1,1,ne0) !next parent
              IF(ne0.NE.0)THEN
                n_generation_ne0=NORD(1,ne0)
              ELSE
                n_generation_ne0=n_generation-1
              ENDIF
            ENDDO !while of the same generation
          ENDIF !not stem
C         Add subtended element lengths in same generation
          ne0=NXI(1,1,ne) !child element (only 1 if in same generation)
          IF(ne0.NE.0)THEN !not terminal
            n_generation_ne0=NORD(1,ne0) !child generation
            DO WHILE(n_generation_ne0.EQ.n_generation) !same generation
              length=length+lengths(ne0) !add child length
              ne0=NXI(1,1,ne0) !next child
              IF(ne0.NE.0)THEN
                n_generation_ne0=NORD(1,ne0)
              ELSE
                n_generation_ne0=n_generation-1
              ENDIF
            ENDDO !while of the same generation
          ENDIF !not terminal

C         Calculate the radius as R=L/(L:D_ratio)/2
          nb=NBJ(nj_radius,ne) !basis function for radius
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=0.5d0*length/ld_ratio*scale_f
          ENDDO !nn
        ENDDO !noelem

      ELSE IF(RADIUS_SCHEME(1:6).EQ.'WEIBEL')THEN
C     Radii equal to measured values from Weibel(1963), based on
C       generation
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          n_generation=MIN(18,NORD(1,ne))
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=0.5d0*WEIBEL_DIAM(n_generation)
     &        *scale_f
            IF(nj_alveoli.NE.0) THEN
              nv=NVJE(nn,nb,nj_alveoli,ne)
              XP(1,nv,nj_alveoli,np)=1.d0 !conducting, therefore a/A=1
            ENDIF
          ENDDO !nn
        ENDDO !noelem
        
      ELSE IF(RADIUS_SCHEME(1:9).EQ.'HORSFIELD')THEN
C     Radii equal to measured values from Horsfield & Cumming (1967),
C     based on Horsfield order
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          ne0=NXI(-1,1,ne)
          nb=NBJ(1,ne)
          n_horsfield=NORD(2,ne)
          IF(ne0.NE.0)THEN
            np=NPNE(2,nb,ne0) !node at end of parent
            nv=NVJE(2,nb,nj_radius,ne0) !version of end parent node
            IF(XP(1,nv,nj_radius,
     &        np).LT.HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield))THEN
              RADIUS=XP(1,nv,nj_radius,np)
            ELSE
              RADIUS=HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)*0.5d0
            ENDIF
c              RADIUS=HORSFIELD_AIRWAY_DIAM(n_horsfield)*0.5d0
          ELSE
            CALL ASSERT(n_horsfield.LT.30,'>>Too many orders!',ERROR,
     &        *9999)
            RADIUS=HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)*0.5d0
c            RADIUS=HORSFIELD_AIRWAY_DIAM(n_horsfield)*0.5d0
          ENDIF
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=RADIUS*scale_f
            IF(nj_alveoli.NE.0) THEN
              nv=NVJE(nn,nb,nj_alveoli,ne)
              XP(1,nv,nj_alveoli,np)=1.d0 !conducting, therefore a/A=1
            ENDIF
          ENDDO !nn
        ENDDO !noelem

      ELSE IF(RADIUS_SCHEME(1:8).EQ.'COMBINED')THEN
C     Radii based on measured values from Horsfield & Cumming (1967),
C     adapted to take into account imaging based measures of airway or
C     vessel size from medical imaging.
        RatioDmajorToDparent=0.86d0
        RatioLengthToDminor=2.7d0
        RatioLengthToDmajor=3.0d0
        RatioDminorToDmajor=0.8d0
        RatioDminorToDparent=RatioDmajorToDparent*RatioDminorToDmajor
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          ne0=NXI(-1,1,ne)
          IF(ne0.EQ.0) ne0=ne
          nb=NBJ(1,ne)
          n_horsfield=NORD(2,ne) !Horsfield order of current branch
          np0=NPNE(1,nb,ne0)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          length=0.d0
          DO nj=1,NJT
            length=length+(XP(1,1,nj,np1)-XP(1,1,nj,np2))**2
            vector0(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np0)
            vector1(nj)=XP(1,1,nj,np2)-XP(1,1,nj,np1)
          ENDDO !nj
          length=DSQRT(length)
          CALL NORMALISE(3,vector0,ERROR,*9999)
          CALL NORMALISE(3,vector1,ERROR,*9999)
          angle1=DCOS(DOT_PROD(vector0,vector1))
C.........Establish whether major or minor.
C.........Assume major if of higher order (not necessarily true)            
          n_horsfield_p=NORD(2,ne0) !Horsfield order of parent branch
          ne1=NXI(1,1,ne0) !1st child
          IF(NXI(1,0,ne0).EQ.1)THEN
            MAJOR_CHILD=.TRUE.
          ELSE
            IF(ne1.EQ.ne) ne1=NXI(1,2,ne0) !2nd child
            n_horsfield_s=NORD(2,ne1)
            IF(n_horsfield.GT.n_horsfield_s)THEN
              MAJOR_CHILD=.TRUE.  
            ELSE IF(n_horsfield.EQ.n_horsfield_s)THEN
              np3=NPNE(2,nb,ne1)
              length2=0.d0
              DO nj=1,NJT
                length2=length2+(XP(1,1,nj,np1)-XP(1,1,nj,np3))**2
                vector2(nj)=XP(1,1,nj,np3)-XP(1,1,nj,np1)
              ENDDO !nj
              length2=DSQRT(length2)
              CALL NORMALISE(3,vector2,ERROR,*9999)
              angle2=DCOS(DOT_PROD(vector0,vector2))
              IF(length.GT.length2)THEN
                MAJOR_CHILD=.FALSE.
              ELSE
                MAJOR_CHILD=.TRUE.
              ENDIF
c            IF(angle1.GT.angle2)THEN
c              MAJOR_CHILD=.FALSE.
c            ELSE
c              MAJOR_CHILD=.TRUE.
c            ENDIF
            ELSE
              MAJOR_CHILD=.FALSE.
            ENDIF
          ENDIF
          IF(n_horsfield.LT.30)THEN
            np=NPNE(2,nb,ne0)            !node at end of parent
            nv=NVJE(2,nb,nj_radius,ne0)  !version of end parent node
            weight1=0.8d0
            weight2=1.d0-weight1
            IF(MAJOR_CHILD)THEN
              RADIUS=weight1*XP(1,nv,nj_radius,np)*RatioDmajorToDparent
     &          +weight2*HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)
c              RADIUS=0.333d0*XP(1,nv,nj_radius,np)*RatioDmajorToDparent
c     &              +0.333d0*HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)
c     &              +0.333d0*0.5d0*length/RatioLengthToDmajor
c              RADIUS=XP(1,nv,nj_radius,np)*RatioDmajorToDparent
            ELSE
              RADIUS=weight1*XP(1,nv,nj_radius,np)*RatioDminorToDparent
     &          +weight2*HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)
c              RADIUS=0.25d0*XP(1,nv,nj_radius,np)*RatioDminorToDparent
c     &              +0.5d0*HORSFIELD_AIRWAY_DIAM_ADJ(n_horsfield)
c     &          +0.25d0*0.5d0*length/RatioLengthToDminor
c              RADIUS=XP(1,nv,nj_radius,np)*RatioDminorToDparent
            ENDIF
          ELSE
            IF(MAJOR_CHILD)THEN
              RADIUS=XP(1,nv,nj_radius,np)*RatioDmajorToDparent
            ELSE
              RADIUS=XP(1,nv,nj_radius,np)*RatioDminorToDparent
            ENDIF
          ENDIF
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=RADIUS*scale_f
            IF(nj_alveoli.NE.0) THEN
              nv=NVJE(nn,nb,nj_alveoli,ne)
              XP(1,nv,nj_alveoli,np)=1.d0 !conducting, therefore a/A=1
            ENDIF
          ENDDO !nn
        ENDDO !noelem

      ELSE IF(RADIUS_SCHEME(1:8).EQ.'ARTERIES'
     &    .OR.RADIUS_SCHEME(1:5).EQ.'VEINS') THEN

C     Radii equal to measured values from Horsfield (1978) for arteries,
C       or Horsfield (1981) for veins, based on Strahler ordering.
        IF(UNSTRAINED) THEN
          IF(RADIUS_SCHEME(1:8).EQ.'ARTERIES') THEN
            compliance_vessel=0.261d0 !/kPa (from Yen:1981)
            Ptm=0.7365d0 !kPa (from Singhal:1973 - pressure difference)
          ELSE
            compliance_vessel=0.071d0 !/kPa (from Yen:1980)
            Ptm=1.325d0 !kPa (from Horsfield:1981)
          ENDIF
          unstrained_factor=(1.d0+compliance_vessel*Ptm) !scales to unstrained radius values
        ELSE
          unstrained_factor=1.d0 !no scaling done - use values as they are
        ENDIF
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          ne0=NXI(-1,1,ne)
          nb=NBJ(1,ne)
          n_strahler=NORD(3,ne) !Strahler order
          IF(ne0.NE.0)THEN
            np=NPNE(2,nb,ne0) !node at end of parent
            nv=NVJE(2,nb,nj_radius,ne0) !version of end parent node
            IF(RADIUS_SCHEME(1:8).EQ.'ARTERIES') THEN
              IF(XP(1,nv,nj_radius,np)
     &          .LT.HORSFIELD_ARTERY_DIAM(n_strahler))THEN
                RADIUS=XP(1,nv,nj_radius,np)
              ELSE
                RADIUS=HORSFIELD_ARTERY_DIAM(n_strahler)
     &            /unstrained_factor*0.5d0
              ENDIF
            ELSE
              RADIUS=HORSFIELD_ARTERY_DIAM(n_strahler)
     &          /unstrained_factor*0.5d0
            ENDIF
          ELSE IF(RADIUS_SCHEME(1:5).EQ.'VEINS') THEN
            IF(XP(1,nv,nj_radius,
     &        np).LT.HORSFIELD_VEIN_DIAM(n_strahler))THEN
              RADIUS=XP(1,nv,nj_radius,np)
            ELSE
              RADIUS=HORSFIELD_VEIN_DIAM(n_strahler)/unstrained_factor
     &          *0.5d0
            ENDIF
          ELSE
            RADIUS=HORSFIELD_VEIN_DIAM(n_strahler)/unstrained_factor
     &        *0.5d0
          ENDIF
          DO nn=1,2
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=RADIUS
            IF(nj_alveoli.NE.0) THEN
              nv=NVJE(nn,nb,nj_alveoli,ne)
              XP(1,nv,nj_alveoli,np)=1.d0 !conducting, therefore a/A=1
            ENDIF
          ENDDO !nn
        ENDDO !noelem

      ELSE IF(RADIUS_SCHEME(1:10).EQ.'RDSTRAHLER'.OR.
     &    RADIUS_SCHEME(1:11).EQ.'RDHORSFIELD')THEN
C     Radii calculated using the Strahler or Horsfield diameter ratio.

        IF(RADIUS_SCHEME(1:10).EQ.'RDSTRAHLER')THEN
          nindex = 3 !for Strahler ordering
        ELSE IF(RADIUS_SCHEME(1:11).EQ.'RDHORSFIELD')THEN
          nindex = 2 !for Horsfield ordering
        ENDIF

        IF(CUBE_ROOT)THEN
          DO N=1,GENM
            NTALLY(N)=0
          ENDDO !N
          
          NMAX_ORD=0
          DO noelem=1,NEELEM(0)
            ne=NEELEM(noelem)
            n_order=NORD(nindex,ne) !Horsfield or Strahler order
            IF(n_order.GT.NMAX_ORD) NMAX_ORD=n_order !highest order
C           ne0=NXI(-1,1,ne)
c           IF(ne0.NE.0)THEN
cc              IF(n_order.NE.n_order_0)THEN
c                NTALLY(n_order)=NTALLY(n_order)+1
cc              ENDIF
c           ELSE
            !Don't want to sum elements (i.e. no bifurcation) within a branch
            ne0=NXI(1,0,ne)
            IF(ne0.NE.1) THEN 
              NTALLY(n_order)=NTALLY(n_order)+1
            ENDIF
          ENDDO !noelem

          
          DO N=1,NMAX_ORD
            X(N)=N
            YREGRESS(N)=DLOG10(DBLE(NTALLY(N)))
          ENDDO !N
          CALL LINREGRESS(NMAX_ORD,R_SQ,slope,X,YREGRESS,ERROR,*9999)
          RATIO_DIAMETER=(10.d0**DABS(slope))**(1.d0/3.d0) !R_d=R_b^(1/3)
        ELSE
          NMAX_ORD=NORD(nindex,ne_start) !highest order from order evaluation
        ENDIF
        ne=ne_start
        np=NPNE(1,NBJ(1,ne),ne)
        diameter=XP(1,1,nj_radius,np)*2.d0

        CALL ASSERT(diameter.GT.0.d0,
     &    '>>Must specify radius of stem branch first',ERROR,*9999)       
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          n_order=NORD(nindex,ne)
          radius=(10.d0**(DLOG10(RATIO_DIAMETER)*DBLE(n_order
     &      -NMAX_ORD)+DLOG10(diameter)))*0.5d0
          nb=NBJ(nj_radius,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=radius*scale_f
          ENDDO !nn
        ENDDO !noelem
        
      ELSE IF(RADIUS_SCHEME(1:3).EQ.'HBW')THEN
C     Radii equal to measured values from Haefeli-Bleuer & Weibel
C     (1988), based on generation. Only for respiratory airways
        ne0=NXI(-1,1,NELIST(1)) !valid only if 1st in list is at top
        IF(ne0.NE.0)THEN
          NGEN_OFFSET=NORD(1,ne0)
        ELSE
          NGEN_OFFSET=0
        ENDIF
        DO noelem=1,NELIST(0) !alter dimensions for terminal airspaces
          ne=NELIST(noelem)
          n_generation=NORD(1,ne)-NGEN_OFFSET !generation
          IF(n_generation.GT.0)THEN !respiratory generation
            nb=NBJ(1,ne)
            np1=NPNE(1,nb,ne)
            np2=NPNE(2,nb,ne)
            
c            BRANCH_L=HBW_LENGTH(n_generation+1)
c            DO nj=1,NJT
c              XP(1,1,nj,np2)=XP(1,1,nj,np1)+BRANCH_L*XP(1,2,nj,np2)
c            ENDDO !nj
            
            BRANCH_D=HBW_DIAM(n_generation+1)
            TERMINAL=.FALSE.
            CHECK=.FALSE.
            ne2=ne
            DO WHILE(.NOT.CHECK)
              IF(NXI(1,0,ne2).EQ.0)THEN !is terminal
                TERMINAL=.TRUE.
                CHECK=.TRUE.
              ELSE
                ne2=NXI(1,1,ne2)
                IF(NORD(1,ne2).NE.NORD(1,ne))THEN
                  CHECK=.TRUE.
                ELSE
                  ne2=NXI(1,1,ne2)
                ENDIF
              ENDIF
            ENDDO
            
            IF(TERMINAL)THEN
              BRANCH_D=SAC_DIAM_MEAN(n_generation+1)
              AA_RATIO=AREA_RATIO_SAC(n_generation+1)
            ELSE
              AA_RATIO=AREA_RATIO_DUCT(n_generation)
            ENDIF
            
            DO nn=1,2
              np=NPNE(nn,nb,ne)
              nb=NBJ(nj_radius,ne)
              nv=NVJE(nn,nb,nj_radius,ne)
              !R_alveolar=SQRT(r_duct**2/(a_A ratio))
              XP(1,nv,nj_radius,np)=BRANCH_D*0.5d0/DSQRT(AA_RATIO)
     &          *scale_f
              IF(nj_alveoli.NE.0) THEN
                nb=NBJ(nj_alveoli,ne)
                nv=NVJE(nn,nb,nj_alveoli,ne)
                XP(1,nv,nj_alveoli,np)=AA_RATIO !alveolar x-sec area
              ENDIF
            ENDDO !nn

          ENDIF
        ENDDO !noelem0

      ELSE IF(RADIUS_SCHEME.EQ.'CE_FIELD')THEN
C Update the radii and alveolar ratio from stored values in CE
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            IF(SCALED)THEN
              XP(1,nv,nj_radius,np)=CE(1,ne)*DSQRT(XAB(1,ne))
            ELSE
              XP(1,nv,nj_radius,np)=CE(1,ne)
            ENDIF
            IF(nj_alveoli.NE.0) THEN
              XP(1,nv,nj_alveoli,np)=CE(2,ne)
            ENDIF
          ENDDO !nn
        ENDDO !noelem

      ELSE IF(RADIUS_SCHEME.EQ.'LUMPED')THEN
C     Calculate an effective radius for the terminal elements that
C     represent the 'lumped parameter' models. 
        radius=DSQRT(ACINUS_VOLUME/(PI*ACINUS_LENGTH))
        DO noelem=1,NELIST(0)
          ne=NELIST(noelem)
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb)
            np=NPNE(nn,nb,ne)
            nv=NVJE(nn,nb,nj_radius,ne)
            XP(1,nv,nj_radius,np)=radius*scale_f
            XP(1,nv,nj_alveoli,np)=AA_R
          ENDDO !nn
          DO nj=1,NJT
            XP(1,1,nj,NPNE(2,nb,ne))=XP(1,1,nj,NPNE(1,nb,ne))
     &        +ACINUS_LENGTH*XP(1,2,nj,NPNE(1,nb,ne))
          ENDDO !nj
          
        ENDDO !noelem

      ENDIF !RADIUS_SCHEME

      
      CALL EXITS('UPMESH_GEOMETRY')
      RETURN
 9999 CALL ERRORS('UPMESH_GEOMETRY',ERROR)
      CALL EXITS('UPMESH_GEOMETRY')
      RETURN 1
      END



