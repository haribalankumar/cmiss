      SUBROUTINE GNARTRY(MIN_ORDER,nb,NBJ,NEELEM,NENP,NKJ,NKJE,NORD,
     '  NP_INTERFACE,NPNE,NPNODE,nr_airway,NRE,nr_artery,NVJE,NVJP,NXI,
     '  CE,SE,XP,ADD_SUPER,ERROR,*)
      
C#### Subroutine: GNARTRY
C###  Description:
C###    GNARTRY generates a pulmonary arterial circulation from an
C###    airway tree.
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'
      
!     Parameter List
      INTEGER MIN_ORDER,nb,NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NORD(5,NE_R_M),NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr_airway,NRE(NEM),nr_artery,
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),
     '  XP(NKM,NVM,NJM,NPM)
      LOGICAL ADD_SUPER
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER gen_next,MAP_AIR_ARTERY(NP_R_M),ne,ne0,
     '  ne_artery,ngen,nj,noelem,noelem_artery,nonode_artery,np_artery,
     '  np,np1,np2
      REAL*8 diameter,length,length_V,SUM_LENGTH,V(NJT)
      
      CALL ENTERS('GNARTRY',*9999)

c      IF(.NOT.DIMENSIONS) THEN
C       MAP_AIR_ARTERY(NPNE(1,nb,NEELEM(1,nr_airway)))=0 !1st airway node not used
        noelem_artery=0 
        ne_artery=NET(0) 
        nonode_artery=0 
        np_artery=NPT(0)         
C... KSB: 02/2003 redone, now similar to GNVEIN
C... Creates arterial nodes      
        DO noelem=1,NEELEM(0,nr_airway)
          ne=NEELEM(noelem,nr_airway)
          ngen=NORD(1,ne)
C          IF(ngen.GE.1) THEN
            np=NPNE(2,nb,ne) !end node of element
            IF(CE(4,ne).NE.0.d0) THEN !CE(4,ne)=radius
              diameter=CE(4,ne)*2.d0
            ELSE !if diameters not defined use length/diam ratio=3
C... Calculate total airway branch length          
              SUM_LENGTH=CE(1,ne)
              IF(NXI(1,0,ne).NE.0) THEN
                ne0=NXI(1,1,ne)
                gen_next=NORD(1,ne0) !generation
                ne0=ne
                DO WHILE(gen_next.EQ.NORD(1,ne))
                  ne0=NXI(1,1,ne0) !neighbouring element
                  SUM_LENGTH=SUM_LENGTH+CE(1,ne0)
                  IF(NXI(1,1,ne0).NE.0)THEN
                    gen_next=NORD(1,NXI(1,1,ne0))
                  ELSE
                    gen_next=NORD(1,ne)+1
                  ENDIF
                ENDDO
              ENDIF
              IF(NXI(-1,0,ne).NE.0)THEN
                ne0=NXI(-1,1,ne)
                gen_next=NORD(1,ne0)
                ne0=ne
                DO WHILE(gen_next.EQ.NORD(1,ne))
                  ne0=NXI(-1,1,ne0) !neighbouring element
                  SUM_LENGTH=SUM_LENGTH+CE(1,ne0)
                  gen_next=NORD(1,NXI(-1,1,ne0))
                ENDDO
              ENDIF
              diameter=SUM_LENGTH/3.d0
            ENDIF !CE(4,ne).NE.0.d0
            IF(ngen.EQ.1) THEN !add pulmonary trunk
              nonode_artery=nonode_artery+1
              np_artery=np_artery+1
              CALL ASSERT(nonode_artery.LE.NP_R_M,'>>Increase NP_R_M',
     '          ERROR,*9999)
              CALL ASSERT(np_artery.LE.NPM,'>>Increase NPM',ERROR,*9999)
              NPNODE(nonode_artery,nr_artery)=np_artery
              np1=NPNE(1,nb,ne)
              length_V=0.d0
              DO nj=1,NJT
                V(nj)=XP(1,1,nj,np1)-XP(1,1,nj,np) !vector from np1-np
                length_V=length_V+V(nj)**2.d0
              ENDDO
              length_V=DSQRT(length_V)
              XP(1,1,1,np_artery)=XP(1,1,1,np1)+diameter
              XP(1,1,2,np_artery)=XP(1,1,2,np1)+diameter
              XP(1,1,3,np_artery)=XP(1,1,3,np)+V(3)
     '          *(HORSFIELD_ARTERY_LENGTH(12)/length_V)
              MAP_AIR_ARTERY(np1)=np_artery !***
            ENDIF
            nonode_artery=nonode_artery+1
            np_artery=np_artery+1
            CALL ASSERT(nonode_artery.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np_artery.LE.NPM,'>>Increase NPM',ERROR,*9999)
            NPNODE(nonode_artery,nr_artery)=np_artery
            XP(1,1,1,np_artery)=XP(1,1,1,np)+diameter
            XP(1,1,2,np_artery)=XP(1,1,2,np)+diameter
            XP(1,1,3,np_artery)=XP(1,1,3,np)
            MAP_AIR_ARTERY(np)=np_artery
C           ENDIF
        ENDDO !noelem
C... Creates arterial elements      
        DO noelem=1,NEELEM(0,nr_airway)
          ne=NEELEM(noelem,nr_airway)
          np1=NPNE(1,nb,ne)
          np2=NPNE(2,nb,ne)
          IF(MAP_AIR_ARTERY(np1).GT.0) THEN
            noelem_artery=noelem_artery+1
            ne_artery=ne_artery+1
            NEELEM(noelem_artery,nr_artery)=ne_artery
            NPNE(1,nb,ne_artery)=MAP_AIR_ARTERY(np1)
            NPNE(2,nb,ne_artery)=MAP_AIR_ARTERY(np2)
            NORD(1,ne_artery)=NORD(1,ne)
            CALL GN1DNEJ(nb,NBJ,ne_artery,NKJ,NKJE,NPNE(2,nb,ne_artery),
     '        NPNE(1,nb,ne_artery),nr_artery,NRE,NVJE,
     '        NVJP,SE,ERROR,*9999)
            length=0.d0 
            DO nj=1,NJT 
              length=length+(XP(1,1,nj,NPNE(2,nb,ne_artery))-
     '          XP(1,1,nj,NPNE(1,nb,ne_artery)))**2.d0
            ENDDO
            length=DSQRT(length)
            CE(1,ne_artery)=length !stores segment length
          ENDIF
        ENDDO !noelem
        
        IF(ADD_SUPER) THEN !add supernumarary branches to network
          NPT(nr_artery)=np_artery
          NET(nr_artery)=ne_artery
          NPT(0)=NPT(nr_artery)
          NET(0)=NET(nr_artery)
          NEELEM(0,nr_artery)=noelem_artery
C... Need to call this here to determine diameters, called at end also
          CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr_artery,ERROR,*9999)
          CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr_artery,NXI,ERROR,*9999)
          CALL BRANCH_ORD(NEELEM(0,nr_artery),NORD,NXI,ERROR,*9999)
          CALL DIAM_DEF_STRAHLER(MIN_ORDER,NEELEM,NORD,nr_artery,NXI,
     '      CE,HORSFIELD_ARTERY_DIAM,ERROR,*9999)
C... Adding supernumerary side branches at right angle to main stem
C          CALL SUPERNUMERARY(nb,NBJ,NEELEM,NELIST,NELIST2,NENP,NKJ,
C     &      NKJE,NORD,NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,NXI,CE,
C     &      D_RATIO,DIAM_RATIO,D_MAX_ORDER,L_BRANCH,SE,SV_FREQ,XP,
C     &      ZD,VEINS,ERROR,*9999)
        ELSE !this done in SUPERNUMERARY if ADD_SUPER=.TRUE.
          NPNODE(0,nr_artery)=nonode_artery
          NEELEM(0,nr_artery)=noelem_artery
          NPT(nr_artery)=nonode_artery
          NET(nr_artery)=noelem_artery
          NPT(0)=np_artery
          NET(0)=ne_artery
          NPNODE(0,0)=np_artery
          NEELEM(0,0)=ne_artery
          CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
          CALL CALC_NENP_1D(nb,NEELEM,NENP,NPNE,nr_artery,ERROR,*9999)
          CALL NENXI_1D(nb,NEELEM,NENP,NPNE,nr_artery,NXI,ERROR,*9999)
        ENDIF !ADD_SUPER
c      ENDIF !DIMENSION
      
C... Determines branch generation numbers & allocates vessel diameters
      DO noelem=1,NEELEM(0,nr_artery)
        ne=NEELEM(noelem,nr_artery)
        ne0=NXI(-1,1,ne)
        IF(ne0.NE.0)THEN
          NORD(1,ne)=NORD(1,ne0)+1
        ELSE
          NORD(1,ne)=1
        ENDIF
      ENDDO
      
      CALL BRANCH_ORD(NEELEM(0,nr_artery),NORD,NXI,ERROR,*9999)
      CALL DIAM_DEF_STRAHLER(MIN_ORDER,NEELEM,NORD,nr_artery,NXI,CE,
     '  HORSFIELD_ARTERY_DIAM,ERROR,*9999)
      
      CALL EXITS('GNARTRY')
      RETURN
 9999 CALL ERRORS('GNARTRY',ERROR)
      CALL EXITS('GNARTRY')
      RETURN 1
      END


      
