      SUBROUTINE GNVEIN(NBJ,NEELEM,NENP,NKJ,NKJE,
     '  NP_INTERFACE,NPLIST,NPNE,NPNODE,nr_host,NRE,nr_vein,NVJE,NVJP,
     &  NXI,CE,SE,XP,ERROR,*)

C#### Subroutine: GNVEIN
C###  Description:
C###    GNVEIN generates a pulmonary venous circulation from a
C###    host tree (either airways or arteries).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'lung00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     '  NP_INTERFACE(0:NPM,0:3),NPLIST(0:NPM),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),nr_host,NRE(NEM),nr_vein,
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM),
     '  NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 CE(NMM,NEM),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb_host,nb_vein,HOST(0:NE_R_M),I,MAP_AIR_VEIN(NP_R_M),ne,
     &  NE_CLOSEST(NE_R_M),ne_host,ne_vein,ne_vein1,nj,noelem,
     &  noelem_next,noelem_vein,noelem_vein1,nonode,nonode_vein,np,
     &  np_next,np_tmp1,np_tmp2,np_vein,np_vein1,np1,np2,TMP(0:NE_R_M)
      REAL*8 dist,length,min_dist,offset
      LOGICAL HALF,USED(NP_R_M)
      
      CALL ENTERS('GNVEIN',*9999)
      
      nonode_vein=NPNODE(0,nr_vein)
      np_vein=NPT(0)
      offset=1.d0
      nb_host=NBJ(1,NEELEM(1,nr_host))
      MAP_AIR_VEIN(NPNE(1,nb_host,NEELEM(1,nr_host)))=0
      noelem_vein=NEELEM(0,nr_vein)
      ne_vein=NET(0)
      IF(NEELEM(0,nr_vein).ne.0) THEN
        DO nonode=1,NP_R_M 
          USED(nonode)=.FALSE. !initialise
        ENDDO
 !there are starting elements -> must find closest matching host element
        nb_vein=NBJ(1,NEELEM(1,nr_vein))
        NPLIST(0)=0
        DO noelem=1,NEELEM(0,nr_vein)
          ne=NEELEM(noelem,nr_vein)
          !First find all terminal nodes in this region
          np2=NPNE(2,nb_vein,ne)
          IF(NENP(np2,0,nr_vein).EQ.1) THEN
            NPLIST(0)=NPLIST(0)+1 !terminal node
            NPLIST(NPLIST(0))=np2
          ENDIF
        ENDDO !noelem
        DO nonode_vein=1,NPLIST(0) !for each terminal node
          np=NPLIST(nonode_vein)
          np_tmp1=np
          np_tmp2=1
          min_dist=1.d6 !initialise
          DO noelem=1,NEELEM(0,nr_host)
            ne=NEELEM(noelem,nr_host)
 !loop through all host nodes to find closest matching node
            np2=NPNE(2,nb_host,ne)
            dist=0.d0
            DO nj=1,NJT
              dist=dist+(XP(1,1,nj,np2)-XP(1,1,nj,np))**2.d0
            ENDDO
            dist=DSQRT(dist)
            IF(dist.LT.min_dist.AND..NOT.USED(np2)) THEN
              min_dist=dist
              NE_CLOSEST(np)=ne
              USED(np2)=.TRUE. !keeps track of nodes already allocated
              !can't use nodes more than once.
              IF(np.EQ.np_tmp1.AND.np_tmp2.NE.np2) USED(np_tmp2)=.FALSE.
              np_tmp1=np
              np_tmp2=np2
            ENDIF
          ENDDO !nonode_vein
        ENDDO !nonode
 !Again loop over each terminal venous element and create venous mesh,
 !until terminal host element is reached.
        DO nonode=1,NPLIST(0)
          np=NPLIST(nonode)
          ne_host=NE_CLOSEST(np)
          MAP_AIR_VEIN(NPNE(1,nb_host,ne_host))=np
          HOST(0)=1
          HOST(HOST(0))=ne_host
          DO WHILE(HOST(0).GT.0)
 !do until terminal node reached
            TMP(0)=0
            DO noelem_next=1,HOST(0)
              ne_host=HOST(noelem_next)
              np_next=NPNE(2,nb_host,ne_host) !2nd node of element
              IF(NXI(1,0,ne_host).EQ.2)THEN
                nonode_vein=nonode_vein+1
                np_vein=np_vein+1
                CALL ASSERT(nonode_vein.LE.NP_R_M,'>>Increase NP_R_M',
     '            ERROR,*9999)
                CALL ASSERT(np_vein.LE.NPM,'>>Increase NPM',
     '            ERROR,*9999)              
                NPNODE(nonode_vein,nr_vein)=np_vein
                np1=NPNE(2,nb_host,NXI(1,1,ne_host))
                np2=NPNE(2,nb_host,NXI(1,2,ne_host))
                TMP(0)=TMP(0)+1
                TMP(TMP(0))=NXI(1,1,ne_host)
                TMP(0)=TMP(0)+1
                TMP(TMP(0))=NXI(1,2,ne_host)
                DO nj=1,NJT
                  XP(1,1,nj,np_vein)=(XP(1,1,nj,np1)*0.5d0
     '              +XP(1,1,nj,np2)*0.5d0+XP(1,1,nj,np_next)*0.d0)
                ENDDO
                MAP_AIR_VEIN(np_next)=np_vein
              ELSE IF(NXI(1,0,ne_host).EQ.1)THEN
                WRITE(OP_STRING,
     '            '('' Only single daughter, from parent: '',I5)')
     &            ne_host
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                nonode_vein=nonode_vein+1
                np_vein=np_vein+1
                CALL ASSERT(nonode_vein.LE.NP_R_M,'>>Increase NP_R_M',
     '            ERROR,*9999)
                CALL ASSERT(np_vein.LE.NPM,'>>Increase NPM',
     '            ERROR,*9999)              
                NPNODE(nonode_vein,nr_vein)=np_vein
                TMP(0)=TMP(0)+1
                TMP(TMP(0))=NXI(1,1,ne_host)
                DO nj=1,NJT
                  XP(1,1,nj,np_vein)=XP(1,1,nj,np_next)+offset
                ENDDO
                MAP_AIR_VEIN(np_next)=np_vein
!               ELSE IF(NXI(1,0,ne).EQ.0) THEN 
              ENDIF
C... create elements on the fly              
              IF(NXI(1,0,ne_host).NE.0) THEN
                noelem_vein=noelem_vein+1
                ne_vein=ne_vein+1
                CALL ASSERT(noelem_vein.LE.NE_R_M,'>>Increase NP_R_M',
     '            ERROR,*9999)
                CALL ASSERT(ne_vein.LE.NEM,'>>Increase NPM',
     '            ERROR,*9999)    
                NEELEM(noelem_vein,nr_vein)=ne_vein
                NPNE(1,nb_vein,ne_vein)=
     &            MAP_AIR_VEIN(NPNE(1,nb_host,ne_host))
                NPNE(2,nb_vein,ne_vein)=np_vein
                CALL GN1DNEJ(nb_vein,NBJ,ne_vein,NKJ,NKJE,
     &            NPNE(2,nb_vein,ne_vein),NPNE(1,nb_vein,ne_vein),
     &            nr_vein,NRE,NVJE,NVJP,SE,ERROR,*9999)
              ENDIF
            ENDDO !ne_host
            HOST(0)=0
            DO i=1,TMP(0)
              HOST(0)=HOST(0)+1
              HOST(HOST(0))=TMP(i)
            ENDDO !i
          ENDDO !WHILE
        ENDDO !nonode_vein
      ELSE
        DO noelem=1,NEELEM(0,nr_host)
          ne=NEELEM(noelem,nr_host)
          np=NPNE(2,nb_host,ne) !2nd node of element
          IF(NXI(1,0,ne).EQ.2)THEN
            nonode_vein=nonode_vein+1
            np_vein=np_vein+1
            CALL ASSERT(nonode_vein.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np_vein.LE.NPM,'>>Increase NPM',
     '        ERROR,*9999)              
            NPNODE(nonode_vein,nr_vein)=np_vein
            np1=NPNE(2,nb_host,NXI(1,1,ne))
            np2=NPNE(2,nb_host,NXI(1,2,ne))
            DO nj=1,NJT
              XP(1,1,nj,np_vein)=(XP(1,1,nj,np1)*0.5d0
     '          +XP(1,1,nj,np2)*0.5d0+XP(1,1,nj,np)*0.d0)
            ENDDO
            MAP_AIR_VEIN(np)=np_vein
          ELSE IF(NXI(1,0,ne).EQ.1)THEN
            WRITE(OP_STRING,
     '        '('' Only single daughter, from parent: '',I5)') ne
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            nonode_vein=nonode_vein+1
            np_vein=np_vein+1
            CALL ASSERT(nonode_vein.LE.NP_R_M,'>>Increase NP_R_M',
     '        ERROR,*9999)
            CALL ASSERT(np_vein.LE.NPM,'>>Increase NPM',
     '        ERROR,*9999)              
            NPNODE(nonode_vein,nr_vein)=np_vein
            DO nj=1,NJT
              XP(1,1,nj,np_vein)=XP(1,1,nj,np)+offset
            ENDDO
            MAP_AIR_VEIN(np)=np_vein
!           ELSE IF(NXI(1,0,ne).EQ.0) THEN 
          ENDIF
        ENDDO !noelem
C      NPNODE(0,nr_vein)=nonode_vein
C... Now create venous elements      
        DO noelem=1,NEELEM(0,nr_host)
          ne=NEELEM(noelem,nr_host)
          np1=NPNE(1,nb_host,ne)
          np2=NPNE(2,nb_host,ne)
          IF(MAP_AIR_VEIN(np1).GT.0.AND.MAP_AIR_VEIN(np2).GT.0.AND.
     '      NXI(1,0,ne).NE.0)THEN
            noelem_vein=noelem_vein+1
            ne_vein=ne_vein+1
            NEELEM(noelem_vein,nr_vein)=ne_vein
            NPNE(1,nb_vein,ne_vein)=MAP_AIR_VEIN(np1)
            NPNE(2,nb_vein,ne_vein)=MAP_AIR_VEIN(np2)
            CALL GN1DNEJ(nb_vein,NBJ,ne_vein,NKJ,NKJE,NPNE(2,nb_vein,
     &        ne_vein),NPNE(1,nb_vein,ne_vein),nr_vein,NRE,NVJE,
     '        NVJP,SE,ERROR,*9999)
            length=0.d0
            DO nj=1,NJT 
              length=length+(XP(1,1,nj,NPNE(2,nb_vein,ne_vein))-
     '          XP(1,1,nj,NPNE(1,nb_vein,ne_vein)))**2.d0
            ENDDO
            length=DSQRT(length)
            CE(1,ne_vein)=length !stores segment length
          ENDIF
        ENDDO !noelem        
      ENDIF
      
      NEELEM(0,nr_vein)=noelem_vein
      NPNODE(0,0)=nonode_vein
      NEELEM(0,0)=noelem_vein 
      NPT(nr_vein)=np_vein
      NET(nr_vein)=ne_vein
      NPT(0)=NPT(nr_vein)
      NET(0)=NET(nr_vein)
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb_vein,NEELEM,NENP,NPNE,nr_vein,ERROR,*9999)
      CALL NENXI_1D(nb_vein,NEELEM,NENP,NPNE,nr_vein,NXI,ERROR,*9999)

C... Create new venous vessels if there is only 2 branches at a
C... bifucation        
      DO noelem=1,NEELEM(0,nr_host)
        ne=NEELEM(noelem,nr_host)
        np1=NPNE(1,nb_host,ne)
        np2=NPNE(2,nb_host,ne)
        IF(MAP_AIR_VEIN(np1).GT.0) THEN
          np_vein1=MAP_AIR_VEIN(np1)
          IF(NENP(np_vein1,0,nr_vein).EQ.2) THEN !only 2 ne's at bifurcation 
            noelem_vein1=0
            DO WHILE(NENP(np_vein1,0,nr_vein).LT.3.AND.noelem_vein1
     '        .LT.NENP(np_vein1,0,nr_vein))
              noelem_vein1=noelem_vein1+1
              ne_vein1=NENP(np_vein1,noelem_vein1,nr_vein)
              IF(NXI(1,0,ne_vein1).EQ.1.AND.NXI(1,0,ne).GT.1) THEN
C... Create new venous vessel in direction of one of the host
C... counterpart branches
                HALF=.FALSE.
                CALL VEIN_ADD(nb_vein,NBJ,noelem_vein,nonode_vein,ne,
     &            NEELEM,NENP,ne_vein,ne_vein1,NKJ,NKJE,np1,np2,NPNE,
     &            NPNODE,np_vein,np_vein1,NRE,nr_vein,NVJE,NVJP,NXI,CE,
     &            SE,XP,HALF,ERROR,*9999)
              ENDIF
            ENDDO !WHILE
          ENDIF
        ENDIF
        IF(NENP(np2,0,nr_host).GT.2) THEN !a bifurcation (at least)
 !If host node is a bifurcation the venous node should also be bifurcation
          IF(MAP_AIR_VEIN(np2).GT.0) THEN
            np_vein1=MAP_AIR_VEIN(np2)
            IF(NENP(np_vein1,0,nr_vein).LE.2) THEN !must add new node and element
              noelem_vein1=0
              HALF=.FALSE.
              DO WHILE(NENP(np_vein1,0,nr_vein).LT.3.AND.noelem_vein1
     '          .LT.NENP(np_vein1,0,nr_vein))
                noelem_vein1=noelem_vein1+1
                ne_vein1=NENP(np_vein1,noelem_vein1,nr_vein)
                CALL VEIN_ADD(nb_vein,NBJ,noelem_vein,nonode_vein,ne,
     &            NEELEM,NENP,ne_vein,ne_vein1,NKJ,NKJE,np1,np2,NPNE,
     &            NPNODE,np_vein,np_vein1,NRE,nr_vein,NVJE,NVJP,NXI,CE,
     &            SE,XP,HALF,ERROR,*9999)
                HALF=.TRUE.
              ENDDO 
            ENDIF
          ENDIF !MAP_AIR_VEIN(np).GT.0
        ENDIF !NENP(np,0,nr_host)
      ENDDO !noelem
      NPNODE(0,nr_vein)=nonode_vein
      NEELEM(0,nr_vein)=noelem_vein
      NPNODE(0,0)=np_vein
      NEELEM(0,0)=ne_vein 
      NPT(nr_vein)=np_vein
      NET(nr_vein)=ne_vein
      NPT(0)=NPT(nr_vein)
      NET(0)=NET(nr_vein)
      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP_1D(nb_vein,NEELEM,NENP,NPNE,nr_vein,ERROR,*9999)
      CALL NENXI_1D(nb_vein,NEELEM,NENP,NPNE,nr_vein,NXI,ERROR,*9999)

      
      CALL EXITS('GNVEIN')
      RETURN
 9999 CALL ERRORS('GNVEIN',ERROR)
      CALL EXITS('GNVEIN')
      RETURN 1
      END

