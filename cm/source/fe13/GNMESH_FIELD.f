      SUBROUTINE GNMESH_FIELD(NBJ,NEELEM,NENP,NPNODE,NKJ,NKJE,
     &  NOVERSIONS,NPNE,nr,NUM_FIELD,NVJE,NVJP,SE,XP,ERROR,*)
      
      
C#### Subroutine: GNMESH_FIELD
C###  Description:
C###  GNMESH_FIELD sets up field structure for pulmonary meshes. Need to
C###  then call UPMESH to update the mesh geometry and input field
C###  values.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'call00.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     &  NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),NKJE(NKM,NNM,NJM,NEM),
     &  NOVERSIONS(0:NJM),NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),nr,
     &  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)      
!     Local variables
      INTEGER i,nb,FREELIST(NJ_LOC_MX),ne,ne2,nj,njj,nk,nn,
     &  noelem,nonode,np,np1,np2,nrr,ntemp,numf,NUM_FIELD,NUMFREE,nv
      REAL*8 length
      CHARACTER STRING*255
      LOGICAL INLIST
      
      CALL ENTERS('GNMESH_FIELD',*9999)

C     Set up fields for length, radii, alveoli, etc.
C     Check to see whether already set up!!
      IF(.NOT.ADD)THEN !set up fields
        NUMFREE=NUM_FIELD
        DO numf=1,NUMFREE
          FREELIST(numf)=4+numf
          nj=FREELIST(numf)
          NJ_LOC(NJL_FIEL,numf,nr)=nj
          NJ_TYPE(nj,1)=NJL_FIEL
          NJ_TYPE(nj,2)=numf
          IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
          IF(nj.GT.NJ_LOC(NJL_FIEL,0,0)) NJ_LOC(NJL_FIEL,0,0)=nj
        ENDDO

        NJ_LOC(NJL_FIEL,0,nr) = NUMFREE
        
        DO nrr=1,NRT
          IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     '      NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
        ENDDO !nrr
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO numf=1,NUMFREE
            nj=FREELIST(numf)
            IF(INLIST(numf,NOVERSIONS(1),NOVERSIONS(0),ntemp))THEN
              NVJP(nj,np)=1
            ELSE
              NVJP(nj,np)=NENP(np,0,nr) !# of versions = # of adjacent elements
            ENDIF
            NKJ(nj,np)=1 !# of derivatives = 0
          ENDDO !numf
        ENDDO !np
        
        CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).LE.NJ_LOC_MX,
     '    '>>Increase dimension of PROMPT_NV',ERROR,*9999)

      ELSE IF(ADD)THEN
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO njj=1,NUM_FIELD
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            IF(INLIST(njj,NOVERSIONS(1),NOVERSIONS(0),ntemp))THEN
              NVJP(nj,np)=1
            ELSE
              NVJP(nj,np)=NENP(np,0,nr) !# of versions = # of adjacent elements
            ENDIF
            NKJ(nj,np)=1 !# of derivatives = 0
          ENDDO
        ENDDO !np
        
        CALL ASSERT(NJ_LOC(NJL_FIEL,0,nr).LE.NJ_LOC_MX,
     '    '>>Increase dimension of PROMPT_NV',ERROR,*9999)

      ENDIF !ADD

C     Always do for all elements, whether new or old      
      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
C       Set up direction vector, stored in nv=2
        np1=NPNE(1,nb,ne)
        np2=NPNE(2,nb,ne)
        length=0.d0
        DO nj=1,NJT
          length=length+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2
        ENDDO !nj
        length=DSQRT(length)
        IF(NVM.LT.3) THEN
          WRITE(ERROR,'(''>>Increase NVM to 3'')')
          GO TO 9999
        ENDIF
        DO nj=1,NJT
          XP(1,2,nj,np2)=(XP(1,1,nj,np2)-XP(1,1,nj,np1))/length
        ENDDO !nj
      ENDDO !noelem

      DO noelem=1,NEELEM(0,nr)
        ne=NEELEM(noelem,nr)
        nb=NBJ(1,ne)
        DO njj=1,NUM_FIELD
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          NBJ(nj,ne)=nb !basis function the same as for geometry
        ENDDO
        DO nn=1,NNT(nb)
          np=NPNE(nn,nb,ne)
C         Set up the number of versions for the field types
          DO i=1,NENP(np,0,nr)
            ne2=NENP(np,i,nr)
            IF(ne2.EQ.ne)THEN
              DO njj=1,NUM_FIELD
                nj=NJ_LOC(NJL_FIEL,njj,nr)
                IF(INLIST(njj,NOVERSIONS(1),NOVERSIONS(0),ntemp))THEN
                  NVJE(nn,nb,nj,ne)=1
                ELSE
                  NVJE(nn,nb,nj,ne)=i !the version of node at nn'th position
                ENDIF
              ENDDO
            ENDIF
          ENDDO !i
          DO njj=1,NUM_FIELD
            nj=NJ_LOC(NJL_FIEL,njj,nr)
            DO nk=1,NKT(nn,nb)
              NKJE(nk,nn,nj,ne)=1
            ENDDO !nk
          ENDDO
          SE(nn,nb,ne)=1.d0
        ENDDO !nn
      ENDDO !noelem
      CALL_FIEL=.TRUE.

      CALL EXITS('GNMESH_FIELD')
      RETURN
 9999 CALL ERRORS('GNMESH_FIELD',ERROR)
      CALL EXITS('GNMESH_FIELD')
      RETURN 1
      END

      
