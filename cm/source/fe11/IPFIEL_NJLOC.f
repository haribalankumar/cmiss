      SUBROUTINE IPFIEL_NJLOC(NJJ_START,NKJ,NPNODE,nr,NUM_FIELD,NVJP,
     &  ERROR,*)

C#### Subroutine: IPFIEL_NJLOC
C###  Description:
C###    IPFIEL_NJLOC sets up NJ_LOC for a new field

      IMPLICIT NONE
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NJJ_START,NKJ(NJM,NPM),NPNODE(0:NP_R_M,0:NRM),nr,
     &  NUM_FIELD,NVJP(NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER FREELIST(NJ_LOC_MX),index,nj,njj,njj1,njj2,nonode,np,nrr,
     &  numf,NUMFREE

      CALL ENTERS('IPFIEL_NJLOC',*9999)

C     Set up NJ_LOC(nj,njj)
C     Store free nj's in FREELIST & check enough room
      IF(ADD)THEN
        nj=NJ_LOC(NJL_FIEL,0,nr)
      ELSE
        nj=0
      ENDIF

      NUMFREE=0
      DO WHILE((NUMFREE.NE.NUM_FIELD).AND.(nj.LE.NJM))
        nj=nj+1
        CALL ASSERT(nj.LE.3*NJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '    //'loc00.cmn',ERROR,*9999)
        IF(NJ_TYPE(nj,1).EQ.0.OR.
     '    (.NOT.ADD.AND.NJ_TYPE(nj,1).EQ.NJL_FIEL)) THEN
C         Empty space or old field info in that nj location
          NUMFREE=NUMFREE+1
          CALL ASSERT(NUMFREE.LE.NJ_LOC_MX,'>>Increase NJ_LOC_MX in '
     '      //'loc00.cmn',ERROR,*9999)
          FREELIST(NUMFREE)=nj
        ENDIF
      ENDDO
      CALL ASSERT(nj.LE.NJM,' >>Increase NJM',ERROR,*9999)
      
C     Clear any existing fields unless doing an add
      IF(.NOT.ADD) THEN
        DO njj1=1,NJ_LOC(NJL_FIEL,0,nr)
          nj=NJ_LOC(NJL_FIEL,njj1,nr)
          NJ_TYPE(nj,1)=0
          NJ_TYPE(nj,2)=0
          NJ_LOC(NJL_FIEL,njj1,nr)=0
        ENDDO
        NJ_LOC(NJL_FIEL,0,nr)=0
        NJ_LOC(0,0,nr)=0
        DO njj1=1,3
          DO njj2=1,NJ_LOC(njj1,0,nr)
            nj=NJ_LOC(njj1,njj2,nr)
            IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
          ENDDO
        ENDDO
      ENDIF !NOT ADD
C     Store field in free space
C KSB 2004 - index enables adding of fields      
      index=NJJ_START
      DO numf=1,NUMFREE
        nj=FREELIST(numf)
C       NJ_LOC(NJL_FIEL,numf,nr)=nj
        NJ_LOC(NJL_FIEL,index,nr)=nj
        NJ_TYPE(nj,1)=NJL_FIEL
C        NJ_TYPE(nj,2)=numf
        NJ_TYPE(nj,2)=index
        IF(nj.GT.NJ_LOC(0,0,nr)) NJ_LOC(0,0,nr)=nj
        IF(nj.GT.NJ_LOC(NJL_FIEL,0,0)) NJ_LOC(NJL_FIEL,0,0)=nj
        index=index+1
      ENDDO
      !KSB: this enables addition of field variables in same region
      NJ_LOC(NJL_FIEL,0,nr)=index-1 !NUM_FIELD
      DO nrr=1,NRT
        IF(NJ_LOC(0,0,nrr).GT.NJ_LOC(0,0,0))
     '    NJ_LOC(0,0,0)=NJ_LOC(0,0,nrr)
      ENDDO !nrr
C KAT 21May99: Initialize NKJ in case nodes are not set and CALC_NUNK
C       crashes.  NVJP is initialized in FEMINI.  Check this in case the
C       node has been defined in another region.
      DO nonode=1,NPNODE(0,nr)
        np=NPNODE(nonode,nr)
        DO njj=1,NJ_LOC(NJL_FIEL,0,nr)
          nj=NJ_LOC(NJL_FIEL,njj,nr)
          IF(NVJP(nj,np).EQ.0) NKJ(nj,np)=0
        ENDDO !nj
      ENDDO !np

      CALL EXITS('IPFIEL_NJLOC')
      RETURN
 9999 CALL ERRORS('IPFIEL_NJLOC',ERROR)
      DO njj1=NJJ_START,NJ_LOC(NJL_FIEL,0,nr)
        NJ_LOC(NJL_FIEL,njj1,nr)=0
      ENDDO
      NJ_LOC(NJL_FIEL,0,nr)=NJJ_START-1
      DO njj1=1,3
        DO njj2=1,NJ_LOC(njj1,0,nr)
          IF(NJ_LOC(njj1,njj2,nr).GT.NJ_LOC(0,0,nr))
     '      NJ_LOC(0,0,nr)=NJ_LOC(njj1,njj2,nr)
        ENDDO !njj2
      ENDDO !njj1
      NJ_LOC(0,0,0)=0
      NJ_LOC(NJL_FIEL,0,0)=0
      DO nr=1,NRT
        IF(NJ_LOC(0,0,nr).GT.NJ_LOC(0,0,0)) NJ_LOC(0,0,0)=NJ_LOC(0,0,nr)
        IF(NJ_LOC(NJL_FIEL,0,nr).GT.NJ_LOC(NJL_FIEL,0,0))
     '    NJ_LOC(NJL_FIEL,0,0)=NJ_LOC(NJL_FIEL,0,nr)
      ENDDO !nr
      CALL EXITS('IPFIEL_NJLOC')
      RETURN 1
      END


