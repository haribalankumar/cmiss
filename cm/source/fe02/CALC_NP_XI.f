      SUBROUTINE CALC_NP_XI(IBT,INP,NBJ,ne,NENP,np,NPNE,NRLIST,XI,
     '  ERROR,*)

C#### Subroutine: CALC_NP_XI
C###  Description:
C###    CALC_NP_XI calculates the xi location and element number
C###    of a node.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),INP(NNM,NIM,NBFM),NBJ(NJM,NEM),ne,
     '  NENP(NPM,0:NEPM,0:NRM),np,NPNE(NNM,NBFM,NEM),NRLIST(0:NRM)
      REAL*8 XI(*)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ni,nii,nn,nnnp,nonrlist,nr
      LOGICAL COLLAPSEDNODE,COLLAPSEDXI,ISATCOLLAPSE

      CALL ENTERS('CALC_NP_XI',*9999)

C*** Try and find the node within the element
      ne=0
      nonrlist=1
      DO WHILE((ne.EQ.0).AND.(nonrlist.LE.NRLIST(0)))
        nr=NRLIST(nonrlist)
        ne=NENP(np,1,nr)
        nonrlist=nonrlist+1
      ENDDO

      CALL ASSERT(ne.GT.0,'>>No element was found',ERROR,*9999)

      nb=NBJ(1,ne)
      nnnp=0
      DO nn=1,NNT(nb)
        IF(np.EQ.NPNE(nn,nb,ne)) THEN
          nnnp=nn
        ENDIF
      ENDDO !nn
      CALL ASSERT(nnnp.NE.0,'>>Could not find local node',
     '  ERROR,*9999)
      COLLAPSEDNODE=ISATCOLLAPSE(IBT(1,1,nb),INP(1,1,nb),nb,nnnp)
      IF(.NOT.COLLAPSEDNODE) THEN
        DO ni=1,NIT(nb)
          IF(IBT(1,ni,nb).EQ.1) THEN !Lagrange Tensor product
            XI(ni)=DBLE(INP(nnnp,ni,nb)-1)/DBLE(IBT(2,ni,nb))
          ELSE IF(IBT(1,ni,nb).EQ.2) THEN !Hermite Tensor product
            XI(ni)=DBLE(INP(nnnp,ni,nb)-1)
          ELSE IF(IBT(1,ni,nb).EQ.3) THEN !Simplex element
            IF(NKT(1,nb).EQ.1) THEN !Apex node 1
              IF(nnnp.EQ.1) THEN
                XI(ni)=0.0d0
              ELSE IF(nnnp.EQ.2) THEN
                XI(ni)=DBLE(ni-1)
              ELSE
                XI(ni)=1.0d0
              ENDIF
            ELSE IF(NKT(3,nb).EQ.1) THEN !Apex node 3
              IF(nnnp.EQ.1) THEN
                XI(ni)=0.0d0
              ELSE IF(nnnp.EQ.2) THEN
                XI(ni)=DBLE(2-ni)
              ELSE
                XI(ni)=DBLE(ni-1)
              ENDIF
            ENDIF
          ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN !Sector
            IF(IBT(2,ni,nb).EQ.4) THEN
              XI(ni)=DBLE(INP(nnnp,ni,nb)-1)
            ELSE
              XI(ni)=DBLE(INP(nnnp,ni,nb)-1)/DBLE(IBT(2,ni,nb))
            ENDIF
          ELSE
            ERROR='>>Unknown basis type'
            GOTO 9999
          ENDIF
        ENDDO !ni
      ELSE
        DO ni=1,NIT(nb)
          COLLAPSEDXI=.FALSE.
          DO nii=1,NIT(nb)
CC AJPs 191297
c            IF((IBT(1,nii,nb).EQ.5.OR.IBT(1,nii,nb).EQ.
c     '        6).AND.(IBT(3,nii,nb).EQ.ni)) COLLAPSEDXI=.TRUE.
            IF((IBT(1,nii,nb).EQ.5.OR.IBT(1,nii,nb).EQ.
     '        6).AND.(IBT(3,nii,nb).NE.ni)) COLLAPSEDXI=.TRUE.
CC AJPe
          ENDDO ! nii
          IF(COLLAPSEDXI) THEN
            XI(ni)=0.0d0
          ELSE
            IF(IBT(1,ni,nb).EQ.1) THEN !Lagrange tensor product
              XI(ni)=DBLE(INP(nnnp,ni,nb)-1)/DBLE(IBT(2,ni,nb))
            ELSE IF(IBT(1,ni,nb).EQ.2) THEN !Hermite tensor product
              XI(ni)=DBLE(INP(nnnp,ni,nb)-1)
            ELSE IF(IBT(1,ni,nb).EQ.5.OR.IBT(1,ni,nb).EQ.6) THEN !Sector
              IF(IBT(2,ni,nb).EQ.4) THEN
                XI(ni)=DBLE(INP(nnnp,ni,nb)-1)
              ELSE
                XI(ni)=DBLE(INP(nnnp,ni,nb)-1)/DBLE(IBT(2,ni,nb))
              ENDIF
            ELSE
              ERROR='>>Unknown basis type'
              GOTO 9999
            ENDIF
          ENDIF
        ENDDO !ni
      ENDIF

      CALL EXITS('CALC_NP_XI')
      RETURN
 9999 CALL ERRORS('CALC_NP_XI',ERROR)
      CALL EXITS('CALC_NP_XI')
      RETURN 1
      END


