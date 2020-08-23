      SUBROUTINE IPMESH3(IBT,NBJ,NEELEM,NENP,NKJE,NKJ,
     '  NP_INTERFACE,NPNE,NPNODE,NRE,NVJE,NVJP,SE,XP,ERROR,*)

C#### Subroutine: IPMESH3
C###  Description:
C###    IPMESH3 defines mesh parameters for a spherical mesh (3d) or
C###    circular mesh (2d).

C**** Dependent variable information is not set up (although material
C**** parameters are requested they need to be reentered in define
C**** material).
C**** FOR 3D:
C**** Also determines the analytic solution for a simple dipole at the
C**** centre of the spheres.  To become more general we need to include
C**** offsets for each sphere and also the offset for a dipole(s) inside
C**** the sphere.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh03.cmn'

!     Parameter List
      INTEGER IBT(3,NIM,NBFM),NBJ(NJM,NEM),
     '  NEELEM(0:NE_R_M,0:NRM),NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NP_INTERFACE(0:NPM,0:3),NPNE(NNM,NBFM,NEM),
     '  NPNODE(0:NP_R_M,0:NRM),NRE(NEM),NVJE(NNM,NBFM,NJM,NEM),
     '  NVJP(NJM,NPM)
      REAL*8 SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER i,ICHAR,INFO,m1,m2,nb,nbb,nb1,nb2,ne,ni,nj,nk,ncirc,nn,
     '  nspher,nonr,ns,np,noelem,nonode,NOQUES,nocircle,nosphere,
     '  NP_INCR,NPSTART,nr,REGLIST(0:2,NSPHEREMX-1),
     '  SPHERELIST(0:2,NSPHEREMX)
      REAL*8 dRdXI,PHI,RbydRdXI,RdTHETAdXI,THETA,S1_DIR(3),S2_DIR(3)
      CHARACTER CHAR1*1,CHAR1B*1
      LOGICAL BEMREGION,CALC_DERIV,FEMREGION,FEMCIRCLE,FILEIP,NKTOBIG

      CALL ENTERS('IPMESH3',*9999)

      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='($,'' Enter the number of regions for the mesh '
     '  //'[1]: '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NR
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,99,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) MESH3_NR=IDATA(1)
      NRT=MESH3_NR
      CALL ASSERT(NRT.LE.NRM,'>>Increase NRM',ERROR,*9999)
      CALL ASSERT(NRT.LE.NSPHEREMX-1,
     '  '>>Increase NSPHEREMX in mesh03.cmn',ERROR,*9999)
      BEMREGION=.FALSE.
      FEMREGION=.FALSE.
      MESH3_BEMTYPE=2
      DO nr=1,NRT
        CALL INIT_NJ_LOC(NJT,NJL_GEOM,nr,ERROR,*9999)
C GMH 14/2/97 Destroy this region and recreate -
C             necessary for cmgui locking
        IF(NPNODE(0,nr).GT.0) THEN
          CALL REGION_DESTROY(nr,ERROR,*9999)
        ENDIF
        IF(nr.EQ.1) THEN
          IDEFLT(1)=1
        ELSE
          IDEFLT(1)=MESH3_TYPE(nr-1)
        ENDIF
        WRITE(CHAR1,'(I1)') nr
        WRITE(CHAR1B,'(I1)') IDEFLT(1)
        FORMAT='('' The mesh type for region '//CHAR1//
     '    ' is ['//CHAR1B//']: '''//
     '    '/''   (1) Boundary Element (BE) mesh'''//
     '    '/''   (2) Finite Element (FE) mesh'''//
     '    '/$,''    '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=MESH3_TYPE(nr)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) MESH3_TYPE(nr)=IDATA(1)
        IF(MESH3_TYPE(nr).EQ.1) THEN
          BEMREGION=.TRUE.
        ELSE IF(MESH3_TYPE(nr).EQ.2) THEN
          FEMREGION=.TRUE.
        ENDIF
        IF(nr.EQ.1.AND.MESH3_TYPE(1).EQ.1) THEN
          IF(NJT.EQ.2) THEN
            FORMAT='('' The type of BE mesh for region 1 is [1]: '''
     '        //'/''   (1) Circular'''
     '        //'/''   (2) Annular'''
     '        //'/$,''    '',I1)'
          ELSE
            FORMAT='('' The type of BE mesh for region 1 is [1]: '''
     '        //'/''   (1) Spherical'''
     '        //'/''   (2) Annular'''
     '        //'/$,''    '',I1)'
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_BEMTYPE
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_BEMTYPE=IDATA(1)
        ENDIF
      ENDDO !nr

      NSPHERES=NRT+(MESH3_BEMTYPE-1)

C***  REGLIST(0:2,nr) is the list of spheres/circles in region nr
C***  SPHERELIST(0:2,nosphere) is the list of regions that the sphere/
C***    circle is in.
      DO nr=1,NRT
        DO i=0,2
          REGLIST(i,nr)=0
        ENDDO !i
      ENDDO !nr
      DO nosphere=1,NSPHERES
        DO i=0,2
          SPHERELIST(i,nosphere)=0
        ENDDO !i
      ENDDO
      nosphere=0
      DO nr=1,NRT
        IF(nr.EQ.1) THEN
          IF(MESH3_TYPE(1).EQ.1) THEN
            IF(MESH3_BEMTYPE.EQ.1) THEN
              nosphere=1
              SPHERELIST(0,1)=1
              SPHERELIST(1,1)=1
              REGLIST(0,1)=1
              REGLIST(1,1)=1
            ELSE IF(MESH3_BEMTYPE.EQ.2) THEN
              nosphere=2
              SPHERELIST(0,1)=1
              SPHERELIST(1,1)=1
              SPHERELIST(0,2)=1
              SPHERELIST(1,2)=1
              REGLIST(0,1)=2
              REGLIST(1,1)=1
              REGLIST(2,1)=2
            ENDIF
          ELSE IF(MESH3_TYPE(1).EQ.2) THEN
            nosphere=2
            SPHERELIST(0,1)=1
            SPHERELIST(1,1)=1
            SPHERELIST(0,2)=1
            SPHERELIST(1,2)=1
            REGLIST(0,1)=2
            REGLIST(1,1)=1
            REGLIST(2,1)=2
          ENDIF
        ELSE
          nosphere=nosphere+1
          SPHERELIST(0,nosphere)=1
          SPHERELIST(1,nosphere)=nr
          SPHERELIST(0,nosphere-1)=SPHERELIST(0,nosphere-1)+1
          SPHERELIST(SPHERELIST(0,nosphere-1),nosphere-1)=nr
          REGLIST(0,nr)=2
          REGLIST(1,nr)=nosphere-1
          REGLIST(2,nr)=nosphere
        ENDIF
      ENDDO !nr

      IF(DOP) THEN
        WRITE(OP_STRING,'('' Number of regions: '',I2)') NRT
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nr=1,NRT
          IF(MESH3_TYPE(nr).EQ.1) THEN
            IF(nr.EQ.1) THEN
              IF(MESH3_BEMTYPE.EQ.1) THEN
                IF(NJT.EQ.2) THEN
                  WRITE(OP_STRING,'('' Region '',I2,'' is a '
     '              //'circular BE region'')') nr
                ELSE IF(NJT.EQ.3) THEN
                  WRITE(OP_STRING,'('' Region '',I2,'' is a '
     '              //'spherical BE region'')') nr
                ENDIF
              ELSE IF(MESH3_BEMTYPE.EQ.2) THEN
                WRITE(OP_STRING,'('' Region '',I2,'' is an '
     '              //'annular BE region'')') nr
              ENDIF
            ELSE
              WRITE(OP_STRING,'('' Region '',I2,'' is a BE region'')')
     '          nr
            ENDIF
          ELSE IF(MESH3_TYPE(nr).EQ.2) THEN
            WRITE(OP_STRING,'('' Region '',I2,'' is a FE region'')') nr
          ENDIF
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nr
        IF(NJT.EQ.2) THEN
          WRITE(OP_STRING,'('' Number of circles: '',I3)') nspheres
        ELSE IF(NJT.EQ.3) THEN
          WRITE(OP_STRING,'('' Number of spheres: '',I3)') nspheres
        ENDIF
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nr=1,NRT
          WRITE(OP_STRING,'('' REGLIST(0..,'',I2,''): '',3I2)')
     '      nr,(REGLIST(i,nr),i=0,REGLIST(0,nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nr
        DO nosphere=1,NSPHERES
          WRITE(OP_STRING,'('' SPHERELIST(0..,'',I3,''): '',3I2)')
     '      nosphere,(SPHERELIST(i,nosphere),i=0,SPHERELIST(0,nosphere))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !nosphere
      ENDIF

      IF(NJT.EQ.2)THEN
c cpb 9/8/95 changing mesh structures
C        FORMAT='(/$,'' Enter basis function number of the '
C     '    //'hermite BE basis [1]: '',I1)'
        IF(BEMREGION) THEN
          FORMAT='(/$,'' Enter basis function number for the '
     '      //'BE basis [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(1,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(1,1)=IDATA(1)
          nb1=MESH3_NB(1,1)
          CALL ASSERT(NBC(nb1).EQ.5,'>>Not a BE basis',ERROR,*9999)
          CALL ASSERT(NIT(nb1).EQ.1,'>>Must only have 1 xi direction',
     '      ERROR,*9999)
          CALL ASSERT(IBT(1,1,nb1).EQ.1.AND.IBT(2,1,nb1).EQ.1.OR.
     '      IBT(1,1,nb1).EQ.2,'>>Not implemented for this basis type',
     '      ERROR,*9999)
        ENDIF
        IF(FEMREGION) THEN
          FORMAT='(/$,'' Enter basis function number for the '
     '      //'FE basis [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(1,2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(1,2)=IDATA(1)
          nb2=MESH3_NB(1,2)
          CALL ASSERT(NBC(nb2).EQ.1,'>>Not a FE basis',ERROR,*9999)
          CALL ASSERT(NIT(nb2).EQ.2,
     '      '>>Must only have 2 xi directions',ERROR,*9999)
          IF(BEMREGION) THEN
            CALL ASSERT(IBT(1,1,nb1).EQ.IBT(1,1,nb2).AND.
     '        IBT(2,1,nb1).EQ.IBT(2,1,nb2),
     '        '>>Inconsistent bases across the interface',ERROR,*9999)
          ENDIF
          CALL ASSERT((IBT(1,1,nb2).EQ.1.AND.IBT(2,1,nb2).EQ.1.OR.
     '      IBT(1,1,nb2).EQ.2).AND.(IBT(1,2,nb2).EQ.1.AND.IBT(2,2,nb2)
     '      .EQ.1.OR.IBT(1,2,nb2).EQ.2),
     '      '>>Not implemented for this basis type',ERROR,*9999)
        ENDIF
C        FORMAT='(/$,'' Enter the number of circles [2]: '',I1)'
C        IDEFLT(1)=2
C        IF(IOTYPE.EQ.3) IDATA(1)=NSPHERES
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,10,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) NSPHERES=IDATA(1)

        DO nocircle=1,NSPHERES
          WRITE(CHAR1,'(I1)') nocircle
          FORMAT='(/$,'' Enter the radius of circle '//CHAR1
     '      //' ['//CHAR1//'.0]: '',D12.4)'
          RDEFLT(1)=DBLE(nocircle)
          IF(IOTYPE.EQ.3) RDATA(1)=MESH3_RAD(nocircle)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_RAD(nocircle)=RDATA(1)

          IF(nocircle.EQ.1) THEN
            IDEFLT(1)=4
          ELSE
            IDEFLT(1)=MESH3_S(nocircle-1,1)
          ENDIF
          WRITE(CHAR1B,'(I1)') IDEFLT(1)
          FORMAT='(/$,'' Enter the number of elements around circle '
     '     //CHAR1//' ['//CHAR1B//']: '',I4)'
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_S(nocircle,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_S(nocircle,1)=IDATA(1)

c cpb 4/9/95 Should use define material
C          FORMAT='(/$,'' Enter the conductivity of circle '
C     '     //CHAR1//' [1.0]: '',D12.4)'
C          RDEFLT(1)=1.0d0
C          IF(IOTYPE.EQ.3) RDATA(1)=SIGMA(nocircle)
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,20,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) SIGMA(nocircle)=RDATA(1)
        ENDDO !End of input for all spheres
        DO nr=1,NRT
          NPNODE(0,nr)=0
        ENDDO !nr
        NKTOBIG=.FALSE.
        np=1
C***    Loop over the circles
        DO nocircle=1,NSPHERES
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Generating nodes for nocircle='','
     '        //'I3)') nocircle
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          DO m1=0,MESH3_S(nocircle,1)-1
C***        Calculate nodal coordinates
            THETA=DBLE(m1)/DBLE(MESH3_S(nocircle,1))*2.0d0*PI
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' m1='',I2,'' np='',I5,'
     '          //''' theta='',D12.4)') m1,np,THETA
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            IF(np.LE.NPM) THEN
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
              NJT=2
              NKJ(1,np)=1
              NKJ(2,np)=1
              IF(ITYP10(1).EQ.1) THEN
                XP(1,1,1,np)=MESH3_RAD(nocircle)*DCOS(THETA)
                XP(1,1,2,np)=MESH3_RAD(nocircle)*DSIN(THETA)
              ELSE IF(ITYP10(1).EQ.2) THEN
                XP(1,1,1,np)=MESH3_RAD(nocircle)
                XP(1,1,2,np)=THETA
              ENDIF
              FEMCIRCLE=.FALSE.
              DO i=1,SPHERELIST(0,nocircle)
                IF(MESH3_TYPE(SPHERELIST(i,nocircle)).EQ.2) THEN
                  nb=nb2
                  FEMCIRCLE=.TRUE.
                ELSE
                  nb=nb1
                ENDIF
              ENDDO
              IF(IBT(1,1,nb).EQ.2) THEN !Hermite
C***            Calculate derivatives wrt s1
                IF(ITYP10(1).EQ.1) THEN
                  IF(NBI(nb).EQ.1) THEN !unit scale factors
                    RdTHETAdXI=MESH3_RAD(nocircle)*(-2)*PI
     '                                     /DBLE(MESH3_S(nocircle,1))
                    S1_DIR(1)=-RdTHETAdXI*DSIN(THETA)
                    S1_DIR(2)=RdTHETAdXI*DCOS(THETA)
                  ELSE
                    S1_DIR(1)=DSIN(THETA)
                    S1_DIR(2)=-DCOS(THETA)
                  ENDIF
                  IF(nocircle.EQ.1.AND.MESH3_TYPE(1).EQ.1.AND.
     '              MESH3_BEMTYPE.EQ.2) THEN
C***                s1 is in the opposite direction
                    DO nj=1,NJT
                      S1_DIR(nj)=-S1_DIR(nj)
                    ENDDO
                  ENDIF
                  NKJ(1,np)=2
                  NKJ(2,np)=2
                  IF(NKJ(1,np).LE.NKM) THEN
                    XP(2,1,1,np)=S1_DIR(1)
                    XP(2,1,2,np)=S1_DIR(2)
                  ELSE
                    NKTOBIG=.TRUE.
                  ENDIF
                ELSE IF(ITYP10(1).EQ.2) THEN
                  ERROR='>>Hermite theta interpolation not implemented'
                  GO TO 9999
C                  NKJ(1,np)=2
C                  NKJ(2,np)=2
C                  IF(NKJ(1,np).LE.NKM) THEN
CC!!! KAT 17Apr98: Something's wrong here.
C                    XP(2,1,1,np)=1.0d0
C                    XP(2,1,1,np)=0.0d0
C                  ELSE
C                    NKTOBIG=.TRUE.
C                  ENDIF
                ENDIF
              ENDIF
              IF(FEMCIRCLE) THEN
                IF(IBT(1,2,nb2).EQ.2) THEN
                  IF(NBI(nb).EQ.1) THEN !unit scale factors
                    IF(nocircle.EQ.1) THEN
                      dRdXI=MESH3_RAD(2)-MESH3_RAD(1)
                    ELSEIF(nocircle.EQ.NSPHERES) THEN
                      dRdXI=MESH3_RAD(NSPHERES)-MESH3_RAD(NSPHERES-1)
                    ELSE
C                     Probably not ideal but a simple geometric
C                     average for different sized elements.
                      dRdXI=(MESH3_RAD(nocircle+1)
     '                       -MESH3_RAD(nocircle-1))/2
                    ENDIF
                  ELSE
                    dRdXI=1.0d0
                  ENDIF
                  IF(ITYP10(1).EQ.1) THEN
                    NKJ(1,np)=NKJ(1,np)*2
                    NKJ(2,np)=NKJ(2,np)*2
                    IF(NKJ(1,np).LE.NKM) THEN
                      RbydRdXI=MESH3_RAD(nocircle)/dRdXI
                      IF(NKJ(1,np).GT.2) THEN
                        XP(3,1,1,np)=XP(1,1,1,np)/RbydRdXI
                        XP(3,1,2,np)=XP(1,1,2,np)/RbydRdXI
C KAT 17Apr98: Second derivatives are not zero.
                        XP(4,1,1,np)=XP(2,1,1,np)/RbydRdXI
                        XP(4,1,2,np)=XP(2,1,2,np)/RbydRdXI
                      ELSE
                        XP(2,1,1,np)=XP(1,1,1,np)/RbydRdXI
                        XP(2,1,2,np)=XP(1,1,2,np)/RbydRdXI
                      ENDIF
                    ELSE
                      NKTOBIG=.TRUE.
                    ENDIF
                  ELSE IF(ITYP10(1).EQ.2) THEN
                    NKJ(1,np)=NKJ(1,np)*2
                    NKJ(2,np)=NKJ(2,np)*2
                    IF(NKJ(1,np).LE.NKM) THEN
                      IF(NKJ(1,np).GT.2) THEN
                        XP(3,1,1,np)=0.0d0
                        XP(3,1,2,np)=dRdXI
                        XP(4,1,1,np)=0.0d0
                        XP(4,1,2,np)=0.0d0
                      ELSE
                        XP(2,1,1,np)=0.0d0
                        XP(2,1,2,np)=dRdXI
                      ENDIF
                    ELSE
                      NKTOBIG=.TRUE.
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
              IF(DOP) THEN
                WRITE(OP_STRING,'('' NKJ(1..,np)= '',2I2)') NKJ(1,np),
     '            NKJ(2,np)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' XP(nk,1,1,np):'',4D12.4)')
     '            (XP(nk,1,1,np),nk=1,NKJ(1,np))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,'('' XP(nk,1,2,np):'',4D12.4)')
     '            (XP(nk,1,2,np),nk=1,NKJ(2,np))
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDIF
            np=np+1
          ENDDO !m1
          np=np-1
          CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
          CALL ASSERT(.NOT.NKTOBIG,'>>Increase NKM',ERROR,*9999)
C***      Loop over regions associated with the circle
          DO nonr=1,SPHERELIST(0,nocircle)
            nr=SPHERELIST(nonr,nocircle)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' nr='',I2)') nr
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NPT(nr)=np
            IF(nr.EQ.1) THEN
              NPNODE(0,nr)=NPT(nr)
              IF(NPNODE(0,nr).LE.NP_R_M) THEN
                DO nonode=1,NPT(nr)
                  NPNODE(nonode,nr)=nonode
                ENDDO !nonode
              ENDIF
            ELSE
              NPSTART=0
              DO ncirc=1,nocircle-1
                NPSTART=NPSTART+MESH3_S(ncirc,1)
              ENDDO
              DO nonode=1,NPT(nr)-NPSTART
                NPNODE(nonode+NPNODE(0,nr),nr)=NPSTART+nonode
              ENDDO
              NPNODE(0,nr)=NPNODE(0,nr)+NPT(nr)-NPSTART
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' NPT(nr)='',I5)') NPT(nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NPNODE(0,nr)='',I5)') NPNODE(0,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NPNODE(1..,nr)='',10(1X,I5),'
     '          //'/:(10(1X,I5)))') (NPNODE(nonode,nr),nonode=1,
     '          NPNODE(0,nr))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nonr (nr)
          np=np+1
        ENDDO !nocircle - End of node construction
        np=np-1
        CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
        NPT(0)=np
        NPNODE(0,0)=np
        CALL ASSERT(NPNODE(0,0).LE.NP_R_M,'>>Increase NP_R_M',
     '    ERROR,*9999)
C***    Construct elements and npne array
        ne=1
        DO nr=1,NRT
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Generating elements for region '','
     '        //'I2,'':'')') nr
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(MESH3_TYPE(nr).EQ.1) THEN !BEM
            DO ncirc=1,REGLIST(0,nr)
              nocircle=REGLIST(ncirc,nr)
              NPSTART=0
              DO i=1,nocircle-1
                NPSTART=NPSTART+MESH3_S(i,1)
              ENDDO !i
              IF(DOP) THEN
                WRITE(OP_STRING,'('' ncirc='',I1,'' nocircle='',I3,'
     '            //''' npstart='',I5)') ncirc,nocircle,npstart
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              DO m1=1,MESH3_S(nocircle,1) !theta steps
                np=m1+NPSTART
                IF(ne.LE.NEM) THEN
C LKC for AJP 191297
CC AJPs
                  NJ_LOC(NJL_GEOM,0,nr)=2
CC AJPe

                  NBJ(1,ne)=nb1
                  NBJ(2,ne)=nb1
                  IF(m1.EQ.MESH3_S(nocircle,1))THEN
                    NP_INCR=1-MESH3_S(nocircle,1)
                  ELSE
                    NP_INCR=1
                  ENDIF
                  IF(nocircle.EQ.1.AND.MESH3_BEMTYPE.EQ.2) THEN
C***                Need to have s1 direction reversed.
C***                Number anticlockwise in this case.
                    NPNE(1,nb1,ne)=np
                    NPNE(2,nb1,ne)=np+NP_INCR
                  ELSE
C***                Number clockwise in this case.
                    NPNE(1,nb1,ne)=np+NP_INCR
                    NPNE(2,nb1,ne)=np
                  ENDIF
                  IF(DOP) THEN
                    WRITE(OP_STRING,'('' ne='',I5,'' m1='',I3,'' np='','
     '                //'I5,'' np_incr='',I5)') ne,m1,np,np_incr
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne):'',2I6)')
     '                nb1,(NPNE(nn,nb1,ne),nn=1,2)
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDIF
                ne=ne+1


              ENDDO !end m1
            ENDDO !ncirc
          ELSE IF(MESH3_TYPE(nr).EQ.2) THEN !FEM
            nocircle=REGLIST(1,nr)
            NPSTART=0
            DO i=1,nocircle-1
              NPSTART=NPSTART+MESH3_S(i,1)
            ENDDO !i
            IF(DOP) THEN
              WRITE(OP_STRING,'('' nocircle='',I3,'' npstart='',I5)')
     '          nocircle,npstart
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO m1=1,MESH3_S(nocircle,1) !theta steps
              np=m1+NPSTART
              IF(ne.LE.NEM) THEN
C LKC for AJP 191297
CC AJPs
                NJ_LOC(NJL_GEOM,0,nr)=2
C                  NJE(ne)=2
CC AJPe
                NBJ(1,ne)=nb2
                NBJ(2,ne)=nb2
                IF(m1.EQ.MESH3_S(nocircle,1))THEN
                  NP_INCR=1-MESH3_S(nocircle,1)
                ELSE
                  NP_INCR=1
                ENDIF
                IF(ITYP10(nr).EQ.1) THEN
C***              Number clockwise in this case.
                  NPNE(1,nb2,ne)=np+NP_INCR
                  NPNE(2,nb2,ne)=np
                  NPNE(3,nb2,ne)=np+NP_INCR+MESH3_S(nocircle,1)
                  NPNE(4,nb2,ne)=np+MESH3_S(nocircle,1)
                ELSE IF(ITYP10(nr).EQ.2) THEN
C***              Number anticlockwise in this case.
                  NPNE(1,nb2,ne)=np
                  NPNE(2,nb2,ne)=np+NP_INCR
                  NPNE(3,nb2,ne)=np+MESH3_S(nocircle,1)
                  NPNE(4,nb2,ne)=np+NP_INCR+MESH3_S(nocircle,1)
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ne='',I5,'' m1='',I3,'' np='','
     '              //'I5,'' np_incr='',I5)') ne,m1,np,np_incr
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne):'',4I6)')
     '              nb2,(NPNE(nn,nb2,ne),nn=1,4)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              ne=ne+1
            ENDDO !end m1
          ENDIF
          ne=ne-1
          CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
          NET(nr)=ne
          IF(nr.EQ.1) THEN
            NEELEM(0,nr)=NET(nr)
            IF(NEELEM(0,nr).LE.NE_R_M) THEN
              DO noelem=1,NET(nr)
                NEELEM(noelem,nr)=noelem
                NRE(noelem)=1
              ENDDO !noelem
            ENDIF
          ELSE
            IF(NET(nr)-NET(nr-1).LE.NE_R_M) THEN
              DO noelem=1,NET(nr)-NET(nr-1)
                NEELEM(noelem,nr)=NET(nr-1)+noelem
                NRE(NET(nr-1)+noelem)=nr
              ENDDO
            ENDIF
            NEELEM(0,nr)=NET(nr)-NET(nr-1)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' NET(nr)='',I5)') NET(nr)
            WRITE(OP_STRING,'('' NEELEM(0,nr)='',I5)') NEELEM(0,nr)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NEELEM(1..,nr)='',10(1X,I5),'
     '        //'/:(10(1X,I5)))') (NEELEM(noelem,nr),noelem=1,
     '        NEELEM(0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NRE(ne)='',30(1X,I1),'
     '        //'/:(30(1X,I1)))') (NRE(NEELEM(noelem,nr)),noelem=1,
     '        NEELEM(0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          ne=ne+1
        ENDDO !nr
        ne=ne-1
        CALL ASSERT(ne.LE.NEM,'>>Increase NEM',ERROR,*9999)
        NET(0)=ne
        NEELEM(0,0)=ne !number of elements in all regions
        CALL ASSERT(NEELEM(0,0).LE.NE_R_M,'>>Increase NE_R_M',
     '    ERROR,*9999)

      ELSE IF(NJT.EQ.3)THEN
c cpb 10/8/95 Adding FE spheres etc

C cpb 2/11/96 Adding BE cubic and linear sectors
C CPB 10/8/95 Only implementing spheres for bicubic Hermite (BE)
C and bicubic Hermite Linear (FE) bases.

C        CALL ASSERT(NBFT.GT.2,'>>Define three basis functions',
C     '    ERROR,*9999)

        IF(BEMREGION) THEN
C          FORMAT='(/$,'' Enter basis number of the bicubic '
C     '      //'hermite BE basis [1]: '',I1)'
C         At least 3 boundary element basis functions should have
C         already been set up  - 2 hermite simplexes and one bicubic
C         hermite.
          FORMAT='(/$,'' Enter basis number of the main (central) '
     '      //'BE basis [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(1,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(1,1)=IDATA(1)
          nb1=MESH3_NB(1,1)
          CALL ASSERT(NBC(nb1).EQ.5,'>>Not a BE basis',
     '      ERROR,*9999)
          CALL ASSERT(NIT(nb1).EQ.2,'>>Must have 2 xi directions',
     '      ERROR,*9999)
          CALL ASSERT((IBT(1,1,nb1).EQ.1.AND.IBT(2,1,nb1).EQ.1.AND.
     '      IBT(1,2,nb1).EQ.1.AND.IBT(2,2,nb1).EQ.1).OR.
     '      (IBT(1,1,nb1).EQ.2.AND.IBT(2,1,nb1).EQ.1.AND.
     '      IBT(1,2,nb1).EQ.2.AND.IBT(2,2,nb1).EQ.1),
     '      '>>Not implemented for this basis type',ERROR,*9999)
C          FORMAT='(/$,'' Enter basis number of the hermite '
C     '      //'simplex BE basis (apex at node 1) [2]: '',I1)'
          FORMAT='(/$,'' Enter basis number of the bottom '
     '      //'BE sector basis (apex at node 1) [2]: '',I2)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(2,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(2,1)=IDATA(1)
          CALL ASSERT(NBC(MESH3_NB(2,1)).EQ.6,'>>Not a BE sector basis',
     '      ERROR,*9999)

C          FORMAT='(/$,'' Enter basis number of the hermite '
C     '      //'simplex BE basis (apex at node 3) [3]: '',I1)'
          FORMAT='(/$,'' Enter basis number of the top '
     '      //'BE sector basis (apex at node 3) [3]: '',I2)'
          IDEFLT(1)=3
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(3,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(3,1)=IDATA(1)
          CALL ASSERT(NBC(MESH3_NB(3,1)).EQ.6,'>>Not a BE sector basis',
     '      ERROR,*9999)
        ENDIF
        IF(FEMREGION) THEN
          FORMAT='(/$,'' Enter basis number of the main (central) '
     '      //'FE basis [1]: '',I2)'
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(1,2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IONE,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(1,2)=IDATA(1)
          nb2=MESH3_NB(1,2)
          CALL ASSERT(NBC(nb2).EQ.1,'>>Not a FE basis',
     '      ERROR,*9999)
          CALL ASSERT(NIT(nb2).EQ.3,'>>Must have 3 xi directions',
     '      ERROR,*9999)
          CALL ASSERT(IBT(1,1,nb2).EQ.2.AND.IBT(2,1,nb2).EQ.1.AND.
     '      IBT(1,2,nb2).EQ.2.AND.IBT(2,2,nb2).EQ.1.AND.
     '      IBT(1,3,nb2).EQ.1.AND.IBT(2,3,nb2).EQ.1,
     '      '>>Not implemented for this basis type',ERROR,*9999)

          FORMAT='(/$,'' Enter basis number of the bottom FE sector '
     '      //'basis (apex node 1) [2]: '',I2)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(2,2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(2,2)=IDATA(1)
          nb=MESH3_NB(2,2)
          CALL ASSERT(NBC(nb).EQ.2,'>>Not a FE sector basis',
     '      ERROR,*9999)
          CALL ASSERT(NIT(nb).EQ.3,'>>Must have 3 xi directions',
     '      ERROR,*9999)
          CALL ASSERT(IBT(1,1,nb).EQ.5.AND.IBT(2,1,nb).EQ.4.AND.
     '      IBT(3,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2.AND.
     '      IBT(2,2,nb).EQ.2.AND.IBT(1,3,nb).EQ.1.AND.
     '      IBT(2,3,nb).EQ.1,
     '      '>>Not implemented for this basis type',ERROR,*9999)

          FORMAT='(/$,'' Enter basis number of the top FE sector '
     '      //'FE basis (apex node 3) [3]: '',I2)'
          IDEFLT(1)=3
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_NB(3,2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NBT,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_NB(3,2)=IDATA(1)
          nb=MESH3_NB(3,2)
          CALL ASSERT(NBC(nb).EQ.2,'>>Not a FE sector basis',
     '      ERROR,*9999)
          CALL ASSERT(NIT(nb).EQ.3,'>>Must have 3 xi directions',
     '      ERROR,*9999)
          CALL ASSERT(IBT(1,1,nb).EQ.6.AND.IBT(2,1,nb).EQ.4.AND.
     '      IBT(3,1,nb).EQ.2.AND.IBT(1,2,nb).EQ.2.AND.
     '      IBT(2,2,nb).EQ.3.AND.IBT(1,3,nb).EQ.1.AND.
     '      IBT(2,3,nb).EQ.1,
     '      '>>Not implemented for this basis type',ERROR,*9999)
        ENDIF

C CPB 2/11/96 Check if there are any derivatives to calculate. In the
C future will be able to have cubics and/or linears in each direction
C but for now assume bilinear or bicubic

        IF(FEMREGION) THEN
          nb=MESH3_NB(1,2)
        ELSE
          nb=MESH3_NB(1,1)
        ENDIF
        CALC_DERIV=.FALSE.
        DO ni=1,NIT(nb)
          IF(IBT(1,ni,nb).EQ.2.OR.IBT(2,ni,nb).EQ.4) CALC_DERIV=.TRUE.
        ENDDO !ni

C        FORMAT='(/$,'' Enter the number of spheres [2]: '',I1)'
C        IDEFLT(1)=2
C        IF(IOTYPE.EQ.3) IDATA(1)=NSPHERES
C        CALL GINOUT(IOTYPE,3,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
C     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C        IF(IOTYPE.NE.3) NSPHERES=IDATA(1)

        DO nosphere=1,NSPHERES
          WRITE(CHAR1,'(I1)') nosphere
          FORMAT='(/$,'' Enter the radius of sphere '//CHAR1
     '      //' ['//CHAR1//'.0]: '',D12.4)'
          RDEFLT(1)=DBLE(nosphere)
          IF(IOTYPE.EQ.3) RDATA(1)=MESH3_RAD(nosphere)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_RAD(nosphere)=RDATA(1)

          FORMAT='(/$,'' Enter the number of elements around sphere '
     '     //CHAR1//' (theta dir) [4]: '',I4)'
          IDEFLT(1)=4
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_S(nosphere,1)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_S(nosphere,1)=IDATA(1)

          FORMAT='(/$,'' Enter the number of elements in the '
     '    //'azimuthal direction of sphere '//CHAR1//' [2]: '',I4)'
          IDEFLT(1)=2
          IF(IOTYPE.EQ.3) IDATA(1)=MESH3_S(nosphere,2)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,2,9999,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) MESH3_S(nosphere,2)=IDATA(1)

c cpb 4/9/95 Should use define material
C          FORMAT='(/$,'' Enter the conductivity of sphere '
C     '     //CHAR1//' [1.0]: '',D12.4)'
C          RDEFLT(1)=1.0d0
C          IF(IOTYPE.EQ.3) RDATA(1)=SIGMA(nosphere)
C          CALL GINOUT(IOTYPE,5,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
C     '      ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,20,
C     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
C          IF(IOTYPE.NE.3) SIGMA(nosphere)=RDATA(1)
        ENDDO !End of input for all spheres

        DO nr=1,NRT
          NPNODE(0,nr)=0
        ENDDO !nr
        IF(CALC_DERIV) CALL ASSERT(NKM.GE.4,'>>Increase NKM',
     '    ERROR,*9999)
        np=1
C***    Loop over the spheres
        DO nosphere=1,NSPHERES
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Generating nodes for nosphere='','
     '        //'I3)') nosphere
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
C***      Calculate nodal coordinates
C***      Calculate the bottom node first
          IF(np.LE.NPM) THEN
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
            NJT=3
            NKJ(1,np)=1
            NKJ(2,np)=1
            NKJ(3,np)=1
            XP(1,1,1,np)=0.0d0
            XP(1,1,2,np)=0.0d0
            XP(1,1,3,np)=-MESH3_RAD(nosphere)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Bottom node. np='',I5)') np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NKJ(1..,np)= '',3I2)') NKJ(1,np),
     '          NKJ(2,np),NKJ(3,np)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,1,np):'',4D12.4)')
     '          (XP(nk,1,1,np),nk=1,NKJ(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,2,np):'',4D12.4)')
     '          (XP(nk,1,2,np),nk=1,NKJ(2,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,3,np):'',4D12.4)')
     '          (XP(nk,1,3,np),nk=1,NKJ(3,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
          np=np+1
C***      Calculate the middle nodes
          DO m2=1,MESH3_S(nosphere,2)-1
            DO m1=0,MESH3_S(nosphere,1)-1
              PHI=PI-DBLE(m2)/DBLE(MESH3_S(nosphere,2))*PI
              THETA=DBLE(m1)/DBLE(MESH3_S(nosphere,1))*2.0d0*PI
              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' m1='',I2,'' m2='',I2,'' np='',I5,'
     '            //''' theta='',D12.4,'' phi='',D12.4)') m1,m2,np,
     '            THETA,PHI
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(np.LE.NPM) THEN
C GMH 8/1/97 Update cmgui link
                CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
                XP(1,1,1,np)=MESH3_RAD(nosphere)*DSIN(PHI)*DCOS(THETA)
                XP(1,1,2,np)=MESH3_RAD(nosphere)*DSIN(PHI)*DSIN(THETA)
                XP(1,1,3,np)=MESH3_RAD(nosphere)*DCOS(PHI)
                NJT=3
                IF(CALC_DERIV) THEN
C***            Calculate derivatives wrt s1
                  S1_DIR(1)=-DSIN(THETA)
                  S1_DIR(2)=DCOS(THETA)
                  S1_DIR(3)=0.0d0
C***            Calculate derivatives wrt s2
                  S2_DIR(1)=-DCOS(PHI)*DCOS(THETA)
                  S2_DIR(2)=-DCOS(PHI)*DSIN(THETA)
                  S2_DIR(3)=DSIN(PHI)
                  IF(nosphere.EQ.1.AND.NSPHERES.GT.1.AND.
     '              MESH3_TYPE(1).EQ.1.AND.MESH3_BEMTYPE.EQ.2) THEN
C***              s2 is in the opposite direction
                    S2_DIR(1)=-S2_DIR(1)
                    S2_DIR(2)=-S2_DIR(2)
                    S2_DIR(3)=-S2_DIR(3)
                  ENDIF
                  XP(2,1,1,np)=S1_DIR(1)
                  XP(3,1,1,np)=S2_DIR(1)
                  XP(4,1,1,np)=0.0d0
                  XP(2,1,2,np)=S1_DIR(2)
                  XP(3,1,2,np)=S2_DIR(2)
                  XP(4,1,2,np)=0.0d0
                  XP(2,1,3,np)=S1_DIR(3)
                  XP(3,1,3,np)=S2_DIR(3)
                  XP(4,1,3,np)=0.0d0
                  NKJ(1,np)=4
                  NKJ(2,np)=4
                  NKJ(3,np)=4
                ELSE
                  NKJ(1,np)=1
                  NKJ(2,np)=1
                  NKJ(3,np)=1
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' NKJ(1..,np)= '',3I2)') NKJ(1,np),
     '              NKJ(2,np),NKJ(3,np)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' XP(nk,1,1,np):'',4D12.4)')
     '              (XP(nk,1,1,np),nk=1,NKJ(1,np))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' XP(nk,1,2,np):'',4D12.4)')
     '              (XP(nk,1,2,np),nk=1,NKJ(2,np))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' XP(nk,1,3,np):'',4D12.4)')
     '              (XP(nk,1,3,np),nk=1,NKJ(3,np))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
              ENDIF
              np=np+1
            ENDDO !m1
          ENDDO !m2
C***      Calculate the top node
          IF(np.LE.NPM) THEN
C GMH 8/1/97 Update cmgui link
            CALL NODE_CHANGE(np,.TRUE.,ERROR,*9999)
            NJT=3
            NKJ(1,np)=1
            NKJ(2,np)=1
            NKJ(3,np)=1
            XP(1,1,1,np)=0.0d0
            XP(1,1,2,np)=0.0d0
            XP(1,1,3,np)=MESH3_RAD(nosphere) !Top node
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' Top node. np='',I5)') np
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NKJ(1..,np)= '',3I2)') NKJ(1,np),
     '          NKJ(2,np),NKJ(3,np)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,1,np):'',4D12.4)')
     '          (XP(nk,1,1,np),nk=1,NKJ(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,2,np):'',4D12.4)')
     '          (XP(nk,1,2,np),nk=1,NKJ(2,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' XP(nk,1,3,np):'',4D12.4)')
     '          (XP(nk,1,3,np),nk=1,NKJ(3,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF
C***      Loop over regions associated with the circle
          DO nonr=1,SPHERELIST(0,nosphere)
            nr=SPHERELIST(nonr,nosphere)
            IF(DOP) THEN
              WRITE(OP_STRING,'(/'' nr='',I2)') nr
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            NPT(nr)=np
            IF(nr.EQ.1) THEN
              NPNODE(0,nr)=NPT(nr)
              IF(NPNODE(0,nr).LE.NP_R_M) THEN
                DO nonode=1,NPT(nr)
                  NPNODE(nonode,nr)=nonode
                ENDDO !nonode
              ENDIF
            ELSE
              NPSTART=0
              DO nspher=1,nosphere-1
                NPSTART=NPSTART+(MESH3_S(nspher,2)-1)*
     '            MESH3_S(nspher,1)+2
              ENDDO !nspher
              IF(NPNODE(0,nr)+NPT(nr)-NPSTART.LE.NP_R_M) THEN
                DO nonode=1,NPT(nr)-NPSTART
                  NPNODE(nonode+NPNODE(0,nr),nr)=NPSTART+nonode
                ENDDO
              ENDIF
              NPNODE(0,nr)=NPNODE(0,nr)+NPT(nr)-NPSTART
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' NPT(nr)='',I5)') NPT(nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NPNODE(0,nr)='',I5)') NPNODE(0,nr)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' NPNODE(1..,nr)='',10(1X,I5),'
     '          //'/:(10(1X,I5)))') (NPNODE(nonode,nr),nonode=1,
     '          NPNODE(0,nr))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO !nonr (nr)
          np=np+1
        ENDDO !nosphere - End of node construction
        np=np-1
        CALL ASSERT(np.LE.NPM,'>>Increase NPM',ERROR,*9999)
        NPT(0)=np
        NPNODE(0,0)=np
        CALL ASSERT(NPNODE(0,0).LE.NP_R_M,'>>Increase NP_R_M',
     '    ERROR,*9999)
C***    Construct elements and npne array
        ne=1
        DO nr=1,NRT
          IF(DOP) THEN
            WRITE(OP_STRING,'(/'' Generating elements for region '','
     '        //'I2,'':'')') nr
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          IF(MESH3_TYPE(nr).EQ.1) THEN
            DO nspher=1,REGLIST(0,nr)
              nosphere=REGLIST(nspher,nr)
              NPSTART=0
              DO i=1,nosphere-1
                NPSTART=NPSTART+(MESH3_S(i,2)-1)*MESH3_S(i,1)+2
              ENDDO !i
              NP_INCR=(MESH3_S(nosphere,2)-1)*MESH3_S(nosphere,1)+2
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nspher='',I1,'' nosphere='',I3,'
     '            //''' npstart='',I5,'' np_incr='',I5)') nspher,
     '            nosphere,npstart,np_incr
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              np=NPSTART+1
              IF(DOP) THEN
                WRITE(OP_STRING,'('' nspher='',I1,'' nosphere='',I3,'
     '            //''' npstart='',I5)') nspher,nosphere,npstart
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF
              IF(nosphere.EQ.1.AND.MESH3_BEMTYPE.EQ.2) THEN
C***            Reverse the direction of s2
                DO m2=1,MESH3_S(nosphere,2) !phi steps
                  DO m1=1,MESH3_S(nosphere,1) !theta steps
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' ne='',I5,'' m1='',I3,'
     '                  //''' m2='',I3,'' np='',I5)') ne,m1,m2,np
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
C LKC for AJP
CC AJPs 191297
C                    NJE(ne)=3
                    NJ_LOC(NJL_GEOM,0,nr)=3
CC AJPe
                    IF(m2.EQ.1) THEN
C***                  Top Elements
                      NBJ(1,ne)=MESH3_NB(3,1)
                      NBJ(2,ne)=MESH3_NB(3,1)
                      NBJ(3,ne)=MESH3_NB(3,1)
                      NPNE(1,MESH3_NB(3,1),ne)=np+m1
                      NPNE(2,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(3,1),ne)=np
                      NPNE(1,MESH3_NB(1,1),ne)=np+m1
                      NPNE(2,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(1,1),ne)=np
                      NPNE(4,MESH3_NB(1,1),ne)=np
                      NPNE(1,MESH3_NB(2,1),ne)=np+m1
                      NPNE(2,MESH3_NB(2,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(2,1),ne)=np
                    ELSE IF(m2.EQ.MESH3_S(nosphere,2)) THEN
C***                  Bottom Elements
                      NBJ(1,ne)=MESH3_NB(2,1)
                      NBJ(2,ne)=MESH3_NB(2,1)
                      NBJ(3,ne)=MESH3_NB(2,1)
                      NPNE(1,MESH3_NB(2,1),ne)=npstart+np_incr
                      NPNE(2,MESH3_NB(2,1),ne)=np+m1
                      NPNE(3,MESH3_NB(2,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(1,1),ne)=npstart+np_incr
                      NPNE(2,MESH3_NB(1,1),ne)=npstart+np_incr
                      NPNE(3,MESH3_NB(1,1),ne)=np+m1
                      NPNE(4,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(3,1),ne)=npstart+np_incr
                      NPNE(2,MESH3_NB(3,1),ne)=np+m1
                      NPNE(3,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                    ELSE
C***                  Middle Elements
                      NBJ(1,ne)=MESH3_NB(1,1)
                      NBJ(2,ne)=MESH3_NB(1,1)
                      NBJ(3,ne)=MESH3_NB(1,1)
                      NPNE(1,MESH3_NB(1,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                      NPNE(2,MESH3_NB(1,1),ne)=np+MESH3_S(nosphere,1)+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(1,1),ne)=np+m1
                      NPNE(4,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(2,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                      NPNE(2,MESH3_NB(2,1),ne)=np+
     '                  MESH3_S(nosphere,1)+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(2,1),ne)=np+m1
                      NPNE(1,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(2,MESH3_NB(3,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                      NPNE(3,MESH3_NB(3,1),ne)=np+m1
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NBJ(1..,ne):'',3(1X,I2))')
     '                  (NBJ(nj,ne),nj=1,3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne):'','
     '                  //'4(1X,I5))') NBJ(1,ne),(NPNE(nn,NBJ(1,ne),ne),
     '                  nn=1,NNT(NBJ(1,ne)))
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
                    ne=ne+1
                  ENDDO !m1
                  IF(m2.NE.1) THEN !Not the Bottom element
                    np=np+MESH3_S(nosphere,1)
                  ENDIF
                ENDDO !m2
              ELSE
                DO m2=1,MESH3_S(nosphere,2) !phi steps
                  DO m1=1,MESH3_S(nosphere,1) !theta steps
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' ne='',I5,'' m1='',I3,'
     '                  //''' m2='',I3,'' np='',I5)') ne,m1,m2,np
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF
C LKC for AJP
CC AJPs 191297
                    NJ_LOC(NJL_GEOM,0,nr)=3
C                    NJE(ne)=3
CC AJPe
                    IF(m2.EQ.1) THEN
C***                  Bottom Elements
                      NBJ(1,ne)=MESH3_NB(2,1)
                      NBJ(2,ne)=MESH3_NB(2,1)
                      NBJ(3,ne)=MESH3_NB(2,1)
                      NPNE(1,MESH3_NB(2,1),ne)=np
                      NPNE(2,MESH3_NB(2,1),ne)=np+m1
                      NPNE(3,MESH3_NB(2,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(1,1),ne)=np
                      NPNE(2,MESH3_NB(1,1),ne)=np
                      NPNE(3,MESH3_NB(1,1),ne)=np+m1
                      NPNE(4,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(3,1),ne)=np
                      NPNE(2,MESH3_NB(3,1),ne)=np+m1
                      NPNE(3,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                    ELSE IF(m2.EQ.MESH3_S(nosphere,2)) THEN
C***                  Top Elements
                      NBJ(1,ne)=MESH3_NB(3,1)
                      NBJ(2,ne)=MESH3_NB(3,1)
                      NBJ(3,ne)=MESH3_NB(3,1)
                      NPNE(1,MESH3_NB(3,1),ne)=np+m1
                      NPNE(2,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(3,1),ne)=npstart+np_incr
                      NPNE(1,MESH3_NB(1,1),ne)=np+m1
                      NPNE(2,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(1,1),ne)=npstart+np_incr
                      NPNE(4,MESH3_NB(1,1),ne)=npstart+np_incr
                      NPNE(1,MESH3_NB(2,1),ne)=np+m1
                      NPNE(2,MESH3_NB(2,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(2,1),ne)=npstart+np_incr
                    ELSE
C***                  Middle Elements
                      NBJ(1,ne)=MESH3_NB(1,1)
                      NBJ(2,ne)=MESH3_NB(1,1)
                      NBJ(3,ne)=MESH3_NB(1,1)
                      NPNE(1,MESH3_NB(1,1),ne)=np+m1
                      NPNE(2,MESH3_NB(1,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(1,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                      NPNE(4,MESH3_NB(1,1),ne)=np+MESH3_S(nosphere,1)+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(2,1),ne)=np+m1
                      NPNE(2,MESH3_NB(2,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                      NPNE(3,MESH3_NB(2,1),ne)=np+
     '                  MESH3_S(nosphere,1)+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(1,MESH3_NB(3,1),ne)=np+m1
                      NPNE(2,MESH3_NB(3,1),ne)=np+
     '                  MOD(m1,MESH3_S(nosphere,1))+1
                      NPNE(3,MESH3_NB(3,1),ne)=np+m1+
     '                  MESH3_S(nosphere,1)
                    ENDIF
                    IF(DOP) THEN
                      WRITE(OP_STRING,'('' NBJ(1..,ne):'',3(1X,I2))')
     '                  (NBJ(nj,ne),nj=1,3)
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                      WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne):'','
     '                  //'4(1X,I5))') NBJ(1,ne),(NPNE(nn,NBJ(1,ne),ne),
     '                  nn=1,NNT(NBJ(1,ne)))
                      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                    ENDIF

C CPB 23-jan-01 Needs to be NEM+1
CC LKC 13-MAR-2000 New assert
C                    CALL ASSERT(NEM.GE.ne+1,'>> Increase NEM',
C     '                ERROR,*9999)
                    CALL ASSERT(NEM+1.GE.ne+1,'>> Increase NEM',
     '                ERROR,*9999)
                    ne=ne+1
                  ENDDO !m1
                  IF(m2.NE.1) THEN !Not the Bottom element
                    np=np+MESH3_S(nosphere,1)
                  ENDIF
                ENDDO !m2
              ENDIF
            ENDDO !nspher
          ELSE IF(MESH3_TYPE(nr).EQ.2) THEN
            nosphere=REGLIST(1,nr)
            NPSTART=0
            DO i=1,nosphere-1
              NPSTART=NPSTART+(MESH3_S(i,2)-1)*MESH3_S(i,1)+2
            ENDDO !i
            NP_INCR=(MESH3_S(nosphere,2)-1)*MESH3_S(nosphere,1)+2
            IF(DOP) THEN
              WRITE(OP_STRING,'('' nosphere='',I3,'' npstart='',I5,'
     '          //''' np_incr='',I5)') nosphere,npstart,np_incr
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            np=NPSTART+1
            DO m2=1,MESH3_S(nosphere,2) !phi steps
              DO m1=1,MESH3_S(nosphere,1) !theta steps
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' ne='',I5,'' m1='',I3,'' m2='','
     '              //'I3,'' np='',I5)') ne,m1,m2,np
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                IF(m2.EQ.1) THEN
C***              Bottom Elements
                  NBJ(1,ne)=MESH3_NB(2,2)
                  NBJ(2,ne)=MESH3_NB(2,2)
                  NBJ(3,ne)=MESH3_NB(2,2)
                  NPNE(1,MESH3_NB(2,2),ne)=np
                  NPNE(2,MESH3_NB(2,2),ne)=np+m1
                  NPNE(3,MESH3_NB(2,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(4,MESH3_NB(2,2),ne)=np+np_incr
                  NPNE(5,MESH3_NB(2,2),ne)=np+m1+np_incr
                  NPNE(6,MESH3_NB(2,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(1,MESH3_NB(1,2),ne)=np
                  NPNE(2,MESH3_NB(1,2),ne)=np
                  NPNE(3,MESH3_NB(1,2),ne)=np+m1
                  NPNE(4,MESH3_NB(1,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(5,MESH3_NB(1,2),ne)=np+np_incr
                  NPNE(6,MESH3_NB(1,2),ne)=np+np_incr
                  NPNE(7,MESH3_NB(1,2),ne)=np+m1+np_incr
                  NPNE(8,MESH3_NB(1,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(1,MESH3_NB(3,2),ne)=np
                  NPNE(2,MESH3_NB(3,2),ne)=np+m1
                  NPNE(3,MESH3_NB(3,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(4,MESH3_NB(3,2),ne)=np+np_incr
                  NPNE(5,MESH3_NB(3,2),ne)=np+m1+np_incr
                  NPNE(6,MESH3_NB(3,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                ELSE IF(m2.EQ.MESH3_S(nosphere,2)) THEN
C***              Top Elements
                  NBJ(1,ne)=MESH3_NB(3,2)
                  NBJ(2,ne)=MESH3_NB(3,2)
                  NBJ(3,ne)=MESH3_NB(3,2)
                  NPNE(1,MESH3_NB(3,2),ne)=np+m1
                  NPNE(2,MESH3_NB(3,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(3,MESH3_NB(3,2),ne)=npstart+np_incr
                  NPNE(4,MESH3_NB(3,2),ne)=np+np_incr+m1
                  NPNE(5,MESH3_NB(3,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(6,MESH3_NB(3,2),ne)=npstart+2*np_incr
                  NPNE(1,MESH3_NB(1,2),ne)=np+m1
                  NPNE(2,MESH3_NB(1,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(3,MESH3_NB(1,2),ne)=npstart+np_incr
                  NPNE(4,MESH3_NB(1,2),ne)=npstart+np_incr
                  NPNE(5,MESH3_NB(1,2),ne)=np+np_incr+m1
                  NPNE(6,MESH3_NB(1,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(7,MESH3_NB(1,2),ne)=npstart+2*np_incr
                  NPNE(8,MESH3_NB(1,2),ne)=npstart+2*np_incr
                  NPNE(1,MESH3_NB(2,2),ne)=np+m1
                  NPNE(2,MESH3_NB(2,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(3,MESH3_NB(2,2),ne)=npstart+np_incr
                  NPNE(4,MESH3_NB(2,2),ne)=np+np_incr+m1
                  NPNE(5,MESH3_NB(2,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(6,MESH3_NB(2,2),ne)=npstart+2*np_incr
               ELSE
C***              Middle Elements
                  NBJ(1,ne)=MESH3_NB(1,2)
                  NBJ(2,ne)=MESH3_NB(1,2)
                  NBJ(3,ne)=MESH3_NB(1,2)
                  NPNE(1,MESH3_NB(1,2),ne)=np+m1
                  NPNE(2,MESH3_NB(1,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(3,MESH3_NB(1,2),ne)=np+m1+MESH3_S(nosphere,1)
                  NPNE(4,MESH3_NB(1,2),ne)=np+MESH3_S(nosphere,1)+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(5,MESH3_NB(1,2),ne)=np+m1+np_incr
                  NPNE(6,MESH3_NB(1,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(7,MESH3_NB(1,2),ne)=np+m1+MESH3_S(nosphere,1)+
     '              np_incr
                  NPNE(8,MESH3_NB(1,2),ne)=np+MESH3_S(nosphere,1)+
     '              MOD(m1,MESH3_S(nosphere,1))+1+np_incr
                  NPNE(1,MESH3_NB(2,2),ne)=np+m1
                  NPNE(2,MESH3_NB(2,2),ne)=np+m1+MESH3_S(nosphere,1)
                  NPNE(3,MESH3_NB(2,2),ne)=np+MESH3_S(nosphere,1)+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(4,MESH3_NB(2,2),ne)=np+m1+np_incr
                  NPNE(5,MESH3_NB(2,2),ne)=np+m1+MESH3_S(nosphere,1)+
     '              np_incr
                  NPNE(6,MESH3_NB(2,2),ne)=np+MESH3_S(nosphere,1)+
     '              MOD(m1,MESH3_S(nosphere,1))+1+np_incr
                  NPNE(1,MESH3_NB(3,2),ne)=np+m1
                  NPNE(2,MESH3_NB(3,2),ne)=np+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(3,MESH3_NB(3,2),ne)=np+m1+MESH3_S(nosphere,1)
                  NPNE(4,MESH3_NB(3,2),ne)=np+m1+np_incr
                  NPNE(5,MESH3_NB(3,2),ne)=np+np_incr+
     '              MOD(m1,MESH3_S(nosphere,1))+1
                  NPNE(6,MESH3_NB(3,2),ne)=np+m1+MESH3_S(nosphere,1)+
     '              np_incr
                ENDIF
                IF(DOP) THEN
                  WRITE(OP_STRING,'('' NBJ(1..,ne):'',3(1X,I2))')
     '              (NBJ(nj,ne),nj=1,3)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' NPNE(nn,'',I2,'',ne):'','
     '              //'8(1X,I5))') NBJ(1,ne),(NPNE(nn,NBJ(1,ne),ne),
     '              nn=1,NNT(NBJ(1,ne)))
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                ENDIF
                ne=ne+1
              ENDDO !m1
              IF(m2.NE.1) THEN !Not the Bottom element
                np=np+MESH3_S(nosphere,1)
              ENDIF
            ENDDO !m2
          ENDIF
          ne=ne-1
          NET(nr)=ne
          IF(nr.EQ.1) THEN
            NEELEM(0,nr)=NET(nr)
            DO noelem=1,NET(nr)
              NEELEM(noelem,nr)=noelem
              NRE(noelem)=1
            ENDDO !noelem
          ELSE
            DO noelem=1,NET(nr)-NET(nr-1)
              NEELEM(noelem,nr)=NET(nr-1)+noelem
              NRE(NET(nr-1)+noelem)=nr
            ENDDO
            NEELEM(0,nr)=NET(nr)-NET(nr-1)
          ENDIF
          IF(DOP) THEN
            WRITE(OP_STRING,'('' NET(nr)='',I5)') NET(nr)
            WRITE(OP_STRING,'('' NEELEM(0,nr)='',I5)') NEELEM(0,nr)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NEELEM(1..,nr)='',10(1X,I5),'
     '        //'/:(10(1X,I5)))') (NEELEM(noelem,nr),noelem=1,
     '        NEELEM(0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,'('' NRE(ne)='',30(1X,I1),'
     '        //'/:(30(1X,I1)))') (NRE(NEELEM(noelem,nr)),noelem=1,
     '        NEELEM(0,nr))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
          ne=ne+1
        ENDDO !nr
        NET(0)=ne-1
        NEELEM(0,0)=ne-1 !number of elements in all regions

c cpb 12/1/96 Huh ??
C        IF(NSPHERES.EQ.1)THEN
C          NRT=1
C        ELSE
C          NRT=NSPHERES-1
C        ENDIF

      ENDIF !End of nj loop

      DO nr=1,NRT
        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
CC AJPs
          NJ_LOC(NJL_GEOM,0,nr)=NJT
CC AJPe
          DO nj=1,NJT
            NVJP(nj,np)=1
          ENDDO !nj
        ENDDO !nonode (np)
        DO noelem=1,NEELEM(0,nr)
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
CC AJPs
          NJ_LOC(NJL_GEOM,0,nr)=NJT
CC AJPe
          DO nbb=1,NBFT
            IF(NBC(nb).EQ.NBC(nbb).AND.NIT(nb).EQ.NIT(nbb)) THEN
              DO ns=1,NST(nb)+NAT(nb)
                SE(ns,nb,ne)=1.0d0
              ENDDO
              DO nn=1,NNT(nbb)
                NPNE(nn,nbb,ne)=NPNE(nn,nb,ne)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  NVJE(nn,nb,nj,ne)=1 !vers. one of nn,nj in elem ne
                ENDDO !nj
              ENDDO !nn
            ENDIF
          ENDDO !nbb
          DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
            nbb=NBJ(nj,ne)
            DO nn=1,NNT(nbb)
              DO nk=1,NKT(nn,nbb)
                NKJE(nk,nn,nj,ne)=nk
              ENDDO !nk
            ENDDO !nn
          ENDDO !nj
        ENDDO !ne
C GMH 14/2/97 create region
        IF(NPNODE(0,nr).GT.0) THEN
          CALL REGION_CREATE(nr,ERROR,*9999)
        ENDIF
      ENDDO !nr

      CALL INTERFACE(NP_INTERFACE,NPNODE,ERROR,*9999)
      CALL CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*9999)

      CALL EXITS('IPMESH3')
      RETURN
 9999 CALL ERRORS('IPMESH3',ERROR)
      CALL EXITS('IPMESH3')
      RETURN 1
      END


