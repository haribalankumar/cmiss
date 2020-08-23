      SUBROUTINE COEFF(INP,ISC_GK,ISC_GKK,ISR_GK,ISR_GKK,NBJ,NHP,NKJE,
     '  NONY,NPF,NPNE,NPNODE,nr,nr_gkk,NRE,NVHP,NVJE,NW,nx,NYNP,CONY,
     '  GK,GKK,PG,SE,XA,XE,XG,XP,ERROR,*)

C#### Subroutine: COEFF
C###  Description:
C###    COEFF calculates the integrations of the singularity in
C###    the fundamental solution for the diagonal element of GK.
C###    The diagonal entry in GK is calculated using row sums for
C###    linear elasticity and Laplace.  This assumes a finite domain.
C###    The code will need to be modified if a problem in a
C###    semi-infinite domain is being solved. For linear elasticity
C###    this assumes there are no body forces.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'b10.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ktyp90.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'tol00.cmn'
      INCLUDE 'solv00.cmn'
!     Parameter List
      INTEGER INP(NNM,NIM,NBFM),ISC_GK(NISC_GKM),ISC_GKK(NISC_GKKM),
     '  ISR_GK(NISR_GKM),ISR_GKK(NISR_GKKM),NBJ(NJM,NEM),NHP(NPM),
     '  NKJE(NKM,NNM,NJM,NEM),NONY(0:NOYM,NYM,NRCM),NPF(9,NFM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M),nr,nr_gkk,NRE(NEM),
     '  NVHP(NHM,NPM,NCM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 CONY(0:NOYM,NYM,NRCM),GK(NZ_GK_M),GKK(NZ_GKK_M),
     '  PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),
     '  XE(NSM,NJM),XG(NJM,NUM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER loop,nb,ne,ng,NG_CLOSE(3),nh,nhx,nh2,nhx2,nj,nn,
     '  no1,no2,noy1,noy2,nonode,nonode2,np,NP2,ns,nu,nv,ny1,ny2,nz,nzz
      REAL*8 A,B,C,co1,co2,CP,DIST,DIST1,SOLID_ANGLE,SUM(6),SUMXG,
     '  TN(3,2,3),XN1(3),XN2(3),XN3(3)
      LOGICAL INTERFACE

      CALL ENTERS('COEFF',*9999)

C     Write out warning if ktyp92 ge 3
      CALL ASSERT(KTYP92.LT.3,' >>Must use conventional BIE',
     '  ERROR,*9999)
      nv=1 !Temporary GMH
      DO nonode=1,NPNODE(0) !Loop over nodes surrounding each BE region.
        np=NPNODE(nonode)
        DO nhx=1,NHP(np)
          nh=NH_LOC(nhx,nx)

CC AJPs 10/10/07 - 191297 - rgb
c          IF((IGREN(nr).LE.2).OR.
c     '      (IGREN(nr).GE.7.AND.IGREN(nr).LE.8).OR.
c     '      (IGREN(nr).GE.13.AND.IGREN(nr).LE.15)) THEN
          IF((IGREN(nr).LE.2).OR.
     '      (IGREN(nr).GE.7.AND.IGREN(nr).LE.10).OR.
     '      (IGREN(nr).GE.13.AND.IGREN(nr).LE.15)) THEN
CC AJPe
            IF(CALC_GLOBAL(nr,nx)) THEN
C             We need to do rowsums adding to 0
              ny1=NYNP(1,1,nh,np,1,1,nr)
C             Zero the sum
              DO nhx2=1,NHP(np)
                nh2=NH_LOC(nhx2,nx)
                SUM(nh2)=0.0D0
              ENDDO !End of nh2 loop
C             Sum up the terms on this row - sep sums for each nh
              DO nonode2=1,NPNODE(0)
                np2=NPNODE(nonode2)
                DO nhx2=1,NHP(np2)
                  nh2=NH_LOC(nhx2,nx)
C                 Miss out the middle
                  IF(np.NE.np2) THEN
                    ny2=NYNP(1,1,nh2,np2,2,1,nr)
                    CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '                NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                    IF(nz.NE.0) SUM(nh2) = SUM(nh2)+GK(nz)
                  ENDIF
                ENDDO !End of nh2 loop
              ENDDO !End of np2 loop
C             Change GK
              DO nhx2=1,NHP(np)
                nh2=NH_LOC(nhx2,nx)
                ny2=NYNP(1,1,nh2,np,2,1,nr)
                CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '            NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                IF(nz.NE.0) GK(nz) = -SUM(nh2)
              ENDDO !End of nh2 loop
            ELSE
C             We need to do rowsums adding to 0
              ny1=NYNP(1,1,nh,np,1,1,nr)
              ny2=NYNP(1,1,nh,np,0,1,nr)
              IF(NONY(0,ny2,2).EQ.0) THEN
C               Zero the sum
                DO nhx2=1,NHP(np)
                  nh2=NH_LOC(nhx2,nx)
                  SUM(nh2)=0.0d0
                ENDDO !End of nh2 loop
C               Sum up the terms on this row - sep sums for each nh
                DO nonode2=1,NPNODE(0)
                  np2=NPNODE(nonode2)
                  DO nhx2=1,NHP(np2)
                    nh2=NH_LOC(nhx2,nx)
C                   Miss out the middle
                    IF(np.NE.np2) THEN
                      ny2=NYNP(1,1,nh2,np2,2,1,nr)
                      CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '                  NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                      IF(nz.NE.0) SUM(nh2) = SUM(nh2)+GK(nz)
                    ENDIF
                  ENDDO !End of nh2 loop
                ENDDO !End of np2 loop
C               Change GK
                DO nhx2=1,NHP(np)
                  nh2=NH_LOC(nhx2,nx)
                  ny2=NYNP(1,1,nh2,np,2,1,nr)
                  CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '              NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
                  IF(nz.NE.0) GK(nz) = -SUM(nh2)
                ENDDO !End of nh2 loop
              ELSE
                DO noy1=1,NONY(0,ny1,1)
                  no1=NONY(noy1,ny1,1)
                  co1=CONY(noy1,ny1,1)
                  IF(DABS(co1).GT.ZERO_TOL) THEN
C                   Zero the sum
                    DO nhx2=1,NHP(np)
                      nh2=NH_LOC(nhx2,nx)
                      SUM(nh2)=0.0d0
                    ENDDO !End of nh2 loop
C                   Sum up the terms on this row - sep sums for each nh
                    DO nonode2=1,NPNODE(0)
                      np2=NPNODE(nonode2)
                      DO nhx2=1,NHP(np2)
                        nh2=NH_LOC(nhx2,nx)
C                       Miss out the middle
                        IF(np.NE.np2) THEN
                          ny2=NYNP(1,1,nh2,np2,0,1,nr)
                          DO noy2=1,NONY(0,ny2,2)
                            no2=NONY(noy2,ny2,2)
                            co2=CONY(noy2,ny2,2)
                            IF(DABS(co1*co2).GT.ZERO_TOL) THEN
                              CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),
     '                          nzz,NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,
     '                          ISR_GKK,SPARSEGKK(nx),ERROR,*9999)
                              SUM(nh2) = SUM(nh2)+GKK(nzz)/(co1*co2)
                            ENDIF
                          ENDDO !noy2
                          IF(NONY(0,ny2,2).EQ.0) THEN
                            ny2=NYNP(1,1,nh2,np2,2,1,nr)
                            CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '                        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,
     '                        ERROR,*9999)
                            IF(nz.NE.0) SUM(nh2) = SUM(nh2)+GK(nz)
                          ENDIF
                        ENDIF
                      ENDDO !End of nh2 loop
                    ENDDO !End of np2 loop
C                   Change GK
                    DO nhx2=1,NHP(np)
                      nh2=NH_LOC(nhx2,nx)
                      ny2=NYNP(1,1,nh2,np,0,1,nr)
                      DO noy2=1,NONY(0,ny2,2)
                        no2=NONY(noy2,ny2,2)
                        co2=CONY(noy2,ny2,2)
                        CALL SPARSE(no1,no2,NOT(1,1,nr_gkk,nx),nzz,
     '                    NZ_GKK_M,NZZT(1,nr_gkk,nx),ISC_GKK,ISR_GKK,
     '                    SPARSEGKK(nx),ERROR,*9999)
C!!! This is assuming that multiple mesh degrees of freedom in different
C!!! regions do not map to the same diagonal solution degree-of-freedom
                        GKK(nzz) = -SUM(nh2)
                      ENDDO !noy2
                    ENDDO !End of nh2 loop
                  ENDIF
                ENDDO !noy1
              ENDIF
            ENDIF
          ELSE
C           Find ny1 and ny2 for the following code.  This is only
C           based on nh for row - ie 1st dependent var along row.
            ny1=NYNP(1,1,nh,np,1,1,nr)
            ny2=NYNP(1,1,nh,np,2,1,nr)
            CALL SPARSE(ny1,ny2,NYT(1,1,nx),nz,NZ_GK_M,
     '        NZT(1,nx),ISC_GK,ISR_GK,KTYP24,ERROR,*9999)
            IF(NVHP(1,np,2).LE.1) THEN !no corner
              IF(COMPLEX) THEN
              ELSE
                GK(nz)=GK(nz)+0.5D0
              ENDIF
            ELSE IF(NVHP(1,np,2).EQ.2) THEN !Only 2d corner or 3d edge

!warning that this needs checking (and in next loop)
C ***       Find angle between outward normals or elements at the
C ***       corner.
C ***       Use this to calculate the corner contribution.
C ***       Calculate normal at first Gauss point of the high order
C ***       scheme.
              IF(NIT(1).EQ.1) THEN !1d integral
                DO loop=1,2
                  IF(loop.EQ.1) THEN
! NEEDS FIXING.
!Search to find the ne's containing node np once we know we have more than
!one version at the current np.
!                  ne=? ! One of the corner node elements
                  ELSE
!                    ne=?
                  ENDIF
                  ne=1 !temporary cpb 9/7/95
                  nb=NBJ(1,ne)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
                  IF(np.EQ.NPNE(1,nb,ne)) THEN !node is 1st node of corn
                    ng=1 !use first Gauss point
                  ELSE ! node must be last node of corner element
                    ng=NGT(nb) !use last Gauss point
                  ENDIF
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    DO nu=1,NUT(nb)
                      SUMXG=0.0D0
                      DO ns=1,NST(nb)
                        SUMXG=SUMXG+PG(ns,nu,ng,nb)*XE(ns,nj)
                      ENDDO
                      XG(nj,nu)=SUMXG
                    ENDDO
                  ENDDO
                  INTERFACE=.FALSE. !May need changing AJP 21-1-93
                  !Find unit out. normal
                  CALL NORMAL(ne,nr,NW,XG,XN2,INTERFACE,ERROR,*9999)
                  IF(loop.EQ.1) THEN
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      XN1(nj)=XN2(nj)
                    ENDDO
                  ENDIF
C ***             Find a tangent vector and adjust it so that it lies in
C ***             the element
                  nn=1
                  DO WHILE(np.NE.NPNE(nn,nb,ne))
                    nn=nn+1 !Find local node number of np
                  ENDDO
                  IF(INP(nn,1,nb).GT.1) THEN !Want negative of the tangent
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,1,loop)=-XG(nj,2)
                    ENDDO
                  ELSE
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,1,loop)=XG(nj,2)
                    ENDDO
                  ENDIF
                ENDDO !End of loop
!                ne=? need to search (or have an array nenv?)
C ***           Find angle between xn1,xn2
                CALL ANGLE(TN(1,1,1),TN(1,1,2),A,XN1,XN2,
     '                     ERROR,*9999)
                CP=A/(2.0D0*PI)
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(COEFF_1)
                  WRITE(OP_STRING,'('' Corner coefficient for node '','
     '              //'I5,''='',D12.4)') np,CP
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(COEFF_1)
                ENDIF
C ***         2d integrals
              ELSE IF(NIT(1).EQ.2) THEN !2d integrals (i.e. 3d edge)
                DO loop=1,2
                  IF(loop.EQ.1) THEN
!                   ne=? need to search (or have an array nenv?)
                  ELSE
!                   ne=? need to search (or have an array nenv?)
                  ENDIF
                  nb=NBJ(1,ne)
                  CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '              NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '              SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C ***             Find Gauss point closest to the current node
                  DIST1=RMAX
                  DO ng=1,NGT(nb)
                    DIST=0.0D0
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      SUMXG=0.0D0
                      DO ns=1,NST(nb)
                        SUMXG=SUMXG+PG(ns,1,ng,nb)*XE(ns,nj)
                      ENDDO
                      DIST=DIST+(SUMXG-XP(1,nv,nj,np))**2
                    ENDDO
                    IF(DIST.LE.DIST1) THEN
                      DIST1=DIST
                      NG_CLOSE(loop)=ng
                    ENDIF
                  ENDDO
C ***             Find normal derivative at this point
                  ng=NG_CLOSE(loop)
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    DO nu=1,NUT(nb)
                      SUMXG=0.0D0
                      DO ns=1,NST(nb)
                        SUMXG=SUMXG+PG(ns,nu,ng,nb)*XE(ns,nj)
                      ENDDO
                      XG(nj,nu)=SUMXG
                    ENDDO
                  ENDDO
                  INTERFACE=.FALSE. !May need changing AJP 21-1-93
                  CALL NORMAL(ne,nr,NW,XG,XN2,INTERFACE,ERROR,*9999)
                  !Find unit out. normal.
                  IF(loop.EQ.1) THEN
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      XN1(nj)=XN2(nj)
                    ENDDO
                  ENDIF
C ***             Find two tangent vectors in the plane and adjust them
C ***             so that they lie in the element
                  nn=1
                  DO WHILE(np.NE.NPNE(nn,nb,ne))
                    nn=nn+1 !Find local node number of np
                  ENDDO
                  IF(INP(nn,1,nb).GT.1) THEN !Want negative of the tangent
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,1,loop)=-XG(nj,2)
                    ENDDO
                  ELSE
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,1,loop)=XG(nj,2)
                    ENDDO
                  ENDIF
                  IF(INP(nn,2,nb).GT.1) THEN !Want negative of the tangent
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,2,loop)=-XG(nj,4)
                    ENDDO
                  ELSE
                    DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                      TN(nj,2,loop)=XG(nj,4)
                    ENDDO
                  ENDIF
                ENDDO !End of loop
!Search to find the ne's containing node np once we know we have more than
!one version at the current np.
!                    ne=? ! One of the corner node elements
C ***           Find angle between xn1,xn2
                CALL ANGLE(TN(1,1,1),TN(1,1,2),A,XN1,XN2,
     '                     ERROR,*9999)
                CP=(2.0D0*A)/(4.0D0*PI) !Solid angle for an edge
                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(COEFF_2)
                  WRITE(OP_STRING,'('' Edge coefficient for node '','
     '              //'I5,''='',D12.4)') np,CP
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(COEFF_2)
                ENDIF
              ELSE
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(COEFF_3)
                WRITE(OP_STRING,'('' >>At an edge and neither 1 or '
     '            //'2d integral ??'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(COEFF_3)
                GOTO 9999
              ENDIF
              IF(COMPLEX) THEN
C                GKC(nz)=GKC(nz)+DCMPLX(CP,0.0D0)
              ELSE
                GK(nz)=GK(nz)+CP
              ENDIF
            ELSE !3d corner
C ***         Find solid angle at the corner
              DO loop=1,3
                IF(loop.EQ.1) THEN
!                    ne=? ! One of the corner node elements
                ELSE IF(loop.EQ.2) THEN
!                    ne=? ! Second of the corner node elements
                ELSE
!                    ne=? ! third corner node element
                ENDIF
                nb=NBJ(1,ne)
                CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '            NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '            SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
C ***           Find Gauss point closest to the current node
                DIST1=RMAX
                DO ng=1,NGT(nb)
                  DIST=0.0D0
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    SUMXG=0.0D0
                    DO ns=1,NST(nb)
                      SUMXG=SUMXG+PG(ns,1,ng,nb)*XE(ns,nj)
                    ENDDO
                    DIST=DIST+(SUMXG-XP(1,nv,nj,np))**2
                  ENDDO
                  IF(DIST.LE.DIST1) THEN
                    DIST1=DIST
                    NG_CLOSE(loop)=ng
                  ENDIF
                ENDDO
C ***           Find normal derivative at this point
                ng=NG_CLOSE(loop)
                DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                  DO nu=1,NUT(nb)
                    SUMXG=0.0D0
                    DO ns=1,NST(nb)
                      SUMXG=SUMXG+PG(ns,nu,ng,nb)*XE(ns,nj)
                    ENDDO
                    XG(nj,nu)=SUMXG
                  ENDDO
                ENDDO
                INTERFACE=.FALSE. !May need changing AJP 21-1-93
                CALL NORMAL(ne,nr,NW,XG,XN3,INTERFACE,ERROR,*9999)
                !Find unit outward normal.
                IF(loop.EQ.1) THEN
                  DO nj=1,3
                    XN1(nj)=XN3(nj)
                  ENDDO
                ELSE IF(loop.EQ.2) THEN
                  DO nj=1,3
                    XN2(nj)=XN3(nj)
                  ENDDO
                ENDIF
C ***           Find two tangent vectors in the plane and adjust them
C ***           so that they lie in the element
                nn=1
                DO WHILE(np.NE.NPNE(nn,nb,ne))
                  nn=nn+1 !Find local node number of np
                ENDDO
                IF(INP(nn,1,nb).GT.1) THEN !Want negative of the tangent
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,1,loop)=-XG(nj,2)
                  ENDDO
                ELSE
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,1,loop)=XG(nj,2)
                  ENDDO
                ENDIF
                IF(INP(nn,2,nb).GT.1) THEN !Want negative of the tangent
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,2,loop)=-XG(nj,4)
                  ENDDO
                ELSE
                  DO nj=1,NJ_LOC(NJL_GEOM,0,nr)
                    TN(nj,2,loop)=XG(nj,4)
                  ENDDO
                ENDIF
              ENDDO !End of loop
!                ne=? ! One of the corner node elements
C ***         Find angles between xn1,xn2  xn1,xn3  and  xn2,xn3
              CALL ANGLE(TN(1,1,1),TN(1,1,2),A,XN1,XN2,
     '                   ERROR,*9999)
              CALL ANGLE(TN(1,1,1),TN(1,1,3),B,XN1,XN3,
     '                   ERROR,*9999)
              CALL ANGLE(TN(1,1,2),TN(1,1,3),C,XN2,XN3,
     '                   ERROR,*9999)
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(COEFF_4)
              WRITE(OP_STRING,'('' A='',D12.4,'', B='',D12.4,'','
     '          //' C='',D12.4)') A,B,C
              CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(COEFF_4)
              SOLID_ANGLE=(A+B+C-PI)/(4.0D0*PI)
              IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(COEFF_5)
                WRITE(OP_STRING,'('' Solid angle at '',I5,'' = '','
     '            //'D12.5)') np,SOLID_ANGLE
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(COEFF_5)
              ENDIF
              IF(COMPLEX) THEN
C                GKC(nz)=GKC(nz)+DCMPLX(SOLID_ANGLE,0.0D0)
              ELSE
                GK(nz)=GK(nz)+SOLID_ANGLE
              ENDIF
            ENDIF !End of 3d corner
          ENDIF !Not laplace or linear elasticity
        ENDDO !End of nh loop
      ENDDO !End of np loop
      CALL EXITS('COEFF')
      RETURN
9999  CALL ERRORS('COEFF',ERROR)
      CALL EXITS('COEFF')
      RETURN 1
      END



