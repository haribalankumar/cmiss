      SUBROUTINE UPGRID_CONN(NBJ,NEELEM,NKJE,NPF,NPNE,NQXI,
     '  nr1,nr2,NVJE,NXQ,NQS,NQNE,
     '  SE,XA,XE,XP,XQ,ERROR,*)

C#### Subroutine: UPGRID_CONN
C###  Description:
C###    UPGRID_CONN updates grid point connectivity in region nr1 by
C       cleavage planes defined as elements in region nr2.
C       For Bidomain, a cleavage plane only changes connectivities in the
C       intracellular domain, so a break is represented by placing a
C       negative sign in front of grid numbers (in NXQ) which are not
C       connected in this domain.
C       Created by Darren Hooks 2001

C       Note that this routine only works for 3D geometries (elements
C       containing bilinear cleavage plane elements)
C       It is required that the element bounding the cleavage planes
C       is of orthogonal (box) geometry, and has xi directions aligning
C       with the xyz axes (ie. xi1 is in x direction). Orientation of the
C       cleavage planes within the bounding element is arbitrary.

      IMPLICIT NONE
      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='UPGRID_CONN')
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NQXI(0:NIM,NQSCM),nr1,nr2,
     '  NVJE(NNM,NBFM,NJM,NEM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  NQS(NEQM),NQNE(NEQM,NQEM)
      REAL*8 SE(NSM,NBFM,NEM),
     '  XA(NAM,NJM,NEM),XE(NSM,NJM),XP(NKM,NVM,NJM,NPM),XQ(NJM,NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ii,STEP1(3),STEP2(3),
     '  DIR1(3),DIR2(3),noelem_nr1,noelem_nr2,ne1,ne2,
     '  ij,gridscheme,nq_dirn1,nq_dirn2,nq_dirn3,count1,
     '  count2,count3,prev_nq
      REAL*8 DIST,PREV_DIST,
     '  a,b,c,d,aa,bb,cc,dd,alpha,
     '  beta,betadelta,gamma,SBG,SBA,BBAG,XI1,XI2,bcXI1,
     '  bbccXI1,daXI1,ddaaXI1
      LOGICAL haveXI1,haveXI2
!     Functions

      CALL ENTERS(ROUTINENAME,*9999)

C set up the directions of travel (this requires bounding element is
C aligned with xyz axes.
        STEP1(1)=3
        STEP1(2)=3
        STEP1(3)=2
        STEP2(1)=2
        STEP2(2)=1
        STEP2(3)=1

        DIR1(1)=2
        DIR1(2)=3
        DIR1(3)=1
        DIR2(1)=3
        DIR2(2)=1
        DIR2(3)=2

        DO noelem_nr1=1,NEELEM(0,nr1)!loop over elements in reg. 1
          ne1=NEELEM(noelem_nr1,nr1)
          gridscheme=NQS(ne1)

          DO noelem_nr2=1,NEELEM(0,nr2)!loop over elems in region nr2:
C            these are the cleavage plane elements
            ne2=NEELEM(noelem_nr2,nr2)

            !get info about each cleavage plane in region nr2
            CALL XPXE(NBJ(1,ne2),NKJE(1,1,1,ne2),NPF(1,1),NPNE(1,1,ne2),
     '        nr2,NVJE(1,1,1,ne2),SE(1,1,ne2),XA(1,1,ne2),XE,XP,ERROR,
     '        *9999)

            DO ij=1,NQXI(0,gridscheme)!loop over xi dirns of grd scheme
              nq_dirn1=NQNE(ne1,1)!get the global grid point number of
C             the first grid point in element ne1
              count1=1
              DO WHILE (count1.LE.NQXI(STEP1(ij),gridscheme))!count1
                count2=1
                nq_dirn2=nq_dirn1
                DO WHILE (count2.LE.NQXI(STEP2(ij),gridscheme))!count2
C                  calculate the parameters a,b,c,d,A,B,C,D,alpha,beta,
C                 gamma for the grid point and cleavage plane
                  a=XE(2,DIR1(ij))-XE(1,DIR1(ij))
                  b=XE(3,DIR1(ij))-XE(1,DIR1(ij))
                  c=XE(1,DIR1(ij))-XE(2,DIR1(ij))-XE(3,DIR1(ij))+
     '              XE(4,DIR1(ij))
                  d=XQ(DIR1(ij),nq_dirn2)-XE(1,DIR1(ij))
                  aa=XE(2,DIR2(ij))-XE(1,DIR2(ij))
                  bb=XE(3,DIR2(ij))-XE(1,DIR2(ij))
                  cc=XE(1,DIR2(ij))-XE(2,DIR2(ij))-XE(3,DIR2(ij))+
     '              XE(4,DIR2(ij))
                  dd=XQ(DIR2(ij),nq_dirn2)-XE(1,DIR2(ij))
                  alpha=aa*c-a*cc
                  beta=d*cc-dd*c+aa*b-a*bb
                  gamma=bb*d-b*dd

                  BBAG=beta*beta-4*alpha*gamma

                  IF (BBAG.GT.0.d0) THEN !solution
                    !find value of Xi1 if it exists (0,1)
                    haveXI1=.FALSE.
                    betadelta=DABS(beta)+SQRT(BBAG)
                    SBG=-2.0d0*SIGN(1.d0,beta)*gamma
                    IF(betadelta.GE.SBG.AND.SBG.GE.0.0d0)THEN
                      haveXI1=.TRUE.
                      XI1=SBG/betadelta
                    ELSE
                      SBA=-2.0d0*SIGN(1.d0,beta)*alpha
                      IF(SBA.GE.betadelta.AND.SBA.NE.0.0d0)THEN
                        haveXI1=.TRUE.
                        XI1=betadelta/SBA
                      ENDIF
                    ENDIF
                    IF (haveXI1) THEN !in plane in xi1
                      !find value of Xi2 if it exists (0,1)
                      haveXI2=.FALSE.
                      bcXI1=b+c*XI1
                      bbccXI1=bb+cc*XI1
                      IF(DABS(bbccXI1).GT.DABS(bcXI1))THEN
                        ddaaXI1=dd-aa*XI1
                        IF(DABS(bbccXI1).GE.DABS(ddaaXI1)
     '                    .AND.bbccXI1.NE.0.0d0)THEN
                          XI2=ddaaXI1/bbccXI1
                          IF(XI2.GE.0.0d0) haveXI2=.TRUE.
                        ENDIF
                      ELSE
                        daXI1=d-a*XI1
                        IF(DABS(bcXI1).GE.DABS(daXI1)
     '                    .AND.bcXI1.NE.0.0d0)THEN
                          XI2=daXI1/bcXI1
                          IF(XI2.GE.0.0d0) haveXI2=.TRUE.
                        ENDIF
                      ENDIF
                      IF (haveXI2) THEN!inside plane in both xi1, xi2
                        PREV_DIST=0.d0
                        nq_dirn3=nq_dirn2
                        count3=1
                        prev_nq=0
                        DO WHILE (count3.LE.NQXI(ij,gridscheme))
C                         calculate the distance in ij direction
C                         from grid pt nq_dirn3 to its projected pt.
                          DIST=XQ(ij,nq_dirn3)-((1-XI1)*(1-XI2)*
     '                      XE(1,ij)+XI1*(1-XI2)*XE(2,ij)+XI2*(1-XI1)*
     '                      XE(3,ij)+XI1*XI2*XE(4,ij))
                          IF (PREV_DIST*DIST.LT.0.d0) THEN !we have
C                           traversed the plane

C                       For Bidomain, a cleavage plane only changes
C                       connectivities in the intracellular domain, so
C                       a break is represented by placing a negative
C                       sign in front of grid numbers which are not
C                       connected in this domain.
                            DO ii=1,NXQ(-ij,0,nq_dirn3,1)
                              NXQ(-ij,ii,nq_dirn3,1)=
     '                          SIGN(NXQ(-ij,ii,nq_dirn3,1),-1)
                            ENDDO
                            DO ii=1,NXQ(ij,0,prev_nq,1)
                              NXQ(ij,ii,prev_nq,1)=
     '                          SIGN(NXQ(ij,ii,prev_nq,1),-1)
                            ENDDO
                            count3=NQXI(ij,gridscheme)+1!exit loop
                          ELSE
                            prev_nq=nq_dirn3
C                           step to next grid pt. in ij direction
                            nq_dirn3=IABS(NXQ(ij,1,nq_dirn3,1))
                            PREV_DIST=DIST
                            count3=count3+1
                          ENDIF
                        ENDDO !count3
                      ENDIF !inside plane in xi2
                    ENDIF !inside plane in xi1
                  ENDIF !real xi1 solution
                  count2=count2+1
                  nq_dirn2=IABS(NXQ(STEP2(ij),1,nq_dirn2,1))
                ENDDO!count2
                nq_dirn1=IABS(NXQ(STEP1(ij),1,nq_dirn1,1))
                count1=count1+1
              ENDDO!count1
            ENDDO!ij
          ENDDO !loop over elements in region nr=2
        ENDDO !loop over elements in region nr=1
C      ENDIF !more than one element in region nr=1

C--type---------------------------------------------------------------
C*** PJH June 1999


      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END


