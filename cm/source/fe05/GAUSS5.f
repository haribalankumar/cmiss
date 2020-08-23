      SUBROUTINE GAUSS5(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,
     '  DET,PG,WG,XIG,ERROR,*)

C#### Subroutine: GAUSS5
C###  Description:
C###    GAUSS5 defines parameters for a BE family of basis functions.
C###    For each family member a call is made to Gauss10 to set up the
C###    XIG, WG and PG arrays.

C**** The basis number, nb, passed to this routine is interpreted as
C**** the parent basis number of the boundary element family.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),nb,NBTOP,NDET(NBFM,0:NNM),NGAP(NIM,NBM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,ISP,nbbem,nb_incr,ng,NGHIGH,
     '  NGHIGH2,NGLOW,NGLOW2,ni,nk,nn,ns,nu
      REAL*8 ALIM(2,99),BLIM(2,99),PSI1,TEST(100,5)

      CALL ENTERS('GAUSS5',*9999)

c cpb 11/7/95 Initialise NDET and DET
      DO nn=0,NNM
        NDET(nb,nn)=1
        DO ng=1,NGM
          DO i=1,6
            DET(nb,nn,ng,i)=1.0d0
          ENDDO !i
        ENDDO !ng
      ENDDO !nn
      NGLOW=NGLIMITS(1,nb,1)
      NGHIGH=NGLIMITS(1,nb,2)
      IF(NIT(nb).EQ.1) THEN !1D integrals
        IF(IBT(1,1,nb).EQ.1) THEN !Lagrange basis function
          IF(IBT(2,1,nb).EQ.1) THEN !linear basis function
C cpb 4/11/96 Adding logarithmic Gaussian quadrature. For the moment
C just create another basis function after the element splitting
C schemes. This is wasteful for non Laplace type Greens functions
C but then all of these basis functions should be calculated as they
C are needed.
C            NBASEF(nb,0)=6
            NBASEF(nb,0)=8
            NFBASE(1,nb)=nb  !family number of global basis number nb
            NFBASE(2,nb)=1   !parent is first child
            !Require a low, mid and 4 high order schemes.
            !First extra high order scheme is for np at nn=1
            !Second extra high order scheme is for np at nn=2
            !Last extra high order scheme is for adaptive int.
            NGAP(1,NBTOP+1)=NGHIGH
            NGAP(1,NBTOP+2)=NGHIGH
C           + another two high order log schemes (1 for each split)
            NGAP(1,NBTOP+3)=NGHIGH !log scheme
            NGAP(1,NBTOP+4)=NGHIGH !log scheme
C            NGAP(1,NBTOP+3)=(NGHIGH+NGLOW)/2
C            NGAP(1,NBTOP+4)=NGLOW
C            NGAP(1,NBTOP+5)=NGHIGH !Adaptive scheme
            NGAP(1,NBTOP+5)=(NGHIGH+NGLOW)/2
            NGAP(1,NBTOP+6)=NGLOW
            NGAP(1,NBTOP+7)=NGHIGH !Adaptive scheme
            ALIM(1,nb)=0.0d0
            BLIM(1,nb)=1.0d0
            DO nbbem=NBTOP+1,NBTOP+NBASEF(nb,0)-1
              ALIM(1,nbbem)=0.0d0
              BLIM(1,nbbem)=1.0d0
            ENDDO
          ELSE IF(IBT(2,1,nb).EQ.2) THEN !quadratic basis function
            NBASEF(nb,0)=8
            NFBASE(1,nb)=nb  !family number of global basis number nb
            NFBASE(2,nb)=1   !parent is first child
            !First extra high scheme is for np at nn=1
            !Next 2 extra high schemes are for np at nn=2
            !Second to last extra high scheme is for np at nn=3
            !Last extra scheme is for adaptive integration.
            NGAP(1,NBTOP+1)=NGHIGH
            NGAP(1,NBTOP+2)=NGHIGH !Integral from 0 to 1/2
            NGAP(1,NBTOP+3)=NGHIGH !Integral from 1/2 to 1
            NGAP(1,NBTOP+4)=NGHIGH
            NGAP(1,NBTOP+5)=(NGHIGH+NGLOW)/2
            NGAP(1,NBTOP+6)=NGLOW
            NGAP(1,NBTOP+7)=NGHIGH
            ALIM(1,nb)=0.0d0
            BLIM(1,nb)=1.0d0
            ALIM(1,NBTOP+1)=0.0d0 !np at nn=1 scheme
            BLIM(1,NBTOP+1)=1.0d0
            ALIM(1,NBTOP+2)=0.0d0 !np at nn=2 scheme
            BLIM(1,NBTOP+2)=0.5d0
            ALIM(1,NBTOP+3)=0.5d0 !np at nn=2 scheme
            BLIM(1,NBTOP+3)=1.0d0
            DO nbbem=NBTOP+4,NBTOP+NBASEF(nb,0)-1
              ALIM(1,nbbem)=0.0d0
              BLIM(1,nbbem)=1.0d0
            ENDDO
            NDET(nb,2)=2
            !For midside node split element into two parts.
            !Determinant automatically handled in basis
            !function setup call to Gauss10.
          ELSE IF(IBT(2,1,nb).EQ.3) THEN !cubic Lagrange basis function
            ERROR='>>Not implemented'
            GOTO 9999
          ENDIF
        ELSE IF(IBT(1,1,nb).EQ.2) THEN !Hermite basis function
          IF(IBT(2,1,1).EQ.1) THEN !Cubic Hermite basis function
C cpb 4/11/96 Adding logarithmic Gaussian quadrature. For the moment
C just create another basis function after the element splitting
C schemes. This is wasteful for non Laplace type Greens functions
C but then all of these basis functions should be calculated as they
C are needed.
            !Generate the same number of basis functions as in the
            !linear Lagrange case (since linear interpolation is
            !to be used for the normal derivative)
C            NBASEF(nb,0)=6
            NBASEF(nb,0)=8
            NFBASE(1,nb)=nb  !family number of global basis number nb
            NFBASE(2,nb)=1   !parent is first child
            !Require a low, mid and 4 high order schemes.
            !First extra high order scheme is for np at nn=1
            !Second extra high order scheme is for np at nn=2
            !Last extra high order scheme is for adaptive int.

C LKC 19-MAR-1998 Assert on NGAP when accessing NBTOP+7 below
            CALL ASSERT(NBM.GE.NBTOP+7,'>>Increase NBM',ERROR,*9999)

            NGAP(1,NBTOP+1)=NGHIGH
            NGAP(1,NBTOP+2)=NGHIGH
C           + another two high order log schemes (1 for each split)
            NGAP(1,NBTOP+3)=NGHIGH !log scheme
            NGAP(1,NBTOP+4)=NGHIGH !log scheme
C            NGAP(1,NBTOP+3)=(NGHIGH+NGLOW)/2
C            NGAP(1,NBTOP+4)=NGLOW
C            NGAP(1,NBTOP+5)=NGHIGH !Adaptive scheme
            NGAP(1,NBTOP+5)=(NGHIGH+NGLOW)/2
            NGAP(1,NBTOP+6)=NGLOW
            NGAP(1,NBTOP+7)=NGHIGH !Adaptive scheme

            ALIM(1,nb)=0.0d0
            BLIM(1,nb)=1.0d0
            DO nbbem=NBTOP+1,NBTOP+NBASEF(nb,0)-1
              ALIM(1,nbbem)=0.0d0
              BLIM(1,nbbem)=1.0d0
            ENDDO
          ENDIF
        ENDIF
      ELSE IF(NIT(nb).EQ.2) THEN !2D integrals
        IF(IBT(2,1,nb).EQ.1.AND.IBT(2,2,nb).EQ.1) THEN
          !Bilinear, bicubic hermite or some combination of
          !linear and cubic hermite
C LKC 3-AUG-97 Added assert line
          CALL ASSERT(NBM.GE.NBTOP+11,'>>Increase NBM to NBTOP+11',
     '      ERROR,*9999)
          NGLOW2=NGLIMITS(2,nb,1)
          !Low order scheme in second Xi direction
          NGHIGH2=NGLIMITS(2,nb,2)
          !High order scheme in second Xi direction
          NBASEF(nb,0)=12
          NFBASE(1,nb)=nb  !family number of global basis number nb
          NFBASE(2,nb)=1   !parent is first child
          !Require a low, mid, 9 high order schemes and a general
          !adaptive scheme.  There is 1 general high order scheme
          !and two high order schemes for each location of the
          !singularity in the element (square element is split
          !into 2 triangular elements based on the singularity
          !location)
          DO nb_incr=1,8
            NGAP(1,NBTOP+nb_incr)=NGHIGH
            NGAP(2,NBTOP+nb_incr)=NGHIGH2
          ENDDO
          NGAP(1,NBTOP+9)=(NGHIGH+NGLOW)/2
          NGAP(2,NBTOP+9)=(NGHIGH2+NGLOW2)/2
          NGAP(1,NBTOP+10)=NGLOW
          NGAP(2,NBTOP+10)=NGLOW2
          NGAP(1,NBTOP+11)=NGHIGH !adaptive scheme
          NGAP(2,NBTOP+11)=NGHIGH2
          DO ni=1,2
            ALIM(ni,nb)=0.0d0
            BLIM(ni,nb)=1.0d0
            DO nb_incr=1,11
              ALIM(ni,NBTOP+nb_incr)=0.0d0
              BLIM(ni,NBTOP+nb_incr)=1.0d0
            ENDDO
          ENDDO
        ELSE
          ERROR='>>This 2d BE basis function is not implemented'
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          GOTO 9999
        ENDIF
      ENDIF !End of 1D or 2D loop
C Need to check total number of basis functions here
      DO nbbem=1,NBASEF(nb,0)
        NBASEF(nb,nbbem)=NBTOP+nbbem-1
      ENDDO
      NBASEF(nb,1)=nb
      CALL GAUSS10(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),
     '  nb,NGAP(1,nb),ALIM(1,nb),BLIM(1,nb),PG(1,1,1,nb),WG(1,nb),
     '  XIG(1,1,nb),ERROR,*9999)
      NGT(nb)=1
      DO ni=1,NIT(nb)
        NGT(nb)=NGT(nb)*NGAP(ni,nb)
      ENDDO
      CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
     '   ERROR,*9999)
      IF(NIT(nb).EQ.2) THEN !2d Integrals
        !Need to calculate the preimage of the Gauss points in
        !the (psi1',psi2') coordinate system.  The mapping from
        !a split element in the (psi1,psi2) space to the
        !(psi1',psi2') space depends on the location of the
        !node at which the singularity is placed.  The
        !determinant of this map for the ith part of the
        !split element is stored in DET(nb,nn,ng,i).
        !Once the preimage of the Gauss points has been found then
        !evaluate the basis function at these preimage points.
        IF(IBT(2,1,nb).EQ.1.AND.IBT(2,2,nb).EQ.1) THEN
          !Bilinear, bicubic hermite or some combination
          !linear and cubic hermite.
          DO nn=1,NNT(nb)
            NDET(nb,nn)=2
            DO ng=1,NGT(nb)
              !Find Jacobian of the 2 transformations
              !NOTE Assumes standard global to local coords used.
              DET(nb,nn,ng,1)=XIG(1,ng,nb) !psi1'
              DET(nb,nn,ng,2)=XIG(2,ng,nb) !psi2'
            ENDDO
          ENDDO
          !Find preimage of Gauss points.
          DO ng=1,NGT(nb)
            !Singularity at nn=1, first part of split element
            XIG(1,ng,NBTOP+1)=XIG(1,ng,nb) !psi1'
            XIG(2,ng,NBTOP+1)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
            !Singularity at nn=1, second part of split element
            XIG(1,ng,NBTOP+2)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
            XIG(2,ng,NBTOP+2)=XIG(2,ng,nb) !psi2'
            !Singularity at nn=2, first part of split element
            XIG(1,ng,NBTOP+3)=1.0D0-XIG(1,ng,nb)
            XIG(2,ng,NBTOP+3)=XIG(1,ng,nb)*(1.0D0-XIG(2,ng,nb))
            !Singularity at nn=2, second part of split element
            XIG(1,ng,NBTOP+4)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            XIG(2,ng,NBTOP+4)=XIG(2,ng,nb)
            !Singularity at nn=3, first part of split element
            XIG(1,ng,NBTOP+5)=XIG(1,ng,nb)
            XIG(2,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            !Singularity at nn=3, second part of split element
            XIG(1,ng,NBTOP+6)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
            XIG(2,ng,NBTOP+6)=1.0D0-XIG(2,ng,nb)
            !Singularity at nn=4, first part of split element
            XIG(1,ng,NBTOP+7)=1.0D0-XIG(1,ng,nb)
            XIG(2,ng,NBTOP+7)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            !Singularity at nn=4, second part of split element
            XIG(1,ng,NBTOP+8)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            XIG(2,ng,NBTOP+8)=1.0D0-XIG(2,ng,nb)
          ENDDO
          !Evaluate basis functions at each of these points.
          DO nb_incr=1,8
            DO ng=1,NGT(nb)
              ns=0
              DO nn=1,NNT(nb)
                DO nk=1,NKT(nn,nb)
                  ns=ns+1
                  DO nu=1,NUT(nb)
                    PG(ns,nu,ng,NBTOP+nb_incr) =
     '                PSI1(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,nu,
     '                nk,nn,XIG(1,ng,NBTOP+nb_incr))
                  ENDDO !End nu loop
                ENDDO !End nn loop
              ENDDO !End nk loop
              NGT(NBTOP+nb_incr)=NGT(nb)
              WG(ng,NBTOP+nb_incr)=WG(ng,nb)
            ENDDO !End of ng loop
          ENDDO !End NB_INCR loop
        ENDIF !End of Bilinear loop
      ENDIF !End of 2d integrals
      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        TEST(nb,1)=0.0d0
        TEST(nb,2)=0.0d0
        DO nb_incr=0,8
          TEST(NBTOP+nb_incr,1)=0.0d0
          TEST(NBTOP+nb_incr,2)=0.0d0
        ENDDO
        DO ng=1,NGT(nb)
          TEST(nb,1)=TEST(nb,1)+WG(ng,nb)*XIG(1,ng,nb)
          TEST(nb,2)=TEST(nb,2)+WG(ng,nb)*XIG(2,ng,nb)
        ENDDO
        WRITE(OP_STRING,*)' Basis number',nb
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' Integral of xi1 over element=',
     '     TEST(nb,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,*)' Integral of xi2 over element=',
     '     TEST(nb,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        IF(NIT(nb).EQ.2) THEN !2d Integrals
          DO ng=1,NGT(nb)
            DO nb_incr=1,4
              TEST(NBTOP+nb_incr,1)=TEST(NBTOP+nb_incr,1)+
     '            WG(ng,NBTOP+nb_incr)*DET(nb,nb_incr,ng,1)
              TEST(NBTOP+nb_incr,2)=TEST(NBTOP+nb_incr,2)+
     '            WG(ng,NBTOP+nb_incr)*DET(nb,nb_incr,ng,2)
            ENDDO
            ISP=0
            DO I=1,NDET(nb,1)
              !Integral of 1/(r*r) for each part of split element
              !with singularity at node 1
              TEST(NBTOP+5,I)=TEST(NBTOP+5,I)+
     '           WG(ng,NBTOP+I)*DET(nb,1,ng,I)/
     '           DSQRT(XIG(1,ng,NBTOP+I)**2+
     '           XIG(2,ng,NBTOP+I)**2)
            ENDDO
            ISP=ISP+NDET(nb,1)
            DO I=1,NDET(nb,2)
              !Integral of 1/(r*r) for each part of split element
              !with singularity at node 2
              TEST(NBTOP+6,I)=TEST(NBTOP+6,I)+
     '            WG(ng,NBTOP+I+ISP)*DET(nb,2,ng,I)/
     '            DSQRT((1.0D0-XIG(1,ng,NBTOP+I+ISP))**2+
     '            XIG(2,ng,NBTOP+I+ISP)**2)
            ENDDO
            ISP=ISP+NDET(nb,2)
            DO I=1,NDET(nb,3)
              !Integral of 1/(r*r) for each part of split element
              !with singularity at node 3
              TEST(NBTOP+7,I)=TEST(NBTOP+7,I)+
     '            WG(ng,NBTOP+I+ISP)*DET(nb,3,ng,I)/
     '            DSQRT((1.0D0-XIG(2,ng,NBTOP+I+ISP))**2+
     '            XIG(1,ng,NBTOP+I+ISP)**2)
            ENDDO
            ISP=ISP+NDET(nb,3)
            DO I=1,NDET(nb,4)
              !Integral of 1/(r*r) for each part of split element
              !with singularity at node 4
              TEST(NBTOP+8,I)=TEST(NBTOP+8,I)+
     '            WG(ng,NBTOP+I+ISP)*DET(nb,4,ng,I)/
     '            DSQRT((1.0D0-XIG(1,ng,NBTOP+I+ISP))**2+
     '            (1.0D0-XIG(2,ng,NBTOP+I+ISP))**2)
            ENDDO
          ENDDO !End of ng loop
          DO nb_incr=1,4
            WRITE(OP_STRING,*)
     '         ' Area of 1st part of split element   =',
     '       TEST(NBTOP+nb_incr,    1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '         ' Area of 2nd part of split element =',
     '         TEST(NBTOP+nb_incr,    2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
          DO nn=1,4
            WRITE(OP_STRING,*)' Centre of singularity at node ',nn
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '         ' Integral of 1/(r*r) in first part of element = ',
     '         TEST(NBTOP+4+nn,1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            WRITE(OP_STRING,*)
     '         ' Integral of 1/(r*r) in second part of element = ',
     '         TEST(NBTOP+4+nn,2)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF !End of 2d integrals tests
        DO ng=1,NGT(nb)
          WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''): '',3E11.3)')
     '      ng,nb,(XIG(ni,ng,nb),ni=1,NIT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO
CC$      call mp_unsetlock()
      ENDIF
      DO nbbem=NBTOP+1,NBTOP+NBASEF(nb,0)-1
        NFBASE(1,nbbem)=nb  !family number of global basis number nbbem
        NFBASE(2,nbbem)=nbbem-NBTOP+1 !child number in family

        IF(NIT(nb).EQ.1) THEN !1D integrals
          IF(((IBT(1,1,nb).EQ.1.AND.IBT(2,1,nb).EQ.1).OR.
     '      (IBT(1,1,nb).EQ.2.AND.IBT(2,1,nb).EQ.1)).AND.
     '      (nbbem.EQ.NBTOP+3.OR.nbbem.EQ.NBTOP+4)) THEN !Log scheme
            IF(nbbem.EQ.NBTOP+3) THEN
              CALL GAUSSLOG(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,1,
     '          NGAP(1,nbbem),PG(1,1,1,nbbem),WG(1,nbbem),
     '          XIG(1,1,nbbem),ERROR,*9999)
            ELSE IF(nbbem.EQ.NBTOP+4) THEN
              CALL GAUSSLOG(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,2,
     '          NGAP(1,nbbem),PG(1,1,1,nbbem),WG(1,nbbem),
     '          XIG(1,1,nbbem),ERROR,*9999)
            ENDIF
          ELSE
            CALL GAUSS10(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '        NGAP(1,nbbem),ALIM(1,nbbem),BLIM(1,nbbem),PG(1,1,1,nbbem),
     '        WG(1,nbbem),XIG(1,1,nbbem),ERROR,*9999)
          ENDIF
          NGT(nbbem)=1
          DO ni=1,NIT(nb)
            NGT(nbbem)=NGT(nbbem)*NGAP(ni,nbbem)
          ENDDO
        ELSE IF(NIT(nb).EQ.2) THEN !2D integrals
          IF(IBT(2,1,nb).EQ.1.AND.IBT(2,2,nb).EQ.1) THEN
            !Bilinear, bicubic hermite or some combination of
            !linear and cubic hermite
            IF(nbbem-NBTOP.GT.8) THEN
              !The basis functions for NBBEM-NBTOP=1-8 have already
              !been evaluated above.
              CALL GAUSS10(IBT(1,1,nb),IDO(1,1,0,nb),
     '           INP(1,1,nb),nb,NGAP(1,nbbem),ALIM(1,nbbem),
     '           BLIM(1,nbbem),PG(1,1,1,nbbem),WG(1,nbbem),
     '           XIG(1,1,nbbem),ERROR,*9999)
              NGT(nbbem)=1
              DO ni=1,NIT(nb)
                NGT(nbbem)=NGT(nbbem)*NGAP(ni,nbbem)
              ENDDO
            ENDIF
          ELSE !Other 2d options e.g. biquadratic
            WRITE(OP_STRING,*)
     '         ' This 2d BE basis function not yet implemented'
            CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
            GOTO 9999
          ENDIF !End of  2d options.
        ENDIF !End of 1D or 2D integral loop
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          DO ng=1,NGT(nbbem)
            WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''): '','
     '        //' 3E11.3)') ng,nbbem,(XIG(ni,ng,nbbem),
     '         ni=1,NIT(nbbem))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO
CC$        call mp_unsetlock()
        ENDIF
      ENDDO !End of nbbem loop

      CALL EXITS('GAUSS5')
      RETURN
 9999 CALL ERRORS('GAUSS5',ERROR)
      CALL EXITS('GAUSS5')
      RETURN 1
      END


