      SUBROUTINE GAUSS6(IBT,IDO,INP,nb,NBTOP,NDET,NGAP,DET,PG,WG,XIG,
     '  ERROR,*)

C#### Subroutine: GAUSS6
C###  Description:
C###    GAUSS6 defines parameters for a BE family of hermite simplex
C###    basis functions. The basis number, nb, passed to this routine
C###    is interpreted as the parent basis number of the boundary
C###    element family.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'bem000.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
     '  nb,NBTOP,NDET(NBFM,0:NNM),NGAP(NIM,NBM)
      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,nb1,nbbem,nb_incr,ng,NGHIGH,NGHIGH2,NGLOW,NGLOW2,
     '  ni,nk,nn,ns,nu,POSITION
      REAL*8 ALIM(2,99),BLIM(2,99),PSI2_HERMITE,PSI5,TEST(20,5)
      LOGICAL HERMSECTOR

      CALL ENTERS('GAUSS6',*9999)

c cpb 30/11/96 Generalising this for arbitrary sectors.

      CALL ASSERT(NNT(nb).EQ.3,'>>Not implemented',ERROR,*9999)

      HERMSECTOR=IBT(1,1,nb).EQ.3.AND.IBT(2,1,nb).EQ.4

c cpb 11/7/95 Initialise NDET and DET
      DO nn=0,NNM
        NDET(nb,nn)=1
        DO ng=1,NGM
          DO i=1,6
            DET(nb,nn,ng,i)=1.0d0
          ENDDO !i
        ENDDO !ng
      ENDDO !nn
      NGLOW=NGLIMITS(1,nb,1) !Low order scheme in first Xi direction
      NGHIGH=NGLIMITS(1,nb,2)!High order scheme in first Xi direction
      NGLOW2=NGLIMITS(2,nb,1) !Low order scheme in second Xi direction
      NGHIGH2=NGLIMITS(2,nb,2) !High order scheme in second Xi direction

c cpb 30/11/96 Work out number of basis functions in the family. In
c general require a low order scheme, a mid order schem, a general high
c high order scheme, 2 high order schemes for each location of the
c singularity within the element (since the 'square' element is split
c into 2 triangular elements based on the singularity location) and
c a general adaptive scheme.

C Check to see if we have enough basis functions
      CALL ASSERT(NBTOP+3+2*NNT(nb).LE.NBM,'>>Increase NBM',ERROR,
     '  *9999)
C Set the number of Gauss points for the high order schemes
      DO nb_incr=1,2*NNT(nb)
        NGAP(1,NBTOP+nb_incr)=NGHIGH
        NGAP(2,NBTOP+nb_incr)=NGHIGH2
      ENDDO !nb_incr
C Set the number of Gauss points for the mid order schemes
      NGAP(1,NBTOP+2*NNT(nb)+1)=(NGHIGH+NGLOW)/2
      NGAP(2,NBTOP+2*NNT(nb)+1)=(NGHIGH2+NGLOW2)/2
C Set the number of Gauss points for the low order schemes
      NGAP(1,NBTOP+2*NNT(nb)+2)=NGLOW
      NGAP(2,NBTOP+2*NNT(nb)+2)=NGLOW2
C Set the number of Gauss points for the adaptive scheme
      NGAP(1,NBTOP+2*NNT(nb)+3)=NGHIGH !adaptive scheme
      NGAP(2,NBTOP+2*NNT(nb)+3)=NGHIGH2
C Check to see if we have enough Gauss points
      NGT(nb)=1
      DO ni=1,NIT(nb)
        NGT(nb)=NGT(nb)*NGAP(ni,nb)
      ENDDO !ni
      CALL ASSERT(NGT(nb).LE.NGM,'>>Increase NGM',ERROR,*9999)
C Set the gauss point locations and limits of integration
      ALIM(1,nb)=0.0d0
      BLIM(1,nb)=1.0d0
      ALIM(2,nb)=0.0d0
      BLIM(2,nb)=1.0d0
      DO nb_incr=1,3+2*NNT(nb)
        ALIM(1,NBTOP+nb_incr)=0.0d0
        BLIM(1,NBTOP+nb_incr)=1.0d0
        ALIM(2,NBTOP+nb_incr)=0.0d0
        BLIM(2,NBTOP+nb_incr)=1.0d0
      ENDDO !nb_incr

      CALL GAUSS12(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,NGAP(1,nb),
     '  ALIM(1,nb),BLIM(1,nb),PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),
     '  HERMSECTOR,ERROR,*9999)

C Now calculate the preimage of the Gauss points in the (psi1',psi2')
C coordinate system.  The mapping from a split element in the
C (psi1,psi2) space to the (psi1',psi2') space depends on the location
C of the node at which the singularity is placed.  The determinant of
C this map for the ith part of the split element is stored in
C DET(nb,nn,ng,i). Once the preimage of the Gauss points have been
C found then evaluate the basis function at these preimage points.

C Set the number of splits of each element (i.e. 2) and the determinents
      DO nn=1,NNT(nb)
        NDET(nb,nn)=2
        DO ng=1,NGT(nb)
          DET(nb,nn,ng,1)=XIG(1,ng,nb) !psi1'
          DET(nb,nn,ng,2)=XIG(2,ng,nb) !psi2'
        ENDDO !ng
      ENDDO !nn

C For the collapsed nodes just use the transformation from the first
C element split.

C cpb 30/11/96 Generalising for sectors
C      !Find preimage of Gauss points.
C      IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
C        DO ng=1,NGT(nb)
C          !Singularity at nn=1, first part of split element
C          XIG(1,ng,NBTOP+1)=XIG(1,ng,nb) !psi1'
C          XIG(2,ng,NBTOP+1)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          !Singularity at nn=1, second part of split element
C          XIG(1,ng,NBTOP+2)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          XIG(2,ng,NBTOP+2)=XIG(2,ng,nb) !psi2'
C          !Singularity at nn=2, first part of split element
C          XIG(1,ng,NBTOP+3)=XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+3)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=2, second part of split element
C          XIG(1,ng,NBTOP+4)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
C          XIG(2,ng,NBTOP+4)=1.0D0-XIG(2,ng,nb)
C          !Singularity at nn=3, first part of split element
C          XIG(1,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=3, second part of split element
C          XIG(1,ng,NBTOP+6)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          XIG(2,ng,NBTOP+6)=1.0D0-XIG(2,ng,nb)
C        ENDDO !ng
C      ELSE IF(NKT(3,nb).EQ.1)THEN !Apex at node 3
C        DO ng=1,NGT(nb)
C          !Singularity at nn=1, first part of split element
C          XIG(1,ng,NBTOP+1)=XIG(1,ng,nb) !psi1'
C          XIG(2,ng,NBTOP+1)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          !Singularity at nn=1, second part of split element
C          XIG(1,ng,NBTOP+2)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          XIG(2,ng,NBTOP+2)=XIG(2,ng,nb) !psi2'
C          !Singularity at nn=2, first part of split element
C          XIG(1,ng,NBTOP+3)=1.0D0-XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+3)=XIG(1,ng,nb)*(1.0D0-XIG(2,ng,nb))
C          !Singularity at nn=2, second part of split element
C          XIG(1,ng,NBTOP+4)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          XIG(2,ng,NBTOP+4)=XIG(2,ng,nb)
C          !Singularity at nn=3, first part of split element
C          XIG(1,ng,NBTOP+5)=XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=3, second part of split element
C          XIG(1,ng,NBTOP+6)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
C          XIG(2,ng,NBTOP+6)=1.0D0-XIG(2,ng,nb)
C        ENDDO !ng
C      ENDIF

      DO nn=1,NNT(nb)
        nbbem=NBTOP+2*(nn-1)+1
        DO ng=1,NGT(nb)
          POSITION=INP(nn,1,nb)+2*(INP(nn,2,nb)-1)
          IF(POSITION.EQ.1) THEN !singularity at 'nn=1'
C           First part of split element, xi1=s, xi2=st
            XIG(1,ng,nbbem)=XIG(1,ng,nb)
            XIG(2,ng,nbbem)=XIG(1,ng,nb)*XIG(2,ng,nb)
C           Second part of split element, xi1=st, xi2=t
            XIG(1,ng,nbbem+1)=XIG(1,ng,nb)*XIG(2,ng,nb)
            XIG(2,ng,nbbem+1)=XIG(2,ng,nb)
          ELSE IF(POSITION.EQ.2) THEN !singularity at 'nn=2'
C           First part of split element, xi1=1-s, xi2=s(1-t)
            XIG(1,ng,nbbem)=1.0D0-XIG(1,ng,nb)
            XIG(2,ng,nbbem)=XIG(1,ng,nb)*(1.0D0-XIG(2,ng,nb))
C           Second part of split element, xi1=1-st, xi2=t
            XIG(1,ng,nbbem+1)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            XIG(2,ng,nbbem+1)=XIG(2,ng,nb)
          ELSE IF(POSITION.EQ.3) THEN
C           First part of split element, xi1=s, xi2=1-st
            XIG(1,ng,nbbem)=XIG(1,ng,nb)
            XIG(2,ng,nbbem)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C           Second part of split element, xi1=t(1-s), xi2=1-t
            XIG(1,ng,nbbem+1)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
            XIG(2,ng,nbbem+1)=1.0D0-XIG(2,ng,nb)
          ELSE IF(POSITION.EQ.4) THEN
C           First part of split element, xi1=1-s,xi2=1-st
            XIG(1,ng,nbbem)=1.0D0-XIG(1,ng,nb)
            XIG(2,ng,nbbem)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C           Second part of split element, xi1=1-st, xi2=1-t
            XIG(1,ng,nbbem+1)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
            XIG(2,ng,nbbem+1)=1.0D0-XIG(2,ng,nb)
          ELSE
            ERROR='>>Invalid position'
            GOTO 9999
          ENDIF
          WG(ng,nbbem)=WG(ng,nb)
          WG(ng,nbbem+1)=WG(ng,nb)
        ENDDO !ng
        NGT(nbbem)=NGT(nb)
        NGT(nbbem+1)=NGT(nb)
      ENDDO !nn

      DO nb1=1,NNT(nb)
        nbbem=NBTOP+2*(nb1-1)+1
        ns=0
        DO nn=1,NNT(nb)
          DO nk=1,NKT(nn,nb)
            ns=ns+1
            DO nu=1,NUT(nb)
              IF(HERMSECTOR) THEN
                DO ng=1,NGT(nb)
                  PG(ns,nu,ng,nbbem)=PSI2_HERMITE(IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,nk,nn,XIG(1,ng,nbbem))
                  PG(ns,nu,ng,nbbem+1)=PSI2_HERMITE(IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,nk,nn,XIG(1,ng,nbbem+1))
                ENDDO !ng
              ELSE !Sector
                DO ng=1,NGT(nb)
                  PG(ns,nu,ng,nbbem)=PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,nk,nn,XIG(1,ng,nbbem))
                  PG(ns,nu,ng,nbbem+1)=PSI5(IBT(1,1,nb),IDO(1,1,0,nb),
     '              INP(1,1,nb),nb,nu,nk,nn,XIG(1,ng,nbbem+1))
                ENDDO !ng
              ENDIF
            ENDDO !nu
          ENDDO !nk
        ENDDO !nn
      ENDDO !nb1

C Find the Gauss points for the low and mid order schemes

      DO nbbem=2*NNT(nb)+1,2*NNT(nb)+3
        nb1=NBTOP+nbbem
        CALL GAUSS12(IBT(1,1,nb),IDO(1,1,0,nb),INP(1,1,nb),nb,
     '    NGAP(1,nb1),ALIM(1,nb1),BLIM(1,nb1),PG(1,1,1,nb1),
     '    WG(1,nb1),XIG(1,1,nb1),HERMSECTOR,ERROR,*9999)
        NGT(NBTOP+nbbem)=1
        DO ni=1,NIT(nb)
          NGT(nb1)=NGT(nb1)*NGAP(ni,nb1)
        ENDDO !ni
        CALL ASSERT(NGT(nb1).LE.NGM,'>>Increase NGM',ERROR,*9999)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          DO ng=1,NGT(nbbem)
            WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''):'','
     '        //'3(1X,D12.4))') ng,nb1,(XIG(ni,ng,nb1),ni=1,NIT(nbbem))
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !ng
CC$        call mp_unsetlock()
        ENDIF
      ENDDO !nbbem

C Set the number of schemes
      NBASEF(nb,0)=4+2*NNT(nb) !number of schemes/basis functions
C Set the child and parent numbers
      NFBASE(1,nb)=nb !family number of global basis number nb
      NFBASE(2,nb)=1 !parent is first child
      DO nbbem=1,NBASEF(nb,0)-1
        NFBASE(1,NBTOP+nbbem)=nb !family number of global basis
        NFBASE(2,NBTOP+nbbem)=nbbem+1 !child number in family
        NBASEF(nb,nbbem)=NBTOP+nbbem-1 !global basis number of child
      ENDDO !nbbem
      NBASEF(nb,NBASEF(nb,0))=NBTOP+NBASEF(nb,0)-1
      NBASEF(nb,1)=nb !global basis of parent

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        TEST(1,1)=0.0d0
        TEST(1,2)=0.0d0
        DO ng=1,NGT(nb)
          TEST(1,1)=TEST(1,1)+WG(ng,nb)*XIG(1,ng,nb)
          TEST(1,2)=TEST(1,2)+WG(ng,nb)*XIG(2,ng,nb)
        ENDDO !ng
        WRITE(OP_STRING,'('' Basis number '',I2)') nb
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Integral of xi1 over element = '',D16.8)')
     '    TEST(1,1)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' Integral of xi2 over element = '',D16.8)')
     '    TEST(1,2)
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        DO nn=1,NNT(nb)
          TEST(nn,1)=0.0d0
          TEST(nn,2)=0.0d0
        ENDDO !nn
        DO nn=1,NNT(nb)
          WRITE(OP_STRING,'('' Singularity at node nn = '',I2)') nn
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO i=1,NDET(nb,nn) !Loop over number of element splits
            nb1=NBTOP+2*(nn-1)+i
            DO ng=1,NGT(nb)
              TEST(nn,i)=TEST(nn,i)+WG(ng,nb1)*DET(nb,nn,ng,i)
            ENDDO !ng
            WRITE(OP_STRING,'(''   Area of part '',I1,'' of the split '
     '        //'element = '',D16.8)') i,TEST(nn,i)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !i
        ENDDO !nn
        DO nn=1,NNT(nb)
          TEST(nn,1)=0.0d0
          TEST(nn,2)=0.0d0
        ENDDO !nn
        DO nn=1,NNT(nb)
          WRITE(OP_STRING,'('' Singularity at node '',I2)') nn
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          POSITION=INP(nn,1,nb)+2*(INP(nn,2,nb)-1)
          DO i=1,NDET(nb,nn) !Loop over number of element splits
            nb1=NBTOP+2*(nn-1)+i
C           Integrate 1/r for each part of the split element
            DO ng=1,NGT(nb)
              IF(POSITION.EQ.1) THEN !singularity at 'nn=1'
                TEST(nn,i)=TEST(nn,i)+WG(ng,nb1)*DET(nb,nn,ng,i)/
     '            DSQRT(XIG(1,ng,nb1)**2+XIG(2,ng,nb1)**2)
              ELSE IF(POSITION.EQ.2) THEN !singularity at 'nn=2'
                TEST(nn,i)=TEST(nn,i)+WG(ng,nb1)*DET(nb,nn,ng,i)/
     '            DSQRT((1.0d0-XIG(1,ng,nb1))**2+XIG(2,ng,nb1)**2)
              ELSE IF(POSITION.EQ.3) THEN !singularity at 'nn=3'
                TEST(nn,i)=TEST(nn,i)+WG(ng,nb1)*DET(nb,nn,ng,i)/
     '            DSQRT(XIG(1,ng,nb1)**2+(1.0d0-XIG(2,ng,nb1))**2)
              ELSE IF(POSITION.EQ.4) THEN !singularity at 'nn=4'
                TEST(nn,i)=TEST(nn,i)+WG(ng,nb1)*DET(nb,nn,ng,i)/
     '            DSQRT((1.0d0-XIG(1,ng,nb1))**2+
     '            (1.0d0-XIG(2,ng,nb1))**2)
              ELSE
                ERROR='>>Invalid position'
                GOTO 9999
              ENDIF
            ENDDO !ng
            WRITE(OP_STRING,'(''   Integral of 1/r in part '',I1,'' of '
     '      //'the split element = '',D16.8)') i,TEST(nn,i)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDDO !i
        ENDDO !nn
        DO ng=1,NGT(nb)
          WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''): '','
     '      //' 3E11.3)') ng,nb,(XIG(ni,ng,nb),ni=1,NIT(nb))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDDO !ng
CC$      call mp_unsetlock()
      ENDIF

      CALL EXITS('GAUSS6')
      RETURN
 9999 CALL ERRORS('GAUSS6',ERROR)
      CALL EXITS('GAUSS6')
      RETURN 1
      END

c cpb 1/11/96 Old GAUSS6 and GAUSS6_HERMITE


C      SUBROUTINE GAUSS6(IDO,INP,nb,NBTOP,NDET,NGAP,
C     '  DET,PG,WG,XIG,ERROR,*)
C
CC#### Subroutine: GAUSS6
CC###  Description:
CC###    GAUSS6 defines parameters for a BE family of hermite simplex
CC###    basis functions.
C
CC**** The basis number, nb, passed to this routine is interpreted as
CC**** the parent basis number of the boundary element family.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:b01.cmn'
C      INCLUDE 'cmiss$reference:bem000.cmn'
C      INCLUDE 'cmiss$reference:cbdi02.cmn'
C      INCLUDE 'cmiss$reference:geom00.cmn'
C!     Parameter List
C      INTEGER IDO(NKM,NNM,0:NIM,NBFM),INP(NNM,NIM,NBFM),
C     '  nb,NBTOP,NDET(NBFM,0:NNM),NGAP(NIM,NBM)
C      REAL*8 DET(NBFM,0:NNM,NGM,6),PG(NSM,NUM,NGM,NBM),
C     '  WG(NGM,NBM),XIG(NIM,NGM,NBM)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER I,ISP,nbbem,nb_incr,ng,NGHIGH,
C     '  NGHIGH2,NGLOW,NGLOW2,ni,nk,nn,ns,nu
C      REAL*8 PSI2_HERMITE,TEST(100,5)
C
C      CALL ENTERS('GAUSS6',*9999)
C
Cc cpb 11/7/95 Initialise NDET and DET
C      DO nn=0,NNM
C        NDET(nb,nn)=1
C        DO ng=1,NGM
C          DO i=1,6
C            DET(nb,nn,ng,i)=1.0d0
C          ENDDO !i
C        ENDDO !ng
C      ENDDO !nn
C      NGLOW=NGLIMITS(1,nb,1)
C      NGHIGH=NGLIMITS(1,nb,2)
C      NGLOW2=NGLIMITS(2,nb,1)
C      !Low order scheme in second Xi direction
C      NGHIGH2=NGLIMITS(2,nb,2)
C      !High order scheme in second Xi direction
C      NBASEF(nb,0)=10
C      NFBASE(1,nb)=nb  !family number of global basis number nb
C      NFBASE(2,nb)=1   !parent is first child
C      !Require a low, mid, 7 high order schemes and a general
C      !adaptive scheme.  There is 1 general high order scheme
C      !and two high order schemes for each location of the
C      !singularity in the element (square element is split
C      !into 2 triangular elements based on the singularity
C      !location)
C      DO nb_incr=1,6
C        NGAP(1,NBTOP+nb_incr)=NGHIGH
C        NGAP(2,NBTOP+nb_incr)=NGHIGH2
C      ENDDO !nb_incr
C      NGAP(1,NBTOP+7)=(NGHIGH+NGLOW)/2
C      NGAP(2,NBTOP+7)=(NGHIGH2+NGLOW2)/2
C      NGAP(1,NBTOP+8)=NGLOW
C      NGAP(2,NBTOP+8)=NGLOW2
C      NGAP(1,NBTOP+9)=NGHIGH !adaptive scheme
C      NGAP(2,NBTOP+9)=NGHIGH2
CC Need to check the total number of basis functions?
C      DO nbbem=1,NBASEF(nb,0)
C        NBASEF(nb,nbbem)=NBTOP+nbbem-1
C      ENDDO !nbbem
C      NBASEF(nb,1)=nb
C      CALL GAUSS6_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),
C     '  nb,NGAP(1,nb),PG(1,1,1,nb),WG(1,nb),XIG(1,1,nb),ERROR,*9999)
C      NGT(nb)=1
C      DO ni=1,NIT(nb)
C        NGT(nb)=NGT(nb)*NGAP(ni,nb)
C      ENDDO !ni
C      CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
C     '   ERROR,*9999)
C      !Need to calculate the preimage of the Gauss points in
C      !the (psi1',psi2') coordinate system.  The mapping from
C      !a split element in the (psi1,psi2) space to the
C      !(psi1',psi2') space depends on the location of the
C      !node at which the singularity is placed.  The
C      !determinant of this map for the ith part of the
C      !split element is stored in DET(nb,nn,ng,i).
C      !Once the preimage of the Gauss points has been found then
C      !evaluate the basis function at these preimage points.
C      DO nn=1,NNT(nb)
C        NDET(nb,nn)=2
C        DO ng=1,NGT(nb)
C          DET(nb,nn,ng,1)=XIG(1,ng,nb) !psi1'
C          DET(nb,nn,ng,2)=XIG(2,ng,nb) !psi2'
C        ENDDO !ng
C      ENDDO !nn
C      !Find preimage of Gauss points.
C      IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
C        DO ng=1,NGT(nb)
C          !Singularity at nn=1, first part of split element
C          XIG(1,ng,NBTOP+1)=XIG(1,ng,nb) !psi1'
C          XIG(2,ng,NBTOP+1)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          !Singularity at nn=1, second part of split element
C          XIG(1,ng,NBTOP+2)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          XIG(2,ng,NBTOP+2)=XIG(2,ng,nb) !psi2'
C          !Singularity at nn=2, first part of split element
C          XIG(1,ng,NBTOP+3)=XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+3)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=2, second part of split element
C          XIG(1,ng,NBTOP+4)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
C          XIG(2,ng,NBTOP+4)=1.0D0-XIG(2,ng,nb)
C          !Singularity at nn=3, first part of split element
C          XIG(1,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=3, second part of split element
C          XIG(1,ng,NBTOP+6)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          XIG(2,ng,NBTOP+6)=1.0D0-XIG(2,ng,nb)
C        ENDDO !ng
C      ELSE IF(NKT(3,nb).EQ.1)THEN !Apex at node 3
C        DO ng=1,NGT(nb)
C          !Singularity at nn=1, first part of split element
C          XIG(1,ng,NBTOP+1)=XIG(1,ng,nb) !psi1'
C          XIG(2,ng,NBTOP+1)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          !Singularity at nn=1, second part of split element
C          XIG(1,ng,NBTOP+2)=XIG(1,ng,nb)*XIG(2,ng,nb) !psi1'*psi2'
C          XIG(2,ng,NBTOP+2)=XIG(2,ng,nb) !psi2'
C          !Singularity at nn=2, first part of split element
C          XIG(1,ng,NBTOP+3)=1.0D0-XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+3)=XIG(1,ng,nb)*(1.0D0-XIG(2,ng,nb))
C          !Singularity at nn=2, second part of split element
C          XIG(1,ng,NBTOP+4)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          XIG(2,ng,NBTOP+4)=XIG(2,ng,nb)
C          !Singularity at nn=3, first part of split element
C          XIG(1,ng,NBTOP+5)=XIG(1,ng,nb)
C          XIG(2,ng,NBTOP+5)=1.0D0-XIG(1,ng,nb)*XIG(2,ng,nb)
C          !Singularity at nn=3, second part of split element
C          XIG(1,ng,NBTOP+6)=XIG(2,ng,nb)*(1.0D0-XIG(1,ng,nb))
C          XIG(2,ng,NBTOP+6)=1.0D0-XIG(2,ng,nb)
C        ENDDO !ng
C      ENDIF
C      !Evaluate basis functions at each of these points.
C      DO nb_incr=1,6
C        DO ng=1,NGT(nb)
C          ns=0
C          DO nn=1,NNT(nb)
C            DO nk=1,NKT(nn,nb)
C              ns=ns+1
C              DO nu=1,NUT(nb)
C                PG(ns,nu,ng,NBTOP+nb_incr) =PSI2_HERMITE(IDO(1,1,0,nb),
C     '            INP(1,1,nb),nb,nu,nk,nn,XIG(1,ng,NBTOP+nb_incr))
C              ENDDO !nu
C            ENDDO !nn
C          ENDDO !nk
C          NGT(NBTOP+nb_incr)=NGT(nb)
C          WG(ng,NBTOP+nb_incr)=WG(ng,nb)
C        ENDDO !ng
C      ENDDO !nb_incr
C      IF(DOP) THEN
C        TEST(nb,1)=0.0d0
C        TEST(nb,2)=0.0d0
C        DO nb_incr=0,8
C          TEST(NBTOP+nb_incr,1)=0.0d0
C          TEST(NBTOP+nb_incr,2)=0.0d0
C        ENDDO !nb_incr
C        DO ng=1,NGT(nb)
C          TEST(nb,1)=TEST(nb,1)+WG(ng,nb)*XIG(1,ng,nb)
C          TEST(nb,2)=TEST(nb,2)+WG(ng,nb)*XIG(2,ng,nb)
C        ENDDO !ng
C        WRITE(OP_STRING,*) ' Basis number',nb
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,*) ' Integral of xi1 over element=',
C     '     TEST(nb,1)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        WRITE(OP_STRING,*) ' Integral of xi2 over element=',
C     '     TEST(nb,2)
C        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        DO ng=1,NGT(nb)
C          DO nb_incr=1,3
C            TEST(NBTOP+nb_incr,1)=TEST(NBTOP+nb_incr,1)+
C     '        WG(ng,NBTOP+nb_incr)*DET(nb,nb_incr,ng,1)
C            TEST(NBTOP+nb_incr,2)=TEST(NBTOP+nb_incr,2)+
C     '        WG(ng,NBTOP+nb_incr)*DET(nb,nb_incr,ng,2)
C          ENDDO !nb_incr
C          ISP=0
C          DO I=1,NDET(nb,1)
C !Integral of 1/(r*r) for each part of split element
C !with singularity at node 1
C            TEST(NBTOP+4,I)=TEST(NBTOP+4,I)+
C     '        WG(ng,NBTOP+I)*DET(nb,1,ng,I)/
C     '        DSQRT(XIG(1,ng,NBTOP+I)**2+
C     '        XIG(2,ng,NBTOP+I)**2)
C          ENDDO !I
C          ISP=ISP+NDET(nb,1)
C          DO I=1,NDET(nb,2)
C !Integral of 1/(r*r) for each part of split element
C !with singularity at node 2
C            IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
C              TEST(NBTOP+5,I)=TEST(NBTOP+5,I)+
C     '          WG(ng,NBTOP+I+ISP)*DET(nb,2,ng,I)/
C     '          DSQRT((1.0D0-XIG(2,ng,NBTOP+I+ISP))**2+
C     '          XIG(1,ng,NBTOP+I+ISP)**2)
C            ELSE IF(NKT(3,nb).EQ.1)THEN !Apex at node 3
C              TEST(NBTOP+5,I)=TEST(NBTOP+5,I)+
C     '          WG(ng,NBTOP+I+ISP)*DET(nb,2,ng,I)/
C     '          DSQRT((1.0D0-XIG(1,ng,NBTOP+I+ISP))**2+
C     '          XIG(2,ng,NBTOP+I+ISP)**2)
C            ENDIF
C          ENDDO !I
C          ISP=ISP+NDET(nb,2)
C          DO I=1,NDET(nb,3)
C !Integral of 1/(r*r) for each part of split element
C !with singularity at node 3
C            IF(NKT(1,nb).EQ.1)THEN !Apex at node 1
C              TEST(NBTOP+6,I)=TEST(NBTOP+6,I)+
C     '          WG(ng,NBTOP+I+ISP)*DET(nb,3,ng,I)/
C     '          DSQRT((1.0D0-XIG(1,ng,NBTOP+I+ISP))**2+
C     '          (1.0D0-XIG(2,ng,NBTOP+I+ISP))**2)
C            ELSE IF(NKT(3,nb).EQ.1)THEN !Apex at node 3
C              TEST(NBTOP+6,I)=TEST(NBTOP+6,I)+
C     '          WG(ng,NBTOP+I+ISP)*DET(nb,3,ng,I)/
C     '          DSQRT((1.0D0-XIG(2,ng,NBTOP+I+ISP))**2+
C     '          XIG(1,ng,NBTOP+I+ISP)**2)
C            ENDIF
C          ENDDO !I
C        ENDDO !ng
C        DO nb_incr=1,3
C          WRITE(OP_STRING,*)
C     '      ' Area of 1st part of split element   =',
C     '      TEST(NBTOP+nb_incr,1)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,*)
C     '      ' Area of 2nd part of split element =',
C     '      TEST(NBTOP+nb_incr,2)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDDO !nb_incr
C        DO nn=1,3
C          WRITE(OP_STRING,*)' Centre of singularity at node ',nn
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,*)
C     '      ' Integral of 1/(r*r) in first part of element = ',
C     '      TEST(NBTOP+3+nn,1)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          WRITE(OP_STRING,*)
C     '      ' Integral of 1/(r*r) in second part of element = ',
C     '      TEST(NBTOP+3+nn,2)
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDDO !nn
C        DO ng=1,NGT(nb)
C          WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''): '','
C     '      //' 3E11.3)') ng,nb,(XIG(ni,ng,nb),ni=1,NIT(nb))
C          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C        ENDDO !ng
C      ENDIF
C      DO nbbem=NBTOP+1,NBTOP+NBASEF(nb,0)-1
C        NBASEF(nb,nbbem-NBTOP+1)=nbbem
C        NFBASE(1,nbbem)=nb !family number of global basis number nbbem
C        NFBASE(2,nbbem)=nbbem-NBTOP+1 !child number in family
C        IF(nbbem-NBTOP.GT.6) THEN
CC         The basis functions for NBBEM-NBTOP=1-6 have already been
CC         evaluated above.
C          CALL GAUSS6_HERMITE(IDO(1,1,0,nb),INP(1,1,nb),
C     '      nb,NGAP(1,nbbem),PG(1,1,1,nbbem),WG(1,nbbem),
C     '      XIG(1,1,nbbem),ERROR,*9999)
C          NGT(nbbem)=1
C          DO ni=1,NIT(nb)
C            NGT(nbbem)=NGT(nbbem)*NGAP(ni,nbbem)
C          ENDDO !ni
C        ENDIF
C        IF(DOP) THEN
C          DO ng=1,NGT(nbbem)
C            WRITE(OP_STRING,'('' XIG(ni,'',I3,'','',I2,''): '','
C     '        //' 3E11.3)') ng,nbbem,(XIG(ni,ng,nbbem),
C     '        ni=1,NIT(nbbem))
C            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
C          ENDDO !ng
C        ENDIF
C      ENDDO !nbbem
C
C      CALL EXITS('GAUSS6')
C      RETURN
C 9999 CALL ERRORS('GAUSS6',ERROR)
C      CALL EXITS('GAUSS6')
C      RETURN 1
C      END
C
C
C      SUBROUTINE GAUSS6_HERMITE(IDO,INP,nb,NGAP,PG,WG,XIG,
C     '  ERROR,*)
C
CC#### Subroutine: GAUSS6_HERMITE
CC###  Description:
CC###    GAUSS6_HERMITE defines the Gaussian quadrature coords XIG and
CC###    weights WG and evaluates the basis function Gauss point array
CC###    PG for hermite sector boundary elements.  The basis functions
CC###    are almost tensor products of 1d hermite and special quadratic
CC###    basis functions.
C
CC**** Area coordinates are not used.
C
C      IMPLICIT NONE
C      INCLUDE 'cmiss$reference:geom00.cmn'
C!     Parameter List
C      INTEGER IDO(NKM,NNM,0:NIM),INP(NNM,NIM),nb,NGAP(NIM)
C      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
C      CHARACTER ERROR*(*)
C!     Local Variables
C      INTEGER I,J,K,ng,ng1,ng2,ng3,ni,nk,nn,ns,nu
C      REAL*8 D(7,7),PSI2_HERMITE,W(7,7),XI(3),XIGG(7,7,7,3)
C
C      DATA D/7*0.0D0,
C     '      -0.2886751345948130D0,         0.2886751345948130D0,5*0.0D0,
C     '      -0.3872983346207410D0, 0.0D0,  0.3872983346207410D0,4*0.0D0,
C     '      -0.4305681557970260D0,        -0.1699905217924280D0,
C     '       0.1699905217924280D0,         0.4305681557970260D0,3*0.0D0,
C     '      -0.4530899229693320D0,        -0.2692346550528410D0,  0.0D0,
C     '       0.2692346550528410D0,         0.4530899229693320D0,2*0.0D0,
C     '      -0.4662347571015760D0,        -0.3306046932331330D0,
C     '      -0.1193095930415990D0,         0.1193095930415990D0,
C     '       0.3306046932331330D0,         0.4662347571015760D0,  0.0D0,
C     '      -0.4745539561713800D0,        -0.3707655927996970D0,
C     '      -0.2029225756886990D0, 0.0D0,  0.2029225756886990D0,
C     '       0.3707655927996970D0,         0.4745539561713800D0/
C      DATA W/1.0D0,6*0.0D0,2*0.50D0,5*0.0D0,0.2777777777777780D0,
C     '       0.4444444444444440D0,         0.2777777777777780D0,4*0.0D0,
C     '       0.1739274225687270D0,         0.3260725774312730D0,
C     '       0.3260725774312730D0,         0.1739274225687270D0,3*0.0D0,
C     '       0.1184634425280940D0,         0.2393143352496830D0,
C     '       0.2844444444444440D0,         0.2393143352496830D0,
C     '       0.1184634425280940D0, 2*0.0D0,
C     '       0.0856622461895850D0,         0.1803807865240700D0,
C     '       0.2339569672863460D0,         0.2339569672863460D0,
C     '       0.1803807865240700D0,         0.0856622461895850D0,  0.0D0,
C     '       0.0647424830844350D0,         0.1398526957446390D0,
C     '       0.1909150252525600D0,         0.2089795918367350D0,
C     '       0.1909150252525600D0,         0.1398526957446390D0,
C     '       0.064742483084435D0/
C
C      CALL ENTERS('GAUSS6_HERMITE',*9999)
C
C      ng1=NGAP(1)
C      ng2=1
C      ng3=1
C      IF(NIT(nb).GT.1) ng2=NGAP(2)
C      IF(NIT(nb).GT.2) ng3=NGAP(3)
C      DO k=1,ng3
C        DO j=1,ng2
C          DO i=1,ng1
C            XIGG(i,j,k,1)=0.50D0+D(i,ng1)
C            XIGG(i,j,k,2)=0.50D0+D(j,ng2)
C            XIGG(i,j,k,3)=0.50D0+D(k,ng3)
C            ng=i+(j-1+(k-1)*ng2)*ng1
C            WG(ng)=W(i,ng1)*W(j,ng2)*W(k,ng3)
C            DO ni=1,NIT(nb)
C              XI(ni)=XIGG(I,J,K,ni)
C              XIG(ni,ng)=XI(ni)
C            ENDDO
C            ns=0
C            DO nn=1,NNT(nb)
C              DO nk=1,NKT(nn,nb)
C                ns=ns+1
C                DO nu=1,NUT(nb)
C                  PG(ns,nu,ng)=PSI2_HERMITE(IDO,INP,nb,nu,nk,nn,XI)
C                ENDDO
C              ENDDO
C            ENDDO
C          ENDDO !i
C        ENDDO !j
C      ENDDO !k
C
C      CALL EXITS('GAUSS6_HERMITE')
C      RETURN
C 9999 CALL ERRORS('GAUSS6_HERMITE',ERROR)
C      CALL EXITS('GAUSS6_HERMITE')
C      RETURN 1
C      END


