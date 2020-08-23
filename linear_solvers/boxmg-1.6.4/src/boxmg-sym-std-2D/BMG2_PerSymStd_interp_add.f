      SUBROUTINE BMG2_PerSymStd_interp_add(
     &                KC, KF, Q, QC, SOR, CI,
     &                IIC, JJC, IIF, JJF, NSORv, IBC 
     &                )

C
C***BEGIN PROLOGUE  BMG2_PerSymStd_interp_add
C***SUBSIDIARY
C***PURPOSE  BMG2_PerSymStd_interp_add interpolates Q from the coarse mesh
C            KC to the fine mesh KF and adds result to Q on fine mesh.
C
C***LIBRARY   SLATEC
C***AUTHOR  DENDY, J. E. JR.
C             LOS ALAMOS NATIONAL LABORATORY
C           VOYTKO, M. H.
C             LAWRENCE LIVERMORE NATIONAL LABORATORY
C
C          Copyright, 1988. The Regents of the University of California.
C          This software was produced under a U.S. Government contract
C          (W-7405-ENG-36) by Los Alamso National Laboratory, which is
C          operated by the University of California for the U.S. Depart-
C          ment of Energy. The U.S. Government is licensed to use,
C          reproduce, and distribute this software. Permission is
C          granted to the public to copy and use this software without
C          charge, provided that this Notice and any statement of
C          authorship are reproduced on all copies. Neither the
C          Government nor the University makes any warranty, expressed
C          or implied, or assumes any liability or responsibility for
C          the use of this software.
C
C***PARAMETERS
C***INPUT
C   KC        KC is grid number for coarse grid.
C   KF        KF is grid number for fine grid.
C   QC        Q for coarse grid.
C   SOR       Refer to BOXMG.
C   CI        Refer to BOXMG.
C   IIC       Number of grid points in x direction on coarse grid,
C             including two fictitious points.
C   JJC       Number of grid points in y direction on coarse grid,
C             including two fictitious points.
C   IIF       Number of grid points in x direction on fine grid,
C             including two fictitious points.
C   JJF       Number of grid points in y direction on fine grid,
C             including two fictitious points.
C   IPN       Refer to BOXMG.
C***INPUT/OUTPUT
C   Q         Refer to BOXMG.
C
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830925  DATE WRITTEN
C   900627  modified to conform to the 4/10/90 "Guide to the SLATEC
C           Common Mathematical Library" by Victor A. Bandy
C
C***END PROLOGUE  MGIADP
C
      IMPLICIT NONE

      INCLUDE 'BMG_constants.h'
      INCLUDE 'BMG_stencils.h'

C     CALLING ARGUMENTS
      integer IIC, IIF, JJC, JJF, IBC, KC, KF, NSORv
      real*8 CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC),
     &       SOR(IIF,JJF,NSORv) 

C     LOCAL ARGUMENTS
      integer ICC, IFF, IICF, IICF1, IIC1, IIF1, JCC, JFF, JJCF, JJCF1, 
     & JJC1, JJF1, IPN
      real*8 A, AQ
C
C   interpolate answers from coarse to fine mesh and add
C   to answers on fine mesh.
C
C***FIRST EXECUTABLE STATEMENT  MGINADP

C
      IPN=IABS(IBC)
      JJF1=JJF-1
      IIF1=IIF-1
      IIC1=IIC-1
      JJC1=JJC-1
      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
      JFF=2
      IFF=2
      Q(2,JFF) = Q(2,JFF) + QC(2,2)

      DO ICC = 3,IICF1

         IFF = IFF + 2
         Q(IFF,JFF) = Q(IFF,JFF) + QC(ICC,2)

         A = CI(ICC,2,LR)*QC( ICC ,2) 
     &     + CI(ICC,2,LL)*QC(ICC-1,2)
         Q(IFF-1,JFF) = Q(IFF-1,JFF) + A + SOR(IFF-1,JFF,MTOT)

      ENDDO

      DO JCC = 3,JJCF1

         JFF = JFF+2
         IFF = 2
         Q(2,JFF) = Q(2,JFF) + QC(2,JCC)

         AQ = CI(2,JCC,LA)*QC(2,JCC) 
     &      + CI(2,JCC,LB)*QC(2,JCC-1)
         Q(2,JFF-1) = Q(2,JFF-1) + AQ + SOR(2,JFF-1,MTOT)

         DO ICC = 3,IICF1

            IFF = IFF + 2
            Q(IFF,JFF) = Q(IFF,JFF) + QC(ICC,JCC)

            A = CI(ICC,JCC,LR)*QC(ICC,JCC) 
     &        + CI(ICC,JCC,LL)*QC(ICC-1,JCC)
            Q(IFF-1,JFF) = Q(IFF-1,JFF) + A + SOR(IFF-1,JFF,MTOT)

            AQ = CI(ICC,JCC,LA)*QC(ICC,JCC) 
     &         + CI(ICC,JCC,LB)*QC(ICC,JCC-1)
            Q(IFF,JFF-1) = Q(IFF,JFF-1) + AQ + SOR(IFF,JFF-1,MTOT)

            A = CI(ICC-1,JCC-1,LSW)*QC(ICC-1,JCC-1)
     &        + CI(ICC-1,JCC-1,LNW)*QC(ICC-1, JCC )
     &        + CI(ICC-1,JCC-1,LNE)*QC( ICC , JCC )
     &        + CI(ICC-1,JCC-1,LSE)*QC( ICC ,JCC-1)
            Q(IFF-1,JFF-1) = Q(IFF-1,JFF-1) + A + SOR(IFF-1,JFF-1,MTOT)

         ENDDO
      ENDDO

      IF( IABS(IBC).NE.1 .AND. IABS(IBC).NE.3 ) GO TO 30
      DO IFF = 1,IIF
         Q(IFF, 1 ) = Q(IFF,JJF1)
         Q(IFF,JJF) = Q(IFF,  2 )
      ENDDO
   30 IF( IABS(IBC).NE.2 .AND. IABS(IBC).NE.3 ) RETURN
      DO JFF = 1,JJF
         Q( 1 ,JFF) = Q(IIF1,JFF)
         Q(IIF,JFF) = Q(  2 ,JFF)
      ENDDO

      RETURN
      END
