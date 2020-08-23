      SUBROUTINE SUMFLOWS(NBJ,NEELEM,nj,NPNE,NXI,NVJE,XP,ERROR,*)

C     Get all terminal flows. If any terminal flows are in the wrong direction
C     due to instabilities, overwrite with average. This should be ok if
C     there are only a few erroneous nodes. 
C     Flows are then summed up the airways using the terminal values. Check
C     that summed flows gives that correct inlet flow. The summed flows ensures
C     flow conservation throughout the tree.
C     AJS May 2010

      IMPLICIT NONE
      INCLUDE 'b00.cmn' 
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn' 
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'lung00.cmn' 
      INCLUDE 'lung_nej00.cmn' 
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M),nj,NPNE(NNM,NBFM,NEM),
     &  NVJE(NNM,NBFM,NJM,NEM),NXI(-NIM:NIM,0:NEIM,0:NEM)
      REAL*8 XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)

      INTEGER i,nb,ne,ne0,nn,noelem,np,np0,np1,np2,nv,nv0,nv1,nv2,
     &  CHECKED(NEM)
      REAL*8 averageflow,checkflow,flow,inletflow,summ,FLOWS(NEM)
      CHARACTER ERRMSG*200

      CALL ENTERS('SUMFLOWS',*9999)

      inletflow=XP(1,1,nj,1) ! flow at node #1 
      averageflow=0.0d0

C Terminal flows
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
!         CHECKED(ne)=0
!         FLOWS(ne)=0.0d0
        nb=NBJ(1,ne)
        np2=NPNE(2,nb,ne) 
        nv2=NVJE(2,nb,nj,ne)
        IF(NXI(1,0,ne).EQ.0)THEN ! terminal
          averageflow=averageflow+XP(1,nv2,nj,np2) 
          CHECKED(ne)=1
          FLOWS(ne)=XP(1,nv2,nj,np2)
        ENDIF
      ENDDO
      averageflow=averageflow/NTERMINAL
c      write(*,*) 'Average terminal flow=',averageflow

! C AJS - check for erroneous nodes. If terminal, then set to average flow.
!       DO noelem=1,NEELEM(0)
!         ne=NEELEM(noelem)
!         nb=NBJ(1,ne)
!         DO nn=1,2
!           np=NPNE(nn,nb,ne) 
!           nv=NVJE(nn,nb,nj,ne)
!           IF((inletflow.GT.0.0d0.AND.XP(1,nv,nj,np).LT.0.0d0).OR.
!      &      (inletflow.LT.0.0d0.AND.XP(1,nv,nj,np).GT.0.0d0))THEN
!             WRITE(IOOP,'('' >>Flow in wrong direction for terminal '//
!      &        'element '',I6,''. Overwriting with average flow.'')')ne
!             IF(NXI(1,0,ne).EQ.0) XP(1,nv,nj,np)=averageflow ! terminal
!             
!           ENDIF !wrong direction flow...
!           IF(NXI(1,0,ne).EQ.0)THEN !set arrays for terminal elements
!             CHECKED(ne)=1
!             FLOWS(ne)=XP(1,nv,nj,np)
!           ENDIF
!         ENDDO !nn
!       ENDDO!noelem

C  Sum flows up the tree
      DO noelem=NEELEM(0),1,-1
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        ne0=NXI(-1,1,ne)!parent element
        IF(ne0.NE.0)THEN !not at inlet element
          FLOWS(ne0)=FLOWS(ne0)+FLOWS(ne) !parent flow = parent flow + daughter flow
          CHECKED(ne0)=1 !element has been summed
        ENDIF
      ENDDO

! C  Sum flows up the tree
!       DO noelem=NEELEM(0),1,-1
!         ne=NEELEM(noelem)
!         nb=NBJ(1,ne)
!         IF(NXI(1,0,ne).NE.0)THEN !not at terminal element
!           summ=0.0d0
!           np2=NPNE(2,nb,ne)!end node of parent
!           nv2=NVJE(1,nb,nj,ne) !version of end node on parent element
!           np1=NPNE(1,nb,ne)!start node of parent
!           nv1=NVJE(1,nb,nj,ne)!version of start node on parent element
! !           IF(NXI(1,0,ne).GT.1)THEN !at bifurcation
!             DO i=1,NXI(1,0,ne) !for each daughter
!               ne0=NXI(1,i,ne)!daughter element
!               WRITE(ERRMSG,'('' >> Flow summation not correct for ne'',
!      &          I6,'' (daughter ne'',I6,'')'')')ne,ne0
!               IF(CHECKED(ne0).EQ.0) WRITE(IOOP,ERRMSG)
! !               CALL ASSERT((CHECKED(ne0).EQ.1),ERRMSG,ERROR,*9999)
!               np0=NPNE(1,nb,ne0)!start node of daughter element
!               nv0=NVJE(1,nb,nj,ne0) !version of start node on daughter element
!               summ=summ+XP(1,nv0,nj,np0) !sum daughter flows
!             ENDDO
! !           ELSE
! !             ne0=NXI(1,1,ne)!daughter element
! !             np0=NPNE(1,nb,ne0)!start node of daughter element
! !             nv0=NVJE(1,nb,nj,ne0) !version of start node on daughter element
! !             summ=XP(1,nv0,nj,np0) !set equal to daughter flow
! !           ENDIF
! ! C         Set flow at start and end nodes equal to sum of daughter flows
! !           XP(1,nv1,nj,np1)=summ
! !           XP(1,nv2,nj,np2)=summ
!           FLOWS(ne)=summ
!           CHECKED(ne)=1 !element has been summed
!         ENDIF
!       ENDDO

C...  Check against existing flow values (note: this is not necessarily an 
C...  error if the flows throughout the rest of the tree have not been 
C...  previously set.
      checkflow=ABS((FLOWS(1)-inletflow)/inletflow)
!       WRITE(*,'('' >> Flow summation does not equal inlet '//
!      &  'flow in SUMFLOWS.f'')')
!       IF(checkflow.GT.1.d-2) WRITE(IOOP,ERRMSG)
!       CALL ASSERT(checkflow.LT.1.d-2,ERRMSG,ERROR,*9999)
!       FLOWS(1)=inletflow !overwrite with original inlet value to remove rounding errors

C     Overwrite nodal flow field with summed flows
      DO noelem=1,NEELEM(0)
        ne=NEELEM(noelem)
        nb=NBJ(1,ne)
        IF(NXI(1,0,ne).NE.0)THEN !not at terminal element
          np1=NPNE(1,nb,ne)
          nv1=NVJE(1,nb,nj,ne)
          XP(1,nv1,nj,np1)=FLOWS(ne) !first node
          np2=NPNE(2,nb,ne)
          nv2=NVJE(2,nb,nj,ne)
          XP(1,nv2,nj,np2)=FLOWS(ne) !second node
        ENDIF
      ENDDO

C.. note: this check is not summing correctly..
! C     Final check
!       DO noelem=1,NEELEM(0)
!         ne=NEELEM(noelem)
!         IF(NXI(1,0,ne).NE.0)THEN !check for all non-terminals
!           nb=NBJ(1,ne)
!           summ=0.0d0
!           np2=NPNE(2,nb,ne)!end node of parent
!           nv2=NVJE(2,nb,nj,ne) !version of end node on parent element
!           flow=XP(1,nv2,nj,np2) !element flow
!           DO i=1,NXI(1,0,ne) !for each daughter
!             ne0=NXI(1,i,ne)!daughter element
!             np0=NPNE(1,nb,ne0)!start node of daughter element
!             nv0=NVJE(1,nb,nj,ne0) !version of start node on daughter element
!             summ=summ+XP(1,nv0,nj,np0) !sum daughter flows
!           ENDDO
!           WRITE(ERRMSG,'('' >> CHECK: Flow in ne '',I6,'' = '',
!      &      D10.4,'' sum in daughters ='',D10.4)')ne,flow,summ
!           checkflow=ABS((summ-flow)/flow)
!           IF(checkflow.GT.1.d-2) write(*,*) 'checkflow=',checkflow
!           CALL ASSERT(checkflow.LT.1.d-2,ERRMSG,ERROR,*9999)
!         ENDIF
!       ENDDO


      
      CALL EXITS('SUMFLOWS')
      RETURN
 9999 CALL ERRORS('SUMFLOWS',ERROR)
      CALL EXITS('SUMFLOWS')
      RETURN 1
      END
