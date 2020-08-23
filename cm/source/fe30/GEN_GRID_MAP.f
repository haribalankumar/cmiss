      SUBROUTINE GEN_GRID_MAP(IBT,IDO,INP,NQET,NQSCNB,NQXI,PGNQE,XIG,
     '  ERROR,*)

C#### Subroutine: GEN_GRID_MAP
C###  Description:
C###    GEN_GRID_MAP is used to create for each scheme, a grid point
C###    to gauss point mapping array and a gauss point to grid point
C###    mapping array. This requires a quadratic basis function to be
C###    defined with the same number of Xi coordinates as the element
C###    and the grid scheme. This is used for local quadratic
C###    interpolation.
C###
C***  NOTE: gauss point to grid point array has not been set up,
C***    is it needed? When do we have a solution at gauss points and
C***    not at the nodes.
C***  Created by Martin Buist, July 1997

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'

!     Parameter list
      INTEGER IBT(3,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  INP(NNM,NIM,NBFM),NQET(NQSCM),NQSCNB(NQSCM),NQXI(0:NIM,NQSCM)
      REAL*8 PGNQE(NGM,NQEM,NQSCM),XIG(NIM,NGM,NBM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER IJF,IJS,IKF,IKS,nb,nbq,ng,ni,ni1,ni2,ni3,NITB,nn,nq,
     '  NQENG(27),nqsc
      REAL*8 DS,DSH,PSI1,XI(3),XIQ(3)
      LOGICAL QUADBASIS

      CALL ENTERS('GEN_GRID_MAP',*9999)

! Loop over all the grid schemes
      DO nqsc=1,NQSCT
! Find a (one/bi/tri) quadratic basis function for the interpolation
        QUADBASIS=.FALSE.
        DO nb=1,NBT
          NITB=NIT(nb)
          IF(NITB.EQ.NQXI(0,nqsc)) THEN
            DO ni=1,NITB
              IF(IBT(1,ni,nb).EQ.1) THEN
                IF(IBT(2,ni,nb).EQ.2) THEN
                  nbq=nb
                  QUADBASIS=.TRUE.
                ELSE
                  QUADBASIS=.FALSE.
                ENDIF
              ENDIF
            ENDDO !ni
          ENDIF
        ENDDO !nb
        CALL ASSERT(QUADBASIS,'>>Must have a quadratic basis defined',
     '    ERROR,*9999)

! nbq is the local quadratic basis function
! nb is the fe basis function
        NITB=NIT(nbq)
        nb=NQSCNB(nqsc)
        DO ni=1,NITB
          XI(ni)=0.0d0
        ENDDO !ni
        DO ni=1,3**NITB
          NQENG(ni)=0
        ENDDO !ni
! Loop over all gauss points from basis nb
        DO ng=1,NGT(nb)
          DO nq=1,NQET(nqsc)
            PGNQE(ng,nq,nqsc)=0.0d0
          ENDDO !nq
! Find the grid point which is closest to the current gauss point
! in xi space
          DSH=RMAX
          IJF=1
          IKF=1
          IJS=1
          IKS=1
          IF(NITB.GE.2) THEN
            IJF=NQXI(2,nqsc)-1
            IJS=2
          ENDIF
          IF(NITB.EQ.3) THEN
            IKF=NQXI(3,nqsc)-1
            IKS=2
          ENDIF
! Loop over the internal grid points, if we restrict ourselves to
! the internal points in each scheme then we can always make a
! local quadratic element without going over global element boundaries.
          DO ni3=IKS,IKF
            DO ni2=IJS,IJF
              DO ni1=2,NQXI(1,nqsc)-1
                XI(1)=DBLE(ni1-1)/DBLE(NQXI(1,nqsc)-1)
                IF(IJS.NE.1) XI(2)=DBLE(ni2-1)/DBLE(NQXI(2,nqsc)-1)
                IF(IKS.NE.1) XI(3)=DBLE(ni3-1)/DBLE(NQXI(3,nqsc)-1)
                DS=0.0d0
                DO ni=1,NITB
                  DS=DS+DABS(XIG(ni,ng,nb)-XI(ni))**2
                ENDDO !ni
                DS=DSQRT(DS)
                IF(DS.LT.DSH) THEN
                  DSH=DS
                  nq=ni1+((ni2-1)*NQXI(1,nqsc))+((ni3-1)*NQXI(1,nqsc)*
     '              NQXI(2,nqsc))  !store local grid point number
                  DO ni=1,NITB
                    XIQ(ni)=XI(ni) !store xi position of nq
                  ENDDO !ni
                ENDIF
              ENDDO !ni1
            ENDDO !ni2
          ENDDO !ni3
          CALL ASSERT(DSH.LT.RMAX,'>>No grid points found',ERROR,*9999)

          IF(DOP) THEN
            WRITE(OP_STRING,'(''Grid: '',I6,'' Gauss: '',I6)') nq,ng
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

! Calculate the xi position of ng relative to the local quad. element
          DO ni=1,NITB
            XI(ni)=0.5d0+((XIG(ni,ng,nb)-XIQ(ni))*(
     '        DBLE(NQXI(ni,nqsc)-1)/2.0d0))
          ENDDO !ni
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Gauss pt: '',I6,'' Xi: '',3F8.6)') ng,
     '        XI(1),XI(2),XI(3)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

! Calulate which local grid point numbers correspond to the local
! quadratic element nodes
          IF(NITB.EQ.1) THEN
            NQENG(1)=nq-1
            NQENG(2)=nq
            NQENG(3)=nq+1
          ELSE IF(NITB.EQ.2) THEN
            NQENG(1)=nq-NQXI(1,nqsc)-1
            NQENG(2)=nq-NQXI(1,nqsc)
            NQENG(3)=nq-NQXI(1,nqsc)+1
            NQENG(4)=nq-1
            NQENG(5)=nq
            NQENG(6)=nq+1
            NQENG(7)=nq+NQXI(1,nqsc)-1
            NQENG(8)=nq+NQXI(1,nqsc)
            NQENG(9)=nq+NQXI(1,nqsc)+1
          ELSE IF(NITB.EQ.3) THEN
            NQENG(1)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(2)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(3)=nq-NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
            NQENG(4)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(5)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(6)=nq-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
            NQENG(7)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(8)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(9)=nq+NQXI(1,nqsc)-(NQXI(1,nqsc)*NQXI(2,nqsc))+1
            NQENG(10)=nq-NQXI(1,nqsc)-1
            NQENG(11)=nq-NQXI(1,nqsc)
            NQENG(12)=nq-NQXI(1,nqsc)+1
            NQENG(13)=nq-1
            NQENG(14)=nq
            NQENG(15)=nq+1
            NQENG(16)=nq+NQXI(1,nqsc)-1
            NQENG(17)=nq+NQXI(1,nqsc)
            NQENG(18)=nq+NQXI(1,nqsc)+1
            NQENG(19)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(20)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(21)=nq-NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))+1
            NQENG(22)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(23)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(24)=nq+(NQXI(1,nqsc)*NQXI(2,nqsc))+1
            NQENG(25)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))-1
            NQENG(26)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))
            NQENG(27)=nq+NQXI(1,nqsc)+(NQXI(1,nqsc)*NQXI(2,nqsc))+1
          ENDIF

          IF(DOP) THEN
            DO ni=1,3**NITB
              WRITE(OP_STRING,'(''Local grid numbers'',I6)') NQENG(ni)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !ni
          ENDIF

! Evaluate the basis functions for the local quadratic element at
! XI which is the position of the gauss point relative to the local
! quadratic element
          DO nn=1,NNT(nbq)
            nq=NQENG(nn)
            PGNQE(ng,nq,nqsc)=PSI1(IBT(1,1,nbq),IDO(1,1,0,nbq),
     '        INP(1,1,nbq),nbq,1,1,nn,XI)
          ENDDO !nn
          IF(DOP) THEN
            DO nq=1,NQET(nqsc)
              WRITE(OP_STRING,'(''Grid weights'',3I6,F12.6)')
     '          ng,nq,nqsc,PGNQE(ng,nq,nqsc)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO !nq
          ENDIF

        ENDDO !gauss point
      ENDDO !scheme

      CALL EXITS('GEN_GRID_MAP')
      RETURN
 9999 CALL ERRORS('GEN_GRID_MAP',ERROR)
      CALL EXITS('GEN_GRID_MAP')
      RETURN 1
      END



