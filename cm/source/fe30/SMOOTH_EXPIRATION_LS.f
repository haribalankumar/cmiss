      SUBROUTINE SMOOTH_EXPIRATION_LS(NBJ,nh,NLIST,NPNE,NVJE,NYNP,XP,YP,
     '  ERROR,*)

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),nh,NLIST(0:300),NPNE(NNM,NBFM,NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NYNP(NKM,NVM,NHM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM),YP(NYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER N,nb,ne,np,ny,ny1
      REAL*8 DISTANCE(300),intercept,slope,XSUM,XXSUM,XYSUM,YSUM
      REAL*8 LENGTH_1D
      
      CALL ENTERS('SMOOTH_EXPIRATION_LS',*9999)


      YSUM=0.d0
      XSUM=0.d0
      XXSUM=0.d0
      XYSUM=0.d0

      !for first node
      ne=NLIST(1)
      np=NPNE(1,nb,ne)
      ny=NYNP(1,1,nh,np)
      YSUM=YP(ny)

      DO N=1,NLIST(0)
        ne=NLIST(N) !don't include first node
        np=NPNE(2,nb,ne)
        ny=NYNP(1,1,nh,np)
        
        IF(N.EQ.1)THEN
          DISTANCE(N)=LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
        ELSE
          DISTANCE(N)=DISTANCE(N-1)+LENGTH_1D(NBJ,ne,NPNE,NVJE,XP)
        ENDIF
        YSUM=YSUM+YP(ny)
        XSUM=XSUM+DISTANCE(N)
        XYSUM=XYSUM+DISTANCE(N)*YP(ny)
        XXSUM=XXSUM+DISTANCE(N)*DISTANCE(N)
      ENDDO !N
      !calculate least squares estimate of straight line thru solution
      slope=(XYSUM-XSUM*YSUM/NLIST(0))/(XXSUM-XSUM*XSUM/(NLIST(0)+1))
      intercept=YSUM/NLIST(0)-slope*XSUM/(NLIST(0)+1)
      
      !update solution values for the branch
      ne=NLIST(1)
      np=NPNE(1,nb,ne)
      ny=NYNP(1,1,nh,np)
      ny1=ny

      intercept=YP(ny)
      ne=NLIST(NLIST(0))
      np=NPNE(2,nb,ne)
      ny=NYNP(1,1,nh,np)

      slope=(YP(ny)-intercept)/DISTANCE(NLIST(0))

      YP(ny)=intercept
      
      DO N=1,NLIST(0)
        ne=NLIST(N)
        np=NPNE(2,nb,ne)
        ny=NYNP(1,1,nh,np)
        YP(ny)=intercept+slope*DISTANCE(N)
      ENDDO !N

      CALL EXITS('SMOOTH_EXPIRATION_LS')
      RETURN
 9999 CALL ERRORS('SMOOTH_EXPIRATION_LS',ERROR)
      CALL EXITS('SMOOTH_EXPIRATION_LS')
      RETURN 1
      END

      
