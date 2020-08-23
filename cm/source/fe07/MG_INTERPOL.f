      SUBROUTINE MG_INTERPOL(na,niq,NAQ,NLQ,NXQ,YQ,ADAPTIVE,ERROR,*)

C#### Subroutine: MG_INTERPOL
C###  Description:
C###    MG_INTERPOL performs multigrid interpolation (prolongation)
C###    of YQ(nq,niq,na) from coarse grid na+1 to fine grid na.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER na,niq,NAQ(NQM,NAM),NLQ(NQM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
      LOGICAL ADAPTIVE
!     Local Variables
      INTEGER INTERP,nq

      CALL ENTERS('MG_INTERPOL',*9999)
      INTERP=10*na+na+1 !is index for interpolating grid na from na+1

      DO nq=1,NQT !loop over fine grid points (level na)
        IF(NAQ(nq,na).EQ.0                                 !nq is in grid na
     '    .AND.(.NOT.ADAPTIVE
     '           .OR.ADAPTIVE.AND.NLQ(nq).EQ.INTERP)) THEN !& needs interpolating

          IF(NAQ(nq,na+1).EQ.0) THEN      !straight copy
            YQ(nq,niq,na)=YQ(nq,niq,na+1)
          ELSE IF(NAQ(nq,na+1).EQ.1) THEN !2 pt interpolation in Xi(1)
            YQ(nq,niq,na)=0.50D0*(YQ(NXQ(-1,1,nq,na),niq,na+1)
     '                           +YQ(NXQ( 1,1,nq,na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.2) THEN !2 pt interpolation in Xi(2)
            YQ(nq,niq,na)=0.50D0*(YQ(NXQ(-2,1,nq,na),niq,na+1)
     '                           +YQ(NXQ( 2,1,nq,na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.3) THEN !2 pt interpolation in Xi(3)
            YQ(nq,niq,na)=0.50D0*(YQ(NXQ(-3,1,nq,na),niq,na+1)
     '                           +YQ(NXQ( 3,1,nq,na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.4) THEN !4 pt interp.n in Xi1,2 plane
            YQ(nq,niq,na)=0.25D0*
     '        (YQ(NXQ(-1,1,NXQ(-2,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ(-1,1,NXQ( 2,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ(-2,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ( 2,1,nq,na),na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.5) THEN !4 pt interp.n in Xi2,3 plane
            YQ(nq,niq,na)=0.25D0*
     '        (YQ(NXQ(-2,1,NXQ(-3,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ(-2,1,NXQ( 3,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 2,1,NXQ(-3,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 2,1,NXQ( 3,1,nq,na),na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.6) THEN !4 pt interp.n in Xi3,1 plane
            YQ(nq,niq,na)=0.25D0*
     '        (YQ(NXQ(-3,1,NXQ(-1,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ(-3,1,NXQ( 1,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 3,1,NXQ(-1,1,nq,na),na),niq,na+1)
     '        +YQ(NXQ( 3,1,NXQ( 1,1,nq,na),na),niq,na+1))
          ELSE IF(NAQ(nq,na+1).EQ.7) THEN !8 pt interpol.n in 1,2,3 space
            YQ(nq,niq,na)=0.125D0*
     '        (YQ(NXQ(-1,1,NXQ(-2,1,NXQ(-3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ(-1,1,NXQ(-2,1,NXQ( 3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ(-1,1,NXQ( 2,1,NXQ(-3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ(-1,1,NXQ( 2,1,NXQ( 3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ(-2,1,NXQ(-3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ(-2,1,NXQ( 3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ( 2,1,NXQ(-3,1,nq,na),na),na),niq,na+1)
     '        +YQ(NXQ( 1,1,NXQ( 2,1,NXQ( 3,1,nq,na),na),na),niq,na+1))
          ENDIF !NAQ
          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_INTERPOL_1)
            WRITE(OP_STRING,
     '       '('' Interpol YQ('',I6,'','',I1,'','',I1,'')='',D12.4)')
     '        nq,niq,na,YQ(nq,niq,na)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_INTERPOL_1)
          ENDIF !DOP
        ENDIF !NAQ
      ENDDO !nq

      CALL EXITS('MG_INTERPOL')
      RETURN
 9999 CALL ERRORS('MG_INTERPOL',ERROR)
      CALL EXITS('MG_INTERPOL')
      RETURN 1
      END


