      SUBROUTINE MG_RESTRICT(na,niq,NITB,NAQ,NWQ,NXQ,YQ,ERROR,*)

C#### Subroutine: MG_RESTRICT
C###  Description:
C###    MG_RESTRICT is the multigrid restriction operator for
C###    projecting solution (niq=1) or residuals (niq>1) from fine
C###    grid to coarse grid.

C**** YQ(nq,niq,na  ) is fine grid solution(niq=1) or residual(niq>1)
C**** YQ(nq,niq,na+1) is coarse "     "            "     "

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER na,niq,NITB,NAQ(NQM,NAM),NWQ(8,0:NQM,NAM),
     '  NXQ(-NIM:NIM,0:4,0:NQM,NAM)
      REAL*8 YQ(NYQM,NIQM,NAM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ij,ik,MQ,nii,nij,nik,nq
      REAL*8 RQ(-1:1,-1:1,-1:1),SUM,WEIGHT(3)
      DATA WEIGHT/0.5D0,0.25D0,0.125D0/
      DATA RQ/0.125D0,0.250D0,0.125D0, 0.250D0,0.500D0,0.250D0,
     '        0.125D0,0.250D0,0.125D0,
     '        0.250D0,0.500D0,0.250D0, 0.500D0,1.000D0,0.500D0,
     '        0.250D0,0.500D0,0.250D0,
     '        0.125D0,0.250D0,0.125D0, 0.250D0,0.500D0,0.250D0,
     '        0.125D0,0.250D0,0.125D0/

      CALL ENTERS('MG_RESTRICT',*9999)
      DO nq=1,NQT
        IF(NAQ(nq,na+1).EQ.0) THEN !nq is coarse g.p.
          IF(NWQ(1,nq,na+1).ne.0) THEN !nq is boundary point
            YQ(nq,niq,na+1)=YQ(nq,niq,na)
          ELSE                    !nq is internal point
            ik=MAX(0,NITB-2) !zero for 1,2D, one for 3D
            ij=MIN(NITB-1,1) !zero for 1D, one for 2,3D
            SUM=0.0D0
            DO nik=-ik,ik
              DO nij=-ij,ij
                DO nii=-1,1
                  mq=NXQ(nii,1,NXQ(nij*2,1,NXQ(nik*3,1,nq,na),na),na)
!                                                              !adj g.p.
                  IF(mq.GT.0) THEN
                    SUM=SUM+WEIGHT(NITB)*RQ(nii,nij,nik)*YQ(mq,niq,na)
                  ENDIF
                ENDDO !nii
              ENDDO !nij
            ENDDO !nik
            YQ(nq,niq,na+1)=SUM
          ENDIF !Internal point

          IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP CRITICAL(MG_RESTRICT_1)
            WRITE(OP_STRING,
     '       '('' Restrict YQ('',I6,'','',I1,'','',I1,'')='',D12.4)')
     '        nq,niq,na+1,YQ(nq,niq,na+1)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP END CRITICAL(MG_RESTRICT_1)
          ENDIF
        ENDIF !nq in coarse grid
      ENDDO !Global grid point nq

      CALL EXITS('MG_RESTRICT')
      RETURN
 9999 CALL ERRORS('MG_RESTRICT',ERROR)
      CALL EXITS('MG_RESTRICT')
      RETURN 1
      END


