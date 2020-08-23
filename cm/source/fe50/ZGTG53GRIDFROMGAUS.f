      SUBROUTINE ZGTG53GRIDFROMGAUS(AZL,ICQS_SPATIAL,IRCQS_SPATIAL,
     ' ne,ng,NQLIST,NQNE,RCQS_SPATIAL,RET_ERROR,*)

C#### Subroutine: ZGTG53GRIDFROMGAUS
C###  Description:
C###    ZGTG53GRIDFROMGAUS updates 1 grid point at a gauss point with 
C###    green strain components


      IMPLICIT NONE

      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='ZGTG53GRIDFROMGAUS')

      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'

!     Parameter List
      INTEGER ICQS_SPATIAL(NQISVM,NQM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     &  NQLIST(0:NQM),NQNE(NEQM,NQEM)
      REAL*8 AZL(3,3),RCQS_SPATIAL(NQRSVM,NQM)
      CHARACTER RET_ERROR*(*)
!     Local Variables
      INTEGER CELL_VARIANT_NUM,ncmpt,ne,ng,niqs,nq,nqrsv
      REAL*8 DEFORMATION_VALUE
      CHARACTER ERROR*(ERRSTRLEN)
      LOGICAL FOUND
      CALL ENTERS(ROUTINENAME,*9999)

      nq=NQNE(ne,ng) !obtain global grid point number
      CELL_VARIANT_NUM=ICQS_SPATIAL(1,nq) !obtain cell variant for nq
C currently assuming that all cell variants store Tdevij in the 
C same indices of RCQS 3,6,8,4,5,7
      NQLIST(0)=6
      NQLIST(1)=3
      NQLIST(2)=6
      NQLIST(3)=8
      NQLIST(4)=4
      NQLIST(5)=5
      NQLIST(6)=7

C     map RCQS indices to appropriate row in RCQS_SPATIAL array
      DO ncmpt=1,6
        FOUND=.FALSE.
        DO nqrsv=1,IRCQS_SPATIAL(0,CELL_VARIANT_NUM)
          IF(NQLIST(ncmpt).EQ.
     &        IRCQS_SPATIAL(nqrsv,CELL_VARIANT_NUM)) THEN
            NQLIST(ncmpt)=nqrsv
            FOUND=.TRUE.
          ENDIF
        ENDDO !nqrsv
      ENDDO !ncmpt
      CALL ASSERT(FOUND,
     &  '>>A green strain component is not spatially varying '
     &  //'or the wrong RCQS index has been given',
     &  ERROR,*9999)
      DO ncmpt=1,6 !6 green strain components
        niqs=NQLIST(ncmpt)
        DEFORMATION_VALUE=0.0d0
C calculate deformation value (green strain) for grid point
        IF(ncmpt.LE.3) THEN !axial component
          DEFORMATION_VALUE=0.5d0*(AZL(ncmpt,ncmpt)-1)
        ELSE IF(ncmpt.EQ.4) THEN
          DEFORMATION_VALUE=0.5d0*AZL(1,2)
        ELSE IF(ncmpt.EQ.5) THEN
          DEFORMATION_VALUE=0.5d0*AZL(1,3)
        ELSE IF(ncmpt.EQ.6) THEN
          DEFORMATION_VALUE=0.5d0*AZL(2,3)
        ENDIF
        RCQS_SPATIAL(niqs,nq)=DEFORMATION_VALUE
      ENDDO

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      RET_ERROR=ERROR
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END
