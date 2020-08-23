      SUBROUTINE EGDW_CELLML(EG,ICQS_SPATIAL,IRCQS_SPATIAL,
     '  ne,ng,NQNE,RCQS_SPATIAL,YQS,DW,ERROR,*)

C#### Subroutine: EGDW_CELLML
C###  Description:
C###    EGDW_CELLML evaluates the derivates of the strain energy function
C###    for a given strain at a grid point.
C###    Inputs:
C###      EG - contains Green strain components, it is assumedto be symmetrical
C###      Other parameters are as set up for the problem
C###    Outputs:
C###      DW - contains derivatives of the strain energy function (deviatoric stress components)

      IMPLICIT NONE

      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='EGDW_CELLML')

      INCLUDE 'error0.inc'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'cell02.cmn'
      INCLUDE 'cellml.cmn'
      INCLUDE 'cell_reserved.inc'
      INCLUDE 'ptr00.cmn'

!     Parameter List
      INTEGER ICQS_SPATIAL(NQISVM,NQM),IRCQS_SPATIAL(0:NQRSVM,NQVM),
     &  ne,ng,NQNE(NEQM,NQEM),niqq
      REAL*8 EG(3,3),RCQS_SPATIAL(NQRSVM,NQM),YQS(NIQSM,NQM),DW(6)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nq,CELL_VARIANT_NUM,ncmpt,niqs,nqrsv,
     &  IDXI(6),IDXJ(6),NQLIST(6),i,j
      LOGICAL FOUND

C     Map arrays to find EG indices i,j from DW index n
      DATA IDXI /1,2,3,1,1,2/
      DATA IDXJ /1,2,3,2,3,3/

      CALL ENTERS(ROUTINENAME,*9999)

C     obtain cell variant for nq
      nq=NQNE(ne,ng)
      CELL_VARIANT_NUM=ICQS_SPATIAL(1,nq) 
C     assuming that all cell variants store Tdevij in the
C     same indices of RCQS 3,6,8,4,5,7
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
        RCQS_SPATIAL(niqs,nq) = EG(IDXI(ncmpt),IDXJ(ncmpt))

!         print *,ROUTINENAME,' RCQS_SPATIAL(niqs,nq) ',
!      &    RCQS_SPATIAL(niqs,nq)
      ENDDO

      CALL ZGTG53EVALCELL(%VAL(CELL_ICQS_VALUE_PTR),
     &      %VAL(CELL_RCQS_VALUE_PTR),
     &      ICQS_SPATIAL,
     &      %VAL(IICQS_SPATIAL_PTR),
     &      IRCQS_SPATIAL,
     &      ne,ng,NQNE,
     &      RCQS_SPATIAL,
     &      YQS,ERROR,*9999)

!       NQLIST(1)=2
!       NQLIST(2)=3
!       NQLIST(3)=4
!       NQLIST(4)=5
!       NQLIST(5)=6
!       NQLIST(6)=7
C     convert from 11 12 13 22 23 33 to 11 22 33 12 13 23 format
      NQLIST(1)=2
      NQLIST(2)=5
      NQLIST(3)=7
      NQLIST(4)=3
      NQLIST(5)=4
      NQLIST(6)=6

C     put stresses in DW
      DO niqq=1,6
        niqs=NQLIST(niqq)
        DW(niqq)=YQS(niqs,NQNE(ne,ng))
      ENDDO !niqq (niqs)

!       print *,ROUTINENAME,' DW ',(DW(i), i=1,6)
!       print *,ROUTINENAME,' EG ',((EG(i,j), i=1,3),j=1,3)

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END
