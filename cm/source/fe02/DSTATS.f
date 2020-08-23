      SUBROUTINE DSTATS(LD,NDDL,NDLT,ERROR,*)

C#### Subroutine: DSTATS
C###  Description:
C###    DSTATS calculates data point book-keeping variables:

C**** NDDL(nd,ne)
C**** NDLT(ne)
C**** from LD(nd).

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER LD(NDM),NDDL(NEM,NDEM),NDLT(NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nd,ne

      CALL ENTERS('DSTATS',*9999)

      CALL ASSERT(USE_DATA.EQ.1,'>> Set USE_DATA to 1',ERROR,*9999)
C LKC 19-OCT-98 looping bounds incorrect
C      DO ne=1,NEM
C        NDLT(ne)=0
C      ENDDO
      DO ne=1,NET(0)
        NDLT(ne)=0
      ENDDO

      DO nd=1,NDT
        ne=LD(nd)
        IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$        call mp_setlock()
          WRITE(OP_STRING,'('' nd='',I4,'' LD(nd)='',I4)') nd,LD(nd)
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$        call mp_unsetlock()
        ENDIF

C LKC 19-OCT-98
C        IF((ne.GT.0).AND.(ne.LE.NEM)) THEN
        IF(ne.EQ.0) THEN
          WRITE(OP_STRING,'(
     '      ''>> WARNING: No element associated with data point'',I5)')
     '      nd
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ELSE IF((ne.GT.0).AND.(ne.LE.NET(0))) THEN
          NDLT(ne)=NDLT(ne)+1
          IF(NDLT(ne).LE.NDEM) THEN   !IJ Le G  6th Jan 1992
            NDDL(ne,NDLT(ne))=nd
          ELSE
            ERROR=' Too many data points in element: NDEM too small'
            GO TO 9999
          ENDIF
        ELSE
          WRITE(OP_STRING,
     '      '('' WARNING : Invalid elem associated with nd'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO

      CALL EXITS('DSTATS')
      RETURN
 9999 CALL ERRORS('DSTATS',ERROR)
      CALL EXITS('DSTATS')
      RETURN 1
      END



