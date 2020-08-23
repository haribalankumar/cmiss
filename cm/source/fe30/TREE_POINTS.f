      SUBROUTINE TREE_POINTS(ELEM,ELEMLIST,GEOMFRAC,MAX_BRANCH,
     '  MAX_RAN_PTS,NOELEMS,NRANDOM,ntt,POINTS,PWEIGHT,RANDOM_COORD,
     '  XJPOWER,WEIGHTS,ERROR,*)

C#### Subroutine: TREE_POINTS
C###  Description:
C###    TREE_POINTS generates the random points for trees generated
C###    by IPMESH8.
C***  Created by Martin Buist December 1996

      IMPLICIT NONE

      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'

!     Parameter list
      INTEGER MAX_BRANCH,MAX_RAN_PTS,NOELEMS(0:1),ntt
      REAL*8 GEOMFRAC(MAX_BRANCH,2,3),RANDOM_COORD(MAX_RAN_PTS,NJT),
     '  XJPOWER(MAX_BRANCH,3),WEIGHTS(MAX_BRANCH)
      CHARACTER ERROR*(*)
      LOGICAL ELEM,POINTS(MAX_BRANCH)
!     Local variables
      INTEGER ELEMLIST(MAX_BRANCH,4),IBEG,ICHAR,IEND,INFO,ISEED,ne,
     '  nee,nez,nj,NOQUES,np,nrc,NRANDOM
      REAL*8 DX(3),noele,PWEIGHT(MAX_RAN_PTS),CM_RANDOM_NUMBER,XI(3)
      CHARACTER CHAR1*10
      LOGICAL FILEIP

      CALL ENTERS('TREE_POINTS',*9999)

      FILEIP=.FALSE.
      ISEED=2
      NRANDOM=100
      NOQUES=0
      ICHAR=999
      DO nj=1,3
        XI(nj)=0.0d0
      ENDDO !nj

      WRITE(CHAR1,'(I1)') ntt
      CALL STRING_TRIM(CHAR1,IBEG,IEND)
      IDEFLT(1)=100
      FORMAT='($,'' Enter the number of random points '
     '  //'for tree '//CHAR1(IBEG:IEND)//' [100]: '',I6)'
      IF(IOTYPE.EQ.3) IDATA(1)=NRANDOM
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) NRANDOM=IDATA(1)
      CALL ASSERT(NRANDOM.LE.MAX_RAN_PTS,
     '  '>>ERROR need to increase parameter MAX_RAN_PTS in ipmesh8 ',
     '  ERROR,*9999)
      CALL ASSERT(NRANDOM.GT.0,
     '  '>>ERROR must have random points to grow trees ',ERROR,*9999)

      IF(ELEM) THEN
        noele=0.0d0
        DO ne=1,NOELEMS(0)
          IF(POINTS(ne)) THEN
            noele=noele+1.0d0
          ENDIF
        ENDDO !ne
        np=0
        nee=INT(NRANDOM/noele)
        DO nez=1,NOELEMS(0)      !Loop over elements
          IF(POINTS(nez)) THEN
            ne=ELEMLIST(nez,1)
            IF(DOP) THEN
              WRITE(OP_STRING,*) ' '
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,'('' Element number '',I4)') ne
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
            DO nj=1,NJT
              DX(nj)=GEOMFRAC(nez,2,nj)-GEOMFRAC(nez,1,nj)
            ENDDO !nj
            DO nrc=np+1,np+nee
              DO nj=1,NJT
                RANDOM_COORD(nrc,nj)=CM_RANDOM_NUMBER(ISEED)**
     '            (1.0d0/(XJPOWER(ne,nj)+1.0d0))
                XI(nj)=(RANDOM_COORD(nrc,nj)*DX(nj))+GEOMFRAC(nez,1,nj)
                RANDOM_COORD(nrc,nj)=XI(nj)+ELEMLIST(nez,nj+1)
              ENDDO !nj
              PWEIGHT(nrc)=WEIGHTS(nez)
            ENDDO !nrc
            IF(DOP) THEN
              DO nrc=np+1,np+nee
                WRITE(OP_STRING,'(I5,F8.3,F8.3,F8.3)') nrc,
     '            RANDOM_COORD(nrc,1),
     '            RANDOM_COORD(nrc,2),RANDOM_COORD(nrc,3)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO !nrc
            ENDIF
            np=np+nee
          ENDIF
        ENDDO !nez

      ELSE !not elem

        DO nj=1,NJT
          DX(nj)=GEOMFRAC(1,2,nj)-GEOMFRAC(1,1,nj)
        ENDDO !nj
        DO nrc=1,NRANDOM
          DO nj=1,NJT
            RANDOM_COORD(nrc,nj)=CM_RANDOM_NUMBER(ISEED)**(
     '        1.0d0/(XJPOWER(1,nj)+1.0d0))
            RANDOM_COORD(nrc,nj)=RANDOM_COORD(nrc,nj)*DX(nj)
     '        +GEOMFRAC(1,1,nj)
          ENDDO !nj
          PWEIGHT(nrc)=1.0d0
        ENDDO !nrc

      ENDIF !elem

      CALL EXITS('TREE_POINTS')
      RETURN
 9999 CALL ERRORS('TREE_POINTS',ERROR)
      CALL EXITS('TREE_POINTS')
      RETURN 1
      END



