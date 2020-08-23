      SUBROUTINE ConstructNAQ(NITB,NAQ,NXQ,RelativePos,ERROR,*)

C#### Subroutine: ConstructNAQ
C###  Description:
C###    ConstructNAQ constructs the NAQ matrix from the fine grid
C###    NXQ matrix.

      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NITB,NAQ(NQM,NAM),NXQ(-NIM:NIM,0:4,0:NQM,NAM),
     '  RelativePos(0:5,0:NQM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nq,i,j,IX(3),neighbour
      LOGICAL define

      CALL ENTERS('ConstructNAQ',*9999)

! Initialize the RelativePos matrix so that the first column is all -1.
! Also initialize the NAQ matrix so that all entries are -1.
C$OMP PARALLEL DO
C$OMP&PRIVATE(nq,na)
C$OMP&SHARED(NQT,RelativePos,NAQ)
      DO nq=1,NQT
        RelativePos(0,nq)=-1
        DO na=1,NMGT
          NAQ(nq,na)=-1
        ENDDO
      ENDDO
C$OMP END PARALLEL DO

      IF(NMGT.EQ.1) THEN
! Only one multigrid level -> all points belong to this level
        DO nq=1,NQT
          NAQ(nq,1)=0
        ENDDO

      ELSE !more than one multigrid level

! Identify a point that will be a member of the coarsest grid to act as
! a reference and set it up as the first point to be scanned.
        nq=1
        DO i=0,NITB
          RelativePos(i,nq)=0
        ENDDO

! This is the main loop which defines the NVT array.
        DO WHILE(nq.LE.NQT)

! Define the values for the current point nq based on the relative
! position for grid level na.
          DO na=1,NMGT
            DO i=1,3
              IF(i.LE.NITB) THEN
                IX(i) = MOD(ABS(RelativePos(i,nq)),2**(na-1))
              ELSE
                IX(i)=0
              ENDIF
            ENDDO

            define=.TRUE. !Decide whether or not to define NAQ(nq,na)
            if(na.GE.3) then
! PJH 9Jan96  if(NAQ(nq,na-1).GT.0) then
              if(NAQ(nq,na-1).NE.0) then
                define=.FALSE. !i.e. leave NAQ(nq,na) as -1
              endif
            endif

            IF(define) then
              IF((IX(1).EQ.0).AND.(IX(2).eq.0).AND.
     '          (IX(3).EQ.0)) THEN
                NAQ(nq,na)=0
              ELSE IF((IX(1).NE.0).AND.(IX(2).EQ.0).AND.
     '            (IX(3).EQ.0)) THEN
                NAQ(nq,na)=1
              ELSE IF((IX(1).EQ.0).AND.(IX(2).NE.0).AND.
     '            (IX(3).EQ.0)) THEN
                NAQ(nq,na)=2
              ELSE IF((IX(1).EQ.0).AND.(IX(2).EQ.0).AND.
     '            (IX(3).NE.0)) THEN
                NAQ(nq,na)=3
              ELSE IF((IX(1).NE.0).AND.(IX(2).NE.0).AND.
     '            (IX(3).EQ.0)) THEN
                NAQ(nq,na)=4
              ELSE IF((IX(1).EQ.0).AND.(IX(2).NE.0).AND.
     '            (IX(3).NE.0)) THEN
                NAQ(nq,na)=5
              ELSE IF((IX(1).NE.0).AND.(IX(2).eq.0).AND.
     '            (IX(3).NE.0)) THEN
                NAQ(nq,na)=6
              ELSE IF((IX(1).NE.0).AND.(IX(2).NE.0).AND.
     '            (IX(3).NE.0)) THEN
                NAQ(nq,na)=7
              ENDIF
            ENDIF
            IF(DOP) THEN
              WRITE(OP_STRING,'('' NAQ(nq='',I5,'',na='',I2,'')='',I3)')
     '          nq,na,NAQ(nq,na)
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ENDDO

! Find all of the neighbours relative to node nq (on the finest grid),
! set them as possible candidates to be scanned, and calculate their
! positions.
          DO i=-NITB,NITB
            neighbour=NXQ(i,1,nq,1)
            IF((neighbour.NE.0).and.
     '        (RelativePos(0,neighbour).EQ.-1)) THEN
              DO j=0,NITB
                RelativePos(j,neighbour)=RelativePos(j,nq)
              ENDDO
              RelativePos(ABS(i),neighbour)=
     '          RelativePos(ABS(i),neighbour)+i/ABS(i)
            ENDIF
          ENDDO
          RelativePos(0,nq)=1 !We have finished with node nq

! Choose another node to operate upon
          nq=1
          DO WHILE((nq.LE.NQT).AND.(RelativePos(0,nq).NE.0))
            nq=nq+1
          ENDDO
        ENDDO !While (nq.lt.nqt)
      ENDIF !NGMT=1

      CALL EXITS('ConstructNAQ')
      RETURN
 9999 CALL ERRORS('ConstructNAQ',ERROR)
      CALL EXITS('ConstructNAQ')
      RETURN 1
      END


