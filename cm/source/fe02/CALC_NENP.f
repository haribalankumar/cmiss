      SUBROUTINE CALC_NENP(NBJ,NEELEM,NENP,NPNE,ERROR,*)

C#### Subroutine: CALC_NENP
C###  Description:
C###    CALC_NENP calculates the list of elements surrounding a node.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),NPNE(NNM,NBFM,NEM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb,ne,nn,noelem,np,nps,npsort,nr

      CALL ENTERS('CALC_NENP',*9999)

C*** Initialise array

      DO np=1,NPT(0)
        DO nr=0,NRT
          NENP(np,0,nr)=0
        ENDDO
      ENDDO

      DO nr=1,NRT
        DO noelem=1,NEELEM(0,nr) !Loop over elements
          ne=NEELEM(noelem,nr)
          nb=NBJ(1,ne)
          DO nn=1,NNT(nb) !Loop over element nodes
            np=NPNE(nn,nb,ne)
            npsort=NENP(np,0,nr)
            DO WHILE (NENP(np,npsort,nr).GT.ne.AND.npsort.GT.0)
C             Insert new element in correct place
              npsort=npsort-1
            ENDDO
            IF(NENP(np,npsort,nr).NE.ne.OR.npsort.EQ.0)THEN
C             Element is not already in list
              NENP(np,0,nr)=NENP(np,0,nr)+1
              CALL ASSERT(NENP(np,0,nr).LE.NEPM,
     '          '>>Increase NEPM',ERROR,*9999)
              DO nps=NENP(np,0,nr),npsort+2,-1
                NENP(np,nps,nr)=NENP(np,nps-1,nr)
              ENDDO
              NENP(np,npsort+1,nr)=ne
            ENDIF
            npsort=NENP(np,0,0)
            DO WHILE (NENP(np,npsort,0).GT.ne.AND.npsort.GT.0)
C             Insert new element in correct place
              npsort=npsort-1
            ENDDO
            IF(NENP(np,npsort,0).NE.ne.OR.npsort.EQ.0)THEN
C             Element is not already in list
              NENP(np,0,0)=NENP(np,0,0)+1
              CALL ASSERT(NENP(np,0,0).LE.NEPM,
     '          '>>Increase NEPM',ERROR,*9999)
              DO nps=NENP(np,0,0),npsort+2,-1
                NENP(np,nps,0)=NENP(np,nps-1,0)
              ENDDO
              NENP(np,npsort+1,0)=ne
            ENDIF
          ENDDO !end nn
        ENDDO !end ne
        IF(DOP)THEN
          WRITE(OP_STRING,'('' NENP('',I4,'',0,'',I1,'')='',I2,'
     '      //''', NENP(np,1..,'',I1,''): '',10(I4,X),'
     '      //':(/32X,10(I4,X)))') np,nr,NENP(np,0,nr),nr,
     '      (NENP(np,ne,nr),ne=1,NENP(np,0,nr))
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO !nr

      CALL EXITS('CALC_NENP')
      RETURN
 9999 CALL ERRORS('CALC_NENP',ERROR)
      CALL EXITS('CALC_NENP')
      RETURN 1
      END


