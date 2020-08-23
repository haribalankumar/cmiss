      SUBROUTINE NIQ_LOC(ACTION,TYPE,niq,SUBTYPE,ERROR,*)

C#### Subroutine: NIQ_LOC
C###  Description:
C###    <HTML>
C###    NIQ_LOC returns niq number.
C###    <PRE>
C###    TYPE is  NIQ_MULTIGRID for multigrid applications
C###             NIQ_GRID for miscellanous grid point information
C###             NIQ_ION for general ionic current equations
C###             NIQ_FHN for Fitzhugh-Nagumo ionic equations
C###             NIQ_VCD for VanCapelle-Durrer ionic equations
C###             NIQ_BR for Beeler-Reuter ionic equations
C###             NIQ_LR for Luo-Rudy ionic equations
C###             NIQ_DN for diFrancesco-Noble ionic equations
C###
C###    SUBTYPE is NIQ_SOLUTION is current solution (multigrid)
C###             NIQ_RESID1 is residual or restriction (multigrid)
C###             NIQ_RESID2 is residual or restriction (multigrid)
C###             NIQ_RESTRIC is restriction (multigrid)
C###             NIQ_SOURCE is source term (multigrid)
C###             NIQ_OLDSOLN is soln at previous time step for
C###                        dynamic eqn (multigrid)
C###             NIQ_TRANSPOT is transmembrane potential (ion)
C###             NIQ_EXTPOT is extracellular potential (ion)
C###             NIQ_RECOV is recovery variable (if needed) (ion)
C###             NIQ_CALCIUM is calcium level (FHN,VCD,BR,LR)
C###             NIQ_SAC is SAC current (FHN)
C###             NIQ_DRECOV is change in recovery variable (VCD)
C###             NIQ_DTRANSPOT is change in transmembrane pot. (VCD)
C###             NIQ_X (BR,LR)
C###             NIQ_M (BR,LR,DN)
C###             NIQ_H (BR,LR,DN)
C###             NIQ_J (BR,LR)
C###             NIQ_D (BR,LR)
C###             NIQ_F (BR,LR)
C###             NIQ_TGRIDACT is activation time (in ms) of grid point
C###               (grid)
C###             NIQ_TRANSAPOTTEMP used for temporary storage of phi(m)
C###               (grid)
C###             NIQ_EXTPOTTEMP used for temporary storage of phi(e)
C###               (grid)
C###             NIQ_MAXDPOT is abs largest value of change in potential
C###               (grid)
C###    ACTION is  NIQ_ALLOCATE to allocate an niq# for a problem type
C###             NIQ_ALLOCATE_AND_LOCK to allocate and lock an niq#
C###                        for a problem type
C###             NIQ_FREE is to free a locked allocated niq#
C###             NIQ_INQUIRE to inquire about the niq# for type/subtype
C###    If ACTION is NIQ_ALLOCATE_AND_LOCK the niq# is
C###    locked until it is freed with another call to NIQ_LOC with an
C###    ACTION of NIQ_FREE.
C###    </PRE>
C###    <HTML>
C***  Martin Buist 15 July 1997

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'nqloc00.cmn'
      INCLUDE 'nqloc00.inc'
!     Parameter List
      INTEGER ACTION,niq,SUBTYPE,TYPE
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER LOCK,niqq
      LOGICAL FOUND

      CALL ENTERS('NIQ_LOC',*9999)

      IF(ACTION.EQ.NIQ_ALLOCATE) THEN
        LOCK=NIQ_NOLOCK
      ELSE IF(ACTION.EQ.NIQ_ALLOCATE_AND_LOCK) THEN
        LOCK=NIQ_LOCK
      ENDIF

      IF(TYPE.EQ.NIQ_MULTIGRID.OR.TYPE.EQ.
     '  NIQ_GRID.OR.TYPE.EQ.NIQ_ION.OR.TYPE.EQ.NIQ_FHN.OR.TYPE.EQ.
     '  NIQ_VCD.OR.TYPE.EQ.NIQ_BR.OR.TYPE.EQ.NIQ_LR.OR.TYPE.EQ.
     '  NIQ_DN.OR.TYPE.EQ.NIQ_JRW) THEN

        IF(ACTION.EQ.NIQ_ALLOCATE.OR.ACTION.EQ.
     '    NIQ_ALLOCATE_AND_LOCK) THEN

          FOUND=.FALSE.
          niq=0
          DO WHILE(niq.LT.NIQM.AND.(.NOT.FOUND))
            niq=niq+1
            IF((NIQ_TYPE(niq).EQ.0).OR.(NIQ_LOCKS(niq).EQ.NIQ_NOLOCK))
     '        FOUND=.TRUE.
          ENDDO
          IF(FOUND) THEN
            NIQ_SUBTYPE(niq)=SUBTYPE
            NIQ_TYPE(niq)=TYPE
            NIQ_LOCKS(niq)=LOCK
            NIQ_LIST(0)=NIQ_LIST(0)+1
            NIQ_LIST(NIQ_LIST(0))=niq
            IF(DOP) THEN
              WRITE(OP_STRING,'('' Alloc. niq= '',I2,'' for type= '','
     '          //'I2,'', subtype= '',I2)') niq,TYPE,SUBTYPE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDIF
          ELSE
            niq=0
            ERROR='>>Increase NIQM'
            GOTO 9999
          ENDIF

        ELSE IF(ACTION.EQ.NIQ_FREE) THEN

          FOUND=.FALSE.
          niqq=0
          DO WHILE(niqq.LT.NIQ_LIST(0).AND.(.NOT.FOUND))
            niqq=niqq+1
            niq=NIQ_LIST(niqq)
            IF((NIQ_TYPE(niq).EQ.TYPE).AND.NIQ_SUBTYPE(niq).EQ.SUBTYPE)
     '        FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) THEN
            ERROR='>>Cannot free niq, none allocated for the '
     '        //'problem type'
            GOTO 9999
          ELSE
            NIQ_LIST(0)=NIQ_LIST(0)-1
            DO niqq=niq,NIQ_LIST(0)
              NIQ_LIST(niqq)=NIQ_LIST(niqq+1)
            ENDDO
            NIQ_LIST(NIQ_LIST(0)+1)=0
            NIQ_SUBTYPE(niq)=0
            NIQ_LOCKS(niq)=NIQ_NOLOCK
            NIQ_TYPE(niq)=0
          ENDIF

        ELSE IF(ACTION.EQ.NIQ_INQUIRE) THEN

          FOUND=.FALSE.
          niqq=0
          DO WHILE(niqq.LT.NIQ_LIST(0).AND.(.NOT.FOUND))
            niqq=niqq+1
            niq=NIQ_LIST(niqq)
            IF((NIQ_TYPE(niq).EQ.TYPE).AND.NIQ_SUBTYPE(niq).EQ.SUBTYPE)
     '        FOUND=.TRUE.
          ENDDO
          IF(.NOT.FOUND) niq=0
          IF(DOP) THEN
            WRITE(OP_STRING,'('' Inquire: niq= '',I2,'' for type= '','
     '        //'I2,'', subtype= '',I2)') niq,TYPE,SUBTYPE
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF

        ELSE
          niq=0
          ERROR='>>Invalid niq action'
          GOTO 9999
        ENDIF
      ELSE
        niq=0
        ERROR='>>Invalid niq type'
        GOTO 9999
      ENDIF

      CALL EXITS('NIQ_LOC')
      RETURN
 9999 CALL ERRORS('NIQ_LOC',ERROR)
      CALL EXITS('NIQ_LOC')
      RETURN 1
      END


