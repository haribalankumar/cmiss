      SUBROUTINE OPLEAD(ERROR,*)

C#### Subroutine: OPLEAD
C###  Description:
C###    OPLEAD outputs electrocardiographic leads.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'lead00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nlead,nleadelec

      CALL ENTERS('OPLEAD',*9999)

      WRITE(OP_STRING,
     '  '(/'' The number of leads is '',I2)') NUMLEADS
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      DO nlead=1,NUMLEADS
        WRITE(OP_STRING,'(/'' Lead '',I2,'': '')') nlead
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The lead title is : '',A)')
     '    LEADTITLE(nlead)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The number of electrodes in the lead '
     '    //'is '',I2)') LEADELECS(0,nlead)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The electrodes are       :'',10I6)')
     '    (LEADELECS(nleadelec,nlead),nleadelec=1,
     '    LEADELECS(0,nlead))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The additive constant is :'',F6.2)')
     '    LEADCOUP(0,nlead)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'('' The coupling coeffs. are :'',10F6.2)')
     '    (LEADCOUP(nleadelec,nlead),nleadelec=1,
     '    LEADELECS(0,nlead))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDDO !nlead

C LKC 24-JUL-98 Changed display lead
C      WRITE(OP_STRING,'(/'' The electrode signal filename is : '',A)')
C     '   SIGFNAME
C      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPLEAD')

      RETURN
 9999 CALL ERRORS('OPLEAD',ERROR)
      CALL EXITS('OPLEAD')
      RETURN 1
      END


