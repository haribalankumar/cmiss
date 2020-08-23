      SUBROUTINE DIMCHK(IUNIT,nx,ERROR,*)

C#### Subroutine: DIMCHK
C###  Description:
C###    DIMCHK checks array dimensions and stops program if any
C###    are exceeded.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IUNIT,nx
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nb

      CALL ENTERS('DIMCHK',*9999)

      ERROR=' '
      IF((NBM.LT.NBT).OR.(NBT.LT.0)) THEN
        WRITE(OP_STRING,'('' NBM='',I4,'' NBT='',I4)') NBM,NBT
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NDM.LT.NDT).OR.(NDT.LT.0)) THEN
        WRITE(OP_STRING,'('' NDM='',I4,'' NDT='',I4)') NDM,NDT
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NEM.LT.NET(1)).OR.(NET(1).LT.0)) THEN
        WRITE(OP_STRING,'('' NEM='',I4,'' NET(1)='',I4)') NEM,NET(1)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NFM.LT.NFT).OR.(NFT.LT.0)) THEN
        WRITE(OP_STRING,'('' NFM='',I4,'' NFT='',I4)') NFM,NFT
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NLM.LT.NLT).OR.(NLT.LT.0)) THEN
        WRITE(OP_STRING,'('' NLM='',I4,'' NLT='',I4)') NLM,NLT
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NOM.LT.NOT(1,1,0,nx)).OR.(NOT(1,1,0,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NOM='',I4,'' NOT(1,1,0,nx)='',I4)') NOM,
     '    NOT(1,1,0,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NOM.LT.NOT(2,1,0,nx)).OR.(NOT(2,1,0,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NOM='',I4,'' NOT(2,1,0,nx)='',I4)') NOM,
     '    NOT(2,1,0,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NPM.LT.NPT(1)).OR.(NPT(1).LT.0)) THEN
        WRITE(OP_STRING,'('' NPM='',I4,'' NPT(1)='',I4)') NPM,NPT(1)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NQM.LT.NQT).OR.(NQT.LT.0)) THEN
        WRITE(OP_STRING,'('' NQM='',I4,'' NQT='',I4)') NQM,NQT
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NYM.LT.NYT(1,1,nx)).OR.(NYT(1,1,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NYM='',I4,'' NYT='',I4)') NYM,NYT(1,1,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NYM.LT.NYT(2,1,nx)).OR.(NYT(2,1,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NYM='',I4,'' NYT='',I4)') NYM,NYT(2,1,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NYM.LT.NYT(1,2,nx)).OR.(NYT(1,2,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NYM='',I4,'' NYT='',I4)') NYM,NYT(1,2,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
      IF((NYM.LT.NYT(2,2,nx)).OR.(NYT(2,2,nx).LT.0)) THEN
        WRITE(OP_STRING,'('' NYM='',I4,'' NYT='',I4)') NYM,NYT(2,2,nx)
        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
        ERROR=' Invalid array dimensions'
      ENDIF
C PJH 14/6/98
C      IF((NZM.LT.NZT(1,nx)).OR.(NZT(1,nx).LT.0)) THEN
C        WRITE(OP_STRING,'('' NZM='',I4,'' NZT(1,nx)='',I6)') NZM,
C     '    NZT(1,nx)
C        CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
C        ERROR=' Invalid array dimensions'
C      ENDIF
      DO nb=1,MIN(NBM,NBFT)
        IF((NAM.LT.NAT(nb)).OR.(NAT(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NAM='',I4,'' NAT('',I2,'')='',I4)')
     '          NAM,nb,NAT(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NGM.LT.NGT(nb)).OR.(NGT(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NGM='',I4,'' NGT('',I2,'')='',I4)')
     '          NGM,nb,NGT(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NIM.LT.NIT(nb)).OR.(NIT(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NIM='',I4,'' NIT('',I2,'')='',I4)')
     '          NIM,nb,NIT(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NKM.LT.NKT(0,nb)).OR.(NKT(0,nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NKM='',I4,'' NKT(0,'',I2,'')='',I4)')
     '          NKM,nb,NKT(0,nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NNM.LT.NNT(nb)).OR.(NNT(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NNM='',I4,'' NNT('',I2,'')='',I4)')
     '          NNM,nb,NNT(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NSM.LT.NST(nb)).OR.(NST(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NSM='',I4,'' NST('',I2,'')='',I4)')
     '          NSM,nb,NST(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
        IF((NUM.LT.NUT(nb)).OR.(NUT(nb).LT.0)) THEN
          WRITE(OP_STRING,'('' NUM='',I4,'' NUT('',I2,'')='',I4)')
     '          NUM,nb,NUT(nb)
          CALL WRITES(IUNIT,OP_STRING,ERROR,*9999)
          ERROR=' Invalid array dimensions'
        ENDIF
      ENDDO
      IF(error.NE.' ') THEN
        GO TO 9999
      ENDIF

      CALL EXITS('DIMCHK')
      RETURN
 9999 CALL EXITS('DIMCHK')
      RETURN 1
      END


