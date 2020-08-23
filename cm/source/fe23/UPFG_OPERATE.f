      SUBROUTINE UPFG_OPERATE(ATPOINTS,CONST,CP,FROM,INDICES,NBH,nc,
     &  NEELEM,NELISTL,NKH,NKJ,NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,NYNP,
     &  OPERATION,SCALE,TO,XAB,XP,YG,YP,ZD,ERROR,*)

C#### Subroutine: UPFG-OPERATE
C###  Description:
C###    This subroutine excecutes the operation (substitute,add,subtract
C###    multiply and divide)


      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp30.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER INDICES(10,2),NBH(NHM,NCM,NEM),nc,NEELEM(0:NE_R_M,0:NRM),
     '  NELISTL(0:NEM),NKH(NHM,NPM,NCM,0:NRM),NKJ(NJM,NPM),
     &  NPLIST_LOCAL(0:1000),NUMVALUES,NVHP(NHM,NPM,NCM,0:NRM),
     '  NVJP(NJM,NPM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8  CONST,CP(NMM,NPM,NXM),SCALE,XAB(NORM,NEM),
     &  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM),YP(NYM,NIYM,NXM),
     &  ZD(NJM,NDM)
      CHARACTER ATPOINTS*16,FROM*16,OPERATION*16,TO*16,ERROR*(*)
!     Local Variables
      INTEGER FGNK(2,0:NUMVALUES),FGNV(2,0:NUMVALUES),nk,np,nv
      REAL*8 F(NKM,NVM,NUMVALUES),FG(NKM,NVM,NUMVALUES),
     '  G(NKM,NVM,NUMVALUES)
      LOGICAL MISMATCH

      CALL ENTERS('UPFG_OPERATE',*9999)

      FGNK(1,0)=0
      FGNV(1,0)=0
      MISMATCH=.FALSE.
      IF(ATPOINTS(1:5).EQ.'NODES')THEN
        ! Copy TO array in F and FROM array into G
        
        CALL UPFGNODES(TO,'FG',F,FGNK,FGNV,CONST,CP,INDICES,nc,NKH,NKJ,
     '    NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,NYNP,1,1.0D0,XP,YP,
     '    ERROR,*9999)

        CALL UPFGNODES(FROM,'FG',G,FGNK,FGNV,CONST,CP,INDICES,nc,NKH,
     &    NKJ,NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,NYNP,2,SCALE,XP,YP,
     '    ERROR,*9999)

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' Copied TO and FROM arrays to F and G '
     '           //'array'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' F Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(1,np)
              WRITE(OP_STRING,'('' F(np='',I6,'',nv='',I3,'',nk=1..) '
     &          //''',8E10.3)') np,nv,(F(nk,nv,np),
     &          nk=1,FGNK(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
          WRITE(OP_STRING,'(/'' G Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(2,np)
              WRITE(OP_STRING,'('' G(np='',I6,'',nv='',I3,'',nk=1..) '
     '          //''',8E10.3)') np,nv,(G(nk,nv,np),
     &          nk=1,FGNK(2,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF !DOP

        ! Operate on FG, F and G, eg. FG=F+G for addition
        IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
          CALL UPFGSUBSTITUTE(FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
          CALL UPFGADD(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
          CALL UPFGSUBTRACT(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
          CALL UPFGMULTIPLY(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:6).EQ.'DIVIDE')THEN
          CALL UPFGDIVIDE(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ENDIF !OPERATION

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' Operated on F and G to give FG ' 
     '          //'array'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' FG Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(1,np)
              WRITE(OP_STRING,'('' FG(np='',I6,'',nv='',I3,'',nk=1..) '
     '          //''',8E10.3)') np,nv,(FG(nk,nv,np),
     &          nk=1,FGNK(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF !DOP

        ! Copy FG into TO array
        CALL UPFGNODES('FG',TO,FG,FGNK,FGNV,CONST,CP,INDICES,nc,NKH,NKJ,
     '    NPLIST_LOCAL,NUMVALUES,NVHP,NVJP,NYNP,1,1.0D0,XP,YP,
     '    ERROR,*9999)

        ! Checks for mismatch
        IF(MISMATCH.AND.DOP)THEN
          WRITE(OP_STRING,'('' >>WARNING: Mismatched number of '
     '      //'derivatives. The FROM array has less derivatives than '
     '      //'the TO array. May need to update derivatives.'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' Copied FG array into the TO array'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF !DOP

      ELSEIF(ATPOINTS(1:8).EQ.'ELEMENTS')THEN
         ! Copy TO array in F and FROM array into G

        CALL UPFGELEMS(TO,'FG',F,FGNK,FGNV,CONST,INDICES,
     '    NELISTL,NUMVALUES,1,1.0D0,XAB,ERROR,*9999)

        CALL UPFGELEMS(FROM,'FG',G,FGNK,FGNV,CONST,INDICES,
     '    NELISTL,NUMVALUES,2,SCALE,XAB,ERROR,*9999)
        ! Operate on FG, F and G, eg. FG=F+G for addition
        IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
          CALL UPFGSUBSTITUTE(FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
          CALL UPFGADD(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
          CALL UPFGSUBTRACT(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
          CALL UPFGMULTIPLY(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:6).EQ.'DIVIDE')THEN
          CALL UPFGDIVIDE(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ENDIF !OPERATION

        ! Copy FG into TO array
        CALL UPFGELEMS('FG',TO,FG,FGNK,FGNV,CONST,INDICES,
     '    NELISTL,NUMVALUES,1,1.0D0,XAB,ERROR,*9999)

      ELSEIF(ATPOINTS(1:8).EQ.'DATA_PTS')THEN
        ! Copy TO array in F and FROM array into G
        CALL UPFGDATA(TO,'FG',F,FGNK,FGNV,CONST,INDICES,NUMVALUES,1,
     &    1.0D0,ZD,ERROR,*9999)
        CALL UPFGDATA(FROM,'FG',G,FGNK,FGNV,CONST,INDICES,NUMVALUES,2,
     &    SCALE,ZD,ERROR,*9999)

        ! Operate on FG, F and G, eg. FG=F+G for addition
        IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
          CALL UPFGSUBSTITUTE(FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
          CALL UPFGADD(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
          CALL UPFGSUBTRACT(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
          CALL UPFGMULTIPLY(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:6).EQ.'DIVIDE')THEN
          CALL UPFGDIVIDE(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ENDIF !OPERATION

        ! Copy FG into TO array
        CALL UPFGDATA('FG',TO,FG,FGNK,FGNV,CONST,INDICES,NUMVALUES,1,
     &    1.0d0,ZD,ERROR,*9999)

      ELSEIF(ATPOINTS(1:9).EQ.'GAUSS_PTS')THEN
        ! Copy TO array in F and FROM array into G
        CALL UPFGGAUSS(TO,'FG',F,FGNK,FGNV,CONST,INDICES,
     '      NBH,NEELEM,NUMVALUES,1,1.0D0,YG,ERROR,*9999)
        CALL UPFGGAUSS(FROM,'FG',G,FGNK,FGNV,CONST,INDICES,
     '      NBH,NEELEM,NUMVALUES,2,SCALE,YG,ERROR,*9999)

        ! Operate on FG, F and G, eg. FG=F+G for addition
        IF(OPERATION(1:10).EQ.'SUBSTITUTE')THEN
          CALL UPFGSUBSTITUTE(FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:3).EQ.'ADD')THEN
          CALL UPFGADD(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'SUBTRACT')THEN
          CALL UPFGSUBTRACT(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:8).EQ.'MULTIPLY')THEN
          CALL UPFGMULTIPLY(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ELSEIF(OPERATION(1:6).EQ.'DIVIDE')THEN
          CALL UPFGDIVIDE(F,FG,FGNK,FGNV,G,MISMATCH,NUMVALUES,
     &      ERROR,*9999)
        ENDIF !OPERATION

        ! Copy FG into TO array
        CALL UPFGGAUSS('FG',TO,FG,FGNK,FGNV,CONST,INDICES,
     '      NBH,NEELEM,NUMVALUES,1,1.0D0,YG,ERROR,*9999)

        ! Checks for mismatch
        IF(MISMATCH.AND.DOP)THEN
          WRITE(OP_STRING,'('' >>WARNING: Mismatched number of gauss '
     '      //'point per element between TO and FROM arrays.'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
        ENDIF

        IF(DOP)THEN
          WRITE(OP_STRING,'(/'' F Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(1,np)
              WRITE(OP_STRING,'('' F(np='',I6,'',nv='',I3,'',nk=1..) '
     '          //''',8E10.3)') np,nv,(F(nk,nv,np),
     &          nk=1,FGNK(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
          WRITE(OP_STRING,'(/'' G Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(2,np)
              WRITE(OP_STRING,'('' G(np='',I6,'',nv='',I3,'',nk=1..) '
     '          //''',8E10.3)') np,nv,(G(nk,nv,np),
     &          nk=1,FGNK(2,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
          WRITE(OP_STRING,'(/'' FG Array :'')')
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          DO np=1,FGNK(1,0)
            DO nv=1,FGNV(1,np)
              WRITE(OP_STRING,'('' FG(np='',I6,'',nv='',I3,'',nk=1..) '
     '          //''',8E10.3)') np,nv,(FG(nk,nv,np),
     &          nk=1,FGNK(1,np))
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
            ENDDO
          ENDDO
        ENDIF !DOP

      ELSEIF(ATPOINTS(1:7).EQ.'GRIDPTS')THEN
        ! Not implemented
      ENDIF

      CALL EXITS('UPFG_OPERATE')
      RETURN
 9999 CALL ERRORS('UPFG_OPERATE',ERROR)
      CALL EXITS('UPFG_OPERATE')
      RETURN 1
      END


