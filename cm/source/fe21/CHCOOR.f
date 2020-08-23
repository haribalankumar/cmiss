      SUBROUTINE CHCOOR(NBJ,NENP,NKJ,NPNE,NPNODE,NUNK,
     '  NVJE,NVJP,XG,XP,STRING,ERROR,*)

C#### Subroutine: CHCOOR
C###  Description:
C###    CHCOOR allows user to redefine coordinate system.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'call00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NENP(NPM,0:NEPM,0:NRM),NKJ(NJM,NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),NUNK(NKM,NJM,NPM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 XP(NKM,NVM,NJM,NPM),XG(NJM,NUM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      REAL*8 XGRC(NJM,NUM)
      INTEGER IBEG,ICOORD,IEND,IFROMC,N3CO,nb,nbref,ne,nep,nj,njnvref,
     '  nk,NKMAX,nn,nonode,np,nr,nu,nv,NVMAX
      LOGICAL CBBREV,CHANGENV

      CALL ENTERS('CHCOOR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change coordinates
C###  Parameter:      <to COORDINATE#[1]>
C###  Description:
C###    Change coordinate system to another predefined system.

        OP_STRING(1)=STRING(1:IEND)//' <to COORDINATE#[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHCOOR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'TO',1,noco+1,NTCO,N3CO)) THEN
          ICOORD=IFROMC(CO(N3CO+1))
        ELSE
          ICOORD=1
        ENDIF
        IF(ICOORD.EQ.1) THEN !change to r.c.
          DO nr=1,NRT
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)
C GMH 8/1/97 Update cmgui link
              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
C             Try do deal with inconsistent nv, nk
              CHANGENV=.FALSE.
              NKMAX=NKJ(1,np)
              NVMAX=NVJP(1,np)
              njnvref=1
              DO nj=2,NJT
                IF(NKJ(nj,np).NE.NKMAX.AND.NKMAX.NE.0) THEN
                  WRITE(OP_STRING,
     '              '('' Differing numbers of derivatives at node '''
     '              //',I2)') np
                  CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                  NKMAX=0
                ENDIF
                IF(NVJP(nj,np).NE.NVMAX) THEN
                  IF(.NOT.CHANGENV) THEN
                    WRITE(OP_STRING,
     '                '('' Differing numbers of versions at node '''
     '                //',I2)') np
                    CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
                    CHANGENV=.TRUE.
                  ENDIF
                  IF(NVJP(nj,np).GT.NVMAX) THEN
                    NVMAX=NVJP(nj,np)
                    njnvref=nj
                  ENDIF
                ENDIF
              ENDDO
              DO nv=1,NVMAX
                DO nj=1,NJT
                  IF(nv.LE.NVJP(nj,np)) THEN
                    DO nk=1,NKMAX
                      nu=NUNK(nk,nj,np)
                      XG(nj,nu)=XP(nk,nv,nj,np)
                    ENDDO
                  ENDIF
                ENDDO
                DO nk=1,NKMAX
                  nu=NUNK(nk,1,np)
C LKC 21-JAN-2003 Changing ZG to a local variable
C                 with compatible dimensions with XG
C                  CALL XZ_DERIV(ITYP10(nr),nu,XG,ZG)
                  CALL XZ_DERIV(ITYP10(nr),nu,XG,XGRC)
                  DO nj=1,NJT
                    XP(nk,nv,nj,np)=XGRC(nj,nu)
                  ENDDO
                ENDDO
              ENDDO !nv
              IF(CHANGENV) THEN
C               Add new versions
                DO nj=1,NJT
                  NVJP(nj,np)=NVMAX
                ENDDO
                IF(CALL_ELEM) THEN
                  DO nep=1,NENP(np,0,nr)
                    ne=NENP(np,nep,nr)
                    nbref=NBJ(njnvref,ne)
C!!!                only works for consistent nn
                    DO nn=1,NNT(nbref)
                      IF(NPNE(nn,nbref,ne).EQ.np) THEN
                        nv=NVJE(nn,nbref,njnvref,ne)
                        DO nj=1,NJT
                          nb=NBJ(nj,ne)
                          IF(NPNE(nn,nb,ne).EQ.np) NVJE(nn,nb,nj,ne)=nv
                        ENDDO
                      ENDIF
                    ENDDO
                  ENDDO !nep
                ENDIF !CALL_ELEM
              ENDIF !change nv
            ENDDO !nonode
          ENDDO !nr
C CPB 26/1/00 adding change to prolate spheroidal from rc. Note: only
C implemented for nodal values, not derivatives.
        ELSE IF(ICOORD.EQ.4) THEN !change to p.s.
          WRITE(OP_STRING,
     '      '(''>>WARNING: Only implmented for nk=1 and nv=1'')')
          CALL WRITES(IOOP,OP_STRING,ERROR,*9999)
          DO nr=1,NRT
            CALL ASSERT(ITYP10(nr).EQ.1,'Can only convert to '
     '        //'prolate spheroidal from rectangular cartesian '
     '        //'coordinates',ERROR,*9999)
            DO nonode=1,NPNODE(0,nr)
              np=NPNODE(nonode,nr)

              CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
              DO nj=1,NJT
C LKC 21-JAN-2003 Change ZG to local variable XGRC with
C               appropriate dimensions
C               ZG(nj,1)=XP(1,1,nj,np)
                XGRC(nj,1)=XP(1,1,nj,np)
              ENDDO !nj
              CALL ZX(ICOORD,XGRC(1,1),XG(1,1))
              DO nj=1,NJT
                XP(1,1,nj,np)=XG(nj,1)
              ENDDO !nj
            ENDDO !nonode
          ENDDO !nr
        ELSE
          ERROR='>>This coordinate change is not yet implemented'
          GO TO 9999
        ENDIF
        ITYP10(1)=ICOORD
      ENDIF

      CALL EXITS('CHCOOR')
      RETURN
 9999 CALL ERRORS('CHCOOR',ERROR)
      CALL EXITS('CHCOOR')
      RETURN 1
      END


