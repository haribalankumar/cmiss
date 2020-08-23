      SUBROUTINE CHDISP(NHP,NKJ,NPLIST,NPNODE,NRLIST,NVJP,NYNP,PAOPTI,XP
     &     ,YP,STRING,FIX,ERROR,*)
 
C#### Subroutine: CHDISP
C###  Description:
C###    CHDISP allows user to change the displacement value for essential 
C###    (Dirichlet) boundary points based on rotations and translations 
C###    based on the current values

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'ktyp50.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'map000.cmn'
!     Parameter List
      INTEGER NKJ(NJM,NPM),NPLIST(0:NPM),NPNODE(0:NP_R_M,0:NRM)
     &     ,NRLIST(0:NRM),NVJP(NJM,NPM),NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM
     &     ,NRM),nh,NHP(NPM,0:NRM,NXM)
      REAL*8 PAOPTI(*),XP(NKM,NVM,NJM,NPM),YP(NYM,NIYM,NXM)
      CHARACTER ERROR*(*),STRING*(MXCH)
      LOGICAL FIX(NYM,NIYFIXM,NXM)
 
!     Local Variables
      INTEGER I,IBEG,IEND,i1,i2,J,N3CO,nc,nj,nk,nolist,np,nr,nrl,NTRL,nv
     &     ,NXLIST(0:NXM),nx,ny,MAX_VERSIONS
      REAL*8 AMOUNT,ANGLE,AXIS(3),RFROMC,RL1(3),RL2(3),RL3(3),
     &     TRANS(3,4),TRANS2(3,4),YY(3,NVM),Z(3),T(12)
      LOGICAL ABBREV,ALL_REGIONS,CBBREV,ROTATE,TRANSL,UNQUAL,TMATRIX
      
      CALL ENTERS('CHDISP',*9999)
C
C     OR 15/06/2005:  This subroutine is to introduce two new commands 
C     which are capable of modifying the displacement values. The 
C     core strucure this soubroutine code was copied over from CHNODS.f.
C
      nc=1 !temporary (displacement bc,NCT(nr,nx))
      TRANSL=.FALSE.
      ROTATE=.FALSE.
      TMATRIX=.FALSE.
      
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)
        
C---------------------------------------------------------------------

C#### Command: FEM change displacement_bc initialise
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the fixed (Dirichlet, essential) boundary points, which
C###    should be initialised ith zeros
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the fixed (Dirichlet, essential) boundary points,
C###    which are initialized with 0.
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    Initialises displacement values with zeros

        OP_STRING(1)=STRING(1:IEND)//' initialise '
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------
        
C#### Command: FEM change displacement_bc rotate by ANGLE
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the fixed (Dirichlet, essential) boundary points,
C###    which are modified by a rigid body rotation
C###  Parameter:      <about AX#[0.0],AY#[0.0],AZ#[0.0]>
C###    Specify a point though which the axis of rotation passes
C###  Parameter:      <axis BX#[0.0],BY#[0.0],BZ#[1.0]>
C###    Specify the direction vector of the axis of rotation
C###  Parameter:      <region (#s/all)[1]>
C###    Specify the regions to rotate
C###  Description:
C###    Rotates specified node values by ANGLE about axis
C###    through point (AX,AY,AZ) with direction vector (BX,BY,BZ).
C###    and computes the displacements from its starting values.
C###    The displacement values are put in the displacement field.

        OP_STRING(1)=STRING(1:IEND)//' rotate by ANGLE'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<about AX#[0.0],AY#[0.0],AZ#[0.0]>'
        OP_STRING(4)=BLANK(1:15)//'<axis BX#[0.0],BY#[0.0],BZ#[1.0]>'
        OP_STRING(5)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)


C---------------------------------------------------------------------

C#### Command: FEM change displacement_bc translate by DX#,DY#,DZ#
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the fixed (Dirichlet, essential) boundary points,
C###    which are modified by a rigid body translation
C###  Parameter:      <region (#s/all)[1]>
C###  Description:
C###    Uses DX,DY,DZ as displacement values for nodes. The option all
C###    refers to all essential (Dirichlet) boundary points

        OP_STRING(1)=STRING(1:IEND)//' translate by DX#,DY#,DZ#'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

C#### Command: FEM change dispplacement_bc tmatrix <3x4 Matrix entries>
C###  Parameter:      <node (#s/all)[all]>
C###    Specify the fixed (Dirichlet, essential) boundary points,
C###    which are modified by a rigid body translation
C###  Parameter:      <region (#s/all)[1]>
C###  Description: 
C###        Applys to the fixed (Dirichlet, essential) boundary points
C###        within the specified region the specified transformation matrix
C###        in three dimendions. The entries for the transformation matrix
C###        are in the following order:
C###        r11, r12, r13, t1, r21, r22, r23, t2, r31, r32, r33, t3
C###        where r stands for the entries of the rotation matrix
C###        and t for the translation entries.
        
        OP_STRING(1)=STRING(1:IEND)//' tmatrix'
        OP_STRING(2)=BLANK(1:15)//'<node (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<region (#s/all)[1]>'
        OP_STRING(4)=BLANK(1:15)//'<matrix (12 comma-separated '
     &           //'numbers) >'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)
C
C---------------------------------------------------------------------
C       
C        
      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHDISP',ERROR,*9999)
      ELSE
        UNQUAL=.TRUE.
        CALL PARSE_CLASS(noco,NTCO,NXLIST,CO,ERROR,*9999)
        nx=NXLIST(1)
        IF(ABBREV(CO(noco+1),'INITIALISE',1)) THEN
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     &         ERROR,*9999)
          DO nrl=1,NRLIST(0)
            nr=NRLIST(nrl)
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              DO nh=1,NHP(np,nr,nx)
                DO nv=1,NVJP(nh,np)
                  DO nk=1,NKJ(nh,np)           
                    ny=NYNP(nk,nv,nh,np,0,nc,nr) 
                    IF (ny.GT.0) THEN
                      IF (FIX(ny,1,nc)) THEN
                        YP(ny,2,nx)=0.0D0
                      ENDIF     ! (FIX...)
                    ENDIF       ! ny>0
                  ENDDO         ! nv
                ENDDO           ! nk
              ENDDO             ! nj=1,NJT
            ENDDO               ! nolist
          ENDDO                 ! nrl
        ENDIF
C
        IF(UNQUAL.AND.ABBREV(CO(noco+1),'ROTATE',1).
     $       OR.ABBREV(CO(noco+1),'TRANSLATE',2).
     &       OR.ABBREV(CO(noco+1),'TMATRIX',2)) THEN
C
          CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,
     &         *9999)
          CALL PARSE_NODES(NPNODE,NPLIST,noco,NRLIST,NTCO,CO,
     &         ERROR,*9999)
C
          CALL RESET(TRANS)
          CALL RESET(TRANS2)
          TRANSL=.FALSE.
          ROTATE=.FALSE.
          IF(ABBREV(CO(noco+1),'TRANSLATE',2)) THEN
            TRANSL=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL1,ERROR,*9999)
            ENDIF
            CALL ASSERT(NTRL.EQ.NJT,'>> # increments entered must '
     $           //'equal # global coordinates',ERROR,*9999)
            AMOUNT=0.d0
            DO nj=1,NJT
              AMOUNT=AMOUNT+RL1(nj)**2
            ENDDO               !nj
            DO nj=NJT+1,3
              RL1(nj)=0.d0
            ENDDO               !nj
            AMOUNT=DSQRT(AMOUNT)
            CALL SHIFT(RL1,AMOUNT,TRANS,ERROR,*9999)
          ELSE IF(ABBREV(CO(noco+1),'ROTATE',2)) THEN
            ROTATE=.TRUE.
            IF(CBBREV(CO,'BY',1,noco+1,NTCO,N3CO)) THEN
              ANGLE=RFROMC(CO(N3CO+1))*PI/180.d0
            ELSE IF(CBBREV(CO,'FROM',2,noco+1,NTCO,N3CO)) THEN
              ANGLE=PAOPTI(1)*PI/180.d0
            ELSE
              ANGLE=0.d0
            ENDIF
            IF(DOP) THEN
C     KAT 14May01: Can't branch out of critical section.
C     Critical section is not essential.
CC$            call mp_setlock()
              WRITE(OP_STRING,'('' Angle='',D12.4)') ANGLE
              CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$            call mp_unsetlock()
            ENDIF
            IF(CBBREV(CO,'ABOUT',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL2,ERROR,*9999)
            ELSE
              RL2(1)=0.d0
              RL2(2)=0.d0
              RL2(3)=0.d0
            ENDIF
            IF(CBBREV(CO,'AXIS',2,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),3,NTRL,RL3,ERROR,*9999)
            ELSE
              RL3(1)=0.d0
              RL3(2)=0.d0
              RL3(3)=1.d0
            ENDIF
            AMOUNT=DSQRT(RL2(1)**2+RL2(2)**2+RL2(3)**2)
            IF(AMOUNT.GT.1.0D-6) THEN
              CALL SHIFT(RL2,-AMOUNT,TRANS,ERROR,*9999)
            ENDIF
            AXIS(1)=RL3(1)
            AXIS(2)=RL3(2)
            AXIS(3)=RL3(3)
            CALL TWIST(AXIS,ANGLE,TRANS,ERROR,*9999)
            CALL TWIST(AXIS,ANGLE,TRANS2,ERROR,*9999)
            IF(AMOUNT.GT.1.0D-6) THEN
              CALL SHIFT(RL2,AMOUNT,TRANS,ERROR,*9999)
            ENDIF
          ELSE IF (ABBREV(CO(noco+1),'TMATRIX',2)) THEN
            TMATRIX=.TRUE.
            DO I=1,12 
              T(I)=0.0d0 
            ENDDO
            IF(CBBREV(CO,'MATRIX',1,noco+1,NTCO,N3CO)) THEN
              CALL PARSRL(CO(N3CO+1),12,NTRL,T,ERROR,*9999)
              IF(NTRL.NE.12) THEN
                WRITE(OP_STRING,'('' >>WARNING: Less then 12 entries ' 
     &               //'for the transformation matrix are specified! ' 
     '               //'The missing numbers get initialised by 0.'')')
                CALL WRITES(IOER,OP_STRING,ERROR,*9999)
              ENDIF
              DO I=1,3
                DO J=1,3
                  TRANS(I,J)=T(4*(I-1)+J)
                  TRANS2(I,J)=T(4*(I-1)+J)
                ENDDO
                TRANS(I,J)=T(4*(I-1)+4)
              ENDDO
            ELSE
              CALL RESET(TRANS)
              CALL RESET(TRANS2)
              WRITE(OP_STRING,'('' >>WARNING: No transformation '
     &             // 'matrix specified, reset transformation '
     &             // 'matrix to identity. '')')
              CALL WRITES(IOER,OP_STRING,ERROR,*9999)
            ENDIF
          ENDIF 
C          
          IF(TRANSL.OR.ROTATE.OR.TMATRIX) THEN
            IF(DOP) THEN
C     KAT 14May01: Can't branch out of critical section.
C     Critical section is not essential.
CC$            call mp_setlock()
              DO I=1,3
                WRITE(OP_STRING,'('' TRANS('',I1,'',j):'',4E12.3)')
     $               I, (TRANS(I,J),J=1,4)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDDO
CC$            call mp_unsetlock()
            ENDIF
C            
            MAX_VERSIONS=1
            DO nolist=1,NPLIST(0)
              np=NPLIST(nolist)
              DO nj=1,NJT
                IF(NVJP(nj,np).GT.MAX_VERSIONS) THEN
                  MAX_VERSIONS=NVJP(nj,np)
                ENDIF
              ENDDO             !nj
            ENDDO               !np


            DO nrl=1,NRLIST(0)
              nr=NRLIST(nrl)
              DO nolist=1,NPLIST(0)
                np=NPLIST(nolist)
                CALL NODE_CHANGE(np,.FALSE.,ERROR,*9999)
                DO nj=1,NJT-1
                  IF(NVJP(nj,np).NE.NVJP(nj+1,np)) THEN
                    WRITE(OP_STRING,'(''Warning: version '
     '                   //'numbers inconsistent for'
     '                   //'transforms'')')
                    CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  ENDIF
                ENDDO              
                CALL ASSERT(MAX_VERSIONS.LE.NVM,
     '               'Increase dim. of YY in CHDISP',
     '               ERROR,*9999) 
                nk=1            ! First we treat its function values, later 
                                ! on, we treat the dirvatives
                DO nv=1,MAX_VERSIONS
                  DO nj=1,NJT
                    IF (nv.LE.NVJP(nj,np)) THEN
                      ny=NYNP(nk,nv,nj,np,0,nc,nr) 
                      YY(nj,nv)=YP(ny,1,nx)+YP(ny,2,nx)!XP(nk,nv,nj,np)!+YP(ny,2,nx)
                    ELSE
                      ny=NYNP(nk,nv,nj,np,0,nc,nr) 
                      YY(nj,nv)=YP(ny,1,nx)+YP(ny,2,nx)!XP(nk,1,nj,np)!+YP(ny,2,nx)
                    ENDIF
                  ENDDO         ! nj             
                  Z(1)=YY(1,nv)
                  Z(2)=YY(2,nv)
                  Z(3)=YY(3,nv)
                  DO i2=1,NJT
                    YY(i2,nv)=0.0D0
                    DO i1=1,NJT
                      YY(i2,nv)=YY(i2,nv)+Z(i1)*TRANS2(i2,i1)
                    ENDDO
                    YY(i2,nv)=YY(i2,nv)+TRANS(i2,4)
                  ENDDO         ! i2
                  DO nj=1,NJT
                    nk=1        
                    ny=NYNP(nk,nv,nj,np,0,nc,nr)  
                    IF (KTYP58(nr).EQ.1) THEN ! then geometric coordinates
C     YY(nj,nv)=YY(nj,nv)-XP(nk,nv,nj,np)-YP(ny,2,nx)
                      YY(nj,nv)=YY(nj,nv)-YP(ny,1,nx)-YP(ny,2,nx)
                    ENDIF
                    IF (FIX(ny,1,nc)) THEN !Update only essential (Dirichlet) 
                                ! boundary points
                      YP(ny,2,nx)=YP(ny,2,nx)+YY(nj,nv)
                    ENDIF
                  ENDDO         ! nj
                ENDDO           ! nv
                
C     Next, we treat all the derivative values. Note, x_(new) = Rot*x_(old) + Trans,
C     then, dx_(new)/ds_1 = Rot*dx_old/ds_1
C
                IF(NKJ(1,np).GT.1) THEN 
                  DO nk=2,NKJ(1,np)
                    DO nv=1,MAX_VERSIONS
                      DO nj=1,NJT
                        IF(nv.LE.NVJP(nj,np)) THEN ! an nv version of XP(nk) exists
                          ny=NYNP(nk,nv,nj,np,0,nc,nr)
                          YY(nj,nv)=XP(nk,nv,nj,np)
                          YY(nj,nv)=YP(ny,1,nx)+YP(ny,2,nx)!XP(nk,nv,nj,np)!+YP(ny,2,nx)
                        ELSE    ! an nv version of XP does not exist for this nk..
                          ny=NYNP(nk,nv,nj,np,0,nc,nr) 
                          YY(nj,nv)=YP(ny,1,nx)+YP(ny,2,nx)!XP(nk,1,nj,np)!+YP(ny,2,nx)
                                !so use 1st version
                        ENDIF
                      ENDDO     ! nj
                      Z(1)=YY(1,nv)
                      Z(2)=YY(2,nv)
                      Z(3)=YY(3,nv)
                      DO i2=1,NJT
                        YY(i2,nv)=0.0D0
                        DO i1=1,NJT
                          YY(i2,nv)=YY(i2,nv)+Z(i1)*TRANS2(i2,i1)
                        ENDDO
                      ENDDO     ! i2
                      DO nj=1,NJT
                        ny=NYNP(nk,nv,nj,np,0,nc,nr)         
                        IF (KTYP58(nr).EQ.1) THEN ! then geometric coordinates
C     YY(nj,nv)=YY(nj,nv)-XP(nk,nv,nj,np)-YP(ny,2,nx)
                          YY(nj,nv)=YY(nj,nv)-YP(ny,1,nx)-YP(ny,2,nx)
                        ENDIF
                        IF (FIX(ny,1,nc)) THEN !Update only essential (Dirichlet) 
                                ! boundary points
                          YP(ny,2,nx)=YP(ny,2,nx)+YY(nj,nv)
                        ENDIF
                      ENDDO     ! nj
                    ENDDO       ! nv
                  ENDDO         ! nk
                ENDIF           ! if NKJ(1,np)>1 ...
              ENDDO             ! np
            ENDDO               ! nrlist
          ENDIF  
        ENDIF                   ! unqualified
      ENDIF
      
      CALL EXITS('CHDISP')
      RETURN
 9999 CALL ERRORS('CHDISP',ERROR)
      CALL EXITS('CHDISP')
      RETURN 1
      END
