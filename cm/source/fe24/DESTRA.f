      SUBROUTINE DESTRA(IBT,IDO,INP,LD,NAN,NBJ,NKJE,NPF,NPNE,NRE,
     '  NVJE,SE,XA,XE,XID,XP,ZP,STRING,ERROR,*)

C#### Subroutine: DESTRA
C###  Description:
C###    DESTRA defines principal strains for vector plots

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
!     Parameter List
      INTEGER IBT(3,NIM,NBFM),INP(NNM,NIM,NBFM),IDO(NKM,NNM,0:NIM,NBFM),
     '  LD(NDM),NAN(NIM,NAM,NBFM),NBJ(NJM,NEM),NKJE(NKM,NNM,NJM,NEM),
     '  NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRE(NEM),NVJE(NNM,NBFM,NJM,NEM)
      REAL*8 SE(NSM,NBFM,NEM),XA(NAM,NJM,NEM),XE(NSM,NJM),XID(NIM,NDM),
     '  XP(NKM,NVM,NJM,NPM),ZP(NKM,NVM,NHM,NPM,NCM)
      CHARACTER ERROR*(*),STRING*(*)
!     Local Variables
      INTEGER i,IBEG,IEND,j,k,N3CO,nd,ne
      REAL*8 CM(2,2),CQP(2,2),DBX(2,2),DET,DSX(2,2),DSXINV(2,2),E11,
     '  E11M,E22,E22M,EQP(3,3),EULAN,EVAL(3),FM(3,3),FQP(3,3),
     '  FQPINV(3,3),PFXI,PQP(2,2),PQPM(2,2),RQP(3,3),
     '  THETA,TQP(3,3),UIQP(3,3),UM(3,3),UQP(3,3),XI(3)
      CHARACTER TYPE*8
      LOGICAL ABBREV,CALCU,CBBREV,INVERSE

      CALL ENTERS('DESTRA',*9999)
 1    IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM define strain;c spamm
C###  Parameter:      <(inverse/no-inverse)[no-inverse]>
C###  Description:
C###    Compute strain parameters using knowledge of XP and ZP.


        OP_STRING(1)=STRING(1:IEND)//';c spamm'
        OP_STRING(2)=BLANK(1:15)//'<(inverse/no-inverse)[no-inverse]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe24','doc','DESTRA',ERROR,*9999)
      ELSE
        CALL CHECKQ('C',noco,1,CO,COQU,STRING,*1)
        CALCU =.FALSE.
        IF(ABBREV(COQU(noco,1),'C',1)) THEN
          CALCU=.TRUE.
        ENDIF

        IF(CALCU) THEN
          IF(CBBREV(CO,'SPAMM',2,noco+1,NTCO,N3CO)) THEN
            TYPE='SPAMM'
          ENDIF

          IF(CBBREV(CO,'INVERSE',2,noco+1,NTCO,N3CO)) THEN
            INVERSE=.TRUE.
          ELSE
            INVERSE=.FALSE.
          ENDIF
        ENDIF

        IF(CALCU) THEN
          IF(TYPE(1:5).EQ.'SPAMM') THEN  !Dave Rahdert 29/10/91

C        CALL STRING_TRIM(FILE02,IBEG,IEND)
C        CALL OPENF(9,'DISK',FILE02(IBEG:IEND)//'.fourier','NEW',
C     '    'SEQUEN','FORMATTED',132,ERROR,*9999)
C     It is meant to compute strain parameters using knowledge of
C     XP and ZP.

            DO nd=1,NDT
              ne=LD(nd)
              XI(1)=XID(1,nd)
              XI(2)=XID(2,nd)

              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
              DSX(1,1)=PFXI(IBT,IDO,INP,NAN,1,2,XE(1,1),XI)  !d-X/d-Xi1
              DSX(2,1)=PFXI(IBT,IDO,INP,NAN,1,2,XE(1,2),XI)  !d-Y/d-Xi1
              DSX(1,2)=PFXI(IBT,IDO,INP,NAN,1,4,XE(1,1),XI)  !d-X/d-Xi2
              DSX(2,2)=PFXI(IBT,IDO,INP,NAN,1,4,XE(1,2),XI)  !d-Y/d-Xi2

              CALL XPXE(NBJ(1,ne),NKJE(1,1,1,ne),NPF(1,1),
     '          NPNE(1,1,ne),NRE(ne),NVJE(1,1,1,ne),
     '          SE(1,1,ne),XA(1,1,ne),XE,ZP,ERROR,*9999)
              DBX(1,1)=PFXI(IBT,IDO,INP,NAN,1,2,XE(1,1),XI)  !d-x/d-Xi1
              DBX(2,1)=PFXI(IBT,IDO,INP,NAN,1,2,XE(1,2),XI)  !d-y/d-Xi1
              DBX(1,2)=PFXI(IBT,IDO,INP,NAN,1,4,XE(1,1),XI)  !d-x/d-Xi2
              DBX(2,2)=PFXI(IBT,IDO,INP,NAN,1,4,XE(1,2),XI)  !d-y/d-Xi2

              DET=DSX(1,1)*DSX(2,2)-DSX(1,2)*DSX(2,1)
              DSXINV(1,1)= DSX(2,2)/DET
              DSXINV(1,2)=-DSX(1,2)/DET
              DSXINV(2,1)=-DSX(2,1)/DET
              DSXINV(2,2)= DSX(1,1)/DET       !DBXINV= inverse(DBX)

              DO i=1,2
                DO j=1,2
                  FQP(i,j)=0.0D0
                  DO k=1,2
                    FQP(i,j)=FQP(i,j)+DBX(i,k)*DSXINV(k,j)  !Strain gradient tensor
                  ENDDO
                ENDDO
              ENDDO

C             INVERSE=.TRUE. means that F should be inverted
              IF(INVERSE) THEN
                DET=FQP(1,1)*FQP(2,2)-FQP(1,2)*FQP(2,1)
                FQPINV(1,1)=FQP(2,2)/DET
                FQPINV(2,2)=FQP(1,1)/DET
                FQPINV(1,2)=-FQP(1,2)/DET
                FQPINV(2,1)=-FQP(2,1)/DET

C ...           Overwrite old FQP with FQPINV
                FQP(1,1)=FQPINV(1,1)
                FQP(1,2)=FQPINV(1,2)
                FQP(2,1)=FQPINV(2,1)
                FQP(2,2)=FQPINV(2,2)
              ENDIF

              CQP(1,1)=FQP(1,1)*FQP(1,1)+FQP(2,1)*FQP(2,1)
              CQP(2,1)=FQP(1,2)*FQP(1,1)+FQP(2,2)*FQP(2,1)
              CQP(1,2)=CQP(2,1)
              CQP(2,2)=FQP(1,2)*FQP(1,2)+FQP(2,2)*FQP(2,2)
              CALL EIGEN1(2,2,CQP,EVAL,PQP,ERROR,*9999)
              E11=DSQRT(DABS(EVAL(1)))           !2 principal extensions
              E22=DSQRT(DABS(EVAL(2)))
              UQP(1,1)=E11*PQP(1,1)*PQP(1,1)+E22*PQP(1,2)*PQP(1,2)
              UQP(1,2)=E11*PQP(1,1)*PQP(2,1)+E22*PQP(1,2)*PQP(2,2)
              UQP(2,1)=UQP(1,2)
              UQP(2,2)=E11*PQP(2,1)*PQP(2,1)+E22*PQP(2,2)*PQP(2,2)
              DET=UQP(1,1)*UQP(2,2)-UQP(1,2)*UQP(2,1)
              UIQP(1,1)= UQP(2,2)/DET
              UIQP(1,2)=-UQP(1,2)/DET
              UIQP(2,1)=-UQP(2,1)/DET
              UIQP(2,2)= UQP(1,1)/DET
              RQP(1,1)=FQP(1,1)*UIQP(1,1)+FQP(1,2)*UIQP(2,1)
              RQP(1,2)=FQP(1,1)*UIQP(1,2)+FQP(1,2)*UIQP(2,2)
              RQP(2,1)=FQP(2,1)*UIQP(1,1)+FQP(2,2)*UIQP(2,1)
              RQP(2,2)=FQP(2,1)*UIQP(1,2)+FQP(2,2)*UIQP(2,2)
              EULAN=180.0D0*DATAN2(RQP(2,1),RQP(1,1))/PI  !Euler angle
              EQP(1,1)=0.5D0*(CQP(1,1)-1.0D0)
              EQP(1,2)=0.5D0*CQP(1,2)
              EQP(2,1)=0.5D0*CQP(2,1)
              EQP(2,2)=0.5D0*(CQP(2,2)-1.0D0)

C ...         Now repeat above formulations with F which is referenced
C             to material axes (ie. XI1 base vector is one base.  The
C             second base vector is perpendicular to the first b.v.)
C             All parameters ending with M (such EM) are the same as the
C             corresponding parameter without the M except that it is
C             referenced to the material axes.

              THETA=DATAN(DBX(2,1)/DBX(1,1)) !Angle between XI1 and X axes
              TQP(1,1)=DCOS(THETA)
              TQP(1,2)=-DSIN(THETA)
              TQP(2,1)=DSIN(THETA)
              TQP(2,2)=DCOS(THETA)

              CALL TRAN(1,FQP,FM,TQP)

              CM(1,1)=FM(1,1)*FM(1,1)+FM(2,1)*FM(2,1)
              CM(2,1)=FM(1,2)*FM(1,1)+FM(2,2)*FM(2,1)
              CM(1,2)=CM(2,1)
              CM(2,2)=FM(1,2)*FM(1,2)+FM(2,2)*FM(2,2)

              CALL EIGEN1(2,2,CM,EVAL,PQPM,ERROR,*9999)
              E11M=DSQRT(DABS(EVAL(1)))                   !2 principal extensions
              E22M=DSQRT(DABS(EVAL(2)))
              UM(1,1)=E11M*PQPM(1,1)*PQPM(1,1)+E22M*PQPM(1,2)*PQPM(1,2)
              UM(1,2)=E11M*PQPM(1,1)*PQPM(2,1)+E22M*PQPM(1,2)*PQPM(2,2)
              UM(2,1)=UM(1,2)
              UM(2,2)=E11M*PQPM(2,1)*PQPM(2,1)+E22M*PQPM(2,2)*PQPM(2,2)
              DET=UM(1,1)*UM(2,2)-UM(1,2)*UM(2,1)

C!!! set but never used
C              UIM(1,1)= UM(2,2)/DET
C              UIM(1,2)=-UM(1,2)/DET
C              UIM(2,1)=-UM(2,1)/DET
C              UIM(2,2)= UM(1,1)/DET
C              RM(1,1)=FM(1,1)*UIM(1,1)+FM(1,2)*UIM(2,1)
C              RM(1,2)=FM(1,1)*UIM(1,2)+FM(1,2)*UIM(2,2)
C              RM(2,1)=FM(2,1)*UIM(1,1)+FM(2,2)*UIM(2,1)
C              RM(2,2)=FM(2,1)*UIM(1,2)+FM(2,2)*UIM(2,2)
C              EM(1,1)=0.5D0*(CM(1,1)-1.0D0)
C              EM(1,2)=0.5D0*CM(1,2)
C              EM(2,1)=0.5D0*CM(2,1)
C              EM(2,2)=0.5D0*(CM(2,2)-1.0D0)

              OPEN(IOFILE2,FILE='strain.dat',STATUS='NEW')
              WRITE(IOFILE2,*)'STRAIN PARAMETERS REFERENCED TO ',
     '          'GLOBAL AXES'
              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '        'The deformation gradient F, referenced to global cds is:'
              WRITE(IOFILE2,'(2(F10.5,1X))')FQP(1,1),FQP(1,2)
              WRITE(IOFILE2,'(2(F10.5,1X))')FQP(2,1),FQP(2,2)


              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'The rotation matrix R, referenced to global cds is:'
              WRITE(IOFILE2,'(2(F10.5,1X))')RQP(1,1),RQP(1,2)
              WRITE(IOFILE2,'(2(F10.5,1X))')RQP(2,1),RQP(2,2)

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'The Euler angle  is (degrees):'
              WRITE(IOFILE2,'(A30,F10.5)')' ',EULAN

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'Matrix U, referenced to global cds is:'
              WRITE(IOFILE2,'(2(F10.5,1X))')UQP(1,1),UQP(1,2)
              WRITE(IOFILE2,'(2(F10.5,1X))')UQP(2,1),UQP(2,2)

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'The eigenvectors of U are:'
              WRITE(IOFILE2,'(A30,2(F10.5,1X))')
     '          '(for the first  e-value):',PQP(1,1),PQP(2,1)
              WRITE(IOFILE2,'(A30,2(F10.5,1X))')
     '          '(for the second e-value):',PQP(1,2),PQP(2,2)

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'The principal extensions are:'
              WRITE(IOFILE2,'(A30,2(F10.5,1X))')
     '          '1st,2nd =',E11,E22

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'Matrix C, referenced to global cds is:'
              WRITE(IOFILE2,'(2(F10.5,1X))')CQP(1,1),CQP(1,2)
              WRITE(IOFILE2,'(2(F10.5,1X))')CQP(2,1),CQP(2,2)

              WRITE(IOFILE2,*)' '
              WRITE(IOFILE2,*)
     '          'The strain matrix E, referenced to global cds is:'
              WRITE(IOFILE2,'(2(F10.5,1X))')EQP(1,1),EQP(1,2)
              WRITE(IOFILE2,'(2(F10.5,1X))')EQP(2,1),EQP(2,2)
              CALL CLOSEF(IOFILE2,ERROR,*9999)
C1001          FORMAT(2(F10.5,1X))
C1002          FORMAT(A30,2(F10.5,1X))
C1003          FORMAT(A30,F10.5)

            ENDDO
          ENDIF
        ENDIF

      ENDIF

      CALL EXITS('DESTRA')
      RETURN
 9999 CALL ERRORS('DESTRA',ERROR)
      CALL EXITS('DESTRA')
      RETURN 1
      END


