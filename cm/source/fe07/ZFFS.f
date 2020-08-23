      SUBROUTINE ZFFS(IDOXFT,JDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,
     '  NORMALSIGN1,nr,nx,CG,FS,PG,RDF1,RDF2,XDF,YGF,ZDF,ERROR,*)

C#### Subroutine: ZFFS
C###  Description:
C###    <HTML>
C###    ZFFS calculates face tangent stiffness matrix FS from
C###    current dependent variable array ZE.
C###    <PRE>
C###    KTYP1D=1 : calculation is done algebraically.
C###    KTYP1D=2 : calculation is done by one-sided finite differences.
C###    KTYP1D=3 : calculation is done by central finite differences.
C###    </PRE>
C###    All scale factor work is assumed to be done above this routine.
C###    </HTML>

C#### Variable: FS(ns,idoxf,ns,jdoxf,neax,nh)
C###  Type: REAL*8
C###  Set_up: ZFFS
C###  Description:
C###    FS(ns,idoxf,ns,jdoxf,neax,nh) are the contributions to the tangent
C###    stiffness matrix from one face.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IDOXFT(NHM),JDOXFT,NBHF(NHM),NBJF(NJM,NFM),NEAXT,NHF,
     '  NIEF(0:2),NORMALSIGN1,nr,nx
      REAL*8 CG(NMM,NGM),FS(NSFM,2,NSFM,2,2,NHM),PG(NSM,NUM,NGM,NBM),
     '  RDF1(NSFM,2,NHM),RDF2(NSFM,2,NHM),
     '  XDF(NSFM,2,2,NJM),YGF(NIYGFM,NGFM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER idoxf,jdoxf,ms,nb,NEAFST,neax,nh,nhx,ns,NSTH
      REAL*8 DELTA,HALFDELTA,ZF_STORE
      LOGICAL FLUX

      CALL ENTERS('ZFFS',*9999)

C      IF(.true..or.KTYP1D.NE.1) THEN
      IF(KTYP1D.NE.1) THEN
C***    Finite difference calculation of tangent stiffness matrix

C        IF(.true..or.KTYP1D.EQ.2) THEN !one-sided finite differences
        IF(KTYP1D.EQ.2) THEN !one-sided finite differences
C         Residual is probably close to linear so a large perturbation
C           can be used (to minimize numerical error).
          DELTA=1.0d-4
C         Calculate reference residual
          CALL ZFRF30(IDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,NORMALSIGN1,nr,nx,
     '      CG,PG,RDF1,XDF,YGF,ZDF,ERROR,*9999)
        ELSE !central differences
C         Residual is probably not so linear so a small perturbation
C           must be used.
          HALFDELTA=1.0d-8
          DELTA=2*HALFDELTA
        ENDIF

C!!!  need to check order of stabilizing derivs
        FLUX=ITYP15(nr,nx).NE.3.OR..TRUE.

        DO nhx=1,NHF
          nh=NH_LOC(nhx,nx)
          nb=NBHF(nh)
          NSTH=NST(nb)
          IF(FLUX) THEN
            NEAFST=NEAXT
          ELSE
            NEAFST=1 !derivs are the same in each neighbouring element
          ENDIF
          DO neax=1,NEAFST
            DO jdoxf=1,JDOXFT
              DO ns=1,NSTH
C!!!            Only implemented for non-isochoric interpolation

C ***           Store current solution
                ZF_STORE=ZDF(ns,jdoxf,neax,nhx)

C              IF(.false..and.KTYP1D.NE.2) THEN !central differences
                IF(KTYP1D.NE.2) THEN !central differences
C ***             Perturb solution in one direction
                  ZDF(ns,jdoxf,neax,nhx)=ZF_STORE-HALFDELTA
C ***             Evaluate perturbed residual
                  CALL ZFRF30(IDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,
     '              NORMALSIGN1,nr,nx,CG,PG,RDF1,XDF,YGF,ZDF,
     '              ERROR,*9999)
                ENDIF

C ***           Perturb solution (in the other direction for cent. diffs.).
                ZDF(ns,jdoxf,neax,nhx)=ZDF(ns,jdoxf,neax,nhx)+DELTA

C ***           Evaluate perturbed residual
                CALL ZFRF30(IDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,NORMALSIGN1,
     '            nr,nx,CG,PG,RDF2,XDF,YGF,ZDF,ERROR,*9999)

                IF(DOP) THEN
C KAT 20Mar01: Can't branch out of critical section.
C              Critical section is not essential.
CC$OMP             CRITICAL(ZFFS_1)
                  WRITE(OP_STRING,'(/'' ****** nhx='',I2,'//
     '              ''' jdoxf='',I1,'' ns='',I2,'' ******'')')
     '              nhx,jdoxf,ns
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,
     '              '('' ****** ZDF(ns,jdoxf,neax,nhx)='',D12.5,
     '              '//''' Perturbation='',D12.5)')
     '              ZDF(ns,jdoxf,neax,nhx),DELTA
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
                  WRITE(OP_STRING,'('' ****** Diagonal term: '','//
     '              '''RDF1(ns,1,nhx)='',D20.10,'//
     '              ''' RDF2(ns,1,nhx)='',D20.10)')
     '              RDF1(ns,1,nhx),RDF2(ns,1,nhx)
                  CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$OMP             END CRITICAL(ZFFS_1)
                ENDIF

C ***           Return ZDF to original value
                ZDF(ns,jdoxf,neax,nhx)=ZF_STORE

C ***           Assemble element stiffness matrix
                DO idoxf=1,IDOXFT(nhx)
                  DO ms=1,NSTH
                    FS(ms,idoxf,ns,jdoxf,neax,nhx)=
     '                (RDF2(ms,idoxf,nhx)-RDF1(ms,idoxf,nhx))/DELTA
                  ENDDO
                ENDDO

              ENDDO !ns
            ENDDO !jdoxf
          ENDDO !neax
        ENDDO !nhx

      ELSE !IF(KTYP1D.EQ.1) THEN
C***    Algebraic calculation of tangent stiffness matrix
        CALL ZFFS30(IDOXFT,JDOXFT,NBHF,NBJF,NEAXT,NHF,NIEF,NORMALSIGN1,
     '    nr,nx,CG,FS,PG,XDF,YGF,ZDF,ERROR,*9999)
      ENDIF

      CALL EXITS('ZFFS')
      RETURN
 9999 CALL ERRORS('ZFFS',ERROR)
      CALL EXITS('ZFFS')
      RETURN 1
      END


