      SUBROUTINE ARCDER(NPL,DERIV,DL,XI,XN_LOCAL,ERROR,*)

C#### Subroutine: ARCDER
C###  Description:
C###    ARCDER calculates arc length derivatives at XI coordinates.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'b14.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER NPL(5,0:3)
      REAL*8 DERIV,DL(3),XI,XN_LOCAL(2,3,4)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER I,N,K,nj
      REAL*8 PH3,PL1,PL2,PL3,PL2S1,PL2S3,SL,SM,XA_LOCAL(3,3)

      CALL ENTERS('ARCDER',*9999)
      DO nj=1,3
        DO K=1,3
          XA_LOCAL(K,nj)=0.0D0
        ENDDO
      ENDDO
      DO nj=1,NJT
        IF(NPL(1,nj).EQ.1) THEN
          DO K=1,3
            XA_LOCAL(K,nj)=0.0D0
            DO N=1,2
              XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL1(N,K,XI)*XN_LOCAL(1,nj,N)
            ENDDO
          ENDDO
        ELSE IF(NPL(1,nj).EQ.2) THEN
          DO K=1,3
            XA_LOCAL(K,nj)=0.0D0
            DO N=1,3
              XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL2(N,K,XI)*XN_LOCAL(1,nj,N)
            ENDDO
          ENDDO
        ELSE IF(NPL(1,nj).EQ.3) THEN
          DO K=1,3
            XA_LOCAL(K,nj)=0.0D0
            DO N=1,4
              XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+PL3(N,K,XI)*XN_LOCAL(1,nj,N)
            ENDDO
          ENDDO
        ELSE IF(NPL(1,nj).EQ.4) THEN
          DO K=1,3
            XA_LOCAL(K,nj)=0.0D0
            DO N=1,2
              XA_LOCAL(K,nj)=XA_LOCAL(K,nj)+
     '           PH3(N,1,K,XI)*XN_LOCAL(1,nj,N)
     '          +PH3(N,2,K,XI)*XN_LOCAL(2,nj,N)*DL(N)
              IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$              call mp_setlock()
                WRITE(OP_STRING,'(1X,3E11.4,3I2,3E11.4)')
     '            (DL(I),I=1,3),nj,K,N,
     '            (XN_LOCAL(I,nj,N),I=1,2),XA_LOCAL(K,nj)
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$              call mp_unsetlock()
              ENDIF
            ENDDO
          ENDDO

        ELSE IF(NPL(1,nj).EQ.6) THEN
!         !Added AJP 2-6-93
!         !Special Hermite simplex line,Apex at node 1
          DO K=1,3
            XA_LOCAL(K,nj)=PL2S1(1,1,K,XI)*XN_LOCAL(1,nj,1)
     '              +PL2S1(2,1,K,XI)*XN_LOCAL(1,nj,2)
     '              +PL2S1(2,2,K,XI)*XN_LOCAL(2,nj,2)*DL(2)
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3E11.4)')
     '        (XA_LOCAL(K,nj),K=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ELSE IF(NPL(1,nj).EQ.7) THEN
!         !Added AJP 2-6-93
!         !Special Hermite simplex line, apex at node 3
          DO K=1,3
            XA_LOCAL(K,nj)=PL2S3(1,1,K,XI)*XN_LOCAL(1,nj,1)
     '              +PL2S3(1,2,K,XI)*XN_LOCAL(2,nj,1)*DL(1)
     '              +PL2S3(2,1,K,XI)*XN_LOCAL(1,nj,2)
          ENDDO
          IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$          call mp_setlock()
            WRITE(OP_STRING,'('' XA_LOCAL(k,nj)='',3E11.4)')
     '        (XA_LOCAL(K,nj),K=1,4)
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$          call mp_unsetlock()
          ENDIF
        ENDIF
      ENDDO

      IF(ITYP10(1).EQ.1) THEN
        DERIV=XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2+XA_LOCAL(2,3)**2
      ELSE IF(ITYP10(1).EQ.2) THEN
        DERIV=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1))**2+
     '    XA_LOCAL(2,3)**2
      ELSE IF(ITYP10(1).EQ.3) THEN
        DERIV=XA_LOCAL(2,1)**2+(XA_LOCAL(2,2)*XA_LOCAL(1,1)*
     '      DCOS(XA_LOCAL(1,3)))**2
     '      +(XA_LOCAL(2,3)*XA_LOCAL(1,1))**2
      ELSE IF(ITYP10(1).EQ.4) THEN
        SL=DSINH(XA_LOCAL(1,1))
        SM=DSIN(XA_LOCAL(1,2))
        DERIV=FOCUS*FOCUS*((SL*SL+SM*SM)*
     '                    (XA_LOCAL(2,1)**2+XA_LOCAL(2,2)**2)
     '                   +(SL*SM*XA_LOCAL(2,3))**2)
      ELSE IF(ITYP10(1).EQ.5) THEN
      ENDIF
      DERIV=DSQRT(DERIV)

      CALL EXITS('ARCDER')
      RETURN
 9999 CALL ERRORS('ARCDER',ERROR)
      CALL EXITS('ARCDER')
      RETURN 1
      END


