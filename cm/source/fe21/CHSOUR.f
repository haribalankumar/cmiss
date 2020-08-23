      SUBROUTINE CHSOUR(DIPOLE_DIR_NTIME,DIPOLE_DIR,DIPOLE_CEN,
     '   NDIPOLES,STRING,ERROR,*)

C#### Subroutine: CHSOUR
C###  Description:
C###    CHSOUR allows user to change the source locations.

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'back00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'file00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'grid00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'sour00.cmn'

!     Parameter List
      CHARACTER ERROR*(*),STRING*(MXCH)
      INTEGER DIPOLE_DIR_NTIME(NDIPOLEM,NRM),
     '    NDIPOLES(NRM),nr


      REAL*8 DIPOLE_CEN(4,0:NDIPTIMM,NDIPOLEM,NRM),
     '     DIPOLE_DIR(4,0:NDIPTIMM,NDIPOLEM,NRM)
      
!     Local Variables
      INTEGER IBEG,IEND,NTRL,dip,npts,N3CO,IFROMC
      LOGICAL CBBREV,ABBREV,TRANSL,SCALE
      REAL*8 TXYZ(3),RL1(3)

      CALL ENTERS('CHSOUR',*9999)
      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM change source
C###  Parameter:     <translate by DX#,DY#,DZ#>
C###    Translates the source cnetre points in the cartesian coordinate system  C###  Parameter:      <scale by SX#,SY#,SZ#>
C###    Scale the dipole vector in cartesian
C###    coordinates
C###  Parameter:     <class #[1]>
c###    Specify the class number  
C###  Description:
C###    Change the source locations

        OP_STRING(1)=BLANK(1:15)//'<translate by DX#,DY#,DZ#>'
        OP_STRING(2)=BLANK(1:15)//'<scale by SX#,SY#,SZ#>'
        OP_STRING(3)=BLANK(1:15)//'<class #[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','CHSOUR',ERROR,*9999)
      ELSE
        IF(CBBREV(CO,'TRANSLATE',1,noco+1,NTCO,N3CO)) THEN
           IF (ABBREV(CO(N3CO+1),'BY',1)) THEN

              CALL PARSRL(CO(N3CO+2),3,NTRL,TXYZ,ERROR,*9999)
              TRANSL=.TRUE.
             
           ENDIF
        ELSE
           TRANSL=.FALSE.
        ENDIF

        IF(CBBREV(CO,'SCALE',1,noco+1,NTCO,N3CO)) THEN
           IF(CBBREV(CO,'BY',1,noco+2,noco+3,N3CO)) THEN
            CALL PARSRL(CO(N3CO+1),3,NTRL,RL1,ERROR,*9999)
            SCALE=.TRUE.
           ENDIF
        ELSE
            SCALE=.FALSE.
        ENDIF



        IF(CBBREV(CO,'CLASS',3,noco+1,NTCO,N3CO)) THEN
           nr=IFROMC(CO(N3CO+1))
        ELSE
           nr=1
        ENDIF

        IF(TRANSL) THEN
           DO dip=1,NDIPOLES(nr)
              DO npts=0,DIPOLE_DIR_NTIME(dip,nr)
                 DIPOLE_CEN(1,npts,dip,nr)=TXYZ(1)+
     '                DIPOLE_CEN(1,npts,dip,nr)
                 DIPOLE_CEN(2,npts,dip,nr)=TXYZ(2)+
     '             DIPOLE_CEN(2,npts,dip,nr)
                 DIPOLE_CEN(3,npts,dip,nr)=TXYZ(3)+
     '             DIPOLE_CEN(3,npts,dip,nr)
              ENDDO
           ENDDO

        ENDIF


        IF(SCALE) THEN
           DO dip=1,NDIPOLES(nr)
              DO npts=0,DIPOLE_DIR_NTIME(dip,nr)
                 DIPOLE_DIR(1,npts,dip,nr)=RL1(1)*
     '                DIPOLE_DIR(1,npts,dip,nr)
                 DIPOLE_DIR(2,npts,dip,nr)=RL1(2)*
     '             DIPOLE_DIR(2,npts,dip,nr)
                 DIPOLE_DIR(3,npts,dip,nr)=RL1(3)*
     '             DIPOLE_DIR(3,npts,dip,nr)
              ENDDO
           ENDDO

        ENDIF

      ENDIF

      CALL EXITS('CHSOUR')
      RETURN
 9999 CALL ERRORS('CHSOUR',ERROR)
      CALL EXITS('CHSOUR')
      RETURN 1
      END


