      SUBROUTINE ZPZF(IDOXFT,JDOXFT,NBHF,NEAXT,NHF,nr,NSFE,NYNS_E,
     '  nx,SN_E,YP,ZDF,ERROR,*)

C#### Subroutine: ZPZF
C###  Description:
C###    ZPZF transfers global parameters ZP to face parameters ZDF for
C###    interpolation within the face.
C***  Uses arrays set up in FACE_INT_LGE.
C#### Variable: ZDF(ns_f,1:IDOXFT/JDOXFT,neax,nhx)
C###  Type: REAL*8
C###  Set_up: ZPZF
C###  Description:
C###    ZDF(ns_f,idoxf/jdoxf,neax,nhx) are the parameters ns_f for
C###    interpolation over a face using adjacent element neax to obtain
C###    either cross face derivative jdoxf or derivative discontinuity
C###    term jdoxf of dependent variable nhx.  Only neax=1 is used for
C###    derivative discontinuity terms.

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER IDOXFT(NHM),JDOXFT,NBHF(NHM),NEAXT,NHF,nr,
     '  NSFE(NSM,2,2,NHM),NYNS_E(0:NSM,2,2,NHM),nx
      REAL*8 SN_E(NSM,2,2,NHM),YP(NYM,NIYM),ZDF(NSFM,2,2,NHM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ijdoxf,IJDOXFT,jdoxf,nh,nhx,nb_f,neax,ns_e,ns_f,NSTB,ny
      LOGICAL FLUX

      CALL ENTERS('ZPZF',*9999)

C!!!  need to check order of stabilizing derivs
      FLUX=ITYP15(nr,nx).NE.3.OR..TRUE.

      DO nhx=1,NHF
        nh=NH_LOC(nhx,nx)
        nb_f=NBHF(nh)
        NSTB=NST(nb_f)
        IF(FLUX) THEN !Flux Term
          IJDOXFT=JDOXFT !ZDF,SN_E depend on interp func deriv
        ELSE !High order derivative discontinuity term
          IJDOXFT=IDOXFT(nhx) !ZDF,SN_E depends on weight func deriv
        ENDIF
        DO neax=1,NEAXT
          DO ijdoxf=1,IJDOXFT !cross face derivs of interp/weight func
            IF(FLUX) THEN !Flux Term
              jdoxf=ijdoxf !ny depends on interp func deriv
            ELSE !High order derivative discontinuity term
              jdoxf=1 !ny does not depends on weight func deriv
            ENDIF
            DO ns_f=1,NSTB
              ZDF(ns_f,ijdoxf,neax,nhx)=0.0d0
            ENDDO
            DO ns_e=1,NYNS_E(0,jdoxf,neax,nhx)
              ny=NYNS_E(ns_e,jdoxf,neax,nhx)
              ns_f=NSFE(ns_e,jdoxf,neax,nhx) !corresponding component of FS
              IF(ny.NE.0) THEN
                ZDF(ns_f,ijdoxf,neax,nhx)=ZDF(ns_f,ijdoxf,neax,nhx)+
     '            YP(ny,1)*SN_E(ns_e,ijdoxf,neax,nhx)
              ELSE
                ZDF(ns_f,ijdoxf,neax,nhx)=0.0d0
              ENDIF
            ENDDO !ns
          ENDDO !ijdoxf
        ENDDO !neax
        IF(.NOT.FLUX) THEN !High order deriv. discon. term
C         Collect adjacent element terms together
          DO ijdoxf=1,IJDOXFT !cross face derivs of interp/weight func
            jdoxf=1
            DO ns_f=1,NSTB
              ZDF(ns_f,ijdoxf,1,nhx)=
     '          ZDF(ns_f,ijdoxf,1,nhx)+ZDF(ns_f,ijdoxf,2,nhx)
            ENDDO !ns
          ENDDO !ijdoxf
        ENDIF
      ENDDO !nh

      CALL EXITS('ZPZF')
      RETURN
 9999 CALL ERRORS('ZPZF',ERROR)
      CALL EXITS('ZPZF')
      RETURN 1
      END

