      SUBROUTINE OPFIT(NBJ,NEELEM,NKH,NMNO,NPNODE,nr,NVHP,nx,NYNP,
     '  WU,FIX,ERROR,*)

C#### Subroutine: OPFIT
C###  Description:
C###    Outputs fit for region nr.

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'fit000.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'ityp00.cmn'
      INCLUDE 'ktyp00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),NKH(NHM,NPM,NCM),
     '  NMNO(1:2,0:NOPM),NPNODE(0:NP_R_M,0:NRM),nr,
     '  NVHP(NHM,NPM,NCM),nx,
     '  NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 WU(0:NUM+1,NEM)
      CHARACTER ERROR*(*)
      LOGICAL FIX(NYM,NIYFIXM)
!     Local Variables
      INTEGER i,nh,nhj,nhx,nj,njj,nk,nonode,NOTOT,np,nv,ny,NYTOT

      CALL ENTERS('OPFIT',*9999)

      IF(KTYP8.EQ.1.AND.ITYP6(nr,nx).EQ.1) THEN
        WRITE(OP_STRING,'('' Linear geometric fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.1.AND.ITYP6(nr,nx).EQ.2) THEN
        WRITE(OP_STRING,'('' Non-linear geom fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.2.AND.NJ_LOC(njl_fibr,0,nr).EQ.1) THEN
        WRITE(OP_STRING,'('' Fibre angle fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.2.AND.NJ_LOC(njl_fibr,0,nr).EQ.2) THEN
        WRITE(OP_STRING,'('' Fibre (imbrication) angle fitting '','
     '    //'''problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.2.AND.NJ_LOC(njl_fibr,0,nr).EQ.3) THEN
        WRITE(OP_STRING,'('' Sheet angle fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.3) THEN
        WRITE(OP_STRING,'('' Field fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.4) THEN
        WRITE(OP_STRING,'('' Signal fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.5) THEN
        WRITE(OP_STRING,'('' Motion fitting problem with Fourier '','
     '    //'''basis'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.6) THEN
        WRITE(OP_STRING,'('' Data fitting problem by optimisation '')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP8.EQ.7) THEN
        WRITE(OP_STRING,'('' Material parameter fitting problem'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF
      IF(KTYP6.EQ.1) THEN
        WRITE(OP_STRING,'('' Data is specified at Gauss points'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ENDIF

      IF(KTYP12.EQ.0) THEN
        WRITE(OP_STRING,'('' No smoothing constraints defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE IF(KTYP12.EQ.1.OR.KTYP12.EQ.2) THEN
        IF(KTYP12.EQ.1) THEN
          WRITE(OP_STRING,
     '      '('' Sobolev smoothing is applied to the field'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          WRITE(OP_STRING,'('' Sobolev smoothing is applied to '
     '      //'deviation from initial field'')')
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,'('' Smoothing scaling for element 1 is :'','
     '    //'D11.3)') WU(0,NEELEM(1,nr))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP8.EQ.5) THEN
          WRITE(OP_STRING,'('' Smoothing weights for element 1 are:'','
     '      //'9D11.3)') (WU(i,1),i=2,NUM)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.1) THEN
            WRITE(OP_STRING,'('' Smoothing weights for element 1 '
     '        //'are:'',2D11.3)') (WU(i,NEELEM(1,nr)),i=2,3)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.2) THEN
            WRITE(OP_STRING,'('' Smoothing weights for element 1 '
     '        //'are:'',5D11.3)') (WU(i,NEELEM(1,nr)),i=2,6)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(NIT(NBJ(1,NEELEM(1,nr))).EQ.3) THEN
            WRITE(OP_STRING,'('' Smoothing weights for element 1 '
     '        //'are:'',9D11.3)') (WU(I,NEELEM(1,nr)),i=2,11)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF !ktyp8
      ENDIF !ktyp12

      WRITE(OP_STRING,'(/'' Number of fitting problems = '',I1)')
     '  NUM_FIT(0)
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      IF(KTYP8.EQ.7) THEN !material fitting
        WRITE(OP_STRING,'(/'' Number of residual sets  = '',I1)')
     '    KTYP28
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE                !all other fitting
      ENDIF !ktyp8

      NYTOT=0
      NOTOT=0
      DO njj=1,NUM_FIT(0)
        WRITE(OP_STRING,'(/'' Fit Problem: '',I1)') njj
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Number of variables = '',I1)')
     '    NUM_FIT(njj)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(KTYP8.EQ.2) THEN !fibre/sheet fitting
          FORMAT='('' Variable '',I1,'' is stored in '','
     '      //'''fibre variable '',I1,'' and nhx '',I1)'
        ELSE IF(KTYP8.EQ.6) THEN !fitting by optimisation
          FORMAT='('' Variable '',I1,'' is stored in '','
     '      //'''geometric variable '',I1,'' and nhx '',I1)'
        ELSE
          FORMAT='('' Variable '',I1,'' is stored in '','
     '      //'''field variable '',I1,'' and nhx '',I1)'
        ENDIF
        DO nhj=1,NUM_FIT(njj)
          WRITE(OP_STRING,FORMAT) nhj,NLH_FIT(nhj,2,njj),
     '      NLH_FIT(nhj,3,njj)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO !nhj

        IF(KTYP8.EQ.7) THEN !material parameter fitting
          WRITE(OP_STRING,'('' Material params in fit: '','
     '      //'12I3)') (NMNO(1,nhj),nhj=1,NMNO(1,0))
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF !ktyp8=7

        DO nhj=1,NUM_FIT(njj)
          IF(KTYP6.EQ.0) THEN
            WRITE(OP_STRING,'('' Data variable number to be fitted '','
     '        //'''for variable '',I1,'' is '',I1)') nhj,NJ_FIT(nhj,njj)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ELSE IF(KTYP6.EQ.1) THEN
            WRITE(OP_STRING,'('' Gauss variable number to be fitted '','
     '        //'''for variable '',I1,'' is '',I1)') nhj,NG_FIT(nhj,njj)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDDO !nhj

        IF(KTYP8.EQ.4) THEN
          WRITE(OP_STRING,'('' Data region to be fitted '',I1)')
     '      DATA_REGION
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Signal sample skip number '',I5)')
     '      SIGNAL_SKIPSAMPLE
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Input signal filename '',A)')
     '      INFILENAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'('' Output history filename '',A)')
     '      OUTFILENAME
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF

        DO nonode=1,NPNODE(0,nr)
          np=NPNODE(nonode,nr)
          DO nhj=1,NUM_FIT(njj)
            nj=NLH_FIT(nhj,1,njj)
            nhx=NLH_FIT(nhj,3,njj)
            nh=NH_LOC(nhx,nx)
            DO nv=1,NVHP(nh,np,1)
              DO nk=1,NKH(nh,np,1)
                NYTOT=NYTOT+1
                ny=NYNP(nk,nv,nh,np,0,1,nr)
                IF(.NOT.FIX(ny,1)) THEN
                  NOTOT=NOTOT+1
                  WRITE(OP_STRING,'('' ny='',I5,'' : np='',I5,'
     '              //''', nj='',I2,'', nh='',I2,'', nv='',I2,'
     '              //''', nk='',I1,'' - In fit'')') ny,np,nj,nh,nv,nk
                ELSE
                  WRITE(OP_STRING,'('' ny='',I5,'' : np='',I5,'
     '              //''', nj='',I2,'', nh='',I2,'', nv='',I2,'
     '              //''', nk='',I1,'' - Not In fit'')')
     '              ny,np,nj,nh,nv,nk
                ENDIF
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDDO !nk
            ENDDO !nv
          ENDDO !nhj
        ENDDO !nonode (np)
      ENDDO !njj

      WRITE(OP_STRING,'(/,'' Total number of geometric variables '','
     '  //'''- NYT : '',I7)') NYTOT
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      WRITE(OP_STRING,'('' Total number of geometric variables '','
     '   //'''in the fit - NOT : '',I7)') NOTOT
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

      CALL EXITS('OPFIT')
      RETURN
 9999 CALL ERRORS('OPFIT',ERROR)
      CALL EXITS('OPFIT')
      RETURN 1
      END


