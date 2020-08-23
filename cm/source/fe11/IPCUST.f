      SUBROUTINE IPCUST(ERROR,*)

C#### Subroutine: IPCUST
C###  Description:
C###    IPCUST defines parameters for customising a mesh.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'chmesh0.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,n,ndat,NOQUES,nores
      LOGICAL FILEIP
      CHARACTER CHAR1*2,CHAR3*3,CHAR4*4

      CALL ENTERS('IPCUST',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      FORMAT='('' Enter type of customisation: '''//
     '      '/''   (1) Use Polynomial and Cosine functions '''//
     '      '/''   (2) Use 2-D measurements for simple scaling'''//
     '      '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=CUSTOMISATION_TYPE-4
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CUSTOMISATION_TYPE=IDATA(1)+4

      IF(CUSTOMISATION_TYPE.EQ.5) THEN
        FORMAT='($,'' Enter the number of polynomial coeffs '
     '    //'[4]: '',I2)'
        IDEFLT(1)=4
        IF(IOTYPE.EQ.3) IDATA(1)=NPC
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NPC=IDATA(1)

        DO n=1,NPC
          WRITE(CHAR3,'(I3)') n
          FORMAT='($,'' Enter coeff value '//CHAR3//' [1.0]:'',D12.4)'
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) RDATA(1)=POLY_COEFFs(n)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) POLY_COEFFS(n)=RDATA(1)
        ENDDO !n
        FORMAT='($,'' Enter the number of cosine coeffs [2]: '',I2)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=NCC
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NCC=IDATA(1)

        DO n=1,NCC
          WRITE(CHAR3,'(I3)') n
          FORMAT='($,'' Enter coeff value '//CHAR3//' [1.0]:'',D12.4)'
          IF(IOTYPE.EQ.3) RDATA(1)=COS_COEFFS(n)
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) COS_COEFFS(n)=RDATA(1)
        ENDDO !n

      ELSEIF(CUSTOMISATION_TYPE.EQ.6) THEN
        FORMAT='('' Enter the number of 2-d measurements '
     '    //'[3]: '',I2)'
        IF(IOTYPE.EQ.3) THEN
          IDATA(1)=NTDM
        ENDIF
        IDEFLT(1)=3
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NMEASUREMX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) THEN
          NTDM=IDATA(1)
        ENDIF
        DO nores=1,NTDM
          WRITE(CHAR1,'(I2)') nores
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=WDMEASURE(1,nores)
          ENDIF
          FORMAT='($,'' Enter the 2d measurement in x direction '
     '      //CHAR1(1:2)//' [1.0]: '',D12.4))'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WDMEASURE(1,nores)=RDATA(1)
          ENDIF
          RDEFLT(1)=1.0d0
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=WDMEASURE(2,nores)
          ENDIF
          FORMAT='($,'' Enter the 2d measurement in y direction '
     '      //CHAR1(1:2)//' [1.0]: '',D12.4))'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WDMEASURE(2,nores)=RDATA(1)
          ENDIF
          RDEFLT(1)=DBLE(nores)/DBLE(NTDM)
          WRITE(CHAR4,'(F4.1)') RDEFLT(1)
          IF(IOTYPE.EQ.3) THEN
            RDATA(1)=WDMEASURE(3,nores)
          ENDIF
          FORMAT='($,'' Enter the height of the measurements '
     '      //CHAR1(1:2)//' ['//CHAR4//']: '',D12.4))'
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,1,0,NOQUES,FILEIP,FORMAT,
     '      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            WDMEASURE(3,nores)=RDATA(1)
          ENDIF
        ENDDO
       ENDIF

      FORMAT='('' Enter type of length scaling: '''//
     '      '/''   (1) Simple scaling '''//
     '      '/''   (2) Variable length scaling'''//
     '      '/$,''    '',I1)'
      IDEFLT(1)=1
      IF(IOTYPE.EQ.3) IDATA(1)=SCLTYPE
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) SCLTYPE=IDATA(1)


      IF(SCLTYPE.EQ.1) THEN
        FORMAT='($,'' Enter the total length '
     '    //'[1.0]: '',D12.4)'
        RDEFLT(1)=1.0d0
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) RDATA(1)=TORSO_LENGTHS(0,1)
        CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) TORSO_LENGTHS(0,1)=RDATA(1)

      ELSEIF(SCLTYPE.EQ.2) THEN
        FORMAT='($,'' Enter the number of length measurements '
     '    //'[2]: '',I2)'
        IDEFLT(1)=2
        IF(IOTYPE.EQ.3) IDATA(1)=NTL
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,2,99,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(IOTYPE.NE.3) NTL=IDATA(1)
        DO n=1,NTL
          WRITE(CHAR3,'(I3)') n
          FORMAT='($,'' Enter height on generic model and corresponding'
     '      //' actual height [1.0,1.0]: '',D12.4,D12.4)'
          IF(IOTYPE.EQ.3) THEN
            DO ndat=1,2
              RDATA(ndat)=TORSO_LENGTHS(ndat,n)
            ENDDO !ndat
          ENDIF
          CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      2,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,NPM,
     '      LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
          IF(IOTYPE.NE.3) THEN
            DO ndat=1,2
             TORSO_LENGTHS(ndat,n)=RDATA(ndat)
            ENDDO !ndat
          ENDIF
        ENDDO !n
      ENDIF
      CALL EXITS('IPCUST')
      RETURN
 9999 CALL ERRORS('IPCUST',ERROR)
      CALL EXITS('IPCUST')
      RETURN 1
      END


