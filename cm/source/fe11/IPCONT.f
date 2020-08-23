      SUBROUTINE IPCONT(ERROR,*)

C#### Subroutine: IPCONT
C###  Description:
C###    IPCONT  inputs parameters used for contact mechanics.
C###    Set up in subroutine IPCONT.
      
      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'nonl00.cmn'

!     Parameter List
      CHARACTER ERROR*(*)
      
!     Local Variables
      INTEGER ICHAR,INFO,NOQUES
      LOGICAL FILEIP
      
      CALL ENTERS('IPCONT',*9999)
      FILEIP=.FALSE.
      NOQUES=0
      ICHAR=999

      IDEFLT(1)=1
      FORMAT='('' Specify frictionless contact method [1]: '''
     &  //'/''   (1) Penalty method'''
     &  //'/''   (2) Cross-constraint method'''
     &  //'/$,''    '',I1)'
      IF(IOTYPE.EQ.3) IDATA(1)=PENALTY
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,2,
     &  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) PENALTY=IDATA(1)

      RDEFLT(1)=1.0d0
      FORMAT='($,'' Enter contact stiffness parameter (frictionless)'
     &  //' [1.0]: '',D12.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=CONT_STIFF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) CONT_STIFF=RDATA(1)

      RDEFLT(1)=1.0d0
      FORMAT='($,'' Enter tied stiffness parameter (tied contact)'
     &  //' [1.0]: '',D12.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=TIED_STIFF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) TIED_STIFF=RDATA(1)

      RDEFLT(1)=1.0d0
      FORMAT='($,'' Enter tangent stiffness parameter (frictional)'
     &  //' [1.0]: '',D12.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=FRIC_STIFF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) FRIC_STIFF=RDATA(1)

      RDEFLT(1)=1.0d0
      FORMAT='($,'' Enter frictional coefficient (frictional)'
     &  //' [1.0]: '',D12.4)'
      IF(IOTYPE.EQ.3) RDATA(1)=FRIC_COEFF
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,0.0d0,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) FRIC_COEFF=RDATA(1)

      RDEFLT(1)=0.0d0
      FORMAT='($,'' Enter gap tolerance'
     &  //' [0.0]: '',D14.6)'
      IF(IOTYPE.EQ.3) RDATA(1)=GAP_TOL
      CALL GINOUT(IOTYPE,IPREAL,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,-RMAX,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) GAP_TOL=RDATA(1)

      IDEFLT(1)=0
      FORMAT='($,'' Enter number of augmentations to perform'
     &  //' [0]: '',I5)'
      IF(IOTYPE.EQ.3) IDATA(1)=AUGMENT
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) AUGMENT=IDATA(1)

      FORMAT='($,'' Do you want to check for all contact points'
     &  //' during Newton steps [N]? '',A)'
      ADEFLT(1)='N'
      IF(IOTYPE.EQ.3) ADATA(1)='N'
      CALL GINOUT(IOTYPE,IPANSW,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &  1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,IMIN,
     &  IMAX,LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,
     &  *9999)
      IF(ADATA(1).EQ.'Y') THEN
        CHECK_CONTACT=.TRUE.
      ELSE
        CHECK_CONTACT=.FALSE.
      ENDIF
C*** 25/09/08 JHC Add user-option to specify when to add the geometric 
C***              term to the contact stiffness matrix
      IDEFLT(1)=0
      FORMAT='($,'' At which Newton step do you want to add'
     &  //' the geometric term to contact stiffness matrix [0]? '',I5)'
      IF(IOTYPE.EQ.3) IDATA(1)=ADD_GEOM
      CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     &  ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     &  LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
      IF(IOTYPE.NE.3) ADD_GEOM=IDATA(1)

      CALL EXITS('IPCONT')
      RETURN
 9999 CALL ERRORS('IPCONT',ERROR)
      CALL EXITS('IPCONT')
      RETURN 1
      END
