      SUBROUTINE BASIS8(IBTYP,NAN,nb,NGAP,PG,WG,XIG,ERROR,*)

C#### Subroutine: BASIS8
C###  Description:
C###    BASIS8 supplements basis function array
C###    PG(ns,nu,ng,nb),ns=1,.NST(nb), with auxiliary functions
C###    ns=NST(nb) + na, where na=1,..NAT(nb), for use with element
C###    parameters ZA(na,nc,nh,ne).
C MPN 7/4/93 - Separating out NST(nb) and NAT(nb)
!old C**** NST(nb) is then replaced by NST(nb)+NAT(nb).

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
!     Parameter List
      INTEGER IBTYP,NAN(NIM,NAM),nb,NGAP(NIM,NBM)
      REAL*8 PG(NSM,NUM,NGM),WG(NGM),XIG(NIM,NGM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER ICHAR,INFO,na,ni,NOQUES,IEND
      CHARACTER CHAR1*100,CHAR2*1
CC JPP 18Mar2003 New string... note the mess with CHAR2 below...      
      CHARACTER CHAR3*12
      LOGICAL FILEIP

      CALL ENTERS('BASIS8',*9999)
      FILEIP=.FALSE.

      NOQUES=0
      INFO=2
      ICHAR=999

      IF(IBTYP.EQ.0) THEN
        IF(nb.GT.1) THEN
          IDEFLT(1)=NIT(nb-1)
        ELSE
          IDEFLT(1)=1
        ENDIF
        WRITE(CHAR1,'(I1)') IDEFLT(1)
        FORMAT='($,'' Enter the number of Xi-coordinates ['//
     '    CHAR1(1:1)//']: '',I1)'
        IF(IOTYPE.EQ.3) IDATA(1)=NIT(nb)
        CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,3,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
        IF(iotype.ne.3) NIT(nb)=IDATA(1)

        NUT(nb)=NIT(nb)*NIT(nb)+2
        CALL ASSERT(NUT(nb).LE.NUM,'>>Increase NUM',ERROR,*9999)
        NGT(nb)=1
        DO ni=1,NIT(nb)
          NGAP(ni,nb)=1
        ENDDO
        DO ni=1,NIT(nb)
          IF(nb.GT.1) THEN
            IF(NIT(nb).EQ.NIT(nb-1)) THEN
              IDEFLT(1)=NGAP(ni,nb-1)
            ENDIF
          ELSE
            IDEFLT(1)=2
          ENDIF
          WRITE(CHAR1,'(I1)') ni
          WRITE(CHAR2,'(I1)') IDEFLT(1)
          FORMAT='($,'' Enter the number of Gauss points in the Xi('//
     '      CHAR1(1:1)//') direction ['//CHAR2(1:1)//']: '',I1)'
          IF(IOTYPE.EQ.3) IDATA(1)=NGAP(ni,nb)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,1,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) THEN
            NGAP(ni,nb)=IDATA(1)
            NGT(nb)=NGT(nb)*NGAP(ni,nb)
            CALL ASSERT(NGT(nb).LE.NGM,'>>Need to increase NGM',
     '        ERROR,*9999)
          ENDIF
        ENDDO
      ENDIF
      DO na=1,NAT(nb)
CC JPP 18Mar2003 Generalizing to variable-length integer
        IEND=0
        CALL APPENDI(IEND,na,CHAR3)
        FORMAT='('' For auxiliary basis function parameter '
     '    //CHAR3(1:IEND)//': '')'
C JPP 18Mar2003 These were the original statements...
C        WRITE(CHAR2,'(I1)') na
C        FORMAT='('' For auxiliary basis function parameter '
C     '    //CHAR2(1:1)//': '')'
        CALL GINOUT(IOTYPE,IPMESS,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,1,
     '    ADATA,ANO,CDATA,CDEFLT,ICHAR,IDATA,IDEFLT,0,IMAX,
     '    LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)

        DO ni=1,NIT(nb)
          WRITE(CHAR1,'(I1)') ni
          IF(NABTYP(nb).EQ.1) THEN      !Legendre auxiliary basis
            FORMAT='($,'' The polynomial degree in the Xi('
     '        //CHAR1(1:1)//') direction is [0]: '',I2)'
          ELSE IF(NABTYP(nb).EQ.2) THEN !Fourier  auxiliary basis
            FORMAT='($,'' Enter n, where 2n is the wavenumber in'
     '        //' the Xi('//CHAR1(1:1)//') direction [0]: '',I2)'
          ELSE IF(NABTYP(nb).EQ.3) THEN !Pressure auxiliary basis
            FORMAT='($,'' The polynomial degree in the Xi('
     '        //CHAR1(1:1)//') direction is [0]: '',I2)'
          ENDIF
          IF(IOTYPE.EQ.3) IDATA(1)=NAN(ni,na)
          CALL GINOUT(IOTYPE,IPINTE,IVDU,IFILE,0,0,NOQUES,FILEIP,FORMAT,
     &      1,ADATA,ADEFLT,CDATA,CDEFLT,ICHAR,IDATA,IZERO,-2,9,
     '      LDATA,LDEFLT,RDATA,RDEFLT,RMIN,RMAX,INFO,ERROR,*9999)
          IF(iotype.ne.3) NAN(ni,na)=IDATA(1)
        ENDDO
      ENDDO
C MPN 7/4/93 - Separating out NST(nb) and NAT(nb)
c      NST(nb)=NST(nb)+NAT(nb)
      CALL GAUSS8(NAN,nb,NGAP(1,nb),PG,WG,XIG,ERROR,*9999)

      CALL EXITS('BASIS8')
      RETURN
 9999 CALL ERRORS('BASIS8',ERROR)
      CALL EXITS('BASIS8')
      RETURN 1
      END


