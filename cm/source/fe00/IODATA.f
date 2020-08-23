      SUBROUTINE IODATA(COMAND,ANGLE_TYPE,TYPE,IUNIT,NDP,NJTT,
     '  TITLE,WD,Z_CONT,ZD,ERROR,*)

C#### Subroutine: IODATA
C###  Description:
C###    <HTML>
C###    IODATA handles I/O of measured field data for optimisation
C###    and plotting.
C###    <BR>
C###    COMAND determines whether data is to be read 'READ' from or
C###    written 'WRITE' to IUNIT.
C###    </HTML>

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'defn00.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ioda00.cmn'
!     Parameter List
      INTEGER IUNIT,NDP(NDM),NJTT
      REAL*8 WD(NJM,NDM),Z_CONT(NDM,2,67),ZD(NJM,NDM)
      CHARACTER ANGLE_TYPE*(*),COMAND*(*),ERROR*(*),TITLE*80,TYPE*(*)

C*** LKC 2-JUL-1999 Note: These parameters (partially) determine the
C format for the output of data. A new smaller format
C LENR_WD has been added to  reduce the output for the weights.
!     Parameters
      INTEGER LENI,LENR,LENR_CONT,LENR_WD,LENS
C GRC 25-JUN-2001 LENI is set to 5 for backward compatibility
C     == minimum space used for writing integers
      PARAMETER (LENI=5,LENR=10,LENR_CONT=25,LENR_WD=6,LENS=1)

!     Local Variables
      INTEGER i0,i1,IBEG,IEND,IOSTAT,nd,NDP_LEN,NDP_MIN,NDP_MAX,nj,
     '  WARN_NDP
      CHARACTER ACCESS*10,CHAR*(LENR),
     '  CI*(LENI),FMTI*(LENI+3),FMTR*(LENR+3),
     '  FMTR_WD*(LENR+3),FORM*11,FORMAT*110



      LOGICAL OPENED

C cpb 7/4/94 this needs to be generalised for nj_loc
      CALL ENTERS('IODATA',*9999)

      IF(DOP) THEN
C KAT 14May01: Can't branch out of critical section.
C              Critical section is not essential.
CC$      call mp_setlock()
        WRITE(OP_STRING,'(A)') ' >IODATA'
        CALL WRITES(IODI,OP_STRING,ERROR,*9999)
CC$      call mp_unsetlock()
      ENDIF

      INQUIRE(UNIT=IUNIT,OPENED=OPENED)
      IF(OPENED) THEN
        INQUIRE(UNIT=IUNIT,FORM=FORM)
        IF(FORM.EQ.'FORMATTED') THEN
          INQUIRE(UNIT=IUNIT,ACCESS=ACCESS)
          IF(ACCESS.EQ.'SEQUENTIAL') THEN
            IF(COMAND.EQ.'READ') THEN
              READ(UNIT=IUNIT,FMT='(A80)') TITLE(1:80)
              IF(ADD) THEN
                NDTOLD=NDT
              ELSE
                NDTOLD=0
              ENDIF
              NDT=NDTOLD

              DO nd=NDTOLD+1,NDM
C GMH 14/2/97 Notify CMGUI
                CALL DATA_CHANGE(nd,.TRUE.,ERROR,*9999)
C CPB 8/4/94 This needs to be generalised for NJ_LOC
                IF(TYPE(1:6).EQ.'FIBRES'.OR.TYPE(1:6).EQ.'SHEETS') THEN
                  READ(UNIT=IUNIT,FMT=*,IOSTAT=IOSTAT,END=9998)
     '              NDP(nd),
     '              (ZD(nj,nd),NJ=1,NJT),ZD(NJTT,nd),
     '              (WD(nj,nd),NJ=1,NJT),WD(NJTT,nd)
C*** 10/10/08 JHC Reading in contact pressure data
                ELSE IF(TYPE(1:7).EQ.'CONTACT') THEN
                  READ(UNIT=IUNIT,FMT=*,IOSTAT=IOSTAT,END=9998)
     &              NDP(nd),(Z_CONT(nd,1,15+nj*2),nj=1,NJTT),
     &              (WD(15+nj*2,nd),nj=1,NJTT)
                  Z_CONT(nd,2,17)=Z_CONT(nd,1,17)
                  Z_CONT(nd,2,19)=Z_CONT(nd,1,19)
                  Z_CONT(nd,2,21)=Z_CONT(nd,1,21)
                ELSE
                  READ(UNIT=IUNIT,FMT=*,IOSTAT=IOSTAT,END=9998)
     '              NDP(nd),(ZD(nj,nd),nj=1,NJTT),(WD(nj,nd),nj=1,NJTT)
                ENDIF
                IF(ANGLE_TYPE(1:7).EQ.'DEGREES') THEN
                  !degrees entered in file - storage is always radians
                  ZD(NJTT,nd)=ZD(NJTT,nd)*PI/180.0d0
                ENDIF
                IF(IOSTAT.GT.0) THEN
                  ERROR=' Read error'
                  GOTO 9999
                ENDIF
                NDT=NDT+1
              ENDDO
              ERROR=' Increase NDM as #data pts in file >NDM'
              GO TO 9999
            ELSE IF(COMAND.EQ.'WRITE') THEN
c             WRITE(UNIT=IUNIT,FMT='(A80)') TITLE(1:80)
              WRITE(UNIT=IUNIT,FMT='('' Data file'')')

C***  LKC 2-JUL-1999 Formats for data and weights
C*** 10/10/08 JHC more decimal pts for contact pressure data
              IF(TYPE(1:7).EQ.'CONTACT') THEN
                FMTR='(G25.17)'
              ELSE
                FMTR='(F10.4)'
              ENDIF
              FMTR_WD='(F6.2)'

              WRITE(CHAR,'(I2)') (NJTT*(LENR+LENS)-3)
              FORMAT='('' NDP: ZD:'','//CHAR//'X,''WD:'')'
C              FORMAT='('' NDP: ZD:'','//
C     '          CFROMI((NJTT*(LENR+LENS)-3),'(I2)')//'X,''WD:'')'
c              WRITE(UNIT=IUNIT,FMT=FORMAT)

C LKC 13-MAY-1999 if NDP is not setup
C              IF(NDP(1).LE.0.AND.NDP(2).LE.0) THEN
C                WRITE(OP_STRING,
C     '            '(/'' >>WARNING: NDP not setup, doing it now'')')
C                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C                DO nd=1,NDT
C                  NDP(nd)=nd
C                ENDDO
C              ENDIF

C GRC 25-JUN-2001 Get largest text length required for any NDP,
C     starting with LENI=5 for backward compatibility.
              NDP_LEN = LENI
C     NDP_MIN and NDP_MAX store the minimum and maximum values that
C     can be written with the current NDP_LEN, assumed to start at 5.
              NDP_MIN = -9999
              NDP_MAX = 99999

              DO nd=1,NDT
                DO WHILE ((NDP(nd).GT.NDP_MAX).OR.
     '            (NDP(nd).LT.NDP_MIN))
                  NDP_MIN = NDP_MIN*10 - 9
                  NDP_MAX = NDP_MAX*10 + 9
                  NDP_LEN = NDP_LEN + 1
                ENDDO
              ENDDO

              WRITE(CI,'(I2)') NDP_LEN
              FMTI='(I'//CI//')'

C LKC 9-OCT-2002 intialise NDP if it has not already been
C     setup.
C             DO nd=1,NDT
              WARN_NDP=0
              DO nd=1,NDT
                IF(NDP(nd).LE.0.OR.NDP(nd).GT.100000000) THEN
                  WARN_NDP=WARN_NDP+1
                  NDP(nd)=nd
                ENDIF

                FORMAT=' '
                i0=1
                i1=NDP_LEN
                WRITE(FORMAT(i0:i1),FMTI) NDP(nd)
C                FORMAT(i0:i1)=CFROMI(NDP(nd),FMTI)
C*** 10/10/08 JHC Added contact pressure data points
                IF(TYPE(1:7).EQ.'CONTACT') THEN 
                  DO nj=1,NJTT
                    i0=i1+LENS+1
                    i1=i0+LENR_CONT-1
                    WRITE(FORMAT(i0:i1),FMTR) Z_CONT(nd,1,15+2*nj)
                  ENDDO
                ELSE 
                  DO nj=1,NJTT
                    i0=i1+LENS+1
                    i1=i0+LENR-1
                    WRITE(FORMAT(i0:i1),FMTR) ZD(nj,nd)
C                    FORMAT(i0:i1)=CFROMR(ZD(nj,nd),FMTR)
                  ENDDO
                ENDIF
                DO nj=1,NJTT
                  i0=i1+LENS+1
                  i1=i0+LENR_WD-1
                  WRITE(FORMAT(i0:i1),FMTR_WD) WD(nj,nd)
C                  FORMAT(i0:i1)=CFROMR(WD(nj,nd),FMTR)
                ENDDO
                WRITE(UNIT=IUNIT,FMT='('''//FORMAT//''')',
     '            IOSTAT=IOSTAT)
                IF(IOSTAT.NE.0) THEN
                  ERROR=' Write error'
                  GOTO 9999
                ENDIF
              ENDDO

C LKC 9-OCT-2002 new warning
              IF(WARN_NDP.NE.0) THEN
                WRITE(OP_STRING,'('' >> Warning: NDP not setup, '','
     '            //'''altering '',I7,'' of '',I7'
     '            //''' reference data num.'')') WARN_NDP,NDT
                CALL WRITES(IODI,OP_STRING,ERROR,*9999)
              ENDIF

            ELSE
              ERROR=' Command error: COMAND='//COMAND
              GOTO 9999
            ENDIF
          ELSE
            WRITE(CI,'(I3)') IUNIT
C            CI=CFROMI(IUNIT,'(I3)')
            CALL STRING_TRIM(CI,IBEG,IEND)
            ERROR=' File '//
     '        CI(IBEG:IEND)//' is not connected for sequential access'
            GOTO 9999
          ENDIF
        ELSE
          WRITE(CI,'(I3)') IUNIT
C          CI=CFROMI(IUNIT,'(I3)')
          CALL STRING_TRIM(CI,IBEG,IEND)
          ERROR=' File '//
     '      CI(IBEG:IEND)//' is not connected for formatted i/o'
          GOTO 9999
        ENDIF
      ELSE
        WRITE(CI,'(I3)') IUNIT
C        CI=CFROMI(IUNIT,'(I3)')
        CALL STRING_TRIM(CI,IBEG,IEND)
        ERROR=' File '//CI(IBEG:IEND)//' is not opened'
        GOTO 9999
      ENDIF

 9998 ERROR=' '
      CALL EXITS('IODATA')
      RETURN
 9999 CALL ERRORS('IODATA',ERROR)
      CALL EXITS('IODATA')
      RETURN 1
      END

