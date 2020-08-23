      SUBROUTINE OPTG50(COORDS,IPOINTTYP,STRESSTYPE,NBH,ng,nr,nx,
     '  PHI,PST,RM,TC,TG,TN,TNA,XG,XI,ZG,FULL,ERROR,*)

C#### Subroutine: OPTG50
C###  Description:
C###    OPTG50 outputs stress fields at current solution wrt
C###    'Fibre' or 'Reference' coordinates as specified in COORDS.

C**** Since base vectors of theta coords are not orthonormal, cpts
C**** of Cauchy and Nominal stress are converted to 'physical' values.
C**** XG   are undeformed theta coords and derivs wrt Xi
C**** ZG   are deformed theta coords and derivs wrt undeformed coords
C**** TG   are tensor cpts of 2nd Piola-Kirchhoff stresses
C**** TN   are physical cpts of Nominal stresses
C**** TNA  is the physical (Cauchy) active stress component
C**** TC   are physical cpts of Cauchy stress
C**** PST  are principal stresses
C**** RM   is the modal matrix whose cols=eigenvectors assoc with PST
C**** PHI  are the Euler angles of the principal stresses

      IMPLICIT NONE
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'opti00.cmn'
      INCLUDE 'stra00.cmn'
!     Parameter List
      INTEGER IPOINTTYP,NBH(NHM),ng,nr,nx
      REAL*8 PHI(3),PST(3),
     '  RM(3,3),TC(3,3),TG(3,3),TN(3,3),TNA,XG(NJM,NUM),
     '  XI(3),ZG(NHM,NUM)
      CHARACTER COORDS*(*),ERROR*(*),STRESSTYPE*(*)
      LOGICAL FULL
!     Local Variables
      INTEGER i0,i1,IBEG,IBEG2,IEND,IEND2,IOSTAT,LENI,LENR,LENS,mi,
     '  nb,nhx,ni,njj1,njj2,NITB
      PARAMETER (LENI=4,LENR=12,LENS=1)
      CHARACTER CI*(LENI),CR*(LENI),
     '  FMTI*(LENI+3),FMTR*(LENI+3),POINTLABEL*3

      CALL ENTERS('OPTG50',*9999)
      nb=NBH(NH_LOC(1,nx))
      NITB=NIT(nb)
      IF(TABLE) THEN
        WRITE(CI,'(I2)') LENI
        WRITE(CR,'(I2)') LENR
        FMTI='(I'//CI//')'
C Changing to fixed format
C        FMTR='(*'//CR//')'
        FMTR='(F10.4)'
        FORMAT=' '
        i0=1
        i1=LENI
        WRITE(FORMAT(i0:i1),FMTI) ng
        DO ni=1,NITB
          DO mi=1,ni
            i0=i1+LENS+1
            i1=i0+LENR-1
C            FORMAT(i0:i1)=CFROMR(TC(mi,ni),FMTR)
            WRITE(FORMAT(i0:i1),FMTR) TC(mi,ni)
          ENDDO
        ENDDO
        DO ni=1,NITB
          i0=i1+LENS+1
          i1=i0+LENR-1
C          FORMAT(i0:i1)=CFROMR(PST(ni),FMTR)
          WRITE(FORMAT(i0:i1),FMTR) PST(ni)
        ENDDO
        DO ni=1,NITB
          i0=i1+LENS+1
          i1=i0+LENR-1
C          FORMAT(i0:i1)=CFROMR(PHI(ni),FMTR)
          WRITE(FORMAT(i0:i1),FMTR) PHI(ni)
        ENDDO
        WRITE(OP_STRING,FMT='('''//FORMAT//''')',IOSTAT=IOSTAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(IOSTAT.NE.0) THEN
          ERROR=' Write error'
          GOTO 9999
        ENDIF
      ELSE
C MPN 25May2000: various changes below to tidy up format of output
C MPN 5Aug2014: add format for data output (nd instead of ng)
        IF(IPOINTTYP.EQ.1) THEN
          POINTLABEL='ng='
        ELSE IF(IPOINTTYP.EQ.2) THEN
          POINTLABEL='nd='
        ELSE
          POINTLABEL='pt.'
        ENDIF
        FORMAT='(/'' '//POINTLABEL//''',I9,''  X(j,0):'',6D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) ng,((XG(NJ_LOC(njj1,njj2,nr),1),
     '    njj2=1,NJ_LOC(njj1,0,nr)),njj1=1,2)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( '' ======'', '' Z(j,0):'',6D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) (ZG(NH_LOC(nhx,nx),1),
     '    nhx=1,NH_LOC(0,nx))
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        FORMAT='( ''       '', ''  Xi(i):'',3D12.4)'
        WRITE(OP_STRING,FMT=FORMAT) (XI(ni),ni=1,NITB)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        CALL STRING_TRIM(COORDS,IBEG,IEND)
        CALL STRING_TRIM(STRESSTYPE,IBEG2,IEND2)
        FORMAT='(7X,''2nd Piola-Kirchhoff, Cauchy & Nominal '//
     '    STRESSTYPE(IBEG2:IEND2)//' stresses wrt '//
     '    COORDS(IBEG:IEND)//' coords:'')'
        WRITE(OP_STRING,FMT=FORMAT)
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO mi=1,NITB
          IF(mi.EQ.1) THEN
            FORMAT='(7X,''TG(1,i):'', D12.4,23X,''TC(1,i):'','//
     '        ' D12.4,23X,''TN(1,i):'',3D12.4)'
          ELSE IF(mi.EQ.2) THEN
            FORMAT='(7X,''   2    '',2D12.4,11X,''   2    '','//
     '        '2D12.4,11X,''   2    '',3D12.4)'
          ELSE IF(mi.EQ.3) THEN
            FORMAT='(7X,''   3    '',3D12.4,     ''  3    '','//
     '        '3D12.4,     ''  3    '',3D12.4)'
          ENDIF
          WRITE(OP_STRING,FMT=FORMAT)
     '      (TG(mi,ni),ni=1,mi),(TC(mi,ni),ni=1,mi),
     '      (TN(mi,ni),ni=1,NITB)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO

        IF(FULL) THEN
          CALL STRING_TRIM(COORDS,IBEG,IEND)
          CALL STRING_TRIM(STRESSTYPE,IBEG2,IEND2)
C !!! It would be helpful to say whether these are Cauchy, Nominal, or 2PK.
C         "Principal directions" may be better than "eigenvectors".
          FORMAT='(7X,''Principal '//STRESSTYPE(IBEG2:IEND2)//
     '      ' stresses, Euler angles (deg) and eigenvectors wrt '//
     '      COORDS(IBEG:IEND)//' coords:'')'
          WRITE(OP_STRING,FMT=FORMAT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          DO mi=1,NITB
            IF(mi.EQ.1) THEN
              FORMAT='(7X,''PST(1): '',D12.4,''  PHI(1):'','//
     '          'D12.4,''  RM(1,i):'',3D12.4)'
            ELSE IF(mi.EQ.2) THEN
              FORMAT='(7X,''    2   '',D12.4,''      2  '','//
     '          'D12.4,''     2    '',3D12.4)'
            ELSE IF(mi.EQ.3) THEN
              FORMAT='(7X,''    3   '',D12.4,''      3  '','//
     '          'D12.4,''     3    '',3D12.4)'
            ENDIF
            WRITE(OP_STRING,FMT=FORMAT)
     '        PST(mi),PHI(mi),(RM(mi,ni),ni=1,NITB)
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          ENDDO
        ENDIF

        IF(STRESSTYPE.EQ.'Active') THEN
          FORMAT='( ''       Active Fibre tension (kPa):'',D12.4)'
          WRITE(OP_STRING,FMT=FORMAT) TNA
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDIF

      CALL EXITS('OPTG50')
      RETURN
 9999 CALL ERRORS('OPTG50',ERROR)
      CALL EXITS('OPTG50')
      RETURN 1
      END


