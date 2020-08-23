      SUBROUTINE CAP_NE(NBJ,NELIST,NENP,NPNE,nr,CE,XP,ERROR,*)

C#### Subroutine: CAP_NE
C###  Description:
C###    CAP_NE inputs parameter values for capillary elements
C###    into CE and defines element arrays via GN1DNEJ.

C*** Created by KSB, 7th March, 2002.

C... dimensions defined are CE(nm_C0,ne): the length of the membrane
C... portion of capillary & CE(nm_a0,ne): length of tissue portion of
C... capillary at zero Ptp and Ptm. From this CE(nm_a,ne) & CE(nm_b,ne),
C... the elliptical approximations, are calculated.

C*** all dimensions should be in mm

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'pulm00.cmn'
!     Parameter list
      INTEGER NBJ(NJM,NEM),NELIST(0:NEM),
     '  NENP(NPM,0:NEPM,0:NRM),NPNE(NNM,NBFM,NEM),nr
      REAL*8 CE(NMM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local variables
      INTEGER nb,ISEED,ne,noelem,np1,np2,nj
      REAL*8 CM_RANDOM_NUMBER,DIA_CAPS(2),DIAMETER,LEN_CAP,length,
     '  perimeter

      CALL ENTERS('CAP_NE',*9999)

      nb=NBJ(1,NELIST(1))
      DIA_CAPS(1)=0.012d0 !(mm) from Huang (2001)
      DIA_CAPS(2)=0.0165d0 !(mm) from Huang (2001)      
      length=0.d0
      ISEED=1
      DIAMETER=0.d0
      DO noelem=1,NELIST(0)
        ne=NELIST(noelem)
C...  defines vessel dimensions
        CE(nm_C0,ne)=CM_RANDOM_NUMBER(ISEED)*(DIA_CAPS(2)-DIA_CAPS(1))+
     &    DIA_CAPS(1)
        CE(nm_a0,ne)=0.006d0
        perimeter=CE(nm_C0,ne)+CE(nm_a0,ne)
        CE(nm_b,ne)=CE(nm_a0,ne)/2.d0
        CE(nm_a,ne)=DSQRT(DABS(perimeter**2.d0/(2.d0*PI**2.d0)-
     '    CE(nm_b,ne)**2.d0))
        CE(nm_Dh,ne)=4.d0*PI*CE(nm_a,ne)*CE(nm_b,ne)/perimeter !(mm) 4Ac/perimeter
        diameter=diameter+CE(nm_Dh,ne) !calculating average hydraulic diameter of mesh
        LEN_CAP=0.d0 !initialise
        np1=NPNE(1,nb,ne) !length of capillary segment (LEN_CAP) =
        np2=NPNE(2,nb,ne) !distance between the 2 nodes.
C... If node not connected flow solution can't be calculated
        IF(ne.NE.CAP_INLET.AND.ne.NE.CAP_OUTLET) THEN
          IF(NENP(np1,0,nr).LE.1) THEN
            WRITE(OP_STRING,'('' WARNING:node not connected '//
     '        '(increase # boundary points):'',I5)') np1
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ELSE IF(NENP(np2,0,nr).LE.1) THEN
            WRITE(OP_STRING,'('' WARNING:node not connected '//
     '        '(increase # boundary points):'',I5)') np2
            CALL WRITES(IODI,OP_STRING,ERROR,*9999)
          ENDIF
        ENDIF
        DO nj=1,NJT
          LEN_CAP=LEN_CAP+(XP(1,1,nj,np2)-XP(1,1,nj,np1))**2.d0
        ENDDO !nj
        CE(nm_length,ne)=DSQRT(LEN_CAP)*scale_factor  !mm:
        IF(ne.NE.CAP_INLET.AND.ne.NE.CAP_OUTLET) THEN
          length=length+CE(nm_length,ne)
        ENDIF
      ENDDO !noelem
      length=length/NELIST(0) !avge length of capillary segments
      WRITE(OP_STRING,'('' Average length of capillary segments in '//
     '  'whole mesh (mm):'',D12.4)') length
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)
      diameter=diameter/NELIST(0) !avge diameter of capillary segments
      WRITE(OP_STRING,'('' Average diameter of capillary segments in '//
     '  'whole mesh (mm):'',D12.4)') diameter
      CALL WRITES(IODI,OP_STRING,ERROR,*9999)

      CALL EXITS('CAP_NE')
      RETURN
 9999 CALL ERRORS('CAP_NE',ERROR)
      CALL EXITS('CAP_NE')
      RETURN 1
      END

