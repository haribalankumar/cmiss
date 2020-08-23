      SUBROUTINE EVERROR(NEELEM,NELIST,NKJE,NPF,NPNE,NRLIST,
     '  NVJE,NW,CG,NEERR,PG,SE,WG,XA,XE,XG,XP,
     '  YG,STRING,ERROR,*)

C#### Subroutine: EVERROR
C###  Description:
C###     EVERROR evaluates the strain energy in an element from fields
C###     of the stress components (fields 1,2 and 3), strain energy
C###     error norm, and the relative energy norm error.
C**** Created by Carey Stevens 23 June 1997

      IMPLICIT NONE
      INCLUDE 'mxch.inc'
      INCLUDE 'b01.cmn'
      INCLUDE 'b13.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'cbdi10.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh00.cmn'
!     Parameter List
      INTEGER NEELEM(0:NE_R_M,0:NRM),NELIST(0:NEM),
     '  NKJE(NKM,NNM,NJM,NEM),NPF(9,NFM),NPNE(NNM,NBFM,NEM),
     '  NRLIST(0:NRM),NVJE(NNM,NBFM,NJM,NEM),NW(NEM,3,NXM)
      REAL*8 CG(NMM),NEERR(NEM,3),PG(NSM,NUM,NGM,NBM),SE(NSM,NBFM,NEM),
     '  WG(NGM,NBM),XA(NAM,NJM,NEM),XE(NSM,NJM),XG(NJM,NUM),
     '  XP(NKM,NVM,NJM,NPM),YG(NIYGM,NGM,NEM)
      CHARACTER ERROR*(*),STRING*(MXCH)
!     Local Variables
      INTEGER IBEG,IEND,IFROMC,N3CO,nb,NBJ_temp(12),ne,ng,nj,njj1,njj2,
     '  no_nelist,no_nrlist,nr,nx
      REAL*8 EG_dif(3,3),EG_disp(3,3),EG_field(3,3),ENERGY_ELEMENT_dif,
     '  ENERGY_ELEMENT_disp,ENERGY_ELEMENT_field,ENERGY_GP_field,
     '  ENERGY_GP_dif,ENERGY_GP_disp,ENERGY_NORM_ELEMENT,
     '  ERR_NORM_ELEMENT,POIS,
     '  RELATIVE_ERR,TG_dif(3,3),TG_disp(3,3),TG_field(3,3),YMOD
      LOGICAL ALL_REGIONS,CBBREV,OPFILE
      CHARACTER FILE*(MXCH)

      CALL ENTERS('EVERROR',*9999)

      IF(CO(noco+1).EQ.'?') THEN
        CALL STRING_TRIM(STRING,IBEG,IEND)

C---------------------------------------------------------------------

C#### Command: FEM evaluate error<;FILENAME>
C###  Parameter:        <element (#s/all)[all]>
C###    Define what elements are to be included.
C###  Parameter:        <basis #[1]>
C###  Parameter:        <region (#s/all)[1]>
C###    Specify the region numbers to evaluate.
C###  Description:
C###    Evaluates the strain energy in an element from fields
C###    of the stress components (fields 1,2 and 3), strain energy
C###    error norm, and the relative energy norm error.
C###    Writes to screen or file FILENAME.operr if qualifier present.

        OP_STRING(1)=STRING(1:IEND)//'<;FILENAME>'
        OP_STRING(2)=BLANK(1:15)//'<element (#s/all)[all]>'
        OP_STRING(3)=BLANK(1:15)//'<basis #[1]>'
        OP_STRING(4)=BLANK(1:15)//'<region (#s/all)[1]>'
        CALL WRITES(IOH1,OP_STRING,ERROR,*9999)

C---------------------------------------------------------------------

      ELSE IF(CO(noco+1).EQ.'??') THEN
        CALL DOCUM('fe21','doc','EVERROR',ERROR,*9999)
      ELSE
        IF(NTCOQU(noco).GT.0) THEN
          FILE=COQU(noco,1)
          CALL STRING_TRIM(FILE,IBEG,IEND)
          OPFILE=.TRUE.
          IOFI=IOFILE1
          CALL OPENF(IOFI,'DISK',FILE(IBEG:IEND)//'.operr','NEW',
     '      'SEQUEN','FORMATTED',132,ERROR,*9999)
        ELSE
          OPFILE=.FALSE.
          IOFI=IOOP
        ENDIF

        nx=1 ! MPN 14Jun2000  may need generalising?

        CALL PARSE_REGIONS(NRLIST,noco,NTCO,CO,ALL_REGIONS,ERROR,*9999)
        CALL PARSE_ELEMENTS(NEELEM,NELIST,noco,NRLIST,NTCO,CO,
     '    ERROR,*9999)

        IF(CBBREV(CO,'BASIS',2,noco+1,NTCO,N3CO)) THEN
           nb=IFROMC(CO(N3CO+1))
        ELSE
           nb=1
        ENDIF

        YMOD=CG(1)  !is Young's modulus
        POIS=CG(2)  !is Poisson's ratio
        ENERGY_NORM_TOTAL=0.0d0
        ERR_NORM_TOTAL=0.0d0
        DO no_nrlist=1,NRLIST(0)
          nr=NRLIST(no_nrlist)

          DO njj1=1,3   !geom/fibres/field
            DO njj2=1,NJ_LOC(njj1,0,nr)
              nj=NJ_LOC(njj1,njj2,nr)
              CALL ASSERT(nj.LE.12,'>>ERROR: increase size '
     '          //'of NBJ_temp to NJM',ERROR,*9999)
                NBJ_temp(nj)=nb
            ENDDO !njj2
          ENDDO !njj1

          DO no_nelist=1,NELIST(0)
            ne=NELIST(no_nelist)

            CALL ASSERT((NW(ne,1,nx).EQ.11).OR.(NW(ne,1,nx).EQ.12),
     '        '>>Only plane stress case implemented',ERROR,*9999)
            CALL ASSERT(IMT(NW(ne,1,nx)).EQ.1,
     '        '>>Only isotropic case implemented',ERROR,*9999)

            WRITE(OP_STRING,
     '        '(/'' Element '',I3,'' :'')') ne
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

            ENERGY_NORM_ELEMENT=0.0d0
            ERR_NORM_ELEMENT=0.0d0
            ENERGY_ELEMENT_field=0.0d0
            ENERGY_ELEMENT_disp=0.0d0
            ENERGY_ELEMENT_dif=0.0d0
            CALL XPXE(NBJ_temp,NKJE(1,1,1,ne),NPF(1,1),
     '        NPNE(1,1,ne),nr,NVJE(1,1,1,ne),
     '        SE(1,1,ne),XA(1,1,ne),XE,XP,ERROR,*9999)
            DO ng=1,NGT(nb)
              CALL XEXG(NBJ_temp,ng,nr,PG,XE,XG,ERROR,*9999)

              TG_field(1,1)=XG(NJ_LOC(NJL_FIEL,1,nr),1)
              TG_field(2,2)=XG(NJ_LOC(NJL_FIEL,2,nr),1)
              TG_field(1,2)=XG(NJ_LOC(NJL_FIEL,3,nr),1)

              IF(NW(ne,1,nx).EQ.11) THEN !plane stress
                EG_field(1,1)=(TG_field(1,1)/YMOD)-
     '            (POIS/YMOD)*TG_field(2,2)
                EG_field(2,2)=(TG_field(2,2)/YMOD)-
     '            (POIS/YMOD)*TG_field(1,1)
                EG_field(1,2)=TG_field(1,2)*(2.0d0*(POIS+1.0d0)/YMOD)
              ELSE IF(NW(ne,1,nx).EQ.12) THEN !plane strain
                EG_field(1,1)=(TG_field(1,1)*((1-POIS**2)/YMOD))-
     '            TG_field(2,2)*((POIS*(POIS+1))/YMOD)
                EG_field(2,2)=(TG_field(2,2)*((1-POIS**2)/YMOD))-
     '            TG_field(1,1)*((POIS*(POIS+1))/YMOD)
                EG_field(1,2)=TG_field(1,2)*((2.0d0*(POIS+1.0d0))/YMOD)
              ENDIF

              TG_disp(1,1)=YG(1,ng,ne)
              TG_disp(2,2)=YG(2,ng,ne)
              TG_disp(1,2)=YG(3,ng,ne)

              IF(NW(ne,1,nx).EQ.11) THEN !plane stress
                EG_disp(1,1)=(TG_disp(1,1)/YMOD)-
     '            (POIS/YMOD)*TG_disp(2,2)
                EG_disp(2,2)=(TG_disp(2,2)/YMOD)-
     '            (POIS/YMOD)*TG_disp(1,1)
                EG_disp(1,2)=TG_disp(1,2)*(2.0d0*(POIS+1.0d0)/YMOD)
              ELSE IF(NW(ne,1,nx).EQ.12) THEN !plane strain
                EG_disp(1,1)=(TG_disp(1,1)*((1-POIS**2)/YMOD))-
     '            TG_disp(2,2)*((POIS*(POIS+1))/YMOD)
                EG_disp(2,2)=(TG_disp(2,2)*((1-POIS**2)/YMOD))-
     '            TG_disp(1,1)*((POIS*(POIS+1))/YMOD)
                EG_disp(1,2)=TG_disp(1,2)*((2.0d0*(POIS+1.0d0))/YMOD)
              ENDIF

              TG_dif(1,1)=TG_disp(1,1)-TG_field(1,1)
              TG_dif(2,2)=TG_disp(2,2)-TG_field(2,2)
              TG_dif(1,2)=TG_disp(1,2)-TG_field(1,2)

              EG_dif(1,1)=EG_disp(1,1)-EG_field(1,1)
              EG_dif(2,2)=EG_disp(2,2)-EG_field(2,2)
              EG_dif(1,2)=EG_disp(1,2)-EG_field(1,2)

              ENERGY_GP_field=0.5d0*(TG_field(1,1)*EG_field(1,1)+
     '          TG_field(2,2)*EG_field(2,2)+
     '          TG_field(1,2)*EG_field(1,2))
              ENERGY_GP_disp=0.5d0*(TG_disp(1,1)*EG_disp(1,1)+
     '          TG_disp(2,2)*EG_disp(2,2)+
     '          TG_disp(1,2)*EG_disp(1,2))
              ENERGY_GP_dif=0.5d0*(TG_dif(1,1)*EG_dif(1,1)+
     '          TG_dif(2,2)*EG_dif(2,2)+TG_dif(1,2)*EG_dif(1,2))

              ENERGY_ELEMENT_field=ENERGY_ELEMENT_field+
     '          ENERGY_GP_field*WG(ng,nb)
              ENERGY_ELEMENT_disp=ENERGY_ELEMENT_disp+
     '          ENERGY_GP_disp*WG(ng,nb)
              ENERGY_ELEMENT_dif=ENERGY_ELEMENT_dif+
     '          ENERGY_GP_dif*WG(ng,nb)

              IF(DOP) THEN
                WRITE(OP_STRING,'(/'' Gauss point '',I2)') ng
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' EG_field(1,1)='',D12.4,'' EG_field(2,2)='',
     '            D12.4,'' EG_field(1,2)='',D12.4)') EG_field(1,1),
     '            EG_field(2,2),EG_field(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' TG_field(1,1)='',D12.4,'' TG_field(2,2)='',
     '            D12.4,'' TG_field(1,2)='',D12.4,'' kN/m'')')
     '            TG_field(1,1),TG_field(2,2),TG_field(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

                WRITE(OP_STRING,
     '            '(/'' EG_disp(1,1)='',D12.4,'' EG_disp(2,2)='',D12.4,'
     '            //''' EG_disp(1,2)='',D12.4)')
     '            EG_disp(1,1),EG_disp(2,2),EG_disp(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' TG_disp(1,1)='',D12.4,'' TG_disp(2,2)='',D12.4,'
     '            //''' TG_disp(1,2)='',D12.4,'' kN/m'')')
     '            TG_disp(1,1),TG_disp(2,2),TG_disp(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

                WRITE(OP_STRING,
     '            '(/'' EG_dif(1,1)='',D12.4,'' EG_dif(2,2)='',D12.4,'
     '            //''' EG_dif(1,2)='',D12.4)')
     '            EG_dif(1,1),EG_dif(2,2),EG_dif(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' TG_dif(1,1)='',D12.4,'' TG_dif(2,2)='',D12.4,'
     '            //''' TG_dif(1,2)='',D12.4,'' kN/m'')')
     '            TG_dif(1,1),TG_dif(2,2),TG_dif(1,2)
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

                WRITE(OP_STRING,
     '            '('' Strain energy at gauss pt '
     '             //'(field) = '',D12.4,'' kJ/m^2'')')
     '            ENERGY_GP_field
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' Strain energy at gauss pt '
     '            //'(disp.) = '',D12.4,'' kJ/m^2'')')
     '            ENERGY_GP_disp
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
                WRITE(OP_STRING,
     '            '('' Strain energy at gauss pt '
     '            //'(diff.) = '',D12.4,'' kJ/m^2'')')
     '            ENERGY_GP_dif
                CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              ENDIF
            ENDDO !ng

            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' Total S.E. for element (field.) '
     '          //'= '',D12.4,'' kJ/m^2'')')
     '          ENERGY_ELEMENT_field
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

              WRITE(OP_STRING,
     '          '('' Total S.E. for element (disp.) '
     '          //'= '',D12.4,'' kJ/m^2'')')
     '          ENERGY_ELEMENT_disp
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

              WRITE(OP_STRING,
     '          '('' Total S.E. for element (diff.) '
     '          //'= '',D12.4,'' kJ/m^2'')')
     '          ENERGY_ELEMENT_dif
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF

            ERR_NORM_ELEMENT=DSQRT(DABS(ENERGY_ELEMENT_dif))

            ENERGY_NORM_ELEMENT=DSQRT(DABS(ENERGY_ELEMENT_field))


            RELATIVE_ERR=ERR_NORM_ELEMENT/ENERGY_NORM_ELEMENT*100d0
            IF(DOP) THEN
              WRITE(OP_STRING,
     '          '('' S.E. norm for element = '',D12.4,'' kJ/m^2'')')
     '          ENERGY_NORM_ELEMENT
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING,
     '          '('' S.E. error norm for element = '',D12.4,
     '          '' kJ/m^2'')') ERR_NORM_ELEMENT
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENDIF
            WRITE(OP_STRING,
     '        '('' Relative energy norm error = '',F8.2,'' %'')')
     '        RELATIVE_ERR
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ENERGY_NORM_TOTAL=ENERGY_NORM_TOTAL+ENERGY_NORM_ELEMENT**2
            ERR_NORM_TOTAL=ERR_NORM_TOTAL+ERR_NORM_ELEMENT**2

            NEERR(ne,1)=RELATIVE_ERR !Percentage
            NEERR(ne,2)=ERR_NORM_ELEMENT

          ENDDO !no_nelist
        ENDDO !no_nrlist

        ERR_NORM_TOTAL=DSQRT(ERR_NORM_TOTAL)
        ENERGY_NORM_TOTAL=DSQRT(ENERGY_NORM_TOTAL)
        RELATIVE_ERR=ERR_NORM_TOTAL/ENERGY_NORM_TOTAL*100.0d0
        IF(DOP) THEN
          WRITE(OP_STRING,'(/'' S.E. norm for listed elements '
     '      //'= '',D12.4,'' kJ/m^2'')') ENERGY_NORM_TOTAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
          WRITE(OP_STRING,'(/'' S.E. error norm for listed elements '
     '      //'= '',D12.4,'' kJ/m^2'')') ERR_NORM_TOTAL
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDIF
        WRITE(OP_STRING,
     '    '('' Relative energy norm error for listed elements '
     '    //'= '',F8.2,'' %'')')
     '    RELATIVE_ERR
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

        IF(OPFILE) THEN
          CALL CLOSEF(IOFI,ERROR,*9999)
          IOFI=IOOP
        ENDIF
      ENDIF
      CALL EXITS('EVERROR')
      RETURN
 9999 CALL ERRORS('EVERROR',ERROR)
      CALL EXITS('EVERROR')
      RETURN 1
      END


