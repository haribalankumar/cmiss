      SUBROUTINE OPANAL(NPNODE,nr,nx,NYNP,YP,ACTIVATION,ERROR,*)

C#### Subroutine: OPANAL
C###  Description:
C###    OPANAL outputs analytic formula parameters.

      IMPLICIT NONE
      INCLUDE 'anal00.cmn'
      INCLUDE 'b00.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'ityp00.cmn'
!     Parameter List
      INTEGER nr,nx,NPNODE(0:NP_R_M,0:NRM),
     '   NYNP(NKM,NVM,NHM,NPM,0:NRCM,NCM,NRM)
      REAL*8 YP(NYM,NIYM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER nonode,np,ny1,ny2,ny3,ny4
      CHARACTER ANAL2DLAPLACE(11)*56,ANAL3DLAPLACE(15)*35,
     '  ANAL2DPOISSON(3)*22,ANAL3DPOISSON(7)*21

CC AJPs
      INTEGER nj,n_site
      LOGICAL ACTIVATION
CC AJPe

      DATA ANAL2DLAPLACE
     '  /'K(x-y)                                                  ',
     '   'K(x^2-y^2)                                              ',
     '   'K(x^2+2xy-y^2)                                          ',
     '   'K(a.r^n.cos(n.t)+b.r^m.sin(m.t)+c.r.cos(t)+d.r.sin(t)+e)',
     '   'Centre dipole (single circle)                           ',
     '   'Centre dipole (mulitple circles)                        ',
     '   'Eccentric dipole (single circle)                        ',
     '   'Eccentric dipole (mulitple circles)                     ',
     '   'Anisotropic annulus                                     ',
     '   'Anisotropic plate, f(z)=k.e^z                           ',
     '   'Anisotropic plate, f(z)=a.sin(z)+b.cos(z)               '/
      DATA ANAL3DLAPLACE
     '  /'K(x-y)                             ',
     '   'K(z)                               ',
     '   'K(x^2-y^2)                         ',
     '   'K(x^2-z^2)                         ',
     '   'K(y^2-z^2)                         ',
     '   'K(x^2+y^2-2z^2)                    ',
     '   'K(x^2-2y^2+z^2)                    ',
     '   'K(-2x^2+y^2+z^2)                   ',
     '   'K(x^2+2xy-y^2)                     ',
     '   'K(x^2+2xz-z^2)                     ',
     '   'K(y^2+2yz-z^2)                     ',
     '   'Centre dipole (single sphere)      ',
     '   'Centre dipole (mulitple spheres)   ',
     '   'Eccentric dipole (single sphere)   ',
     '   'Eccentric dipole (mulitple spheres)'/

      DATA ANAL2DPOISSON
     '  /'A*r+N*log(r)+C        ',
     '   'Bidomain transfer test',
     '   'Double circle problem '/

      DATA ANAL3DPOISSON
     '  /'Case 1b              ',
     '   'Case 3a              ',
     '   'Case 2b              ',
     '   'Case 1a              ',
     '   'Case 2a              ',
     '   'Case 3b              ',
     '   'Double sphere problem'/

      CALL ENTERS('OPANAL',*9999)

      WRITE(OP_STRING,'(/'' Region '',I2,'' (nx='',I1,'')'')') nr,nx
      CALL WRITES(IOFI,OP_STRING,ERROR,*9999)

CC AJPs
      IF(ACTIVATION) THEN
        WRITE(OP_STRING,'(/'' Analytic activation sequence defined'')')
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        WRITE(OP_STRING,'(/'' Number of activation sites is '',I2)')
     '    N_SITES
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        DO n_site=1,N_SITES
          WRITE(OP_STRING,'('' Coordinates of activation site '',I2,'
     '    //''' are '',3D12.4)')n_site,(ACTVN_SITE(nj,n_site),nj=1,NJT)
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ENDDO
        WRITE(OP_STRING,'(/'' Activation time interval is  '',D12.4)')
     '    ACTVN_TIME
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        IF(ANAL_CHOICE(nr).EQ.1) THEN
          WRITE(OP_STRING,'('' Activation times calc at nodes'')')
        ELSEIF(ANAL_CHOICE(nr).EQ.2) THEN
          WRITE(OP_STRING,'('' Activation times calc at grid '
     '      //'points'')')
        ENDIF
        CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
      ELSE
        IF(ITYP5(nr,nx).EQ.1) THEN !static analysis
          IF(ITYP2(nr,nx).EQ.3) THEN !Laplaces equation
            WRITE(OP_STRING,'(/'' Equation is Laplace''''s equation'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NJT.EQ.2) THEN
              WRITE(OP_STRING,'(/'' Analytic solution is '',A)')
     '          ANAL2DLAPLACE(ANAL_CHOICE(nr))
            ELSE IF(NJT.EQ.3) THEN
              WRITE(OP_STRING,'(/'' Analytic solution is '',A)')
     '          ANAL3DLAPLACE(ANAL_CHOICE(nr))
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(NJT.EQ.2) THEN
              IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.3) THEN
                WRITE(OP_STRING(1),'('' Value of K = '',D12.4)') ANAL_K
              ELSE IF(ANAL_CHOICE(nr).EQ.4) THEN
                WRITE(OP_STRING(1),'('' Value of K = '',D12.4)') ANAL_K
                WRITE(OP_STRING(2),'('' Value of n = '',I12)')   ANAL_N
                WRITE(OP_STRING(3),'('' Value of m = '',I12)')   ANAL_M
                WRITE(OP_STRING(4),'('' Value of a = '',D12.4)') ANAL_A
                WRITE(OP_STRING(5),'('' Value of b = '',D12.4)') ANAL_B
                WRITE(OP_STRING(6),'('' Value of c = '',D12.4)') ANAL_C
                WRITE(OP_STRING(7),'('' Value of d = '',D12.4)') ANAL_D
                WRITE(OP_STRING(8),'('' Value of e = '',D12.4)') ANAL_E
              ELSE IF(ANAL_CHOICE(nr).EQ.9) THEN
                WRITE(OP_STRING(1),'('' Value of K     = '',D12.4)')
     '            ANISO_K
                WRITE(OP_STRING(2),'('' Value of m     = '',I12)')
     '            ANISO_M
              ELSE IF(ANAL_CHOICE(nr).EQ.10) THEN
                WRITE(OP_STRING(1),'('' Value of K     = '',D12.4)')
     '            ANISO_K
                WRITE(OP_STRING(3),'('' Value of L     = '',D12.4)')
     '            ANISO_L
                WRITE(OP_STRING(4),'('' Value of H     = '',D12.4)')
     '            ANISO_H
                WRITE(OP_STRING(5),'('' Value of theta = '',D12.4,'
     '            //''' deg'')') ANISO_THETA*180.0d0/PI
              ELSE IF(ANAL_CHOICE(nr).EQ.11) THEN
                WRITE(OP_STRING(1),'('' Value of a     = '',D12.4)')
     '            ANAL_A
                WRITE(OP_STRING(2),'('' Value of b     = '',D12.4)')
     '            ANAL_B
                WRITE(OP_STRING(4),'('' Value of L     = '',D12.4)')
     '            ANISO_L
                WRITE(OP_STRING(5),'('' Value of H     = '',D12.4)')
     '            ANISO_H
                WRITE(OP_STRING(6),'('' Value of theta = '',D12.4,'
     '            //''' deg'')') ANISO_THETA*180.0d0/PI
              ENDIF
            ELSE IF(NJT.EQ.3) THEN
              IF(ANAL_CHOICE(nr).GE.1.AND.ANAL_CHOICE(nr).LE.11) THEN
                WRITE(OP_STRING(1),'('' Value of K = '',D12.4)') ANAL_K
              ENDIF
            ENDIF
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
C CPB 29/11/00 Adding Poissons equation analytic solutions
          ELSE IF(ITYP2(nr,nx).EQ.3) THEN !Poissons equation
            WRITE(OP_STRING,'(/'' Equation is Poisson''''s equation'')')
            CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            IF(ITYP3(nr,nx).EQ.2) THEN !special rhs term
              WRITE(OP_STRING,'(/'' Equation has a special rhs term'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(NJT.EQ.2) THEN
                WRITE(OP_STRING,'(/'' Analytic solution is '',A)')
     '            ANAL2DPOISSON(ANAL_CHOICE(nr))
              ELSE IF(NJT.EQ.3) THEN
                WRITE(OP_STRING,'(/'' Analytic solution is '',A)')
     '            ANAL3DPOISSON(ANAL_CHOICE(nr))
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              IF(NJT.EQ.2) THEN
                IF(ANAL_CHOICE(nr).EQ.1) THEN
                  WRITE(OP_STRING(1),'('' Value of A = '',D12.4)')
     '              ANAL_A
                  WRITE(OP_STRING(2),'('' Value of B = '',D12.4)')
     '              ANAL_B
                  WRITE(OP_STRING(3),'('' Value of C = '',D12.4)')
     '              ANAL_C
                ELSE IF(ANAL_CHOICE(nr).EQ.3) THEN
                  WRITE(OP_STRING(1),
     '              '('' Value of C_1       = '',D12.4)') ANAL_C_1
                  WRITE(OP_STRING(3),
     '              '('' Value of D_1       = '',D12.4)') ANAL_D_1
                  WRITE(OP_STRING(3),
     '              '('' Value of H_1       = '',D12.4)') ANAL_H_1
                  WRITE(OP_STRING(3),
     '              '('' Value of reference = '',D12.4)') ANAL_REF
                ENDIF
              ELSE IF(NJT.EQ.3) THEN
                IF(ANAL_CHOICE(nr).EQ.7) THEN
                  WRITE(OP_STRING(1),
     '              '('' Value of C_1       = '',D12.4)') ANAL_C_1
                  WRITE(OP_STRING(2),
     '              '('' Value of C_11      = '',D12.4)') ANAL_C_11
                  WRITE(OP_STRING(3),
     '              '('' Value of D_11      = '',D12.4)') ANAL_D_11
                  WRITE(OP_STRING(3),
     '              '('' Value of H_11      = '',D12.4)') ANAL_H_11
                  WRITE(OP_STRING(3),
     '              '('' Value of reference = '',D12.4)') ANAL_REF
                ENDIF
              ENDIF
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
            ELSE
              ERROR='>>Equation type not implemented'
              GOTO 9999
            ENDIF
          ELSE
            ERROR='>>Equation type not implemented'
            GOTO 9999
          ENDIF
        ELSE IF(ITYP5(nr,nx).EQ.2) THEN !time integration
          IF(ITYP2(nr,nx).EQ.3) THEN !advection-diffusion
C DMAL 22-MAY-2002
            IF(ANAL_CHOICE(nr).EQ.1) THEN
              WRITE(OP_STRING(1),
     '          '(''The equation is the diffusion equation'//
     '          ' on a rectangular plate (width=A, height=B) with'//
     '          ' the temperature on three sides fixed at zero and'//
     '          ' the other having a half sine wave profile.'')')
              WRITE(OP_STRING(2),
     '          '(''The differnce between the numerical and analytic'//
     '          ' solution is saved in YP(iy=5).'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING(1),
     '          '(''Width, A of plate (Xi=1) = '//
     '          ''',D12.4)') ANAL_A
              WRITE(OP_STRING(2),
     '          '(''Height, B of plate (Xi=2) = '//
     '          ''',D12.4)') ANAL_B
              WRITE(OP_STRING(3),
     '          '(''Diffusion constant, D = '//
     '          ''',D12.4)') ANAL_D
              WRITE(OP_STRING(4),
     '          '(''Analytical solution calculated at time'//
     '          ' = '',D12.4)') ANAL_TIME
              ANAL_E=0.0d0
C calculating rms error at nodes and saving difference to YP(iy=5)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny1=NYNP(1,1,1,np,0,1,nr)
                ny2=NYNP(2,1,1,np,0,1,nr)
                ny3=NYNP(3,1,1,np,0,1,nr)
                ny4=NYNP(4,1,1,np,0,1,nr)
                ANAL_E=ANAL_E+(YP(ny1,1)-YP(ny1,7))**2
                YP(ny1,5)=YP(ny1,1)-YP(ny1,7)
                YP(ny2,5)=YP(ny2,1)-YP(ny2,7)
                YP(ny3,5)=YP(ny3,1)-YP(ny3,7)
                YP(ny4,5)=YP(ny4,1)-YP(ny4,7)
              ENDDO
              ANAL_E=SQRT(ANAL_E/NPNODE(0,nr))
              WRITE(OP_STRING(6),
     '          '(''RMS error'//
     '          ' = '',D12.4)') ANAL_E
            ELSE IF(ANAL_CHOICE(nr).EQ.2) THEN
              WRITE(OP_STRING(1),
     '          '(''The equation is the diffusion equation'//
     '          ' on a rod (Length=L) with'//
     '          ' T(0,t)=T1 and T(L,t)=T2'')')
              WRITE(OP_STRING(2),
     '          '(''The differnce between the numerical and analytic'//
     '          ' solution is saved in YP(iy=5).'')')
              CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
              WRITE(OP_STRING(1),
     '          '(''Length,L of rod (Xi=1) = '//
     '          ''',D12.4)') ANAL_A
              WRITE(OP_STRING(2),
     '          '(''Temperature at x=0 = '//
     '          ''',D12.4)') ANAL_B
              WRITE(OP_STRING(3),
     '          '(''Temperature at x=L = '//
     '          ''',D12.4)') ANAL_C
              WRITE(OP_STRING(4),
     '          '(''Diffusion constant, D = '//
     '          ''',D12.4)') ANAL_D
              WRITE(OP_STRING(4),
     '          '(''Analytical solution calculated at time'//
     '          ' = '',D12.4)') ANAL_TIME
              ANAL_E=0.0d0
C calculating rms error at nodes and saving difference to YP(iy=5)
              DO nonode=1,NPNODE(0,nr)
                np=NPNODE(nonode,nr)
                ny1=NYNP(1,1,1,np,0,1,nr)
                ny2=NYNP(2,1,1,np,0,1,nr)
                ANAL_E=ANAL_E+(YP(ny1,1)-YP(ny1,7))**2
                YP(ny1,5)=YP(ny1,1)-YP(ny1,7)
                YP(ny2,5)=YP(ny2,1)-YP(ny2,7)
              ENDDO
              ANAL_E=SQRT(ANAL_E/NPNODE(0,nr))
              WRITE(OP_STRING(6),
     '          '(''RMS error'//
     '          ' = '',D12.4)') ANAL_E
            ELSE
              ERROR='>>Analysis type not implemented'
              GOTO 9999
            ENDIF
          ELSE
            ERROR='>>Analysis type not implemented'
            GOTO 9999
          ENDIF
          CALL WRITES(IOFI,OP_STRING,ERROR,*9999)
        ELSE
          ERROR='>>Analysis type not implemented'
          GOTO 9999
        ENDIF
      ENDIF !activation

      CALL EXITS('OPANAL')
      RETURN
 9999 CALL ERRORS('OPANAL',ERROR)
      CALL EXITS('OPANAL')
      RETURN 1
      END
CC AJPe

