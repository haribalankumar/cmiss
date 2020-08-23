      SUBROUTINE STRAIN_ENERGY(nr,ne,ng,EG,gse,ERROR,*)
C#### Subroutine: STRAIN_ENERGY
C###  Description:
C###    STRAIN_ENERGY evaluates the strain energy at a gauss point for 3D finite
C###      elasticity problems with grid points at gauss points and the constitutive
C###      law defined using cellml.
C###    Inputs:
C###      EG - contains Green strain components, it is assumed to be symmetrical.
C###      Other parameters are as set up for the problem.
C###    Outputs:
C###      gse - the strain energy for this gauss/grid point.

! Maxima code to do the same given green strain in E, dW/dE in S(E)
! and the deformation gradient in terms of k in Fc
! StrainEnergy(x):=block(
!   [Fev,k,E,se,G],
!   k:x,
!   Fev:ev(Fc),
!   E:(1/2)*(transpose(Fev).Fev-ident(dim)),
!   se:0,
!
!   G:zeromatrix(dim,dim),
!   for i:1 thru dim do for j:1 thru dim do (
!       G[i,j]:z,
!       rb:ev(romberg( S(G)[i,j],z,0,E[i,j])),
!       G[i,j]:E[i,j],
!       se:se+rb
!   ),
!   se
! );

      IMPLICIT NONE

      CHARACTER*(*) ROUTINENAME
      PARAMETER(ROUTINENAME='STRAIN_ENERGY')

      INCLUDE 'error0.inc'
      INCLUDE 'sten50.cmn'

!     Parameter List
      INTEGER ne,ng,nr
      REAL*8 EG(3,3),
     &  gse ! gauss point strain energy density

      INTEGER WORKARRLEN
      PARAMETER(WORKARRLEN=200)

      ! locals
      INTEGER
     &  neval,ier,limit,lenw,last,iwork(WORKARRLEN),i,j
      REAL*8
     &  work(4*WORKARRLEN),epsabs,epsrel,energy,abserr
      CHARACTER ERROR*(*)

      EXTERNAL EVAL_S11
      EXTERNAL EVAL_S22
      EXTERNAL EVAL_S33
      EXTERNAL EVAL_S12
      EXTERNAL EVAL_S13
      EXTERNAL EVAL_S23

      CALL ENTERS(ROUTINENAME,*9999)

!       PRINT *,ROUTINENAME, ' enter EG ',((EG(i,j),i=1,3),j=1,3)

      ! initialise variables
      gse=0.0d0
      energy=0.0d0

      ! copy parameters to common variables so they are available to EVAL_S*
      cnr=nr
      cne=ne
      cng=ng
      cerror=ERROR
      cerr_ptr=%LOC(ERROR)
      DO i=1,3
        DO j=1,3
          cEG(i,j)=0.0d0
        ENDDO
      ENDDO

      epsabs=1.0d-10
      epsrel=1.0d-10
      limit=WORKARRLEN
      lenw=4*limit

      CALL dqags(EVAL_S11,0.0d0,EG(1,1),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF
      cEG(1,1)=EG(1,1)

      CALL dqags(EVAL_S22,0.0d0,EG(2,2),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF
      cEG(2,2)=EG(2,2)

      CALL dqags(EVAL_S33,0.0d0,EG(3,3),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF
      cEG(3,3)=EG(3,3)

      CALL dqags(EVAL_S12,0.0d0,EG(1,2),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+2.0d0*energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF
      cEG(1,2)=EG(1,2)
      cEG(2,1)=EG(1,2)

      CALL dqags(EVAL_S13,0.0d0,EG(1,3),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+2.0d0*energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF
      cEG(1,3)=EG(1,3)
      cEG(3,1)=EG(1,3)

      CALL dqags(EVAL_S23,0.0d0,EG(2,3),
     &  epsabs,epsrel,energy,
     &  abserr,neval,ier,limit,lenw,last,iwork,work)
      gse=gse+2.0d0*energy
      IF(ier.NE.0) THEN
        WRITE(ERROR,'(''>>dqags failed ( '',I1,'')'')') ier
        GOTO 9999
      ENDIF

!       PRINT *,ROUTINENAME,' i!=j cEG=',
!      &((cEG(i,j),i=1,3),j=1,3)
!       PRINT *,ROUTINENAME, ' gse tot=', gse

      CALL EXITS(ROUTINENAME)
      RETURN
 9999 CALL ERRORS(ROUTINENAME,ERROR)
      CALL EXITS(ROUTINENAME)
      RETURN 1
      END

      REAL*8 FUNCTION EVAL_S11(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,1,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S11=dWdE
 9999   END

      REAL*8 FUNCTION EVAL_S22(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,2,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S22=dWdE
 9999   END

      REAL*8 FUNCTION EVAL_S33(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,3,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S33=dWdE
 9999   END

      REAL*8 FUNCTION EVAL_S12(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,4,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S12=dWdE
 9999   END

      REAL*8 FUNCTION EVAL_S13(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,5,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S13=dWdE
 9999   END

      REAL*8 FUNCTION EVAL_S23(x)
        IMPLICIT NONE
        INCLUDE 'mxch.inc'
        INCLUDE 'ptr00.cmn'
        INCLUDE 'sten50.cmn'
        REAL*8 x,dWdE

        CALL EVALDW(x,cEG,6,cnr,cne,cng,%VAL(NQNE_PTR),dWdE,
     &    cerror,*9999)
        EVAL_S23=dWdE
 9999   END
