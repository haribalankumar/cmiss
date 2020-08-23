      SUBROUTINE INTEGRATE_LINE(CHXE,INTEGRAL,ERROR,*)

C#### Subroutine: INTEGRATE_LINE
C###  Description:
C###    
C**** Written by Duane Malcolm, 06 November 2003

      IMPLICIT NONE
      INCLUDE 'geom00.cmn'

!     Parameter List
      REAL*8 INTEGRAL,CHXE(NSM,4)
      CHARACTER ERROR*(*)

!     Local Variables
      INTEGER ng
      REAL*8 DPHI(4),DXDXI,DYDXI,DZDXI,FVALUE,GP(3),JACOBIAN,
     &  PHI(4),WEIGHTS(3)
C      LOGICAL 

      DATA DPHI/0.0D0,
     &          0.0D0,
     &          0.0D0,
     &          0.0D0/
      DATA PHI /0.0D0,
     &          0.0D0,
     &          0.0D0,
     &          0.0D0/
      DATA GP/0.112701665D0,
     &        0.5D0,
     &        0.887298335D0/
      DATA WEIGHTS /0.277777778D0,
     &              0.444444444D0,
     &              0.277777778D0/

      CALL ENTERS('INTEGRATE_LINE',*9999)
      
      INTEGRAL=0.0D0
      DO ng=1,3
        PHI(1)=1.0D0-3.0D0*GP(ng)**2.0D0+2.0D0*GP(ng)**3.0D0
        PHI(2)=GP(ng)*(GP(ng)-1.0D0)**2.0D0
        PHI(3)=GP(ng)**2.0D0*(3.0D0-2.0D0*GP(ng))
        PHI(4)=GP(ng)**2.0D0*(GP(ng)-1.0D0)
        
        DPHI(1)=-6.0D0*GP(ng)+6.0D0*GP(ng)**2
        DPHI(2)=3.0D0*GP(ng)**2-4.0D0*GP(ng)+1
        DPHI(3)=6.0D0*GP(ng)-6.0D0*GP(ng)**2
        DPHI(4)=3.0D0*GP(ng)**2-2.0D0*GP(ng)
        
        DXDXI=DPHI(1)*CHXE(1,1)+DPHI(2)*CHXE(2,1)+
     &        DPHI(3)*CHXE(3,1)+DPHI(4)*CHXE(4,1)
        DYDXI=DPHI(1)*CHXE(1,2)+DPHI(2)*CHXE(2,2)+
     &        DPHI(3)*CHXE(3,2)+DPHI(4)*CHXE(4,2)
        DZDXI=DPHI(1)*CHXE(1,3)+DPHI(2)*CHXE(2,3)+
     &        DPHI(3)*CHXE(3,3)+DPHI(4)*CHXE(4,3)
        JACOBIAN=SQRT(DXDXI**2+DYDXI**2+DZDXI**2)
        
        FVALUE=PHI(1)*CHXE(1,4)+PHI(2)*CHXE(2,4)+
     &         PHI(3)*CHXE(3,4)+PHI(4)*CHXE(4,4)
        
        INTEGRAL=INTEGRAL+WEIGHTS(ng)*FVALUE*JACOBIAN
      ENDDO
      
      CALL EXITS('INTEGRATE_LINE')
      RETURN
 9999 CALL ERRORS('INTEGRATE_LINE',ERROR)
      CALL EXITS('INTEGRATE_LINE')
      RETURN 1
      END


