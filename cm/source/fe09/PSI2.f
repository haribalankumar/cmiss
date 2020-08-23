      REAL*8 FUNCTION PSI2(IBT,INP,nb,nu,nn,XL)

C#### Function: PSI2
C###  Type: REAL*8
C###  Description:
C###    PSI2 evaluates Simplex basis functions at area coordinates XL.

      IMPLICIT NONE
      INCLUDE 'b01.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'cbdi02.cmn'
!     Parameter List
      INTEGER IBT(3,NIM),INP(NNM,NIM),nb,nn,nu
      REAL*8 XL(3)
!     Local Variables
      INTEGER IPU(15,4),L(4),li
      REAL*8 PSIM
      CHARACTER ERROR *10

      DATA IPU/1,2,3,1,1,2,1,1,2,1,1,1,2,1,1,
     '         1,1,1,2,3,2,1,1,1,2,1,1,1,2,1,
     '         1,1,1,1,1,1,2,3,2,2,1,1,1,1,2,
     '         1,1,1,1,1,1,1,1,1,1,2,3,2,2,2/

      L(2)=INP(nn,1)-1
      L(3)=INP(nn,2)-1
      L(4)=INP(nn,3)-1
      L(1)=IBT(2,1)-L(2)-L(3)
      PSI2=1.0D0
      DO li=1,NIT(nb)+1
        PSI2=PSI2*PSIM(L(li),IPU(nu,li),IBT(2,1),XL(li))
        IF(DOP) THEN
          WRITE(OP_STRING,'('' nu='',I2,'' nn='',I2,'' L(1)='',I2,'
     '      //''' L(2)='',I2,'' L(3)='',I2,'' L(4)='',I2,'' LI='','
     '      //'I2,'' XL='',E10.3,'' PSI2='',E10.3)')
     '      nu,nn,L(1),L(2),L(3),L(4),li,XL(li),PSI2
          CALL WRITES(IODI,OP_STRING,ERROR,*9999)
        ENDIF
      ENDDO

 9999 RETURN
      END


