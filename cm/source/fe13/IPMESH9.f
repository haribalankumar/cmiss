      SUBROUTINE IPMESH9(NBJ,NEELEM,NENP,NKJE,NKJ,
     '  NP_INTERFACE,NPNE,NPNODE,nr,NRE,NVJE,NVJP,MRNA,RNA,
     '  SE,XP,ERROR,*)

C#### Subroutine: IPMESH9
C###  Description:
C###    IPMESH9 defines mesh parameters for DNA mesh.

      IMPLICIT NONE
      INCLUDE 'b00.cmn'
      INCLUDE 'b01.cmn'
      INCLUDE 'cbdi02.cmn'
      INCLUDE 'geom00.cmn'
      INCLUDE 'inout00.cmn'
      INCLUDE 'loc00.cmn'
      INCLUDE 'loc00.inc'
      INCLUDE 'mesh01.cmn'
      INCLUDE 'mxch.inc'

!     Parameter List
      INTEGER NBJ(NJM,NEM),NEELEM(0:NE_R_M,0:NRM),
     '  NENP(NPM,0:NEPM,0:NRM),
     '  NKJE(NKM,NNM,NJM,NEM),NKJ(NJM,NPM),
     '  NPNE(NNM,NBFM,NEM),NPNODE(0:NP_R_M,0:NRM),
     '  NP_INTERFACE(0:NPM,0:3),nr,NRE(NEM),
     '  NVJE(NNM,NBFM,NJM,NEM),NVJP(NJM,NPM)
      REAL*8 MRNA(3),RNA(3),SE(NSM,NBFM,NEM),XP(NKM,NVM,NJM,NPM)
      CHARACTER ERROR*(*)
!     Local Variables
      INTEGER base,base_pairs,first_mover,PAIR(40),
     '  i,IBEG2,IBEG3,IEND2,IEND3,ICHAR,INFO,last_mover,M,m1,m2,m3,
     '  MRNA_start,MRNA_nodes,
     '  node_num_start,node_num_start_1,node_num_start_2,
     '  N,N12,nb,nb1,nc,ne,ne1,NE2,NE3,ni,nj,nk,nn,nnp,ns,np,
     '  np2,noelem,nonode,NOQUES
      REAL*8 angle,angle_inc,chain,
     '  BASE_1(400),BASE_2(400),delta_x,DNA_radius,head,space,
     '  radius,RDEFAULT(3,8),spacer,
     '  RNA_angle,RNA_radius,RNA_x, RNA_y, RNA_z,
     '  SIDE(3,3),SIDE_COARSE(3,3),SIDE_FINE(3,3),
     '  SIDE_LENGTH(3),SIDE_TOT(3,3),SIDE_X(3),tail,x_inc,x,y,z
      CHARACTER CDEFAULT1(2)*1,CDEFAULT2(4)*3,CDEFAULT3(8)*5,
     '  CHAR1*1,CHAR2*2,CHAR3*5
      LOGICAL FILEIP,SAMETYPE,SAMENUMXI
      DATA RDEFAULT/0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,
     '              1.0d0, 0.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0,
     '              1.0d0, 1.0d0, 0.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0,
     '              1.0d0, 1.0d0, 1.0d0/
      DATA CDEFAULT1/'0','1'/,
     '     CDEFAULT2/'0,0','1,0','0,1','1,1'/,
     '     CDEFAULT3/'0,0,0','1,0,0','0,1,0','1,1,0',
     '               '0,0,1','1,0,1','0,1,1','1,1,1'/

      CALL ENTERS('IPMESH9',*9999)

C!!!! CS 2/9/2003 removed for now....

      CALL EXITS('IPMESH9')
      RETURN
 9999 CALL ERRORS('IPMESH9',ERROR)
      CALL EXITS('IPMESH9')
      RETURN 1
      END

