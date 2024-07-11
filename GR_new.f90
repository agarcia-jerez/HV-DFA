
SUBROUTINE GR(G1,IMG3, INDEX,OFFSET_R)
  
! Contribution of one mode (TR,CR) to the Rayleigh-wave Green's function
use marc,only:X,G_NX,VALUES_R
use globales,only:long_float,long_cmplx,I,PI
IMPLICIT NONE
INTEGER,INTENT(IN)::INDEX,OFFSET_R
REAL(long_float),intent(out)::G1(G_NX),IMG3(G_NX)
REAL,POINTER::XMOD(:),SLOW(:)
REAL(long_float)DSR,DW ! Doble precission Rayleigh Slowness and circular frequency
INTEGER J
REAL(long_float) I1,I2,II1,ID2,vg,AR
COMPLEX(long_cmplx)u0_w0

! Initializing...
G1(1:OFFSET_R)=0.d0
IMG3(1:OFFSET_R)=0.d0
XMOD=>X(OFFSET_R+1:G_NX)
SLOW=>VALUES_R((INDEX-1)*G_NX+OFFSET_R+1:INDEX*G_NX)

!!$OMP PARALLEL DO PRIVATE(u0_w0,I1,I2,II1,ID2,VG,AR,DSR,DW)
DO j=1,G_NX-OFFSET_R
  DSR=REAL(SLOW(J),KIND=long_float)
  DW =REAL(XMOD(J),KIND=long_float)
  CALL WANGMETHOD_R(u0_w0,I1,I2,II1,ID2,REAL(XMOD(J),KIND=long_float),REAL(SLOW(J),KIND=long_float))
  vg=0.5D0*DSR*(  ((1.D0/DSR)**2*(I1+I2)-1.D0/(DSR*DW)**2*ID2)  +II1)/(I1+I2)
  SLOW(j)=1/vg
  AR=DSR/(2.D0*vg*(I1+I2))
  G1(OFFSET_R+j)=-0.5d0*AR*AIMAG(u0_w0)*AIMAG(u0_w0)
  IMG3(OFFSET_R+j)=-0.5d0*AR
ENDDO
!!$OMP END PARALLEL DO
END SUBROUTINE
