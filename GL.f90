
SUBROUTINE GL(IMG2,INDEX,OFFSET_L)

! Contribution of one mode to the Love-wave Green's function
USE marc,only:X,G_NX,VALUES_L
USE globales,only:long_float,long_cmplx,I,PI
IMPLICIT NONE
INTEGER INDEX,J,OFFSET_L
REAL(long_float)AL
COMPLEX(long_cmplx) Integral,Integralratio
REAL(long_float),INTENT(OUT)::IMG2(G_NX)
! Initializing
IMG2=0.D0
!$OMP PARALLEL DO PRIVATE(Integral,Integralratio,AL)      
DO J=1,G_NX-OFFSET_L
  CALL WangMethod_L(Integral,Integralratio,REAL(X(OFFSET_L+J),KIND=long_float),&
                    REAL(VALUES_L((INDEX-1)*G_NX+OFFSET_L+J),KIND=long_float))  
  AL=1.D0/(2.D0*Integralratio*Integral)
  VALUES_L((INDEX-1)*G_NX+OFFSET_L+J)=1.D0/(REAL(VALUES_L((INDEX-1)*G_NX+OFFSET_L+J),KIND=long_float)*Integralratio)
  IMG2(OFFSET_L+J)=0.5d0*AL
ENDDO
!$OMP END PARALLEL DO 
END SUBROUTINE

