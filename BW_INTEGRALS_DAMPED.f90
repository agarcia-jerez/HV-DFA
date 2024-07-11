SUBROUTINE BWR(SUMV,SUMPSV,SUMSH, Nks,W)

! Variables
USE TYPES
USE GLOBALES,ONLY:I,twopi,SHDAMP,PSVDAMP,EmatrixNorma,GenLmatrix,GenLmatrixInv
USE MARC,ONLY:NCAPAS,ALFA,BTA,H,MU
IMPLICIT NONE
INTEGER,INTENT(IN)::Nks ! number of positive ks for integration
REAL(LONG_FLOAT),INTENT(OUT)::SUMV,SUMPSV,SUMSH ! G11 and G33 for body waves
REAL(LONG_FLOAT)dk,k,K2,HPSV(Nks+1),V(Nks+1),SH(Nks+1)
REAL W
REAL(LONG_FLOAT)expo,norma,AUX1
INTEGER ik,il
COMPLEX(LONG_CMPLX)YSH(2,1),DSH(2,1),YPSV(4,2),DPSV(4,2),LSH(2,2),LSHINV(2,2),LPSV(4,4),LPSVINV(4,4),&
                   QPSV(2,2),EPSVINV(4,4),ESHINV(2,2),nua,nubsh,nubpsv
COMPLEX(LONG_CMPLX)CWSH,CCSH,CWPSV,CCPSV ! COMPLEX OMEGA, COMPLEX C

CWSH =CMPLX(W,-W*SHDAMP, KIND=LONG_CMPLX)
CWPSV=CMPLX(W,-W*PSVDAMP,KIND=LONG_CMPLX)

! Go!
dk=DBLE(w)/bta(NCAPAS)/Nks

! WANGMETHOD FOR BODY WAVES:

!DO ik=1,Nks+1
  IK=NKS+2
  K=DBLE(W)/BTA(NCAPAS)+DK  
1 CONTINUE  ! starting point of the loop
  K=K-DK
2 CONTINUE  !Insertion point for the last iteration
  IK=IK-1
  K2=K*K
  CCSH=CWSH/K;CCPSV=CWPSV/K;  
  ! Find out the signs of w, c
  nua=conjg(csqrt(cmplx(K2-(CWPSV/alfa(NCAPAS))**2)))
  nubPSV=conjg(csqrt(cmplx(K2-(CWPSV/bta(NCAPAS))**2)))
  nubSH =conjg(csqrt(cmplx(K2-(CWSH/ bta(NCAPAS))**2)))  
  CALL GenLmatrix(LPSV, alfa(NCAPAS),bta(NCAPAS),mu(NCAPAS),-cwPSV,-ccPSV,nua,nubPSV)
  CALL LmatrixSH(LSH, mu(NCAPAS),nubSH) !LSH=[1 1;-1*((nub(j))*mu(j)) (nub(j))*mu(j)];
  DPSV=RESHAPE((/1.,0.,0.,0., 0.,1.,0.,0./),(/4,2/))
  DSH(:,1)=(/1.,0./)
  YPSV=MATMUL(LPSV,DPSV)
  YSH=MATMUL(LSH,DSH)
  DO il=NCAPAS-1,1,-1
    nua=conjg(csqrt(cmplx(K2-(CWPSV/alfa(il))**2)))
    nubPSV=conjg(csqrt(cmplx(K2-(CWPSV/bta(il))**2)))
    nubSH =conjg(csqrt(cmplx(K2-(CWSH/ bta(il))**2)))    
    CALL GenLmatrix(LPSV, alfa(il),bta(il),mu(il),-cwPSV,-ccPSV,nua,nubPSV)!LPSV=LM(alfa(j),bta(j),mu(j),-w,-c,nua(j),nub(j));  
    CALL LmatrixSH(LSH, mu(il),nubSH)!LSH=[1 1;-1*((nub(j))*mu(j)) (nub(j))*mu(j)]
    CALL GenLmatrixInv(LPSVINV, alfa(il),bta(il),mu(il),-cwPSV,-ccPSV,nua,nubPSV)! LmatrixInv(salida,alfa,bta,mu,w,c,gam,nu)
    CALL LmatrixInv2(LSHINV, mu(il),nubSH)
    DPSV=MATMUL(LPSVINV,YPSV)      
    DSH=MATMUL(LSHINV,YSH)   

    ! Transforming and scaling D
    norma=1/sqrt(real(dot_product(DPSV(:,1),DPSV(:,1)))*real(dot_product(DPSV(:,2),DPSV(:,2))))
    QPSV=RESHAPE((/DPSV(2,2),-DPSV(2,1),-DPSV(1,2),DPSV(1,1)/)*norma,(/2,2/))
    DPSV=MATMUL(DPSV,QPSV) 
    DPSV(1,2)=0
    DPSV(2,1)=0
    DSH=DSH*(1/sqrt(real(dot_product(DSH(:,1),DSH(:,1)))))
  
    !Propagating upwards
    CALL EmatrixNorma(EPSVINV,expo,-H(il),nua,nubPSV,3)
    CALL EmatrixNorma2(ESHINV,expo,-H(il),nubSH) 
    YSH=MATMUL(MATMUL(LSH,ESHINV),DSH)
    YPSV=MATMUL(MATMUL(LPSV,EPSVINV),DPSV)
  ENDDO

  SH(ik)=REAL(I*k*YSH(1,1)/YSH(2,1))
  YPSV(:,1)=YPSV(:,1)/MAXVAL(ABS(YPSV(:,1)))
  YPSV(:,2)=YPSV(:,2)/MAXVAL(ABS(YPSV(:,2)))
  V(ik)=REAL(k*I*(YPSV(2,2)*YPSV(3,1)-YPSV(2,1)*YPSV(3,2))/(-YPSV(3,2)*YPSV(4,1)+YPSV(3,1)*YPSV(4,2)))
  HPSV(ik)=REAL(I*k*(-YPSV(1,2)*YPSV(4,1)+YPSV(1,1)*YPSV(4,2))/&
                        (-YPSV(3,2)*YPSV(4,1)+YPSV(3,1)*YPSV(4,2)))

!!!!! ENDDO
IF(IK.NE.2.AND.IK.NE.1)GOTO 1
!K=DK/10;IF(IK.EQ.2)GOTO 2
IF(IK.EQ.2)THEN
  K=DK*0.1d0
  GOTO 2
ENDIF

! Integrals
AUX1=-1/TWOPI*dk
SUMV  = AUX1*sum(V)
SUMPSV= AUX1*0.5*sum(HPSV)
SUMSH = AUX1*0.5*sum(SH)

END SUBROUTINE