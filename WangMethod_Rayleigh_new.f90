
subroutine WangMethod_R(elip,Integral1,Integral2,IntegralI1,IntegralD2, W,SLOW)
! Implements the improved Wang's method for Rayleigh waves
use marc,only:G_SLOWS,G_SLOWP,ALFA,BTA,H,RHO,MU,NCAPAS
use globales,only:LONG_FLOAT,LONG_CMPLX,I,PI,DIAG,MATMUL1X1,EMATRIXNORMA
implicit none
! Inputs
real(LONG_FLOAT),intent(in)::W,SLOW
! Outputs
COMPLEX(LONG_CMPLX),intent(out)::elip
REAl(LONG_FLOAT),intent(out)::IntegralI1,IntegralD2,Integral1,Integral2
! Functions
complex(LONG_CMPLX)gam_val_from_slow,nu_val_from_slow
! Local variables
real(LONG_FLOAT)C,expo_prov1,expo_prov2,A_by_B_fac,exposEinvInt(3,ncapas),exposQacum(ncapas),inv_abs_y0_2,&
                norma,exps_v2(2,2)
complex(LONG_CMPLX)y0(4),L(4,4,ncapas),Einv1(4,4),Einv2(4,4),AuxMtrx(4,4,6,ncapas),Q(2,2),&
                   Qacum(2,2,ncapas),gam,nu,D(4,2),facD(4,4),deno(4,4),A_by_B,v2(2,2),p1,p2,p3,Linv(4,4)  

complex(LONG_CMPLX)Y(4,2)!Y(4,2,ncapas)

                   
! NO TENGO CLARO SI p1 Y p2 PUEDEN SER COMPLEJOS. p3 sí que puede
integer indx,index,indxx,IntegralStyle(ncapas)
complex(LONG_CMPLX)Integral(2),int_parcial(2),int_parcial2(2) 
! Body
L=0
AuxMtrx=0
exposEinvInt=0
exposQacum=0
expo_prov1=0
expo_prov2=0
IntegralStyle=0
Integral=0
IntegralD2=0
IntegralI1=0
do index=1,ncapas
  Qacum(1,1:2,index)=(/1.,0./)
  Qacum(2,1:2,index)=(/0.,1./)
enddo
!     Go!
c=1.d0/slow
index=ncapas;
gam=gam_val_from_slow(slow,w,REAL(G_SLOWP(index),KIND=long_float))
nu=nu_val_from_slow(slow,w,REAL(G_SLOWS(index),KIND=long_float))
call Lmatrix(L(:,:,index),alfa(index),bta(index),mu(index),w,c,gam,nu)
D=reshape((/(1.,0.),(0.,0.),(0.,0.),(0.,0.),(0.,0.),(1.,0.),(0.,0.),(0.,0.)/),(/4,2/))
Y=matmul(L(:,:,index),D)
AuxMtrx(1,1,1,index)=-1/(-gam+conjg(-gam))
AuxMtrx(2,2,2,index)=-1/(-nu+conjg(-nu))
AuxMtrx(2,1,3,index)=-1/(-nu+conjg(-gam))
AuxMtrx(1,1,4,index)=AuxMtrx(1,1,1,index)*(-gam*conjg(-gam));
AuxMtrx(2,2,5,index)=AuxMtrx(2,2,2,index)*(-nu*conjg(-nu));
AuxMtrx(2,1,6,index)=AuxMtrx(2,1,3,index)*(-nu*conjg(-gam));
exposEinvInt(:,index)=(/0,0,0/);
AuxMtrx(:,:,1,index)=matmul1x1(D(:,1),D(:,1))*AuxMtrx(:,:,1,index)
AuxMtrx(:,:,2,index)=matmul1x1(D(:,2),D(:,2))*AuxMtrx(:,:,2,index)
AuxMtrx(:,:,3,index)=matmul1x1(D(:,2),D(:,1))*AuxMtrx(:,:,3,index)
AuxMtrx(:,:,4,index)=matmul1x1(D(:,1),D(:,1))*AuxMtrx(:,:,4,index)
AuxMtrx(:,:,5,index)=matmul1x1(D(:,2),D(:,2))*AuxMtrx(:,:,5,index)
AuxMtrx(:,:,6,index)=matmul1x1(D(:,2),D(:,1))*AuxMtrx(:,:,6,index)
do index=ncapas-1,1,-1
  gam=gam_val_from_slow(slow,w,REAL(G_SLOWP(index),KIND=long_float))
  nu=nu_val_from_slow(slow,w,REAL(G_SLOWS(index),KIND=long_float))
  call Lmatrix(L(:,:,index),alfa(index),bta(index),mu(index),w,c,gam,nu)
  call LmatrixInv(Linv,alfa(index),bta(index),mu(index),w,c,gam,nu)
  D=matmul(Linv,Y)
  norma=1/sqrt(real(dot_product(D(:,1),D(:,1)))*real(dot_product(D(:,2),D(:,2))))
  Q=reshape((/D(2,2),-D(2,1),-D(1,2),D(1,1)/)*norma,(/2,2/))
  D=matmul(D,Q)
  D(1,2)=0
  D(2,1)=0
  Q(1,:)=Q(1,:)*exp(expo_prov2-expo_prov1) 
  do indx=ncapas,index+1,-1
    Qacum(:,:,indx)=matmul(Qacum(:,:,indx),Q)
    exposQacum(indx)=exposQacum(indx)-expo_prov2        
  enddo
  call EmatrixNorma(Einv1,expo_prov1,-H(index),gam,nu,1)
  Y(:,1)=matmul(matmul(L(:,:,index),Einv1),D(:,1))
  call EmatrixNorma(Einv2,expo_prov2,-H(index),gam,nu,2)
  Einv2(1,1)=0
  Y(:,2)=matmul(matmul(L(:,:,index),Einv2),D(:,2))
          
! For integrals:
  AuxMtrx(:,:,1,index)=matmul1x1(diag(Einv1),conjg(diag(Einv1)))-exp(-2*expo_prov1)
  AuxMtrx(:,:,2,index)=matmul1x1((/CMPLX(0,KIND=LONG_CMPLX),Einv2(2,2),Einv2(3,3),Einv2(4,4)/),&
  (/CMPLX(0,KIND=LONG_CMPLX),conjg(Einv2(2,2)),conjg(Einv2(3,3)),conjg(Einv2(4,4))/))-exp(-2*expo_prov2)
  AuxMtrx(:,:,3,index)=matmul1x1((/CMPLX(0,KIND=LONG_CMPLX),Einv2(2,2),Einv2(3,3),Einv2(4,4)/),conjg(diag(Einv1)))-&
  exp(-expo_prov1-expo_prov2)
  exposEinvInt(:,index)=(/2*expo_prov1,2*expo_prov2,expo_prov1+expo_prov2/)     
  
! deno=reshape((/-gam+conjg(-gam),-gam+conjg(-nu),-gam+conjg(gam),-gam+conjg(nu),-nu+conjg(-gam),-nu+conjg(-nu),&
!      -nu+conjg(gam),-nu+conjg(nu),gam+conjg(-gam),gam+conjg(-nu),gam+conjg(gam),gam+conjg(nu),nu+conjg(-gam),nu+conjg(-nu),&
!       nu+conjg(gam),nu+conjg(nu)/),(/4,4/),ORDER=(/2,1/))

  deno(1,1)=-2*real(gam);    deno(1,2)=-gam-conjg(nu);deno(1,3)=-2*aimag(gam)*I;  deno(1,4)=-gam+conjg(nu)
  deno(2,1)=conjg(deno(1,2));deno(2,2)=-2*real(nu);   deno(2,3)=-conjg(deno(1,4));deno(2,4)=-2*aimag(nu)*I
  deno(3,1)=-deno(1,3);      deno(3,2)=-deno(1,4);    deno(3,3)=-deno(1,1);       deno(3,4)=-deno(1,2)
  deno(4,1)=-deno(2,3);      deno(4,2)=-deno(2,4);    deno(4,3)=-deno(2,1);       deno(4,4)=-deno(2,2)
        
!facD=reshape((/gam*conjg(gam),gam*conjg(nu),-gam*conjg(gam),-gam*conjg(nu),nu*conjg(gam),nu*conjg(nu),&
!     -nu*conjg(gam),-nu*conjg(nu),-gam*conjg(gam),-gam*conjg(nu),gam*conjg(gam),gam*conjg(nu),-nu*conjg(gam),-nu*conjg(nu),&
!      nu*conjg(gam),nu*conjg(nu)/),(/4,4/),ORDER=(/2,1/))       

 facD(1,1)=gam*conjg(gam);  facD(1,2)=gam*conjg(nu);facD(1,3)=-facD(1,1);facD(1,4)=-facD(1,2);
 facD(2,1)=conjg(facD(1,2));facD(2,2)=nu*conjg(nu); facD(2,3)=-facD(2,1);facD(2,4)=-facD(2,2);
 facD(3,1)=facD(1,3);       facD(3,2)=facD(1,4);    facD(3,3)=facD(1,1); facD(3,4)=facD(1,2);
 facD(4,1)=facD(2,3);       facD(4,2)=facD(2,4);    facD(4,3)=facD(2,1); facD(4,4)=facD(2,2);

  do indx=1,4
    do indxx=1,4
      if (deno(indx,indxx)/=0) then
        ! Relaccionado con la integraci¢n de los cuadrados
        AuxMtrx(indx,indxx,1,index)= AuxMtrx(indx,indxx,1,index)/(-deno(indx,indxx))
        AuxMtrx(indx,indxx,2,index)= AuxMtrx(indx,indxx,2,index)/(-deno(indx,indxx))
        AuxMtrx(indx,indxx,3,index)= AuxMtrx(indx,indxx,3,index)/(-deno(indx,indxx))
      else
        ! Relaccionado con la integraci¢n de los cuadrados
        AuxMtrx(indx,indxx,1,index)=H(index)*exp(-exposEinvInt(1,index))
        AuxMtrx(indx,indxx,2,index)=H(index)*exp(-exposEinvInt(2,index))
        AuxMtrx(indx,indxx,3,index)=H(index)*exp(-exposEinvInt(3,index))
      endif
      ! Relaccionado con la integraci¢n de derivadas al cuadrado:
      AuxMtrx(indx,indxx,4,index)=AuxMtrx(indx,indxx,1,index)*facD(indx,indxx)
      AuxMtrx(indx,indxx,5,index)=AuxMtrx(indx,indxx,2,index)*facD(indx,indxx)
      AuxMtrx(indx,indxx,6,index)=AuxMtrx(indx,indxx,3,index)*facD(indx,indxx)!cruzado
    enddo
  enddo
  AuxMtrx(:,:,1,index)=matmul1x1(D(:,1),conjg(D(:,1)))*AuxMtrx(:,:,1,index)
  AuxMtrx(:,:,2,index)=matmul1x1(D(:,2),conjg(D(:,2)))*AuxMtrx(:,:,2,index)
  AuxMtrx(:,:,3,index)=matmul1x1(D(:,2),conjg(D(:,1)))*AuxMtrx(:,:,3,index)
  AuxMtrx(:,:,4,index)=matmul1x1(D(:,1),conjg(D(:,1)))*AuxMtrx(:,:,4,index)
  AuxMtrx(:,:,5,index)=matmul1x1(D(:,2),conjg(D(:,2)))*AuxMtrx(:,:,5,index)
  AuxMtrx(:,:,6,index)=matmul1x1(D(:,2),conjg(D(:,1)))*AuxMtrx(:,:,6,index)
enddo

A_by_B=-Y(3,2)/Y(3,1)
y0=A_by_B*Y(:,1)+Y(:,2)
inv_abs_y0_2=1/real(y0(2)*conjg(y0(2)))
elip=conjg(y0(1)/(i*y0(2)))
A_by_B_fac=exp(expo_prov2-expo_prov1)

!Integrals  
v2=matmul1X1((/A_by_B,cmplx(1,KIND=LONG_CMPLX)/),(/conjg(A_by_B),cmplx(1,KIND=LONG_CMPLX)/))
exps_v2=reshape((/2*(expo_prov2-expo_prov1),(expo_prov2-expo_prov1),(expo_prov2-expo_prov1),0.d0/),(/2,2/),ORDER=(/2,1/))
do index=1,ncapas
  if(index==1)then
    p1=v2(1,1)*exp(exps_v2(1,1)+exposEinvInt(1,index)-2*expo_prov2)
    p2=v2(2,2)*exp(exps_v2(2,2)+exposEinvInt(2,index)-2*expo_prov2)
    p3=v2(2,1)*exp(exps_v2(2,1)+exposEinvInt(3,index)-2*expo_prov2) 
  else
    p1=dot_product(Qacum(1,:,index),&
                  matmul(Qacum(1,:,index),v2*exp(exps_v2+exposEinvInt(1,index)-2*expo_prov2+(2*exposQacum(index)))))
    p2=dot_product(Qacum(2,:,index),&
                  matmul(Qacum(2,:,index),v2*exp(exps_v2+exposEinvInt(2,index)-2*expo_prov2+(2*exposQacum(index)))))
    p3=dot_product(Qacum(1,:,index),&
                  matmul(Qacum(2,:,index),v2*exp(exps_v2+exposEinvInt(3,index)-2*expo_prov2+(2*exposQacum(index)))))
  endif
  ! Kinetic energy integration conjg(L*AuxMtrx) 
  do indx=1,2
    int_parcial(indx)=&
      dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,1,index)))*p1+&
      dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,2,index)))*p2+&
      2*real(dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,3,index)))*p3)
    int_parcial2(indx)=&
      dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,4,index)))*p1+&
      dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,5,index)))*p2+&
      2*real(dot_product(L(indx,:,index),matmul(L(indx,:,index),AuxMtrx(:,:,6,index)))*p3)
  enddo
  Integral=Integral+(rho(index)*int_parcial)
  ! Integral I1
  IntegralI1=IntegralI1+mu(index)*(alfa(index)*alfa(index)/(bta(index)*bta(index))*int_parcial(1)+int_parcial(2))
  IntegralD2=IntegralD2+mu(index)*(alfa(index)*alfa(index)/(bta(index)*bta(index))*int_parcial2(2)+int_parcial2(1))
enddo
Integral=Integral*inv_abs_y0_2
Integral1=Integral(1)
Integral2=Integral(2)
IntegralI1=IntegralI1*inv_abs_y0_2
IntegralD2=IntegralD2*inv_abs_y0_2
end subroutine
