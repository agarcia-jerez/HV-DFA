
      subroutine WangMethod_L(Integral,Integralratio, W,SLOW)
      use marc,only:G_SLOWS,H,RHO,MU,NCAPAS
      use globales,only:long_float,long_cmplx,I,PI,DIAG,MATMUL1X1
      implicit none
      ! Inputs
      real(long_float),intent(in)::W,SLOW
      ! Outputs
      complex(long_cmplx),intent(out)::Integral,Integralratio
      ! Functions
      complex(long_cmplx) nu_val_from_slow
      ! Local variables
      integer interfaz_fuente,index,indx,indxx
      real(long_float)::expo,exposEinvInt(ncapas),exposQacum(ncapas),Q,Qacum(ncapas),inv_abs_y0_2,p1
      complex(long_cmplx)::Y(2),L(2,2,ncapas),AuxMtrx(2,2,ncapas),D(2),Linv(2,2),Einv1(2,2),Integralmu,&
                           deno(2,2),int_parcial,nu
      ! Body
      Qacum=1
      exposQacum=0
      AuxMtrx=0
      expo=0
!     L=0
      interfaz_fuente=1
      Integral=0
      Integralmu=0
      index=ncapas
!      IF SLOW.NE.REAL(G_SLOWS(index),KIND=long_float)
        nu=nu_val_from_slow(SLOW,W,REAL(G_SLOWS(index),KIND=long_float))
!        write(*,*)real(nu),aimag(nu)
!      ELSE
        
!      ENDIF        
!     IF(REAL(NU)<=0)THEN
!       WRITE(*,*)INDEX,W,SLOW,REAL(G_SLOWS(index),KIND=long_float),NU
!       WRITE(*,*)SLOW.Eq.REAL(G_SLOWS(index),KIND=long_float)
!       READ(*,*)
!     ENDIF
      exposEinvInt(index)=0            
      call LmatrixSH(L(:,:,index), mu(index),nu)
      D=(/1,0/)
      Y=matmul(L(:,:,index),D)
      AuxMtrx(1,1,index)=-1/(-nu+conjg(-nu))
      !AuxMtrx(:,:,index)=matmul1x1(D,conjg(D))*AuxMtrx(:,:,index)
      do index=ncapas-1,1,-1
        nu=nu_val_from_slow(SLOW,W,REAL(G_SLOWS(index),KIND=long_float))
        call LmatrixSH(L(:,:,index), mu(index),nu)
        call LmatrixInv2(Linv, mu(index),nu)
        D=matmul(Linv,Y)
        Q=1/sqrt(dot_product(D,D))
        do indx=ncapas,index+1,-1
          Qacum(indx)=Qacum(indx)*Q
          exposQacum(indx)=exposQacum(indx)-expo
        enddo
        D=D*Q
        call EmatrixNorma2(Einv1,expo, -H(index),nu)
        Y=matmul(matmul(L(:,:,index),Einv1),D)
        ! Related with integration:
        AuxMtrx(:,:,index)=matmul1x1(diag(Einv1),conjg(diag(Einv1)))-exp(-2*expo)
        exposEinvInt(index)=2*expo
        deno(2,1)= 2*I*aimag(nu);deno(2,2)=2*real(nu);deno(1,1)=-deno(2,2);deno(1,2)=-deno(2,1);        
        do indx=1,2
          do indxx=1,2
            if (deno(indx,indxx)/=0)then
              AuxMtrx(indx,indxx,index)=AuxMtrx(indx,indxx,index)/(-deno(indx,indxx))
            else
              AuxMtrx(indx,indxx,index)=H(index)*exp(-exposEinvInt(index))
            endif
          enddo
        enddo
        AuxMtrx(:,:,index)=matmul1x1(D,conjg(D))*AuxMtrx(:,:,index)
      enddo
      inv_abs_y0_2=1/(Y(1)*conjg(Y(1)));
      do index=1,ncapas
        ! Integrations
        if (index==1)then
          p1=exp(exposEinvInt(index)-2*expo)
        else          
          p1=Qacum(index)*(exp(exposEinvInt(index)-2*expo+2*exposQacum(index)))*Qacum(index)
        endif
        int_parcial=dot_product(L(1,:,index),matmul(L(1,:,index),AuxMtrx(:,:,index)))*p1
        Integral=Integral+rho(index)*int_parcial
        Integralmu=Integralmu+mu(index)*int_parcial
      enddo
      Integral=Integral*inv_abs_y0_2
      Integralmu=Integralmu*inv_abs_y0_2
      Integralratio=Integralmu/Integral
      end subroutine

