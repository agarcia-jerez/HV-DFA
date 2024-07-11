
      REAL FUNCTION Y(SLOW)
      USE Marc,ONLY:ISRAYLEIGH
      USE Globales,ONLY:long_float,long_cmplx
      REAL,INTENT(IN)::SLOW
      REAL(long_float)YL,YR
      IF(ISRAYLEIGH)THEN
        Y=real(YR(SLOW),kind=kind(Y))
      ELSE
        Y=real(YL(SLOW),kind=kind(Y))
      ENDIF
      END FUNCTION

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!                R A Y L E I G H                 !!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function yR(slow)
      use Marc,ONLY:G_OMEGA,ALFA,BTA,H,MU,NCAPAS
      use Globales,ONLY:i,pi,twopi,det2x2,EmatrixNorma
      Use Types 
      IMPLICIT NONE                  
      real(long_float)yR
      REAL,INTENT(IN)::SLOW
      complex(long_cmplx)gam_val,nu_val
      ! Local variables
      real(long_float)w,c,exposYb,expo_prov1,expo_prov2,norma
      complex(long_cmplx)Yb(4,2),L(4,4),Linv(4,4),Einv1(4,4),Einv2(4,4),Q(2,2),Qacum(2,2),gam,nu,D(4,2)
      integer index
      ! Initializing
      expo_prov1=0;expo_prov2=0
      ! Body
      w=DBLE(G_OMEGA)
      c = 1.D0 / DBLE(slow)
      exposYb=0
      Qacum = reshape((/1.,0.,0.,1./),(/2,2/))
      ! Go!
      index=ncapas;
      gam=gam_val(c,w,alfa(index));
      nu=nu_val(c,w,bta(index));
      call Lmatrix(L,alfa(index),bta(index),mu(index),w,c,gam,nu)
      D = reshape((/1.,0.,0.,0., 0.,1.,0.,0./),(/4,2/))
      Yb=matmul(L,D)
      do index=ncapas-1,1,-1
        gam=gam_val(c,w,alfa(index))
        nu=nu_val(c,w,bta(index))
        call Lmatrix(L,alfa(index),bta(index),mu(index),w,c,gam,nu)
        call LmatrixInv(Linv,alfa(index),bta(index),mu(index),w,c,gam,nu)
        D=matmul(Linv,Yb)
        norma=1/sqrt(real(dot_product(D(:,1),D(:,1)))*real(dot_product(D(:,2),D(:,2))))
        Q=reshape((/D(2,2),-D(2,1),-D(1,2),D(1,1)/)*norma,(/2,2/))
        D=matmul(D,Q);D(1,2)=0;D(2,1)=0
        Qacum=matmul(Qacum,Q)
        call EmatrixNorma(Einv1,expo_prov1,-H(index),gam,nu,1)
        Yb(:,1)=matmul(matmul(L,Einv1),D(:,1))
        call EmatrixNorma(Einv2,expo_prov2,-H(index),gam,nu,2)
        Einv2(1,1)=0
        Yb(:,2)=matmul(matmul(L,Einv2),D(:,2))
      enddo
      yR=real(det2x2(matmul(Yb(3:4,:),conjg(Qacum))))      
      return
      end function


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!                   L O V E                      !!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION YL(SLOW)
      USE Marc,ONLY:G_OMEGA,BTA,H,MU,NCAPAS
      USE Globales,ONLY:long_float,long_cmplx,i,pi,twopi
      implicit none
      real(long_float)yl
      real,intent(in)::slow
      complex(long_cmplx) nu_val
      ! Local variables
      integer index
      real(long_float)::w,c,expo_prov,Q!,Qacum,exposQacum,exposYb
      complex(long_cmplx)::Yb(2,1),L(2,2),Linv(2,2),nu,D(2,1),Einv1(2,2)
      ! Body
      w=G_OMEGA
      c=1/slow
      index=ncapas
      nu=nu_val(c,w,bta(index));
      CALL LmatrixSH(L,mu(index),nu)
      D(:,1)=(/1.,0./)
      Yb=matmul(L,D)
      do index=ncapas-1,1,-1
        nu=nu_val(c,w,bta(index));
        CALL LmatrixSH(L,mu(index),nu)
        CALL LmatrixInv2(Linv,mu(index),nu)
        D=matmul(Linv,Yb)
        Q=1/sqrt(real(D(1,1)*conjg(D(1,1))+D(2,1)*conjg(D(2,1))))
        D=D*Q
        CALL EmatrixNorma2(Einv1,expo_prov,-H(index),nu)
        Yb=matmul(matmul(L,Einv1),D)
      enddo
      yL=real(Yb(2,1))
      RETURN
      END FUNCTION

      
