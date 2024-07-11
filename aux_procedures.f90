
      function matmul1x1(A,B) ! interface required
      ! Matrix multiplication of vectors: (nrow x 1) (1 x ncol)
      USE TYPES
      implicit none
      complex(long_cmplx),dimension(:),intent(in)::A,B
      complex(long_cmplx),dimension(size(A),size(B))::matmul1x1
      integer indxx,indx
      do indxx=1,size(A)
       do indx=1,size(B)
         matmul1x1(indxx,indx)=A(indxx)*B(indx)
       enddo
      enddo
      return
      end function

      function diag(A) ! Interface required
      ! extracts the diagonal of a complex square matriz
      USE TYPES
      implicit none
      complex(long_cmplx),dimension(:,:),intent(in)::A
      complex(long_cmplx),dimension(size(A,1))::diag
      integer index
      do index=1,size(diag)
        diag(index)=A(index,index)
      enddo
      return
      end function

      function gam_val(c,w,alfa)
      use globales,only:long_float,long_cmplx,i
      implicit none      
      complex(long_cmplx)gam_val            
      real(long_float),intent(in)::c,w,alfa
      gam_val=-i*w/alfa*csqrt(cmplx(1-alfa*alfa/(c*c)))
      end function

      function gam_val_from_slow(slow,w,slowP)
      use globales,only:long_float,long_cmplx,i
      implicit none      
      complex(long_cmplx)gam_val_from_slow           
      real(long_float),intent(in)::slow,w,slowP
      IF(slow.NE.slowP)THEN
        gam_val_from_slow=-i*w*slowP*sqrt(cmplx(1-slow*slow/(slowP*slowP),KIND=long_cmplx))!double precission
        !gam_val_from_slow=-i*w*slowP*csqrt(cmplx(1-slow*slow/(slowP*slowP)))!single precission
      ELSE
        gam_val_from_slow=-i*w*slowP*csqrt(cmplx(-tiny(0.)))
      ENDIF
      end function

      function nu_val(c,w,bta)
      use globales,only:long_float,long_cmplx,i
      implicit none
      complex(long_cmplx)nu_val      
      real(long_float),intent(in)::c,w,bta
      nu_val=-i*w/bta*csqrt(cmplx(1-bta*bta/(c*c)))
      end function   
      
      function nu_val_from_slow(slow,w,slowS)
      use globales,only:long_float,long_cmplx,i
      implicit none
      complex(long_cmplx)nu_val_from_slow     
      real(long_float),intent(in)::slow,w,slowS
      IF(slow.NE.slowS)THEN
        nu_val_from_slow=-i*w*slowS*sqrt(cmplx(1-slow*slow/(slowS*slowS),KIND=long_cmplx)) !double precission
        !nu_val_from_slow=-i*w*slowS*csqrt(cmplx(1-slow*slow/(slowS*slowS))) ! single precission
      ELSE
        nu_val_from_slow=-i*w*slowS*csqrt(cmplx(-tiny(0.)))
      ENDIF
      end function

      function det2x2(A)
      use globales,only:long_cmplx
      complex(long_cmplx)det2x2
      complex(long_cmplx),intent(in)::A(2,2)
      det2x2=A(1,1)*A(2,2)-A(2,1)*A(1,2)
      end function

!     Rayleigh

      subroutine EmatrixNorma(salida,expo, z,gam,nu,tipo_norma)
!     tipo_norma=1 normalize by the diverging gamma-type term
!     tipo_norma=2 normalize by the diverging nu-type term
!     tipo_norma=3 normaliza by the fastest diverging term
      use globales,only:long_float,long_cmplx,i
      implicit none
      real(long_float),intent(in)::z
      complex(long_cmplx),intent(in)::gam,nu
      integer,optional::tipo_norma
      complex(long_cmplx),dimension(4,4),intent(out)::salida
      real(long_float),intent(out)::expo
      ! local
      real(long_float)rg,ig,rnu,inu
      complex(long_cmplx)aux1
      rg=real(gam)*z
      ig=aimag(gam)*z
      rnu=real(nu)*z
      inu=aimag(nu)*z
      salida=0
      if(.not.present(tipo_norma).or.tipo_norma==3)then
       if (rg>0) then
        if (rg > rnu) then
         aux1=exp(i*inu)
         salida(3,3)=exp(i*ig)
         salida(1,1)=conjg(salida(3,3))*exp(-2*rg)
         salida(2,2)=conjg(aux1)*exp(-(rnu+rg))        
         salida(4,4)=aux1*exp((rnu-rg))
         expo=rg
        else
         aux1=exp(i*ig)      
         salida(4,4)=exp(i*inu)          
         salida(1,1)=conjg(aux1)*exp(-(rg+rnu))
         salida(2,2)=conjg(salida(4,4))*exp(-2*rnu)
         salida(3,3)=aux1*exp((rg-rnu))
         expo=rnu
        endif
       else
        if (rg < rnu) then
         aux1=exp(i*inu)
         salida(1,1)=exp(-i*ig)
         salida(2,2)=conjg(aux1)*exp(-rnu+rg)
         salida(3,3)=conjg(salida(1,1))*exp(2*rg)
         salida(4,4)=aux1*exp(rnu+rg)
         expo=-rg
        else
         aux1=exp(i*ig)
         salida(1,1)=conjg(aux1)*exp(-rg+rnu)
         salida(2,2)=exp(-i*inu)
         salida(3,3)=aux1*exp(rg+rnu)
         salida(4,4)=conjg(salida(2,2))*exp(2*rnu)
         expo=-rnu
        endif
       endif
      else
       if (rg>0)then
        if (tipo_norma==1)then! all by exp(rg)
          aux1=exp(i*inu)
          salida(3,3)=exp(i*ig)         
          salida(1,1)=conjg(salida(3,3))*exp(-2*rg)
          salida(2,2)=conjg(aux1)*exp(-(rnu+rg))
          salida(4,4)=aux1*exp((rnu-rg))
          expo=rg
        else ! all by exp(rnu)
          aux1=exp(i*ig)
          salida(1,1)=conjg(aux1)*exp(-(rg+rnu))
          salida(4,4)=exp(i*inu)
          salida(2,2)=conjg(salida(4,4))*exp(-2*rnu)
          salida(3,3)=aux1*exp((rg-rnu))         
          expo=rnu
        endif
       else
        if(tipo_norma==1)then!% all by exp(-rg)
          aux1=exp(i*inu)
          salida(1,1)=exp(-i*ig)
          salida(2,2)=conjg(aux1)*exp(-rnu+rg)
          salida(3,3)=conjg(salida(1,1))*exp(2*rg)
          salida(4,4)=aux1*exp(rnu+rg)
          expo=-rg
        else ! % all by exp(-rnu)
          aux1=exp(i*ig)
          salida(1,1)=conjg(aux1)*exp(-rg+rnu)
          salida(2,2)=exp(-i*inu)
          salida(3,3)=aux1*exp(rg+rnu)
          salida(4,4)=conjg(salida(2,2))*exp(2*rnu)
          expo=-rnu
        endif
       endif
      endif
      end subroutine

      subroutine LmatrixInv(salida,alfa,bta,mu,w,c,gam,nu)
      use types
      implicit none
      real(long_float),intent(in)::alfa,bta,mu,w,c
      complex(long_cmplx),intent(in)::gam,nu
      complex(long_cmplx),dimension(4,4),intent(out)::salida
      ! local
      real(long_float)k
      complex(long_cmplx)aux1
      real(long_float)aux2,aux3
      ! Go      
      k=w/c
      aux1=0.5*bta*(k*k+nu*nu)
      aux2=bta*bta/alfa
      aux3=1/w
      salida(1,1)=aux2*k
      salida(3,2)=bta*aux1/(alfa*gam)
      salida(1,2)=-salida(3,2)
      salida(1,4)=0.5*aux2/mu
      salida(1,3)=-salida(1,4)*(k/gam)  !-0.5*bta*bta*k/(alfa*gam*mu)      
      salida(2,1)=-aux1/nu
      salida(2,2)=+bta*k
      salida(2,3)=+0.5*bta/mu
      salida(2,4)=-salida(2,3)*(k/nu)
      salida(3,1)=salida(1,1)
      salida(3,3)=0.5*salida(1,1)/(gam*mu) !+0.5*bta*bta*k/(alfa*gam*mu)
      salida(3,4)=salida(1,4)
      salida(4,1)=salida(2,1)
      salida(4,2)=-salida(2,2)!-bta*k
      salida(4,3)=-salida(2,3)
      salida(4,4)=salida(2,4)
      salida=salida*aux3
      end subroutine
      
      subroutine CmplxWLmatrixInv(salida,alfa,bta,mu,w,c,gam,nu)
      use types
      implicit none
      real(long_float),intent(in)::alfa,bta,mu
      complex(long_cmplx),intent(in)::w,c,gam,nu
      complex(long_cmplx),dimension(4,4),intent(out)::salida
      ! local
      real(long_float)k
      complex(long_cmplx)aux1,aux3
      real(long_float)aux2
      ! Go      
      k=real(w,KIND=long_float)/real(c,KIND=long_float)
      aux1=0.5*bta*(k*k+nu*nu)
      aux2=bta*bta/alfa
      aux3=1/w
      salida(1,1)=aux2*k
      salida(3,2)=bta*aux1/(alfa*gam)
      salida(1,2)=-salida(3,2)
      salida(1,4)=0.5*aux2/mu
      salida(1,3)=-salida(1,4)*(k/gam)  !-0.5*bta*bta*k/(alfa*gam*mu)      
      salida(2,1)=-aux1/nu
      salida(2,2)=+bta*k
      salida(2,3)=+0.5*bta/mu
      salida(2,4)=-salida(2,3)*(k/nu)
      salida(3,1)=salida(1,1)
      salida(3,3)=0.5*salida(1,1)/(gam*mu) !+0.5*bta*bta*k/(alfa*gam*mu)
      salida(3,4)=salida(1,4)
      salida(4,1)=salida(2,1)
      salida(4,2)=-salida(2,2)!-bta*k
      salida(4,3)=-salida(2,3)
      salida(4,4)=salida(2,4)
      salida=salida*aux3
      end subroutine

      subroutine Lmatrix(salida, alfa,bta,mu,w,c,gam,nu)
      use types      
      implicit none
      real(long_float),intent(in)::alfa,bta,mu,w,c
      complex(long_cmplx),intent(in)::gam,nu
      complex(long_cmplx),dimension(4,4),intent(out)::salida
      complex(long_cmplx)aux1,aux2
      ! local
      real(long_float) k
      ! Go
      k=w/c
      aux1=mu*(k**2+nu**2)/w
      aux2=2*mu*k          
      salida(1,1)=alfa*k/w
      salida(1,2)=bta*nu/w
      salida(1,3)=salida(1,1)
      salida(1,4)=salida(1,2)
      salida(2,1)=alfa*gam/w
      salida(2,2)=bta*k/w
      salida(2,3)=-salida(2,1)
      salida(2,4)=-salida(2,2)
      salida(3,3)=aux2*salida(2,1) !2*alfa*mu*k*gam/w
      salida(3,1)=-salida(3,3)
      salida(3,4)=bta*aux1      
      salida(3,2)=-salida(3,4)    
      salida(4,1)=-alfa*aux1
      salida(4,2)=-aux2*salida(1,2)
      salida(4,3)=salida(4,1)
      salida(4,4)=salida(4,2)!-2*bta*mu*k*nu/w      
      end subroutine

      subroutine CmplxWLmatrix(salida, alfa,bta,mu,w,c,gam,nu)
      use types      
      implicit none
      real(long_float),intent(in)::alfa,bta,mu
      complex(long_cmplx),intent(in)::w,c,gam,nu
      complex(long_cmplx),dimension(4,4),intent(out)::salida
      complex(long_cmplx)aux1,aux2
      ! local
      real(long_float) k
      ! Go
      k=real(w,KIND=long_float)/real(c,KIND=long_float)
      aux1=mu*(k**2+nu**2)/w
      aux2=2*mu*k          
      salida(1,1)=alfa*k/w
      salida(1,2)=bta*nu/w
      salida(1,3)=salida(1,1)
      salida(1,4)=salida(1,2)
      salida(2,1)=alfa*gam/w
      salida(2,2)=bta*k/w
      salida(2,3)=-salida(2,1)
      salida(2,4)=-salida(2,2)
      salida(3,3)=aux2*salida(2,1) !2*alfa*mu*k*gam/w
      salida(3,1)=-salida(3,3)
      salida(3,4)=bta*aux1      
      salida(3,2)=-salida(3,4)    
      salida(4,1)=-alfa*aux1
      salida(4,2)=-aux2*salida(1,2)
      salida(4,3)=salida(4,1)
      salida(4,4)=salida(4,2)!-2*bta*mu*k*nu/w      
      end subroutine
      
!     Love Waves

      SUBROUTINE EmatrixNorma2(salida,expo, z,nu)
      USE Globales,ONLY:long_float,long_cmplx,i
      IMPLICIT NONE
      COMPLEX(long_cmplx),INTENT(OUT)::SALIDA(2,2)
      REAL(long_float),INTENT(OUT)::EXPO
      COMPLEX(long_cmplx),INTENT(IN)::NU
      REAL(long_float),INTENT(IN)::Z
      REAL(long_float)::RNU,INU
      RNU=REAL(NU)*Z;
      INU=AIMAG(NU)*Z;
      IF(RNU>0)THEN
       SALIDA(2,2)=EXP(I*INU)
       SALIDA(1,1)=CONJG(SALIDA(2,2))*EXP(-2*RNU)
       SALIDA(2,1)=0
       SALIDA(1,2)=0
       EXPO=RNU
      ELSE
       SALIDA(1,1)=EXP(-I*INU)
       SALIDA(2,1)=0
       SALIDA(1,2)=0
       SALIDA(2,2)=CONJG(SALIDA(1,1))*EXP(2*RNU)
       EXPO=-RNU
      ENDIF
      END SUBROUTINE
      
      SUBROUTINE LmatrixSH(salida, mu,nu)
      USE Globales,ONLY:long_float,long_cmplx
      IMPLICIT NONE
      COMPLEX(long_cmplx),INTENT(OUT)::SALIDA(2,2)
      COMPLEX(long_cmplx),INTENT(IN)::NU
      REAL(long_float),INTENT(IN)::MU
      SALIDA(1,1)=CMPLX(1);
      SALIDA(2,2)=nu*mu
      SALIDA(2,1)=-SALIDA(2,2)
      SALIDA(1,2)=CMPLX(1)
      END SUBROUTINE

      SUBROUTINE LmatrixInv2(salida,mu,nu)
      USE Globales,ONLY:long_float,long_cmplx
      IMPLICIT NONE
      COMPLEX(long_cmplx),INTENT(OUT)::SALIDA(2,2)
      COMPLEX(long_cmplx),INTENT(IN)::NU
      REAL(long_float),INTENT(IN)::MU
      SALIDA(1,1)=CMPLX(0.5);
      SALIDA(2,2)=1./(2.*nu*mu)
      SALIDA(1,2)=-SALIDA(2,2)
      SALIDA(2,1)=CMPLX(0.5)
      END SUBROUTINE
      
      FUNCTION HALFSPACE_RAYLEIGH(iLayer)
      USE Marc,ONLY:G_SLOWP,G_SLOWS
      USE Types
      IMPLICIT NONE
      REAL(long_float)::HALFSPACE_RAYLEIGH
      INTEGER,INTENT(IN)::ILAYER
      real(long_float)::x
      real(long_float)::boa,a3,a4,sqx,y,dy
!     Equation : (2-sq_cob)-4*sqrt(1-sq_boa*sq_cob)*sqrt(1-sq_cob)=0
!     The sqrt can be elevated to square, interiors are positive: [(2-x)]=16* (1-sq_boa*x) * (1-x), where x=sq_cob
!     Calculating all terms: x * [x-8*x+(24-16*sq_boa) *x+(16*sq_boa-16)]=0
!     x=0 is not a valid solution Thus we are looking for the roots of: x-8*x+a3*x+a4 where, a4=16*sq_boa-16 a3=8-a4
!     3rd degree poly which solutions are calculated by newton raphson using the derivative : 3*x-16*x+a3
!     (c/beta)=0.95 solution for a=infinity, this is a maximum bound value for the Rayleigh velocity */
      X=0.9025
      boa=G_slowP(iLayer)/G_slowS(iLayer)
      a4=16*boa*boa-16;
      a3=8-a4;      
      sqx=x*x
      y=x*sqx-8*sqx+a3*x+a4
      dy=3*sqx-16*x+a3
      x=x-y/dy
      do while(abs(y)>1e-10)! This condition over y achieves a precision over x of about 1e-15
       sqx=x*x
       y=x*sqx-8*sqx+a3*x+a4
       dy=3*sqx-16*x+a3
       x=x-y/dy
      enddo
      HALFSPACE_RAYLEIGH=G_SLOWS(iLayer)/sqrt(x)
      end function
