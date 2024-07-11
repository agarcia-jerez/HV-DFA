      
      module TYPES
      integer,parameter::long_float=kind(0.0d0)
      integer,parameter::long_cmplx=kind((0.0d0,0.0d0))
      end module
      
      module globales
      USE TYPES
      complex(long_cmplx),parameter::i=(0,1)
      real(long_float),parameter::pi=3.141592654;
      real(long_float),parameter::twopi=2.d0*3.141592654;
      !logical SHOW_PH,SHOW_GR,SHOW_HV
      !integer TO_STD ! which output goes to the standard output; 0=report 1=phase 2=group 3=hv
      INTEGER UNIT_RDPGS(5) ! Units for: Report, Debug, Phase-vel., Group-vel., Spectra
      REAL SHDAMP,PSVDAMP
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! I N T E R F A C E S
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      interface matmul1x1
       ! multiplica vector (nx1)*(1xm) para dar matriz (nxm)
       function matmul1x1(A,B)
        USE TYPES
        complex(long_cmplx),intent(in):: A(:),B(:)
        complex(long_cmplx)matmul1x1(size(A),size(B))
       end function matmul1x1
      end interface
      interface diag
       ! Extracts the diagonal from a square matrix as a vector
       function diag(A)
        USE TYPES
        complex(long_cmplx),dimension(:,:),intent(in)::A
        complex(long_cmplx),dimension(size(A,1))::diag
       end function diag
      end interface
      interface det2x2
       ! Determinant of the 2x2 square matrix
        function det2x2(A)
        USE TYPES
        complex(long_cmplx)det2x2
        complex(long_cmplx),dimension(2,2),intent(in)::A
       end function det2x2
      end interface
      interface EmatrixNorma
       ! last argument is optional
       subroutine EmatrixNorma(salida,expo, z,gam,nu,tipo_norma)
        USE TYPES
        real(long_float),intent(in)::z
        complex(long_cmplx),intent(in)::gam,nu
        integer,optional::tipo_norma
        complex(long_cmplx),dimension(4,4),intent(out)::salida
        real(long_float),intent(out)::expo
       end subroutine
      end interface
      interface GenLmatrix
       ! c and w can be either real or complex
       subroutine Lmatrix(salida, alfa,bta,mu,w,c,gam,nu)        
        USE TYPES
        real(long_float),intent(in)::alfa,bta,mu,w,c
        complex(long_cmplx),intent(in)::gam,nu
        complex(long_cmplx),dimension(4,4),intent(out)::salida
       end subroutine
       subroutine CmplxWLmatrix(salida, alfa,bta,mu,w,c,gam,nu)        
        USE TYPES
        real(long_float),intent(in)::alfa,bta,mu
        complex(long_cmplx),intent(in)::w,c,gam,nu
        complex(long_cmplx),dimension(4,4),intent(out)::salida
       end subroutine        
      end interface
      interface GenLmatrixInv
       ! c and w can be either real or complex
       subroutine LmatrixInv(salida,alfa,bta,mu,w,c,gam,nu)
        use types
        real(long_float),intent(in)::alfa,bta,mu,w,c
        complex(long_cmplx),intent(in)::gam,nu
        complex(long_cmplx),dimension(4,4),intent(out)::salida
       end subroutine
       subroutine CmplxWLmatrixInv(salida,alfa,bta,mu,w,c,gam,nu)
        use types
        real(long_float),intent(in)::alfa,bta,mu
        complex(long_cmplx),intent(in)::w,c,gam,nu
        complex(long_cmplx),dimension(4,4),intent(out)::salida
       end subroutine        
      end interface       
      end module

      MODULE Marc
      USE TYPES
      REAL::G_x1,G_x2,G_y1,G_y2,G_POLARITY ! TO BRACKET THE ROOT AND STORE THE RESULT OF Y
      REAL G_OMEGA
      REAL::G_PRECISION ! MARC WAS USING 1.E-7
      REAL,ALLOCATABLE,TARGET::X(:) ! VECTOR OF CIRCULAR FREQUENCIES (OMEGA), GROWING
      REAL,ALLOCATABLE,TARGET::VALUES_L(:),VALUES_R(:) ! SLOWNESS
      LOGICAL,ALLOCATABLE,TARGET::VALID_L(:),VALID_R(:)
      INTEGER G_NX! NUMBER OF OMEGA VALUES.
      INTEGER G_NMODES! NUMBER OF MODES
      REAL::G_DX=0.1! SLOWNESS INCREMENT (DEFAULT 10%)
      LOGICAL::G_DXTYPE=.FALSE. ! TYPE OF INCREMENT .TRUE.=ABSOLUTE,.FALSE.=RELATIVE (DEFAULT)
      LOGICAL GO_DOWN ! .TRUE.= DECREMENTING SLOWNESS, .FALSE.=INCREMENTING SLOWNESS
      LOGICAL forbitVelocityInversion ! WILL BE SET TO .TRUE. IF THERE IS ANY INVERSION IN THE S WAVE VELOCITY LAYER MODEL
      REAL,ALLOCATABLE::G_SLOWS(:),G_SLOWP(:)
      REAL(LONG_FLOAT),ALLOCATABLE::ALFA(:),BTA(:),H(:),RHO(:),MU(:)
      INTEGER NCAPAS,G_VELOCITYINVERSION
      INTEGER,PARAMETER::ROOTSOLVER_MAX_ITERATIONS=500
      REAL SOLVER_RELATIVESTEP,SOLVER_SEARCHSTEP
      REAL G_MAXRAYLEIGHSLOWNESS,G_SLOWSMIN,G_SLOWSMAX
      INTEGER ILAYER                  
      REAL,PARAMETER::PI=3.141592654D0
      LOGICAL::ISRAYLEIGH ! TURN TO .FALSE. FOR LOVE WAVES
      ! GLOBALLY AVAILABLE FUNCTIONS:
      INTERFACE Y
        REAL FUNCTION Y(SLOWNESS) !WHOSE ROOTS REPRESENT DISPERSION CURVES
          REAL,INTENT(IN)::SLOWNESS
          ! IT ALSO USES G_OMEGA AND THE GROUND MODEL PARAMETERS
        END FUNCTION
      END INTERFACE
      INTERFACE HALFSPACE_RAYLEIGH
        FUNCTION HALFSPACE_RAYLEIGH(ILAYER)
        USE TYPES
        REAL(long_float)HALFSPACE_RAYLEIGH 
        INTEGER,INTENT(IN)::ILAYER
        END FUNCTION
      END INTERFACE   
      END MODULE
