
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                  !
!                          H    H   V     V         DDDD     FFFFF   AAAAAA                                        !
!                          H    H   V     V         D   D    F       A    A                                        !
!                          HHHHHH    V   V    ===   D    D   FFFFF   AAAAAA                                        !
!                          H    H     V V           D   D    F       A    A                                        !
!                          H    H      V            DDDD     F       A    A                                        !
!                                                                                                                  ! 
!                                            version 1.0                                                           !
!                                                                                                                  !
!       This is an implementation of the method by Sanchez-Sesma et al (2011, Geoph. J. Int.) ,                    !
!       using a propagator matrix scheme based on Wang (1999, Bull, Seism. Soc. Am.)and contour                    !
!       integration in the complex-wavenumber plane.                                                               !  
!                                                                                                                  !                                   
!                            HV-INV Project Team (UNAM - University of Almeria)                                    !
!                                                                                                                  !
!                              (c) Antonio Garcia Jerez, Jose Piña Flores                                          !
!                                                                                                                  !
!       Support: A.G.J. (agarcia-jerez@ual.es)                                                                     !
!                                                                                                                  !
!       Download compatible inversion software (HV-Inv) and updated versions of HV-DFA at                          !
!                                                                                                                  !
!                                http://www.ual.es/GruposInv/hv-inv/                                               !
!                                                                                                                  !
!                                                                                                                  !                                                                                                               !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                                                  !
!      Compilation: gfortran HV.f90 -o HV.exe [-static] [-fopenmp]    see readme.txt                               !
!      Some other compilers are supported, requiring modification of MISCELLANEA section (see a few lines below)   !
!                                                                                                                  !
!      Usage (example): HV.exe -f model.txt -nf 100 -fmin 0.1 -fmax 10 -logsam -nmr 20 -nml 20 -ph -hv  > HV.dat   !
!                                                                                                                  !
!      Type 'HV.exe -h' for help                                                                                   !
!                                                                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! REQUIRED MODULES AND GENERAL PROCEDURES
      INCLUDE 'modules.f90'                  ! Three modules, including common variables and definitions of KIND of variables      
      INCLUDE 'aux_procedures.f90'           ! Related with matrix manipulation and forward solving SH and PSV plane wave propagation 

      ! PROCEDURES REQUIRED FOR COMPUTATION OF DISPERSION CURVES
      INCLUDE 'Dispersion.f90'                ! Root search procedure (translated into FORTRAN and adapted from Geopsy)
      INCLUDE 'RootSolverTemplate.f90'        ! Aux. procedures of the root solver (translated into FORTRAN and adapted from Geopsy)  
      INCLUDE 'y.f90'                         ! Frequency equation, which vanishes on the dispersion curves (new)
              
      ! PROCEDURES FOR SW PART OF THE GREEN'S FUNCTIONS (COMPUTATION OF DISPERSION CURVES EXCLUDED)       
      INCLUDE 'GR_new.f90'
      INCLUDE 'GL.f90'
      INCLUDE 'WangMethod_Rayleigh_new.f90'
      INCLUDE 'WangMethod_Love.f90'
      
      ! PROCEDURES FOR BW PART OF THE GREEN'S FUNCTIONS
      INCLUDE 'BW_INTEGRALS_DAMPED.f90'
            
      ! MISCELLANEA
      !Load one of these files, depending on your compiler:
      !INCLUDE 'read_command_line.f90'       ! NON STANDARD F90. It runs under gfortran, g95?
      !INCLUDE 'read_command_line_MDSv5.f90' ! NON STANDARD F90. It runs under MS Developer Studio
      INCLUDE 'read_command_line.f03'        ! STANDARD F03.     It runs under gfortran, ftn95, and probably any fortran-2003 compiler


      !!!!!!!!!!!!!!!!!!!!!!!!   MAIN PROGRAM   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      PROGRAM HV           
      USE TYPES !defines long_float and long_cmplx
      USE Marc,ONLY:VALUES_L,VALUES_R,VALID_L,VALID_R,G_NX,X,ALFA,BTA,H,RHO,MU,NCAPAS,G_NMODES,ISRAYLEIGH,G_DX,G_PRECISION
      USE Globales,ONLY:UNIT_RDPGS,PI
      IMPLICIT NONE
      LOGICAL,PARAMETER::NORMALIZE=.FALSE.! Enable (True) /Disable (False) normalization of model parameter
      INTEGER INDX,INDXX ! General purpouse counters
      INTEGER NMODESRL(2),NKS ! Requested number of modes (read from the command line) and number of points for numerical integration 
      INTEGER::NM_R=0,NM_L=0,NM_RL=0 ! Number of modes found (Raylegh, Love, max of both)     
      LOGICAL,POINTER::SLOWNESSVALID(:) ! Slowness found for a particular mode (pointer to some elements of VALID)  PASS TO MODULE?
      REAL,   POINTER::SLOWNESSVALUES(:)! Slowness found for a particular mode (pointer to some elements of VALUES) PASS TO MODULE? 
      LOGICAL SALIDA ! true while everithing is going ok with the dispersion curves
      REAL PREC ! Precision in % of slowness read from the command line
      REAL(LONG_FLOAT),ALLOCATABLE::G1(:,:),IMG2(:,:),IMG3(:,:)! Outputs
      REAL(LONG_FLOAT),ALLOCATABLE::IMG11_pihalf(:),IMG33(:),IMVV(:),IMHPSV(:),IMHSH(:) !
      REAL(LONG_FLOAT)NORMH,NORMV,NORMG,NORMD,NORMM ! Used for normalizing
      REAL NORMT                                    ! Used for normalizing
      INTEGER,ALLOCATABLE::OFFSET_R(:),OFFSET_L(:) ! Number of frequencies below the cutoff frequency of each mode 
      INTERFACE DISPERSION
        LOGICAL FUNCTION DISPERSION(VALUES,VALID) RESULT(RETORNO)
        REAL,INTENT(INOUT),DIMENSION(:),TARGET::VALUES
        LOGICAL,INTENT(INOUT),DIMENSION(:),TARGET:: VALID
        END FUNCTION        
      END INTERFACE
      
      ! Read from the command line   
      CALL READ_COMMAND_LINE(NMODESRL,PREC,NKS)! Also updates: NCAPAS,G_NX,X,ALFA,BTA,H,RHO,MU,UNIT_RDPGS,SHDAMP,PSVDAMP  and allocates: G_SLOWS,G_SLOWP

      !Normalizing (optional)
      IF(NORMALIZE)THEN ! Normalize
        IF(NCAPAS>1)THEN
          NORMT=SUM(H/BTA(1:NCAPAS-1))/NCAPAS
          NORMH=SUM(H)/NCAPAS        
          NORMV=NORMH/NORMT
        ELSE
          NORMT=1 ! 1Hz is the reference frequency
          NORMV=BTA(1)
          NORMH=BTA(1)! Wavelength at 1Hz
        ENDIF
        NORMD=SUM(RHO)/NCAPAS
        NORMG=1/(NORMD*NORMV*NORMV*NORMH)! = NORMT^2/(NORMD*NORMH**3); T^2 M^-1
        NORMM=NORMV*NORMV*NORMD
        ! Do Normalization
        ALFA=ALFA/NORMV
        BTA=BTA/NORMV
        H=H/NORMH
        RHO=RHO/NORMD
        MU=MU/NORMM
        X=X*NORMT
      ELSE
        NORMT=1.D0;NORMH=1.D0;NORMV=1.D0;NORMD=1.D0;NORMG=1.D0;NORMT=1.D0;
      ENDIF
      
      !Open formatted error report      
      IF(UNIT_RDPGS(1)>6)THEN ! -1 = NO REPORT; 6 = TO STD OUTPUT
        OPEN(UNIT_RDPGS(1),file='report.dat',form='formatted')
        REWIND UNIT_RDPGS(1)
      ENDIF
      
      ! Compute some limits for the SW velocities from the model
      CALL SET_MODEL_PARAMETERS ! sets G_VELOCITYINVERSION,G_SLOWSMIN,G_SLOWSMAX,HALFSPACE_RAYLEIGH

      !Compute Rayleigh curves     
      IF(UNIT_RDPGS(5)>0.OR.(UNIT_RDPGS(3)>0.AND.NMODESRL(1)>0).OR.(UNIT_RDPGS(4)>0.AND.NMODESRL(1)>0))THEN
        G_NMODES=NMODESRL(1)
        G_PRECISION=PREC*1.E-2;
        G_DX=0.1
        ! allocate resulting dispersion curves and quality control vectors        
        ALLOCATE(VALUES_R(G_NX*G_NMODES),VALID_R(G_NX*G_NMODES))             
        VALUES_R=0;VALID_R=.FALSE. !Reset
        ISRAYLEIGH=.TRUE.
        SALIDA=DISPERSION(VALUES_R,VALID_R)
        ! Finding out the cutoff frequencies and the actual number of modes  
        ! Rayleigh
        ALLOCATE(OFFSET_R(G_NMODES))
        OFFSET_R=-1 !Default for not existing mode
        NM_R=G_NMODES ! Default: all the tried modes are found
        INDX=1
        DO INDXX=1,G_NMODES
          SLOWNESSVALUES=>VALUES_R((INDXX-1)*G_NX+1:INDXX*G_NX)
          SLOWNESSVALID=>VALID_R((INDXX-1)*G_NX+1:INDXX*G_NX)
          DO WHILE((.NOT.SLOWNESSVALID(INDX)).AND.(INDX<G_NX))  
            INDX=INDX+1
          ENDDO
          IF(SLOWNESSVALID(INDX))THEN          
            OFFSET_R(INDXX)=INDX-1
          ELSE
            NM_R=INDXX-1
            EXIT
          ENDIF        
        ENDDO
        IF(UNIT_RDPGS(3)>0)THEN
          IF(UNIT_RDPGS(3)>6)THEN
            ! save Rayleigh curves to file (ascii)
            OPEN(UNIT_RDPGS(3),file='Rph.dat',form='formatted')
            REWIND UNIT_RDPGS(3)
          ENDIF
          WRITE(UNIT_RDPGS(3),*)G_NX,G_NMODES
          WRITE(UNIT_RDPGS(3),*)(VALUES_R(INDXX)/NORMV,INDXX=1,G_NX*G_NMODES) ! Slowness
          WRITE(UNIT_RDPGS(3),*)(VALID_R(INDXX) ,INDXX=1,G_NX*G_NMODES)
          IF(UNIT_RDPGS(3)>6)CLOSE(UNIT_RDPGS(3))
        ENDIF
      ELSE
        SALIDA=.TRUE. ! If no Rayleigh mode is required, there is no problem with Rayleigh modes.   
      ENDIF

      ! Compute Love curves
      IF((UNIT_RDPGS(5)>0.AND.NMODESRL(2)>0).OR.(UNIT_RDPGS(3)>0.AND.NMODESRL(2)>0).OR.(UNIT_RDPGS(4)>0.AND.NMODESRL(2)>0))THEN           
        G_NMODES=NMODESRL(2)
        G_PRECISION=PREC*1.E-2
        G_DX=0.1
        ! ALLOCATE RESULTING DISPERSION CURVES AND QUALITY CONTROL VECTORS
        ALLOCATE(VALUES_L(G_NX*G_NMODES),VALID_L(G_NX*G_NMODES))      
        VALUES_L=0;VALID_L=.FALSE. !Reset
        ISRAYLEIGH=.FALSE.
        !SALIDA=SALIDA .AND. DISPERSION(VALUES_L,VALID_L)! Written in this order, Love wave are not calculated if Rayleigh waves fail
        SALIDA=DISPERSION(VALUES_L,VALID_L) .AND. SALIDA! Written in this order, Love wave are calculated even though Rayleigh waves fail        
        ! Actual number of Love modes and their cutoff frequencies:    
        ALLOCATE(OFFSET_L(G_NMODES))
        OFFSET_L=-1 !Default for not existing mode
        NM_L=G_NMODES ! Default: all the tried modes are found
        INDX=1
        DO INDXX=1,G_NMODES
          SLOWNESSVALUES=>VALUES_L((INDXX-1)*G_NX+1:INDXX*G_NX)
          SLOWNESSVALID=>VALID_L((INDXX-1)*G_NX+1:INDXX*G_NX)
          DO WHILE((.NOT.SLOWNESSVALID(INDX)).AND.(INDX<G_NX))  
            INDX=INDX+1
          ENDDO
          IF(SLOWNESSVALID(INDX))THEN          
            OFFSET_L(INDXX)=INDX-1
          ELSE
            NM_L=INDXX-1
            EXIT
          ENDIF        
        ENDDO              
        IF(UNIT_RDPGS(3)>0)THEN
          IF(UNIT_RDPGS(3)>6)THEN
            ! save Love curves to file (ascii)
            OPEN(UNIT_RDPGS(3),file='Lph.dat',form='formatted')
            REWIND UNIT_RDPGS(3)
          ENDIF 
          WRITE(UNIT_RDPGS(3),*)G_NX,G_NMODES
          WRITE(UNIT_RDPGS(3),*)(VALUES_L(INDXX)/NORMV,INDXX=1,G_NX*G_NMODES)
          WRITE(UNIT_RDPGS(3),*)(VALID_L(INDXX) ,INDXX=1,G_NX*G_NMODES)
          IF(UNIT_RDPGS(3)>6)CLOSE(UNIT_RDPGS(3))
        ENDIF
      ENDIF
      IF(.NOT.SALIDA.AND.UNIT_RDPGS(1)>0)THEN
        WRITE(UNIT_RDPGS(1),*)'Problems with dispersion curves calculations'
        CLOSE(UNIT_RDPGS(1))
      ENDIF

      ! Prepare calculation of SW Green functions
      NM_RL=MAX(NM_R,NM_L) 
      ALLOCATE(G1(G_NX,NM_RL),IMG2(G_NX,NM_RL),IMG3(G_NX,NM_RL))
      ALLOCATE(IMG11_pihalf(G_NX),IMG33(G_NX))      
      G1=0.;IMG2=0.;IMG3=0.;
      IMG11_pihalf=0;IMG33=0;  
      
      !Call the calculation of the SW part of Green's functions
      IF(UNIT_RDPGS(5)>0.OR.(UNIT_RDPGS(4)>0.AND.NM_R>0))THEN        
        DO INDXX=1,NM_R
          CALL GR(G1(:,INDXX),IMG3(:,INDXX),INDXX,OFFSET_R(INDXX)) 
          IMG11_pihalf(OFFSET_R(INDXX)+1:G_NX)=IMG11_pihalf(OFFSET_R(INDXX)+1:G_NX)+0.5*REAL(G1(OFFSET_R(INDXX)+1:G_NX,INDXX)) 
          IMG33(OFFSET_R(INDXX)+1:G_NX)=IMG33(OFFSET_R(INDXX)+1:G_NX)+IMG3(OFFSET_R(INDXX)+1:G_NX,INDXX)
        ENDDO
      ENDIF
      IF(UNIT_RDPGS(5)>0.OR.(UNIT_RDPGS(4)>0.AND.NM_L>0))THEN        
        DO INDXX=1,NM_L
          CALL GL(IMG2(:,INDXX),INDXX,OFFSET_L(INDXX))
          IMG11_pihalf(OFFSET_L(INDXX)+1:G_NX)=IMG11_pihalf(OFFSET_L(INDXX)+1:G_NX)-0.5*IMG2(OFFSET_L(INDXX)+1:G_NX,INDXX)  
        ENDDO
      ENDIF
      
      ! Calculation of BW integrals
      ALLOCATE(IMVV(G_NX),IMHPSV(G_NX),IMHSH(G_NX))      
      IF(NKS>0)THEN    
        !$OMP PARALLEL DO        
        DO INDXX=1,G_NX
          CALL BWR(IMVV(INDXX),IMHPSV(INDXX),IMHSH(INDXX), NKS,X(INDXX))
        ENDDO
        !$OMP END PARALLEL DO             
      ELSE
        IMVV=0.
        IMHPSV=0.
        IMHSH=0.
      ENDIF

      !Print out group velocities 
      IF(UNIT_RDPGS(4)>0.AND.NM_R>0)THEN
        IF(UNIT_RDPGS(4)>6)THEN
          OPEN(UNIT_RDPGS(4),file='Rgr.dat',form='formatted')
          REWIND UNIT_RDPGS(4)
        ENDIF
        WRITE(UNIT_RDPGS(4),*)G_NX,NM_R
        WRITE(UNIT_RDPGS(4),*)(VALUES_R(INDXX)/NORMV,INDXX=1,G_NX*NM_R)
        WRITE(UNIT_RDPGS(4),*)(VALID_R(INDXX) ,INDXX=1,G_NX*NM_R)
        IF(UNIT_RDPGS(4)>6)CLOSE(UNIT_RDPGS(4))
      ENDIF
      IF(UNIT_RDPGS(4)>0.AND.NM_L>0)THEN
        IF(UNIT_RDPGS(4)>6)THEN
          OPEN(UNIT_RDPGS(4),file='Lgr.dat',form='formatted')
          REWIND UNIT_RDPGS(4)
        ENDIF
        WRITE(UNIT_RDPGS(4),*)G_NX,NM_L
        WRITE(UNIT_RDPGS(4),*)(VALUES_L(INDXX)/NORMV,INDXX=1,G_NX*NM_L)
        WRITE(UNIT_RDPGS(4),*)(VALID_L(INDXX) ,INDXX=1,G_NX*NM_L)
        IF(UNIT_RDPGS(4)>6)CLOSE(UNIT_RDPGS(4))
      ENDIF
                    
      !Print out modal contributions
      IF(UNIT_RDPGS(2)>0)THEN
        IF(UNIT_RDPGS(2)>6)OPEN(UNIT_RDPGS(2),file='G123.dat',form='formatted')
        DO INDXX=1,NM_RL
          WRITE(UNIT_RDPGS(2),*)'      FREQ.           ReG1             ImG2             ImG3'
          DO INDX=1,G_NX
            WRITE(UNIT_RDPGS(2),*)X(INDX)/(2*PI*NORMT),0.5*G1(INDX,INDXX)*NORMG,-0.5*IMG2(INDX,INDXX)*NORMG,IMG3(INDX,INDXX)*NORMG
          ENDDO
        ENDDO
        IF(UNIT_RDPGS(2)>6)THEN
          CLOSE(UNIT_RDPGS(2))
          OPEN(UNIT_RDPGS(2),file='GBW.dat',form='formatted')
        ENDIF
        WRITE(UNIT_RDPGS(2),*)'         FREQ.                   IMHPSV                  IMHSH                   IMVPSV'
        DO INDX=1,G_NX
          WRITE(UNIT_RDPGS(2),*)X(INDX)/(2*PI*NORMT),IMHPSV(INDX)*NORMG,IMHSH(INDX)*NORMG,IMVV(INDX)*NORMG
        ENDDO
        IF(UNIT_RDPGS(2)>6)CLOSE(UNIT_RDPGS(2))        
      ENDIF
      
      ! H/V output
      IF(UNIT_RDPGS(5)>0)THEN
        IF(UNIT_RDPGS(5)>6)OPEN(UNIT_RDPGS(5),file='HV.dat',form='formatted')
        WRITE(UNIT_RDPGS(5),*)(X(INDX)/(2.D0*PI*NORMT),&
                               SQRT(2.D0*(IMG11_pihalf(INDX)+IMHPSV(INDX)+IMHSH(INDX))/(IMG33(INDX)+IMVV(INDX))),INDX=1,G_NX)                   
        IF(UNIT_RDPGS(5)>6)CLOSE(UNIT_RDPGS(5))
      ENDIF
      
      END PROGRAM HV

      SUBROUTINE SET_MODEL_PARAMETERS()
      USE Marc,ONLY:NCAPAS,G_MAXRAYLEIGHSLOWNESS,G_SLOWS,G_SLOWP,ALFA,BTA,G_VELOCITYINVERSION,G_SLOWSMIN,G_SLOWSMAX,&
                    HALFSPACE_RAYLEIGH
      IMPLICIT NONE
      INTEGER I,IMAXSLOWS
      G_SLOWS=1/REAL(BTA)
      G_SLOWP=1/REAL(ALFA)
      IMAXSLOWS=MAXLOC(G_SLOWS,1)
      G_MAXRAYLEIGHSLOWNESS=REAL(HALFSPACE_RAYLEIGH(IMAXSLOWS))
      G_SLOWSMAX=G_SLOWS(IMAXSLOWS)
      G_SLOWSMIN=G_SLOWS(NCAPAS) ! MINIMUM SLOWNESS FOR THE HIGHER MODES
      DO I=2,NCAPAS
        IF(G_SLOWS(I)>G_SLOWS(I-1).OR.G_SLOWP(I)>G_SLOWP(I-1))THEN
          G_VELOCITYINVERSION=I;RETURN; ! -1,2,3,4...NCAPAS
        ENDIF
      ENDDO
      G_VELOCITYINVERSION= -1;RETURN;
      END SUBROUTINE
     
 