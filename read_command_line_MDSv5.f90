      
      ! This subroutine uses old GNU methods for reading the command line instead of
      ! the procedures of the FORTRAN 2003 standard (used in read_command_line.F95)
      !  
      ! IARGC  replaces COMMAND_ARGUMENT_COUNT
      ! GETARG replaces GET_COMMAND_ARGUMENT
      !
      ! It also loads DFPORT and DFLIB modules from Microsoft Developement Studio

      SUBROUTINE READ_COMMAND_LINE(NMODESRL,PREC,NKS)
      USE DFPORT
	  USE DFLIB 
      USE Marc,ONLY:X,NCAPAS,G_NX,ALFA,BTA,H,RHO,MU,G_SLOWS,G_SLOWP
      USE Globales,ONLY:PI,TWOPI,UNIT_RDPGS,SHDAMP,PSVDAMP
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::NMODESRL(2) ! Number of modes indicated in the command line 
      REAL,INTENT(OUT)::PREC ! PRECISION IN %      
      REAL FMIN,FMAX,DISCARD
      INTEGER(2)INDXX
      INTEGER SAMTYPE,INDX,IOSS,LENGTH,STATUS,NARGMNS
      CHARACTER(120)FLAG,VALOR,FREQFILE,MODLFILE !WARNING 
      LOGICAL ISFREQFILE,ISMODLFILE!,SHOW_REPORT
      INTEGER TO_STD,NKS                               
      ! Initializing first to invalid or default values 
      FMIN=-1;FMAX=-1;G_NX=-1;NCAPAS=-1;NMODESRL=0
      SAMTYPE=0;ISFREQFILE=.FALSE.;
      PREC=1.E-4;NKS=0
      SHDAMP=1.0d-5;PSVDAMP=1.0d-5
      ! Starting command line analysis
      NARGMNS=IARGC() !COMMAND_ARGUMENT_COUNT()
      IF(NARGMNS==0)THEN
        PRINT*,'Insufficient arguments. Use flag -h for help'
        STOP
      ENDIF      
      INDXX=0      
      DO WHILE(INDXX<NARGMNS)
        INDXX=INDXX+1
        CALL GETARG(INDXX,FLAG);LENGTH=LEN_TRIM(FLAG)
        ! CALL LCASE@(FLAG)! no va con gfortran        
        IF(FLAG(1:1)=='-')THEN
          SELECT CASE(FLAG(2:LENGTH))
                       
!         To build the vector of frequencies           
          CASE('fmin')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)FMIN  
          CASE('fmax')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)FMAX
          CASE('nf')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)G_NX
            ALLOCATE(X(G_NX))! allocate circular frequencies vector
          CASE('logsam')!CASE('sam')
!            INDXX=INDXX+1
!            CALL GETARG(INDXX,VALOR)
!            READ(VALOR,*)SAMTYPE
            SAMTYPE=1
          CASE('ff')
            INDXX=INDXX+1
            CALL GETARG(INDXX,FREQFILE)
            ISFREQFILE=.TRUE.  

!         Controlling dispersion curves                       
          CASE('nmr')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)NMODESRL(1)
          CASE('nml')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)NMODESRL(2)            
          CASE('prec')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)PREC

!         Controlling body wave integrals
          CASE('nks')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)NKS            
          CASE('ash')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)SHDAMP
          CASE('apsv')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)PSVDAMP
                        
!         Reading and building the model from the command line            
          CASE('nlyr')
            INDXX=INDXX+1
            CALL GETARG(INDXX,VALOR)
            READ(VALOR,*)NCAPAS
            ALLOCATE(ALFA(NCAPAS),BTA(NCAPAS),G_SLOWS(NCAPAS),G_SLOWP(NCAPAS),RHO(NCAPAS),H(NCAPAS-1))               
          CASE('vs')
            IF(NCAPAS<0)THEN
              WRITE(*,*)'flag -nlyr required before listing model'
              STOP
            ELSE
              DO INDX=1,NCAPAS
                INDXX=INDXX+1
                CALL GETARG(INDXX,VALOR)
                READ(VALOR,*)BTA(INDX)
              ENDDO
            ENDIF
          CASE('vp')
            IF(NCAPAS<0)THEN
              WRITE(*,*)'flag -nlyr required before listing model'
              STOP
            ELSE
              DO INDX=1,NCAPAS
                INDXX=INDXX+1
                CALL GETARG(INDXX,VALOR)
                READ(VALOR,*)ALFA(INDX)
              ENDDO
            ENDIF
          CASE('dens')
            IF(NCAPAS<0)THEN
              WRITE(*,*)'flag -nlyr required before listing model'
              STOP
            ELSE
              DO INDX=1,NCAPAS
                INDXX=INDXX+1
                CALL GETARG(INDXX,VALOR)
                READ(VALOR,*)RHO(INDX)
              ENDDO              
            ENDIF
          CASE('thk')
            IF(NCAPAS<0)THEN
              WRITE(*,*)'flag -nlyr required before listing model'
              STOP
            ELSE
              DO INDX=1,NCAPAS-1
                INDXX=INDXX+1
                CALL GETARG(INDXX,VALOR)
                READ(VALOR,*)H(INDX)
              ENDDO              
            ENDIF

!         Reading and building the model from model file
          CASE('f')
            INDXX=INDXX+1
            CALL GETARG(INDXX,MODLFILE) 
            INQUIRE(FILE=MODLFILE,EXIST=ISMODLFILE)
            IF(ISMODLFILE)THEN            
              OPEN(3,file=MODLFILE,form='formatted')
              REWIND 3
              READ(3,*)NCAPAS
              ALLOCATE(ALFA(NCAPAS),BTA(NCAPAS),G_SLOWS(NCAPAS),G_SLOWP(NCAPAS),RHO(NCAPAS),H(NCAPAS-1),MU(NCAPAS))                      
              DO INDX=1,NCAPAS-1
                READ(3,*) H(INDX),ALFA(INDX),BTA(INDX),RHO(INDX)
              ENDDO
              READ(3,*)DISCARD,ALFA(NCAPAS),BTA(NCAPAS),RHO(NCAPAS)
              CLOSE(3)
            ELSE
              WRITE(*,*)'missing model file'
            ENDIF
          CASE('rep')
            UNIT_RDPGS(1)=8
            UNIT_RDPGS(2)=9            
            TO_STD=1 ! provisionally identifies last output flag read in the commandline            
          CASE('ph')
            UNIT_RDPGS(3)=10            
            TO_STD=3
          CASE('gr')
            UNIT_RDPGS(4)=11            
            TO_STD=4            
          CASE('hv')
            UNIT_RDPGS(5)=12                        
            TO_STD=5            

!         Listing help   
          CASE('h')
            PRINT*
            PRINT*
            PRINT*,'                                 HV-DFA    version 1.0                                                     '
            PRINT*
            PRINT*,'   A software to compute HV under the DFA and the contributions to Im[Gij] at source                       '
            PRINT*
            PRINT*,'This is an implementation of the Sanchez-Sesma et al (2011, Geoph. J. Int.) theory for H/V                 '
            PRINT*,'using a propagator matrix scheme based on Wang (1999, Bull, Seism. Soc. Am.)and contour                    '
            PRINT*,'integration in the complex-wavenumber plane.                                                               '
            PRINT*
            PRINT*,'Root-search algorithm translated into FORTRAN 90 from M. Wathelet software (Geopsy)                        '                       
            PRINT*
            PRINT*,'                               HV-INV Project Team                                                         '      
            PRINT*                                                                                                                                               
            PRINT*,'             Antonio Garcia Jerez                 Jose Pina Flores                                         '
            PRINT*,'             University of Almeria                Instituto de Ingenieria-UNAM                             '
            PRINT*,'             agj574@ual.es                                                                                 '
            PRINT*
            PRINT*,'Download compatible inversion software (HV-INV) and updated versions of HV-DFA at                          '
            PRINT*
            PRINT*,'                        http://www.ual.es/GruposInv/hv-inv/                                                '
            PRINT*
            PRINT*
            PRINT* 
            PRINT*,'Usage: HV [OPTIONS]'
            PRINT*
            PRINT*,'OPTIONS FOR CONTROLING HV:'
            PRINT*
            PRINT*,'-fmin X : Minimum frequency'  
            PRINT*,'-fmax X : Maximum frequency'          
            PRINT*,'-nf N   : Number of frequencies'
            PRINT*,'-nmr  N : Max. of Rayleigh modes to be considered'
            PRINT*,'-nml  N : Max. of Love modes to be considered'            
            PRINT*,'-prec X : Rel. precission in slowness (default: 1E-4 per cent)'            
            PRINT*,'-nks  N : Number of k values for numeric integrals'
            PRINT*,'          Use 0 or missing flag to skip BW calculation'
            PRINT*,'-apsv X : Attenuation to stabilize PSV. w -> w-I*apsv*w (DEFAULT = 0)'
            PRINT*,'-ash  X : Attenuation to stabilize SH.  w -> w-I*ash*w  (DEFAULT = 0)'            
            PRINT*,'-logsam  : Regular sampling in log[frequency] (default is regular in frequency)'
            PRINT*,'-ff <file> frequency list (sorted) from file'
            PRINT*
            PRINT*,'INPUT MODEL PROPERTIES:'
            PRINT*
            PRINT*,'-nlyr N : Number of layer (including halfspace)'            
            PRINT*,'-dens Density'            
            PRINT*,'-thk  Layer thicknesses'
            PRINT*,'-vp   P wave velocity'
            PRINT*,'-vs   S wave velocity'
            PRINT*,'-f <file> Model from file (see below for format)'
            PRINT*
            PRINT*,'OUTPUTS SELECTION:'
            PRINT*
            PRINT*,'-hv     : Outputs ( freq., H/V) pairs -> HV.dat'
            PRINT*,'-ph     : Outputs phase slowness -> Rph.dat,Lph.dat containing:'
            PRINT*,'          Number of frequencies (N), number of modes (NM), list of N x NM'
            PRINT*,'          slowness values for increasing frequencies and mode number, list of'
            PRINT*,'          N x NM logical values. False (F) and 0 slowness indicate that the'
            PRINT*,'          mode does not exist at that frequency or could not be calculated.' 
            PRINT*,'          The inner loop runs on the frequency.'           
            PRINT*,'-gr     : Outputs group velocities -> Rgr.dat,Lgr.dat. Same file structure.'
            PRINT*,'-rep    : Writes a report if computation fails and details of calculations ->'
            PRINT*,'        -> report.dat,GBW.dat (body waves),G123.dat (surface waves),ERR.DAT '
            PRINT*,'The last flag listed of these four will write into the std output instead of'
            PRINT*,'using the specified file'
            PRINT*
            PRINT*,'FORMAT FOR LAYERED MODEL FILES :'
            PRINT*
            PRINT*,'Line 1    <number of layers including halfspace>'
            PRINT*,'Line 2    <thickness(m)> <Vp(m / s)> <Vs(m / s)> <Density(kg / m3)>'
            PRINT*,'...'
            PRINT*,'Line N    0 <Vp(m / s)> <Vs(m / s)> <Density(kg / m3)>'
            PRINT*,'********************************************************************'           
            PRINT*            
            STOP           
          CASE DEFAULT
            WRITE(*,*)'Unknown flag ',FLAG
            STOP
          END SELECT
        ELSE
          WRITE(*,*)'Flag expected, ',FLAG,' found'
          STOP
        ENDIF  
      ENDDO
      
      ! CONTROL AND FINAL TASKS
      ! What variable does go to the std output?
      IF(ALL(UNIT_RDPGS<0))THEN
        PRINT*,'Use -ph -gr and/or -hv to select an output'
        STOP
      ELSE
        UNIT_RDPGS(TO_STD)=6
      ENDIF
      ! Calculation of Mu
      IF(NCAPAS>0)THEN
        mu=bta**2.d0*rho
      ELSE
        PRINT*,'Error: Undefined ground model';STOP
      ENDIF
      ! obtain circular frequency vector
      IF(ISFREQFILE)THEN
        ! READ FREQUENCY FILE
        INQUIRE(FILE=FREQFILE,EXIST=ISFREQFILE)
        IF(ISFREQFILE)THEN
!         print*,'vamos a cargar archivo frecuencias'
          OPEN(3,file=FREQFILE,form='formatted')
          REWIND 3
          IF(G_NX==-1)THEN ! X has not been allocated yet
            G_NX=0
            READ(3,*,IOSTAT=IOSS)            
            DO WHILE(IOSS==0)
              G_NX=G_NX+1
              READ(3,*,IOSTAT=IOSS)
            ENDDO
            REWIND 3
            IF(G_NX>0)THEN 
              ALLOCATE(X(G_NX))
            ELSE
              WRITE(*,*)'Error in frequency file';STOP
            ENDIF
          ENDIF          
          IOSS=0
          INDX=0
          DO WHILE(IOSS==0.AND.INDX<G_NX)
            INDX=INDX+1 
            READ(3,*,IOSTAT=IOSS)X(INDX)
          ENDDO
          IF(INDX<G_NX)G_NX=INDX-1          
          CLOSE(3)
          !WRITE(*,*)'volcado lista frecuencias'
          !WRITE(*,*)(INDXX,X(INDXX),INDXX=1,G_NX)
          !WRITE(*,*)'fin, G_NX=',G_NX
          X=TWOPI*X
        ELSE
          WRITE(*,*)'Missing frequency file'
          STOP
        ENDIF
      ELSE
        ! BUILD UP X VECTOR
        IF(FMIN<0.OR.FMAX<0.OR.FMAX<FMIN.OR.G_NX<0)THEN 
          WRITE(*,*)'Incorrect or insufficient input values. Use flag -h for help';STOP
        ENDIF
        IF(G_NX==1)THEN
          X=2*PI*FMIN   
        ELSEIF(SAMTYPE==0)THEN ! Constant increment      
          DO INDXX=1,G_NX
            X(INDXX)= (FMIN + (INDXX-1)*(FMAX-FMIN)/(G_NX-1))*2*PI
          ENDDO
        ELSEIF(SAMTYPE==1)THEN ! Logarithmic increment  
          FMIN=LOG10(FMIN*2*PI)
          FMAX=LOG10(FMAX*2*PI)
          DO INDXX=1,G_NX
            X(INDXX)=10**(FMIN+(INDXX-1)*(FMAX-FMIN)/(G_NX-1))
          ENDDO
        ELSE
          WRITE(*,*)'Invalid type of freq. increment (flag -sam)';STOP                    
        ENDIF
      ENDIF
      
      IF(ALL(NMODESRL.LE.0).AND.NKS.LE.0)THEN
        PRINT*,'Nothing to compute. Check -nmr -nml -nks flags';STOP
      ELSEIF(UNIT_RDPGS(5).LE.0.AND.NKS.GT.0)THEN
        IF (ALL(NMODESRL.LE.0))THEN
           PRINT*,'Error: BW will be calculated but there is nothing to show (-hv flag missing)';STOP
        ELSE
           PRINT*,'Warning: BW will be calculated but not shown (-hv flag missing or remove -nks)'
        ENDIF     
      ELSEIF(NCAPAS.LE.0)THEN
        PRINT*,'Number of layers has to be positive integer';STOP            
      ENDIF
      END SUBROUTINE
      
