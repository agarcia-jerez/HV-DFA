
! This function is based on free code written in C by Marc Wathelet and included in
! Geopsy: www.geopsy.org

LOGICAL FUNCTION DISPERSION(VALUES,VALID) RESULT(RETORNO)
USE Marc
USE globales,ONLY:UNIT_RDPGS
IMPLICIT NONE
REAL,INTENT(INOUT),DIMENSION(:),TARGET::VALUES
LOGICAL,INTENT(INOUT),DIMENSION(:),TARGET:: VALID
INTEGER IMODE,I,J ! MODE AND GENERAL-PURPOSE COUNTERS
REAL MINSLOW,MAXSLOW,TO,CURSLOW,NEWCURSLOW
REAL LASTSUPBOUND,NEWLASTSUPBOUND
LOGICAL::ERRORYETDETECTED
! FUNCTIONS
LOGICAL SEARCHUP,SEARCHDOWN
! POINTERS
REAL,ALLOCATABLE::CURMODE(:)
REAL,ALLOCATABLE::LASTMODE(:)
REAL,ALLOCATABLE::TMPMODE(:)
LOGICAL,POINTER::SLOWNESSVALID(:)
REAL,POINTER::SLOWNESSVALUES(:)
! ALLOCATING POINTERS
ALLOCATE(CURMODE(G_NX),LASTMODE(G_NX),TMPMODE(G_NX),SLOWNESSVALID(G_NX),SLOWNESSVALUES(G_NX))
LASTMODE=-9999. ! Initializing to absurd value
MINSLOW=G_SLOWSMIN
IF (ISRAYLEIGH)THEN
  LASTSUPBOUND=G_MAXRAYLEIGHSLOWNESS
ELSE
  LASTSUPBOUND=G_SLOWS(1)
ENDIF
FORBITVELOCITYINVERSION=(G_VELOCITYINVERSION==-1)
!  Initialize equation's polarity
G_OMEGA=0.05 ! PROVISIONAL very low value for calculation of G_POLARITY. The sign of polarity is the same at all frequencies
G_POLARITY=Y(LASTSUPBOUND)
! Some more initializations
G_X2=LASTSUPBOUND ! AGJ
!!!!!!!!!!!!!!!!!!!!!!!!!
DO IMODE=0,G_NMODES-1  !!
!!!!!!!!!!!!!!!!!!!!!!!!!
  IF(IMODE==0)THEN
    IF (ISRAYLEIGH)THEN
      MAXSLOW=G_MAXRAYLEIGHSLOWNESS
      G_DX=0.25*(1-G_SLOWS(1)/G_MAXRAYLEIGHSLOWNESS)
    ELSE
      MAXSLOW=G_SLOWSMAX
      G_MAXRAYLEIGHSLOWNESS=REAL(HALFSPACE_RAYLEIGH(MAXLOC(G_SLOWS,1)))
      G_DX=0.25*(1-G_SLOWS(1)/G_MAXRAYLEIGHSLOWNESS)! CHECK for Love waves
    ENDIF
  ELSE
    ! Exchange LASTMODE and CURMODE
    TMPMODE=LASTMODE
    LASTMODE=CURMODE
    CURMODE=TMPMODE
    MAXSLOW=lastMode(G_NX)
    G_DX=0.5*(G_MAXRAYLEIGHSLOWNESS-G_SLOWS(1))/MAXSLOW*(2*PI/x(G_NX)) !REVISED
    G_POLARITY=-G_POLARITY! CHECK for Love waves
  ENDIF
  SLOWNESSVALUES=>VALUES(IMODE*G_NX+1:(IMODE+1)*G_NX)
  SLOWNESSVALID=>  VALID(IMODE*G_NX+1:(IMODE+1)*G_NX)
!!!!!!!!!!!!!!!!!!!!!!!!!
  DO WHILE(.TRUE.)     !! RESTARTS THE INNER BUCLE ON THE FREQUENCIES AFTER ERRORS
!!!!!!!!!!!!!!!!!!!!!!!!!
    ERRORYETDETECTED=.FALSE.
    CURSLOW=MAXSLOW
!!!!!!!!!!!!!!!!!!!!!!!!!
    DO I=G_NX,1,-1   !!
!!!!!!!!!!!!!!!!!!!!!!!!!
      G_OMEGA=X(I)
      IF(iMode > 0)THEN
        maxSlow=lastMode(i)
        IF(curSlow > maxSlow)curSlow=maxSlow
      ENDIF
      IF(searchDown(curSlow, minSlow, maxSlow))THEN
        CALL neville;
        newCurSlow=MIN(G_X1,G_X2)
      ELSE
        IF(iMode==0)THEN
          !For fundamental mode this condition can never occur unless a mode jumping occured
          errorYetDetected=.true.
          IF(UNIT_RDPGS(1)>0)THEN
            WRITE(UNIT_RDPGS(1),*)"** Warning ** : mode jumping for mode",IMODE," (end), reducing step ratio to ",G_DX*0.1
          ENDIF
          EXIT
        ELSE
          !The end of the curve, the mode does not exist at this period (no root)
          DO j=i,1,-1
            SLOWNESSVAlID(J)=.FALSE.
            curMode(j)=minSlow
          ENDDO
          curSlow=minSlow;
          EXIT
        ENDIF
      ENDIF
      newlastSupBound=MAX(G_X1,G_X2)
      IF(i<G_NX.AND.newCurSlow-curSlow>G_precision)THEN
        IF(UNIT_RDPGS(1)>0)THEN
          WRITE(UNIT_RDPGS(1),*)'INVERSION DETECTED'
        ENDIF
        IF(forbitVelocityInversion)THEN
          errorYetDetected=.true. 
          IF(UNIT_RDPGS(1)>0)THEN
            WRITE(UNIT_RDPGS(1),*)"** Warning ** : retrograde dispersion not allowed (mode ",iMode,"), Rayleigh=",ISRAYLEIGH
            WRITE(UNIT_RDPGS(1),*)"At Frequency ",i," omega ",x(i),"freq ",x(i)/(2*pi),"s = ",newCurSlow
            WRITE(UNIT_RDPGS(1),*)"Reducing step ratio to ",G_DX*0.1
            ! Opening a binary diagnosis file with evaluations of y(s) at the current frequency
            IF(UNIT_RDPGS(2)>6)THEN
              OPEN(UNIT_RDPGS(  2),file='ERR.dat',form='unformatted')
              REWIND UNIT_RDPGS(2)
              WRITE(UNIT_RDPGS(2))G_OMEGA/(2*PI)
              WRITE(UNIT_RDPGS(2))G_SLOWSMIN,G_MAXRAYLEIGHSLOWNESS
              WRITE(UNIT_RDPGS(2))1000
              WRITE(UNIT_RDPGS(2))(Y(J*(G_MAXRAYLEIGHSLOWNESS-G_SLOWSMIN)/999+G_SLOWSMIN),J=0,999)
              CLOSE(UNIT_RDPGS(2))
              WRITE(UNIT_RDPGS(1),*)'DIAGNOSIS FILE WRITTEN'
            ENDIF
          ENDIF
          EXIT ! Exits from the frequency loop to the endless loop. Frequency loop will start again
               ! from the beginning. After this exit, G_DX will be decremented due to the setting errorYetDetected=.true.
        ELSE
          ! Inversions are allowed, check at last point that no mode existED at a higher slowness AT THE PREVIOUS FREQUENCY
          G_OMEGA=x(i+1)
          G_DX=G_DX*0.1
          IF(iMode>0)THEN
            TO=lastMode(i+1)
          ELSE
            TO=maxSlow
          ENDIF
          IF(searchUp(lastSupBound,TO))THEN
            !recalculating curve because of mode jumping
            errorYetDetected=.TRUE.
            G_DX=G_DX*10.
            EXIT ! DEJAMOS ANTICIPADAMENTE EL BUCLE DE LAS FRECUENCIAS
          ENDIF
          ! LA CURVA A LA FRECUENCIA PREVIA ESTABA BIEN, RESTORE G_DX:
          G_DX=G_DX*10.
        ENDIF
      ENDIF
      lastSupBound=newLastSupBound
      curSlow=newCurSlow
      curMode(i)=curSlow
      slownessValues(i)=curSlow
      slownessVALID(i)=.TRUE.
    ENDDO ! End of the the loop on the frequencies (index I)
    ! When calculating the dispersion curve, stepRatio is positive to search downwards
    ! (From maxSlow towards minSlow)
    ! Here we might want to search carefully upwards: from curSlow to maxSlow
    ! If at least one root is found, this proves that mode jumping occured
    IF(curSlow > minSlow)THEN
      !The curSlow was calculated using _solve and G_x2 is the value returned (slightly greater than the true root).
      !Here we need the G_x1 value, e.g. slightly less than the true root,in order to be sure that the root eventually
      !found will not be the same as curSlow
      curSlow=MAX(G_X1,G_X2)
    ENDIF
    IF(errorYetDetected.OR.searchUp(curSlow,maxSlow))THEN
      !recalculating curve because of mode jumping
      IF(G_DX<1E-4)THEN                                                     !REVISED
        ! G_DX IS SMALL ENOUGH. NO MORE REFINEMENTS.
        ! THIS MODE AND THE HIGHER MODES FAIL. ABORT
        IF(UNIT_RDPGS(1)>0)THEN
          IF(ISRAYLEIGH)THEN
            WRITE(UNIT_RDPGS(1),*)"Cancel Rayleigh at ",i,"-th Freq. ",x(i)/(2*pi),", Mode ",iMode
          ELSE
            WRITE(UNIT_RDPGS(1),*)"Cancel Love at ",i,"-th Freq. ",x(i)/(2*pi),", Mode ",iMode
          ENDIF
        ENDIF             
        I=IMODE
        DO WHILE(I<G_NMODES)
          SLOWNESSVALID=>VALID(I+1:I+G_NX)
          DO j=G_NX,1,-1
            SLOWNESSVALID(j)=.false.!CANCELL THIS MODE AND HIHGER ONES
          ENDDO
          I=I+1
        ENDDO 
        RETORNO=.FALSE.;RETURN
      ENDIF
      G_PRECISION=G_PRECISION*0.1
      G_DX=G_DX*0.1
      IF(iMode > 0)maxSlow=lastMode(G_NX);
    ELSE
      EXIT
    ENDIF
  ENDDO ! ENDLESS LOOP
ENDDO ! LOOP ON THE MODES
RETORNO=.TRUE.;RETURN;
END FUNCTION