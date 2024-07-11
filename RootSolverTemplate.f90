     
!     This file is an adaptation and translation into FORTRAN of subroutined used in
!     Geopsy (www.geopsy.org) for root bracketing    

      SUBROUTINE halving()
      USE Marc,ONLY:ROOTSOLVER_MAX_ITERATIONS,G_PRECISION,G_X1,G_X2,G_Y1,G_Y2
      USE GLOBALES,ONLY:UNIT_RDPGS
      IMPLICIT NONE
      Real diff,maxDiff,minDiff
      Real x3,y3,tmpAbs1
      INTEGER INDEX
      ! FUNCTIONS
      LOGICAL isFound
      REAL Y
!     First halving
      x3=0.5*(g_x1 + g_x2)
      y3=y(x3);
!     Relative Convergence criteria
      tmpAbs1=g_x1 + g_x2
      if(tmpAbs1 < 0)tmpAbs1=-tmpAbs1
      maxDiff=tmpAbs1 * G_precision
      minDiff=-maxDiff;
      DO index=1,ROOTSOLVER_MAX_ITERATIONS-1
        IF(y3==0.0)THEN
!         just move the X3 by one tenth of actuel range
          x3 = x3 -(g_x2 - g_x1) * 0.1
!         recalculate the Y3, it will always be different from 0
          y3=y(x3)
        ENDIF
        IF(isFound(g_y1, y3))THEN
          g_x2=x3
          g_y2=y3
        ELSE
          g_x1=x3
          g_y1=y3
        ENDIF
!       Check for convergence with a relative criteria, if so exit
        diff=g_x1-g_x2
        IF(minDiff <= diff .and. diff <= maxDiff)RETURN
        x3=0.5*(g_x1 + g_x2)
        y3=y(x3)
      ENDDO
      IF(UNIT_RDPGS(1)>0)WRITE(UNIT_RDPGS(1),*)'Halving maximum iteration reached'
      END SUBROUTINE
            
      LOGICAL FUNCTION ISFOUND(Y1,Y2)
!     True if range (Y1 Y2) contains zero
      REAL Y1,Y2
      ISFOUND=(Y1<0.0.AND.Y2>0.0).OR.(Y1>0.0.AND.Y2<0.0);
      RETURN
      END FUNCTION

      LOGICAL FUNCTION SEARCHDOWN(FROM,MINIMUM,MAXIMUM)
      USE Marc,ONLY:G_dxType,G_dx,GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::FROM,MINIMUM,MAXIMUM
      REAL DX
      ! FUNCTIONS
      LOGICAL SEARCH2
      GO_DOWN=.TRUE.
      IF(G_dxType)THEN !ABSOLUTE
        DX=G_DX
        SEARCHDOWN=SEARCH2(FROM, MINIMUM, MAXIMUM,DX)
        RETURN
      ELSE
        DX=1.0-G_DX
        SEARCHDOWN=SEARCH2(FROM, MINIMUM, MAXIMUM,DX)
        RETURN
      ENDIF
      END FUNCTION

      LOGICAL FUNCTION SEARCHUP(FROM, TOWARD)
      ! RETURNS FALSE IF FROM ALREADY EXCEEDS TOWARD (FROM>TOWARD) OR
      ! IF THE ROOT COULDN'T BE BRACKETED BETWEEN G_X1 AND G_X2
      USE Marc,ONLY:G_dxType,G_dx,GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::FROM,TOWARD
      REAL DX
      ! FUNCTIONS
      LOGICAL SEARCH1
      !STARTING SEARCHUP
      GO_DOWN=.FALSE.
      IF(from>TOWARD)THEN
        SEARCHUP=.FALSE.;return
      ENDIF
      IF(G_dxType)THEN!==Absolute
        DX=G_DX
        SEARCHUP=SEARCH1(FROM,TOWARD,DX);RETURN
      ELSE
        DX=1.0+G_DX
        SEARCHUP=SEARCH1(FROM,TOWARD,DX);return
      ENDIF
      END FUNCTION

      LOGICAL FUNCTION SEARCHUP2(FROM,MINIMUM,MAXIMUM)
      USE Marc,ONLY:G_dxType,G_dx,GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::FROM,MINIMUM,MAXIMUM
      REAL DX
      ! FUNCTIONS
      LOGICAL SEARCH2
      GO_DOWN=.FALSE.
      IF(G_dxType)THEN!ABSOLUTE
        DX=G_DX
        SEARCHUP2=SEARCH2(FROM, MINIMUM, MAXIMUM,DX)
        RETURN
      ELSE
        DX=1.0+G_DX
        SEARCHUP2=SEARCH2(FROM, MINIMUM, MAXIMUM,DX)
        RETURN
      ENDIF
      END FUNCTION

      REAL FUNCTION IT_NEXT(X,DX)
      USE Marc,ONLY: G_DXTYPE,GO_DOWN
      REAL,INTENT(IN)::DX
      IF (G_dxType)THEN !ABSOLUTE
        IF (GO_DOWN)THEN
          IT_NEXT=X-DX
          RETURN
        ELSE
          IT_NEXT=X+DX
          RETURN
        ENDIF
      ELSE
        IT_NEXT=X*DX
        RETURN
      ENDIF
      END FUNCTION
      
      LOGICAL FUNCTION IT_AT_END1(X,TOWARD)
      USE Marc,ONLY:GO_DOWN
      IMPLICIT NONE
      REAL X,TOWARD
      IF (GO_DOWN)THEN
        IT_AT_END1=(X <= TOWARD)
        RETURN
      ELSE
        IT_AT_END1=(X >= TOWARD)
        RETURN
      ENDIF
      END FUNCTION

      LOGICAL FUNCTION IT_AT_END2(X,MINIMUM,MAXIMUM)
      USE Marc,ONLY:GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::X,MINIMUM,MAXIMUM
      IF (GO_DOWN)THEN
        IT_AT_END2=(X<=MINIMUM)
        RETURN
      ELSE
        IT_AT_END2=(X>=MAXIMUM)
        RETURN
      ENDIF
      END FUNCTION
      
      REAL FUNCTION IT_BEGIN(MINIMUM,MAXIMUM)
      USE Marc,ONLY:GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::MINIMUM,MAXIMUM
      IF (GO_DOWN)THEN
        IT_BEGIN=MAXIMUM
        RETURN
      ELSE
        IT_BEGIN=MINIMUM
        RETURN
      ENDIF
      END FUNCTION

      REAL FUNCTION IT_END(MINIMUM,MAXIMUM)
      USE Marc,ONLY:GO_DOWN
      IMPLICIT NONE
      REAL,INTENT(IN)::MINIMUM,MAXIMUM
      IF (GO_DOWN)THEN
       IT_END=MINIMUM
       RETURN
      ELSE
       IT_END=MAXIMUM
       RETURN
      ENDIF
      END FUNCTION

      LOGICAL FUNCTION SEARCH1(FROM,TOWARD,DX)
      !Outputs .TRUE. IF THE ROOT HAS BEEN BRACKETED BETWEEN G_X1 AND G_X2
      USE Marc,ONLY:G_X1,G_Y1,G_X2,G_Y2
      IMPLICIT NONE
      REAL,INTENT(IN)::FROM,TOWARD,DX
      ! FUNCTIONS
      LOGICAL ISFOUND,IT_AT_END1
      REAL Y,IT_NEXT
      G_X1=FROM
      G_Y1=Y(G_X1)
      G_X2=IT_NEXT(G_X1,DX)
      ! G_X2 IS CLOSER TO TOWARD, G_X1 IS FURTHER
      DO WHILE(.NOT.IT_AT_END1(G_X2,TOWARD))
        G_y2=y(G_x2);
        IF(ISFOUND(G_Y1,G_Y2))THEN
          SEARCH1=.TRUE.;RETURN
        ELSE
          G_x1=G_x2;
          G_y1=G_y2;
          G_X2=IT_NEXT(G_X1,DX);
        ENDIF
      ENDDO
!     Reaching the maximum slowness, tests if there's a last root
      IF(.NOT.(G_x1 >= TOWARD))THEN
        ! ATTEMPT OF BRACKETING WITH G_X2=TOWARD
        G_x2=TOWARD;
        G_y2=y(G_x2);
        IF(isFound(G_y1,G_y2))THEN
          SEARCH1=.TRUE.;RETURN
        ENDIF
      ENDIF
      SEARCH1=.FALSE.;RETURN
      END FUNCTION

      LOGICAL FUNCTION SEARCH2(FROM,MINIMUM,MAXIMUM,DX)
      USE Marc,ONLY:G_POLARITY,G_X1,G_X2,G_Y1,G_Y2
      IMPLICIT NONE
      REAL,INTENT(IN)::FROM,MINIMUM,MAXIMUM,DX
      ! FUNCTIONS
      REAL IT_NEXT,IT_BEGIN,IT_END,Y
      LOGICAL ISFOUND,IT_AT_END2
!     Find the search direction by looking at polarity
      G_X1=FROM
      G_Y1=y(G_X1)
      G_Y2=G_POLARITY
      IF(isFound(G_y1,G_y2))THEN
        IF(IT_AT_END2(FROM,MINIMUM,MAXIMUM))THEN
!         If _x1 reaches the min (there is necessarly a solution
!         greater than min because signs are different), there is a lot of
!         chance that a mode jumping occured. The only one solution is to
!         reduce the step size and recalculating all the curve until not
!         getting this condition.
          SEARCH2=.FALSE.;RETURN
        ELSE
          G_x2=G_x1;
          G_y2=G_y1;
          G_X1=IT_BEGIN(MINIMUM,MAXIMUM);
          G_Y1=G_POLARITY
          SEARCH2=.TRUE.;RETURN
        ENDIF
      ELSE
!       Root is situated below x0, between g_x1 and min, thus search by step
        G_x2=it_next(G_x1,DX);
        DO WHILE(.NOT. IT_AT_END2(G_X2,MINIMUM,MAXIMUM))
          G_y2=y(G_x2)
          IF(isFound(G_y1, G_y2))THEN
            SEARCH2=.TRUE.;RETURN
          ELSE
            G_x1=G_x2
            G_y1=G_y2
            G_x2=IT_NEXT(G_x1,DX)
          ENDIF
        ENDDO
      ENDIF
!     reaching the minimum slowness, tests if there's a last root
      IF(.NOT. IT_AT_END2(G_x1,MINIMUM,MAXIMUM))THEN
        G_x2=it_end(MINIMUM, MAXIMUM);
        G_y2=y(G_x2);
        IF(isFound(G_y1,G_y2))THEN
          SEARCH2=.TRUE.;RETURN
        ENDIF
      ENDIF
!     If x1 reaches the min/max (there is maybe a solution
!     greater/less than min/max, signs are equal), there is a lot of
!     chance that a mode jumping occured for fundamental mode. The only one solution is to
!     reduce the step size and recalculating all the curve until not
!     getting this condition. For higher modes this condition can occur.
      SEARCH2=.FALSE.;RETURN
      END FUNCTION

      SUBROUTINE neville()
      USE Marc,ONLY:G_PRECISION,G_X1,G_X2,G_Y1,G_Y2,ROOTSOLVER_MAX_ITERATIONS,Y      
!     Variables for polynomial fit of Neville
      IMPLICIT NONE
      REAL::TMPABS1,TMPABS2,DENOM,X3,Y3,maxDiff,TNY,NUM
      INTEGER::I,J,M
      REAL,DIMENSION(0:19)::p,nevY
      LOGICAL ISFOUND
      m = 0
      TNY=TINY(TMPABS2)
!     First halving
      x3=0.5 * (G_x1 + G_x2);
      y3=y(x3);
!     Relative Convergence criteria
      maxDiff=abs(x3)*G_precision;
      DO i=1,ROOTSOLVER_MAX_ITERATIONS-1
!       define new bounds according to signs of y1 and y3
!       If y1 and y3 have different signs then the root is between
!       y1 and y3, thus forget (x2,y2) else forget (x1,y1).
!       this case is not implemented in Herrmann's code
!       In our case, it is very important to closely bracket the root
!       X1 and X2 must ALWAYS stay on the same sides of the root during all
!       Neville's iterations.
        IF(y3==0.0)THEN
!         Move the x3 towards the furthest limit
          IF(abs(G_x1-x3)>abs(G_x2-x3))THEN
            x3=0.5*(G_x1+x3);
          ELSE
            x3=0.5*(G_x2+x3);
          ENDIF
!         recalculate the Y3, it will always be different from 0
          Y3=Y(X3)
        ENDIF
        IF(isFound(G_y1, y3))THEN
          G_x2=x3
          G_y2=y3
        ELSE
          G_x1=x3
          G_y1=y3
        ENDIF
!       Check for convergence with a relative criteria, if so exit
        IF(abs(G_x1-G_x2)<=maxDiff)THEN
          RETURN
        ENDIF
!       Check the slopes between x1,x3 and x3,x2
!       If they are different, the function is strongly non linear,
!       thus use only halving, not Neville
        IF(isFound(G_y1 - y3, y3 - G_y2))THEN
          x3=0.5* (G_x1 + G_x2)
          y3=y(x3)
          m=0
        ELSE
          IF(M>0)THEN
            p(m + 1)=x3
            nevY(m + 1)=y3
          ELSE
            p(0)=G_x1
            nevY(0)=G_y1
            p(1)=G_x2
            nevY(1)=G_y2
          ENDIF
!         perform Neville iteration. Instead of generating y(x)
!         we interchange the x and y of formula to solve for x(y)
!         when y=0.
!         A new point is added to the pyramidal contruction.
!         The new point is always added at the base, adding a new
!         slice at each iteration, each time larger
!         As pyramide grows, the degree of the polynome increase
!         thus the number of point taken into account increase also.
          DO J=M,0,-1
            denom=nevY(m + 1)-nevY(j)
            tmpAbs1=nevY( m + 1 )
            IF(tmpAbs1 < 0) tmpAbs1=-tmpAbs1
            tmpAbs2=denom
            IF(tmpAbs2 < 0) tmpAbs2=-tmpAbs2

            NUM=(nevY(m+1)*p(j)-nevY(j)*p(j+1))                            
            IF(EXPONENT(NUM)-EXPONENT(TMPABS2)<RANGE(NUM)-1)THEN ! POSSIBLE OVERFLOW
              p(0)=0.5*(G_x1+G_x2);
              m=0
              EXIT
            ELSE
              p(j)=NUM/denom;                                    
            ENDIF
          ENDDO
          x3=p(0);
!         Just move the X3 by one tenth of actual range towards the limit with the highest y
          tmpAbs1=ABS(G_Y1)
          tmpAbs2=ABS(G_Y2)     
          IF(tmpAbs1 > tmpAbs2)THEN
            x3 =x3 - 0.1 * (x3 - G_x1);
          ELSE
            x3 =x3 - 0.1 * (x3 - G_x2);
          ENDIF
          m=m+1
          IF(m > 18)m=1
!         Make sure new estimate is inside the previous values
!         if not perform interval halving
          IF(G_x1 <= G_x2)THEN
            IF(x3 < G_x1 .OR. x3 > G_x2)THEN
              x3=0.5* (G_x1 + G_x2);
              m=1;
            ENDIF
          ELSE
            IF(x3 < G_x2 .OR. x3 > G_x1)THEN
              x3=0.5*(G_x1 + G_x2);
              m=1;
            ENDIF
          ENDIF
          y3=y(x3);
        ENDIF
      ENDDO                                           
      END SUBROUTINE