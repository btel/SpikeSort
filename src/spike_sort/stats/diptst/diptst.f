      SUBROUTINE DIPTST1(X,N,DIP,XL,XU,IFAULT,GCM,LCM,MN,MJ,DDX,DDXSGN)
C
C     ALGORITHM AS 217 APPL. STATIST. (1985) VOL.34, NO.3
C
C     Does the dip calculation for an ordered vector X using the
C     greatest convex minorant and the least concave majorant, skipping
C     through the data using the change points of these distributions.
C     It returns the dip statistic 'DIP' and the modal interval
C     (XL, XU).
C
C     MODIFICATIONS SEP 2 2002 BY F. MECHLER TO FIX PROBLEMS WITH
C	UNIMODAL (INCLUDING MONOTONIC) INPUT 
C
      REAL X(N)
      INTEGER MN(N), MJ(N), LCM(N), GCM(N), HIGH
      REAL ZERO, HALF, ONE
C     NEXT TWO LINES ARE ADDED
      REAL DDX(N)
      INTEGER DDXSGN(N), POSK, NEGK
      DATA ZERO/0.0/, HALF/0.5/, ONE/1.0/
C
      IFAULT = 1
      IF (N .LE. 0) RETURN
      IFAULT = 0
C
C     Check if N = 1
C
      IF (N .EQ. 1) GO TO 4
C
C     Check that X is sorted
C
      IFAULT = 2
      DO 3 K = 2, N
        IF (X(K) .LT. X(K-1)) RETURN
    3 CONTINUE
      IFAULT = 0
C
C     Check for all values of X identical,
C     and for 1 < N < 4.
C
      IF (X(N) .GT. X(1) .AND. N .GE. 4) GO TO 505
    4 XL = X(1)
      XU = X(N)
      DIP = ZERO
      RETURN
C The code amendment below is intended to be inseted above the line marked "5" in the original FORTRAN code
C The amendment checks the condition whether the input X is perfectly unimodal
C Hartigan's original DIPTST algorithm did not check for this condition
C and DIPTST runs into an infinite cycle for a unimodal input
C The condition that the input is unimodal is equivalent to having 
C at most 1 sign change in the second derivative of the input p.d.f.
C In MATLAB syntax, we check the flips in the function xsign=-sign(diff(1./diff(x)))=-sign(diff(diff(x)));
C with DDXSGN=xsign in the fortran code below
505   NEGK=0
      POSK=0
      DO 104 K = 3,N
      DDX(K) = X(K)+X(K-2)-2*X(K-1)
      IF (DDX(K) .LT. 0) DDXSGN(K) = 1
      IF (DDX(K) .EQ. 0) DDXSGN(K) = 0
      IF (DDX(K) .GT. 0) DDXSGN(K) = -1
      IF (DDXSGN(K) .GT. 0) POSK = K
      IF ((DDXSGN(K) .LT. 0) .AND. (NEGK .EQ. 0)) NEGK = K      
104   CONTINUE

C The condition check below examines whether the greatest position with a positive second derivative 
C is smaller than the smallest position with a negative second derivative
C The boolean check gets it right even if 
C the unimodal p.d.f. has its mode in the very first or last point of the input 

      IF ((POSK .GT. NEGK) .AND. (NEGK .GT. 0)) GOTO 5
      XL=X(1)
      XU=X(N)
      DIP=0
      IFAULT=5
      RETURN
C
C     LOW contains the index of the current estimate of the lower end
C     of the modal interval, HIGH contains the index for the upper end.
C
5     FN = FLOAT(N)
      LOW = 1
      HIGH = N
      DIP = ONE / FN
      XL = X(LOW)
      XU = X(HIGH)
C
C     Establish the indices over which combination is necessary for the
C     convex minorant fit.
C
      MN(1) = 1
      DO 28 J = 2, N
        MN(J) = J - 1
   25   MNJ = MN(J)
        MNMNJ = MN(MNJ)
        A = FLOAT(MNJ - MNMNJ)
        B = FLOAT(J - MNJ)
        IF (MNJ .EQ. 1 .OR. (X(J) - X(MNJ))*A .LT. (X(MNJ) - X(MNMNJ))
     +          *B) GO TO 28
        MN(J) = MNMNJ
        GO TO 25
   28 CONTINUE
C
C     Establish the indices over which combination is necessary for the
C     concave majorant fit.
C
      MJ(N) = N
      NA = N - 1
      DO 34 JK = 1, NA
        K = N - JK
        MJ(K) = K + 1
   32   MJK = MJ(K)
        MJMJK = MJ(MJK)
        A = FLOAT(MJK - MJMJK)
        B = FLOAT(K - MJK)
        IF (MJK .EQ. N .OR. (X(K) - X(MJK))*A .LT. (X(MJK) - X(MJMJK))
     +          *B) GO TO 34
        MJ(K) = MJMJK
        GO TO 32
   34 CONTINUE
C
C     Start the cycling.
C     Collect the change points for the GCM from HIGH to LOW.
C
   40 IC = 1
      GCM(1) = HIGH
   42 IGCM1 = GCM(IC)
      IC = IC + 1
      GCM(IC) = MN(IGCM1)
      IF (GCM(IC) .GT. LOW) GO TO 42
      ICX = IC
C
C     Collect the change points for the LCM from LOW to HIGH.
C
      IC = 1
      LCM(1) = LOW
   44 LCM1 = LCM(IC)
      IC = IC + 1
      LCM(IC) = MJ(LCM1)
      IF (LCM(IC) .LT. HIGH) GO TO 44
      ICV = IC
C
C     ICX, IX, IG are counters for the convex minorant,
C     ICV, IV, IH are counters for the concave majorant.
C
      IG = ICX
      IH = ICV
C
C     Find the largest distance greater than 'DIP' between the GCM and
C     the LCM from LOW to HIGH.
C
      IX = ICX - 1
      IV = 2
      D = ZERO
      IF (ICX .NE. 2 .OR. ICV .NE. 2) GO TO 50
      D = ONE / FN
      GO TO 65
   50 IGCMX = GCM(IX)
      LCMIV = LCM(IV)
      IF (IGCMX .GT. LCMIV) GO TO 55
C
C     If the next point of either the GCM or LCM is from the LCM,
C     calculate the distance here.
C
      LCMIV1 = LCM(IV - 1)
      A = FLOAT(LCMIV - LCMIV1)
      B = FLOAT(IGCMX - LCMIV1 - 1)
      DX = (X(IGCMX) - X(LCMIV1))*A / (FN*(X(LCMIV) - X(LCMIV1)))
     +          - B / FN
      IX = IX - 1
      IF (DX .LT. D) GO TO 60
      D = DX
      IG = IX + 1
      IH = IV
      GO TO 60
C
C     If the next point of either the GCM or LCM is from the GCM,
C     calculate the distance here.
C
   55 LCMIV = LCM(IV)
      IGCM = GCM(IX)
      IGCM1 = GCM(IX + 1)
      A = FLOAT(LCMIV - IGCM1 + 1)
      B = FLOAT(IGCM - IGCM1)
      DX = A / FN - ((X(LCMIV) - X(IGCM1))*B) / (FN * (X(IGCM)
     +          - X(IGCM1)))
      IV = IV + 1
      IF (DX .LT. D) GO TO 60
      D = DX
      IG = IX + 1
      IH = IV - 1
   60 IF (IX .LT. 1) IX = 1
      IF (IV .GT. ICV) IV = ICV
      IF (GCM(IX) .NE. LCM(IV)) GO TO 50
   65 IF (D .LT. DIP) GO TO 100
C
C     Calculate the DIPs for the current LOW and HIGH.
C
C     The DIP for the convex minorant.
C
      DL = ZERO
      IF (IG .EQ. ICX) GO TO 80
      ICXA = ICX - 1
      DO 76 J = IG, ICXA
        TEMP = ONE / FN
        JB = GCM(J + 1)
        JE = GCM(J)
        IF (JE - JB .LE. 1) GO TO 74
        IF (X(JE) .EQ. X(JB)) GO TO 74
        A = FLOAT(JE - JB)
        CONST = A / (FN * (X(JE) - X(JB)))
        DO 72 JR = JB, JE
          B = FLOAT(JR - JB + 1)
          T = B / FN - (X(JR) - X(JB))*CONST
          IF (T .GT. TEMP) TEMP = T
   72   CONTINUE
   74   IF (DL .LT. TEMP) DL = TEMP
   76 CONTINUE
C
C     The DIP for the concave majorant.
C
   80 DU = ZERO
      IF (IH .EQ. ICV) GO TO 90
      ICVA = ICV - 1
      DO 88 K = IH, ICVA
        TEMP = ONE / FN
        KB = LCM(K)
        KE = LCM(K + 1)
        IF (KE - KB .LE. 1) GO TO 86
        IF (X(KE) .EQ. X(KB)) GO TO 86
        A = FLOAT(KE - KB)
        CONST = A / (FN * (X(KE) - X(KB)))
        DO 84 KR = KB, KE
          B = FLOAT(KR - KB - 1)
          T = (X(KR) - X(KB))*CONST - B / FN
          IF (T .GT. TEMP) TEMP = T
   84   CONTINUE
   86   IF (DU .LT. TEMP) DU = TEMP
   88 CONTINUE
C
C     Determine the current maximum.
C
   90 DIPNEW = DL
      IF (DU .GT. DL) DIPNEW = DU
      IF (DIP .LT. DIPNEW) DIP = DIPNEW
      LOW = GCM(IG)
      HIGH = LCM(IH)
C
C     Recycle
C
      GO TO 40
C
  100 DIP = HALF * DIP
      XL = X(LOW)
      XU = X(HIGH)
C
      RETURN
      END

