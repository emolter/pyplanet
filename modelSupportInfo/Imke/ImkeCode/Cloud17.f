C
C***********************************************************************
C**                                                                   **
C**                      Program CLOUD17.FOR                          **
C**                                                                   **
C**                          Paul Romani                              **
C**                                                                   **
C***********************************************************************
C**                                                                   **
C**  CLOUD15.FOR computes the cloud densities, temperature            **
C**  structures, and mixing ratio profiles of condensates in the      **
C**  tropospheres of the outer planets. The lapse rate is user        **
C**  controlled and assumed to be adiabatic. It is the dry lapse      **
C**  rate plus LAPSE percent of the wet correction. So the lapse      **
C**  rate can vary from the pure dry, to dry plus some percent of     **
C**  the wet correction, to the appropiate wet adiabat. Condensation  **
C**  occurs at equilibrium saturation vapor pressures. The possible   **
C**  clouds are; NH3-H20 solution of variable concentration,          **
C**  H2O ice, NH3 ice, H2S ice, CH4 ice, Ar ice, and NH4SH solid.     **
C**  H2S is allowed to disolve into the already formed aqueous        **
C**  ammonia cloud by using an empirical expression from Leyko.       **
C**  CLOUD15.FOR allows for ortho-para conversion of H2 by modifying  **
C**  the Cp of H2. The ortho-para equilibration process assumes       **
C**  complete equilibration at each altitude step and all of H2 to    **
C**  be in the vibrational level V = 0. Unlike its ancestor,          **
C**  CLOUD9.FOR, CLOUD15.FOR allows for the concentration of the      **
C**  NH3-H2O solution to come to equilibrium. Unlike its predessor,   **
C**  CLOUD92.FOR, CLOUD15.FOR uses FORTRAN 77 structure, and allows   **
C**  for only thermodynamic equilibrium to control which cloud(s)     **
C**  form. CLOUD11.FOR, compared to CLOUD10.FOR which also did all    **
C**  of the above, fixed a bug in DISOLV the subroutine that          **
C**  disolves H2S into the solution cloud, and uses a faster          **
C**  algorthim to find the concentration of the solution cloud.       **
C**  CLOUD15.FOR, compared to previous versions of this program,      **
C**  uses new subroutines for the cubic spline cuve fitting for the   **
C**  NH3 & H2O vapor pressures above the aqueous ammonia solution.    **
C**  The difference is that the subrotuines were rewritten to get rid **
C**  of the multiple entry and return points in them.
C** 
C** Cloud17 contains Imke's FP ortho-para fraction as a function of z
C**  from PART in RT program THIS SEEMS TO WORK 
C Also. corrected for pressure dependence in CP (eq. 11; but note typo in paper:
C (1+0.01P e^(-T/70)  not T/270)                                    
C** Note that ENRHE = He/H2 ratio; not enrichment factor.                                                   
C                                                           
C***********************************************************************
C
C f77 -O -o cloud Cloud17.f
C ifort -132 Cloud17.f -o cloud.exe
C ./cloud.exe       to run it
C
C
C  Dimensions and Declerations
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION YA(20), CA(20), RTA(20), GTA(20)
	DIMENSION ADX(100), ADY1(100), ADY2(100), ADY3(100)
	DIMENSION AA1(100), AA2(100), AA3(100)
C
      DIMENSION YWLT(18), CWLT(18), RTWLT(18), GTWLT(18)
	DIMENSION WLTDX(100), WLTDY1(100), WLTDY2(100), WLTDY3(100)
	DIMENSION WLTA1(100), WLTA2(100), WLTA3(100)
C
      DIMENSION YWHT(20), CWHT(20), RTWHT(20), GTWHT(20)
	DIMENSION WHTDX(100), WHTDY1(100), WHTDY2(100), WHTDY3(100)
	DIMENSION WHTA1(100), WHTA2(100), WHTA3(100)
C
      DIMENSION C(800)
      DIMENSION P(800), T(800)
      DIMENSION XH2O(800), XNH3(800), XH2S(800), XCH4(800), XAR(800)
      DIMENSION FJ(10), DJ(10)
C
      REAL*8 LAPSE
      INTEGER ICASE, PLANET(2)
C
      COMMON /B1/ A2, A4, A5, CP, R, DZ, T1
      COMMON /B2/ S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12,
     1  S13, S14, S15
      COMMON /B3/ ENL, ENSH2O, ENSNH3, ENSH2S, ENH4SH, ENSCH4, ENSAR
      COMMON /B4/ XO, XN, XS, XC, XA
C
	COMMON /A1SPLN/ NA, YA, CA, RTA, GTA
	COMMON /A2SPLN/ ADX, ADY1, ADY2, ADY3, AA1, AA2, AA3
C
	COMMON /WLT1SP/ NWLT, YWLT, CWLT, RTWLT, GTWLT
	COMMON /WLT2SP/ WLTDX, WLTDY1, WLTDY2, WLTDY3, WLTA1, WLTA2, WLTA3
C
	COMMON /WHT1SP/ NWHT, YWHT, CWHT, RTWHT, GTWHT
	COMMON /WHTSPL/ WHTDX, WHTDY1, WHTDY2, WHTDY3, WHTA1, WHTA2, WHTA3
C
      COMMON /MBLOCK/ M0, M1, M2, M3, M4, M5, M6
      COMMON /OP/ RR, FACT, FJ, DJ
C
C  Open files; input then output
C
      OPEN (UNIT=1, NAME='INPUTSOL.DAT', STATUS='OLD')
      OPEN (UNIT=2, NAME='INPUTP.DAT', STATUS='OLD')
C
      OPEN (UNIT=3, NAME='THERMAL.DAT', STATUS='NEW')
      OPEN (UNIT=4, NAME='MIXRAT.DAT', STATUS='NEW')
      OPEN (UNIT=7, NAME='CLDDEN.DAT', STATUS='NEW')
      OPEN (UNIT=8, NAME='IMKEIN.DAT', STATUS='NEW')
C
C  Conversion factors
C
C  Q1 converts mm of Hg to dynes/cm**2,
C  Q2 converts atmospheres to dynes/cm**2
C  Q3 converts dynes/cm**2 to bars
C  Q4 converts cm to km
C
      Q1 = 1013250.0D0 / 760.0D0
      Q2 = 1013250.0D0
      Q3 = 1.0D-06
      Q4 = 1.0D-05
C
C  R is the univeresal gas constant in erg/Mole K
C
      R = 8.3143D+07
      RR = R
C
C  Read in the interpolation data for the subrotuines STPA, STPWLT,
C  and STPWHT. These subroutines compute the coefficients for the
C  the vapor pressure equations of ammonia and water over aqueous
C  ammonia solutions as a function of concentration.
C
      NA = 20
      NWLT = 18
      NWHT = 20
C
      DO 5 J = 1, NA
    5 READ(1,*) YA(J), CA(J), RTA(J), GTA(J)
      DO 10 J = 1, NWLT
   10 READ(1,*) YWLT(J), CWLT(J), RTWLT(J), GTWLT(J)
      DO 15 J = 1, NWHT
   15 READ(1,*) YWHT(J), CWHT(J), RTWHT(J), GTWHT(J)
C
C  Pass interpolation points to subroutines that calculate the spline
C  coefficients.
C
      CALL SETPA
      CALL SETPWL
      CALL SETPWH
C
C  Calculate the coefficients of the vapor pressure equations
C  to correct at low/high concentrations. See subrotuine PPRES
C  for details
C
      CCA = 0.0D0
      CCW = 1.0D0
      CALL SPLNA( CCA, S7, S8, S9 )
      CALL SPNWLT( CCW, S10, S11, S12 )
      CALL SPNWHT( CCW, S13, S14, S15 )
C
C  Read in initial data.
C
      READ (2,20) PLANET
   20 FORMAT (2A4)
      READ(2,*) ICASE, Z0, T(1), P(1), G, DZ, TTROP, LAPSE, SUPER
      READ(2,*) COSHE, COSCH4, COSH2O, COSNH3, COSH2S, COSAR
      READ(2,*) ENRHE, ENRCH4, ENRH2O, ENRNH3, ENRH2S, ENRAR
      READ(2,*) WH2, WHE, WH2O, WCH4, WNH3, WH2S, WAR, WNH4SH
      READ(2,*) ENL, ENSH2O, ENH4SH, ENSNH3, ENSH2S, ENSCH4, ENSAR
C
C  Calculate the mixing ratios of the gases in one mole of uncondensed
C  gas from their "cosmic" abundances realtive to H2 (COS variables)
C  the planets HE/H2 value, and assumed enrichments relative to their
C  cosmic values (ENR variables).
C
      XH2 = 1.0D0 / ( 1.0D0 + ENRHE  + ENRCH4*COSCH4 + ENRH2O*COSH2O +
     1  ENRNH3*COSNH3 + ENRH2S*COSH2S + ENRAR*COSAR )
      XHE = ENRHE  * XH2
      XCH4(1) = ENRCH4 * COSCH4 * XH2
      XH2O(1) = ENRH2O * COSH2O * XH2
      XNH3(1) = ENRNH3 * COSNH3 * XH2
      XH2S(1) = ENRH2S * COSH2S * XH2
      XAR(1) = ENRAR * COSAR * XH2
C
      XO = XH2O(1)
      XC = XCH4(1)
      XN = XNH3(1)
      XS = XH2S(1)
      XA = XAR(1)
C
      SUMF = XH2 + XHE + XO + XC + XN + XS + XA
      HEH2 = XHE / XH2
C
C  Calculate intial mean molecular weight
C
      W = XH2*WH2 + XHE*WHE + XO*WH2O + XC*WCH4 + XN*WNH3
     1  + XS*WH2S + XA*WAR
C
C  Initialize variables
C
      Z = Z0
      T1 = T(1)
      DTDZ = 0.0D0
      P(1) = P(1) / Q3
      P1 = P(1)
C
      C(1) = 0.0D0
      DC = 0.0D0
      DELC = 0.00001D0
      IWMAX = 99999
      CH2S = 0.0D0
C
      DXH2O = 0.0D0
      DXNH3 = 0.0D0
      DXH2S = 0.0D0
      DXCH4 = 0.0D0
      DXAR = 0.0D0
C
      DSOLN = 0.0D0
      DH2OS = 0.0D0
      DNH4SH = 0.0D0
      DNH3S = 0.0D0
      DH2SS = 0.0D0
      DCH4S = 0.0D0
      DARS = 0.0D0
C
      M0 = 0
      M1 = 0
      M2 = 0
      M3 = 0
      M4 = 0
      M5 = 0
      M6 = 0
      LAPSE = LAPSE / 100.0D0
C
C  Calculate rotational energy levels and degenercies of molecular
C  hydrogen (needed for subrotuine H2CP).
C
      CALL ORTPAR
C
C  Calculate intial mean molar heat capacity
C
      CALL PART(T1,S1)
      fp=S1
      CALL H2CP( T1, CPH2, fp )
      CPHE = 2.0786D+08
      CPH2O = ( 30.359D0 + 9.61D-03*T1 + 1.184D-06*T1**2 ) * 1.0D+07
      CPCH4 = ( 14.146D0 + 7.55D-02*T1 - 1.799D-05*T1**2 ) * 1.0D+07
      CPNH3 = ( 25.895D0 + 3.30D-02*T1 - 3.0460D-06*T1**2 ) * 1.0D+07
      CPH2S = ( 28.719D0 + 1.612D-02*T1 + 3.2840D-06*T1**2 ) * 1.0D+07
      CPAR = 2.0786D+08
      CPH2 = CPH2 * (1. + P1*Q3*1.0E-2*EXP(-T1/70.)) 
      CP = XH2*CPH2 + XHE*CPHE + XO*CPH2O + XC*CPCH4 + XN*CPNH3
     1   + XS*CPH2S + XA*CPAR
C
C  Define some frequently used combinations
C
      A1 = G * W
      A2 = -A1 / CP
      A3 = -A1 / R
      A4 = XH2O(1) / ( 1.0D0 - C(1) )
      A5 = CP * R
C
C  Echo input and write out intialized variables and constants
C  into first output file.
C
      ZZ = Z0 * Q4
      DZZ = DZ * Q4
      PP = P(1) * Q3
      WLAPSE = LAPSE * 100.0D0
C
      WRITE(3,100) PLANET, ICASE
  100 FORMAT( '1', 4X, 'Intial conditions for ', 2A4, '  Case ',
     1  I3, '  CLOUD15 Model Normal H2'/ ' '/ ' '/ )
      WRITE(3,105) ZZ, T(1), PP
  105 FORMAT( 5X, 'Base level; Height = ', F3.1, ' km, Temperature',
     1  ' = ', F6.2, ' K,  Pressure = ', 1PE9.3, ' Bars'/ ' ' )
      WRITE(3,110) G, DZZ, TTROP
  110 FORMAT( 5X, 'Acceleration due to gravity = ', F7.2, ' cm/sec**2',
     1  ',  Delta Z = ', F3.1, ' km,  Tropopause Temperature = ',
     2  F6.2, ' K'/ ' ' )
      WRITE(3,115) WLAPSE
  115 FORMAT(5X, 'The lapse rate is the dry adiabatic lapse rate ',
     1  'plus ', F5.1, '% of the wet corrrection.'/ ' ')
      WRITE(3,120) SUPER
  120 FORMAT(5X, 'The lapse rate is then multiplied by a factor of ', 
     1  F4.2/ ' ')
      WRITE(3,125) COSHE, COSCH4, COSH2O, COSNH3, COSH2S, COSAR
  125 FORMAT( 5X, 'Mixing ratios relative to H2 for a pure solar ',
     1  'atmosphere:'/ 5X, 'He = ', 0PF5.3, '  CH4 = ', 1PE8.2,
     2  '  H2O = ', 1PE8.2, '  NH3 = ', 1PE8.2, '  H2S = ', 1PE8.2,
     3  '  Ar = ', 1PE8.2/ ' ' )
      WRITE(3,130) ENRHE, ENRCH4, ENRH2O, ENRNH3, ENRH2S, ENRAR
  130 FORMAT( 5X, 'Planets He/H2 ratio and enrichments above solar ',
     1  'for the other gases:'/ 5X, 'He/H2 = ', F5.3, '  CH4 - ',
     2  F6.2, '  H2O - ', F6.2, '  NH3 - ', F6.2, '  H2S - ', F6.2, 
     3  '  Ar - ', F6.2/ ' ')
      WRITE(3,140) XH2, XHE, XC, XO, XN, XS, XA
  140 FORMAT( 5X, 'Initial Composition - Mole Fractions (Mixing ',
     1  'ratios)'/ 5X, 'H2 = ', 1PE8.2, '  He = ', 1PE8.2, 
     2  '  CH4 = ', 1PE8.2, '  H2O = ', 1PE8.2, '  NH3 = ', 1PE8.2, 
     3  '  H2S = ', 1PE8.2, '  Ar = ', 1PE8.2/ ' ' )
      WRITE(3,150) WH2, WHE, WH2O, WCH4, WNH3, WH2S, WAR, WNH4SH
  150 FORMAT( 5X, 'Molecular Weights - Grams/Mole'/ 5X, 'H2 = ', F5.3,
     1  '  He = ', F5.3, '  H2O = ', F6.3, '  CH4 = ', F6.3,
     2  '  NH3 = ', F6.3, '  H2S = ', F6.3, '  Ar = ', F6.3,
     3  '  NH4SH = ', F6.3/ ' ' )
      WRITE(3,160) W
  160 FORMAT( 5X, 'Intial Mean Molecular Weight = ', F6.3,
     1  ' Grams/Mole'/ ' ' )
      WRITE(3,170) CPH2, CPHE, CPH2O, CPCH4, CPNH3, CPH2S, CPAR
  170 FORMAT( 5X, 'Intial Specific Heats at Constant Pressure - ',
     1  'Ergs/Mole-K'/ 5X, 'H2 = ', 1PE9.3, '  He = ', 1PE9.3,
     2  '  H2O = ', 1PE9.3, '  CH4 = ', 1PE9.3, '  NH3 = ', 1PE9.3,
     3  '  H2S = ', 1PE9.3, '  Ar = ', 1PE9.3/ ' ' )
      WRITE(3,180) CP
  180 FORMAT(5X, 'Intial Mean Specific Heat at Constant Pressure',
     1  ' = ', 1PE9.3, ' Ergs/Mole-K'/ ' ')
      WRITE(3,190) ENL
  190 FORMAT(5X, 'Heat of Vaporization of NH3 - H2O Solution = ',
     1  1PE9.3, ' Ergs/Mole'/ ' ')
      WRITE(3,200) ENH4SH
  200 FORMAT(5X, 'Heat of Formation of NH4SH = ', 1PE9.3,
     1  ' Ergs/Mole'/ ' ')
      WRITE(3,210) ENSH2O, ENSNH3, ENSH2S, ENSCH4, ENSAR
  210 FORMAT(5X, 'Heats of Sublimation of Ices - Ergs/Mole'/ 5X,
     1  'H20 = ', 1PE10.4, '  NH3 = ', 1PE10.4, '  H2S = ',
     2  1PE10.4, '  CH4 = ', 1PE10.4, '  Ar = ', 1PE10.4)
C
C  Write out column headings for all output files. 
C
      WRITE(3,250) PLANET, ICASE
  250 FORMAT('1', 3X, 2A4, 3X, 'CASE ', I3, 3X, 'CLOUD15 Model ',
     1  'Normal H2'/ ' ')
      WRITE(3,255)
  255 FORMAT(5X, 'Thermal Structure/Solution Cloud'/ ' ')
      WRITE(3,260)
  260 FORMAT(69X, 'NH3 - H2O Solution'/ ' '/ 3X, 'Z', 7X, 'T', 8X, 'P',
     1  6X, 'DT/DZ', 5X, 'Cp H2', 6X,'fp', 5x, 'Cp Atm', 5X, 'W', 8X, 'Conc.',
     2  5X, 'Delta', 3X, 'Conc. H2S'/ 2X, '(km)', 4X, '(K)', 5X,
     3  '(bars)', 3X, '(K/km)', 2X, 7X,'(Erg/M-K)', 2X, '(Erg/M-K)', 2X,
     4  '(g/M)', 6X, '% NH3', 5X, 'Conc.', 3X, '(Mol/Lit)'/ ' ')
      WRITE(4,250) PLANET, ICASE
      WRITE(4,265)
  265 FORMAT(5X, 'Mole Fraction/Mixing Ratio Profiles'/ ' ')
      WRITE(4,270)
  270 FORMAT(3X, 'Z', 7X, 'T', 8X, 'P', 8X, 'H2O', 7X, 'Delta', 6X,
     1  'NH3', 7X, 'Delta', 6X, 'H2S', 7X, 'Delta', 6X, 'CH4', 7X,
     2  'Delta', 6X, 'Ar', 8X, 'Delta'/ 2X, '(km)', 4X, '(K)', 5X,
     3  '(bars)', 16X, 'H2O', 18X, 'NH3', 18X, 'H2S', 18X, 'CH4', 18X,
     4  'Ar'/' ')
      WRITE(7,250) PLANET, ICASE
      WRITE(7,275)
  275 FORMAT(5X, 'Cloud Densities - gm/cm**3'/ ' ')
      WRITE(7,280)
  280 FORMAT(3X, 'Z', 7X, 'T', 8X, 'P', 6X, 'NH3 - H2O', 5X, 'H2O',
     1  7X, 'NH4SH', 7X, 'NH3', 8X, 'H2S', 8X, 'CH4', 8X, 'Ar'/ 2X,
     2  '(km)', 4X, '(K)', 5X, '(bars)', 3X, 'Solution', 6X, 'Ice',
     3  7X, 'Solid', 7X, 'Ice', 8X, 'Ice', 8X, 'Ice', 8X, 'Ice'/ ' ')
C
C  For ease of plotting up cloud densities write out first level
C  values into cloud density file.
C
      WRITE(7,520) ZZ, T(1), PP, DSOLN, DH2OS, DNH4SH, DNH3S, DH2SS,
     1  DCH4S, DARS
C
C  Output file for Imke de Pater
C
      WRITE(8,530) ZZ, T(1), PP, XH2, XHE, XCH4(1), XNH3(1), XH2O(1),
     1  XH2S(1), DSOLN, DNH4SH
C 
C***********************************************************************
C
C  Beginning of iterative loop over height
C
C***********************************************************************
C
      J = 1
  290 J = J + 1
  295 IF ( J .GT. 800 ) GO TO 570
C
      Z = Z + DZ
C
C  Compute lapse rate, DTDZ, and new temperature variables
C
      CALL TEMP(DTDZ, LAPSE)
      DTDZ = DTDZ * SUPER
      DT = DZ * DTDZ
      T2 = T1 + DT
      T(J) = T2
      T1 = T2
C
C  Calculations end when tropopause temperature is reached
C
      IF (T(J) .LE. TTROP) GO TO 550
C
C  Compute new pressure variables
C
      TAVG = T(J-1) + DT*0.5D0
      P(J) = P(J-1) * DEXP(A3*DZ/TAVG)
      P1 = P(J)
      DP = P(J) - P(J-1)
      PAVG = P(J-1) + DP*0.5D0
C
C  Reset Variables
C
      C(J) = C(J-1)
      DC = 0.0D0
      CH2S = 0.0D0
C
      DSOLN = 0.0D0
      DH2OS = 0.0D0
      DNH3S = 0.0D0
      DH2SS = 0.0D0
      DNH4SH = 0.0D0
      DCH4S = 0.0D0
      DARS = 0.0D0
C
      DXH2O = 0.0D0
      DXNH3 = 0.0D0
      DXH2S = 0.0D0
      DXCH4 = 0.0D0
      DXAR = 0.0D0
C
      M0 = 0
      M1 = 0
      M2 = 0
      M3 = 0
      M4 = 0
      M5 = 0
      M6 = 0
C
      XH2O(J) = XH2O(J-1)
      XNH3(J) = XNH3(J-1)
      XH2S(J) = XH2S(J-1)
      XCH4(J) = XCH4(J-1)
      XAR(J) = XAR(J-1)
C
      XO = XH2O(J)
      XC = XCH4(J)
      XN = XNH3(J)
      XS = XH2S(J)
      XA = XAR(J)
C
C  If T > 188 K, the triple point for H2S, then a liquid rather than a
C  solid H2S cloud will form, adjust the latent heat release attributed
C  to the H2S cloud for that
C
      IF ( T(J) .GT. 188.0D0 ) THEN
        ENSH2S = 1.866D+11
      ELSE
        ENSH2S = 2.103D+11
      END IF
C
C  Compute new partial pressures
C
      PPH2O = XH2O(J) * P(J)
      PPNH3 = XNH3(J) * P(J)
      PPH2S = XH2S(J) * P(J)
      PPCH4 = XCH4(J) * P(J)
      PPAR = XAR(J) * P(J)
C
C  Compute the equilibrium saturation vapor pressures of the ices
C
      GOL10T = DLOG10( T(J) )
      TT = T(J) * T(J)
C
C  H2O over H2O ice
C
      VPH2OS = -2445.5646D0/T(J) + 8.2312D0*GOL10T - 0.01677006D0*T(J)
     1  + 1.20514D-05*TT - 6.757169D0
      VPH2OS = (10.0D0**VPH2OS) * Q1
C
C  NH3 over NH3 ice
C
      VPNH3S = -1790.00D0/T(J) - 1.81630D0*GOL10T + 14.97593D0
      VPNH3S = (10.0D0**VPNH3S) * Q1
C
C  Vapor of H2S over H2S liquid, T > 188 K, over H2S ice for T < 188 K
C
      IF ( T(J) .GT. 212.75 ) THEN
        VPH2SS = 7.7547D0 - 976.0D0/T(J) - 0.12058D0*GOL10T 
      ELSE IF ( T(J) .GT. 188.0 ) THEN
        VPH2SS = 15.4859D0 - 1264.3574D0/T(J) - 2.86206*GOL10T
      ELSE
        VPH2SS = 57.19D0 - 2461.84D0/T(J) - 18.443D0*GOL10T
      END IF
      VPH2SS = (10.0D0**VPH2SS) * Q1
C
C  CH4 over CH4 ice
C
      VPCH4S = 4.425070D0 - 453.92414D0/T(J) - 4055.6016D0/TT
     1  + 115352.19D0/T(J)**3 - 1165560.7D0/T(J)**4
      VPCH4S = (10.0D0**VPCH4S) * Q2
C
C  Ar over Ar ice
C
      VPARS = 6.204D0 - 237.0D0/( T(J) - 16.80D0 )
      VPARS = (10.0D0**VPARS) * Q1
C
C  Test for the formation of the aqueous ammonia solution cloud.
C
C  Check for a solultion cloud only if; (1) the temperature is greater 
C  than the lowest possible temperature a NH3-H2O solution can be 
C  liquid at, (2) the temperature is less than the highest possible 
C  temperature a NH3-H2O solution cloud can form at (the crtical 
C  temperature of water), and (3) there is NH3 left in the atmosphere.
C
      IF ( 180.0D0 .LT. T(J) .AND. T(J) .LT. 647.0D0 .AND. 
     1   XNH3(J) .GT. 0.0D0 ) THEN
C
C  Find C1, the lowest possible concentration of the solution, i.e. 
C  start at pure H2O and keep increasing the conecntration of the
C  solution (thus lowering the vapor pressure of H2O) until the water
C  vapor pressrue is less than the water partial pressure.
C
        IW = -1
  300   IW = IW + 1
        IF ( IW .GT. IWMAX ) GO TO 345
        C1 = DELC * DFLOAT( IW )
        IF ( T1 .LT. 300.0D0 ) THEN
          CALL SPNWLT(C1, S4, S5, S6)
        ELSE
          CALL SPNWHT(C1, S4, S5, S6)
        END IF
        CALL PPRESW(C1, VPH2OL)
        IF ( VPH2OL .GT. PPH2O ) GO TO 300
        IF ( T(J) .LT. 273.15D0 .AND. VPH2OL .GT. VPH2OS ) GO TO 300
C
C  C1 has the minimum posible NH3 vapor pressrue above it, only if 
C  this vapor pressure is smaller than the NH3 partial pressure can
C  a solution cloud form.
C
        IF ( IW .EQ. 0 ) THEN
          C1 = DELC
          IW = 1
        END IF
        CALL SPLNA(C1, S1, S2, S3)
        CALL PPRESA(C1, VPNH3L)
        IF ( VPNH3L .LT. PPNH3 ) THEN
C
C  Now find C2, the maximum possible concentration that the solution 
C  can form at, i.e. start at C1 and keep increasing the concentration
C  until the NH3 vapor pressure exceeds the NH3 partial pressure
C
          IA = IW
  315     IA = IA + 1
          C2 = DELC * DFLOAT( IA )
          CALL SPLNA(C2, S1, S2, S3)
          CALL PPRESA(C2, VPNH3L)
          IF ( VPNH3L .LT. PPNH3 ) GO TO 315
          IA = IA - 1
          C2 = C2 - DELC
C
C  If fortously C1 = C2, check to see if it this an equilibrium 
C  solution (the amount of H2O and NH3 removed from the atmosphere
C  creates a solution of the same concentration that triggered
C  condensation).
C
          IF ( DABS( C1 - C2 ) .LE. DELC ) THEN
            CS = ( C1 + C2 ) / 2.0D0
            IF ( T1 .LT. 300.0D0 ) THEN
              CALL SPNWLT(CS, S4, S5, S6)
            ELSE
              CALL SPNWHT(CS, S4, S5, S6)
            END IF
            CALL PPRESW(CS, VPH2OL)
            DH2O = XH2O(J) - ( VPH2OL / P(J) )
            DNH3 = ( CS*DH2O ) / ( 1.0D0 - C1 )
            PPNH3P = ( XNH3(J-1) - DNH3 ) * P(J)
            CS1 = CS - DELC
            CS2 = CS + DELC
            CALL SPLNA(CS1, S1, S2, S3)
            CALL PPRESA(CS1, VPNH31)
            CALL SPLNA(CS2, S1, S2, S3)
            CALL PPRESA(CS2, VPNH32)
            IF ( ( VPNH31 .LE. PPNH3P ) .AND. ( PPNH3P .LE. VPNH32 ) )
     1        GO TO 340
            GO TO 345
          END IF
C
C  If not so lucky, then start at C2 and count down to C1 in steps 
C  of DELC calculating the vapor pressures of NH3 and H2O and 
C  seeing if an equilibrium solution forms.
C
          IS = IA + 1
  320     IS = IS - 1
          IF ( IS .LT. IW ) GO TO 345
          CS = DELC * DFLOAT( IS )
          IF ( T1 .LT. 300.0D0 ) THEN
            CALL SPNWLT(CS, S4, S5, S6)
          ELSE
            CALL SPNWHT(CS, S4, S5, S6)
          END IF
          CALL PPRESW(CS, VPH2OL)
          DH2O = XH2O(J) - ( VPH2OL / P(J) )
          DNH3 = ( CS*DH2O ) / ( 1.0D0 - CS )
          PPNH3P = ( XNH3(J-1) - DNH3 ) * P(J)
          CS1 = CS - DELC
          CS2 = CS + DELC
          CALL SPLNA(CS1, S1, S2, S3)
          CALL PPRESA(CS1, VPNH31)
          CALL SPLNA(CS2, S1, S2, S3)
          CALL PPRESA(CS2, VPNH32)
          IF ( ( VPNH31 .LE. PPNH3P ) .AND. ( PPNH3P .LE. VPNH32 ) )
     1      GO TO 340
          GO TO 320
C
C  HIZZAH!!!!!!! an equilibrium solution cloud formed, make it real.
C
  340     CONTINUE
          C(J) = CS
          DC = C(J) - C(J-1)
C
          PPH2O = VPH2OL
          XH2O(J) = PPH2O / P(J)
          XO = XH2O(J)
          DXH2O = XH2O(J) - XH2O(J-1)
C
          DXNH3 = ( CS*DXH2O ) / ( 1.0D0 - CS )
          XNH3(J) = XNH3(J-1) + DXNH3
          PPNH3 = XNH3(J) * P(J)
          XN = XNH3(J)
C
          A4 = XH2O(J) / ( 1.0D0 - C(J) )
          M0 = 1
C
C  Disolve H2S into the solution cloud
C
          CALL DISOLV(T1, P1, DXH2O, DXNH3, WH2O, WNH3, Q1, XS,
     1      DXH2S, CH2S)
          XH2S(J) = XH2S(J-1) + DXH2S
          XS = XH2S(J)
          PPH2S = XH2S(J) * P(J)
C
C  Finally, compute the density of the solution cloud.
C
          DSOLN = -(WNH3*DXNH3 + WH2O*DXH2O + WH2S*DXH2S) * PAVG / 
     1      (A1*DZ)
        END IF
      END IF
  345 CONTINUE
C
C  If solution cloud did not form, check for either a pure liquid
C  H2O cloud or a pure liquid NH3 cloud. If they form, they will be
C  catogrized as solution clouds, with the appropitate concentrations.
C
C  Pure liquid H2O cloud
C
      IF ( M0 .EQ. 0 .AND. T(J) .GT. 273.15D0 .AND.
     1  T(J) .LT. 647.0D0 ) THEN
        C1 = 0.0D0
        IF ( T(J) .LT. 300.0D0 ) THEN
          CALL SPNWLT(C1, S4, S5, S6)
        ELSE
          CALL SPNWHT(C1, S4, S5, S6)
        END IF
        CALL PPRESW(C1, VPH2OL)
        IF ( VPH2OL .LT. PPH2O) THEN
          C(J) = C1
          DC = C(J) - C(J-1)
          PPH2O = VPH2OL
          XH2O(J) = PPH2O / P(J)
          XO = XH2O(J)
          DXH2O = XH2O(J) - XH2O(J-1)
          A4 = XH2O(J)
          DSOLN = -( WH2O*DXH2O*PAVG ) / ( A1*DZ )
          M0 = 1
        END IF
      END IF
C
C  Pure liquid NH3 cloud
C
      IF ( M0 .EQ. 0 .AND. XNH3(J) .GT. 0.0D0 ) THEN
        C1 = 1.0D0 
        CALL SPLNA(C1, S1, S2, S3)
        CALL PPRESA(C1, VPNH3L)
        IF (VPNH3L .LT. PPNH3) THEN
          C(J) = C1
          DC = C(J) - C(J-1)
          PPNH3 = VPNH3L
          XNH3(J) = PPNH3 / P(J)
          XN = XNH3(J)
          DXNH3 = XNH3(J) - XNH3(J-1)
          A4 = XNH3(J)
          DSOLN = -( WNH3*DXNH3*PAVG ) / ( A1*DZ )
          M0 = 1
        END IF
      END IF
C
C  For each of the remaining clouds; 1) check to see if it forms.
C  If it does form, then; 2) reduce amount of condensate in the
C  atmosphere, 3) compute the cloud density, and 3) set the
C  appropiate lapse rate flag on; M1, M2, etc.
C
C  H2O ice cloud
C
      IF (VPH2OS .LT. PPH2O .AND. T(J) .LT. 273.15D0 ) THEN
        PPH2O = VPH2OS
        XH2O(J) = PPH2O / P(J)
        XO = XH2O(J)
        DXH2O = XH2O(J) - XH2O(J-1)
        DH2OS = -( WH2O*DXH2O*PAVG / (A1*DZ) )
        M1 = 1
      END IF
C
C  If the solution cloud formed, and was not pure liquid NH3, and the 
C  H2O then froze out as H2O ice, then remove the soltuion cloud.
C  Furthermore, if the solution cloud was not pure water, put the 
C  the NH3 and H2S back into the atmosphere before going on.
C
      IF ( M0 .EQ. 1 .AND. C(J) .LT. 1.0D0 .AND. M1 .EQ. 1 ) THEN
        M0 = 0
        DSOLN = 0.0D0
        IF ( C(J) .GT. 0.0D0 ) THEN
          XNH3(J) = XNH3(J-1)
          DXNH3 = 0.0D0
          XN = XNH3(J)
          PPNH3 = XNH3(J) * P(J)
          XH2S(J) = XH2S(J-1)
          DXH2S = 0.0D0
          XS = XH2S(J)
          PPH2S = XH2S(J) * P(J)
        END IF
      END IF
C
C  NH3 ice cloud
C
      IF (VPNH3S .LT. PPNH3) THEN
        PPNH3 = VPNH3S
        XNH3(J) = PPNH3 / P(J)
        XN = XNH3(J)
        DXNH3 = XNH3(J) - XNH3(J - 1)
        DNH3S = -( WNH3 * DXNH3 * PAVG / (A1 * DZ) )
        M2 = 1
      END IF
C
C  If the NH3 ice cloud formed then it is the equilibrium sink for NH3,
C  and the solution cloud must be removed. So, if the solution cloud 
C  formed, and was not pure water, and only the NH3 froze out as ice, 
C  then the solution cloud must now be removed. If necessary put any
C  H2O and H2S back into the atmosphere and recheck for pure H2O 
C  clouds (liquid and ice).
C
      IF ( M0 .EQ. 1 .AND. C(J) .GT. 0.0D0 .AND. M1 .EQ. 0 
     1  .AND. M2 .EQ. 1 ) THEN 
        M0 = 0
        DSOLN = 0.0D0
        IF ( C(J) .LT. 1.0D0 ) THEN
C
          XH2O(J) = XH2O(J-1)
          DXH2O = 0.0D0
          XO = XH2O(J)
          PPH2O = XH2O(J) * P(J)
C
          XH2S(J) = XH2S(J-1)
          DXH2S = 0.0D0
          XS = XH2S(J)
          PPH2S = XH2S(J) * P(J)
C
C  recheck for H2O ice/liquid cloud
C
          IF ( T(J) .LT. 273.15D0 ) THEN
            IF (VPH2OS .LT. PPH2O) THEN
              PPH2O = VPH2OS
              XH2O(J) = PPH2O / P(J)
              XO = XH2O(J)
              DXH2O = XH2O(J) - XH2O(J - 1)
              DH2OS = -(WH2O*DXH2O*PAVG/(A1*DZ))
              M1 = 1
            END IF
          ELSE
            C1 = 0.0D0
            IF ( T(J) .LT. 300.0D0 ) THEN
              CALL SPNWLT(C1, S4, S5, S6)
            ELSE
              CALL SPNWHT(C1, S4, S5, S6)
            END IF
            CALL PPRESW(C1, VPH2OL)
            IF ( VPH2OL .LT. PPH2O) THEN
              C(J) = C1
              DC = C(J) - C(J-1)
              PPH2O = VPH2OL
              XH2O(J) = PPH2O / P(J)
              XO = XH2O(J)
              DXH2O = XH2O(J) - XH2O(J-1)
              A4 = XH2O(J)
              DSOLN = -( WH2O*DXH2O*PAVG ) / ( A1*DZ )
              M0 = 1
            END IF
          END IF
        END IF
      END IF
C
C  H2S ice cloud
C
      IF (VPH2SS .LT. PPH2S) THEN
        PPH2S = VPH2SS
        XH2S(J) = PPH2S / P(J)
        XS = XH2S(J)
        DXH2S = XH2S(J) - XH2S(J-1)
        DH2SS = -(WH2S*DXH2S*PAVG/(A1*DZ))
        M3 = 1
      END IF
C
C  If the H2S ice cloud and solution cloud both formed, and the
C  solution cloud was not either pure H2O or NH3, remove H2S from
C  from solution cloud
C
      IF ( M0 .EQ. 1 .AND. M3 .EQ. 1 .AND. C(J) .GT. 0.0D0
     1   .AND. C(J) .LT. 1.0D0 ) THEN
        DSOLN = -( WNH3*DXNH3 + WH2O*DXH2O ) * PAVG / (A1*DZ)
      END IF
C
C  NH4SH - solid cloud
C
C  If there is no NH3 or H2S left, skip the whole bloddy mess.
C
      IF ( (XH2S(J) .GT. 0.0D0) .AND. (XNH3(J) .GT. 0.0D0) ) THEN
C
        B1 = PPNH3 * PPH2S
        GOLB2 = ( 14.83D0 - (4715.0D0/T(J)) )  
        B2 = (10.0D0**GOLB2) * Q2 * Q2
        IF (B1 .GT. B2) THEN
          B3 = XNH3(J-1) + XH2S(J-1)
          B4 = XNH3(J-1)*XH2S(J-1) - B2/P(J)**2
          B5 = B3**2 - 4.0D0 * B4
          B6 = (B3 - DSQRT(B5)) / 2.0D0
C
C  Insure physically real values of XNH3 and XH2S
C
          IF ( B6 .LE. 0.0D0 ) B6 = (B3 + DSQRT(B5)) / 2.0D0
          IF ( B6 .LE. 0.0D0 ) GO TO 440
          IF ( B6 .GE. XH2S(J-1) ) B6 = XH2S(J-1)
          IF ( B6 .GE. XNH3(J-1) ) B6 = XNH3(J-1)
C
          DXNH3 = -B6
          DXH2S = -B6
          XNH3(J) = XNH3(J-1) + DXNH3
          XH2S(J) = XH2S(J-1) + DXH2S
          PPNH3 = XNH3(J) * P(J)
          PPH2S = XH2S(J) * P(J)
          DNH4SH = (WNH4SH*B6) * PAVG / (A1*DZ)
          XS = XH2S(J)
          XN = XNH3(J)
          M4 = 1
        END IF
      END IF
  440 CONTINUE
C
C  If the NH4SH cloud and the solution cloud formed, and the
C  solution cloud was not pure water, then assume the NH4SH cloud is
C  the equilibrium sink for NH3. First, remove the solution cloud;
C  second, if the solution cloud was not pure NH3, dump the water
C  from the now defunct solution cloud back into the atmosphere and
C  recheck for a pure H2O water solution cloud or H2O ice cloud.
C  
      IF ( M4 .EQ. 1 .AND. M0 .EQ. 1 .AND. C(J) .GT. 0.0D0 ) THEN 
        M0 = 0
        DSOLN = 0.0D0
        IF ( C(J) .LT. 1.0D0 ) THEN
          XH2O(J) = XH2O(J-1)
          DXH2O = 0.0D0
          XO = XH2O(J)
          PPH2O = XH2O(J) * P(J)
          IF ( T(J) .LT. 273.15D0 ) THEN
            IF (VPH2OS .LT. PPH2O) THEN
              PPH2O = VPH2OS
              XH2O(J) = PPH2O / P(J)
              XO = XH2O(J)
              DXH2O = XH2O(J) - XH2O(J - 1)
              DH2OS = -(WH2O*DXH2O*PAVG/(A1*DZ))
              M1 = 1
            END IF
          ELSE
            C1 = 0.0D0
            IF ( T(J) .LT. 300.0D0 ) THEN
              CALL SPNWLT(C1, S4, S5, S6)
            ELSE
              CALL SPNWHT(C1, S4, S5, S6)
            END IF
            CALL PPRESW(C1, VPH2OL)
            IF ( VPH2OL .LT. PPH2O) THEN
              C(J) = C1
              DC = C(J) - C(J-1)
              PPH2O = VPH2OL
              XH2O(J) = PPH2O / P(J)
              XO = XH2O(J)
              DXH2O = XH2O(J) - XH2O(J-1)
              A4 = XH2O(J)
              DSOLN = -( WH2O*DXH2O*PAVG ) / ( A1*DZ )
              M0 = 1
            END IF
          END IF
        END IF
      END IF
C
C  If both the NH4SH cloud and the NH3 ice cloud formed, then sort out
C  where the NH3 goes. If the NH4SH solid and NH3 ice crystals are 
C  at equilibrium, then the abundance of NH3 in the atmosphere will be
C  controlled by the vapor pressure of NH3 The total amount of NH3
C  removed from the atmosphere must be enough so there is enough NH3
C  to form the necessary amount of NH4SH and still have some left over
C  to form NH3 ice crystals. If not, then the NH4SH is the sole
C  equilibrium sink for the NH3 and need to turn off the NH3 ice cloud.
C
      IF ( M2 .EQ. 1 .AND. M4 .EQ. 1 .AND. M3 .EQ. 0 ) THEN
        XNH3T = VPNH3S / P(J)
        DXNH3T = XNH3(J-1) - XNH3T
        XH2ST = B2 / ( P(J)**2 * XNH3T )
        DXH2ST = XH2S(J-1) - XH2ST
        IF ( DXH2ST .LT. 0.0D0 ) DXH2ST = XH2S(J-1)
        DXNH31 = DXH2ST
        DXNH32 = DXNH3T - DXNH31
        IF ( DXNH32 .GT. 0.0D0 ) THEN
          PPNH3 = VPNH3S
          XNH3(J) = XNH3T
          XN = XNH3(J)
          DXNH3 = -DXNH3T
          DNH3S = ( WNH3*DXNH32*PAVG ) / ( A1*DZ )
          XH2S(J) = XH2ST
          DXH2S = -DXH2ST
          PPH2S = XH2S(J) * P(J)
          XS = XH2S(J)
          DNH4SH = ( WNH4SH*DXH2ST*PAVG ) / ( A1*DZ )
        ELSE
          DNH3S = 0.0D0
          M2 = 0
        END IF
      END IF
C
C  If both the NH4SH cloud and the H2S ice cloud formed, then sort out
C  where the H2S goes. If the NH4SH solid and H2S ice crystals are 
C  at equilibrium, then the abundance of H2S in the atmosphere will be
C  controlled by the vapor pressure of H2S. The total amount of H2S
C  removed from the atmosphere must be enough so there is enough H2S
C  to form the necessary amount of NH4SH and still have some left over
C  to form H2S ice crystals. If not, then the NH4SH is the sole 
C  equilibrium sink for the H2S, and need to turn off the H2S ice cloud.
C
      IF ( M3 .EQ. 1 .AND. M4 .EQ. 1 .AND. M2 .EQ. 0) THEN
        XH2ST = VPH2SS / P(J)
        DXH2ST = XH2S(J-1) - XH2ST
        XNH3T = B2 / ( P(J)**2 * XH2ST )
        DXNH3T = XNH3(J-1) - XNH3T
        IF ( DXNH3T .LE. 0.0D0 ) DXNH3T = XNH3(J-1)
        DXH2S1 = DXNH3T
        DXH2S2 = DXH2ST - DXH2S1
        IF ( DXH2S2 .GT. 0.0D0 ) THEN
          PPH2S = VPH2SS
          XH2S(J) = XH2ST
          XS = XH2S(J)
          DXH2S = -DXH2ST
          DH2SS = ( WH2S*DXH2S2*PAVG ) / ( A1*DZ )
          XNH3(J) = XNH3T
          DXNH3 = -DXNH3T
          PPNH3 = XNH3(J) * P(J)
          XN = XNH3(J)
          DNH4SH = ( WNH4SH*DXNH3T*PAVG ) / ( A1*DZ )
        ELSE
          DH2SS = 0.0D0
          M3 = 0
        END IF
      END IF
C
C  If all three clouds formed; NH3 ice, H2S ice and NH4SH solid,
C  then the NH4SH cloud is the equilibrium sink for both NH3 and H2S.
C
      IF ( M2 .EQ. 1 .AND. M3 .EQ. 1 .AND. M4 .EQ. 1 ) THEN
        DNH3S = 0.0D0
        M2 = 0
        DH2SS = 0.0D0
        M3 = 0
      END IF
C
C  CH4 ice cloud
C
      IF (VPCH4S .LT. PPCH4) THEN
        PPCH4 = VPCH4S
        XCH4(J) = PPCH4 / P(J)
        XC = XCH4(J)
        DXCH4 = XCH4(J) - XCH4(J - 1)
        DCH4S = -(WCH4*DXCH4*PAVG) / (A1*DZ)
        M5 = 1
      END IF
C
C  Ar ice cloud
C
      IF (VPARS .LT. PPAR) THEN
        PPAR = VPARS
        XAR(J) = PPAR / P(J)
        XA = XAR(J)
        DXAR = XAR(J) - XAR(J - 1)
        DARS = -(WAR*DXAR*PAVG) / (A1*DZ)
        M6 = 1
      END IF
C
C  Entrain "dry air" into the parcel. As the various clouds form
C  the mixing ratios of the condensates will be reduced in the
C  ascending parcel. To keep the sum of the mixing ratios at its
C  original value (approx 1.0), add a mix of only H2 and He into the
C  parcel with the original H2/He ratio set at the base level.
C
      SUMFF = XH2 + XHE + XO + XC + XN + XS + XA
      DELTAF = SUMF - SUMFF
      DELH2 = DELTAF / (1.0D0 + HEH2)
      DELHE = DELTAF - DELH2
      XH2 = XH2 + DELH2
      XHE = XHE + DELHE
C
C  Calculate new mean molecular weight
C
      W = XH2*WH2 + XHE*WHE + XO*WH2O + XC*WCH4 + XN*WNH3
     1  + XS*WH2S + XA*WAR
C
C  Calculate new mean molar heat capacity
C
      CALL PART(T1,S1)
      fp=S1
      CALL H2CP( T1, CPH2, fp )
      CPH2O = ( 30.359D0 + 9.61D-03*T1 + 1.184D-06*T1**2 ) * 1.0D+07
      CPCH4 = ( 14.146D0 + 7.55D-02*T1 - 1.799D-05*T1**2 ) * 1.0D+07
      CPNH3 = ( 25.895D0 + 3.30D-02*T1 - 3.0460D-06*T1**2 ) * 1.0D+07
      CPH2S = ( 28.719D0 + 1.612D-02*T1 + 3.2840D-06*T1**2 ) * 1.0D+07
      CPH2 = CPH2 * (1. + P1*Q3*1.0E-2*EXP(-T1/70.)) 
      CP = XH2*CPH2 + XHE*CPHE + XO*CPH2O + XC*CPCH4 + XN*CPNH3
     1   + XS*CPH2S + XA*CPAR
C
C  Redefine some frequently used combinations
C
      A1 = G * W
      A2 = -A1 / CP
      A3 = -A1 / R
      A4 = XH2O(J) / ( 1.0D0 - C(J) )
      A5 = CP * R
C
C  Write out results
C
      ZZ = Z * Q4
      DTDZZ = DT / DZZ
      PP = P(J) * Q3
C
      WRITE(3,500) ZZ, T(J), PP, DTDZZ, CPH2, fp, CP, W, C(J), DC, CH2S
 500  FORMAT(1X, F5.1, 2X, F7.3, 2X, 1PE10.4, 2X, 0PF6.3, 2X, 1PE9.3, 2x, 
     1 0pf5.3,2X, 1PE9.3, 2X, 0PF5.2, 5X, 0PF7.5, 2X, 0PF8.5, 2X, 1PE8.2)
      WRITE(4,510) ZZ, T(J), PP, XH2O(J), DXH2O, XNH3(J), DXNH3,
     1  XH2S(J), DXH2S, XCH4(J), DXCH4, XAR(J), DXAR
  510 FORMAT(1X, F5.1, 2X, F6.2, 2X, 1PE9.3, 6(2X, 1PE8.2, 2X, 1PE9.2) )
      WRITE(7,520) ZZ, T(J), PP, DSOLN, DH2OS, DNH4SH, DNH3S, DH2SS,
     1  DCH4S, DARS
  520 FORMAT( 1X, F5.1, 2X, F6.2, 2X, 1PE9.3, 7(2X, 1PE9.3) )
      WRITE(8,530) ZZ, T(J), PP, XH2, XHE, XCH4(J), XNH3(J), XH2O(J),
     1  XH2S(J), DSOLN, DNH4SH
  530 FORMAT(1X, F5.1, 2X, F6.2, 2X, 1PE11.4, 2(2X, 0PF5.3),
     1  6(2X, 1PE8.2) )
C
C***********************************************************************
C
C  The end of the iterative loop over height
C
C***********************************************************************
C
  540 GO TO 290
C
  550 WRITE(3,560)
      WRITE(4,560)
      WRITE(7,560)
  560 FORMAT( ' '/ 5X, 'Next Temperature <= Tropopause Temperature'/
     1  5X, 'Execution Terminated' )

  570 CONTINUE
C
      STOP
      END
C
C  End of main program nothin' but subroutines from now on.
C
C***********************************************************************
C
      SUBROUTINE SETPA
C
C***********************************************************************
C
C  NAASA  6.1.002 SPLINE   FTN   03-16-80     
C  INTERPOLATION BY PIECEWISE CUBIC SPLINES
C
C  COMPUTES COEFFICIENTS FOR SPLINE INTERPOLATION
C  INPUT.. N = NUMBER OF DATA POINTS
C          X(I),Y1(I),Y2(I),Y3(I) = DATA TO BE INTERPOLATED
C          X(1).LT.X(2).LT. ... .LT.X(N)
C  OUTPUT.. PUT INTO COMMON BLOCK FOR USE BY SUBROTUINE SPLNA
C
      IMPLICIT REAL*8(A-H, O-Z)
C
	DIMENSION YA(20), CA(20), RTA(20), GTA(20)
	DIMENSION ADX(100), ADY1(100), ADY2(100), ADY3(100)
	DIMENSION AA1(100), AA2(100), AA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /A1SPLN/ NA, YA, CA, RTA, GTA
	COMMON /A2SPLN/ ADX, ADY1, ADY2, ADY3, AA1, AA2, AA3
C
      N = NA
      N1 = N - 1
C
	DO I = 1, N
	  X(I) = YA(I)
	  Y1(I) = CA(I)
	  Y2(I) = RTA(I)
	  Y3(I) = GTA(I)
	END DO
C
      T(1) = 0.E0
C
C     AA(I) = S''(X(I))/6.
C
      A1(1) = 0.E0
      A2(1) = 0.E0
      A3(1) = 0.E0
      A1(N) = 0.E0
      A2(N) = 0.E0
      A3(N) = 0.E0
C
      DO I = 1, N1
        DX(I) = X(I + 1) - X(I)
        DY1(I) = (Y1(I + 1) - Y1(I)) / DX(I)
        DY2(I) = (Y2(I + 1) - Y2(I)) / DX(I)
        DY3(I) = (Y3(I + 1) - Y3(I)) / DX(I)
      END DO
C
      DO I = 2, N1
        PIV = 2.E0 * (DX(I - 1) + DX(I)) - DX(I - 1) * T(I - 1)
        T(I) = DX(I) / PIV
        A1(I) = (DY1(I) - DY1(I - 1) - DX(I - 1)*A1(I - 1)) / PIV
        A2(I) = (DY2(I) - DY2(I - 1) - DX(I - 1)*A2(I - 1)) / PIV
        A3(I) = (DY3(I) - DY3(I - 1) - DX(I - 1)*A3(I - 1)) / PIV
      END DO
C
      DO IB = 2, N1
        I = N + 1 - IB
        A1(I) = A1(I) - T(I) * A1(I + 1)
        A2(I) = A2(I) - T(I) * A2(I + 1)
        A3(I) = A3(I) - T(I) * A3(I + 1)
	END DO
C
      DO I = 1, N1
	  ADX(I) = DX(I)
	  ADY1(I) = DY1(I)
	  ADY2(I) = DY2(I)
	  ADY3(I) = DY3(I)
	END DO
C
      DO I = 1, N
	  AA1(I) = A1(I)
	  AA2(I) = A2(I)
	  AA3(I) = A3(I)
	END DO
C
      RETURN
	END
C
C***********************************************************************
C
      SUBROUTINE SPLNA(XX, S1, S2, S3 )
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H, O-Z)
C
	DIMENSION YA(20), CA(20), RTA(20), GTA(20)
	DIMENSION ADX(100), ADY1(100), ADY2(100), ADY3(100)
	DIMENSION AA1(100), AA2(100), AA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /A1SPLN/ NA, YA, CA, RTA, GTA
	COMMON /A2SPLN/ ADX, ADY1, ADY2, ADY3, AA1, AA2, AA3
C
      N = NA
	NSAVE = N
      N1 = N - 1
	K = 1
C
	DO I = 1, N
	  X(I) = YA(I)
	  Y1(I) = CA(I)
	  Y2(I) = RTA(I)
	  Y3(I) = GTA(I)
	  A1(I) = AA1(I)
	  A2(I) = AA2(I)
	  A3(I) = AA3(I)
	END DO
C
      DO I = 1, N1
	  DX(I) = ADX(I)
	  DY1(I) = ADY1(I)
	  DY2(I) = ADY2(I)
	  DY3(I) = ADY3(I)
	END DO
C
C  INPUT.. XX = ANY VALUE
C  OUTPUT.. S1,S2,S3 = SPLINE EVALUATED AT XX
C
   40 IF (XX .LT. X(K)) GO TO 60
   50 IF (XX .GE. X(K + 1)) GO TO 70
C
      W = (XX - X(K)) / DX(K)
      V = 1.E0-W
      S1 = W * Y1(K + 1) + V * Y1(K) + DX(K) ** 2 * ((W**3 - W)*A1(K +
     1     1) + (V**3 - V)*A1(K))
      S2 = W * Y2(K + 1) + V * Y2(K) + DX(K) ** 2 * ((W**3 - W)*A2(K +
     1     1) + (V**3 - V)*A2(K))
      S3 = W * Y3(K + 1) + V * Y3(K) + DX(K) ** 2 * ((W**3 - W)*A3(K +
     1     1) + (V**3 - V)*A3(K))
	GO TO 999
C
   60 CONTINUE
      K = K - 1
      IF (K .NE. 0) GO TO 40
      K = 1
      S1 = Y1(1) + (DY1(1) - DX(1)*A1(2)) * (XX - X(1))
      S2 = Y2(1) + (DY2(1) - DX(1)*A2(2)) * (XX - X(1))
      S3 = Y3(1) + (DY3(1) - DX(1)*A3(2)) * (XX - X(1))
      GO TO 999
C
   70 CONTINUE
      K = K + 1
      IF (K .NE. NSAVE) GO TO 50
      K = NSAVE - 1
      S1 = Y1(NSAVE) + (DY1(K) + DX(K)*A1(K)) * (XX - X(NSAVE))
      S2 = Y2(NSAVE) + (DY2(K) + DX(K)*A2(K)) * (XX - X(NSAVE))
      S3 = Y3(NSAVE) + (DY3(K) + DX(K)*A3(K)) * (XX - X(NSAVE))
C
  999 CONTINUE
      RETURN
      END	
C
C***********************************************************************
C
      SUBROUTINE SETPWL
C
C***********************************************************************
C
C  NAASA  6.1.002 SPLINE   FTN   03-16-80     
C  INTERPOLATION BY PIECEWISE CUBIC SPLINES
C
C  COMPUTES COEFFICIENTS FOR SPLINE INTERPOLATION
C  INPUT.. N = NUMBER OF DATA POINTS
C          X(I),Y1(I),Y2(I),Y3(I) = DATA TO BE INTERPOLATED
C          X(1).LT.X(2).LT. ... .LT.X(N)
C  OUTPUT.. PUT INTO COMMON BLOCK FOR USE BY SUBROTUINE SPNWLT
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION YWLT(18), CWLT(18), RTWLT(18), GTWLT(18)
	DIMENSION WLTDX(100), WLTDY1(100), WLTDY2(100), WLTDY3(100)
	DIMENSION WLTA1(100), WLTA2(100), WLTA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /WLT1SP/ NWLT, YWLT, CWLT, RTWLT, GTWLT
	COMMON /WLT2SP/ WLTDX, WLTDY1, WLTDY2, WLTDY3, WLTA1, WLTA2, WLTA3
C
      N = NWLT
      N1 = N - 1
C
	DO I = 1, N
	  X(I) = YWLT(I)
	  Y1(I) = CWLT(I)
	  Y2(I) = RTWLT(I)
	  Y3(I) = GTWLT(I)
	END DO
C
      T(1) = 0.E0
C
C     AA(I) = S''(X(I))/6.
C
      A1(1) = 0.E0
      A2(1) = 0.E0
      A3(1) = 0.E0
      A1(N) = 0.E0
      A2(N) = 0.E0
      A3(N) = 0.E0
C
      DO I = 1, N1
        DX(I) = X(I + 1) - X(I)
        DY1(I) = (Y1(I + 1) - Y1(I)) / DX(I)
        DY2(I) = (Y2(I + 1) - Y2(I)) / DX(I)
        DY3(I) = (Y3(I + 1) - Y3(I)) / DX(I)
      END DO
C
      DO I = 2, N1
        PIV = 2.E0 * (DX(I - 1) + DX(I)) - DX(I - 1) * T(I - 1)
        T(I) = DX(I) / PIV
        A1(I) = (DY1(I) - DY1(I - 1) - DX(I - 1)*A1(I - 1)) / PIV
        A2(I) = (DY2(I) - DY2(I - 1) - DX(I - 1)*A2(I - 1)) / PIV
        A3(I) = (DY3(I) - DY3(I - 1) - DX(I - 1)*A3(I - 1)) / PIV
      END DO
C
      DO IB = 2, N1
        I = N + 1 - IB
        A1(I) = A1(I) - T(I) * A1(I + 1)
        A2(I) = A2(I) - T(I) * A2(I + 1)
        A3(I) = A3(I) - T(I) * A3(I + 1)
	END DO
C
      DO I = 1, N1
	  WLTDX(I) = DX(I)
	  WLTDY1(I) = DY1(I)
	  WLTDY2(I) = DY2(I)
	  WLTDY3(I) = DY3(I)
	END DO
C
      DO I = 1, N
	  WLTA1(I) = A1(I)
	  WLTA2(I) = A2(I)
	  WLTA3(I) = A3(I)
	END DO
C
      RETURN
	END
C
C***********************************************************************
C
      SUBROUTINE SPNWLT(XX, S1, S2, S3 )
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION YWLT(18), CWLT(18), RTWLT(18), GTWLT(18)
	DIMENSION WLTDX(100), WLTDY1(100), WLTDY2(100), WLTDY3(100)
	DIMENSION WLTA1(100), WLTA2(100), WLTA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /WLT1SP/ NWLT, YWLT, CWLT, RTWLT, GTWLT
	COMMON /WLT2SP/ WLTDX, WLTDY1, WLTDY2, WLTDY3, WLTA1, WLTA2, WLTA3
C
      N = NWLT
	NSAVE = N
      N1 = N - 1
	K = 1
C
	DO I = 1, N
	  X(I) = YWLT(I)
	  Y1(I) = CWLT(I)
	  Y2(I) = RTWLT(I)
	  Y3(I) = GTWLT(I)
	  A1(I) = WLTA1(I)
	  A2(I) = WLTA2(I) 
	  A3(I) = WLTA3(I)
	END DO
C
      DO I = 1, N1
	  DX(I) = WLTDX(I)
	  DY1(I) = WLTDY1(I) 
	  DY2(I) = WLTDY2(I)
	  DY3(I) = WLTDY3(I)
	END DO
C
C  INPUT.. XX = ANY VALUE
C  OUTPUT.. S1,S2,S3 = SPLINE EVALUATED AT XX
C
   40 IF (XX .LT. X(K)) GO TO 60
   50 IF (XX .GE. X(K + 1)) GO TO 70
C
      W = (XX - X(K)) / DX(K)
      V = 1.E0-W
      S1 = W * Y1(K + 1) + V * Y1(K) + DX(K) ** 2 * ((W**3 - W)*A1(K +
     1     1) + (V**3 - V)*A1(K))
      S2 = W * Y2(K + 1) + V * Y2(K) + DX(K) ** 2 * ((W**3 - W)*A2(K +
     1     1) + (V**3 - V)*A2(K))
      S3 = W * Y3(K + 1) + V * Y3(K) + DX(K) ** 2 * ((W**3 - W)*A3(K +
     1     1) + (V**3 - V)*A3(K))
	GO TO 999
C
   60 CONTINUE
      K = K - 1
      IF (K .NE. 0) GO TO 40
      K = 1
      S1 = Y1(1) + (DY1(1) - DX(1)*A1(2)) * (XX - X(1))
      S2 = Y2(1) + (DY2(1) - DX(1)*A2(2)) * (XX - X(1))
      S3 = Y3(1) + (DY3(1) - DX(1)*A3(2)) * (XX - X(1))
      GO TO 999
C
   70 CONTINUE
      K = K + 1
      IF (K .NE. NSAVE) GO TO 50
      K = NSAVE - 1
      S1 = Y1(NSAVE) + (DY1(K) + DX(K)*A1(K)) * (XX - X(NSAVE))
      S2 = Y2(NSAVE) + (DY2(K) + DX(K)*A2(K)) * (XX - X(NSAVE))
      S3 = Y3(NSAVE) + (DY3(K) + DX(K)*A3(K)) * (XX - X(NSAVE))
C
  999 CONTINUE
      RETURN
      END	
C
C***********************************************************************
C
      SUBROUTINE SETPWH
C
C***********************************************************************
C
C  NAASA  6.1.002 SPLINE   FTN   03-16-80     
C  INTERPOLATION BY PIECEWISE CUBIC SPLINES
C
C  COMPUTES COEFFICIENTS FOR SPLINE INTERPOLATION
C  INPUT.. N = NUMBER OF DATA POINTS
C          X(I),Y1(I),Y2(I),Y3(I) = DATA TO BE INTERPOLATED
C          X(1).LT.X(2).LT. ... .LT.X(N)
C  OUTPUT.. PUT INTO COMMON BLOCK FOR USE BY SUBROTUINE SPNWHT
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION YWHT(20), CWHT(20), RTWHT(20), GTWHT(20)
	DIMENSION WHTDX(100), WHTDY1(100), WHTDY2(100), WHTDY3(100)
	DIMENSION WHTA1(100), WHTA2(100), WHTA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /WHT1SP/ NWHT, YWHT, CWHT, RTWHT, GTWHT
	COMMON /WHTSPL/ WHTDX, WHTDY1, WHTDY2, WHTDY3, WHTA1, WHTA2, WHTA3
C
      N = NWHT
      N1 = N - 1
C
	DO I = 1, N
	  X(I) = YWHT(I)
	  Y1(I) = CWHT(I)
	  Y2(I) = RTWHT(I)
	  Y3(I) = GTWHT(I)
	END DO
C
      T(1) = 0.E0
C
C     AA(I) = S''(X(I))/6.
C
      A1(1) = 0.E0
      A2(1) = 0.E0
      A3(1) = 0.E0
      A1(N) = 0.E0
      A2(N) = 0.E0
      A3(N) = 0.E0
C
      DO I = 1, N1
        DX(I) = X(I + 1) - X(I)
        DY1(I) = (Y1(I + 1) - Y1(I)) / DX(I)
        DY2(I) = (Y2(I + 1) - Y2(I)) / DX(I)
        DY3(I) = (Y3(I + 1) - Y3(I)) / DX(I)
      END DO
C
      DO I = 2, N1
        PIV = 2.E0 * (DX(I - 1) + DX(I)) - DX(I - 1) * T(I - 1)
        T(I) = DX(I) / PIV
        A1(I) = (DY1(I) - DY1(I - 1) - DX(I - 1)*A1(I - 1)) / PIV
        A2(I) = (DY2(I) - DY2(I - 1) - DX(I - 1)*A2(I - 1)) / PIV
        A3(I) = (DY3(I) - DY3(I - 1) - DX(I - 1)*A3(I - 1)) / PIV
      END DO
C
      DO IB = 2, N1
        I = N + 1 - IB
        A1(I) = A1(I) - T(I) * A1(I + 1)
        A2(I) = A2(I) - T(I) * A2(I + 1)
        A3(I) = A3(I) - T(I) * A3(I + 1)
	END DO
C
      DO I = 1, N1
	  WHTDX(I) = DX(I)
	  WHTDY1(I) = DY1(I)
	  WHTDY2(I) = DY2(I)
	  WHTDY3(I) = DY3(I)
	END DO
C
      DO I = 1, N
	  WHTA1(I) = A1(I)
	  WHTA2(I) = A2(I)
	  WHTA3(I) = A3(I)
	END DO
C
      RETURN
	END
C
C***********************************************************************
C
      SUBROUTINE SPNWHT(XX, S1, S2, S3 )
C
C***********************************************************************
C
      IMPLICIT REAL*8(A-H, O-Z)
C
      DIMENSION YWHT(20), CWHT(20), RTWHT(20), GTWHT(20)
	DIMENSION WHTDX(100), WHTDY1(100), WHTDY2(100), WHTDY3(100)
	DIMENSION WHTA1(100), WHTA2(100), WHTA3(100)
C
      DIMENSION X(20), Y1(20), Y2(20), Y3(20)
      DIMENSION DX(100), DY1(100), DY2(100), DY3(100)
	DIMENSION A1(100), A2(100), A3(100), T(100)
C
	COMMON /WHT1SP/ NWHT, YWHT, CWHT, RTWHT, GTWHT
	COMMON /WHTSPL/ WHTDX, WHTDY1, WHTDY2, WHTDY3, WHTA1, WHTA2, WHTA3
C
      N = NWHT
	NSAVE = N
      N1 = N - 1
	K = 1
C
	DO I = 1, N
	  X(I) = YWHT(I)
	  Y1(I) = CWHT(I)
	  Y2(I) = RTWHT(I)
	  Y3(I) = GTWHT(I)
	  A1(I) = WHTA1(I)
	  A2(I) = WHTA2(I)
	  A3(I) = WHTA3(I)
	END DO
C
      DO I = 1, N1
	  DX(I) = WHTDX(I) 
	  DY1(I) = WHTDY1(I)
	  DY2(I) = WHTDY2(I)
	  DY3(I) = WHTDY3(I) 
	END DO
C
C  INPUT.. XX = ANY VALUE
C  OUTPUT.. S1,S2,S3 = SPLINE EVALUATED AT XX
C
   40 IF (XX .LT. X(K)) GO TO 60
   50 IF (XX .GE. X(K + 1)) GO TO 70
C
      W = (XX - X(K)) / DX(K)
      V = 1.E0-W
      S1 = W * Y1(K + 1) + V * Y1(K) + DX(K) ** 2 * ((W**3 - W)*A1(K +
     1     1) + (V**3 - V)*A1(K))
      S2 = W * Y2(K + 1) + V * Y2(K) + DX(K) ** 2 * ((W**3 - W)*A2(K +
     1     1) + (V**3 - V)*A2(K))
      S3 = W * Y3(K + 1) + V * Y3(K) + DX(K) ** 2 * ((W**3 - W)*A3(K +
     1     1) + (V**3 - V)*A3(K))
	GO TO 999
C
   60 CONTINUE
      K = K - 1
      IF (K .NE. 0) GO TO 40
      K = 1
      S1 = Y1(1) + (DY1(1) - DX(1)*A1(2)) * (XX - X(1))
      S2 = Y2(1) + (DY2(1) - DX(1)*A2(2)) * (XX - X(1))
      S3 = Y3(1) + (DY3(1) - DX(1)*A3(2)) * (XX - X(1))
      GO TO 999
C
   70 CONTINUE
      K = K + 1
      IF (K .NE. NSAVE) GO TO 50
      K = NSAVE - 1
      S1 = Y1(NSAVE) + (DY1(K) + DX(K)*A1(K)) * (XX - X(NSAVE))
      S2 = Y2(NSAVE) + (DY2(K) + DX(K)*A2(K)) * (XX - X(NSAVE))
      S3 = Y3(NSAVE) + (DY3(K) + DX(K)*A3(K)) * (XX - X(NSAVE))
C
  999 CONTINUE
      RETURN
      END	
C
C***********************************************************************
C
      SUBROUTINE PPRESA(C1, VPN)
C
C***********************************************************************
C
C  This subroutine computes the partial equilibrium saturation vapor
C  pressure of NH3 over an aqueous ammonia solution at concentration
C  of C1 in mole fraction of NH3.
C
      IMPLICIT REAL*8(A-H, O-Z)
      COMMON /B1/ A2, A4, A5, CP, R, DZ, T1
      COMMON /B2/ S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12,
     1  S13, S14, S15
C
      GOLT1 = DLOG(T1)
      GOLPN = S1 + S2 / T1 + S3 * GOLT1
      VPN = DEXP(GOLPN)
C
C  When C1 equals zero, the vapor pressure of NH3 should be zero.
C  But, alas, the vapor pressure equation will return a non zero value
C  due to the way the coefficients for the equation are calculated. To
C  correct this, subtact off the fictous non-zero vapor pressure in a
C  linear manner until at the lowest concentration for which there are
C  non-interpolated/extrapolated coefficients.
C
      IF ( C1 .LT. 0.05D0 ) THEN
        GOLP0N = S7 + S8 / T1 + S9 * GOLT1
        P0N = DEXP(GOLP0N)
        X = ( 0.05D0 - C1 ) / 0.05D0
        P0N = X * P0N
        VPN = VPN - P0N
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE PPRESW(C1, VPH)
C
C***********************************************************************
C
C  This subroutine computes the partial equilibrium saturation vapor
C  pressure of H2O over an aqueous ammonia solution at concentration
C  of C1 in mole fraction of NH3.
C
      IMPLICIT REAL*8(A-H, O-Z)
      COMMON /B1/ A2, A4, A5, CP, R, DZ, T1
      COMMON /B2/ S1, S2, S3, S4, S5, S6, S7, S8, S9, S10, S11, S12,
     1  S13, S14, S15
C
      GOLT1 = DLOG(T1)
      GOLPH = S4 + S5 / T1 + S6 * GOLT1
      VPH = DEXP(GOLPH)
C
C  When C1 equals one, the vapor pressure of H2O should be zero.
C  But, alas, the vapor pressure equation will return a non zero value
C  due to the way the coefficients for the equation are calculated. To
C  correct this, subtact off the fictous non-zero vapor pressure in a
C  linear manner until at the highest concentration for which there are
C  non-interpolated/extrapolated coefficients.
C
      IF ( C1 .GT. 0.95D0 ) THEN
        IF ( T1 .LT. 300.0D0 ) THEN
          GOLP0H = S10 + S11/T1 + S12*GOLT1
        ELSE
          GOLP0H = S13 + S14/T1 + S15*GOLT1
        END IF
        P0H = DEXP( GOLP0H )
        X = ( C1 - 0.9500 ) / 0.0500
        P0H = X * P0H
        VPH = VPH - P0H
      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE TEMP(DTDZ, LAPSE)
C
C***********************************************************************
C
C  This subroutine computes the change in temperature, i.e the lapse 
C  rate, using the dry adiabatic lapse rate plus LAPSE percent of the
C  wet corrrection. If LAPSE is 1.0 then the lapse rate is the 
C  appropiate wet adiabatic lapse rate. At each level the subroutine
C  first calculates the dry adiabatic lapse rate, then it tests to
C  see if any clouds formed. If a cloud formed and LAPSE is non-zero
C  it then calculates the wet correction to the dry adiabatic lapse
C  rate and adds LAPSE percent onto it. The wet correction is found by
C  using a Runge-Kutta method to calculate the wet adiabatic lapse
C  rate and subtracting off the dry lapse rate. To determine what/if
C  any clouds formed the subroutine evaluates the M flags which were
C  set in the main program as follows:
C
C          1 = condensation occured, 0 = no condensation
C
C          M0 - NH3 - H20 Solution
C          M1 - H2O Ice
C          M2 - NH3 Ice
C          M3 - H2S Ice
C          M4 - NH4SH Solid
C          M5 - CH4 Ice
C          M6 - Ar Ice
C
      IMPLICIT REAL*8(A-H, O-Z)
      REAL*8 LAPSE
      LOGICAL CLOUDS
      COMMON /B1/ A2, A4, A5, CP, R, DZ, T1
      COMMON /B3/ ENL, ENSH2O, ENSNH3, ENSH2S, ENH4SH, ENSCH4, ENSAR
      COMMON /B4/ XO, XN, XS, XC, XA
      COMMON /MBLOCK/ M0, M1, M2, M3, M4, M5, M6
C
C  The dry adiabtic lapse rate is just A2.
C
      DRY = A2
C
      CLOUDS = M0 .EQ. 1 .OR. M1 .EQ. 1 .OR. M2 .EQ. 1 .OR. M3 .EQ. 1
     1  .OR. M4 .EQ. 1 .OR. M5 .EQ. 1 .OR. M6 .EQ. 1
C
C  IF only the pure dry adiabtic lapse rate is desired or no clouds 
C  formed in the atmosphere THEN - the lapse rate is just the dry
C  adibat; ELSE - calcualte the wet adiabatic lapse rate and then the
C  desired lapse rate (mix of wet and dry).
C
      IF ( LAPSE .EQ. 0.0D0 .OR. .NOT. CLOUDS ) THEN
        DTDZ = DRY
      ELSE
C
C  First sum E1 and E2 over condensates 
C
        E1 = 0.00
        E2 = 0.00
C
        IF ( M0 .EQ. 1 ) THEN
          E1 = E1 + ( A4 * ENL ) / R
          E2 = E2 + ( A4 * ENL**2 ) / A5
        END IF
        IF ( M1 .EQ. 1 ) THEN
          E1 = E1 + ( XO * ENSH2O ) / R
          E2 = E2 + ( XO * ENSH2O**2 ) / A5
        END IF
        IF ( M2 .EQ. 1 ) THEN
          E1 = E1 + ( XN * ENSNH3 )/ R
          E2 = E2 + ( XN * ENSNH3**2 )/ A5
        END IF
        IF ( M3 .EQ. 1 ) THEN
          E1 = E1 + ( XS * ENSH2S ) / R
          E2 = E2 + ( XS * ENSH2S**2 ) / A5
        END IF
        IF ( M4 .EQ. 1 ) THEN
          XSXN = ( XS * XN ) / ( XS + XN )
          E1 = E1 + ( 2.0D0 * XSXN * ENH4SH ) / R
          E2 = E2 + ( XSXN * ENH4SH * R * 10845.2D0 ) / A5
        END IF
        IF ( M5 .EQ. 1 ) THEN
          E1 = E1 + ( XC * ENSCH4 ) / R
          E2 = E2 + ( XC * ENSCH4**2 ) / A5
        END IF
        IF ( M6 .EQ. 1 ) THEN
          E1 = E1 + ( XA * ENSAR ) / R
          E2 = E2 + ( XA * ENSAR**2 ) / A5
        END IF
C
C  Runge-Kutta approximation to the wet adiabatic lapse rate
C
        D0 = A2 * (1. + E1/T1)/(1. + E2/T1**2 )
        T11 = T1 + DZ * D0 / 2.
        D1 = A2 * (1. + E1/T11)/(1. + E2/T11**2 )
        T12 = T1 + DZ * D1 / 2.
        D2 = A2 * (1. + E1/T12)/(1. + E2/T12**2 )
        T13 = T1 + DZ * D2
        D3 = A2 * (1. + E1/T13)/(1. + E2/T13**2 )
        WET = ( D0 + 2.0*D1 + 2.0*D2 + D3 ) / 6.00
C
C  Finally, calculate the desired lapse rate
C
        WETCOR = WET - DRY
        DTDZ = DRY + ( LAPSE * WETCOR )

      END IF
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE DISOLV(TEMP, PRESRE, DFH2O, DFNH3, WTH2O, WTNH3,
     1  CMMHG, FH2S, DFH2S, CH2S )
C
C***********************************************************************
C
C   Paul Romani  June 1980
C
C   This subroutine calculates how much H2S disolves into the already
C   formed aqueous ammonia cloud. It does so by calculating, via an
C   iterative method, an emperical expression that gives the partial
C   pressure of H2S above an H2O-NH3-H2S solution as a function of
C   temperature and the concentrations of NH3 and H2S in moles/liter.
C   The source of the empircal expression is Leyko, Bull. Acad. Polon.
C   Sci. Ser. Chim., 12, 275-276, 1964. The density of the aqueous
C   ammonia cloud is from an equation given by Croft, Lunine, and Kargel
C   in Icarus 1987 (1988?). The density is a function of the mass
C   fraction of ammonia only, any temperature and pressure efects are
C   ignored. If convergence of the concentration of H2S in the solution
C   is not achieved (to 0.01% within 50 iterations), or the amount of
C   H2S in the atmosphere becomes zero or less, no H2S is allowed to
C   disolve into the solution. Aloha.
C
      IMPLICIT REAL*8(A-H, O-Z)
C
C  Intialize variables
C
      CM3LIT = 1000.0D0
      POWER = 1.0D0 / ( 1.130D0 + 1.8953D0 )
      OLDC = 0.0D0
C
C  Find the density of the solution cloud
C
      X = ( DFNH3*WTNH3 ) / ( DFNH3*WTNH3 + DFH2O*WTH2O )
      DENSOL = 0.9991D0 + X*( -0.4336D0 + X*( 0.3303D0 + X*( 0.2833D0
     1  + X*( -1.9716D0 + X*( 2.1396D0 - X*0.7294D0) ) ) ) )
C
C  Calculte the partial pressure of H2S in the atmosphere
C  in mm Hg.
C
      PPH2S = ( FH2S*PRESRE ) / CMMHG
C
C  Find the concentration of the solution cloud in moles
C  of NH3 / liter
C
      CM3SOL = -( DFH2O*WTH2O + DFNH3*WTNH3 ) / DENSOL
      VOLSOL = CM3SOL / CM3LIT
      CONSOL = -DFNH3 / VOLSOL
C
C  Iterate through the equation until a constant
C  CH2S is found. If not do nothing.
C
      FUNTMP = DEXP( 22.221D0 - 5388.0D0/TEMP )
      FUNNH3 = CONSOL**1.8953D0
      F = FUNNH3 / FUNTMP
C
      DO 10 I = 1 , 50
      CH2S = ( PPH2S*F )**POWER
      ECH2S = DABS( (CH2S - OLDC) / CH2S ) * 100.0D0
      IF ( ECH2S .LE. 0.001D0 ) GO TO 30
      OLDC = CH2S
      PPH2S = ( FH2S - (CH2S*VOLSOL) ) * ( PRESRE/CMMHG )
      IF ( PPH2S .LE. 0.0D0 ) GO TO 20
   10 CONTINUE
C
   20 DFH2S = 0.0D0
      CH2S = 0.0D0
      GO TO 40
C
   30 DFH2S = -( CH2S*VOLSOL )
C
   40 CONTINUE
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE ORTPAR
C
C***********************************************************************
C
C  This subroutine calculates the rotational energy levels (J vlaues 0
C  through 9) and the respective degenercies for molecular hydrogen.
C
C  N.B. the molecular hydrogen is assumed to be in the vibrational
C  level v = 0.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FJ(10), DJ(10)
      COMMON /OP/ RR, FACT, FJ, DJ
C
      H = 6.6262D-27
      BK = 1.3806D-16
      C = 2.9979D+10
      FACT = H * C / BK
C
      B0 = 59.322
      D0 = 4.71D-02
C
      DO 50 JC = 1, 10
      J = JC - 1
      FPJ = FLOAT( J )
      IF ( J/2 * 2 .EQ. J ) THEN
        JD = 2 * J + 1
      ELSE
        JD = 3 * ( 2 * J + 1 )
      END IF
      DJ(JC) = DFLOAT( JD )
      FJ(JC) = B0 * FPJ * ( FPJ + 1 ) - D0 * ( FPJ * ( FPJ + 1 ) )**2
   50 CONTINUE
C
      RETURN
      END
C***********************************************************************
C
      SUBROUTINE PART(TEMP,S1)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION CC(20),CC2(20)                                                  
C     SEE TATUM,AP J SUPPL SER,VOL 14,PAGE 21,1967                              
C     FOR THE EXPRESSIONS FOR THE BOLTZMANN STATISTICS                          
C     THAT APPLY TO THE GROUND STATE OF H2.                                     
C     HERZBERG,MOLECULAR SPECT AND MOL STRUCT,VOL 4,CONSTS                      
C     OF DIATOMIC MOLECULES IS USED FOR THE BV,DE VALUES.
C This is a subset from what is in Imke's RT programs                       
      HCK=1.4388                                                                
      BV=60.853-(3.062*0.5)                                                     
      DE=4.71E-2                                                                
      SNUC=0.5                                                                  
      CI1=SNUC/((2.0*SNUC)+1)                                                   
      CI2=(SNUC+1)/((2.0*SNUC)+1)                                               
      CONST=((2.0*0)+1)*((2.0*SNUC)+1)*((2.0*SNUC)+1)                           
      SUM=0.0                                                                   
      SUM2=0.0                                                                  
      L=0                                                                       
      L1=1                                                                      
      LR=5                                                                      
      DO 31 J=1,LR                                                              
      DG1=(2.0*L)+1.0                                                           
      DG2=(2.0*L1)+1.0                                                          
      J1=L*(L+1)                                                                
      J2=L1*(L1+1)                                                              
      F=(BV*J1)-(DE*J1*J1)                                                      
      F2=(BV*J2)-(DE*J2*J2)                                                     
      CALL EH2(L,E)                                                             
      F=E                                                                       
      CALL EH2(L1,E)                                                            
      F2=E                                                                      
      TX=HCK*(F/TEMP)                                                           
      CALL EX(TX,S)                                                             
      CC(J)=S*CONST*CI1*DG1                                                     
      SUM=SUM+CC(J)                                                             
      TX2=HCK*(F2/TEMP)                                                         
      CALL EX(TX2,S2)                                                           
      CC2(J)=S2*CONST*CI2*DG2                                                   
      SUM2=SUM2+CC2(J)                                                          
      L=L+2                                                                     
      L1=L1+2                                                                   
  31  CONTINUE                                                                  
      QROT=SUM+SUM2                                                             
      QR1=SUM                                                                   
      QR2=SUM2                                                                  
      S1=0.0                                                                    
      S2=0.0                                                                    
      S3=0.0                                                                    
      S4=0.0                                                                    
C     S1=PROB,TOTAL PARA,EQUILIBRIUM                                            
C     S2=PROB,TOTAL ORTHO,EQUILIBRIUM                                           
C     S3=PROB,TOTAL PARA,NORMAL                                                 
C     S4=PROB,TOTAL ORTHO,NORMAL                                                
C     P1=PROB,TOTAL ORTHO/PARA,EQUILIBRIUM                                      
C     P2=PROB,ORTHO J=1/PARA J=0,EQUILIBRIUM                                    
C     P3=PROB,TOTAL ORTHO/PARA,NORMAL                                           
C     P4=PROB,ORTHO J=1/PARA J=0,NORMAL                                         
C     P5=PROB,ORTHO J=1,EQUILIBRIUM                                             
C     P6=PROB,PARA J=0,EQUILIBRIUM                                              
C     P7=PROB,ORTHO J=1,NORMAL                                                  
C     P8=PROB,PARA J=0,NORMAL                                                   
      DO 32 J=1,LR                                                              
      S1=S1+(CC(J)/QROT)                                                        
      S2=S2+(CC2(J)/QROT)                                                       
      S3=S3+(CC(J)/QR1)                                                         
      S4=S4+(CC2(J)/QR2)                                                        
  32  CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EX(TX,S)                                                       
      IMPLICIT REAL*8(A-H,O-Z)
C     SINCE THE COMPUTER DOES NOT LIKE TOO LARGE OR SMALL ARGUMENTS             
C     FOR THE EXPONENTIAL FUNCTION, WE SET THE VALUE TO ZERO IF                 
C     NEED BE.                                                                  
      IF (TX .GT. 75.0) S=0.0                                                   
      IF (TX .LT. 75.0) S=EXP(-TX)                                              
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EH2(J,E)                                                       
      IMPLICIT REAL*8(A-H,O-Z)
C     ENERGY VALUE IN WAVENUMBERS(CM-1) FOR THE                                 
C     JTH ROTATIONAL LEVEL OF H2,SEE NBS 80,EQUATION 2.                         
      A1=J                                                                      
      A2=A1+1.00                                                                
      A3=A1*A1                                                                  
      A4=A2*A2                                                                  
      A5=A3*A1                                                                  
      A6=A4*A2                                                                  
      E=(59.3392*A1*A2)-(0.04599*A3*A4)+(5.2E-5*A5*A6)                          
      RETURN                                                                    
      END                                                                       
      SUBROUTINE H2CP( T, CPH2, fp)
C
C***********************************************************************
C
C  This subroutine calculates Cp at a given temperature
C  for the following types of molecular hydrogen;
C
C     CPEQ - equilibrium mixture of ortho and para hydrogen
C     CPORTH - orth H2 only
C     CPPARA - para H2 only
C     CPNORM - normal H2, aka unequilibrated H2, the ortho and para
C              fractions are set at their high temperature values of
C              3/4 and 1/4 respectively
C
C  The user then selects which Cp is desired by assigning it to the
C  variable CPH2 at the end of the subroutine
C
C  N.B. The molecular hydrogen is assumed to be in the vibrational
C  level v = 0, and the molecular constants and degenercies are
C  calculated by subroutine ORTPAR.
C
C
C     modified: returns a Cp corresponding to input fp
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION FJ(10), DJ(10)
      COMMON /OP/ RR, FACT, FJ, DJ
C
      FORTH = 0.75D0
      FPARA = 0.25D0
C
      ZEQ = 0.D0
      SEQ1 = 0.D0
      SEQ2 = 0.D0
C
      ZORTHO = 0.D0
      SO1 = 0.D0
      SO2 = 0.D0
C
      ZPARA = 0.D0
      SP1 = 0.0D0
      SP2 = 0.0D0
C
      DO 50 JJ = 1, 10
C
      J = 11 - JJ
      JQ = J - 1
      X = ( FACT / T ) * FJ(J)
      P = DJ(J) * DEXP( -X )
C
      ZEQ = ZEQ + P
      SEQ1 = SEQ1 + DJ(J) * X**2 * DEXP( -X )
      SEQ2 = SEQ2 + DJ(J) * X * DEXP( -X )
C
      IF ( JQ/2 * 2 .EQ. JQ ) THEN
        ZPARA = ZPARA + P
        SP1 = SP1 + DJ(J) * X**2 * DEXP( -X )
        SP2 = SP2 + DJ(J) * X * DEXP( -X )
      ELSE
        ZORTHO = ZORTHO + P
        SO1 = SO1 + DJ(J) * X**2 * DEXP( -X )
        SO2 = SO2 + DJ(J) * X * DEXP( -X )
      END IF
C
   50 CONTINUE
C
      CPEQ = 5.0D0/2.0D0 + SEQ1/ZEQ - ( SEQ2/ZEQ )**2
      CPEQ = CPEQ * RR
      CPORTH = 5.0D0/2.0D0 + SO1/ZORTHO - ( SO2/ZORTHO )**2
      CPORTH = CPORTH * RR
      CPPARA = 5.0D0/2.0D0 + SP1/ZPARA - ( SP2/ZPARA )**2
      CPPARA = CPPARA * RR
      CPNORM = FORTH*CPORTH + FPARA*CPPARA
C
C  Normal H2    CPH2=CPNORM 
C  Intermediate H2 ?? this may not be intermediate  CPH2 = fp*Cppara + (1.0 - fp)*Cporth
C  Equilibrium H2     CPH2 = CPEQ
      CPH2 = fp*Cppara + (1.0 - fp)*Cporth
C      CPH2 = CPEQ
C       CPH2 = CPNORM
C For the adiabats, CPNORM is needed, but for opacities CPEQ
C in fact, intermediate is best
C So in Paul's program we use CPNORM, and in Imke's we use equilibrium opacities
C
      RETURN
      END
C
