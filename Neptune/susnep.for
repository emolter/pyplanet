   PROGRAM TB 
C  PROGRAM CALCULATES BRIGHTNESS TEMP. OF THE JOVIAN PLANETS 
C  AS A FUNCTION OF WAVELENGTH. 
C 
C  PROGRAM ESPECIALLY FOR Neptune, based upon Sushil/romani's cloud models
C  L AMMONIA 
C f77 -O -o sus susnep_v2003_orton.for
C
C ifort -132 susnep_2012_orton.for -o sus.exe
C ./sus.exe       to run it
C  man ifort to look for compiling options, like -132
C
C  This program has improved NH3 (joined Spilker and Joiner/Steffes
C  This program has H2S de Boer (H2S2, slightly better compared to lab data than
C  my older program, improved H2O, and PH3.
C Note: PH3 abundace is absent right now; need to add saturation curves, etc.
C
C  AS 1,VOL2:JTP.DAT    CONTAINS TEMP PRESS DATA IN UPPER ATM. 
C  AS 2,VOL2:JPAR.DAT       PARTIAL PRESS DATA ARE CALCUL AND PUT IN THIS 
C  AS 3,VOL2:JTPZ.DAT        TEMP PRESS ALT ARE CALC AND PUT HERE 
C  AS 5,VOL2:JAMMON.DAT    INPUT 
C  AS 6,SPR: 
C  AS 7,VOL2:KAKAR.DAT    CONTAINS AMMONIA LINES AND CP FOR H2 
C  RUN WITH JTP.JOB 
C  M5 EQ 0  CALCULATE TABLES; NE 0 ASSUMES TABLES EXIST ALLREADY 
C     VSET ARE WAVENUMBERS (CM-1). 
      DATA CSZA/.9978,.9911,.9798,.9638,.9428,.9165, 
     1  .8844,.8459,.8000,.750,.700,.650,.600,.500, 
     2  .400,.300,.200,.100/ 
      DATA DCSZA/.10,.0668,.0668,.0667,.0667,.0667, 
     1  .0667,.0666,.0640,.0571,.0492,.0429,.0530, 
     2  .0583,.0440,.0316,.0205,.0126/ 
      DATA STR/0.0,.10,.1666,.2333,.30,.3667,.4334, 
     1  .50,.5667,.6307,.6878,.737,.7799,.8329, 
     2  .8912,.9352,.9668,.9873,1.0/ 
C 
      OPEN(1,NAME='neptune.paulCO_cloud21_fletcher_best_dry')
      OPEN(10,NAME='nepcloud_CO.cloud21_fletcher_best_dry')
      Open(2,NAME='npar.dat')
      OPEN(3,NAME='ntpz.dat')
      OPEN(5,NAME='nammon.dat')
C     OPEN(6,NAME='nepco_spectotal_allgases_dry.dat')
      OPEN(6,NAME='test.dat')
      Open(7,NAME='kakar.dat')
      open(11,NAME='orton_H2.tables',form='formatted')
C     TO SUPPRESS THE LENGHTLY PRINT OUT OF INTERMEDIATE VALUES, 
C     SET NOPR AND NOPR2 TO 0 (NOT EQUAL TO 1). 
C     NOPR3 CONTROLS PRINT OUT OF MIXING RATIOS 
C     DEN=INTEGRATED CLOUD DENSITY; ALT=ALTITUDE RANGE OF CLOUD
C     for some reason I thought I calculated H2-CH4 I guess, 
C               but this hasn't been done
      CZB=0.0
      DEN=7.2e6
      NOPR=0 
      NOPR2=0 
      NOPR3=0
      NOPR4=1

      READ(5,779) MII
      READ(5,779) KQQ 
  779 FORMAT(I2) 
      READ(5,773) SIZE
C      WRITE(6,772) SIZE
      SIZE=SIZE/2.0
  772 FORMAT(2X,'WATER DROPLET SIZE',2X,E10.5,/)
  773 FORMAT(F10.4)
C     EQUIL: IQ=2, FP=S1, REB=1.0
C     NORMAL: IQ=0; FP=0.25; REB=0
C     INTERM: IQ=1; REB=1.0.; FP=S1
C     IQ=2: TEMP FOR EQUIL. CASE
      READ (5,778) IQ,REB,FP 
      if(NOPR4.eq.1) WRITE(6,778) IQ,REB,FP 
  778 FORMAT(I2,2F10.4) 
C 
      DO 1357 I=1,18 
          RR(I)=ACOS(CSZA(I)) 
          RR(I)=SIN(RR(I)) 
          if(NOPR4.eq.1) WRITE(6,1345) CSZA(I),DCSZA(I),RR(I) 
 1345     FORMAT(3E15.5) 
 1357 CONTINUE 
      DO 2345 J=1,16 
          DO 2345 K=1,J 
              READ(7,3333) VKAK(J,K) 
 2345 CONTINUE 
	  DO 2346 J=1,16 
	      DO 2346 K=1,J 
              READ(7,2222) DVKAK(J,K) 
              IF(NOPR.EQ.1) WRITE(6,4444) J,K,VKAK(J,K),DVKAK(J,K) 
              VKAK(J,K)=1.0E-3*VKAK(J,K) 
 2346 CONTINUE 
 3333 FORMAT(F12.3) 
 2222 FORMAT(F8.1) 
C 
      READ(7,103) (TC(I),CPH2D(I),I=1,36) 
      READ(7,108) (TC(I),CPH2D(I),I=37,38) 
  108 FORMAT(2(F3.0,F4.3,1X)) 
      READ(7,104) (TC(I),CPH2D(I),I=39,44) 
  103 FORMAT(9(F3.0,F4.3,1X)) 
  104 FORMAT(6(F4.0,F4.3,1X)) 
  105 FORMAT(5(F5.0,1X,F7.3))
C     for equilibrium H2 
      IF(IQ.EQ.2) READ(7,103) (TC(I),CPH2D(I),I=1,24) 
C     WRITE(6,105) (TC(I),CPH2D(I),I=1,44)
	  close(7)
C     read in H2S lines  THIS IS OLD
c	  OPEN(7,NAME='h2s_file.list')
c	  read(7,7465) (F0(i),S0(i),EL(i),i=1,311)
c7465 format(F13.4,9x,F7.4,3x,F9.4,32x)


C 
      READ(5,1041) M5,NL
 1041 FORMAT(I1,1X,I4) 
      if(NOPR4.eq.1) WRITE(6,1046) M5,NL 
 1046 FORMAT(2X,'PARAMETER M5 IS EQUAL TO',/,2I5) 
 1049 FORMAT(F10.4) 
C     G=1130.0 
      G=1140.0 
      GE=G 
      if(NOPR4.eq.1) WRITE(6,1050) G 
 1050 FORMAT(2X,'THE GRAVITY G IS ',/,F12.4) 
C     READ(5,1045) AB1,AB2,AB3,AB4,AB5,AB6 
 1045 FORMAT(6E10.3) 
C     WRITE(6,1047) AB1,AB2,AB3,AB4,AB5,AB6 
 1047 FORMAT(2X,'ABUNDANCES H2,HE,NH3,H2O,CH4,H2S ARE',/,6E12.5) 
      IF(M5.EQ.0) CALL ATMOS(NL,SOL) 
      if(NOPR4.eq.1) WRITE(6,1051) NL
 1051 FORMAT(2X,'THE PARAMETER NL AFTER ATMOS IS EQUAL TO',/,I6)
      NL2=NL-1 
      REWIND 2
      REWIND 3
      do i=1,NL
          READ(3,23) ZT(I),ZDELZ(I),ZSOLN(I),CLH2O(I),CLNH4SH(I),
     +               CLNH3(I),CLH2S(I),CLCH4(I),CLAR(I)
      enddo
C  zpco etc are partial pressures
      do i=1,NL
          READ(2,13) ZTPR(I),ZPH2(I),ZPHE(I),ZPNH3(I),ZPH2O(I),
     +               ZPCH4(I),ZPH2S(I),ZPPH3(i),ZPCO(i),ZPCO13(i),ZPHCN(i)
      enddo
   23 FORMAT(10X,F10.4,8E13.8) 
   13 FORMAT(11E13.8) 
C 
	  READ(5,1060) NWN 
	  if(NOPR4.eq.1) WRITE(6,1060) NWN 
1060  FORMAT(I5)
C     The following could be changed in the future to a read vset(1) and dvset 
C     DO 31 I=1,NWN 
C  31 READ(5,1061) VSET(I) 
 1061 FORMAT(2F10.4)  
C     WILL PERFORM THE CALCULATIONS FOR THE WAVENUMBERS
	  read(5,1061) vset(1),vsetlast
	  dvset=(vset(1)-vsetlast)/nwn
	  do n=2,nwn
          vset(n)=vset(n-1)-dvset
      enddo   
C
      DO 200 N=1,NWN 
          VNU=VSET(N) 
          FREQ=C*VNU*1.0E9 
          VVNU=1.0/VNU
C      WRITE(6,9) VNU,FREQ,VVNU 
   9      FORMAT(/,/,2X,'V(CM-1),FREQ(HZ) ',1P,3(2X,E10.3)) 
C     INITIALIZE THE OPTICAL DEPTH TO ZERO. 
          DH2=0.0 
          DHE=0.0 
          DNH3=0.0 
          DH2O=0.0 
          DCH4=0.0 
	      DH2S=0.0
          dco=0.0
          dco13=0.0
          dhcn=0.0
          dpph3=0.0
          DOPTD=0.0 
C         SUM OVER THE VARIOUS PRESSURE LEVELS. 
          DO 100 M=1,NL2 
C             SPECIFY THE PRESSSURES OF H2,HE,NH3,H2O,AND CH4. 
              T=ZT(M)
              DELZ=ZDELZ(M)
              PH2=ZPH2(M)
              PHE=ZPHE(M)
              PNH3=ZPNH3(M)
              PH2O=ZPH2O(M)
              PCH4=ZPCH4(M)
	          PH2S=ZPH2S(M)
			  PPH3=ZPPH3(M)
              pco=zpco(M)
              pco13=zpco13(M)
              phcn=zphcn(M)
              TPR=ZTPR(M)
              T1=273./T 
              T2=300./T 
              T3=T2**C2 
              T4=T1**C1 
              T5=T1**C2 
              T6=1.0/(T**3.5) 
              T7=4.8/T
              ARGON=tpr*6.0E-14
              wtmol=PH2*2.0158 + PHE*4.0026 + PH2O*18.0153 + PCH4*16.04 + PNH3*17.03
     1            + PH2S*34.08 + ARGON*39.06 + PCO*28.01 + PHCN*27.018+ PPH3*33.9976
			  wtmol=wtmol*1.67e-24/tpr
C             CALL THE VARIOUS SUBROUTINES WHICH CALCULATE THE ABSORPTION COEFFICIENTS.                                                  
              IF(IQ.GE.1) CALL PART(T)                                                  
              IF(IQ.GE.1) FP=S1                                                         
			  IF(TPR.gt.1.0.and.PNH3.NE.0.0.and.vnu.lt.1.0) CALL NH3TOM(VNU,T,ANH3)
	          IF(TPR.le.1.0.and.PNH3.NE.0.0.and.vnu.lt.1.0) CALL NH3(VNU,T,ANH3)
 	          if(vnu.ge.1.0.and.PNH3.ne.0.0)  CALL NH3JOIN(VNU,T,ANH3) 
              IF(PNH3.EQ.0.0) ANH3=0.0
			  CALL H2O(VNU,T,AH2O)
C             AH2O=0.0
C             ANH3=0.0
C
C             H2H2 WITH MASSIE FORMALISM
C             H2H2JS WITH JOINER-STEFFES FORMALISM; not very good
C             H2H2OR  with Orton's tables;  this is best, but note the range inside which it works.
              if(vnu.le.frequ(1)) CALL H2H2(VNU,T,A,B)   
              if(T.le.Tmin) CALL H2H2(VNU,T,A,B)   
              if(T.ge.Tmax) CALL H2H2(VNU,T,A,B)   
              if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     +            CALL H2H2OR(IQ,VNU,T,A,B,CZB)   
	          CALL absbr_H2S2(ah2s,t,tpr,vnu)
C             Ah2S=0.0
c	          CALL absbr_ph3(aph3,t,tpr,vnu,ZPPH3)
C             CALCULATE THE OPTICAL DEPTH TAU 
C             ACO = absco(lco,fr0,freq,tpr,T,wtmol)
C             the line shape to be used is determined in absco routine.
              call absx_co(aco,t,tpr,vnu)
C             statement above is including all CO lines,  VVW only; this is  useful for continuum measurements
              ACO=ACO*pco/tpr
C             ****still an error in CO13. Also when simply running absx_co13 --  assumed to be continuum 
C             I need to test this. Note that fr0 is wrong for CO13 if the CO12 parameters are used
C             ACO13=0.0
C             ACO13 = absco13(lco,fr0,freq,tpr,T,wtmol)
C             call absx_co13(aco13,t,tpr,vnu)
C             ACO13=ACO13*pco13/tpr
c             write(6,1345) aco,pco,tpr
              A8=0.0
c	          goto 9191
              IF(MII.NE.1) GOTO 9191
c             H2O liquid and ice
			  IF(T.LT.273.0) THEN
			      CALL DIEL(2,0,VVNU,T,E1,E2,D1,D2)
                  X1=SIZE
                  X2=X1*2.0*3.1415/VVNU
                  CALL MIE(X2, D1, D2, QABS)
                  X=18*1.67E-24
                  Y=SIZE*SIZE*SIZE/8.0E-24
                  X=X*Y
                  XDEN=CLH2O(M)/X
                  A8=A8+XDEN*QABS*3.1415*DELZ*size*size
              ENDIF
              IF(T.GE.273.0) THEN
                  CALL DIEL(1,0,VVNU,T,E1,E2,D1,D2)
C                 X1 IS RADIUS PARTICLE
                  X1=SIZE
                  X2=X1*2.0*3.1415/VVNU
                  CALL MIE(X2, D1, D2, QABS)
C                 TAKE MOLEC.SIZE FEW ANGSTROMS; OR ITS RADIUS 2.00E-8 
C                 TAKE H2O-NH3 SOLUTION CLOUD WITH H2O-ICE /WATER PROPERTIES
C                 X=MASS GR/CM3; 
                  X=18*1.67E-24
                  Y=SIZE*SIZE*SIZE/8.0E-24
C                 WEIGHT OF ONE DROPLET IS X*Y
                  X=X*Y
C                 NUMBER OF DROPLETS PER CM**3:
                  XDEN=ZSOLN(M)/X
                  A8=A8+XDEN*QABS*3.1415*DELZ*size*size
			  ENDIF
c	          goto 9191
C             NH4SH CLOUD
              X1=SIZE
              X2=X1*2.0*3.1415/VVNU
              D1=1.7
              D2=.05
              CALL MIE(X2, D1, D2, QABS)
              X=51*1.67E-24
              Y=SIZE*SIZE*SIZE/8.0E-24
              X=X*Y
              XDEN=CLNH4SH(M)/X
              A8=A8+XDEN*QABS*3.1415*DELZ*size*size
	          GOTO 9191
C             NH3-ICE
              D1=1.3
              D2=.05
              X1=SIZE
              X2=X1*2.0*3.1415/VVNU
              CALL MIE(X2, D1, D2, QABS)
              X=17*1.67E-24
              Y=SIZE*SIZE*SIZE/8.0E-24
              X=X*Y
			  XDEN=CLNH3(M)/X
              A8=A8+XDEN*QABS*3.1415*DELZ*size*size
C             H2S-ICE
              D1=1.3
              D2=.01
              X1=SIZE
              X2=X1*2.0*3.1415/VVNU
              CALL MIE(X2, D1, D2, QABS)
              X=34*1.67E-24
              Y=SIZE*SIZE*SIZE/8.0E-24
              X=X*Y
              XDEN=CLH2S(M)/X
              A8=A8+XDEN*QABS*3.1415*DELZ*size*size
	          GOTO 9191
C             CALCULATE THE OPTICAL DEPTH TAU                                              
 9191         CONTINUE
C             FOR MASSIE FORMALISM USE FOLLOWING EXPRESSIONS
              DENH2=PH2/T
              DENHE=PHE/T 
              DENCH4=PCH4/T
              A2=DENH2*DENH2*A*5.246E-3*DELZ/2.9979 
              A3=DENH2*DENHE*B*5.246E-3*DELZ/2.9979 
              A12=0.
C             For Orton tables (best to use): amagat is p(bar)*269.6/T with 269.6=273.15/1.013
              if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     +             A2=DENH2*DENH2*273.*273.*delz*A/1.013/1.013
              if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     +             A3=DENH2*DENHE*273.*273.*delz*B/1.013/1.013
              if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     +             A12=DENH2*DENCH4*273.*273.*DELZ*CZB/1.013/1.013
C
              AH2HE=(A2+A3+A12) 
C             WRITE(6,1003) T,PH2,A,B,CZB,A2,A3,A12
C             write(6,1004) nfreq,ntemp
c 1003        FORMAT(8(1X,E10.3))
c 1004        format(2I10)
c             Joiner\Steffes formalism:
c	          AH2HE=AH2H2*DELZ
C             we absorb PH3 absorption into dh2s
              A4=ANH3*DELZ 
              A5=AH2O*DELZ 
	          A6=AH2S*DELZ
	          A9=APH3*DELZ
              DH2=DH2+A2 
              DHE=DHE+A3 
              DCH4=DCH4+A12
              DNH3=DNH3+A4 
              DH2O=DH2O+A5 
	          DH2S=DH2S+A6+A9
              DCO=DCO+ACO*delz
              DHCN=DHCN+AHCN*delz
              DCO13=DCO13+ACO13*delz
              A7=A4+A5+AH2HE+A8+A6+A9+ACO*DELZ
C             NH3, H2O, HCN and 13CO do not make a difference.
              DOPTD=DOPTD+A7 
C             TAU IS THE OPTICAL DEPTH. 
			  TAU(M)=DOPTD 
              GTAU(1,M)=DH2 
              GTAU(2,M)=DHE 
			  GTAU(3,M)=DNH3 
              GTAU(4,M)=DH2O 
              GTAU(5,M)=DH2S 
c             wt=a7/delz*exp(-tau(m))
c             write weighting function
c             WRITE(6,262) T,TPR,DELZ,wt
  262         FORMAT(4E15.5) 
  100     CONTINUE 
C         EVALUATE THE BRIGHTNESS TEMPERATURE OF THE CALCULATED INTENSITY. 
          DBT=0.0 
          DDR=0.0 
          DO 2001 JJ=1,18 
              CS=CSZA(JJ) 
              DC=DCSZA(JJ) 
              CALL BRIGHT(VNU,FREQ,NL,CS,BRI) 
C             CALL BRIGHTor(VNU,FREQ,NL,CS,BRI) 
              STR(JJ)=BRI	
              DBT=DBT+BRI*RR(JJ)*DC 
              DDR=DDR+RR(JJ)*DC 
2001      CONTINUE 
          DBT=DBT/DDR 
          WAVEL=1./VNU 
          GHZ=C*VNU
C         WRITE(6,225) VNU,WAVEL,GHZ 
225       FORMAT(2X,'V(CM-1)',4X,'WAVELENGTH(CM)',3X,'FREQ(GHZ)', 
     1           /,1P,E11.4,3X,E11.4,4X,E11.4) 
C         WRITE(6,226) 
226       FORMAT('       ') 
C         WRITE(6,235) wavel,STR(1),DBT
235       FORMAT(' wavel, wave^-1, center disk, DISK-AVERAGED TB',/,1P,6X,3E11.4) 
C 
          XWAVE(N)=WAVEL
          XTEMP(N)=DBT
          XNADIR(N)=STR(1)
  200 CONTINUE 
      DO 290 N=1, NWN
         x=1./xwave(N)
         GHZ=C*VSET(N)
         WRITE(6,5768) GHZ,XWAVE(N), x, XNADIR(N), XTEMP(N)
  290 CONTINUE
 5768 FORMAT(5F20.8)
      STOP 
      END
CENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND
CENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND
CENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND

      SUBROUTINE BRIGHT(VNU,FREQ,NL,CMU,BRI) 
      DIMENSION ZT(7000),ZDELZ(7000)
      DIMENSION TEM(7000),TOPR(7000) 
      DIMENSION TAU(7000),FINT(7000),GTAU(5,7000) 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4 
      COMMON/V10/TAU,FINT,GTAU 
      COMMON/ZP/ZT,ZDELZ 
C     BR IS THE INTENSITY AT THE TOP OF THE ATMOSPHERE. 
C     CMU IS THE COSINE OF THE SOLAR ZENITH ANGLE. 
C     FINT(I) IS THE FRACTION OF THE INTENSITY CONTRIBUTED BY LEVEL I. 
C     BR IS GIVEN BY EQUATION 90,PAGE 16 OF CHANDRASEKHAR'S 
C     RADIATIVE TRANSFER BOOK.  BR IS DERIVED FROM THE INTENSITY 
C     CALCULATED FROM EQUATION 90. EQUATION 90 ASSUMES THAT 
C     THE SOURCE FUNCTION IS SIMPLY GIVEN BY THE PLANCK FUNCTION. 
      BR=0.0 
      YY=FREQ
      A1=2.0*6.626E-9*YY*YY*YY/(C*C) 
      A2=6.626*1.0E-11*YY/1.3806 
C     WRITE(6,190) 
C190  FORMAT(2X,'LEVEL,H*FREQ/KT,EXP(-TAU),DTAU,B(T),',/, 
C    1  2X,'TAU,DINT,PARTIAL SUM OF I') 
      L1=NL 
      L2=L1-1 
      BR=0.0 
      DO 200 L=2,L2 
          K=L1+1-L 
          K2=K-1 
          S1=ZT(K) 
          S2=ZT(K2) 
          TAVE=0.5*(S1+S2) 
          A4=A2/TAVE 
          A5=(TAU(K)+TAU(K2))/(2.0*CMU) 
          A6=(TAU(K)-TAU(K2))/CMU 
          PL=A1/(EXP(A4)-1.00) 
          CALL EX(A5,A7) 
          DINT=PL*A7*A6 
          BR=BR+DINT 
          FINT(K)=DINT 
C         if(cmu.gt.0.99) WRITE(6,205) K,Tave,A7,A6,PL,TAU(K),DINT,BR 
 205      FORMAT(2X,I4,7(2X,E10.3)) 
 200  CONTINUE 
C     A8 IS THE BRIGHTNESS TEMPERATURE CORRESPONDING TO BR. 
      A7=A1/BR 
      A8=A2/(ALOG(1.00+A7)) 
C     WRITE(6,888) FREQ,CMU,BR,A8 
 888  FORMAT(2X,'FREQ,CMU,BR,T ARE ',/,4E15.5) 
      BRI=A8 
      RETURN 
      END
CENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND
      SUBROUTINE BRIGHTor(VNU,FREQ,NL,CMU,BRI) 
      real*8 BR,DINT,PL,A44,A55,A66,A77,A88,A99,A11,A22,YY
      DIMENSION ZT(7000),ZDELZ(7000)
      DIMENSION TEM(7000),TOPR(7000) 
      DIMENSION TAU(7000),FINT(7000),GTAU(5,7000) 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
      COMMON/V10/TAU,FINT,GTAU 
      COMMON/ZP/ZT,ZDELZ
C     BR IS THE INTENSITY AT THE TOP OF THE ATMOSPHERE. 
C     CMU IS THE COSINE OF THE SOLAR ZENITH ANGLE. 
C     FINT(I) IS THE FRACTION OF THE INTENSITY CONTRIBUTED BY LEVEL I. 
C     BR IS GIVEN BY EQUATION 90,PAGE 16 OF CHANDRASEKHAR'S 
C     RADIATIVE TRANSFER BOOK.  BR IS DERIVED FROM THE INTENSITY 
C     CALCULATED FROM EQUATION 90. EQUATION 90 ASSUMES THAT 
C     THE SOURCE FUNCTION IS SIMPLY GIVEN BY THE PLANCK FUNCTION. 
      BR=0.0 
      YY=FREQ
      A11=2.0*6.626E-9*YY*YY*YY/(C*C) 
      A22=6.626*1.0E-11*YY/1.3806 
C     WRITE(6,190) 
C190  FORMAT(2X,'LEVEL,H*FREQ/KT,EXP(-TAU),DTAU,B(T),',/, 
C    1  2X,'TAU,DINT,PARTIAL SUM OF I') 
      L1=NL 
      L2=L1-1 
      BR=0.0 
      K=NL 
      S1=ZT(K) 
      TAVE=S1
      A44=A22/TAVE 
      PL=A11/(EXP(A44)-1.00) 
C     determine BR at bottom atm
      BR=PL
C     if(cmu.gt.0.99) WRITE(6,204) K,Tave,PL,BR,tau(K)
 204  format(I5,4E15.5)
      DO 200 L=2,L2 
          K=L1+1-L 
          K2=K-1 
          S1=ZT(K) 
          S2=ZT(K2) 
          TAVE=0.5*(S1+S2) 
          A44=A22/TAVE 
C         A5=TAU(K2)/CMU 
          A55=(TAU(K)+TAU(K2))/(2.0*CMU) 
          A66=(TAU(K)-TAU(K2))/CMU 
          PL=A11/(EXP(A44)-1.00) 
          CALL EX(A66,A77) 
C         DINT=PL*A7*A6 
          BR=BR*A77+PL*(1.-A77) 
C         if(cmu.gt.0.99) WRITE(6,205) K,Tave,BR,A77,A66,PL,TAU(K),TAU(K2) 
 205      FORMAT(2X,I4,7(2X,E10.3)) 
 200  CONTINUE 
C     A8 IS THE BRIGHTNESS TEMPERATURE CORRESPONDING TO BR. 
      A99=A11/BR 
      A88=A22/(DLOG(A99+1.000)) 
C     WRITE(6,888) FREQ,CMU,BR,A8 
 888  FORMAT(2X,'FREQ,CMU,BR,T ARE ',/,4E15.5) 
      BRI=A88 
      RETURN 
      END 
CENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDENDEND
      SUBROUTINE ATMOS(NN,SOL) 
      DIMENSION XP(7000),XT(7000), XG(7000), H2S(7000), CH4(7000),
     +          co(7000),co13(7000),hcn(7000)
      DIMENSION TEMPA(7000),TEMPW(7000),TEMPP(7000),TEMPHY(7000)
      DIMENSION TEMPHE(7000),TEMPXG(7000),HELIUM(7000),XH(7000)
	  DIMENSION TEMPT(7000),TMPCH4(7000),TMPH2S(7000),TSOLN(7000),
     +          tmpph3(7000),tmpco(7000),tmphcn(7000)
      DIMENSION ZDELZ(7000),TDELZ(7000),SOLN(7000),fosf(7000)
      DIMENSION AMMON(7000),WATER(7000),DRUK(7000),HYDR(7000),fosfien(7000)
      DIMENSION TC(45),CPH2D(45) 
	  DIMENSION TCLH2O(7000),TCLNH4SH(7000),TCLNH3(7000),TCLH2S(7000)
	  DIMENSION TCLCH4(7000),TCLAR(7000)
	  DIMENSION CLH2O(7000),CLNH4SH(7000),CLNH3(7000),CLH2S(7000)
	  DIMENSION CLCH4(7000),CLAR(7000)
	  COMMON/CL/CLH2O,CLNH4SH,CLNH3,CLH2S,CLCH4,CLAR
	  COMMON/U2/C1,C2,C3,C4 
      COMMON/U5/SP,SL,SLAM,SLH2O,KB,KS 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
	  COMMON/ABUN/AB1,AB2,AB3,AB4,AB5,AB6 
      COMMON/ATM/TC,CPH2D 
C 
C     AMMON,WATER MIXING RATIOS NH3,H2O 
C     THEY ARE READ IN AT PRESSURE LEVELS OF 
C     .01,.02,.03,.05,.07,.1,.2,.3,.5,.7,1,2,3,5,7,10  BARS. 
      TFLAT=1.0E10
      KOUNT=0
      ZZZ=0.0 
      READ(1,19) JMAX
  19  FORMAT(1x,I5)
	  write(6,19) JMAX
      SOL=0.0
      ZD=0.0
      zdelz(1)=0.0
      DO 8079 II=1, JMAX
          READ(1,*) XH(II),XT(II),XP(II),HYDR(II),HELIUM(II),CH4(II),
     +              AMMON(II),WATER(II),H2S(II),SOLN(II)
C         read(1,*) s,XH(II),XP(II),XT(II)
C         Hydr(ii)=0.831
C         helium(ii)=0.149
C         h2s(ii)=0.0
C         water(ii)=2.0E-4
C         ammon(ii)=1.0E-4
C         ch4(ii)=2.0E-2
          co(ii)=1.0E-6
          if(xp(ii).ge.0.1585) co(ii)=0.0
C         ***NOTE: check on pco in lats loop; do 50 M=1,----
          co13(ii)=co(ii)*1.0E-2
          hcn(ii)=0.0
          soln(ii)=0.0
 222      FORMAT(F5.1,2X,F6.2,E13.3,2F7.3,5(E10.3))
 223      FORMAT(I5,F7.2,E11.3,2F7.3,4(E10.3))
c         IF(II.GT.1) ZDELZ(II)=(XH(II-1)-XH(II))*1.0E5*900./1130.
          IF(II.GT.1) ZDELZ(II)=abs(XH(II-1)-XH(II))*1.0E5
c	      SOLN(II)=SOLN(II)/10.0
          IF(II.NE.1) SOL=SOL+SOLN(II)*ZDELZ(II)
          IF(SOLN(II).GT.0.0) ZD=ZD+ZDELZ(II)
C         IF(XT(II).GT.400.0) ammon(ii)=1.94E-4 
	      fosfie=0.0
	      fosfien(ii)=fosfie
C         CALCULATE VAPOR PRESSURE FOR NH3
C         Skip steps if Paul's files are used
	  goto 8079
	  if(ammon(ii).gt.3.0E-6) ammon(ii)=3.0E-6
      IF(XT(II).GT.195.0) GOTO 608
      PNH3=1.342E7*EXP(-3753.6/XT(II)) 
      PRNH3=PNH3/XP(II)
      IF(PRNH3.GE.AMMON(II)) PNH3=AMMON(II)*XP(II) 
      AMMON(II)=PNH3/XP(II)
  608 CONTINUE
C     calculate vapor pressure for water
C     skip step if Paul's profiles are used
      goto 8079
      IF(XT(II).le.273.) goto 9910
      IF(XT(II).gt.273.) goto 9920
 9910 PH2O=1.013*4.327E7*EXP(-6194./XT(II)) 
      PRH2O=PH2O/XP(II)             
      GOTO 9921                     
 9920 PH2O=1.013*1.5350E6*EXP(-5273./XT(II))   
      PRH2O=PH2O/XP(II)          
 9921 IF(PRH2O.GE.water(ii)) PH2O=water(ii)*xp(ii)
      water(ii)=PH2O/xp(ii)  
c     TT=XT(ii)*XT(ii)
c     VPCH4S = 4.425070D0 - 453.92414D0/XT(ii) - 4055.6016D0/TT
c    1   + 115352.19D0/XT(ii)**3 - 1165560.7D0/XT(ii)**4
c      PCH4 = (10.0D0**VPCH4S) * 1.013
c      PRCH4=PCH4/XP(II)
c      if(PRCH4.ge.CH4(ii)) PCH4=ch4(II)*XP(II)
C      if(xt(ii).lt.80.) CH4(ii)=0.0
 8079  continue
C  skip step if Paul's profiles are used
       goto 8069
C      
       xtabun=ammon(jmax)+water(jmax)+H2S(jmax)+ch4(jmax)+co(jmax)
       fhe=helium(jmax)/hydr(jmax)
       fxtab=xtabun/hydr(jmax)
       DO 8089 II=1, JMAX
       xabun=ch4(ii)+ammon(II)+water(ii)+H2S(ii)+co(ii)
       fac=(xabun/xtabun)*fxtab
       hydr(ii)=1.0/(1.0+fhe+fac)
       helium(ii)=hydr(ii)*fhe
 8089  continue
 8069  continue
C
        DO i=1,jmax
            j=jmax+1-i
c           read(10,*) x,y,z,soln(j),clh2o(j),clnh4sh(j),
c    +                 clnh3(j),clh2s(j),clch4(j),clar(j)
            soln(j)=0.
            clh2o(j)=0.
            clnh4sh(j)=0.
            clnh3(j)=0.
            clh2s(j)=0.
            clch4(j)=0.
            clar(j)=0.
        enddo
        ZDELZ(1)=ZDELZ(2)
c       SOL=SOL/ZD
C       WRITE(6,994) SOL,ZD
  994   FORMAT('DENSITY AND ALTITUDE RANGE ARE',2E12.4)
        TT=XT(1)
        PR=XP(1)
        JM=JMAX-1
C       skip all te rest?
C       goto 6904
        DO 601 II=1,JM
            X=ALOG(XT(II+1)/XT(II))
		    Y=ALOG(XP(II+1)/XP(II))
            XG(II)=X/Y
  601   CONTINUE
	    write(6,*) xg(jm)
        TEMPT(1)=XT(1)
        TEMPP(1)=XP(1)
        TEMPHY(1)=HYDR(1)
        TEMPHE(1)=HELIUM(1)
        TMPCH4(1)=CH4(1)
        TEMPA(1)=AMMON(1)
        TEMPW(1)=WATER(1)
        TMPph3(1)=fosfien(1)
        TMPco(1)=co(1)
        tmphcn(1)=hcn(1)
        TMPH2S(1)=H2S(1)
        TEMPXG(1)=XG(1)
        TDELZ(1)=ZDELZ(1)
        TSOLN(1)=SOLN(1)
        TCLH2O(1)=CLH2O(1)
        TCLNH4SH(1)=CLNH4SH(1)
        TCLNH3(1)=CLNH3(1)
        TCLH2S(1)=CLH2S(1)
        TCLCH4(1)=CLCH4(1)
        TCLAR(1)=CLAR(1)
		DO 387 IJ=1,2
	        KMAX=JMAX
			JJ=1
	        write(6,*) jmax,jm
            DO 102 II=2,JMAX
                DO 101 KK=1,3
                    JJ=JJ+1
                    TEMPXG(JJ)=XG(II)
                    TEMPT(JJ)=KK*(XT(II)-XT(II-1))/3.0 + XT(II-1)
				    TEMPP(JJ)=KK*(XP(II)-XP(II-1))/3.0 + XP(II-1)
                    TEMPHY(JJ)=KK*(HYDR(II)-HYDR(II-1))/3.0 + HYDR(II-1)
                    TEMPHE(JJ)=KK*(HELIUM(II)-HELIUM(II-1))/3.0 + HELIUM(II-1)
                    TMPCH4(JJ)=KK*(CH4(II)-CH4(II-1))/3.0 + CH4(II-1)
                    TEMPA(JJ)=KK*(AMMON(II)-AMMON(II-1))/3.0 + AMMON(II-1)
				    TEMPW(JJ)=KK*(WATER(II)-WATER(II-1))/3.0 + WATER(II-1)
				    TMPH2S(JJ)=KK*(H2S(II)-H2S(II-1))/3.0 + H2S(II-1)
                    TMPco(JJ)=KK*(co(II)-co(II-1))/3.0 + co(II-1)
				    TMPHcn(JJ)=KK*(Hcn(II)-hcn(II-1))/3.0 + hcn(II-1)
                    TMPph3(JJ)=KK*(fosfien(II)-fosfien(II-1))/3.0 + fosfien(II-1)
				    TDELZ(JJ)=ZDELZ(II)/3.0
                    TSOLN(JJ)=kk*(SOLN(II)-soln(ii-1))/3.0 + soln(ii-1)
                    if(soln(ii).eq.0..and.soln(ii-1).gt.0.)
     +                  tsoln(jj)=(4-kk)*soln(ii-1)/3.0
				    TCLH2O(JJ)=kk*(CLH2O(II)-clh2o(ii-1))/3.0 + clh2o(ii-1)
                    if(clh2o(ii).eq.0..and.clh2o(ii-1).gt.0.) 
     +                  tclh2o(jj)=(4-kk)*clh2o(ii-1)/3.0
                    TCLNH4SH(JJ)=kk*(CLNH4SH(II)-clnh4sh(ii-1))/3.0 + clnh4sh(ii-1)
                    if(clnh4sh(ii).eq.0..and.clnh4sh(ii-1).gt.0.) 
     +                  tclnh4sh(jj)=(4-kk)*clnh4sh(ii-1)/3.0
                    TCLNH3(JJ)=CLNH3(II)
                    TCLH2S(JJ)=CLH2S(II)
                    TCLCH4(JJ)=CLCH4(II)
                    TCLAR(JJ)=CLAR(II)
  101           CONTINUE
                KMAX=KMAX+2
  376           CONTINUE
  102       CONTINUE
			JMAX=KMAX
            WRITE(6,19) JMAX
            DO 103 II=1,JMAX
                ZDELZ(II)=TDELZ(II)
		        SOLN(II)=TSOLN(II)
                CLH2O(II)=TCLH2O(II)
                CLNH4SH(II)=TCLNH4SH(II)
                CLNH3(II)=TCLNH3(II)
                CLH2S(II)=TCLH2S(II)
                CLCH4(II)=TCLCH4(II)
                CLAR(II)=TCLAR(II)
                XT(II)=TEMPT(II)
                XP(II)=TEMPP(II)
                HYDR(II)=TEMPHY(II)
                HELIUM(II)=TEMPHE(II)
                CH4(II)=TMPCH4(II)
                AMMON(II)=TEMPA(II)
                WATER(II)=TEMPW(II)
                H2S(II)=TMPH2S(II)     
                Hcn(II)=TMPHcn(II)
                co(II)=TMPco(II)
                fosfien(ii)=tmpph3(ii)
                XG(II)=TEMPXG(II)
  388           CONTINUE
  103       CONTINUE
  387   CONTINUE
C  if skipping extension of the adiabat, then (there is an error in the below, and it is better to use Paul's code)
 6904  continue
       KOUNT=0
CHEREIAMHEREIAM
       DO 50 M=1,JMAX
C   DELP IS THE DIFFERENCE IN THE PRESSURE LEVELS OF THE 
C   MODEL ATMOSPHERE, IN LOG10 COORDINATES. 
      PO=PR 
      TO=TT 
      IF(M.EQ.1) GOTO 50 
      PR=XP(M)
      TT=XT(M)
      PCH4=CH4(M)*PR
      H2OMR=WATER(M)
      AMR=AMMON(M)
      PHE=HELIUM(M)*PR
      PH2S=H2S(M)*PR
      PNH3=AMR*PR
      PH2O=H2OMR*PR
      PH2=HYDR(M)*PR 
	PPh3=fosfien(m)*pr
        if(pr.ge.0.1585) co(m)=0.0
        pco=co(m)*pr
        pco13=co13(m)*pr
        phcn=hcn(m)*pr
C 
C  COMPUTE DELZ, DIFFERENCE IN HEIGHT BETWEEN THE TWO 
C  PRESSURE LEVELS PO AND PR 
C 
      DELZ=ZDELZ(M)
      ZZZ=ZZZ+DELZ*1.0E-5 
      TPR=PR
      WRITE(2,7) TPR,PH2,PHE,PNH3,PH2O,PCH4,PH2S,PCO,PCO13,PHCN
      IF(TT.GT.TFLAT) 
     +WRITE(3,8) TPR,TFLAT,DELZ,SOLN(M),CLH2O(M),CLNH4SH(M),
     + CLNH3(M),CLH2S(M),CLCH4(M),CLAR(M)
      IF(TT.LE.TFLAT) 
     +WRITE(3,8) TPR,TT,DELZ,SOLN(M),CLH2O(M),CLNH4SH(M),
     + CLNH3(M),CLH2S(M),CLCH4(M),CLAR(M)
    7 FORMAT(10E13.8) 
    8 FORMAT(E10.5,F10.4,8E13.8) 
       KOUNT=KOUNT+1 
 50   CONTINUE 
 9999 NN=KOUNT 
      M5=0 
 1049 FORMAT(F10.4) 
  110 FORMAT(I1,1X,I4,I4) 
      RETURN 
      END 
	  
	  
	  
	  
	  
      SUBROUTINE MIE(X,A,B,QABS)
      COMPLEX ZM
      ZM=CMPLX(A,-B)
      ZM=4.0*X*((ZM*ZM-1.0)/(ZM*ZM+2.0))
      QABS=-AIMAG(ZM) 
C      A2=2.0*A*B
C      A3=A2*A2
C      A4=(A*A)-(B*B)+2.0
C      A5=A4*A4
C      QABS=4.0*X*((3.0*B)/(A3+A5))
      RETURN
      END
C
C
      SUBROUTINE EX(TX,S)                                                       
C     SINCE THE COMPUTER DOES NOT LIKE TOO LARGE OR SMALL ARGUMENTS             
C     FOR THE EXPONENTIAL FUNCTION, WE SET THE VALUE TO ZERO IF                 
C     NEED BE.                                                                  
      IF (TX .GT. 75.0) S=0.0                                                   
      IF (TX .LT. 75.0) S=EXP(-TX)                                              
      RETURN                                                                    
      END                                                                       
      BLOCK DATA                                                                
C     CONSTANTS USED IN THE PROGRAM ARE SPECIFIED HERE.                         
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/U2/C1,C2,C3,C4                                                     
C     PI,SPEED OF LIGHT,PLANCK'S CONST/TWO PI,TWO PI C,                         
C     BOLTZMANN'S CONST,AW,1 AMAGAT=2.69E19 CM-3.                               
C     AW IS 2*PI*PI/3*HB*C                                                      
      DATA PI,C,HB,TPC,BOLC,AW,DEN/3.14159,2.99792458E1,1.054,                 
     1  1.884,1.380622,2.082,2.687E19/                                
C     BAR IS 1.0E6 DYNES/CM2.                                                   
      DATA BAR/1.0E6/                                                           
C     THESE MOLECULAR CONSTANTS FOR NH3 ARE IN MHZ.                             
C     C1,C2,C3,C4 ARE 13/3,2/3,1/3,30*760*SQRT(293)                             
      DATA C1,C2,C3,C4/4.333,0.666,0.333,3.902E+05/                             
      END 
C
      REMOVED ABSORPTION
C                                                                               
C                                                                               
C                                           C                                                                               
C                                                                               
