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
      dimension Ttab(10),frequ(2428)
      dimension abeh2h2(2428,10),abnh2h2(2428,10)
      dimension abeh2he(2428,10),abnh2he(2428,10)
      dimension abeh2ch4(2428,10),abnh2ch4(2428,10)
      DIMENSION VKAK(16,16),DVKAK(16,16) 
      DIMENSION TC(45),CPH2D(45) 
      DIMENSION RR(18),CSZA(18),DCSZA(18),STR(19)
      DIMENSION ZT(7000),ZDELZ(7000),ZTPR(7000),ZPH2(7000),ZPH2S(7000),
     + ZPHE(7000),ZPNH3(7000),ZPH2O(7000),ZPCH4(7000),ZSOLN(7000),
     + ZPCO(7000),ZPCO13(7000),ZPHCN(7000),ZPPH3(7000)
      DIMENSION XWAVE(30000),XTEMP(30000), XNADIR(30000),VSET(30000) 
      DIMENSION TAU(7000),FINT(7000),GTAU(5,7000) 
	DIMENSION CLH2O(7000),CLNH4SH(7000),CLNH3(7000),CLH2S(7000)
	DIMENSION CLCH4(7000),CLAR(7000)
	real*4  F0(311),S0(311),EL(311),aph3,wgths(311)
	real*4  FP0(728),SP0(728),ELP(728)
	REAL*4 wgtS0(40),wgtFGB(40),wgtSB(40)
	real*4 t3,gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,
     + deltot,PCO,PCO13,PHCN
      COMMON/orton/Ttab,frequ,abeh2h2,abnh2h2,abeh2he,abnh2he,
     + abeh2ch4,abnh2ch4,nfreq,ntemp
       COMMON/ph3/FP0,SP0,ELP,wgtS0,wgtFGB,wgtSB
       COMMON/h2s/F0,S0,EL,wgths
       COMMON/CL/CLH2O,CLNH4SH,CLNH3,CLH2S,CLCH4,CLAR
       COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
      COMMON/U2/C1,C2,C3,C4 
      COMMON/U3/VKAK,DVKAK 
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7 
      COMMON/U5/SP,SL,SLAM,SLH2O,KB,KS 
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
      COMMON/V3/QROT,QR1,QR2,S1,S2,S3,S4,P1,P2,P3,P4 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4 
      COMMON/V6/P5,P6,P7,P8 
      COMMON/V7/PRO1(20),PRO2(20),STOR(20) 
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2 
      COMMON/V10/TAU,FINT,GTAU 
      COMMON/W1/FP,REB,RHO(20) 
      COMMON/ABUN/AB1,AB2,AB3,AB4,AB5,AB6 
      COMMON/ATM/TC,CPH2D 
      COMMON/ZP/ZT,ZDELZ,ZTPR,ZPH2,ZPHE,ZPNH3,ZPH2O,ZPCH4,ZPH2S,
     +ZPCO,ZPCO13,ZPHCN,ZPPH3
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
C      OPEN(6,NAME='nepco_spectotal_allgases_dry.dat')
      OPEN(6,NAME='test.dat')
      Open(7,NAME='kakar.dat')
      open(11,NAME='orton_H2.tables',form='formatted')
C     TO SUPPRESS THE LENGHTLY PRINT OUT OF INTERMEDIATE VALUES, 
C     SET NOPR AND NOPR2 TO 0 (NOT EQUAL TO 1). 
C     NOPR3 CONTROLS PRINT OUT OF MIXING RATIOS 
C  DEN=INTEGRATED CLOUD DENSITY; ALT=ALTITUDE RANGE OF CLOUD
C  for some reason I thought I calculated H2-CH4 I guess, 
C but this hasn't been done
      CZB=0.0
      DEN=7.2e6
      NOPR=0 
      NOPR2=0 
      NOPR3=0
      NOPR4=1
C
      read (11,*) ntemp,Tmax,Tmin
      dltemp=(alog(Tmax)-alog(Tmin))/float(ntemp-1)
      do n=1,ntemp
         Ttab(n)=alog(Tmin)+(n-1)*dltemp
      enddo
C
      read (11,*) nfreq
      read (11,*) (frequ(i),i=1,nfreq)
      do n=1,nfreq
         read (11,*) (abeh2h2(n,i),i=1,ntemp)
      enddo
      do n=1,nfreq
         read (11,*) (abnh2h2(n,i),i=1,ntemp)
      enddo
      do n=1,nfreq
         read (11,*) (abeh2he(n,i),i=1,ntemp)
      enddo
      do n=1,nfreq
         read (11,*) (abnh2he(n,i),i=1,ntemp)
      enddo
      do n=1,nfreq
         read (11,*) (abeh2ch4(n,i),i=1,ntemp)
      enddo
      do n=1,nfreq
         read (11,*) (abnh2ch4(n,i),i=1,ntemp)
      enddo
C 
      READ(5,779) MII
      READ(5,779) KQQ 
  779 FORMAT(I2) 
      READ(5,773) SIZE
C      WRITE(6,772) SIZE
      SIZE=SIZE/2.0
  772 FORMAT(2X,'WATER DROPLET SIZE',2X,E10.5,/)
 773  FORMAT(F10.4)
C  EQUIL: IQ=2, FP=S1, REB=1.0
C  NORMAL: IQ=0; FP=0.25; REB=0
C  INTERM: IQ=1; REB=1.0.; FP=S1
C  IQ=2: TEMP FOR EQUIL. CASE
      READ (5,778) IQ,REB,FP 
      if(NOPR4.eq.1) WRITE(6,778) IQ,REB,FP 
  778 FORMAT(I2,2F10.4) 
C 
      DO 1357 I=1,18 
      RR(I)=ACOS(CSZA(I)) 
      RR(I)=SIN(RR(I)) 
       if(NOPR4.eq.1) WRITE(6,1345) CSZA(I),DCSZA(I),RR(I) 
 1345 FORMAT(3E15.5) 
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
C  foe equilibrium H2 
      IF(IQ.EQ.2) READ(7,103) (TC(I),CPH2D(I),I=1,24) 
C      WRITE(6,105) (TC(I),CPH2D(I),I=1,44)
	close(7)
C   read in H2S lines  THIS IS OLD
c	OPEN(7,NAME='h2s_file.list')
c	read(7,7465) (F0(i),S0(i),EL(i),i=1,311)
c 7465 format(F13.4,9x,F7.4,3x,F9.4,32x)
      OPEN(7,NAME='H2S_deboer.LIN')
	read(7,*) (F0(i),S0(i),EL(i),wgths(i),i=1,311)
	close(7)
      OPEN(7,NAME='ph3_file.list')
	read(7,*) (FP0(i),SP0(i),ELP(i),i=1,728)
	close(7)
      OPEN(7,NAME='ph3wgt.dat')
	read(7,*) (k,wgtS0(i),wgtFGB(i),wgtSB(i),i=1,40)
	close(7)
C 
 4444 FORMAT(I5,I5,2X,F15.4,2X,F15.4) 
C 
      READ(5,1041) M5,NL
 1041 FORMAT(I1,1X,I4) 
      if(NOPR4.eq.1) WRITE(6,1046) M5,NL 
 1046 FORMAT(2X,'PARAMETER M5 IS EQUAL TO',/,2I5) 
 1049 FORMAT(F10.4) 
C      G=1130.0 
      G=1140.0 
      GE=G 
      if(NOPR4.eq.1) WRITE(6,1050) G 
 1050 FORMAT(2X,'THE GRAVITY G IS ',/,F12.4) 
C      READ(5,1045) AB1,AB2,AB3,AB4,AB5,AB6 
 1045 FORMAT(6E10.3) 
C      WRITE(6,1047) AB1,AB2,AB3,AB4,AB5,AB6 
 1047 FORMAT(2X,'ABUNDANCES H2,HE,NH3,H2O,CH4,H2S ARE',/,6E12.5) 
      IF(M5.EQ.0) CALL ATMOS(NL,SOL) 
      if(NOPR4.eq.1) WRITE(6,1051) NL
 1051  FORMAT(2X,'THE PARAMETER NL AFTER ATMOS IS EQUAL TO',/,I6)
      NL2=NL-1 
      REWIND 2
      REWIND 3
      do i=1,NL
      READ(3,23) ZT(I),ZDELZ(I),ZSOLN(I),CLH2O(I),CLNH4SH(I),
     + CLNH3(I),CLH2S(I),CLCH4(I),CLAR(I)
      enddo
C  zpco etc are partial pressures
      do i=1,NL
      READ(2,13) ZTPR(I),ZPH2(I),ZPHE(I),ZPNH3(I),ZPH2O(I),
     +ZPCH4(I),ZPH2S(I),ZPPH3(i),ZPCO(i),ZPCO13(i),ZPHCN(i)
      enddo
 23   FORMAT(10X,F10.4,8E13.8) 
   13 FORMAT(11E13.8) 
C 
       READ(5,1060) NWN 
       if(NOPR4.eq.1) WRITE(6,1060) NWN 
1060   FORMAT(I5)
C The following could be changed in the future to a read vset(1) and dvset 
C       DO 31 I=1,NWN 
C  31  READ(5,1061) VSET(I) 
 1061 FORMAT(2F10.4)  
C  WILL PERFORM THE CALCULATIONS FOR THE WAVENUMBERS
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
   9  FORMAT(/,/,2X,'V(CM-1),FREQ(HZ) ',1P,3(2X,E10.3)) 
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
C     SUM OVER THE VARIOUS PRESSURE LEVELS. 
      DO 100 M=1,NL2 
C     SPECIFY THE PRESSSURES OF H2,HE,NH3,H2O,AND CH4. 
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
     1  + PH2S*34.08 + ARGON*39.06 + PCO*28.01 + PHCN*27.018+ PPH3*33.9976
       wtmol=wtmol*1.67e-24/tpr
C     CALL THE VARIOUS SUBROUTINES WHICH CALCULATE THE                          
C     ABSORPTION COEFFICIENTS.                                                  
      IF(IQ.GE.1) CALL PART(T)                                                  
      IF(IQ.GE.1) FP=S1                                                         
      IF(TPR.gt.1.0.and.PNH3.NE.0.0.and.vnu.lt.1.0) CALL NH3TOM(VNU,T,ANH3)
	IF(TPR.le.1.0.and.PNH3.NE.0.0.and.vnu.lt.1.0) CALL NH3(VNU,T,ANH3)
 	if(vnu.ge.1.0.and.PNH3.ne.0.0)  CALL NH3JOIN(VNU,T,ANH3) 
       IF(PNH3.EQ.0.0) ANH3=0.0
       CALL H2O(VNU,T,AH2O)
C       AH2O=0.0
C       ANH3=0.0
C
C  H2H2 WITH MASSIE FORMALISM
C H2H2JS WITH JOINER-STEFFES FORMALISM; not very good
C  H2H2OR  with Orton's tables;  this is best, but note the range inside which it works.
       if(vnu.le.frequ(1)) CALL H2H2(VNU,T,A,B)   
       if(T.le.Tmin) CALL H2H2(VNU,T,A,B)   
       if(T.ge.Tmax) CALL H2H2(VNU,T,A,B)   
       if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     +  CALL H2H2OR(IQ,VNU,T,A,B,CZB)   
	CALL absbr_H2S2(ah2s,t,tpr,vnu)
C        Ah2S=0.0
c	CALL absbr_ph3(aph3,t,tpr,vnu,ZPPH3)
C  CALCULATE THE OPTICAL DEPTH TAU 
C          ACO = absco(lco,fr0,freq,tpr,T,wtmol)
C  the line shape to be used is determined in absco routine.
         call absx_co(aco,t,tpr,vnu)
C statement above is including all CO lines,  VVW only; this is  useful for continuum measurements
           ACO=ACO*pco/tpr
C****still an error in CO13. Also when simply running absx_co13 --  assumed to be continuum 
C  I need to test this. Note that fr0 is wrong for CO13 if the CO12 parameters are used
C           ACO13=0.0
C          ACO13 = absco13(lco,fr0,freq,tpr,T,wtmol)
C          call absx_co13(aco13,t,tpr,vnu)
C          ACO13=ACO13*pco13/tpr
c       write(6,1345) aco,pco,tpr
        A8=0.0
c	goto 9191
      IF(MII.NE.1) GOTO 9191
c  H2O liquid and ice
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
C  X1 IS RADIUS PARTICLE
C  X1 IS RADIUS PARTICLE
      X1=SIZE
      X2=X1*2.0*3.1415/VVNU
      CALL MIE(X2, D1, D2, QABS)
C  TAKE MOLEC.SIZE FEW ANGSTROMS; OR ITS RADIUS 2.00E-8 
C  TAKE H2O-NH3 SOLUTION CLOUD WITH H2O-ICE /WATER PROPERTIES
C  X=MASS GR/CM3; 
      X=18*1.67E-24
      Y=SIZE*SIZE*SIZE/8.0E-24
C  WEIGHT OF ONE DROPLET IS X*Y
      X=X*Y
C  NUMBER OF DROPLETS PER CM**3:
      XDEN=ZSOLN(M)/X
      A8=A8+XDEN*QABS*3.1415*DELZ*size*size
        ENDIF
c	goto 9191
C  NH4SH CLOUD
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
C  NH3-ICE
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
C  H2S-ICE
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
C  CALCULATE THE OPTICAL DEPTH TAU                                              
C
 9191 CONTINUE
C  FOR MASSIE FORMALISM USE FOLLOWING EXPRESSIONS
C
      DENH2=PH2/T
      DENHE=PHE/T 
      DENCH4=PCH4/T
      A2=DENH2*DENH2*A*5.246E-3*DELZ/2.9979 
      A3=DENH2*DENHE*B*5.246E-3*DELZ/2.9979 
      A12=0.
C For Orton tables (best to use): amagat is p(bar)*269.6/T with 269.6=273.15/1.013
             if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     + A2=DENH2*DENH2*273.*273.*delz*A/1.013/1.013
             if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     + A3=DENH2*DENHE*273.*273.*delz*B/1.013/1.013
             if(T.gt.Tmin.and.T.lt.Tmax.and.vnu.gt.frequ(1)) 
     + A12=DENH2*DENCH4*273.*273.*DELZ*CZB/1.013/1.013
C
      AH2HE=(A2+A3+A12) 
C       WRITE(6,1003) T,PH2,A,B,CZB,A2,A3,A12
C       write(6,1004) nfreq,ntemp
c 1003  FORMAT(8(1X,E10.3))
c 1004  format(2I10)
c  Joiner\Steffes formalism:
c	AH2HE=AH2H2*DELZ
C  we absorb PH3 absorption into dh2s
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
C NH3, H2O, HCN and 13CO do not make a difference.
      DOPTD=DOPTD+A7 
C     TAU IS THE OPTICAL DEPTH. 
      TAU(M)=DOPTD 
      GTAU(1,M)=DH2 
      GTAU(2,M)=DHE 
      GTAU(3,M)=DNH3 
      GTAU(4,M)=DH2O 
      GTAU(5,M)=DH2S 
c       wt=a7/delz*exp(-tau(m))
c  write weighting function
c      WRITE(6,262) T,TPR,DELZ,wt
  262 FORMAT(4E15.5) 
  100 CONTINUE 
C 
C 
C     EVALUATE THE BRIGHTNESS TEMPERATURE OF THE CALCULATED 
C     INTENSITY. 
      DBT=0.0 
      DDR=0.0 
      DO 2001 JJ=1,18 
      CS=CSZA(JJ) 
      DC=DCSZA(JJ) 
        CALL BRIGHT(VNU,FREQ,NL,CS,BRI) 
C       CALL BRIGHTor(VNU,FREQ,NL,CS,BRI) 
      STR(JJ)=BRI	
      DBT=DBT+BRI*RR(JJ)*DC 
      DDR=DDR+RR(JJ)*DC 
2001  CONTINUE 
      DBT=DBT/DDR 
      WAVEL=1./VNU 
      GHZ=C*VNU
C      WRITE(6,225) VNU,WAVEL,GHZ 
225   FORMAT(2X,'V(CM-1)',4X,'WAVELENGTH(CM)',3X,'FREQ(GHZ)', 
     1  /,1P,E11.4,3X,E11.4,4X,E11.4) 
C      WRITE(6,226) 
226   FORMAT('       ') 
C      WRITE(6,235) wavel,STR(1),DBT
235   FORMAT(' wavel, wave^-1, center disk, DISK-AVERAGED TB',/,1P,6X,3E11.4) 
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
 5768   FORMAT(5F20.8)
      STOP 
      END 
      SUBROUTINE BRIGHT(VNU,FREQ,NL,CMU,BRI) 
      DIMENSION ZT(7000),ZDELZ(7000)
      DIMENSION TEM(7000),TOPR(7000) 
      DIMENSION TAU(7000),FINT(7000),GTAU(5,7000) 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4 
      COMMON/V10/TAU,FINT,GTAU 
      COMMON/ZP/ZT,ZDELZ
C
C 
C     BR IS THE INTENSITY AT THE TOP OF THE ATMOSPHERE. 
C     CMU IS THE COSINE OF THE SOLAR ZENITH ANGLE. 
C     FINT(I) IS THE FRACTION OF THE INTENSITY CONTRIBUTED 
C     BY LEVEL I. 
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
C      if(cmu.gt.0.99) WRITE(6,205) K,Tave,A7,A6,PL,TAU(K),DINT,BR 
 205  FORMAT(2X,I4,7(2X,E10.3)) 
 200  CONTINUE 
C     A8 IS THE BRIGHTNESS TEMPERATURE CORRESPONDING TO BR. 
      A7=A1/BR 
      A8=A2/(ALOG(1.00+A7)) 
C     WRITE(6,888) FREQ,CMU,BR,A8 
 888  FORMAT(2X,'FREQ,CMU,BR,T ARE ',/,4E15.5) 
      BRI=A8 
      RETURN 
      END
C
      SUBROUTINE BRIGHTor(VNU,FREQ,NL,CMU,BRI) 
      real*8 BR,DINT,PL,A44,A55,A66,A77,A88,A99,A11,A22,YY
      DIMENSION ZT(7000),ZDELZ(7000)
      DIMENSION TEM(7000),TOPR(7000) 
      DIMENSION TAU(7000),FINT(7000),GTAU(5,7000) 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
      COMMON/V10/TAU,FINT,GTAU 
      COMMON/ZP/ZT,ZDELZ
C
C 
C     BR IS THE INTENSITY AT THE TOP OF THE ATMOSPHERE. 
C     CMU IS THE COSINE OF THE SOLAR ZENITH ANGLE. 
C     FINT(I) IS THE FRACTION OF THE INTENSITY CONTRIBUTED 
C     BY LEVEL I. 
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
C      determine BR at bottom atm
      BR=PL
C      if(cmu.gt.0.99) WRITE(6,204) K,Tave,PL,BR,tau(K)
 204  format(I5,4E15.5)
      DO 200 L=2,L2 
      K=L1+1-L 
      K2=K-1 
      S1=ZT(K) 
      S2=ZT(K2) 
      TAVE=0.5*(S1+S2) 
      A44=A22/TAVE 
C      A5=TAU(K2)/CMU 
      A55=(TAU(K)+TAU(K2))/(2.0*CMU) 
      A66=(TAU(K)-TAU(K2))/CMU 
      PL=A11/(EXP(A44)-1.00) 
      CALL EX(A66,A77) 
C      DINT=PL*A7*A6 
      BR=BR*A77+PL*(1.-A77) 
C      if(cmu.gt.0.99) WRITE(6,205) K,Tave,BR,A77,A66,PL,TAU(K),TAU(K2) 
 205  FORMAT(2X,I4,7(2X,E10.3)) 
 200  CONTINUE 
C     A8 IS THE BRIGHTNESS TEMPERATURE CORRESPONDING TO BR. 
      A99=A11/BR 
      A88=A22/(DLOG(A99+1.000)) 
C      WRITE(6,888) FREQ,CMU,BR,A8 
 888  FORMAT(2X,'FREQ,CMU,BR,T ARE ',/,4E15.5) 
      BRI=A88 
      RETURN 
      END 
C 
-----------------------------------------ATMOS---------------------------------
      TFLAT=1.0E10
      KOUNT=0

      TT=XT(1)
      PR=XP(1)
      DO 387 IJ=1,2
	KMAX=JMAX
      	JJ=1
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
     				+  tsoln(jj)=(4-kk)*soln(ii-1)/3.0
      			TCLH2O(JJ)=kk*(CLH2O(II)-clh2o(ii-1))/3.0 + clh2o(ii-1)
        		if(clh2o(ii).eq.0..and.clh2o(ii-1).gt.0.) 
     				+ tclh2o(jj)=(4-kk)*clh2o(ii-1)/3.0
       			TCLNH4SH(JJ)=kk*(CLNH4SH(II)-clnh4sh(ii-1))/3.0 + 
     				+ clnh4sh(ii-1)
        		if(clnh4sh(ii).eq.0..and.clnh4sh(ii-1).gt.0.) 
     				+ tclnh4sh(jj)=(4-kk)*clnh4sh(ii-1)/3.0
        		TCLNH3(JJ)=CLNH3(II)
        		TCLH2S(JJ)=CLH2S(II)
        		TCLCH4(JJ)=CLCH4(II)
        		TCLAR(JJ)=CLAR(II)
  101 		CONTINUE
      	KMAX=KMAX+2
  102 CONTINUE
      JMAX=KMAX
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
  103   CONTINUE
  387 CONTINUE

       KOUNT=0
C
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
C  COMPUTE DELZ, DIFFERENCE IN HEIGHT BETWEEN THE TWO 
C  PRESSURE LEVELS PO AND PR 
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
       KOUNT=KOUNT+1 
50CONTINUE 

 9999 NN=KOUNT 
      M5=0 
      RETURN 
      END 


---------------------------------END ATMOS-------------------------------

      SUBROUTINE MIE(X,A,B,QABS)
      COMPLEX ZM
      ZM=CMPLX(A,-B)
      ZM=4.0*X*((ZM*ZM-1.0)/(ZM*ZM+2.0))
      QABS=-AIMAG(ZM) 
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
      SUBROUTINE H2O(VNU,T,AH2O)                                                
C     SEE BERGE AND GULKIS PAPER IN GEHRELS JUPITER BOOK,PG 676,1969.           
C     SEE PAGE 118 OF GOODMAN'S THESIS.                                         
C  Updated July 2003 with de Boer's thesis (p. 66), and refs therein
        REAL*4 AON,T,P,F,V
	real*4 gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,PCO,PCO13,PHCN
	REAL*4 QLOG2,SHAPE2
	REAL*4 C,ZH2O,DV,VNU,gamma
        EXTERNAL QLOG2
        REAL*4 DF,DFC,DFD,EL(10),Q0(5),F0(10),FS,FV
        REAL*4 S0(10),S1(10),TCOR,X,X1,Y
	real*4 gh2(10),ghe(10),gh2o(10),xh2(10),xhe(10),xh2o(10)
C  S0 stands for A in de Boer's table 3.5
C        DATA Q0 /2.711, 2.5237, 2.2605, 1.8164,1.3778/ !LOG PF AT 300,225,150,75,37.5K  Values aren't right, I think
        DATA Q0 /2.2507, 2.0645, 1.804, 1.3649,0.9335/ !LOG PF AT 300,225,150,75,37.5K  These are right, I think
	 DATA F0(1),S0(1),EL(1),gh2(1),ghe(1),gh2o(1),xh2(1),xhe(1),xh2o(1)
     &  /22235.15, 1.0, 644.0,  2.395, 0.67, 10.67, 0.90, 0.515, 0.626/
        DATA F0(2),S0(2),EL(2),gh2(2),ghe(2),gh2o(2),xh2(2),xhe(2),xh2o(2)
     &  /183310.12, 41.9, 196.0, 2.400, 0.71, 11.64, 0.95, 0.490, 0.649/
        DATA F0(3),S0(3),EL(3),gh2(3),ghe(3),gh2o(3),xh2(3),xhe(3),xh2o(3)
     &  /323000.0, 334.4, 1850.0, 2.395, 0.67, 9.59, 0.90, 0.515, 0.420 /
        DATA F0(4),S0(4),EL(4),gh2(4),ghe(4),gh2o(4),xh2(4),xhe(4),xh2o(4)  
     &  /325153.8, 115.7, 454.0,  2.395, 0.67, 11.99, 0.90, 0.490, 0.619/
        DATA F0(5),S0(5),EL(5),gh2(5),ghe(5),gh2o(5),xh2(5),xhe(5),xh2o(5)
     &  /380196.8, 651.8, 306.0,  2.390, 0.63, 12.42, 0.85, 0.540, 0.630/
        DATA F0(6),S0(6),EL(6),gh2(6),ghe(6),gh2o(6),xh2(6),xhe(6),xh2o(6)
     &  /390000.0, 127.0, 2199.0,  2.395, 0.67, 9.16, 0.90, 0.515, 0.330/
        DATA F0(7),S0(7),EL(7),gh2(7),ghe(7),gh2o(7),xh2(7),xhe(7),xh2o(7)
     &  /436000.0, 191.4, 1507.0, 2.395, 0.67, 6.32, 0.90, 0.515, 0.290/
        DATA F0(8),S0(8),EL(8),gh2(8),ghe(8),gh2o(8),xh2(8),xhe(8),xh2o(8)
     &  /438000.0, 697.6, 1070.0, 2.395, 0.67, 8.34, 0.90, 0.515, 0.360 /
        DATA F0(9),S0(9),EL(9),gh2(9),ghe(9),gh2o(9),xh2(9),xhe(9),xh2o(9)
     &  /442000.0, 590.2, 1507.0, 2.395, 0.67, 6.52, 0.90, 0.515, 0.332 /
        DATA F0(10),S0(10),EL(10),gh2(10),ghe(10),gh2o(10),xh2(10),xhe(10),xh2o(10)
     &  /448000.8, 973.1, 412.0,  2.395, 0.67, 11.57, 0.90, 0.515, 0.510/
       COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
       COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11 
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                                
C                                                                               
C     CONVERT VNU(CM-1) TO GHZ. C = 29.97   
	NLINES=10                                             
	C=2.99791E1
	F=VNU*C*1.0E3                  !FREQU IN MHZ
	V=VNU*C				!FREQ IN GHZ
	rho=1.0E12*18.0*PH2O/(8.34143E7*T)	!rho in g/m^3
	PA=0.81*PH2 + 0.35*PHe
	AON=0.0
        DO N = 1,NLINES
	flc=f0(n)                       !line center in MHz
	VLINE=F0(N)/1.0E3		!LINE CENTER IN GHZ
        TCOR = S0(N) * exp(-EL(N)/T)
C  THE GROSS LINE SHAPE
        X2=(-v*v+vline*vline)
        X2=X2*X2
	g1=gh2(N)*PH2*(300.0/T)**xh2(N)
	g2=ghe(N)*PHe*(300.0/T)**xhe(N)
	g3=gh2o(N)*PH2O*(300.0/T)**xh2o(N)
	gamma=g1+g2+g3
        X3=gamma*gamma*4.*v*v
        FJK=TCOR*gamma/(X2+X3)
C  Divide FJK by 10 to test Steffes student's idea that DeBoer is off by 10.
	FJK=FJK/10.0
	AON=AON+FJK
C  perhaps multiply by 4/pi  
	enddo
C	AON=AON*4.0/3.14159	!we think that it is all included in deBoer's stuff; It is
	X1=(300.0/T)**2.1
	X1=1.08E-11*X1*PA*v*v*rho
	X2=(300.0/T)**2.5
	X2=2.0*v*v*rho*X2*2.3E-5*AON
C	write(6,*) AON,X1,X2
	AON=X1+X2
C  
C     H2O ABSORPTION COEFFICIENT IN CM-1 ;
       AH2O=AON                                                          
C                                                                               
      RETURN                                                                    
      END
      SUBROUTINE H2Oold(VNU,T,AH2O)
C This is the old H2O routine                                                
C     SEE BERGE AND GULKIS PAPER IN GEHRELS JUPITER BOOK,PG 676,1969.           
C     SEE PAGE 118 OF GOODMAN'S THESIS.                                         
       COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
       COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11  
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                    
C                                                                               
C     CONVERT VNU(CM-1) TO GHZ. C = 29.97                              
      V=VNU*C                                                            
C                                                                               
C     EVALUATE EQUATION 29.                                                     
      A2=V/29.97                                                                
      A3=(A2-0.74)*(A2-0.74)                                                    
      A4=(A2+0.74)*(A2+0.74)                                                    
      A5=9.88E-2*T5*((0.81*PH2)+(0.35*PHE))                                     
      A5=A5/1.013                                                               
      A6=A5*A5                                                                  
      A7=(A5/(A3+A6))+(A5/(A4+A6))                                              
C     H2O ABSORPTION COEFFICIENT IN CM-1                                        
      AH2O=PH2O*T4*V*V*((9.07E-9*A7)+(1.45E-7*A5))                              
C  CORRECTION TERM B(T)/B(293)   SEE GOODMAN,P.117                              
      AH2O=AH2O*1.179                                                           
      AH2O=AH2O/1.013                                                           
C                                                                               
      IF (NOPR .EQ. 1) WRITE(6,10) V,A2,A3,A4,A5,A6,A7,AH2O   
  10  FORMAT(2X,'H2O V(GHZ),A2,A3,A4,A5,A6,A7,AH2O',/,                          
     1  2X,1P,8(2X,E10.3),/)                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE NH3TOM(VNU,T,ANH3)             
C     CALCULATE THE AMMONIA ABSORPTION COEFFICIENT (CM-1)                       
C     USED BY GULKIS AND POYNTER (1972),AND GIVEN IN BERGE                      
C     AND GULKIS, GEHRELS JUPITER BOOK,PAGE 675.                                
      DIMENSION VKAK(16,16),DVKAK(16,16)                                        
       COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
      COMMON/U2/C1,C2,C3,C4                                                     
      COMMON/U3/VKAK,DVKAK                                                      
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11 
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                        
C
C   VKAK(J,K)=CENTER FREQU FOR THE (J,K) TRANSITION ACCORDING TO THE            
C  TABULATED VALUES OF POYNTER AND GULKIS (1975).                               
C  DVKAK(J,K) ARE THE CORRESPONDING SELF BROADENED LINE WIDTHS                  
C  AS TABULATED BY THEM (RESP. GHZ, MHZ/TORR)                                   
C     CHANGE VNU(CM-1) TO GHZ.                                                  
      V=(C*VNU)                                                          
C  CHANGE PRESSURES FROM BARS TO ATM, AND BACK AT END OF SUBROUTINE; MARCH89
	PH2=PH2/1.013
	PHE=PHE/1.013
	PH2O=PH2O/1.013
	PNH3=PNH3/1.013
	PCH4=PCH4/1.013
	TPPR=PH2+PHE+PNH3+PH2O+PCH4
C     CORRECTION FRACTOR IS GIVEN ON PAGE 676 OF THE BERGE AND GULKIS           
C     PAPER, GEHRELS JUPITER BOOK.                                              
      CORBG=1.0075+((0.0308+(0.0552*(PH2/T)))*(PH2/T))               
	CORR=-0.33664+T/110.4-T*T/70600.
	IF(T.GT.320.) CORR=1.1115
	IF(T.LT.180.) then
	x=(180.*180./70600.) - 0.337
	y=(1./110.4) - (180./35300.)
	corr=x + y*T
	endif
C
	CORR=CORR*CORBG
	IF(TPPR.LE.1.0) CORR=1.0
C                                                                               
      SUM=0.0                                                                   
C     CALCULATE THE INVERSION LINE WAVENUMBERS                                  
C     FOR NH3, FROM THE DATA IN POYNTER AND KAKAR                               
C     AP J SUPPL SERIES,VOL 29,P87,1975.                                        
C     J AND K ARE QUANTUM NUMBERS.                                              
C                                                                               
C     NLINE ARE THE NUMBER OF LINES TO BE CALCULATED.                           
      NLINE=16                                                                  
C                                                                               
      DO 10 J=1,NLINE                                                           
      L1=J*(J+1)                                                                
      L3=0                                                                      
      DO 20 K=1,J                                                               
      L2=K*K                                                                    
      L3=L3+1                                                                   
C     SB IS THE SELF BRODENED LINE WIDTH OF NH3 IN MHZ/TORR.                    
C     C4 IS GIVEN IN THE BLOCK DATA STATEMENT.                                  
C     SK IS 1.5 OR 3, DEPENDING UPON THE K VALUE.                               
      SK=1.5                                                                    
      IF (L3 .NE. 3) GO TO 21                                                   
      SK=3.0                                                                    
      L3=0                                                                      
  21  VLINE=VKAK(J,K)                                                           
       SB=DVKAK(J,K)                                                            
      IF (NOPR .EQ. 1) WRITE(6,22) J,K,SK,VLINE,SB              
  22  FORMAT(2X,'J,K,SK,VLINE(GHZ),SB(MHZ/TORR)',/,                             
     1 2X,I3,2X,I3,2X,F4.1,2X,1P,2(2X,E10.3))                                   
C     CALCULATE A(J,K).                                                         
      CALL SAJTOM(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)                                     
      IF(J.EQ.1.AND.K.EQ.1) ZZV=2.0*AJK*EXP(23.3/T)/3.0                         
C     FJK IS THE FREQUENCY DEPENDENT LINE SHAPE OF BEN REUVEN.                  
      CALL SFJTOM(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      SB=AJK*FJK                                                                
      SUM=SUM+(CORR*AJK*FJK)                                                         
c  no correction term in Joiner/Steffes;
c     SUM=SUM+(AJK*FJK)
  20  CONTINUE                                                                  
  10  CONTINUE                                                                  
C                                                                               
C  CALCULATE ROTATION LINE J=0;K=0                                              
C  VLINE=5.73E5 MHZ; CHOSE DVLINE=14. MHZ                                       
C                                                                               
      L1=2                                                                      
      L2=1                                                                      
      VLINE=5.73E2                                                              
      SB=14.0                                                                   
      GM=(2.318*T3*PH2)+(0.79*T3*PHE)+(0.75*T2*SB*PNH3)                         
      GM=GM                                                               
      CALL SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      ZZV1=ZZV*7.0E-2*T*(1.0-EXP(-28.6/T))                                      
      ZZV2=ZZV*5.25E-2*T*EXP(-23.3/T)*(1.0-EXP(-57.2/T))                        
      ZZV3=ZZV*7.50E-2*T*EXP(-28.6/T)*(1.0-EXP(-57.2/T))                        
      SUM=SUM+(ZZV1*FJK)                                                        
      L1=2                                                                      
      L2=1                                                                      
      VLINE=1169.0                                                              
      SB=23.8                                                                   
      CALL SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      SUM=SUM+(ZZV2*FJK)+(ZZV3*FJK)                                             
C                                                                               
C                                                                               
C     THE AMMONIA ABSORPTION COEFFICIENT IS IN CM-1.                            
        ANH3=sum
      IF (NOPR .EQ. 1) WRITE(6,30) CORR,SUM,ANH3                                
  30  FORMAT(2X,'CORR,SUM,ANH3',2X,1P,3(2X,E10.3),/)                            
C                                                                               
	PH2=PH2*1.013
	PHE=PHE*1.013
	PH2O=PH2O*1.013
	PNH3=PNH3*1.013
	PCH4=PCH4*1.013
      RETURN                                                                    
      END                                                                       
C                                                            
      SUBROUTINE NH3(VNU,T,ANH3)                                                
C     CALCULATE THE AMMONIA ABSORPTION COEFFICIENT (CM-1)                       
C     USED BY GULKIS AND POYNTER (1972),AND GIVEN IN BERGE                      
C     AND GULKIS, GEHRELS JUPITER BOOK,PAGE 675.                                
      DIMENSION VKAK(16,16),DVKAK(16,16)                                        
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
      COMMON/U2/C1,C2,C3,C4                                                     
      COMMON/U3/VKAK,DVKAK                                                      
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4           
C
C   VKAK(J,K)=CENTER FREQU FOR THE (J,K) TRANSITION ACCORDING TO THE            
C  TABULATED VALUES OF POYNTER AND GULKIS (1975).                               
C  DVKAK(J,K) ARE THE CORRESPONDING SELF BROADENED LINE WIDTHS                  
C  AS TABULATED BY THEM (RESP. GHZ, MHZ/TORR)                                   
C     CHANGE VNU(CM-1) TO GHZ.                                                  
      V=(C*VNU)                                                          
C  CHANGE PRESSURES FROM BARS TO ATM, AND BACK AT END OF SUBROUTINE; MARCH89
	PH2=PH2/1.013
	PHE=PHE/1.013
	PH2O=PH2O/1.013
	PNH3=PNH3/1.013
	PCH4=PCH4/1.013
C     CORRECTION FRACTOR IS GIVEN ON PAGE 676 OF THE BERGE AND GULKIS           
C     PAPER, GEHRELS JUPITER BOOK.                                              
C                                                                               
      CORR=1.0075+((0.0308+(0.0552*(PH2/T)))*(PH2/T))               
C      CORR=SQRT(PH2/T)
C      CORR=1.0+(0.155*CORR+0.52*SQRT(CORR))*SQRT(CORR)
      IF (NOPR .EQ. 1) WRITE(6,5) V,CORR                                        
   5  FORMAT(2X,'NH3 V(GHZ),CORR',2X,1P,2(2X,E10.3))                            
      IF (NOPR .EQ. 1) WRITE(6,6) C1,C2,C3,C4                                   
   6  FORMAT(2X,'CONSTANTS C1,C2,C3,C4',2X,1P,4(2X,E10.3))                      
C                                                                               
      SUM=0.0                                                                   
C     CALCULATE THE INVERSION LINE WAVENUMBERS                                  
C     FOR NH3, FROM THE DATA IN POYNTER AND KAKAR                               
C     AP J SUPPL SERIES,VOL 29,P87,1975.                                        
C     J AND K ARE QUANTUM NUMBERS.                                              
C                                                                               
C     NLINE ARE THE NUMBER OF LINES TO BE CALCULATED.                           
      NLINE=16                                                                  
C                                                                               
      DO 10 J=1,NLINE                                                           
      L1=J*(J+1)                                                                
      L3=0                                                                      
      DO 20 K=1,J                                                               
      L2=K*K                                                                    
      L3=L3+1                                                                   
C     SB IS THE SELF BRODENED LINE WIDTH OF NH3 IN MHZ/TORR.                    
C     C4 IS GIVEN IN THE BLOCK DATA STATEMENT.                                  
C     SK IS 1.5 OR 3, DEPENDING UPON THE K VALUE.                               
      SK=1.5                                                                    
      IF (L3 .NE. 3) GO TO 21                                                   
      SK=3.0                                                                    
      L3=0                                                                      
  21  VLINE=VKAK(J,K)                                                           
       SB=DVKAK(J,K)                                                            
      IF (NOPR .EQ. 1) WRITE(6,22) J,K,SK,VLINE,SB                              
  22  FORMAT(2X,'J,K,SK,VLINE(GHZ),SB(MHZ/TORR)',/,                             
     1 2X,I3,2X,I3,2X,F4.1,2X,1P,2(2X,E10.3))                                   
C     CALCULATE A(J,K).                                                         
      CALL SAJK(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)                                     
      IF(J.EQ.1.AND.K.EQ.1) ZZV=2.0*AJK*EXP(23.3/T)/3.0                         
C     FJK IS THE FREQUENCY DEPENDENT LINE SHAPE OF BEN REUVEN.                  
      CALL SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      SB=AJK*FJK                                                                
      SUM=SUM+(CORR*AJK*FJK)                                                         
      IF (NOPR .EQ. 1) WRITE(6,23) AJK,FJK,SB                                   
  23  FORMAT(2X,'AJK,FJK,(AJK*FJK)',2X,1P,3(2X,E10.3))                          
  20  CONTINUE                                                                  
  10  CONTINUE                                                                  
C                                                                               
C  CALCULATE ROTATION LINE J=0;K=0                                              
C  VLINE=5.73E5 MHZ; CHOSE DVLINE=14. MHZ                                       
C                                                                               
      L1=2                                                                      
      L2=1                                                                      
      VLINE=5.73E2                                                              
      SB=14.0                                                                   
      GM=(2.318*T3*PH2)+(0.79*T3*PHE)+(0.75*T2*SB*PNH3)                         
      GM=GM                                                               
      CALL SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      ZZV1=ZZV*7.0E-2*T*(1.0-EXP(-28.6/T))                                      
      ZZV2=ZZV*5.25E-2*T*EXP(-23.3/T)*(1.0-EXP(-57.2/T))                        
      ZZV3=ZZV*7.50E-2*T*EXP(-28.6/T)*(1.0-EXP(-57.2/T))                        
      SUM=SUM+(ZZV1*FJK)                                                        
      L1=2                                                                      
      L2=1                                                                      
      VLINE=1169.0                                                              
      SB=23.8                                                                   
      CALL SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      SUM=SUM+(ZZV2*FJK)+(ZZV3*FJK)                                             
C                                                                               
C                                                                               
C     THE AMMONIA ABSORPTION COEFFICIENT IS IN CM-1.                            
      ANH3=SUM                                                             
      IF (NOPR .EQ. 1) WRITE(6,30) CORR,SUM,ANH3                                
  30  FORMAT(2X,'CORR,SUM,ANH3',2X,1P,3(2X,E10.3),/)                            
C                                                                               
	PH2=PH2*1.013
	PHE=PHE*1.013
	PH2O=PH2O*1.013
	PNH3=PNH3*1.013
	PCH4=PCH4*1.013
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE SAJK(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)  
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
C     GM IS GAMMA(J,K) IN GHZ.                                                  
C   2.318 IPV 3.318
      GM=(2.318*T3*PH2)+(0.79*T3*PHE)+(0.75*T2*SB*PNH3)                         
C      GM=(4.8*PH2*(T2**0.50))+(0.79*T3*PHE)+(0.75*T2*SB*PNH3)
C      GM=17.4*T2*PNH3+1.73*T8*PH2+0.40*T8*PHE
      A1=(2.*J+1.)*FLOAT(L2)/FLOAT(L1)                                          
      A5=((2.98*L1)-(1.09*L2))*T7                                               
C     AJK IS A(J,K),SEE EQN 24 OF BERGE AND GULKIS PAPER,GEHRELS                
C     JUPITER BOOK,PAGE 675.                                                    
      CALL EX(A5,A10)                                                           
      AJK=1.23E3*A1*SK*PNH3*T6*A10                                              
      IF (NOPR .EQ. 1) WRITE(6,10) GM,A1,A5,AJK                                 
  10  FORMAT(2X,'GM,A1,A5,AJK',2X,1P,4(2X,E10.3))                               
      RETURN                                                                    
      END                                                                       
      SUBROUTINE SAJTOM(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)                               
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4 
C     GM IS GAMMA(J,K) IN GHZ.                                                  
	PTOT=PH2+PHE+PNH3+PH2O+PCH4
c  spilker uses equations below
	GNH3=0.74
	GHE=0.46+T/3000.
	R=8.79*EXP(-T/83.)
	X=2.122*EXP(-T/116.8)
	Y=EXP(9.024-T/20.3)
	Y=(Y-0.9918+PH2)**R
	GH2=2.34*(1.0-X/Y)
c  Joanna Joiner uses:
c        gh2=1.69
c        ghe=0.75
c        gnh3=0.60
      GM=(GH2*T3*PH2)+(GHE*T3*PHE)+(GNH3*T2*SB*PNH3)                         
      A1=(2.*J+1.)*FLOAT(L2)/FLOAT(L1)                                          
      A5=((2.98*L1)-(1.09*L2))*T7                                               
C     AJK IS A(J,K),SEE EQN 24 OF BERGE AND GULKIS PAPER,GEHRELS                
C     JUPITER BOOK,PAGE 675.                                                    
      CALL EX(A5,A10)                                                           
      AJK=1.23E3*A1*SK*PNH3*T6*A10                                              
      IF (NOPR .EQ. 1) WRITE(6,10) GM,A1,A5,AJK                                 
  10  FORMAT(2X,'GM,A1,A5,AJK',2X,1P,4(2X,E10.3))                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C 
      SUBROUTINE SFJK(T,L1,L2,V,VLINE,GM,GH2,SB,FJK) 
	real*4 t3,gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,
     + deltot,PCO,PCO13,PHCN
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
       COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
C     FJK IS THE BEN REUVEN LINE SHAPE,SEE EQN 25 OF BERGE AND                  
C     GULKIS ARTICLE,PAGE 675 OF GEHRELS JUPITER BOOK.                          
C                                                                               
      X=PH2+PHE+PH2O+PNH3+PCH4                                                  
C  IF X.LT.1 ATM VAN VLECK-WEISKOPF LINE SHAPE CAN BE                           
C  USED; HOWEVER, BEN REUVEN IS STILL BETTER OUT IN THE                         
C  WINGS OF THE AMMONIA LINES. 
c     A2 IS THE COUPLING FACTOR ZETA(J,K) IN GHZ.                               
C	GOTO 500
C  for Gross lineshape, goto 500 (first VVW, then gross
      IF(VLINE.GT.100.0) GOTO 500                                               
      A2=(1.92*PH2*T3)+(0.49*T2*PNH3*SB)+                                       
     1  (0.3*T3*PHE)                                                            
C      A2=(0.49*T2*PNH3*SB)+(0.3*T3*PHE)+
C     +(4.80*(T2**0.5)-0.6*(T2**2.0))*PH2
C      A2=(1.73*T8-0.87*T10)*PH2+
C     +  (0.40*T8-0.21*T9)*PHE+
C     +  (17.4*T2-6.0*T11)*PNH3
C      A2=0.655*GM
C     A4 IS THE PRESSURE SHIFT DEL (GHZ)                                        
	A4=-0.45*PNH3
C     SQUARES OF VARIOUS TERMS.                                                 
      A5=(VLINE+A4)*(VLINE+A4)                                                  
C      A5=VLINE*VLINE
      A6=GM*GM                                                                  
      A7=A2*A2                                                                  
      A8=V*V                                                                    
C     NUMERATOR                                                                 
      A9=((GM-A2)*A8)+((GM+A2)*(A5+A6-A7))                                      
      A10=A8-A5-A6+A7                                                           
C     DENOMINATOR                                                               
      A10=(A10*A10)+(4.0*A8*A6)                                                 
C     THE BEN REUVEN LINE SHAPE.                                                
      FJK=2.0*A8*(A9/A10)                                                       
      IF (NOPR .EQ. 1) WRITE(6,10) A2,A4,A5,A6,A7,A8,A9,                        
     1  A10,FJK                                                                 
      FFF=FJK                                                                   
      GOTO 501                                                                  
  500 CONTINUE                                                                  
C  THE VAN VLECK-WEISKOPF LINE SHAPE                                            
c
      Z1=0.755*PH2+0.231*PHE                                                    
      Z2=FLOAT(L2)/FLOAT(L1)                                                    
      Z1=Z1*(Z2**0.3333)                                                  
      Z2=SQRT(Z2)*PNH3*6.23
      Z=(Z1+Z2)*1.0E3/T                                                         
      X1=V*V*Z                                                                  
      X2=(V-VLINE)*(V-VLINE)                                                    
      X3=(V+VLINE)*(V+VLINE)                                                    
      FJK=X1/(X2+Z*Z)+X1/(X3+Z*Z)                                               
C      GOTO 501                                                                
  502 CONTINUE
C  THE GROSS LINE SHAPE
        X2=(-V*V+VLINE*VLINE)
        X2=X2*X2
       if(vline.lt.1000.) Z1=123.3551*PH2+60.7441*PHE  
	if(vline.ge.1000.) z1=2.885*PH2+1.00*PHE
C      Z2=FLOAT(L2)/FLOAT(L1)                                                    
	Z2=300./T
      Z1=Z1*(Z2**0.3)
	if(vline.lt.1000.) Z=154.1939*PNH3*Z2
	if(vline.ge.1000.) Z=22.989*PNH3*Z2
       Z=Z1+Z
        X3=Z*Z*4.*V*V
        X1=v*VLINE*Z*4.
        X1=X1*V*VLINE
        FJK=X1/(X2+X3)
       FFF=FJK                                                                   
  501 CONTINUE                                                                  
  10  FORMAT(2X,'A2,A4,A5,A6,A7,A8,A9,A10,FJK',/,                               
     1  2X,1P,9(2X,E10.3))                                                      
      RETURN                                                                    
      END                                                                       
C                                                            
      SUBROUTINE SFJTOM(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)
	real*4 t3,gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,
     + deltot,PCO,PCO13,PHCN
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
       COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
C                                                                               
C     FJK IS THE BEN REUVEN LINE SHAPE,SEE EQN 25 OF BERGE AND                  
C     GULKIS ARTICLE,PAGE 675 OF GEHRELS JUPITER BOOK.                          
C                                                                               
      X=PH2+PHE+PH2O+PNH3+PCH4  
C  IF X.LT.1 ATM VAN VLECK-WEISKOPF LINE SHAPE SHOUL BE                           
C  USED; HOWEVER, BEN REUVEN IS STILL BETTER OUT IN THE                         
C  WINGS OF THE AMMONIA LINES. SO OPTION VAN VLECK-WEISKOPF                     
C  IS NEGLECTED BY SETTING IN NEXT LINE:                                        
C  IF(X.LT.0.0) INSTEAD OF IF(X.LT.1.0)                                         
C
C       IF(X.LT.1.) GOTO 500                                                     
C
C     A2 IS THE COUPLING FACTOR ZETA(J,K) IN GHZ.                               
c  spliker's formalism uses equations below
	ZNH3=0.5
	ZHE=0.28-T/1750.
	ZH2=5.7465+GH2*(-7.7644+GH2*(9.1931+GH2*(-5.6816+1.2307*GH2)))
	if(v.le.30.0) goto 600
c  Joanna Joiner's paper gives parameters below
C  These parameters are better at freq. (v in GHz) above 40 GHz
C  see de Boer's thesis, p. 61-64
        znh3=0.20
       zh2=1.35
        zhe=0.30
C 
      IF(VLINE.GT.100.0) GOTO 500                                               
 600   A2=(ZH2*PH2*T3)+(ZNH3*T2*PNH3*SB)+(ZHE*T3*PHE)                                                            
C     A4 IS THE PRESSURE SHIFT DEL (GHZ)                                        
	A4=-0.45*PNH3
C     SQUARES OF VARIOUS TERMS.                                                 
      A5=(VLINE+A4)*(VLINE+A4)                                                  
      A6=GM*GM                                                                  
      A7=A2*A2                                                                  
      A8=V*V                                                                    
C     NUMERATOR                                                                 
      A9=((GM-A2)*A8)+((GM+A2)*(A5+A6-A7))                                      
      A10=A8-A5-A6+A7                                                           
C     DENOMINATOR                                                               
      A10=(A10*A10)+(4.0*A8*A6)                                                 
C     THE BEN REUVEN LINE SHAPE.                                                
      FJK=2.0*A8*(A9/A10)                                                       
      FFF=FJK                                                                   
C      IF(V.LT.29.0) GOTO 501 
c  April 2003: we skip VVW because Spilker's formalism does already 
c  do the right thing, I believe                                                   
	GOTO 501
  500 CONTINUE                                                                  
C  THE VAN VLECK-WEISKOPF LINE SHAPE                                            
      Z1=0.755*PH2+0.231*PHE                                                    
      Z2=FLOAT(L2)/FLOAT(L1)                                                    
      Z1=Z1*(Z2**0.3333)                                                  
      Z2=SQRT(Z2)*PNH3*6.23
      Z=(Z1+Z2)*1.0E3/T                                                         
      X1=V*V*Z                                                                  
      X2=(V-VLINE)*(V-VLINE)                                                    
      X3=(V+VLINE)*(V+VLINE)                                                    
      FJK=X1/(X2+Z*Z)+X1/(X3+Z*Z)                                               
C      GOTO 501 
  502 CONTINUE
C  THE GROSS LINE SHAPE
        X2=(-V*V+VLINE*VLINE)
        X2=X2*X2
       if(vline.lt.1000.) Z1=123.3551*PH2+60.7441*PHE  
	if(vline.ge.1000.) z1=2.885*PH2+1.00*PHE
C      Z2=FLOAT(L2)/FLOAT(L1)                                                    
	Z2=300./T
      Z1=Z1*(Z2**0.3)
	if(vline.lt.1000.) Z=154.1939*PNH3*Z2
	if(vline.ge.1000.) Z=22.989*PNH3*Z2
       Z=Z1+Z
        X3=Z*Z*4.*V*V
        X1=v*VLINE*Z*4.
        X1=X1*V*VLINE
        FJK=X1/(X2+X3)
       FFF=FJK                           
  501 CONTINUE                                                                  
  10  FORMAT(2X,'A2,A4,A5,A6,A7,A8,A9,A10,FJK',/,                               
     1  2X,1P,9(2X,E10.3))                                                      
      RETURN                                                                    
      END                                                                       
C
C                                           
      SUBROUTINE H2H2(VNU,T,A,B)                                                
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V3/QROT,QR1,QR2,S1,S2,S3,S4,P1,P2,P3,P4                            
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
      COMMON/V6/P5,P6,P7,P8                                                     
      COMMON/V7/PRO1(20),PRO2(20),STOR(20)                                      
      COMMON/W1/FP,REB,RHO(20)                                                  
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
C     H2 OPACITIES ARE CALCULATED FOLLOWING THE PAPER "ANALYSIS OF              
C     THE SHAPE OF THE FAR-INFRARED SPECTRA OF H2-H2 AND H2-HE                  
C     COLLISIONS,BY COHEN AND BIRNBAUM,APRIL 1981.                              
C     A RELATED PAPER IS GIVEN BY BIRNBAUM AND COHEN, CANADIAN                  
C     JOURNAL OF PHYSICS,VOL 54,PAGE 593,1976.                                  
C     THIS WORK IS NBSIR 80-2175(R).                                            
C     A AND B COEFFICIENTS IN THE TRAFTON TRADITION ARE CALCULATED.             
C     SEE TRAFTON, ASTROPHYSICAL JOURNAL,V 147,PAGE 765,1967.                   
C                                                                               
C     THE H2-H2 COLLISIONS ARE CONSIDERED HERE FOR THE A VALUES.                
C                                                                               
      TEMP=T                                                                    
C     BOLTZMANN STATISTICS ARE CALCULATED.                                      
      CALL JROT(TEMP)                                                           
C     EVALUATE SQ,TAU1,TAU2 PARAMETERS.                                         
      CALL PARAM(1)                                                             
      IF (NOPR .EQ. 1) WRITE(6,330) SQ,TAU1,TAU2                                
 330  FORMAT(/,2X,'FOR H2-H2 COLLISIONS',/,2X,'SQ,TAU1,TAU2',                   
     1  2X,1P,3(2X,E10.3))                                                      
C     EVALUATE THE A FACTOR.                                                    
      CALL FAC(10,AFAC)                                                         
      IF (NOPR .EQ. 1) WRITE(6,355) AFAC                                        
 355  FORMAT(2X,'AFAC=',2X,1P,E10.3)                                            
C                                                                               
      V=VNU                                                                     
      IF (NOPR .EQ. 1) WRITE(6,310) V                                           
 310  FORMAT(2X,'V(CM-1)= ',1P,E10.3)                                           
C     OM IS TWO PI TIMES C TIMES THE WAVENUMBER.                                
C     FREQ IS SEC-1.                                                            
      OM=TPC*V                                                                  
      FREQ=C*V*1.0E9                
      BETA=1.0/(BOLC*TEMP)                                                      
      IF (NOPR .EQ. 1) WRITE(6,315) OM,FREQ,BETA                                
 315  FORMAT(2X,'OM(SEC-1),FREQ(SEC-1),BETA(ERG-1)',/,                          
     1  2X,1P,3(2X,E10.3))                                                      
C     EVALUATE THE 1-EXP TERM.                                                  
      TX=BETA*HB*OM                                                             
      CALL EX(TX,S)                                                             
C  AW2 SHOULD BE MULTIPLIED BY 1.0E-36
      AW2=AW*OM*(1.0-S)*SQ                                                      
      IF (NOPR .EQ. 1) WRITE(6,340) AW2                                         
 340  FORMAT(2X,'AW2= ',1P,E10.3)                                               
C                                                                               
C     THE SIX TERMS OF EQUATION 7 OF NBS 80 ARE CALCULATED.                     
C     THE FIRST TERM IS EVALUATED.                                              
      OM2=OM                                                                    
      CALL GAMMA(OM2,GAM)                                                       
      IF (NOPR .EQ. 1) WRITE(6,356) TEMP,TAU1,TAU2,OM,BETA,GAM                  
 356  FORMAT(2X,1P,6(2X,E10.3))                                                 
      TERM1=AFAC*GAM                                                            
C         ZERO AND TWO ARE EVALUATED.                                           
      CALL EVAL(0,2,OM,B1,B2)                                                   
      TERM2=(RHO(1)*B1)+(RHO(3)*B2)                                             
C         ONE AND THREE ARE EVALUATED.                                          
      CALL EVAL(1,3,OM,B1,B2)                                                   
      TERM3=(RHO(2)*B1)+(RHO(4)*B2)                                             
      TERM3=(9.0/5.0)*TERM3                                                     
C         TWO AND FOUR ARE EVALUATED.                                           
      CALL EVAL(2,4,OM,B1,B2)                                                   
      TERM4=(RHO(3)*B1)+(RHO(5)*B2)                                             
      TERM4=(18.0/7.0)*TERM4                                                    
C         THREE AND FIVE ARE EVALUATED.                                         
      CALL EVAL(3,5,OM,B1,B2)                                                   
      TERM5=(RHO(4)*B1)+(RHO(6)*B2)                                             
      TERM5=(10.0/3.0)*TERM5                                                    
C                                                                               
C  SUM3 SHOULD BE MULTIPLIED BY 1.0E-14
      SUM3=TERM1+TERM2+TERM3+TERM4+TERM5                                        
      IF (NOPR .EQ. 1) WRITE(6,360) TERM1,TERM2,TERM3,TERM4,                    
     1  TERM5,SUM3                                                              
 360  FORMAT(2X,'TERM1,TERM2,TERM3,TERM4,TERM5,SUM',/,                          
     1  2X,1P,6(2X,E10.3))                                                      
C     NOTE:THE RESULT IS MULTIPLIED BY 1.0E36                                   
      SUM=AW2*SUM3*2.9979E-4
C                                                                               
C       NOTICE FACTOR OF 2.  NEED TO RESOLVE THIS DISCREPANCY                   
C       WITH BIRNBAUM. (FACTOR OF 2 GIVES D(V) VALUES WHICH                     
C       AGREE WITH BIRNBAUM'S GRAPHED VALUES.)                                  
C                                                                               
C     A(T) IN THE TRAFTON STYLE IS PUT INTO A                                   
      A=SUM                                                                     
      IF (NOPR .EQ. 1) WRITE(6,104) V,T,A                                       
 104  FORMAT(2X,'V,T,A = ',2X,1P,3(2X,E10.3))                                   
C                                                                               
C     THE H2-HE COLLISIONS ARE CONSIDERED HERE FOR THE B VALUES.                
C                                                                               
      TEMP=T                                                                    
C     BOLTZMANN STATISTICS ARE CALCULATED.                                      
      CALL JROT(TEMP)                                                           
      IF (NOPR .EQ. 1) WRITE(6,200)                                             
 200  FORMAT(/,2X,'FOR H2-HE COLLISIONS')                                       
C     EVALUATE SQ,TAU1,TAU2 PARAMETERS.                                         
C     FOR THE ISOTROPIC S AND TAU FACTORS.                                      
      CALL PARAM(2)                                                             
      SQA=SQ                                                                    
      TAU1A=TAU1                                                                
      TAU2A=TAU2                                                                
      IF (NOPR .EQ. 1) WRITE(6,210) SQA,TAU1A,TAU2A                             
 210  FORMAT(2X,'ISOTROPIC SQA,TAU1,TAU2A',2X,1P,3(2X,E10.3))                   
C     EVALUATE SQ,TAU1,TAU2 PARAMETERS.                                         
C     FOR THE ANISOTROPIC S AND TAU FACTORS.                                    
      CALL PARAM(3)                                                             
      SQB=SQ                                                                    
      TAU1B=TAU1                                                                
      TAU2B=TAU2                                                                
      IF (NOPR .EQ. 1) WRITE(6,220) SQB,TAU1B,TAU2B                             
 220  FORMAT(2X,'ATISOTROPIC SQB,TAU1B,TAU2B',2X,1P,3(2X,E10.3))                
C     EVALUATE THE A FACTOR.                                                    
      CALL FAC(10,AFAC)                                                         
      IF (NOPR .EQ. 1) WRITE(6,355) AFAC                                        
C                                                                               
      V=VNU                                                                     
      IF (NOPR .EQ. 1) WRITE(6,310) V                                           
C     OM IS TWO PI TIMES C TIMES THE WAVENUMBER.                                
C     FREQ IS SEC-1.                                                            
      OM=TPC*V                                                                  
      FREQ=C*V*1.0E9                                                                  
      BETA=1.0/(BOLC*TEMP)                                                      
      IF (NOPR .EQ. 1) WRITE(6,315) OM,FREQ,BETA                                
C     EVALUATE THE 1-EXP TERM.                                                  
      TX=BETA*HB*OM                                                             
      CALL EX(TX,S)                                                             
C     NOTE THAT THE SQ FACTOR IS NOT IN THE AW2 EXPRESSION.                     
C     FOR H2-HE, THERE IS AN ADDITIONAL FACTOR OF 2.                            
C     THE SQA AND SQB FACTORS WILL BE INTRODUCED LATER.                         
      AW2=2.0*AW*OM*(1.0-S)                                                     
      IF (NOPR .EQ. 1) WRITE(6,340) AW2                                         
C                                                                               
C     THE ISOTROPIC TERM IS CALCULATED HERE.                                    
      SQ=SQA                                                                    
      TAU1=TAU1A                                                                
      TAU2=TAU2A                                                                
      IF (NOPR .EQ. 1) WRITE(6,430) SQ,TAU1,TAU2                                
 430  FORMAT(2X,'ISOTROPIC SQ,TAU1,TAU2',2X,1P,3(2X,E10.3))                     
      OM2=OM                                                                    
      CALL GAMMA(OM2,GAM)                                                       
      IF (NOPR .EQ. 1) WRITE(6,356) TEMP,TAU1,TAU2,OM,BETA,GAM                  
C     NOTE:THE RESULT IS MULTIPLIED BY 1.0E36                                   
      TERM6=SQ*GAM                                                     
      IF (NOPR .EQ. 1) WRITE(6,357) SQ,GAM,TERM6                                
 357  FORMAT(2X,'ISOTROPIC SQ,GAM,TERM6 ',1P,6(2X,E10.3))                       
C                                                                               
C     THE ANISOTROPIC TERMS ARE CALCULATED HERE.                                
C     THE SIX TERMS OF EQUATION 7 OF NBS 80 ARE CALCULATED.                     
      SQ=SQB                                                                    
      TAU1=TAU1B                                                                
      TAU2=TAU2B                                                                
      IF (NOPR .EQ. 1) WRITE(6,431) SQ,TAU1,TAU2                                
 431  FORMAT(2X,'ANISOTROPIC SQ,TAU1,TAU2',2X,1P,3(2X,E10.3))                   
C     THE FIRST TERM IS EVALUATED.                                              
      OM2=OM                                                                    
      CALL GAMMA(OM2,GAM)                                                       
      IF (NOPR .EQ. 1) WRITE(6,356) TEMP,TAU1,TAU2,OM,BETA,GAM                  
      TERM1=AFAC*GAM                                                            
C         ZERO AND TWO ARE EVALUATED.                                           
      CALL EVAL(0,2,OM,B1,B2)                                                   
      TERM2=(RHO(1)*B1)+(RHO(3)*B2)                                             
C         ONE AND THREE ARE EVALUATED.                                          
      CALL EVAL(1,3,OM,B1,B2)                                                   
      TERM3=(RHO(2)*B1)+(RHO(4)*B2)                                             
      TERM3=(9.0/5.0)*TERM3                                                     
C         TWO AND FOUR ARE EVALUATED.                                           
      CALL EVAL(2,4,OM,B1,B2)                                                   
      TERM4=(RHO(3)*B1)+(RHO(5)*B2)                                             
      TERM4=(18.0/7.0)*TERM4                                                    
C         THREE AND FIVE ARE EVALUATED.                                         
      CALL EVAL(3,5,OM,B1,B2)                                                   
      TERM5=(RHO(4)*B1)+(RHO(6)*B2)                                             
      TERM5=(10.0/3.0)*TERM5                                                    
C                                                                               
      SUM3=TERM1+TERM2+TERM3+TERM4+TERM5                                        
C     THE ANISOTROPIC TERMS ARE MULTIPLED BY SQB                                
C     NOTE:THE RESULT IS MULTIPLIED BY 1.0E36                                   
      SUM3=SQ*SUM3                                                     
C     ADD UP THE ISOTROPIC AND ANISOTROPIC TERMS.                               
      SUM3=SUM3+TERM6                                                           
      IF (NOPR .EQ. 1) WRITE(6,460) TERM1,TERM2,TERM3,TERM4,TERM5,              
     1  SUM3                                                                    
 460  FORMAT(2X,'TERM1,TERM2,TERM3,TERM4,TERM5,TOTAL SUM*SQ',/,                 
     1  2X,1P,7(2X,E10.3))                                                      
C     NOTE:THE RESULT IS MULTIPLIED BY 1.0E36                                   
      SUM=AW2*SUM3*2.9979E-4        
C                                                                               
C       NOTICE FACTOR OF 2.  NEED TO RESOLVE THIS DISCREPANCY                   
C       WITH BIRNBAUM. (FACTOR OF 2 GIVES D(V) VALUES WHICH                     
C       AGREE WITH BIRNBAUM'S GRAPHED VALUES.)                                  
C                                                                               
C     B(T) IN THE TRAFTON STYLE IS PUT INTO B                                   
      B=SUM                                                                     
      IF (NOPR .EQ. 1) WRITE(6,404) V,T,B                                       
 404  FORMAT(2X,'V,T,B = ',2X,1P,3(2X,E10.3),/)                                 
C                                                                               
C                                                                               
      RETURN                                                                    
      END                                      
C Subroutine for H2H2 absorption from Joiner and Steffes
	SUBROUTINE H2H2JS(VNU,T,AH2H2)
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V3/QROT,QR1,QR2,S1,S2,S3,S4,P1,P2,P3,P4                            
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
      COMMON/V6/P5,P6,P7,P8                                                     
      COMMON/V7/PRO1(20),PRO2(20),STOR(20)                                      
      COMMON/W1/FP,REB,RHO(20)                                                  
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN                           
C                                                                               
	X1=(273.0/T)**3.12
	X2=(273.0/T)**2.24
	X3=(273.0/T)**3.34
	X=PH2*X1 + 1.382*PHE*X2 + 9.322*PCH4*X3
      FREQ=C*VNU*1.0E9         	
      ZL=1.0/VNU                
	Y=3.557E-11*PH2/(ZL*ZL)	
	AH2H2=Y*X
	RETURN
	END
C Following subroutine is based on Orton's tables
C Following subroutine is based on Orton's tables
	SUBROUTINE H2H2OR(IQ,VNU,T,A,B,CZB)
      dimension Ttab(10),frequ(2428)
      dimension ah2h2e(10),ah2h2n(10),ah2hee(10),ah2hen(10)
      dimension ah2che(10),ah2chn(10)
      dimension abeh2h2(2428,10),abnh2h2(2428,10)
      dimension abeh2he(2428,10),abnh2he(2428,10)
      dimension abeh2ch4(2428,10),abnh2ch4(2428,10)
C      real*8 A,B,CZB,ah2e,ah2n,ahee,ahen,ache,achn
c      integer ntemp,nfreq
      COMMON/orton/Ttab,frequ,abeh2h2,abnh2h2,abeh2he,abnh2he,
     + abeh2ch4,abnh2ch4,nfreq,ntemp
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN                           
      COMMON/W1/FP,REB,RHO(20)                                                  
C  
      algT=alog(T)
C  CHANGE PRESSURES FROM BARS TO ATM, AND BACK AT END OF SUBROUTINE; 
C      
        kp=0
        do 10 ip=1,nfreq
           if(vnu.ge.frequ(ip)) kp=kp+1
           if(vnu.lt.frequ(ip)) goto 11
 10        continue
C  so vnu is in between kp-1 and kp (after next line)
 11         kp=kp+1
            x=vnu-frequ(kp-1)
           do i=1,ntemp
            ah2h2e(i)=abeh2h2(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abeh2h2(kp,i)-abeh2h2(kp-1,i))
            ah2h2n(i)=abnh2h2(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abnh2h2(kp,i)-abnh2h2(kp-1,i))
            ah2hee(i)=abeh2he(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abeh2he(kp,i)-abeh2he(kp-1,i))
            ah2hen(i)=abnh2he(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abnh2he(kp,i)-abnh2he(kp-1,i))
            ah2che(i)=abeh2ch4(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abeh2ch4(kp,i)-abeh2ch4(kp-1,i))
            ah2chn(i)=abnh2ch4(kp-1,i)
     +      +x/(frequ(kp)-frequ(kp-1))*(abnh2ch4(kp,i)-abnh2ch4(kp-1,i))
           enddo
C  
           kp=0
        do 12 ip=1,ntemp
           if(algT.ge.Ttab(ip)) kp=kp+1
           if(algT.lt.Ttab(ip)) goto 13
 12     continue
C  so algT is in between kp-1 and kp
 13            kp=kp+1
               x=algT-Ttab(kp-1)
C
              ah2e=ah2h2e(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2h2e(kp)-ah2h2e(kp-1))
              ah2n=ah2h2n(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2h2n(kp)-ah2h2n(kp-1))
              ahee=ah2hee(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2hee(kp)-ah2hee(kp-1))
              ahen=ah2hen(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2hen(kp)-ah2hen(kp-1))
              ache=ah2che(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2che(kp)-ah2che(kp-1))
              achn=ah2chn(kp-1)
     +          +x/(Ttab(kp)-Ttab(kp-1))*(ah2chn(kp)-ah2chn(kp-1))
C
        if(IQ.ge.1) A=exp(ah2e)
        if(IQ.ge.1) B=exp(ahee)
        if(IQ.ge.1) CZB=exp(ache)
        if(IQ.eq.0) A=exp(ah2n)
        if(IQ.eq.0) B=exp(ahen)
        if(IQ.eq.0) CZB=exp(achn)
C
C        write(6,667) (Ttab(i),i=1,10)
C        write(6,668) vnu,algt,nfreq,ntemp
C        write(6,669) abeh2h2(128,3),frequ(128)
C        write(6,666) T,A,B,CZB,ah2e,ahee,ache
 666    format(7E15.6)
 667    format(10f8.3)
 668    format(2f8.3,2I10)
 669    format(2E15.5)
	RETURN
	END
C                                                                               
C                                                                               
      SUBROUTINE EVAL(L,K,OM,B1,B2)                                             
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
C     CALCULATE THE GAMMA VALUES.                                               
C     SEE NBS 80.                                                               
      CALL EH2(L,E)                                                             
      B1=E                                                                      
      CALL EH2(K,E)                                                             
      B2=E                                                                      
C     B3 IS W(L,K) IN 2 PI RAD PER SEC.                                         
      B3=TPC*(B2-B1)                                                            
      OM2=OM-B3                                                                 
      IF (NOPR .EQ. 1) WRITE(6,10) L,K,OM,B3                                    
  10  FORMAT(2X,I3,2X,1P,I3,2X,1P,2(2X,E10.3))                                  
      CALL GAMMA(OM2,GAM)                                                       
      IF (NOPR .EQ. 1) WRITE(6,15) TEMP,TAU1,TAU2,OM2,BETA,GAM                  
  15  FORMAT(2X,1P,6(2X,E10.3))                                                 
      B1=GAM                                                                    
      OM2=OM+B3                                                                 
      CALL GAMMA(OM2,GAM)                                                       
      IF (NOPR .EQ. 1) WRITE(6,15) TEMP,TAU1,TAU2,OM2,BETA,GAM                  
      B2=GAM                                                                    
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE PARAM(M)                                                       
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
C     SEE NBS 80 FOR THE SQ,TAU1,TAU2 PARAMETER VALUES.                         
      IF (M .EQ. 1) GO TO 10                                                    
      IF (M .EQ. 2) GO TO 20                                                    
      IF (M .EQ. 3) GO TO 30                                                    
C     FOR H2-H2 COLLISIONS.                                                     
  10  A1=TEMP/273.15                                                            
C      A2=142.0*EXP(0.26*ALOG(A1))                                               
C      A3=4.85*EXP(-0.593*ALOG(A1))                                              
C      A4=2.17*EXP(-0.523*ALOG(A1))                                              
      A2=281.0*(A1**0.235)                                                      
      A3=4.68*(A1**(-0.605))                                                    
      A4=2.23*(A1**(-0.607))                                                    
      SQ=1.38*A2                                                    
      TAU1=A3                                                           
      TAU2=A4                                                           
      GO TO 40                                                                  
C     FOR H2-HE COLLISIONS,THE ISOTROPIC PARAMETERS.                            
  20  CONTINUE
C     A1=TEMP/273.15                                                            
C      A2=112.0*EXP(0.93*ALOG(A1))                                               
C      A3=1.74*EXP(-0.54*ALOG(A1))                                               
C      A4=3.4*EXP(-0.30*ALOG(A1))                                                
      A1=77.4/TEMP                                                              
      A2=33.53/A1                                                               
      A3=3.43*SQRT(A1)                                                          
      A4=6.56*SQRT(A1)                                                          
      SQ=1.38*A2                                                    
      TAU1=A3                                                           
      TAU2=A4                                                           
      GO TO 40                                                                  
C     FOR H2-HE COLLISIONS,ANISOTROPIC PARAMETERS.                              
C     VALUES ARE FOR 195 DEGREES KELVIN.                                        
  30  CONTINUE
C     A1=TEMP/273.15                                                            
C      A2=1.914E1                         
C      A3=4.13                                    
C      A4=2.38                                         
      A1=TEMP/77.4                                                              
      A2=12.06*(A1**0.57)                                                       
      A3=3.02*(A1**(-0.30))                                                     
      A4=8.94*(A1**(-0.60))                                                     
      SQ=1.38*A2                                                    
      TAU1=A4                                                           
      TAU2=A3                                                           
      GO TO 40                                                                  
  40  RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE FAC(L,AFAC)                                                    
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V7/PRO1(20),PRO2(20),STOR(20)                                      
      COMMON/W1/FP,REB,RHO(20)                                                  
C     SEE NBS 80 FOR THE AFAC FACTOR,EQUATION 8.                                
      A8=0.0                                                                    
      DO 10 I=1,L                                                               
      J=I-1                                                                     
      A1=J                                                                      
      A2=A1+1.00                                                                
      A3=2.0*A1                                                                 
      A4=A3+1.00                                                                
      A5=A3-1.00                                                                
      A6=A3+3.00                                                                
      A7=(A1*A2*A4*RHO(I))/(A5*A6)                                              
      A8=A8+A7                                                                  
  10  CONTINUE                                                                  
      AFAC=A8                                                                   
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE BESSEL(Z,BES)                                                  
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
C     BES IS Z*K1(Z)                                                            
C     SEE PAGES 378-379 OF ABRAMOWITZ AND STEGUN,HANDBOOK OF                    
C     MATHEMATICAL FUNCTIONS.                                                   
C     EQNS. 9.8.3,9.8.4,9.8.7,9.8.8 ARE USED TO CALCULATE Z*K1(Z)               
      IF (Z .LE. 2.00) GO TO 10                                                 
      IF (Z .GT. 2.00) GO TO 20                                                 
C     A1,A2,..ARE Z/2 TO THE 2,4,..POWER. SEE EQN 9.8.7.                        
C     BESI IS I1(Z).                                                            
  10  CALL BESS2(Z,BESI)                                                        
      A1=(Z*Z)/(2.0*2.0)                                                        
      A2=A1*A1                                                                  
      A3=A2*A1                                                                  
      A4=A3*A1                                                                  
      A5=A4*A1                                                                  
      A6=A5*A1                                                                  
      BES=(Z*ALOG(Z/2.0)*BESI)+(1.000)+(.15443144*A1)                           
      BES=BES-(.67278579*A2)-(.18156897*A3)-(.01919402*A4)                      
      BES=BES-(.00110404*A5)-(.00004686*A6)                                     
      GO TO 30                                                                  
C     A1,A2,.. ARE Z/2 TO THE 1,2,.. POWER. SEE EQN 9.8.8.                      
  20  A1=2.0/Z                                                                  
      A2=A1*A1                                                                  
      A3=A2*A1                                                                  
      A4=A3*A1                                                                  
      A5=A4*A1                                                                  
      A6=A5*A1                                                                  
      BES=1.25331414+(.23498619*A1)-(.03655620*A2)                              
      BES=BES+(.01504268*A3)-(.00780353*A4)                                     
      BES=BES+(.00325614*A5)-(.00068245*A6)                                     
      CALL EX(Z,A1)                                                             
      BES=BES*SQRT(Z)*A1                                                        
  30  RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE BESS2(Z,BESSI)                                                 
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
C     BESSI IS I1(Z).                                                           
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      IF (Z .LE. 3.75) GO TO 10                                                 
      IF (Z .GT. 3.75) GO TO 20                                                 
C     A1,A2,.. IS Z/3.75 TO THE 2,4,.. POWER. SEE EQN. 9.8.3.                   
  10  A1=(Z*Z)/(3.75*3.75)                                                      
      A2=A1*A1                                                                  
      A3=A2*A1                                                                  
      A4=A3*A1                                                                  
      A5=A4*A1                                                                  
      A6=A5*A1                                                                  
      BESSI=0.5000+(.87890594*A1)+(.51498869*A2)                                
      BESSI=BESSI+(.15084934*A3)+(.02658733*A4)                                 
      BESSI=BESSI+(.00301532*A5)+(.00032411*A6)                                 
      BESSI=BESSI*Z                                                             
      GO TO 30                                                                  
C     A1,A2,.. ARE 3.75/Z TO THE 1,2,.. POWER. SEE EQN. 9.8.4.                  
  20  A1=3.75/Z                                                                 
      A2=A1*A1                                                                  
      A3=A2*A1                                                                  
      A4=A3*A1                                                                  
      A5=A4*A1                                                                  
      A6=A5*A1                                                                  
      A7=A6*A1                                                                  
      A8=A7*A1                                                                  
      BESSI=.39894228-(.03988024*A1)-(.00362018*A2)                             
      BESSI=BESSI+(.00163801*A3)-(.01031555*A4)                                 
      BESSI=BESSI+(.02282967*A5)-(.02895312*A6)                                 
      BESSI=BESSI+(.01787654*A7)-(.00420059*A8)                                 
      BESSI=BESSI*(EXP(Z)/SQRT(Z))                                              
      IF ((Z .GT. 150.) .AND. (NOPR .EQ. 1)) WRITE(6,29)                        
  29  FORMAT(2X,'Z .GT. 150.,EXP PROBLEMS')                                     
  30  RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE GAMMA(OM2,GAM)                                                 
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
      COMMON/V9/TEMP,BETA,SQ,TAU1,TAU2                                          
C     SEE EQNS 17 OF BIRNBAUM AND COHEN (CAN J PHYS,54,PG 593,1976)             
C     ALSO EQUATION 5 OF NBS 80.                                                
      A1=TAU1                                                                   
      A2=TAU2                                                                   
      A3=0.5*BETA*HB*OM2                                                        
      A4=1.00+(OM2*OM2*A1*A1*1.0E-6)      
      A5=(A2*A2*1.0E-6)+((BETA*BETA*HB*HB)/4.0)         
      GAM=(A1/(PI*A4))*EXP((A2/A1)+A3)                                          
        A6=GAM                                                                  
      Z=SQRT(A4*A5)*1.0E3/A1                                
      CALL BESSEL(Z,BES)                                                        
C  GAM SHOULD BE MULTIPLIED BY 1.0E-14
      GAM=GAM*BES                                                               
      IF (NOPR .EQ. 1) WRITE(6,10) A1,A2,A3,A4,A5,A6,Z,GAM                      
  10  FORMAT(2X,'TAU1,TAU2,A3,A4,A5,A6,Z,GAM',/,                                
     1  2X,1P,8(2X,E10.3))                                                      
      RETURN                                                                    
      END                                                                       
      SUBROUTINE EH2(J,E)                                                       
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
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
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE JROT(TEMP)                                                     
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4                                          
      COMMON/V7/PRO1(20),PRO2(20),STOR(20)                                      
      COMMON/W1/FP,REB,RHO(20)                                                  
C     THE BOLTZMANN STATISTICS FOR H2 ARE CALCULATED                            
      CALL PART(TEMP)                                                           
C      RHO(L) GIVES THE FRACTION OF H2 IN THE                                   
C     K=L-1 ROTATIONAL LEVEL,DIVIDED BY 2J+1.                                   
C     AND MULTIPLIED BY FP FOR THE NORMAL CASE,AND                              
C     BY 1.00 IN THE EQUIL CASE (SEE SUBROUTINE                                 
C     PART FOR DEFINITION OF STOR).                                             
C     FRACT IS THE FRACTION OF H2 IN THE K=L-1 LEVEL.                           
C     NOTICE SEPARATE EQUIL AND NORMAL DEFINITIONS OF FRACT.                    
C     EQUIL CASE,KL=0 AND A2=1.00                                               
C     NORMAL CASE,KL=10 AND A2=FP.                                              
 345  FORMAT(/,2X,'BOLTZMANN STATISTICS FOR H2',/,                              
     1  2X,'KL IS EQUAL TO  ',I3)                                               
C     FOR NORMAL CASE,FRACTION OF PARA IS 0.25                                  
C       PARA LEVELS J=0,2,4..ARE SPECIFIED.                                     
      DO 20 L=1,5                                                               
C  K=1,3,5,7,9                                                                  
      K=((L-1)*2)+1                                                             
      J=K-1                                                                     
      A1=(2*(K-1))+1                                                            
C                                                                               
      POP=(REB*STOR(K))+((1.0-REB)*FP*STOR(K+10))                               
      RHO(K)=POP/A1                                                             
 350  FORMAT(2X,'J,FRACT,PRO1 ',2X,I3,2X,1P,2(2X,E10.3))                        
  20  CONTINUE                                                                  
C       ORTHO LEVELS J=1,3,5.. ARE SPECIFIED.                                   
      DO 21 L=1,5                                                               
      K=((L-1)*2)+2                                                             
      J=K-1                                                                     
      A1=(2*(K-1))+1                                                            
      POP=(REB*STOR(K))+((1.0-REB)*(1.0-FP)*STOR(K+10))                         
      RHO(K)=POP/A1                                                             
  21  CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE PART(TEMP)                                                     
      COMMON/V3/QROT,QR1,QR2,S1,S2,S3,S4,P1,P2,P3,P4                            
      COMMON/V6/P5,P6,P7,P8                                                     
      COMMON/V7/PRO1(20),PRO2(20),STOR(20)                                      
      DIMENSION CC(20),CC2(20)                                                  
C     SEE TATUM,AP J SUPPL SER,VOL 14,PAGE 21,1967                              
C     FOR THE EXPRESSIONS FOR THE BOLTZMANN STATISTICS                          
C     THAT APPLY TO THE GROUND STATE OF H2.                                     
C     HERZBERG,MOLECULAR SPECT AND MOL STRUCT,VOL 4,CONSTS                      
C     OF DIATOMIC MOLECULES IS USED FOR THE BV,DE VALUES.                       
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
      P1=S2/S1                                                                  
      P2=CC2(1)/CC(1)                                                           
      P3=S4/S3                                                                  
      P4=(CC2(1)/CC(1))*(QR1/QR2)                                               
      P5=CC2(1)/QROT                                                            
      P6=CC(1)/QROT                                                             
      P7=CC2(1)/QR2                                                             
      P8=CC(1)/QR1                                                              
      DO 40 J=1,LR                                                              
      PRO1(J)=CC(J)/QROT                                                        
      PRO2(J)=CC2(J)/QROT                                                       
      PRO1(J+10)=CC(J)/QR1                                                      
      PRO2(J+10)=CC2(J)/QR2                                                     
  40  CONTINUE                                                                  
C     STOR(I) I=1,2,3,4.. FOR JROT=0,1,2,3,..                                   
C     EQUIL DISTRIBUTION FOR I=1 TO 10                                          
C     NORMAL DISTRIBUTION FOR I=11 TO 20                                        
      DO 50 J=1,LR                                                              
      L=((J-1)*2)+1                                                             
      L1=L+1                                                                    
C     AA1=1.00,AA2=1.00 FOR STANDARD BOLTZMANN STATISTICS                       
C     WHERE THE DEGENERACY IS UPPSTAIRS AND DOWNSTAIRS.                         
C     AA1=.25*L,AA2=.75*L1 FOLLOWS EQUATION 9B OF                               
C     BIRNBAUM AND COHEN(WRONG?)                                                
      AA1=1.00                                                                  
      AA2=1.00                                                                  
      STOR(L)=PRO1(J)/AA1                                                       
      STOR(L1)=PRO2(J)/AA2                                                      
      STOR(L+10)=PRO1(J+10)/AA1                                                 
      STOR(L1+10)=PRO2(J+10)/AA2                                                
  50  CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
      SUBROUTINE DIEL(I,NPR,W,TK,E1,E2,D1,D2)
	COMPLEX ZM,ZZ
	REAL*4 LF,LF0,LF1
	DIMENSION EIMAG(9),FR(9)
C  SUBROUTINE TO CALCULATE DIELECTRIC CONST/REFRACTIVE INDICES FOR
C  WATER, BASED UPON ULABY ET AL, 1981.
	F=30.0E9/W
      IF(I.EQ.1) GOTO 100
      IF(I.EQ.2) GOTO 200
C  WATER
  100 TEMP=TK-273.0
      A1=TEMP
      A2=A1*A1
      A3=A2*A1
	RELT=(1.1109E-10)-A1*(3.824E-12)+A2*(6.938E-14)-
     + A3*(5.096E-16)
      E0=88.045-0.4147*A1+A2*(6.295E-4)+A3*(1.075E-5)
	if(E0.LT.0.0) E0=0.0
      EINF=4.9
      E1=(E0-EINF)/(1+(F*RELT)*(F*RELT))
	E1=EINF+E1
	E2=F*RELT*(E0-EINF)
	E2=E2/(1+(F*RELT)*(F*RELT))
	IF(E2.LT.0.0) E2=0.0
      GOTO 900
C  ICE
  200  CONTINUE
	TEMP=TK-273.0
	FR(1)=1.0E8
	FR(2)=3.0E8
	FR(3)=1.0E9
	FR(4)=2.0E9
	FR(5)=3.0E9
	FR(6)=5.0E9
	FR(7)=1.0E10
	FR(8)=3.0E10
	FR(9)=1.0E11
	EIMAG(1)=8.0E-3
	EIMAG(2)=1.5E-3
	EIMAG(3)=8.0E-4
	EIMAG(4)=1.0E-3
	EIMAG(5)=1.2E-3
	EIMAG(6)=1.5E-3
	EIMAG(7)=3.0E-3
	EIMAG(8)=8.0E-3
	EIMAG(9)=2.0E-2
	E1=3.15
	LF=LOG10(F)
	IF(F.LE.FR(1)) THEN
	E2=8.0E-3
	GOTO 900
	ENDIF
	IF(F.GE.FR(9)) THEN
	E2=2.0E-2
        GOTO 900
	ENDIF
	IF(F.GT.FR(1).AND.F.LE.FR(2)) J=1
	IF(F.GT.FR(2).AND.F.LE.FR(3)) J=2
	IF(F.GT.FR(3).AND.F.LE.FR(4)) J=3
	IF(F.GT.FR(4).AND.F.LE.FR(5)) J=4
	IF(F.GT.FR(5).AND.F.LE.FR(6)) J=5
	IF(F.GT.FR(6).AND.F.LE.FR(7)) J=6
	IF(F.GT.FR(7).AND.F.LE.FR(8)) J=7
	IF(F.GT.FR(8).AND.F.LE.FR(9)) J=8
	LF0=ALOG10(FR(J))
	LF1=ALOG10(FR(J+1))	
	DLF=(LF-LF0)/(LF1-LF0)
	X0=ALOG10(EIMAG(J))
	X1=ALOG10(EIMAG(J+1))
	DX=X0+DLF*(X1-X0)
	E2=10**DX
  900   CONTINUE
c	ZM=CMPLX(E1,-E2)
c	ZZ=SQRT(ZM)
c	D1=REAL(ZZ)
c	D2=-AIMAG(ZZ)
      D1=(E1+SQRT((E1*E1)+(E2*E2)))/2.0
      D1=SQRT(D1)
      IF(D1.EQ.0.0) GOTO 910
      D2=E2/(2.0*D1)
      GOTO 9999
  910 D2=0
 9999 RETURN
      END
C
C
C
       FUNCTION SHAPE2(P,T,F,FLC)
        implicit none
        REAL*4 SHAPE2
        real*4 DL,P,T,F,FLC,DFC,DL2
        DL=F-FLC
	DL2=F+FLC
c  line width is taken as water line width
	DFC = 0.08*P*((273/T)**0.666) 
	DFC = DFC * 30.0		!LINEWIDTH IN MHZ
       SHAPE2=((F*F)/(FLC*FLC))*(DFC/(DL*DL+DFC*DFC)
     + +DFC/(DL2*DL2+DFC*DFC))
	SHAPE2=SHAPE2/3.1415926536		!DIVIDE BY PI
       RETURN
       END
C
       REAL*4 FUNCTION QLOG2(Q0,T)
       IMPLICIT NONE
       REAL*4 Q0(5),T
       REAL*4 SLOPE, TLOG
C LOG(T) FOR T=300,225,150,75, 37.5K
       REAL*4 TLOG0(5) /2.47712,2.35218,2.17609,1.8751, 1.5740/
       TLOG=ALOG10(T)
       IF(TLOG.GT.TLOG0(1)                         )GO TO 4
       IF(TLOG.LE.TLOG0(1) .AND. TLOG.GT. TLOG0(2))GO TO 5 
       IF(TLOG.LE.TLOG0(2) .AND. TLOG.GT. TLOG0(3))GO TO 6 
       IF(TLOG.LE.TLOG0(3) .AND. TLOG.GT. TLOG0(4))GO TO 7
       IF(TLOG.LE.TLOG0(4) .AND. TLOG.GT. TLOG0(5))GO TO 8
       IF(TLOG.LE.TLOG0(5)                        )GO TO 9
4      QLOG2=Q0(1)    !USE 300 K VALUE NEEDS TO BE FIXED
       GO TO 10
5      SLOPE=(Q0(2)-Q0(1))/(TLOG0(2)-TLOG0(1))
       QLOG2=Q0(1)+SLOPE*(TLOG-TLOG0(1))
       GO TO 10
6      SLOPE=(Q0(3)-Q0(2))/(TLOG0(3)-TLOG0(2))
       QLOG2=Q0(2)+SLOPE*(TLOG-TLOG0(2)) 
       GO TO 10
7      SLOPE=(Q0(4)-Q0(3))/(TLOG0(4)-TLOG0(3))
       QLOG2=Q0(3)+SLOPE*(TLOG-TLOG0(3)) 
       GO TO 10
 8     SLOPE=(Q0(5)-Q0(4))/(TLOG0(5)-TLOG0(4))
       QLOG2=Q0(4)+SLOPE*(TLOG-TLOG0(4)) 
       GO TO 10
 9     QLOG2=Q0(5)              !USE 37.5 K VALUE -NEEDS TO BE FIXED
10     CONTINUE
       RETURN
       END
C
C
      SUBROUTINE SAJph3(N,T,T3,GM)                               
        REAL*4 T,P,F,V
	real*8 x,ajk,fjk,tcor
	real*4 T3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
	real*4  FP0(728),SP0(728),ELP(728)
	REAL*4 wgtS0(40),wgtFGB(40),wgtSB(40)
	integer*2 N
        COMMON/ph3/FP0,SP0,ELP,wgtS0,wgtFGB,wgtSB
        COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
C     GM IS GAMMA(N) IN GHZ.                                                  
	GH2=3.293
	GHE=1.6803
	Gph3=4.2157
	if(N.gt.40) goto 40
        GM=(GH2*T3*PH2)+(GHE*T3*PHE)
	GM=GM*wgtfgb(N)
	GM=GM+(GPH3*(300.0/T)*PPh3*wgtSB(N))                         
	goto 50
  40   GM=(GH2*T3*PH2)+(GHE*T3*PHE)
	GM=GM+(GPH3*(300.0/T)*PPh3)                         
  50  continue
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE SFJPH3(N,T,V,VLINE,GM,FJK)                                
        REAL*4 T,P,F,V
	real*8 x,ajk,fjk,tcor,y,z
        real*4 t3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
	real*4  FP0(728),SP0(728),ELP(728)
	REAL*4 wgtS0(40),wgtFGB(40),wgtSB(40)
	integer*2 N
        COMMON/ph3/FP0,SP0,ELP,wgtS0,wgtFGB,wgtSB
        COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
C                                                                               
	zeta=0.
	delta=0.0
C     FJK IS THE BEN REUVEN LINE SHAPE, FROM DEBOER AND STEFFES
C                                                                               
	x = (gm+zeta)*((vline+delta)*(vline+delta) + GM*GM - zeta*zeta)
	x=x+(gm-zeta)*v*v
        y = (v*v - (vline+delta)*(vline+delta) - gm*gm + zeta*zeta)
	z=y*y  + 4.0*v*v*GM*GM
        FJK  = 29.97925*2.0*v*v*x/3.14159/z/vline/vline
      RETURN                                                                    
      END                                                                       
C
       SUBROUTINE ABSBR_h2s2(AON,T,P,VNU)
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF H2S, IN UNITS OF 1/cm
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/cm
C	T	TEMPERATURE				      K
C    	P	PRESSUME				BARS
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V
	real*4 t3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
        real*4 DL,FLC,DL2
	REAL*4 C,DV,VNU,coef,hck
	real*8 x,ajk,fjk,tcor
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=311)   !NUMBER OF SPECTRAL LINES
	real*4  F0(311),S0(311),EL(311),Q0(5),wgths(311)
	REAL*4 QLOG2,SHAPE2
        EXTERNAL QLOG2
        DATA Q0 /2.711, 2.5244, 2.2619, 1.8164, 1.3778/ !LOG PF AT 300,225,150,75,37.5K
        COMMON/h2s/F0,S0,EL,wgths
         COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
C
	C=2.997925E1
	coef=7.244E21
	hck=1.438396	!hc/k
	V=VNU*C				!FREQ IN GHZ
	T3=(296./T)**0.6667
	AON=0.0
        DO N = 1,NLINES
	VLINE=F0(N)		!LINE CENTER IN GHZ
        TCOR = EXP(-(hck*EL(N)*(1./T - 3.3784E-03)))
	call sajh2s2(N,T,T3,GM)
        ajk=S0(N)*TCOR
        CALL SFJh2s2(N,T,V,VLINE,GM,FJK)
       AON=AON + ajk*FJK
       ENDDO
	x=(296.0/T)**3.5
	AON=7.244E21*x*AON*Ph2s/296.0		!2.4163E3 is from deboer
C
       RETURN
       END
C
      SUBROUTINE SAJh2s2(N,T,T3,GM)                               
        REAL*4 T,P,F,V
	real*8 x,ajk,fjk,tcor
	real*4 T3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
	real*4  F0(311),S0(311),EL(311),wgths(311)
	integer*2 N
       COMMON/h2s/F0,S0,EL,wgths
        COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
C     GM IS GAMMA(N) IN GHZ.
 	GH2=1.96
	GHE=1.20
	GH2S=wgths(N)
      GM=(GH2*T3*PH2)+(GHE*T3*PHE)+(GH2S*T3*PH2S)                         
      RETURN                                                                    
      END                                                                       
C                                                                               
      SUBROUTINE SFJh2s2(N,T,V,VLINE,GM,FJK)                                
        REAL*4 T,P,F,V
	real*8 x,ajk,fjk,tcor,y,z
        real*4 t3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
	real*4  F0(311),S0(311),EL(311),wgths(311)
	integer*2 N
        COMMON/h2s/F0,S0,EL,wgths
        COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
C                                                                               
	zeta=gm
	delta=1.28*ph2s
C     FJK IS THE BEN REUVEN LINE SHAPE, FROM DEBOER AND STEFFES
C                                                                               
	x = (gm+zeta)*((vline+delta)*(vline+delta) + GM*GM - zeta*zeta)
	x=x+(gm-zeta)*v*v
        y = (v*v - (vline+delta)*(vline+delta) - gm*gm + zeta*zeta)
	z=y*y  + 4.0*v*v*GM*GM
        FJK  = 29.97925*2.0*v*v*x/3.14159/z/vline/vline
      RETURN                                                                    
      END                                                                       
C
      SUBROUTINE NH3JOIN(VNU,T,ANH3)            
C     CALCULATE THE AMMONIA ABSORPTION COEFFICIENT (CM-1)                       
C     USED BY GULKIS AND POYNTER (1972),AND GIVEN IN BERGE                      
C     AND GULKIS, GEHRELS JUPITER BOOK,PAGE 675.                                
      DIMENSION VKAK(16,16),DVKAK(16,16)                                        
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
      COMMON/U2/C1,C2,C3,C4                                                     
      COMMON/U3/VKAK,DVKAK                                                      
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V2/PI,HB,C,TPC,BOLC,AW,DEN,BAR,G                                   
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4           
C
C   VKAK(J,K)=CENTER FREQU FOR THE (J,K) TRANSITION ACCORDING TO THE            
C  TABULATED VALUES OF POYNTER AND GULKIS (1975).                               
C  DVKAK(J,K) ARE THE CORRESPONDING SELF BROADENED LINE WIDTHS                  
C  AS TABULATED BY THEM (RESP. GHZ, MHZ/TORR)                                   
C     CHANGE VNU(CM-1) TO GHZ.                                                  
      V=(C*VNU)                                                          
C  CHANGE PRESSURES FROM BARS TO ATM, AND BACK AT END OF SUBROUTINE; MARCH89
	PH2=PH2/1.013
	PHE=PHE/1.013
	PH2O=PH2O/1.013
	PNH3=PNH3/1.013
	PCH4=PCH4/1.013
C     CORRECTION FRACTOR IS GIVEN ON PAGE 676 OF THE BERGE AND GULKIS           
C     PAPER, GEHRELS JUPITER BOOK.                                              
C                                                                               
	CORR=1.0
      IF (NOPR .EQ. 1) WRITE(6,5) V,CORR                                        
   5  FORMAT(2X,'NH3 V(GHZ),CORR',2X,1P,2(2X,E10.3))                            
      IF (NOPR .EQ. 1) WRITE(6,6) C1,C2,C3,C4                                   
   6  FORMAT(2X,'CONSTANTS C1,C2,C3,C4',2X,1P,4(2X,E10.3))                      
C                                                                               
      SUM=0.0                                                                   
C     CALCULATE THE INVERSION LINE WAVENUMBERS                                  
C     FOR NH3, FROM THE DATA IN POYNTER AND KAKAR                               
C     AP J SUPPL SERIES,VOL 29,P87,1975.                                        
C     J AND K ARE QUANTUM NUMBERS.                                              
C                                                                               
C     NLINE ARE THE NUMBER OF LINES TO BE CALCULATED.                           
      NLINE=16                                                                  
C                                                                               
      DO 10 J=1,NLINE                                                           
      L1=J*(J+1)                                                                
      L3=0                                                                      
      DO 20 K=1,J                                                               
      L2=K*K                                                                    
      L3=L3+1                                                                   
C     SB IS THE SELF BRODENED LINE WIDTH OF NH3 IN MHZ/TORR.                    
C     C4 IS GIVEN IN THE BLOCK DATA STATEMENT.                                  
C     SK IS 1.5 OR 3, DEPENDING UPON THE K VALUE.                               
      SK=1.5                                                                    
      IF (L3 .NE. 3) GO TO 21                                                   
      SK=3.0                                                                    
      L3=0                                                                      
  21  VLINE=VKAK(J,K)                                                           
       SB=DVKAK(J,K)                                                            
      IF (NOPR .EQ. 1) WRITE(6,22) J,K,SK,VLINE,SB                              
  22  FORMAT(2X,'J,K,SK,VLINE(GHZ),SB(MHZ/TORR)',/,                             
     1 2X,I3,2X,I3,2X,F4.1,2X,1P,2(2X,E10.3))                                   
C     CALCULATE A(J,K).                                                         
      CALL SAJKJOIN(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)  
      IF(J.EQ.1.AND.K.EQ.1) ZZV=2.0*AJK*EXP(23.3/T)/3.0                         
C     FJK IS THE FREQUENCY DEPENDENT LINE SHAPE OF BEN REUVEN.                  
      CALL SFJKJOIN(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      SB=AJK*FJK                                                                
      SUM=SUM+(CORR*AJK*FJK)                                                         
      IF (NOPR .EQ. 1) WRITE(6,23) AJK,FJK,SB                                   
  23  FORMAT(2X,'AJK,FJK,(AJK*FJK)',2X,1P,3(2X,E10.3))                          
  20  CONTINUE                                                                  
  10  CONTINUE                                                                  
C                                                                               
C  CALCULATE ROTATION LINE J=0;K=0                                              
C  VLINE=5.73E5 MHZ; CHOSE DVLINE=14. MHZ                                       
C                                                                               
      L1=2                                                                      
      L2=1                                                                      
      VLINE=5.73E2                                                              
      SB=14.0                                                                   
      GM=(1.89*T3*PH2)+(0.75*T3*PHE)+(0.60*T2*SB*PNH3)                         
      CALL SFJKJOIN(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
      ZZV1=ZZV*7.0E-2*T*(1.0-EXP(-28.6/T))                                      
      ZZV2=ZZV*5.25E-2*T*EXP(-23.3/T)*(1.0-EXP(-57.2/T))                        
      ZZV3=ZZV*7.50E-2*T*EXP(-28.6/T)*(1.0-EXP(-57.2/T))                        
      SUM=SUM+(ZZV1*FJK)                                                        
      L1=2                                                                      
      L2=1                                                                      
      VLINE=1169.0                                                              
      SB=23.8                                                                   
      CALL SFJKJOIN(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                      
C   multiply fjk by 4.0 to test sensitivity
C	FJK=4.0*FJK
      SUM=SUM+(ZZV2*FJK)+(ZZV3*FJK)                                             
C                                                                               
C                                                                               
C     THE AMMONIA ABSORPTION COEFFICIENT IS IN CM-1.                            
      ANH3=SUM                                                             
c      IF (NOPR .EQ. 1) WRITE(6,30) CORR,SUM,ANH3                                
  30  FORMAT(2X,'CORR,SUM,ANH3',2X,1P,3(2X,E10.3),/)                            
C                                                                               
	PH2=PH2*1.013
	PHE=PHE*1.013
	PH2O=PH2O*1.013
	PNH3=PNH3*1.013
	PCH4=PCH4*1.013
      RETURN                                                                    
      END                                                                       
C
      SUBROUTINE SAJKJOIN(T,VLINE,J,L1,L2,SK,SB,GM,GH2,AJK)
	real*4 t3,gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
      COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
C     GM IS GAMMA(J,K) IN GHZ.                                                  
      GM=(1.89*T3*PH2)+(0.75*T3*PHE)+(0.60*T2*SB*PNH3)                         
      A1=(2.*J+1.)*FLOAT(L2)/FLOAT(L1)                                          
      A5=((2.98*L1)-(1.09*L2))*T7                                               
C     AJK IS A(J,K),SEE EQN 24 OF BERGE AND GULKIS PAPER,GEHRELS                
C     JUPITER BOOK,PAGE 675.                                                    
      CALL EX(A5,A10)                                                           
      AJK=1.23E3*A1*SK*PNH3*T6*A10                                              
c      IF (NOPR .EQ. 1) WRITE(6,10) GM,A1,A5,AJK                                 
  10  FORMAT(2X,'GM,A1,A5,AJK',2X,1P,4(2X,E10.3))                               
      RETURN                                                                    
      END                                                                       
C
      SUBROUTINE SFJKJOIN(T,L1,L2,V,VLINE,GM,GH2,SB,FJK)                                
	real*4 t3,gm,vline,fjk,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,deltot,PCO,PCO13,PHCN
      COMMON/V1/A1,A2,A3,A4,A5,A6,A7,A8,A9,A10                                  
      COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN
       COMMON/U4/T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11
      COMMON/V4/NOPR,NOPR2,NOPR3,NOPR4
C                                                                               
C     FJK IS THE BEN REUVEN LINE SHAPE,SEE EQN 25 OF BERGE AND                  
C     GULKIS ARTICLE,PAGE 675 OF GEHRELS JUPITER BOOK.                          
C                                                                               
      X=PH2+PHE+PH2O+PNH3+PCH4                                                  
C  IF X.LT.1 ATM VAN VLECK-WEISKOPF LINE SHAPE CAN BE                           
C  USED; HOWEVER, BEN REUVEN IS STILL BETTER OUT IN THE                         
C  WINGS OF THE AMMONIA LINES. 
c     A2 IS THE COUPLING FACTOR ZETA(J,K) IN GHZ.                               
c	GOTO 500
C  for Gross lineshape, goto 500 (first VVW, then gross
      IF(VLINE.GT.100.0) GOTO 500                                               
      A2=(1.35*PH2*T3)+(0.30*T2*PNH3*SB)+                                       
     1  (0.2*T3*PHE)                                                            
C     A4 IS THE PRESSURE SHIFT DEL (GHZ)                                        
	A4=-0.45*PNH3
C     SQUARES OF VARIOUS TERMS.                                                 
      A5=(VLINE+A4)*(VLINE+A4)                                                  
C      A5=VLINE*VLINE
      A6=GM*GM                                                                  
      A7=A2*A2                                                                  
      A8=V*V                                                                    
C     NUMERATOR                                                                 
      A9=((GM-A2)*A8)+((GM+A2)*(A5+A6-A7))                                      
      A10=A8-A5-A6+A7                                                           
C     DENOMINATOR                                                               
      A10=(A10*A10)+(4.0*A8*A6)                                                 
C     THE BEN REUVEN LINE SHAPE.                                                
      FJK=2.0*A8*(A9/A10)                                                       
      IF (NOPR .EQ. 1) WRITE(6,10) A2,A4,A5,A6,A7,A8,A9,                        
     1  A10,FJK                                                                 
	goto 501
  500 CONTINUE                                                                  
C  THE VAN VLECK-WEISKOPF LINE SHAPE                                            
c
      Z1=0.755*PH2+0.231*PHE                                                    
      Z2=FLOAT(L2)/FLOAT(L1)                                                    
      Z1=Z1*(Z2**0.3333)                                                  
      Z2=SQRT(Z2)*PNH3*6.23
      Z=(Z1+Z2)*1.0E3/T                                                         
      X1=V*V*Z                                                                  
      X2=(V-VLINE)*(V-VLINE)                                                    
      X3=(V+VLINE)*(V+VLINE)                                                    
      FJK=X1/(X2+Z*Z)+X1/(X3+Z*Z)   
C      GOTO 501                                                                
C  THE GROSS LINE SHAPE
        X2=(-V*V+VLINE*VLINE)
        X2=X2*X2
       if(vline.lt.1000.) Z1=123.3551*PH2+60.7441*PHE  
	if(vline.ge.1000.) z1=2.885*PH2+1.00*PHE
C      Z2=FLOAT(L2)/FLOAT(L1)                                                    
	Z2=300./T
      Z1=Z1*(Z2**0.3)
	if(vline.lt.1000.) Z=154.1939*PNH3*Z2
	if(vline.ge.1000.) Z=22.989*PNH3*Z2
       Z=Z1+Z
        X3=Z*Z*4.*V*V
        X1=v*VLINE*Z*4.
        X1=X1*V*VLINE
        FJK=X1/(X2+X3)
 501    FFF=FJK                                                                   
  10  FORMAT(2X,'A2,A4,A5,A6,A7,A8,A9,A10,FJK',/,                               
     1  2X,1P,9(2X,E10.3))                                                      
      RETURN                                                                    
      END
c******************************************************************
      function absco(j,fr0,f,p,t,wtmol) 
c
c     this function calculates the absorption coeff of CO
C  f, fr0 in Hz; p in bar, but the routine expects:
C f, fr0 in Hz, and p in microbar
c     normalized by the CO mixing ratio (xco)
c     the formulae were taken from Clancy's PhD thesis
C
      if(j.eq.0) then
          absco = 0.0
          return
      endif
      p=p*1.0E6
c
      pi = 3.14159265308
c
c     offset frequency from line center
      fpm = f - fr0
      fpp = f + fr0
c
c     pressure broadened line width  in Hz
C  this is the same as in the SHAPE function (IdP)
C  note that microbar = dynes/cm^2 
C + this number should be the half-width of CO in H2-broadened gas
C  which is about 0.071 cm-1/atm
C  or 2100 instead of old 3300
C for HCN that would be 2 times larger
C  the following line is pretty consistent with the Menger et al 2000 
C  results for 2-1 transition
       wp = 2100.0 * p * (300.0/t)**0.7

c     doppler width (e-folding ... not FWHM)  in Hz
      wd = fr0 / 3.00e10 * sqrt(2.0*1.38e-16*t/wtmol)
C in his thesis wd is smaller by f. 1.20
c
c     number density of atmosphere (ideal gas p = dkt)
      d = p / 1.38e-16 / t
c
      a = wp/wd
      x = fpm/wd
c
c     use the Voigt function line shape
      shape = voigt(x,a) / wd
c       
c     VanVleck-Weisskopf shape (for future use at large offsets)
      vvwshape = 
     & 1.0 / pi * f*f/fr0/fr0 * (wp/(fpm*fpm+wp*wp)+wp/(fpp*fpp+wp*wp))
c
c     Lorentz Shape
      vl = wp / pi / (fpm*fpm + wp*wp)
c     
c
        shape=shape*1.0E6
        v1=v1*1.0E6
        vvwshape=vvwshape*1.0E6
        fshape=shape
        if(v1.gt.fshape) fshape=v1
C  The following is correct, I think -- 15 April 2011
        if(vvwshape.gt.fshape) fshape=vvwshape
c  I think pin = p/1.013E6 so P in is in atm; it was written as *1.013E6
        pin = p / 1.01325E6
         xshape=fshape
C         xshape=vvwshape
        call absx_coline(j,aon,t,pin,f,xshape)
C  xshape: p in millibar, f in MHz (p/f); I don't think we need to use this
C vvwshape: p in microbar; f in Hz
C We multiply shape by 1E6 to correct answer to come out in bars and MHz
C
         absco = aon
c       write(6,1003) j,aon,absco,shape,f,x,a,wp,wd
c       write(6,1004) p,xshape,vvwshape,v1,shape
c 1003  format(I5,8E12.4)
c 1004  format(5E12.4)
      p=p/1.0E6
C  if other gases are considered, add them here as well I guess; e.g. 
C  SO2, H2SO4 vapor. VVW line shape is probably OK for lower layers
c
      return
      end
c******************************************************************
c******************************************************************
      function absco13(j,fr0,f,p,t,wtmol)     
c
c     this function calculates the absorption coeff of CO
C  f, fr0 in Hz; p in bar, but the routine expects:
C f, fr0 in Hz, and p in microbar
c     normalized by the CO mixing ratio (xco)
c     the formulae were taken from Clancy's PhD thesis
C
      if(j.eq.0) then
          absco13 = 0.0
          return
      endif
      p=p*1.0E6
c
      pi = 3.14159265308
c
c     offset frequency from line center
      fpm = f - fr0
      fpp = f + fr0
c
c     pressure broadened line width; take same parameters as for CO12
      wp = 2100.0 * p * (300.0/t)**0.7

c     doppler width (e-folding ... not FWHM)  in Hz
      wd = fr0 / 3.00e10 * sqrt(2.0*1.38e-16*t/wtmol)
c
c     number density of atmosphere (ideal gas p = dkt)
      d = p / 1.38e-16 / t
c
      a = wp/wd
      x = fpm/wd
c
c     use the Voigt function line shape
      shape = voigt(x,a) / wd
c       
c     VanVleck-Weisskopf shape (for future use at large offsets)
      vvwshape = 
     & 1.0 / pi * f*f/fr0/fr0 * (wp/(fpm*fpm+wp*wp)+wp/(fpp*fpp+wp*wp))
c
c     Lorentz Shape
      vl = wp / pi / (fpm*fpm + wp*wp)
c
        shape=shape*1.0E6
        v1=v1*1.0E6
        vvwshape=vvwshape*1.0E6
        fshape=shape
        if(v1.gt.fshape) fshape=v1
C  The following is correct, I think -- 15 April 2011
        if(vvwshape.gt.fshape) fshape=vvwshape
c  I think pin = p/1.013E6 so P in is in atm; it was written as *1.013E6
        pin = p / 1.01325E6
        xshape=fshape
C         xshape=shape
         call absx_co13line(j,aon,t,pin,f,xshape)
C  xshape: p in millibar, f in MHz (p/f)
C vvwshape: p in microbar; f in Hz
C We multiply shape by 1E6 to correct answer to come out in bars and MHz
       absco13 = aon
c       write(6,1003) j,aon,absco13,shape,f,x,a,wp,wd
c       write(6,1004) p,xshape,vvwshape,v1,shape
c 1003  format(I5,8E12.4)
c 1004  format(5E12.4)
      p=p/1.0E6
C  if other gases are considered, add them here as well I guess; e.g. 
C  SO2, H2SO4 vapor. VVW line shape is probably OK for lower layers
c
      return
      end
c******************************************************************
c******************************************************************
      function voigt(xx,a)
C  Need to know: T, P, But it is all passed along in xx, a
c
c     voigt function subroutine obtained from D. van Blerkom
c
c     This is not the ONE TRUE Voigt funtion ...
c     The normal definition of voigt function is:
c            H(a,b) = a/pi * INTEGRAL ( exp(-y*y) / ((b-y)**2 + a**2) dy )
c            where a = wp / wd
c                  b = (f - fr0) / wd 
c     We compute:
c            voigt = 1/sqrt(pi) * H(a,b)
c
c     Function requires further normalization by doppler width.
c
      dimension c(31)
      data ic/0/,epsiln/2.5e-12/
c
c
C       write(6,1003) xx,a
C 1003  format(2E12.4)
      if(ic.ne.0) goto 20
   10 k=-16
      pie  = 3.1415926535898
      sqp  = sqrt(pie)
      sqpp = 1./sqp
      pie3 = 3.0*pie
      pie32= pie3*pie
      piec = 0.5/(pie*sqp)
      do 15 n=1,31
         k=k+1
   15    c(n)=piec*exp(-float(k*k)/9.)
      ic=1
c
   20 x=abs(xx)
c
      if(a.gt.1000.0)  go to 40
      if(a.lt.0.001)  go to 45
c
      a1=3.0*a
      a2=a*a
      if(a.lt.0.1) go to 25
      x1=-pie3*a
      voigt=0.0
      rz=0.0
      b2=0.0
      if(x1.lt.-675.0) go to 30
      y1=pie3*x
      e1=exp(x1)
      rz=e1*cos(y1)
      b2=-e1*sin(y1)
      go to 30
   25 x1=pie3*x
      y1=pie3*a
      e1=exp(y1)
      e2=1.0/e1
      rz=0.5*(e1+e2)*cos(x1)
      b2=-0.5*(e2-e1)*sin(x1)
      qin=a2-x*x
      voigt=0.0
      if(qin.gt.-675.0) voigt=sqpp*exp(qin)*cos(2.0*a*x)
   30 b1=(1.0-rz)*a*1.5
      s=-8.0-1.5*x
      t=s*s+2.25*a2
      c2=-a*pie32
      do 35 n=1,31
         t=t+s+0.25
         s=s+0.5
         b1=a1-b1
         b2=-b2
         con=c2
         if(t.gt.epsiln) con=(b1+b2*s)/t
   35    voigt=voigt+con*c(n)
      return
c
c     large a limit
c
   40 qin=x*x/a+a
      voigt=1.0/(pie*qin)
      return
c
c     small a limit
c
   45 voigt=0.0
      x2=x*x
      if(x2.lt.675.0) voigt=sqpp*exp(-x2)
      return
      end
c******************************************************************
              SUBROUTINE ABSX_CO13(AON,T,P,vnu)
C
c  Note: This subroutine is for 13CO.
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF CO13, IN UNITS OF 1/KM
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/KM
C	T	TEMPERATURE				      K
C    	P	PRESSUME				MILLIBARS
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C  the Q0, F0, S0, EL data are from Poynter and Picket, JPL catalogue
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V,wtmol
	REAL*4 QLOG2,SHAPE2,SHAPE
	REAL*4 C,ZCO,DV,VNU
        EXTERNAL QLOG2,SHAPE2,SHAPE
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=9)   !NUMBER OF SPECTRAL LINES
        REAL*4 DF,DFC,DFD,EL(NLINES),Q0(5),F0(NLINES),FS,FV
        REAL*4 S0(NLINES),S1(NLINES),TCOR,X,X1,Y
        DATA Q0 /2.0564, 1.9317, 1.7565, 1.4579, 1.1620/ !LOG PF AT 300,225,150,75,37.5K
         DATA F0(1),S0(1),EL(1)
     &  /110201.3672, -5.0614, 0.0000/
        DATA F0(2),S0(2),EL(2)  
     &  /220398.7031, -4.1698, 3.6759/
        DATA F0(3),S0(3),EL(3)  
     &  /330588.0937, -3.6607, 11.0276/
        DATA F0(4),S0(4),EL(4)
     &  /440765.4375, -3.3127, 22.0549/
        DATA F0(5),S0(5),EL(5)
     &  /550926.7500, -3.0564, 36.7572/
        DATA F0(6),S0(6),EL(6)
     &  /661068.0625, -2.8609, 55.1341/
        DATA F0(7),S0(7),EL(7)
     &  /771185.3125, -2.7098, 77.1850/
        DATA F0(8),S0(8),EL(8)
     &  /881274.5625, -2.5933, 102.9090/
        DATA F0(9),S0(9),EL(9)
     &  /991331.8125, -2.5049, 132.3051/
	P=P*1000.0/1.013		!MILLIBAR
	C=2.99791E1
	F=vnu*C*1.0E3                  !FREQU IN MHZ
	AON=0.0
       DO N = 1,NLINES
       TCOR = EXP(1.4388*EL(N)*(3.33356E-03-1./T))
     & *(1.-EXP(-4.7993E-05*F0(N)/T))
     & /(1.-EXP(-1.5998E-07*F0(N)))
       AON=AON+TCOR*SHAPE(P,T,F,F0(N))*(10.**S0(N))
       enddo
       AON=AON*10**(Q0(1)-QLOG2(Q0,T)) !T DEPENDENCE OF PARTITION FUNCTION
       AON=0.72435E5*AON*P/T       !YIELD ABSORPTION COEFFICIENT IN 1/cm
C                                    WHEN MULTIPLIED BY MIXING RATIO
C
	P=P*1.013/1000.0
       RETURN
       END
C
c******************************************************************
              SUBROUTINE ABSX_CO13line(n,AON,T,P,f,xshape)
C
c  Note: This subroutine is for 13CO.
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF CO13, IN UNITS OF 1/KM
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/KM
C	T	TEMPERATURE				      K
C    	P	PRESSUME				MILLIBARS
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C  the Q0, F0, S0, EL data are from Poynter and Picket, JPL catalogue
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V,wtmol
	REAL*4 QLOG2,SHAPE2,xshape
	REAL*4 C,ZCO,DV,VNU
        EXTERNAL QLOG2,SHAPE2
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=9)   !NUMBER OF SPECTRAL LINES
        REAL*4 DF,DFC,DFD,EL(NLINES),Q0(5),F0(NLINES),FS,FV
        REAL*4 S0(NLINES),S1(NLINES),TCOR,X,X1,Y
        DATA Q0 /2.0564, 1.9318, 1.7565, 1.4579, 1.1620/ !LOG PF AT 300,225,150,75, 37.5K
         DATA F0(1),S0(1),EL(1)
     &  /110201.3672, -5.0614, 0.0000/
        DATA F0(2),S0(2),EL(2)  
     &  /220398.7031, -4.1698, 3.6759/
        DATA F0(3),S0(3),EL(3)  
     &  /330588.0937, -3.6607, 11.0276/
        DATA F0(4),S0(4),EL(4)
     &  /440765.4375, -3.3127, 22.0549/
        DATA F0(5),S0(5),EL(5)
     &  /550926.7500, -3.0564, 36.7572/
        DATA F0(6),S0(6),EL(6)
     &  /661068.0625, -2.8609, 55.1341/
        DATA F0(7),S0(7),EL(7)
     &  /771185.3125, -2.7098, 77.1850/
        DATA F0(8),S0(8),EL(8)
     &  /881274.5625, -2.5933, 102.9090/
        DATA F0(9),S0(9),EL(9)
     &  /991331.8125, -2.5049, 132.3051/
	P=P*1000.0/1.013		!MILLIBAR
	C=2.99791E1
        F=F/1.0E6
C	F=vnu*C*1.0E3                  !FREQU IN MHZ
	AON=0.0
C        xshape=SHAPE2(P,T,F,F0(N))
       TCOR = EXP(1.4388*EL(N)*(3.33356E-03-1./T))
     & *(1.-EXP(-4.7993E-05*F0(N)/T))
     & /(1.-EXP(-1.5998E-07*F0(N)))
       AON=AON+TCOR*xshape*(10.**S0(N))
       AON=AON*10**(Q0(1)-QLOG2(Q0,T)) !T DEPENDENCE OF PARTITION FUNCTION
       AON=0.72435E5*AON*P/T       !YIELD ABSORPTION COEFFICIENT IN 1/cm
C                                    WHEN MULTIPLIED BY MIXING RATIO
C
	P=P*1.013/1000.0
        F=F*1.0E6
       RETURN
       END
C
C********************
C
       SUBROUTINE ABSX_CO(AON,T,P,vnu)
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF CO, IN UNITS OF 1/cm
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/cm
C	T	TEMPERATURE				      K
C    	P	PRESSUME				atm
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C  the Q0, F0, S0, EL data are from Poynter and Picket, JPL catalogue
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V,wtmol
	REAL*4 QLOG2,SHAPE2,SHAPE
	REAL*4 C,ZCO,DV,VNU
        EXTERNAL QLOG2,SHAPE2,SHAPE
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=8)   !NUMBER OF SPECTRAL LINES
        REAL*4 DF,DFC,DFD,EL(NLINES),Q0(5),F0(NLINES),FS,FV
        REAL*4 S0(NLINES),S1(NLINES),TCOR,X,X1,Y
        DATA Q0 /2.0369, 1.9123, 1.7370, 1.4386, 1.1429/ !LOG PF AT 300,225,150,75,37.5K
        DATA F0(1),S0(1),EL(1)
     &  /115271.2018, -5.0105, 0.0000/
        DATA F0(2),S0(2),EL(2)  
     &  /230538.0000, -4.1197, 3.8450/
        DATA F0(3),S0(3),EL(3)  
     &  /345795.9899, -3.6118, 11.5350/
        DATA F0(4),S0(4),EL(4)
     &  /461040.7682, -3.2657, 23.0695/ 
        DATA F0(5),S0(5),EL(5)
     &  /576267.9305, -3.0118, 38.4481/
        DATA F0(6),S0(6),EL(6)
     &  /691473.0763, -2.8193, 57.6704/
        DATA F0(7),S0(7),EL(7)
     &  /806651.8060, -2.6715, 80.7354/
        DATA F0(8),S0(8),EL(8)
     &  /921799.7000, -2.5590, 107.6424/
	P=P*1000.0/1.013		!MILLIBAR 
	C=2.99791E1
	F=VNU*C*1.0E3                  !FREQU IN MHZ
        AON=0.0
       DO N = 1,NLINES
       TCOR = EXP(1.4384*EL(N)*(3.3784E-03-1./T))
     & *(1.-EXP(-4.7993E-05*F0(N)/T))
     & /(1.-EXP(-1.5998E-07*F0(N)))
       AON=AON+TCOR*SHAPE(P,T,F,F0(N))*(10.**S0(N))
       enddo
       AON=AON*10**(Q0(1)-QLOG2(Q0,T)) !T DEPENDENCE OF PARTITION FUNCTION
       AON=0.72435E5*AON*P/T       !YIELD ABSORPTION COEFFICIENT IN 1/cm
C                                   WHEN MULTIPLIED BY MIXING RATIO
       P=P*1.013/1000.0
C       write(6,1003) n,aon,xz,xz2
 1003  format(I5,3E12.4)
C
       RETURN
       END
C
C********************
C
       SUBROUTINE ABSX_COline(n,AON,T,P,f,xshape)
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF CO, IN UNITS OF 1/cm
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/cm
C	T	TEMPERATURE				      K
C    	P	PRESSUME				atm
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C  the Q0, F0, S0, EL data are from Poynter and Picket, JPL catalogue
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V,wtmol
	REAL*4 QLOG2,SHAPE2,SHAPE,xshape
	REAL*4 C,ZCO,DV,VNU
        EXTERNAL QLOG2,SHAPE2,SHAPE
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=8)   !NUMBER OF SPECTRAL LINES
        REAL*4 DF,DFC,DFD,EL(NLINES),Q0(5),F0(NLINES),FS,FV
        REAL*4 S0(NLINES),S1(NLINES),TCOR,X,X1,Y
        DATA Q0 /2.0369, 1.9123, 1.7370, 1.4386, 1.1429/ !LOG PF AT 300,225,150,75, 37.5K
        DATA F0(1),S0(1),EL(1)
     &  /115271.2018, -5.0105, 0.0000/
        DATA F0(2),S0(2),EL(2)  
     &  /230538.0000, -4.1197, 3.8450/
        DATA F0(3),S0(3),EL(3)  
     &  /345795.9899, -3.6118, 11.5350/
        DATA F0(4),S0(4),EL(4)
     &  /461040.7682, -3.2657, 23.0695/ 
        DATA F0(5),S0(5),EL(5)
     &  /576267.9305, -3.0118, 38.4481/
        DATA F0(6),S0(6),EL(6)
     &  /691473.0763, -2.8193, 57.6704/
        DATA F0(7),S0(7),EL(7)
     &  /806651.8060, -2.6715, 80.7354/
        DATA F0(8),S0(8),EL(8)
     &  /921799.7000, -2.5590, 107.6424/
	P=P*1000.0/1.013		!MILLIBAR 
	C=2.99791E1
	F=f/1.0E6				!FREQ IN MHz
C	F=VNU*C*1.0E3                  !FREQU IN MHZ
        AON=0.0
C        xshape=SHAPE(P,T,F,F0(N))
       TCOR = EXP(1.4384*EL(N)*(3.3784E-03-1./T))
     & *(1.-EXP(-4.7993E-05*F0(N)/T))
     & /(1.-EXP(-1.5998E-07*F0(N)))
       AON=AON+TCOR*xshape*(10.**S0(N))
       AON=AON*10**(Q0(1)-QLOG2(Q0,T)) !T DEPENDENCE OF PARTITION FUNCTION
       AON=0.72435E5*AON*P/T       !YIELD ABSORPTION COEFFICIENT IN 1/cm
C                                   WHEN MULTIPLIED BY MIXING RATIO
       P=P*1.013/1000.0
C       write(6,1003) n,aon,xz,xz2
C 1003  format(I5,3E12.4)
       f=f*1.0E6
C
       RETURN
       END
C
       SUBROUTINE ABSX_HCN(AON,T,P,vnu)
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF Hcn, THIOFORMALDEHYDE, IN UNITS OF 1/KM
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/CM
C	T	TEMPERATURE				      K
C    	P	PRESSUME				MILLIBARS
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(P,T,F,FLC) 	LINE SHAPE FACTOR
C
C  the Q0, F0, S0, EL data are from Poynter and Picket, JPL catalogue
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V
	REAL*4 QLOG2,SHAPE
	REAL*4 C,ZHCN,DV,VNU
        EXTERNAL QLOG2,SHAPE
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=11)   !NUMBER OF SPECTRAL LINES
        REAL*4 DF,DFC,DFD,EL(NLINES),Q0(5),F0(NLINES),FS,FV
        REAL*4 S0(NLINES),S1(NLINES),TCOR,X,X1,Y
        DATA Q0 /2.6277,2.5031,2.3276,2.0286,1.7317/ !LOG PF AT 300,225,150,75,37.5K
        DATA F0(1),S0(1),EL(1)
     & /177262.0156,  -3.8432,   2.9564/
        DATA F0(2),S0(2),EL(2)
     & /88630.4141,  -2.9588,     0.000/
        DATA F0(3),S0(3),EL(3)
     & /88631.8437,  -2.7370,     0.000/
        DATA F0(4),S0(4),EL(4)
     & /88633.9375,  -3.4359,     0.000/
        DATA F0(5),S0(5),EL(5)
     & /177259.6719,  -2.6671,     2.9564/
        DATA F0(6),S0(6),EL(6)
     & /177259.9219,  -2.5422,     2.9564/
        DATA F0(7),S0(7),EL(7)
     & /177261.1094,  -2.1900,     2.9564/
        DATA F0(8),S0(8),EL(8)
     & /177261.2188,  -1.9189,     2.9564/
        DATA F0(9),S0(9),EL(9)
     & /177263.4375,  -2.6671,     2.9564/
        DATA F0(10),S0(10),EL(10)
     & /265886.1875,  -1.0751,     8.8692/
        DATA F0(11),S0(11),EL(11)
     & /354505.4687,  -0.7220,  17.7382/
 	P=P*1000.0/1.013		!MILLIBAR
	C=2.99791E1
	F=VNU*C*1.0E3                  !FREQU IN MHZ
	AON=0.0
       DO N = 1,NLINES
       TCOR = EXP(1.4388*EL(N)*(3.33356E-03-1./T))
     & *(1.-EXP(-4.7993E-05*F0(N)/T))
     & /(1.-EXP(-1.5998E-07*F0(N)))
       AON=AON + TCOR*SHAPE(P,T,F,F0(N))*(10.**S0(N))
       ENDDO
       AON=AON*10**(Q0(1)-QLOG2(Q0,T)) !T DEPENDENCE OF PARTITION FUNCTION
       AON=0.72435E5*AON*P/T       !YIELD ABSORPTION COEFFICIENT IN 1/cm
C                                    WHEN MULTIPLIED BY MIXING RATIO
C
	P=P*1.013/1000.0
       RETURN
       END
C*************
C
       FUNCTION SHAPE(P,T,F,FLC)
	implicit none
	REAL*4 SHAPE
	real*4 DL,P,T,F,FLC,DFC,DL2
C  HPFW for collisional broadening from Clancy's thesis, for CO
C P in mbar, F and DFC in MHz when called from absx_co routines
C CO 2-1: 2.1*P*((300.0/T)**0.7)
C       DFC=188.57*P/(T**0.75)
       DFC=2.1*P*((300.0/T)**0.7)
       DL=F-FLC
       DL2=F+FLC
       SHAPE=((F*F)/(FLC*FLC))*(DFC/(DL*DL+DFC*DFC)
     + +DFC/(DL2*DL2+DFC*DFC))
	SHAPE=SHAPE/3.1415926536		!DIVIDE BY PI
       RETURN
       END
C
       SUBROUTINE ABSBR_PH3(AON,T,P,VNU,zpph3)
C
C THIS SUBROUTINE RETURNS QUANTITY AON, WHICH WHEN MULTIPLIED
C BY THE VOLUME MIXING RATIO, GIVES THE ABSORPTION COEFFICIENT
C OF PH3, IN UNITS OF 1/cm
C
C DEFINITIONS
C 	AON	ABSORPTION COEFFICIENT/VOLUME MIXING RATIO   1/cm
C	T	TEMPERATURE				      K
C    	P	PRESSUME				BARS
C  	F	FREQUENCY				MHZ
C
C ADDITIONAL FUNCTIONS REQUIRED
C	QLOG	BASE TEN LOGARITHM OF THE RATIO OF PARTITION
C		FUNCTION AT TEMPERATURE T TO PARTITION FUNCTION
C		AT T=300K
C	SHAPE(F,P) 	LINE SHAPE FACTOR
C
C*******************************************************************
        IMPLICIT NONE
        REAL*4 AON,T,P,F,V
	real*4 t3,gm,vline,ph2,phe,ph2s,pnh3,pch4,ph2o,pph3,
     + deltot,PCO,PCO13,PHCN
        real*4 DL,FLC,DL2,zpph3
	REAL*4 C,DV,VNU,coef,hck
	real*8 x,ajk,fjk,tcor
        INTEGER*2 N,NLINES
        PARAMETER (NLINES=320)   !NUMBER OF SPECTRAL LINES
	real*4  FP0(728),SP0(728),ELP(728),Q0(5)
	REAL*4 wgtS0(40),wgtFGB(40),wgtSB(40)
	REAL*4 QLOG2,SHAPE2
        EXTERNAL QLOG2
       DATA Q0 /2.905, 2.722, 2.4606, 2.0132, 1.5705/ !LOG PF AT 300,225,150,75, 37.5K
        COMMON/ph3/FP0,SP0,ELP,wgtS0,wgtFGB,wgtSB
        COMMON/U1/PH2,PHE,PNH3,PH2O,PCH4,PH2S,PPH3,PCO,PCO13,PHCN 
	C=2.997925E1
	coef=7.244E21
	hck=1.438396	!hc/k
	V=VNU*C				!FREQ IN GHZ
	T3=(300/T)**0.6667
	AON=0.0
        DO N = 1,NLINES
	VLINE=FP0(N)		!LINE CENTER IN GHZ
        TCOR = EXP(-(hck*ELP(N)*(1./T - 3.3333E-03)))
	call sajph3(N,T,T3,GM)
 	if (N.gt.40) goto 40
	ajk=SP0(N)*wgtS0(N)*TCOR
	goto 50
  40    ajk=SP0(N)*TCOR
  50   continue
        CALL SFJph3(N,T,V,VLINE,GM,FJK)
       AON=AON + ajk*FJK
       ENDDO
	x=(300.0/T)**3.5
	AON=2.4163E3*x*AON*PPh3/300.0		!2.4163E3 is from deboer
C
       RETURN
       END

C
C                                                                               
