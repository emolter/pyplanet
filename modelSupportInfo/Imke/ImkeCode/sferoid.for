	PROGRAM sferoid
c  program is set up for Saturn; if other planet is desired, change
c  zrho parameters, and gridsize.
c  alpha: negative if ring plane is tilted up; positive if down.
C f77 -O -o sf sferoid.for
C ifort -132 sferoid.for -o sf.exe
C ./sf.exe       to run it
C  man ifort to look for compiling options, like -132
	dimension zfi(1800),zrho(1800),zy(1800),zyprime(1800,13)
	dimension zx(13),ZFIP(13,25),zyp(25),ZYY(1800,13)
	REAL*8 ALPHA,OR,OP,CS,THETA
	OPEN(6,NAME='SPr.DAT')
	OPEN(1,NAME='zmu26.dat',FORM='UNFORMATTED',STATUS='NEW',recl=25)
c increment in x and y is .0889
c  zfi is fiprime
	zfi(1)=-90.0
	zx(1)=0.0
	zyp(1)=-12.0*0.0889
c inclination angle alpha
	alpha=26.0
	ABSA=ABS(ALPHA)
	pi=3.14159265
	pi180=pi/180.0
	alpha=alpha*pi180
	zfi(1)=zfi(1)*pi180
	DO I=2,1800
	ZFI(I)=ZFI(I-1)+PI180/10.0
	ENDDO
	DO I=2,13
	ZX(I)=(I-1)*0.0889
	ENDDO
	do i=2,25
	zyp(i)=(I-1)*0.0889+zyp(1)
	enddo
	DO I=1,1800
	THETA=PI/2.0-ZFI(I)
	CS=COS(THETA)*COS(THETA)
	ZRHO(I)=-0.0075*CS*CS*CS+0.0146*CS*CS-0.0135*CS+1.1075
	ZY(I)=ZRHO(I)*SIN(ZFI(I))
	DO J=1,13
	zsqu=zx(j)*zx(j)+zy(i)*zy(i)
	ORSQU=ZRHO(I)*ZRHO(I)-ZX(J)*ZX(J)
	OR=SQRT(ORSQU)
	OP=ZRHO(I)*ZRHO(I)-zsqu
	if(op.ge.0.0) OP=SQRT(OP)
	if(OP.lt.0.0) OP=0.0
	if(zy(i).ne.0.0) x=(op/zy(i))
	IF(ZY(I).NE.0.0) ORP=atan(x)
	IF(ZY(I).EQ.0.0) ORP=PI/2.0
	orp=orp+ALPHA
	ZYPRIME(I,J)=OR*COS(ORP)
	IF(ZFI(i).lt.ALPHA.AND.ZYPRIME(I,J).GT.0.0) 
     + ZYPRIME(I,J)=-ZYPRIME(I,J)
	IF(zFI(i).GE.ALPHA.AND.ZYPRIME(I,J).LT.0.0) 
     + ZYPRIME(I,J)=-ZYPRIME(I,J)
	ENDDO
	ENDDO
	DO I=1,1800
	ZFI(I)=ZFI(I)/PI180
	ENDDO
	do IY=1,25
	DO IX=1,13
	do i=1,1799
	IF(ZYPRIME(I,IX).GE.Zyp(Iy).AND.ZYPRIME(I,Ix).LT.Zyp(Iy+1)) GOTO 10
	ENDDO
  10    ZFIP(IX,IY)=ZFI(I)
	ENDDO
	ENDDO
C	WRITE(6,205) (Zx(K),K=1,13)
C	do i=1,180
C	write(6,20) zfi(i),zy(i),zrho(i)
C	write(6,30) (zyprime(i,k),k=1,13)
C	enddo
 30     format(13(F6.3))
 20     format(3e12.5)
	WRITE(6,200) (Zyp(K),K=1,12)
	DO K=1,13
	WRITE(6,100) ZX(K),(ZFIP(K,J),J=1,12)
	ENDDO
	WRITE(6,200) (Zyp(K),K=13,24)
	DO K=1,13
	WRITE(6,105) ZX(K),(ZFIP(K,J),J=13,25)
	ENDDO
	do k=1,13
	do j=1,25
	zfip(k,j)=zfip(k,j)*pi180
	enddo
	enddo
	DO K=1,13
	WRITE(1) (ZFIP(K,J),J=1,25)
	ENDDO
 100    FORMAT(F6.2,1x,12(F5.1,1X))
 105    FORMAT(F6.2,1x,13(F5.1,1X))
 110    FORMAT(1x,12(F5.1,1X))
 120    FORMAT(1x,13(F5.1,1X))
 200    FORMAT(5X,12(F6.2))
 205    FORMAT(5X,13(F6.2))
	STOP
	END
