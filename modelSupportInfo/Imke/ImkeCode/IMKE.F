C
C  Paul Romani May 1987
C 
C  Program IMKE.FOR is a reverse transcription program. It
C  takes the output file produced by CLOUD9.FOR for Imke and 
C  reverses the order of the lines. This is because CLOUD9.FOR
C  program runs from high to low pressures, while her program runs
C  low to high pressures.
C 
C   ifort -132 IMKE.F -o imke.exe
C ./imke.exe       to run it 
C First add in #lines in IMKE.F
C
C
      IMPLICIT REAL*8(A-Z)
      INTEGER I, IMAX J, JJ, JMAX
C
      DIMENSION X(11,1200)
C
      OPEN (UNIT=1, NAME='IMKEIN.DAT', STATUS='OLD')
      OPEN (UNIT=2, NAME='IMKEOUT.DAT', STATUS='NEW')
C
      IMAX = 11
C
      READ(1,*) JMAX
      DO 100 J = 1, JMAX
  100 READ(1,*) (X(I,J), I = 1, IMAX)
C
      WRITE(2,1000) JMAX
      DO 200 JJ = 1, JMAX
      J = JMAX - JJ + 1
      WRITE(2,2000) (X(I,J), I = 1, IMAX)
  200 CONTINUE
C
 1000 FORMAT(2X, I3)
 2000 FORMAT(1X, F6.1, 2X, F9.2, 2X, 1PE11.4, 2(2X, 0PF5.3),
     1   6(2X, 1PE8.2) )
C
      STOP
      END
