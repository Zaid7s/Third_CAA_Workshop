PROGRAM MAIN

  USE Precision_Def
  USE FFT1D

  IMPLICIT NONE

  REAL(KIND = r_def), DIMENSION(16):: xCos, xSin
  REAL(KIND = r_def):: pi, omega, ds, s, t, omega2
  INTEGER(KIND = i_def):: m, n

  m = 16_i_def

  pi = 4.0_r_def*ATAN(1.0_r_def)

  omega  = 4.0_r_def*pi
  omega2 = 8.0_r_def*pi

  ds = 1.0/REAL(m, i_def)

  WRITE(6, *) 
  WRITE(6, *) 'index   FFT Cos input    FFT Sin Input '
  WRITE(6, *) '-----   ------------     ------------- '
  WRITE(6, *) 

  DO n = 1, 16
   s = REAL(n-1, r_def)*ds 
   t = REAL(n, r_def)/ds
   xCos(n) = COS(omega*s) + COS(omega2*s)
   xSin(n) = SIN(omega*s) + SIN(omega2*s)
   WRITE(6, *) n, xCos(n), xSin(n)
   WRITE(16, *) n, xCos(n), xSin(n)
  END DO

  CALL real_fft(data  = xCos, &
                n     = m,    &
                isign = 1_i_def)
  CALL real_fft(data  = xSin, &
                n     = m,    &
                isign = 1_i_def)

  WRITE(6, *) 
  WRITE(6, *) 'index   FFT Cos output   FFT Sin Output '
  WRITE(6, *) '-----   ------------     -------------- '
  WRITE(6, *) 

  DO n = 1, 16
   WRITE(6, *) n, xCos(n), xSin(n)
   WRITE(16, *) n, xCos(n), xSin(n)
  END DO

  STOP
END PROGRAM MAIN
