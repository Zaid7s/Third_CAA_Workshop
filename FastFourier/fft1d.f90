MODULE FFT1D

  USE Precision_Def  ! Precision_Def.f90
  IMPLICIT NONE
  PRIVATE
  PUBLIC:: real_fft, one_d_fft

CONTAINS

SUBROUTINE real_fft(data, n, isign)

  REAL(KIND = r_def), INTENT(INOUT), DIMENSION(:):: data
  INTEGER(KIND = i_def), INTENT(IN):: n, isign

! define local variables

  INTEGER(KIND = i_def):: i, i1, i2, i3, i4, n2p3
  REAL(KIND = r_def):: c1, c2, h1i, h1r, h2i, h2r, wis, wrs
  REAL(KIND = 8):: theta, wi, wpi, wpr, wr, wtemp, pi

  CONTINUE  ! computation starts here

  pi = 4.0_8*(ATAN(1.0_8))

  theta = pi/REAL(n/2_i_def, r_def)
  c1 = 0.5_r_def

  IF (isign == 1) THEN
    c2 = -0.5_r_def
    CALL one_d_fft(data  = data,       &
                   nn    = n/2_i_def,  &
                   isign = 1_i_def)
  ELSE
    c2 = 0.5_r_def
    theta = -theta
  END IF

  wpr = -2.0_8*SIN(0.5_8*theta)**2
  wpi = SIN(theta)
  wr  = 1.0_8+wpr
  wi  = wpi
  n2p3 = n+3_i_def

  DO i = 2, n/4_i_def
   i1 = 2_i_def*i-1_i_def
   i2 = i1+1_i_def
   i3 = n2p3-i2
   i4 = i3+1_i_def
   wrs = REAL(wr, r_def)
   wis = REAL(wi, r_def)
   h1r =  c1*(data(i1)+data(i3))
   h1i =  c1*(data(i2)-data(i4))
   h2r = -c2*(data(i2)+data(i4))
   h2i =  c2*(data(i1)-data(i3))
   data(i1) =  h1r+wrs*h2r-wis*h2i
   data(i2) =  h1i+wrs*h2i+wis*h2r
   data(i3) =  h1r-wrs*h2r+wis*h2i
   data(i4) = -h1i+wrs*h2i+wis*h2r
   wtemp = wr
   wr = wr*wpr-wi*wpi+wr
   wi = wi*wpr+wtemp*wpi+wi
  END DO

  IF (isign == 1_i_def) THEN
   h1r = data(1)
   data(1) = h1r+data(2)
   data(2) = h1r-data(2)
  ELSE
   h1r = data(1)
   data(1) = c1*(h1r+data(2))
   data(2) = c1*(h1r-data(2))
   CALL one_d_fft(data  = data,       &
                  nn    = n/2_i_def,  &
                  isign = -1_i_def)
  END IF

  RETURN
END SUBROUTINE real_fft

SUBROUTINE one_d_fft(data, nn, isign)
 
  INTEGER(KIND = i_def), INTENT(IN):: nn, isign
  REAL(KIND = r_def), DIMENSION(:), INTENT(INOUT):: data

! Replaces data(1:2*nn) by its discrete Fourier transform, if isign
!  is input as 1; or replaces data(1:2*nn) by nn times its inverse
!  discrete Fourier transform, if isign is input as-1.
!
! data is a real array of length 2*nn (F77 version), but actually contains
!  a complex array of length nn (stacked in as real, imag, real, imag...).  
!
! nn MUST be an integer power of 2.
!

! local variables

  INTEGER(KIND = i_def):: i, istep, j, m, mmax, n
  REAL(KIND = r_def):: tempi, tempr, nreal
  REAL(KIND = r_def):: theta, wi, wpi, wpr, wr, wtemp, two_pi

  two_pi = 8.0_r_def*ATAN(1.0_r_def)
  n = 2*nn
  j = 1

  DO i = 1, n, 2
   IF (j > i) THEN
    tempr     = data(j)  
    tempi     = data(j+1)  
    data(j)   = data(i)
    data(j+1) = data(i+1)
    data(i)   = tempr
    data(i+1) = tempi
   END IF
   m = nn

 1 CONTINUE
   
   IF ((m >= 2) .AND. (j > m)) THEN
    j = j-m
    m = m/2
    GO TO 1
   END IF
   j = j+m
  END DO
  mmax = 2

 2 CONTINUE

  IF (n > mmax) THEN
   istep = 2*mmax
   theta = two_pi/(isign*mmax)
   wpr = -2.0_r_def*SIN(0.5_r_def*theta)**2
   wpi = SIN(theta)
   wr = 1.0_r_def
   wi = 0.0_r_def
   DO m = 1, mmax, 2
    DO i = m, n, istep
     j = i+mmax
     tempr = REAL(wr, 4)*data(j) - REAL(wi, 4)*data(j+1)
     tempi = REAL(wr, 4)*data(j+1) + REAL(wi, 4)*data(j)
     data(j)   = data(i)   - tempr
     data(j+1) = data(i+1) - tempi
     data(i)   = data(i)   + tempr
     data(i+1) = data(i+1) + tempi
    END DO
    wtemp = wr
    wr = wr*wpr-wi*wpi+wr
    wi = wi*wpr+wtemp*wpi+wi
   END DO
   mmax = istep
   GO TO 2
  END IF

! normalize

  IF (isign == 1) THEN
   nreal = 1.0/REAL(nn)
   DO i = 1, 2*nn
    data(i) = data(i)*nreal
   END DO
  END IF

  RETURN
END SUBROUTINE one_d_fft 

END MODULE FFT1D

