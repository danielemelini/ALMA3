!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
!
! Adapted for ALMA 2.0 from the F90 version by J.P. Moreau, Paris 
! Revised Giorgio Spada May 2010 for g95 porting of ALMA
!
 Subroutine LUDCMP(A,N,NP,INDX,D,CODE)
 USE FMZM
 IMPLICIT NONE 
 INTEGER, PARAMETER :: NMAX=100 
 TYPE(FM) STINY, AMAX, DUM, VV(NMAX)
 TYPE(ZM) A(N,N), XSUM, CDUM
 INTEGER I, J, K, IMAX, CODE, D, N, NP, INDX(N)
!
 D=1
 CODE=0
!
 STINY=TO_FM(1.5e-16)
!
 DO I=1,N
   AMAX=TO_FM('0.0')
   DO J=1,N
     IF (ABS(A(I,J)).GT.AMAX) AMAX=ABS(A(I,J))
   END DO  
   IF(AMAX.LT.STINY) THEN
     CODE = 1
     RETURN
   END IF
   VV(I) = TO_FM('1.0')/ AMAX
 END DO  
!
 DO J=1,N
   DO I=1,J-1
     XSUM = A(I,J)
     DO K=1,I-1
       XSUM = XSUM - A(I,K)*A(K,J) 
     END DO  
     A(I,J) = XSUM
   END DO 
   AMAX = TO_FM('0.0')
   DO I=J,N
     XSUM = A(I,J)
     DO K=1,J-1
       XSUM = XSUM - A(I,K)*A(K,J) 
     END DO  
     A(I,J) = XSUM
     DUM = VV(I)*ABS(XSUM)
     IF(DUM.GE.AMAX) THEN
       IMAX = I
       AMAX = DUM
     END IF
   END DO  
!   
   IF(J.NE.IMAX) THEN
     DO K=1,N
       CDUM = A(IMAX,K)
       A(IMAX,K) = A(J,K)
       A(J,K) = CDUM
     END DO  
     D = -D
     VV(IMAX) = VV(J)
   END IF
!
   INDX(J) = IMAX
   IF(ABS(A(J,J)) < STINY) A(J,J) = STINY
!
   IF(J.NE.N) THEN
     CDUM = TO_ZM('1.0')/ A(J,J)
     DO I=J+1,N
       A(I,J) = A(I,J)*CDUM
     END DO  
   END IF 
 END DO  
!
 RETURN
 END subroutine LUDCMP
!
!
!
!
!
!
!  ******************************************************************
!  * Solves the set of N linear equations A . X = B.  Here A is     *
!  * input, not as the matrix A but rather as its LU decomposition, *
!  * determined by the routine LUDCMP. INDX is input as the permuta-*
!  * tion vector returned by LUDCMP. B is input as the right-hand   *
!  * side vector B, and returns with the solution vector X. A, N and*
!  * INDX are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. This routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
!
! Adapted for ALMA 2.0 from the F90 version by J.P. Moreau, Paris 
! Revised Giorgio Spada May 2010 for g95 porting of ALMA
!
 Subroutine LUBKSB(A,N,NP,INDX,B)
 USE FMZM
 IMPLICIT NONE 
!
 TYPE(ZM) XSUM, A(N,N), B(N)
 INTEGER N, NP, II, I, J, K, LL, INDX(N)
!
 II = 0
!
 DO I=1,N
   LL = INDX(I)
   XSUM = B(LL)
   B(LL) = B(I)
   IF(II.NE.0) THEN
     DO J=II,I-1
       XSUM = XSUM - A(I,J)*B(J)
     END DO  
   ELSE IF(XSUM.NE.0.) THEN
     II = I
   END IF
   B(I) = XSUM
 END DO  
!
 DO I=N,1,-1
   XSUM = B(I)
   IF(I < N) THEN
     DO J=I+1,N
       XSUM = XSUM - A(I,J)*B(J)
     END DO  
   END IF
   B(I) = XSUM / A(I,I)
 END DO 
!
 RETURN
 END subroutine LUBKSB
!
!
!
!
!
!
