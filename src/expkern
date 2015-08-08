! --------------------------------------------------
SUBROUTINE expkern (omega, Kmat, Umat, Dvec, Ksum, Bmat, &
	& nobs, y, nlam, ulam, eps, &
    & maxit, anlam, npass, jerr, alpmat)
! --------------------------------------------------
      IMPLICIT NONE
    ! - - - arg types - - -
      INTEGER :: nobs
      INTEGER :: nlam
      INTEGER :: anlam
      INTEGER :: jerr
      INTEGER :: maxit
      INTEGER :: npass (nlam)
      DOUBLE PRECISION :: eps
      DOUBLE PRECISION :: omega
      DOUBLE PRECISION :: Kmat (nobs, nobs)
      DOUBLE PRECISION :: Umat (nobs, nobs)
      DOUBLE PRECISION :: Dvec (nobs)
      DOUBLE PRECISION :: Ksum (nobs)
      DOUBLE PRECISION :: Bmat (nobs, nobs)
      DOUBLE PRECISION :: y (nobs)
      DOUBLE PRECISION :: ulam (nlam)
      DOUBLE PRECISION :: alpmat (nobs+1, nlam)
    ! - - - local declarations - - -
      INTEGER :: j
      INTEGER :: l
      DOUBLE PRECISION :: K0 (nobs+1, nobs+1)
      DOUBLE PRECISION :: Ki (nobs+1, nobs)
	  DOUBLE PRECISION :: mbd
	  DOUBLE PRECISION :: Ginv
	  DOUBLE PRECISION :: r (nobs)
	  DOUBLE PRECISION :: phi (nobs)
	  DOUBLE PRECISION :: alpvec (nobs+1)
	  DOUBLE PRECISION :: oalpvec (nobs+1)
	  DOUBLE PRECISION :: dif (nobs+1)
	  DOUBLE PRECISION :: Ddiag (nobs, nobs)
	  DOUBLE PRECISION :: Ainv (nobs, nobs)
	  DOUBLE PRECISION :: BAmat (nobs, nobs)
	  DOUBLE PRECISION :: Qinv (nobs, nobs)
	  DOUBLE PRECISION :: QKsum (nobs)
	  DOUBLE PRECISION :: KUinv (nobs+1, nobs+1)
    
    ! - - - begin - - -
    K0 = 0.0D0
    K0(2:(nobs+1),2:(nobs+1)) = Kmat
    Ki = 1.0D0
    Ki(2:(nobs+1),:) = Kmat
	mbd = 2.0D0 * MAX(1.0D0-omega, omega)
	npass = 0
	r = y
	alpmat = 0.0D0
	alpvec = 0.0D0
	DO l = 1,nlam
		dif = 0.0D0
	! - - - computing Ku inverse - - - 
	    Ddiag = 0.0D0
	    DO j = 1, nobs
		    Ddiag(j, j) = 1/(Dvec(j)*Dvec(j) + 2.0D0*REAL(nobs)*ulam(l)*Dvec(j)/mbd)
	    ENDDO
		Ainv = MATMUL(MATMUL(Umat,Ddiag), TRANSPOSE(Umat))
! 		CALL dblepr("Ainv",-1,Ainv,nobs*nobs)
		BAmat = MATMUL(Bmat, Ainv)
! 		CALL dblepr("BAmat",-1,BAmat,nobs*nobs)
		Ginv = 1.0D0
		DO j = 1, nobs
		    Ginv = Ginv + BAmat(j, j)
	    ENDDO
		Ginv = 1 / Ginv
		Qinv = Ainv - Ginv * MATMUL(Ainv, BAmat)
		QKsum = MATMUL(Qinv, Ksum) / REAL(nobs)
		KUinv = 0.0D0
		KUinv(1, 1) = (1.0D0 + DOT_PRODUCT(QKsum, Ksum)) / REAL(nobs)
		KUinv(1, 2:(nobs+1))  = -QKsum
		KUinv(2:(nobs+1), 1)  = -QKsum
		KUinv(2:(nobs+1), 2:(nobs+1))  = Qinv
		! update alpha
		DO
			DO j = 1, nobs
                IF (r(j) <= 0.0D0) THEN
                    phi (j) = 2.0D0 * (1.0D0 - omega) * r(j)
                ELSE
                    phi (j) = 2.0D0 * omega * r(j)
                END IF
            END DO
		    oalpvec = alpvec
		    alpvec =  oalpvec + (2.0D0*REAL(nobs)/mbd) * MATMUL(KUinv, & 
	     	          & (-ulam(l)*MATMUL(K0, oalpvec) + 0.5D0*MATMUL(Ki, phi)/REAL(nobs)))
		    dif = alpvec - oalpvec
		    r = r - MATMUL(dif, Ki)
		    IF (DOT_PRODUCT(dif,dif)/DOT_PRODUCT(oalpvec,oalpvec) < eps) EXIT
		    npass(l) = npass(l) + 1
		    IF (SUM(npass) > maxit) EXIT
		ENDDO
		alpmat(:, l) = alpvec
		IF (SUM(npass) > maxit) THEN
			jerr = -l
			EXIT
		ENDIF
		anlam = l
	ENDDO

END SUBROUTINE expkern