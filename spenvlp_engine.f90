! --------------------------------------------------
SUBROUTINE spenvlp_engine (npass, maxit, nobs, lambda, weight, eps, gamma, a_vec, b2, &
	& b3, A1, A2, A3)
! --------------------------------------------------
	IMPLICIT NONE
! - - - arg declarations - - -
	INTEGER :: npass
	INTEGER :: maxit
	INTEGER :: nobs
	DOUBLE PRECISION :: lambda
	DOUBLE PRECISION :: weight
	DOUBLE PRECISION :: eps
	DOUBLE PRECISION :: gamma
	DOUBLE PRECISION :: a_vec (nobs)
	DOUBLE PRECISION :: b2 (nobs)
	DOUBLE PRECISION :: b3 (nobs)
	DOUBLE PRECISION :: A1 (nobs, nobs)
	DOUBLE PRECISION :: A2 (nobs, nobs)
	DOUBLE PRECISION :: A3 (nobs, nobs)
! - - - local declarations - - -
	DOUBLE PRECISION :: U_norm
	DOUBLE PRECISION :: t
	DOUBLE PRECISION :: olda_vec (nobs)
	DOUBLE PRECISION :: tmp1 (nobs)
	DOUBLE PRECISION :: tmp2 (nobs)
	DOUBLE PRECISION :: tmp3 (nobs)
	DOUBLE PRECISION :: U (nobs)
	DOUBLE PRECISION :: U_working (nobs)
	DOUBLE PRECISION :: dif (nobs)
! - - - program start - - -
	olda_vec = a_vec
	npass = 0		
	DO
		olda_vec = a_vec
		tmp1=MATMUL(A1, olda_vec)
		tmp2=MATMUL(A2, olda_vec + b2)
		tmp3=MATMUL(A3, olda_vec + b3)
		U = 4.0E0*tmp1/(1.0E0+DOT_PRODUCT(olda_vec, tmp1))-	&				
		& 2.0E0*tmp2/(1.0E0+DOT_PRODUCT(olda_vec+b2, tmp2))- &		
		& 2.0E0*tmp3/(1.0E0+DOT_PRODUCT(olda_vec+b3, tmp3))
		U_working = U + gamma * olda_vec
		U_norm = sqrt(DOT_PRODUCT(U_working,U_working))
		t = U_norm - weight * lambda
		IF(t > 0) THEN
			a_vec = U_working * t / (gamma * U_norm)
		ELSE
			a_vec = 0.0E0
		ENDIF		
	    dif = a_vec - olda_vec
	    if(gamma * DOT_PRODUCT(dif, dif) < eps) EXIT
	    npass = npass + 1
	    if(npass > maxit) EXIT
	ENDDO
END SUBROUTINE spenvlp_engine