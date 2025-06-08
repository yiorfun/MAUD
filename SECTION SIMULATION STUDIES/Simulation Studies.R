####################################################################
### SIMULATION STUDIES     						                 ###
####################################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"mvnfast",		
	### mvnfast::rmvn(), fast generate multivariate normal samples
    "matlib", 
	### matlib::symMat(), create a symmetric matrix from a vector
	"matrixcalc",
	### matrixcalc::vech(), create a vector from a symmetric matrix
	"Matrix",
	### Matrix::bdiag(), create a block-diagonal matrix
	"mmrm",
	### mmrm::mmrm(), perform the mmrm model with various covariance
	"ggplot2",
	### ggplot2::ggplot(), create boxplots for the relative loss
	"svglite"
	### svglite::svglite, create a svg file for the plot
)

CHECK_PACKAGES <- lapply(X = REQUIRED_PACKAGES,
					   FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
    }
    require(x, character.only = TRUE)
  }
)

#############################
### Loading the functions ###
#############################

BLOCK_HADAMARD_PRODUCT <- function(A, B, p_vector){
	K <- length(p_vector)
	COL_temp <- c()
	for(k in 1 : K){
		ROW_temp <- c()
		for(kp in 1 : K){
			if(k == kp){
				ROW_temp <- cbind(ROW_temp, A[k, k] * diag(p_vector[k]) + B[k, k] * matrix(1, p_vector[k], p_vector[k]))
			} else{
				ROW_temp <- cbind(ROW_temp, B[k, kp] * matrix(1, p_vector[k], p_vector[kp]))
			}
		}
		COL_temp <- rbind(COL_temp, ROW_temp)
	}
	M <- COL_temp
	return(M)
}

SPLIT_UB_FOR_SUM <- function(MATRIX, p_vector, A, B){
	
	K <- length(p_vector)
	I_PARA <- matrix(0, nrow = K, ncol = K) ### trace(Nkk)
	J_PARA <- matrix(0, nrow = K, ncol = K) ### sum(Nkk')
	for(k in seq(K)){
		for(kp in seq(K)){
			SUB <- MATRIX[(sum(p_vector[0 : (k - 1)]) + 1) : 
						   sum(p_vector[0 : k]), 
						  (sum(p_vector[0 : (kp - 1)]) + 1) : 
						   sum(p_vector[0 : kp])]
			if(kp == k){
				#UNI <- unique(round(as.vector(SUB), 8))
				#if(length(UNI) > 2) stop("more than 2 elements in I")
				I_PARA[k, kp] <- sum(diag(SUB))
				J_PARA[k, kp] <- sum(SUB)
			} else{
				#UNI <- unique(round(as.vector(SUB), 8))
				#if(length(UNI) > 2) stop("more than 1 elements in J")
				J_PARA[k, kp] <- sum(SUB)
			}	
		}
	}
	return(sum(diag(I_PARA * A)) + sum(J_PARA * B))
}

PRODUCT_OF_FOUR_UB_MATRICES <- function(A1, A2, A3, A4, B1, B2, B3, B4, p_vector){
	
	PRODUCT_OF_TWO_UB_MATRICES <- function(A1, A2, B1, B2, p_vector){
		return(list(A = A1 %*% A2, B = A1 %*% B2 + B1 %*% A2 + B1 %*% diag(p_vector) %*% B2))
	}
	
	RES1_2 <- PRODUCT_OF_TWO_UB_MATRICES(A1 = A1, A2 = A2, B1 = B1, B2 = B2, p_vector = p_vector)
	A1_2 <- RES1_2$A
	B1_2 <- RES1_2$B
	
	RES12_3 <- PRODUCT_OF_TWO_UB_MATRICES(A1 = A1_2, A2 = A3, B1 = B1_2, B2 = B3, p_vector = p_vector)
	A12_3 <- RES12_3$A
	B12_3 <- RES12_3$B
	
	RES123_4 <- PRODUCT_OF_TWO_UB_MATRICES(A1 = A12_3, A2 = A4, B1 = B12_3, B2 = B4, p_vector = p_vector)
	A123_4 <- RES123_4$A
	B123_4 <- RES123_4$B
	
	return(list(A = A123_4, B = B123_4))
}

FISHER_INFORMATION_MATRIX_FOR_GAMMA <- function(AU, BU, AS, BS, p_vector, n, PRODUCT_OF_FOUR_UB_MATRICES){
	
	### AU = diag(- gamma_gg)
	### BU = (gamma_gg')
	
	PARTIAL_PRECISION <- function(AU, BU, p_vector, g, gp){
		
		### note that g_temp <= gp_temp
		
		Eggp <- matrix(0, length(p_vector), length(p_vector))
		Eggp[g, gp] <- 1
		R <- diag(p_vector)
		
		if(g == gp){
			PAT_A <- 2 * (1 + BU[g, gp]) * Eggp
			PAT_B <- - 2 * (1 + BU[g, gp]) * Eggp - (Eggp %*% BU + BU %*% Eggp) + Eggp %*% R %*% BU + BU %*% R %*% Eggp
		} else{
			PAT_A <- matrix(0, length(p_vector), length(p_vector))
			PAT_B <- - 2 * (Eggp + t(Eggp)) + (Eggp + t(Eggp)) %*% (AU + R %*% BU) + (AU + BU %*% R) %*% (Eggp + t(Eggp))
		}
		
		return(list(A = PAT_A, B = PAT_B))
	}

	G <- length(p_vector)
	dim_gamma <- G * (G + 1) / 2
	PAT_PRECISION_A_LIST <- list()
	PAT_PRECISION_B_LIST <- list()
	
	for(j in 1 : dim_gamma){
		
		J_VEC <- rep(0, dim_gamma)
		J_VEC[j] <- 1
		G_MAT <- matrix(0, G, G)
		G_MAT[lower.tri(G_MAT, diag = TRUE)] <- J_VEC
		### in this way, the direction is by column in the lower triangle
		### using lower.tri(), the direction is also by column in it
		G_MAT <- t(G_MAT) ### make sure the direction is by ROW
		G_MAT
		if(sum(G_MAT) == 1){
			REDS_ggp <- which(G_MAT == 1, arr.ind = TRUE)
			g_temp <- REDS_ggp[1]
			gp_temp <- REDS_ggp[2]
		} else{
			stop("non-unique location index")
		}
		### note that g_temp <= gp_temp
		
		RES_PAT_PRE_TEMP <- PARTIAL_PRECISION(AU = AU, BU = BU, p_vector = p_vector, g = g_temp, gp = gp_temp)
		PAT_PRECISION_A_LIST[[j]] <- RES_PAT_PRE_TEMP$A
		PAT_PRECISION_B_LIST[[j]] <- RES_PAT_PRE_TEMP$B
	}
	
	PHI_GAMMA_MAT <- matrix(0, dim_gamma, dim_gamma)
	
	for(j in 1 : dim_gamma){
		for(jp in j : dim_gamma){
		
			A1_TEMP <- PAT_PRECISION_A_LIST[[j]]
			A2_TEMP <- AS
			A3_TEMP <- PAT_PRECISION_A_LIST[[jp]]
			A4_TEMP <- AS
			B1_TEMP <- PAT_PRECISION_B_LIST[[j]]
			B2_TEMP <- BS
			B3_TEMP <- PAT_PRECISION_B_LIST[[jp]]
			B4_TEMP <- BS
		
			RES_PROD <- PRODUCT_OF_FOUR_UB_MATRICES(A1 = A1_TEMP, A2 = A2_TEMP, A3 = A3_TEMP, A4 = A4_TEMP, B1 = B1_TEMP, B2 = B2_TEMP, B3 = B3_TEMP, B4 = B4_TEMP, p_vector = p_vector)
			
			PHI_GAMMA_MAT[j, jp] <- sum(p_vector * (diag(RES_PROD$A) + diag(RES_PROD$B))) * (n / 2)
		}
	}
	PHI_GAMMA_MAT_DIAG <- diag(PHI_GAMMA_MAT)
	PHI_GAMMA_MAT <- PHI_GAMMA_MAT + t(PHI_GAMMA_MAT)
	diag(PHI_GAMMA_MAT) <- PHI_GAMMA_MAT_DIAG
	return(PHI_GAMMA_MAT)
}

NEG_SCORE_FUN_FOR_GAMMA <- function(gamma_vec, p_vector, S, n, SPLIT_UB_FOR_SUM, BLOCK_HADAMARD_PRODUCT){
	
	PARTIAL_PRECISION <- function(AU, BU, p_vector, g, gp){
		
		### note that g_temp <= gp_temp
		
		Eggp <- matrix(0, length(p_vector), length(p_vector))
		Eggp[g, gp] <- 1
		R <- diag(p_vector)
		
		if(g == gp){
			PAT_A <- 2 * (1 + BU[g, gp]) * Eggp
			PAT_B <- - 2 * (1 + BU[g, gp]) * Eggp - (Eggp %*% BU + BU %*% Eggp) + Eggp %*% R %*% BU + BU %*% R %*% Eggp
		} else{
			PAT_A <- matrix(0, length(p_vector), length(p_vector))
			PAT_B <- - 2 * (Eggp + t(Eggp)) + (Eggp + t(Eggp)) %*% (AU + R %*% BU) + (AU + BU %*% R) %*% (Eggp + t(Eggp))
		}
		
		return(list(A = PAT_A, B = PAT_B))
	}

	G <- length(p_vector)
	dim_gamma <- G * (G + 1) / 2
	L_MAT <- diag(p_vector)
	BU <- matlib::symMat(gamma_vec, diag = TRUE, byrow = FALSE) 
	### direction by ROW, the function symMat is incorrect
	AU <- diag(- diag(BU))
	AIU <- diag(G) - AU
	BIU <- - BU
	AO <- AIU %*% AIU
	BO <- AIU %*% BIU + BIU %*% AIU + BIU %*% L_MAT %*% BIU
	DO <- AO + BO %*% L_MAT
	AS <- solve(AO)
	BS <- - solve(DO) %*% BO %*% solve(AO)
	if (!isSymmetric(BS)) {
		BS <- (BS + t(BS)) / 2
	}
	
	GRAD_VEC <- rep(0, dim_gamma)
	
	for(j in 1 : dim_gamma){
		
		J_VEC <- rep(0, dim_gamma)
		J_VEC[j] <- 1
		G_MAT <- matrix(0, G, G)
		G_MAT[lower.tri(G_MAT, diag = TRUE)] <- J_VEC
		### in this way, the direction is by column in the lower triangle
		### using lower.tri(), the direction is also by column in it
		G_MAT <- t(G_MAT) ### make sure the direction is by ROW
		if(sum(G_MAT) == 1){
			REDS_ggp <- which(G_MAT == 1, arr.ind = TRUE)
			g_temp <- REDS_ggp[1]
			gp_temp <- REDS_ggp[2]
		} else{
			stop("non-unique location index")
		}
		### note that g_temp <= gp_temp
		
		RES_PAT_PRE_TEMP <- PARTIAL_PRECISION(AU = AU, BU = BU, p_vector = p_vector, g = g_temp, gp = gp_temp)
		PAT_A <- RES_PAT_PRE_TEMP$A
		PAT_B <- RES_PAT_PRE_TEMP$B
		
		SUM_SIGMA <- sum(diag(PAT_A * AS * L_MAT + PAT_A * BS * L_MAT + PAT_B * AS * L_MAT)) + sum(PAT_B * BS * as.matrix(p_vector %*% t(p_vector)))
		
		SUM_S <- SPLIT_UB_FOR_SUM(MATRIX = S, p_vector = p_vector, A = PAT_A, B = PAT_B)
		
		GRAD_VEC[j] <- (n / 2) * (SUM_SIGMA - SUM_S)
	}
	return(- GRAD_VEC)
}

NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA <- function(gamma_vec, p_vector, S, n, SPLIT_UB_FOR_SUM, BLOCK_HADAMARD_PRODUCT, drop_const = FALSE){
	
	BU <- matlib::symMat(gamma_vec, diag = TRUE, byrow = FALSE) 
	### direction by ROW, the function symMat is incorrect
	AU <- diag(- diag(BU))
	G <- length(p_vector)
	L_MAT <- diag(p_vector)
	
	AIU <- diag(G) - AU
	BIU <- - BU
	AO <- AIU %*% AIU
	BO <- AIU %*% BIU + BIU %*% AIU + BIU %*% L_MAT %*% BIU
	if (!isSymmetric(BO)) {
		BO <- (BO + t(BO)) / 2
	}
	DO <- AO + BO %*% L_MAT
	
	if (any(is.complex(eigen(DO)$values))){
		stop("Complex value produces.")
	}
	if (min(eigen(DO)$values) <= 0) {
		stop("NaN value produces.")
	}
	log_det_OMEGA <- sum((p_vector - 1) * log(diag(AO))) + sum(log(eigen(DO)$values))
	trace_S_OMEGA <- SPLIT_UB_FOR_SUM(S, p_vector, AO, BO)
	
	if(drop_const == TRUE){
		log_likelihood <- log_det_OMEGA - trace_S_OMEGA
	} else{
		log_likelihood <- - (n * sum(p_vector)) / 2 * log(2 * pi) + (log_det_OMEGA - trace_S_OMEGA) * (n / 2)
	}

	return(- log_likelihood)
}

SOLVING_SCORE_FUN_FOR_GAMMA_OPTIM <- function(p_vector, S, n, iters_max, gamma_ini_lb, gamma_ini_ub, BLOCK_HADAMARD_PRODUCT, SPLIT_UB_FOR_SUM, NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, NEG_SCORE_FUN_FOR_GAMMA){

	G <- length(p_vector)
	dim_gamma <- G * (G + 1) / 2
	gamma_temp <- runif(dim_gamma, gamma_ini_lb, gamma_ini_ub)
	Neg_Log_Likelihood_temp <- NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA(gamma_temp, p_vector, S, n, SPLIT_UB_FOR_SUM, BLOCK_HADAMARD_PRODUCT)
	iters <- 1
	
	while(iters <= iters_max){
	
		tryCatch({
	
		repeat{
			RES_OPTIM <- optim(
									par = runif(dim_gamma, gamma_ini_lb, gamma_ini_ub), 
									 fn = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
									 gr = NEG_SCORE_FUN_FOR_GAMMA, 
								 method = "L-BFGS-B",
							### control = list(maxit = 1e5, factr = 1e1, pgtol = 0),
							    control = list(maxit = 1e5),
								hessian = FALSE,
						       p_vector = p_vector, 
									  S = S, 
									  n = n,			  
					   SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
                 BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT) 
			if(RES_OPTIM$convergence == 0){
				break
			}
		}
		gamma_update <- RES_OPTIM$par
		Neg_Log_Likelihood_update <- NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA(gamma_update, p_vector, S, n, SPLIT_UB_FOR_SUM, BLOCK_HADAMARD_PRODUCT)
		
		if(Neg_Log_Likelihood_update < Neg_Log_Likelihood_temp){
			Neg_Log_Likelihood_temp <- Neg_Log_Likelihood_update
			gamma_temp <- gamma_update
		} else {
			Neg_Log_Likelihood_temp <- Neg_Log_Likelihood_temp
			gamma_temp <- gamma_temp
		}
		
		### cat(' repetitions:', iters, "\r")
	
		iters <- iters + 1
		
		}, error = function(e){})
	}

	return(gamma_temp = gamma_temp)
}

SOLVING_AMBD_MODEL_FOR_GAMMA <- function(S, p_vector){
	
	K <- length(p_vector)
	R_MAT <- cov2cor(S)
	p_cumsum <- cumsum(p_vector)
	r_EST <- rep(0, K)
	
	for (k in 1 : K){
	
		if(k == 1){
			R_MAT_gg <- R_MAT[1 : p_cumsum[k], 1 : p_cumsum[k]]			  
		} else {
			R_MAT_gg <- R_MAT[(p_cumsum[k - 1] + 1) : p_cumsum[k], 
							  (p_cumsum[k - 1] + 1) : p_cumsum[k]]
		}
	r_EST[k] <- (sum(R_MAT_gg) - p_vector[k]) / (p_vector[k] * (p_vector[k] - 1))
	}

	d_EST <- r_EST / (1 - r_EST)
	
	c_EST <- - (p_vector * d_EST - d_EST + sqrt(p_vector * d_EST + 1) + 1) / (d_EST + (p_vector * d_EST + 1) * (p_vector - 2))
	
	gamma_EST <- as.vector(matrixcalc::vech(diag(- c_EST)))
	
	SIGMA_EST <- matrix(0, 1, 1)
	for(k in 1 : K){
		SIGMA_EST <- bdiag(SIGMA_EST, (1 - c_EST[k]) ^ (- 2) * (diag(p_vector[k]) + d_EST[k] * matrix(1, p_vector[k], p_vector[k])))
	}
	SIGMA_EST <- as.matrix(SIGMA_EST)[- 1, - 1]
	return(list(gamma_EST = gamma_EST, SIGMA_EST = SIGMA_EST))
}

SOLVING_MMRM_MODEL <- function(XT_mat, Y_mat, divided_by, n, p_vector, structure_model){

	R <- sum(p_vector)
	if(R > divided_by){
		if(R - sum(rep(divided_by, R %/% divided_by)) == 0){
			DIVIDED_SEQUENCES <- rep(divided_by, R %/% divided_by)
		} else {
			DIVIDED_SEQUENCES <- c(rep(divided_by, R %/% divided_by), R - sum(rep(divided_by, R %/% divided_by))) 
		### sum = R 
		}
		START <- c(1, (cumsum(DIVIDED_SEQUENCES) + 1)[- length(DIVIDED_SEQUENCES)])
		END <- cumsum(DIVIDED_SEQUENCES)
	} else {
		DIVIDED_SEQUENCES <- R
		START <- 1
		END <- R
	}
	
	
	EST_MMRM_RES <- c()
	VAR_MMRM_RES <- c()
	Sigma_EST_MMRM <- matrix(0, R, R)
	
	for(Seq in 1 : length(DIVIDED_SEQUENCES)){
	
		SELECTED_SAMPLE_SIZE <- n
		SELECTED_SAMPLE_ID <- seq(SELECTED_SAMPLE_SIZE)
		SELECTED_BIOMAKER_SIZE <- DIVIDED_SEQUENCES[Seq] 
		SELECTED_BIOMAKER_ID <- START[Seq] : END[Seq]
		SUBJECT_ID <- seq(n)
		BIOMARKERS_NAME <- seq(R)
	
		TIMES_VEC <- rep(SELECTED_BIOMAKER_SIZE, SELECTED_SAMPLE_SIZE)
		ID <- rep(SUBJECT_ID[SELECTED_SAMPLE_ID], TIMES_VEC)
		Biomarker <- rep(BIOMARKERS_NAME[SELECTED_BIOMAKER_ID], SELECTED_SAMPLE_SIZE)
		Intercept <- rep(XT_mat[SELECTED_SAMPLE_ID, 1], TIMES_VEC)
		Slope <- rep(XT_mat[SELECTED_SAMPLE_ID, 2], TIMES_VEC)
		Measure <- as.vector(Y_mat[SELECTED_BIOMAKER_ID, SELECTED_SAMPLE_ID])
			
		Biomarker <- factor(Biomarker, levels = BIOMARKERS_NAME[SELECTED_BIOMAKER_ID])

		NMR_DF <- data.frame(ID = factor(ID), 
					  Biomarker = factor(Biomarker), 
					  Intercept = Intercept, 
						  Slope = Slope)
			
		if(structure_model == "ar1"){
			
			fit_ar1 <- mmrm(formula = Measure ~ - 1 + Biomarker + Biomarker : Slope + ar1(Biomarker | ID),
							   data = NMR_DF,
							   reml = FALSE)
			### note the build-in p-values are calculated based on normal distribution rather than a t-distribution				   
		}

		EST_MMRM_RES_TEMP <- matrix(fit_ar1$beta_est, nrow = SELECTED_BIOMAKER_SIZE, ncol = p, byrow = FALSE)
		VAR_MMRM_RES_TEMP <- matrix(diag(fit_ar1$beta_vcov), nrow = SELECTED_BIOMAKER_SIZE, ncol = p, byrow = FALSE)
		
		EST_MMRM_RES <- rbind(EST_MMRM_RES, EST_MMRM_RES_TEMP)
		VAR_MMRM_RES <- rbind(VAR_MMRM_RES, VAR_MMRM_RES_TEMP)
		Sigma_EST_MMRM[START[Seq] : END[Seq], START[Seq] : END[Seq]] <- fit_ar1$cov
	}
	
	B_EST_MMRM <- EST_MMRM_RES
	SE_B_MMRM <- sqrt(VAR_MMRM_RES)
	return(list(B_EST_MMRM = B_EST_MMRM, SE_B_MMRM = SE_B_MMRM, Sigma_EST_MMRM = Sigma_EST_MMRM))
}

SOLVING_MLRM_MODEL <- function(XT_mat, YT_mat){

	RES_GLM_REG <- lm(YT_mat ~ - 1 + XT_mat)
	SSCPE <- crossprod(YT_mat - RES_GLM_REG$fitted.values)
	SigmaHat <- SSCPE / (nrow(YT_mat) - nrow(coef(RES_GLM_REG)))
	### in our setting, MLRM assumes a diagonal covariance matrix
	return(diag(diag(SigmaHat)))
}

####################################
### Setting the model parameters ###
####################################

SET_NO <- 1

SET_UP <- matrix(c(
		   100, 3, 10, 1.8,  NA, 0, 1e3, 2025,  ### study 1
		   200, 3, 10, 1.8,  NA, 0, 1e3, 2025,  ### study 1
		   300, 3, 10, 1.8,  NA, 0, 1e3, 2025,  ### study 1
		   100, 5,  5, 1.8,  NA, 0, 1e3, 2025,  ### study 1
		   200, 5,  5, 1.8,  NA, 0, 1e3, 2025,  ### study 1
		   300, 5,  5, 1.8,  NA, 0, 1e3, 2025,  ### study 1
			50, 3, 10, 1.8, 100, 0, 1e3, 2025,  ### study 2
			50, 3, 15, 1.8, 100, 0, 1e3, 2025,  ### study 2
			50, 3, 20, 1.8, 100, 0, 1e3, 2025,  ### study 2	
			50, 4, 10, 1.8, 100, 0, 1e3, 2025,  ### study 2	
			50, 4, 15, 1.8, 100, 0, 1e3, 2025,  ### study 2	
			50, 4, 20, 1.8, 100, 0, 1e3, 2025,  ### study 2	
			50, 3, 10, 1.8, 100, 1, 1e3, 2025   ### study 3
		), nrow = 13, ncol = 8, byrow = TRUE)
### cols: sample size, G, multiple, kappa, sigma_min, eps, reps_max, model_seed

n <- SET_UP[SET_NO, 1]
### sample size, n > max{p = 2, G(G + 1) / 2 = 6, 15, 55}
if (SET_UP[SET_NO, 2] == 3){
	l_vec <- c(3, 3, 4)  
	### vector of group sizes
} else if (SET_UP[SET_NO, 2] == 4){
	l_vec <- c(2, 3, 2, 3)
} else if (SET_UP[SET_NO, 2] == 5){
	l_vec <- c(2, 2, 2, 2, 2)
}
multiple <- SET_UP[SET_NO, 3]
l_vec <- l_vec * multiple				   
effect_size <- SET_UP[SET_NO, 4]
### mean over standard deviation
sigma_min_scale <- SET_UP[SET_NO, 5]
### components of Sigma_0 are between sigma_min_scale and 2 times
eps <- SET_UP[SET_NO, 6] * 1e-2
### noise level for misspecification
reps_max <- SET_UP[SET_NO, 7]
### number of replications for Monte Carlo simulations
model_seed <- SET_UP[SET_NO, 8]
### seed for generating random number

sig_level <- 0.05		 			
### significance level
p <- 2			     		    
### number of covariate variables: 1 and cont.
L_mat <- diag(l_vec) 
G <- length(l_vec)   
### number of groups
R <- sum(l_vec)      
### total number of features, R = 10 * multiple
dim_gamma <- G * (G + 1) / 2
### G = 3, dim_gamma = 6; G = 5, dim_gamma = 15

if(G == 3 && multiple == 10 && is.na(sigma_min_scale)){
	gamma_0 <- c(0.399841830, 0.005572861, -0.511547315, 0.188360052, -0.908824807, -0.637345389)
} else if (G == 5 && multiple == 5 && is.na(sigma_min_scale)){
	gamma_0 <- c(1.724321785, 0.005205265, 0.565925845, 1.972455421, 0.335524385, 1.351378041, 0.677802318, 0.072193136, -0.718428492, 1.699842690, -0.145024523, 1.586879025, 0.578328007, 1.648949168, 1.694381397)
} else if (G == 5 && multiple == 10 && is.na(sigma_min_scale)){
	gamma_0 <- c(0.86879590, -0.57277313, -0.04493324, 0.11241645, 1.19520147, -0.30809017, 0.80493488, 1.43089400, 1.61352353, -0.84699558, 1.21450331, -0.60276354, 1.09865097, 0.82983710, -0.17773762)
} else if (G == 3 && multiple == 10 && sigma_min_scale == 100){
	gamma_0 <- c(0.6299008, 1.7586950, -0.6530756, 1.2432015, 0.6358242, -0.9034923)
} else if(G == 3 && multiple == 15 && sigma_min_scale == 100){
	gamma_0 <- c(1.0530112, 0.4515030, -0.3041058, 1.5767517, 0.8114608,  0.8035750)
} else if(G == 3 && multiple == 20 && sigma_min_scale == 100){
    gamma_0 <- c(1.8718713, 0.9953679, 0.3631638, 1.7293809, -0.8711309,  1.1011621)
} else if(G == 4 && multiple == 10 && sigma_min_scale == 100){
    gamma_0 <- c(-0.2203274, -0.8097071, 1.4816519, 0.4267440, 0.5997492,  0.9599127, 1.9112539, -0.5376684, -0.3865847, 1.6232989)
} else if(G == 4 && multiple == 15 && sigma_min_scale == 100){
    gamma_0 <- c(0.18937762, 0.80057518, -0.34700112, 1.02245594, 0.36014646, 1.10882481, 0.01185789, 1.50743013, -0.21671494, 0.81410686)
} else if(G == 4 && multiple == 20 && sigma_min_scale == 100){
    gamma_0 <- c(0.90113596, 1.94427826, 0.13423611, 1.41693515, 1.41889462, -0.06046755, -0.23999089, 1.74623299, 0.95604457, -0.80917951)
}

gamma_ini_lb <- - 1
### the lower bound for the random initial values for gamma_0
gamma_ini_ub <- 2
### the upper bound for the random initial values for gamma_0

####################################################################
### Setting the mean vector and covariance matrix for dependence ###
####################################################################

Mu_0 <- rep(0, R)

BU_0 <- matlib::symMat(gamma_0, diag = TRUE, byrow = FALSE) 
### direction by ROW, the function symMat is incorrect
AU_0 <- diag(- diag(BU_0))
AIU_0 <- diag(G) - AU_0
BIU_0 <- - BU_0
AO_0 <- AIU_0 %*% AIU_0
BO_0 <- AIU_0 %*% BIU_0 + BIU_0 %*% AIU_0 + BIU_0 %*% L_mat %*% BIU_0
DO_0 <- AO_0 + BO_0 %*% L_mat
AS_0 <- solve(AO_0)
BS_0 <- - solve(DO_0) %*% BO_0 %*% solve(AO_0)
if (!isSymmetric(BS_0)) {
		BS_0 <- (BS_0 + t(BS_0)) / 2
}

if(eps == 0){
	Sigma_0 <- BLOCK_HADAMARD_PRODUCT(AS_0, BS_0, l_vec)
} else {
	set.seed(10000)
	E_3 <- rWishart(n = 1, df = R, Sigma = 3 * eps * diag(R))[,,1]
	E_6 <- rWishart(n = 1, df = R, Sigma = 6 * eps * diag(R))[,,1]
	E_9 <- rWishart(n = 1, df = R, Sigma = 9 * eps * diag(R))[,,1]
	Sigma_0 <- BLOCK_HADAMARD_PRODUCT(AS_0, BS_0, l_vec) 
	Sigma_3 <- Sigma_0 + E_3
	Sigma_6 <- Sigma_0 + E_6
	Sigma_9 <- Sigma_0 + E_9
}

#################################
### Setting the design matrix ###
#################################

prop_non_null <- 0.30
### H_0j: Beta_1j = 0 for j = 1, ..., R
### proportion of non-null/non-zero/non-true Hj
### When generating data, Beta_1j = 0 refers to truth
X <- matrix(rnorm((p - 1) * n), nrow = p - 1, ncol = n)
X <- rbind(t(rep(1, n)), X) 
### p by n design matrix, i-th column = t(xi), first element of xi = 1
X[2, ] <- (X[2, ] - mean(X[2, ])) / sqrt(mean((X[2, ] - mean(X[2, ])) ^ 2))
### such that sum xi * t(xi) = n * I2
### It works for discrete variable, so omit it
inv_x_sum_mat <- diag(p) / n
### p by p, p = 2 only
X_mat <- X 
### p by n, each column = p by 1 of xi
XT_mat <- t(X_mat) 
### n by p

###########################################
### Setting the regression coefficients ###
###########################################

COV_beta_0 <- kronecker(Sigma_0, inv_x_sum_mat)
SE_beta_0 <- sqrt(diag(COV_beta_0))
Beta_0_temp <- matrix(0, nrow = R, ncol = p) 
l_cumsum <- cumsum(l_vec)
for(g in 1 : G){
	if(g == 1){
		Beta_0_temp[1 : (floor(l_vec[1] * prop_non_null) + 1), ] <- rep(1, floor(l_vec[1] * prop_non_null) + 1) %*% t(c(0.95, 1))
	} else{
		Beta_0_temp[(l_cumsum[g - 1] + 1) : (l_cumsum[g - 1] + floor(l_vec[g] * prop_non_null) + 1), ] <- rep(1, (floor(l_vec[g] * prop_non_null) + 1)) %*% t(c(0.95, 1))
	}
}
beta_0 <- as.vector(matrixcalc::vec(t(Beta_0_temp))) * effect_size * SE_beta_0
### (Rp) by 1 true values of beta
Beta_0 <- matrix(beta_0, nrow = R, ncol = p, byrow = TRUE) 
### direction by ROW, R by p true values of Beta


######################
### Other settings ###
######################
OUTPUT_ADDRESS <- "..."
### the default output address
reps <- 1
iters_max <- 1e2
### the number of initial values for each estimate of gamma_0

STUDY_1_NO <- 1 : 6
gamma_MCMC_MAT <- matrix(0, nrow = reps_max, ncol = dim_gamma)
### a matrix contains the Monte Carlo gamma estimates
PHI_gamma_AASY_1 <- PHI_gamma_AASY_2 <- PHI_gamma_AASY_3 <- matrix(0, nrow = reps_max, ncol = dim_gamma * (dim_gamma + 1) / 2)
### a matrix contains the triangle elements of asy. covariance matrix 
gamma_EWCP_MAT_1 <- gamma_EWCP_MAT_2 <- gamma_EWCP_MAT_3 <- matrix(0, nrow = reps_max, ncol = dim_gamma)
### a matrix contains if the 95% empirical Wald-type CI includes the true
beta_est_vec <- rep(0, length(beta_0))
### the estimate of beta

STUDY_2_NO <- 7 : 12
p_adj_method <- "BH"
### how to adjust the p-values
decision_mat <- matrix(0, nrow = R * 1, ncol = 5)
colnames(decision_mat) <- c("TRUE", "MAUD", "AMBD", "MMRM", "MLRM")
loss_mat_F <- loss_mat_S <- matrix(0, nrow = reps_max, ncol = 5)
colnames(loss_mat_F) <- colnames(loss_mat_S) <- c("TRUE", "MAUD", "AMBD", "MMRM", "MLRM")
computing_time_mat <- matrix(0, nrow = reps_max, ncol = 5)

STUDY_3_NO <- 13
mis_loss_mat_F <- mis_loss_mat_S <- matrix(0, nrow = reps_max, ncol = 4)
colnames(mis_loss_mat_F) <- colnames(mis_loss_mat_S) <- c("sig_0", "sig_3", "sig_6", "sig_9")


#################################
### Generating simulated data ###
#################################

start_time <- Sys.time()
set.seed(model_seed)

while(reps <= reps_max){

	tryCatch({	
	
		if(SET_NO %in% c(STUDY_1_NO, STUDY_2_NO)){
	
			E_mat_0 <- t(mvnfast::rmvn(n = n, mu = Mu_0, sigma = Sigma_0))
			### R by n
			Y_mat <- Beta_0 %*% X_mat + E_mat_0
			### R by n
			YT_mat <- t(Y_mat) 
			### n by R
			BetaT_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% YT_mat
			Beta_OLS <- t(BetaT_OLS)							 		   
			### R by p
			E_mat <- Y_mat - Beta_OLS %*% X_mat
			### R by n
			ET_mat <- t(E_mat)   				       
			### n by R
			S <- (n - 1) * cov(ET_mat) / n 					     	   
			### R by R
			beta_OLS <- as.vector(BetaT_OLS)                     
			### by row of Beta_OLS
		
		} else if(SET_NO %in% STUDY_3_NO){
		
			E_mat_0 <- t(mvnfast::rmvn(n = n, mu = Mu_0, sigma = Sigma_0))
			E_mat_3 <- t(mvnfast::rmvn(n = n, mu = Mu_0, sigma = Sigma_3))
			E_mat_6 <- t(mvnfast::rmvn(n = n, mu = Mu_0, sigma = Sigma_6))
			E_mat_9 <- t(mvnfast::rmvn(n = n, mu = Mu_0, sigma = Sigma_9))
			### R by n
			Y_0_mat <- Beta_0 %*% X_mat + E_mat_0
			Y_3_mat <- Beta_0 %*% X_mat + E_mat_3
			Y_6_mat <- Beta_0 %*% X_mat + E_mat_6
			Y_9_mat <- Beta_0 %*% X_mat + E_mat_9
			### R by n
			Y_0_T_mat <- t(Y_0_mat) 
			Y_3_T_mat <- t(Y_3_mat)
			Y_6_T_mat <- t(Y_6_mat)
			Y_9_T_mat <- t(Y_9_mat)
			### n by R
			Beta_0_T_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% Y_0_T_mat
			Beta_3_T_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% Y_3_T_mat
			Beta_6_T_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% Y_6_T_mat
			Beta_9_T_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% Y_9_T_mat
			
			Beta_0_OLS <- t(Beta_0_T_OLS)
			Beta_3_OLS <- t(Beta_3_T_OLS)
			Beta_6_OLS <- t(Beta_6_T_OLS)
			Beta_9_OLS <- t(Beta_9_T_OLS)			
			### R by p
			E_0_mat <- Y_0_mat - Beta_0_OLS %*% X_mat
			E_3_mat <- Y_3_mat - Beta_3_OLS %*% X_mat
			E_6_mat <- Y_6_mat - Beta_6_OLS %*% X_mat
			E_9_mat <- Y_9_mat - Beta_9_OLS %*% X_mat
			### R by n
			E_0_T_mat <- t(E_0_mat)  
			E_3_T_mat <- t(E_3_mat) 
			E_6_T_mat <- t(E_6_mat) 
			E_9_T_mat <- t(E_9_mat) 			
			### n by R
			S_0 <- (n - 1) * cov(E_0_T_mat) / n	
			S_3 <- (n - 1) * cov(E_3_T_mat) / n
			S_6 <- (n - 1) * cov(E_6_T_mat) / n
			S_9 <- (n - 1) * cov(E_9_T_mat) / n			
			### R by R
			beta_0_OLS <- as.vector(Beta_0_T_OLS)
			beta_3_OLS <- as.vector(Beta_3_T_OLS) 
			beta_6_OLS <- as.vector(Beta_6_T_OLS) 
			beta_9_OLS <- as.vector(Beta_9_T_OLS) 
			### by row of B_OLS
		}
		
		
		if(SET_NO %in% STUDY_1_NO){
		
			### 1 MAUD ###
			### 1-1 estimating gamma vector ###
			gamma_EST_MAUD <- SOLVING_SCORE_FUN_FOR_GAMMA_OPTIM(
						   p_vector = l_vec, 
								  S = S, 
								  n = n, 
						  iters_max = iters_max, 
					   gamma_ini_lb = gamma_ini_lb, 
					   gamma_ini_ub = gamma_ini_ub, 
			 BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, 
				   SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
   NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
			NEG_SCORE_FUN_FOR_GAMMA = NEG_SCORE_FUN_FOR_GAMMA)
		gamma_MCMC_MAT[reps, ] <- gamma_EST_MAUD	
		
			### 1-2-1 estimating asy. cov(gamma) matrix using the formula ###
			BU_EST_MAUD <- matlib::symMat(gamma_EST_MAUD, diag = TRUE, byrow = FALSE) 
			### direction by ROW, the function symMat is incorrect
			AU_EST_MAUD <- diag(- diag(BU_EST_MAUD))
			AIU_EST_MAUD <- diag(G) - AU_EST_MAUD
			BIU_EST_MAUD <- - BU_EST_MAUD
			AO_EST_MAUD <- AIU_EST_MAUD %*% AIU_EST_MAUD
			BO_EST_MAUD <- AIU_EST_MAUD %*% BIU_EST_MAUD + BIU_EST_MAUD %*% AIU_EST_MAUD + BIU_EST_MAUD %*% L_mat %*% BIU_EST_MAUD
			DO_EST_MAUD <- AO_EST_MAUD + BO_EST_MAUD %*% L_mat
			AS_EST_MAUD <- solve(AO_EST_MAUD)
			BS_EST_MAUD <- - solve(DO_EST_MAUD) %*% BO_EST_MAUD %*% solve(AO_EST_MAUD)
			if (!isSymmetric(BS_EST_MAUD)) {
				BS_EST_MAUD <- (BS_EST_MAUD + t(BS_EST_MAUD)) / 2
			}
			PHI_gamma_AASY_TEMP <- FISHER_INFORMATION_MATRIX_FOR_GAMMA(
									 AU = AU_EST_MAUD, 
									 BU = BU_EST_MAUD, 
								     AS = AS_EST_MAUD, 
									 BS = BS_EST_MAUD,
							   p_vector = l_vec, 
									  n = n, 
			PRODUCT_OF_FOUR_UB_MATRICES = PRODUCT_OF_FOUR_UB_MATRICES)
			### Phi_gamma
			PHI_gamma_AASY_1[reps, ] <- as.vector(matrixcalc::vech(solve(PHI_gamma_AASY_TEMP)))
			### = cov(gamma_EST) = Phi_gamma^(-1)
			
			### 1-2-2 estimating asy. cov(gamma) matrix using the optim-hessian ###
			PHI_gamma_AASY_OPTIM <- optimHess(
							par = gamma_EST_MAUD, 
							 fn = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
							 gr = NEG_SCORE_FUN_FOR_GAMMA, 
					   p_vector = l_vec, 
							  S = S, 
							  n = n, 
			   SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
		 BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT) 
			### it is already - Hessian when using control = list(fnscale = -1)
			PHI_gamma_AASY_2[reps, ] <- as.vector(matrixcalc::vech(solve(PHI_gamma_AASY_OPTIM)))
			
			### 1-2-3 estimating asy. cov(gamma) matrix using the formula at truth ###
			PHI_gamma_AASY_TRUTH <- FISHER_INFORMATION_MATRIX_FOR_GAMMA(
									 AU = AU_0, 
									 BU = BU_0, 
								     AS = AS_0, 
									 BS = BS_0,
							   p_vector = l_vec, 
									  n = n, 
			PRODUCT_OF_FOUR_UB_MATRICES = PRODUCT_OF_FOUR_UB_MATRICES)
			PHI_gamma_AASY_3[reps, ] <- as.vector(matrixcalc::vech(solve(PHI_gamma_AASY_TRUTH)))
			
			### 1-3 computing the empirical WCP ###
			for(i in 1 : 3){
				if(i == 1){
					SE_gamma_EST_MAUD <- sqrt(diag(solve(PHI_gamma_AASY_TEMP)))
					CI_lower_temp <- gamma_EST_MAUD - qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_upper_temp <- gamma_EST_MAUD + qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_temp <- rbind(CI_lower_temp, CI_upper_temp)
					gamma_EWCP_temp <- 1 * (CI_temp[1, ] < gamma_0 & CI_temp[2, ] > gamma_0)
					gamma_EWCP_MAT_1[reps, ] <- gamma_EWCP_temp
				} else if(i == 2){
					SE_gamma_EST_MAUD <- sqrt(diag(solve(PHI_gamma_AASY_OPTIM)))
					CI_lower_temp <- gamma_EST_MAUD - qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_upper_temp <- gamma_EST_MAUD + qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_temp <- rbind(CI_lower_temp, CI_upper_temp)
					gamma_EWCP_temp <- 1 * (CI_temp[1, ] < gamma_0 & CI_temp[2, ] > gamma_0)
					gamma_EWCP_MAT_2[reps, ] <- gamma_EWCP_temp
				} else if(i == 3){
					SE_gamma_EST_MAUD <- sqrt(diag(solve(PHI_gamma_AASY_TRUTH)))
					CI_lower_temp <- gamma_EST_MAUD - qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_upper_temp <- gamma_EST_MAUD + qnorm(1 - sig_level / 2) * SE_gamma_EST_MAUD
					CI_temp <- rbind(CI_lower_temp, CI_upper_temp)
					gamma_EWCP_temp <- 1 * (CI_temp[1, ] < gamma_0 & CI_temp[2, ] > gamma_0)
					gamma_EWCP_MAT_3[reps, ] <- gamma_EWCP_temp
				}
			}
			
			
			### 1-4 computing the bias for beta
			beta_est_vec <- beta_est_vec + beta_OLS
		
		} else if(SET_NO %in% STUDY_2_NO){

			### 1 TRUE ###
			time_start_TRUE <- Sys.time()
			loss_mat_F[reps, 1] <- 0 
			loss_mat_S[reps, 1] <- 0
			T_Beta_TRUE <- matrix(beta_OLS / sqrt(diag(COV_beta_0)), nrow = R, ncol = p, byrow = TRUE) 
			### direction by ROW, R by p t-values		
			decision_mat[, 1] <- decision_mat[, 1] + as.vector((p.adjust(2 * (1 - pnorm(abs(T_Beta_TRUE[, 2]))), method = p_adj_method) < sig_level) * 1)
			computing_time_mat[reps, 1] <- Sys.time() - time_start_TRUE
			
			### 2 MAUD ###
			time_start_MAUD <- Sys.time()
			gamma_EST_MAUD <- SOLVING_SCORE_FUN_FOR_GAMMA_OPTIM(
						   p_vector = l_vec, 
								  S = S, 
								  n = n, 
						  iters_max = iters_max, 
					   gamma_ini_lb = gamma_ini_lb, 
					   gamma_ini_ub = gamma_ini_ub, 
			 BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, 
				   SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
   NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
			NEG_SCORE_FUN_FOR_GAMMA = NEG_SCORE_FUN_FOR_GAMMA)
			BU_EST_MAUD <- matlib::symMat(gamma_EST_MAUD, diag = TRUE, byrow = FALSE) 
			### direction by ROW, the function symMat is incorrect
			AU_EST_MAUD <- diag(- diag(BU_EST_MAUD))
			AIU_EST_MAUD <- diag(G) - AU_EST_MAUD
			BIU_EST_MAUD <- - BU_EST_MAUD
			AO_EST_MAUD <- AIU_EST_MAUD %*% AIU_EST_MAUD
			BO_EST_MAUD <- AIU_EST_MAUD %*% BIU_EST_MAUD + BIU_EST_MAUD %*% AIU_EST_MAUD + BIU_EST_MAUD %*% L_mat %*% BIU_EST_MAUD
			DO_EST_MAUD <- AO_EST_MAUD + BO_EST_MAUD %*% L_mat
			AS_EST_MAUD <- solve(AO_EST_MAUD)
			BS_EST_MAUD <- - solve(DO_EST_MAUD) %*% BO_EST_MAUD %*% solve(AO_EST_MAUD)
			if (!isSymmetric(BS_EST_MAUD)) {
				BS_EST_MAUD <- (BS_EST_MAUD + t(BS_EST_MAUD)) / 2
			}
			Sigma_EST_MAUD <- BLOCK_HADAMARD_PRODUCT(AS_EST_MAUD, BS_EST_MAUD, l_vec)
			COV_beta_EST_MAUD <- kronecker(Sigma_EST_MAUD, inv_x_sum_mat)
			loss_mat_F[reps, 2] <- norm(COV_beta_EST_MAUD - COV_beta_0, "F") / norm(COV_beta_0, "F")
			loss_mat_S[reps, 2] <- norm(COV_beta_EST_MAUD - COV_beta_0, "2") / norm(COV_beta_0, "2")
			T_Beta_MAUD <- matrix(beta_OLS / sqrt(diag(COV_beta_EST_MAUD)), nrow = R, ncol = p, byrow = TRUE) 
			### direction by ROW, R by p t-values		
			decision_mat[, 2] <- decision_mat[, 2] + as.vector((p.adjust(2 * (1 - pt(abs(T_Beta_MAUD[, 2]), df = n - 1)), method = p_adj_method) < sig_level) * 1)
			computing_time_mat[reps, 2] <- Sys.time() - time_start_MAUD
		
			### 3 AMBD ###
			time_start_AMBD <- Sys.time()
			RES_AMBD <- SOLVING_AMBD_MODEL_FOR_GAMMA(
					S = S, 
			 p_vector = l_vec)
			gamma_EST_AMBD <- RES_AMBD$gamma_EST
			Sigma_EST_AMBD <- RES_AMBD$SIGMA_EST
			COV_beta_EST_AMBD <- kronecker(Sigma_EST_AMBD, inv_x_sum_mat)
			loss_mat_F[reps, 3] <- norm(COV_beta_EST_AMBD - COV_beta_0, "F") / norm(COV_beta_0, "F")
			loss_mat_S[reps, 3] <- norm(COV_beta_EST_AMBD - COV_beta_0, "2") / norm(COV_beta_0, "2")
			T_Beta_AMBD <- matrix(beta_OLS / sqrt(diag(COV_beta_EST_AMBD)), nrow = R, ncol = p, byrow = TRUE) 
			decision_mat[, 3] <- decision_mat[, 3] + as.vector((p.adjust(2 * (1 - pt(abs(T_Beta_AMBD[, 2]), df = n - 1)), method = p_adj_method) < sig_level) * 1)
			computing_time_mat[reps, 3] <- Sys.time() - time_start_AMBD
			
			### 4 MMRM ###
			time_start_MMRM <- Sys.time()
			RES_MMRM <- SOLVING_MMRM_MODEL(
					XT_mat = XT_mat, 
					 Y_mat = Y_mat, 
				divided_by = 50, 
						 n = n, 
				  p_vector = l_vec, 
		   structure_model = "ar1")
			Sigma_EST_MMRM <- RES_MMRM$Sigma_EST_MMRM
			COV_beta_EST_MMRM <- kronecker(Sigma_EST_MMRM, inv_x_sum_mat)
			loss_mat_F[reps, 4] <- norm(COV_beta_EST_MMRM - COV_beta_0, "F") / norm(COV_beta_0, "F")
			loss_mat_S[reps, 4] <- norm(COV_beta_EST_MMRM - COV_beta_0, "2") / norm(COV_beta_0, "2")
			beta_EST_MMRM <- as.vector(matrixcalc::vec(t(RES_MMRM$B_EST_MMRM)))
			SE_beta_MMRM <- as.vector(matrixcalc::vec(t(RES_MMRM$SE_B_MMRM)))
			T_Beta_MMRM <- matrix(beta_EST_MMRM / SE_beta_MMRM, nrow = R, ncol = p, byrow = TRUE) 
			decision_mat[, 4] <- decision_mat[, 4] + as.vector((p.adjust(2 * (1 - pt(abs(T_Beta_MMRM[, 2]), df = n - 1)), method = p_adj_method) < sig_level) * 1)
			computing_time_mat[reps, 4] <- Sys.time() - time_start_MMRM
			
			### 5 MLRM ###
			time_start_MLRM <- Sys.time()
			Sigma_EST_MLRM <- SOLVING_MLRM_MODEL(XT_mat, YT_mat)
			COV_beta_EST_MLRM <- kronecker(Sigma_EST_MLRM, inv_x_sum_mat)
			loss_mat_F[reps, 5] <- norm(COV_beta_EST_MLRM - COV_beta_0, "F") / norm(COV_beta_0, "F")
			loss_mat_S[reps, 5] <- norm(COV_beta_EST_MLRM - COV_beta_0, "2") / norm(COV_beta_0, "2")
			T_Beta_MLRM <- matrix(beta_OLS / sqrt(diag(COV_beta_EST_MLRM)), nrow = R, ncol = p, byrow = TRUE) 
			decision_mat[, 5] <- decision_mat[, 5] + as.vector((p.adjust(2 * (1 - pt(abs(T_Beta_MLRM[, 2]), df = n - 1)), method = p_adj_method) < sig_level) * 1)
			computing_time_mat[reps, 5] <- Sys.time() - time_start_MLRM
			
		} else if(SET_NO %in% STUDY_3_NO){
			
			for(sig in c(0, 3, 6, 9)){
			
				if(sig == 0){
					S <- S_0
					indx <- 1
				} else if(sig == 3){
					S <- S_3
					indx <- 2
				} else if(sig == 6){
					S <- S_6
					indx <- 3
				} else if(sig == 9){
					S <- S_9
					indx <- 4
				}
				
				gamma_EST_MAUD <- SOLVING_SCORE_FUN_FOR_GAMMA_OPTIM(
						   p_vector = l_vec, 
								  S = S, 
								  n = n, 
						  iters_max = iters_max, 
					   gamma_ini_lb = gamma_ini_lb, 
					   gamma_ini_ub = gamma_ini_ub, 
			 BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT, 
				   SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
   NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
			NEG_SCORE_FUN_FOR_GAMMA = NEG_SCORE_FUN_FOR_GAMMA)
				BU_EST_MAUD <- matlib::symMat(gamma_EST_MAUD, diag = TRUE, byrow = FALSE) 
				### direction by ROW, the function symMat is incorrect
				AU_EST_MAUD <- diag(- diag(BU_EST_MAUD))
				AIU_EST_MAUD <- diag(G) - AU_EST_MAUD
				BIU_EST_MAUD <- - BU_EST_MAUD
				AO_EST_MAUD <- AIU_EST_MAUD %*% AIU_EST_MAUD
				BO_EST_MAUD <- AIU_EST_MAUD %*% BIU_EST_MAUD + BIU_EST_MAUD %*% AIU_EST_MAUD + BIU_EST_MAUD %*% L_mat %*% BIU_EST_MAUD
				DO_EST_MAUD <- AO_EST_MAUD + BO_EST_MAUD %*% L_mat
				AS_EST_MAUD <- solve(AO_EST_MAUD)
				BS_EST_MAUD <- - solve(DO_EST_MAUD) %*% BO_EST_MAUD %*% solve(AO_EST_MAUD)
				if (!isSymmetric(BS_EST_MAUD)) {
					BS_EST_MAUD <- (BS_EST_MAUD + t(BS_EST_MAUD)) / 2
				}
				Sigma_EST_MAUD <- BLOCK_HADAMARD_PRODUCT(AS_EST_MAUD, BS_EST_MAUD, l_vec)
				COV_beta_EST_MAUD <- kronecker(Sigma_EST_MAUD, inv_x_sum_mat)
				mis_loss_mat_F[reps, indx] <- norm(COV_beta_EST_MAUD - COV_beta_0, "F") / norm(COV_beta_0, "F")
				mis_loss_mat_S[reps, indx] <- norm(COV_beta_EST_MAUD - COV_beta_0, "2") / norm(COV_beta_0, "2")
			}
			
		}
	
		cat(' repetitions:', reps, "\r")

		reps <- reps + 1

	}, error = function(e){})
}

end_time <- Sys.time() - start_time
end_time

FILENAME <- paste("SIMULATION_RESULTS_SETUPNO", SET_NO, "reps_max", reps_max, sep = "_")
FILENAME <- paste(FILENAME, ".RData", sep = "")
save.image(paste(OUTPUT_ADDRESS, FILENAME, sep = "/"))

if(SET_NO %in% STUDY_1_NO){

	gamma_BIAS <- apply(gamma_MCMC_MAT, 2, mean) - gamma_0
	gamma_MCSD <- apply(gamma_MCMC_MAT, 2, sd)
	### gamma_AASY <- sqrt(diag(matlib::symMat(apply(PHI_gamma_AASY, 2, mean), diag = TRUE, byrow = FALSE)))
	gamma_AASY_1 <- apply(apply(PHI_gamma_AASY_1, 1, function(x) sqrt(diag(matlib::symMat(x, diag = TRUE, byrow = FALSE)))), 1, mean)
	gamma_EWCP_1 <- apply(gamma_EWCP_MAT_1, 2, mean)
	
	gamma_AASY_2 <- apply(apply(PHI_gamma_AASY_2, 1, function(x) sqrt(diag(matlib::symMat(x, diag = TRUE, byrow = FALSE)))), 1, mean, na.rm = TRUE)
	gamma_EWCP_2 <- apply(gamma_EWCP_MAT_2, 2, mean, na.rm = TRUE)
	gamma_AASY_3 <- apply(apply(PHI_gamma_AASY_3, 1, function(x) sqrt(diag(matlib::symMat(x, diag = TRUE, byrow = FALSE)))), 1, mean, na.rm = TRUE)
	gamma_EWCP_3 <- apply(gamma_EWCP_MAT_3, 2, mean, na.rm = TRUE)
	
	print(round(cbind(gamma_BIAS, gamma_MCSD, gamma_AASY_1, gamma_EWCP_1, gamma_AASY_2, gamma_EWCP_2, gamma_AASY_3, gamma_EWCP_3) * 100, 2))
	
	print(norm(beta_est_vec / reps_max - beta_0, "2"))

} else if(SET_NO %in% STUDY_2_NO){

	PLOT_DATA <- as.data.frame(cbind(
	c(loss_mat_F[, 2], loss_mat_F[, 3], loss_mat_F[, 4], loss_mat_F[, 5],
	  loss_mat_S[, 2], loss_mat_S[, 3], loss_mat_S[, 4], loss_mat_S[, 5]), 	
	c(rep("MAUD", reps_max), rep("AMBD", reps_max), rep("MMRM", reps_max), rep("MLRM", reps_max), rep("MAUD", reps_max), rep("AMBD", reps_max), rep("MMRM", reps_max), rep("MLRM", reps_max)),
	c(rep("Frobenius Norm", 4 * reps_max), rep("Spectral Norm", 4 * reps_max))))
	PLOT_DATA[, 1] <- as.numeric(PLOT_DATA[, 1])
	PLOT_DATA[, 2] <- factor(PLOT_DATA[, 2], levels = c("MAUD", "AMBD", "MMRM", "MLRM"))
	colnames(PLOT_DATA) <- c("relative_loss", "method", "norm")

	font_scale <- 20
	plot_size_in <- 6

	SVGNAME <- paste("SIMULATION_RESULTS_SETUPNO", SET_NO, "reps_max", reps_max, sep = "_")
	SVGNAME <- paste(SVGNAME, ".svg", sep = "")
	svglite::svglite(file = paste(OUTPUT_ADDRESS, SVGNAME, sep = "/"), 
					width = plot_size_in, 
				   height = plot_size_in) 
	ggplot2::ggplot(data = PLOT_DATA, 
				 mapping = aes(x = norm, 
							   y = relative_loss, 
							fill = method)) +
	geom_boxplot() +
	theme_bw() + 
	# theme(legend.position.inside = "top", legend.justification = "right") +
	theme(legend.position.inside = c(0.85, 0.85)) +
	theme(axis.title.x = element_blank()) +
	ylab("Relative Loss") + 
	theme(text = element_text(size = font_scale)) 
	dev.off()
	
	print(cbind(Beta_0[, 2], decision_mat / reps_max))
	### first non-zero: empirical powers
	### first zero: empirical sizes
	
	print(apply(computing_time_mat, 2, mean))
	print(apply(computing_time_mat, 2, sd))
    ### TRUE, MAUD, AMBD, MMRM, MLRM
	
} else if(SET_NO %in% STUDY_3_NO){

	MIS_PLOT_DATA <- as.data.frame(cbind(
	c(mis_loss_mat_F[, 1], mis_loss_mat_F[, 2], mis_loss_mat_F[, 3], mis_loss_mat_F[, 4], mis_loss_mat_S[, 1], mis_loss_mat_S[, 2], mis_loss_mat_S[, 3], mis_loss_mat_S[, 4]), 	
	c(rep("0 level", reps_max), rep("3 level", reps_max), rep("6 level", reps_max), rep("9 level", reps_max), rep("0 level", reps_max), rep("3 level", reps_max), rep("6 level", reps_max), rep("9 level", reps_max)),
	c(rep("Frobenius Norm", 4 * reps_max), rep("Spectral Norm", 4 * reps_max))))
	MIS_PLOT_DATA[, 1] <- as.numeric(MIS_PLOT_DATA[, 1])
	MIS_PLOT_DATA[, 2] <- factor(MIS_PLOT_DATA[, 2], levels = c("0 level", "3 level", "6 level", "9 level"))
	colnames(MIS_PLOT_DATA) <- c("relative_loss", "noise", "norm")

	font_scale <- 20
	plot_size_in <- 6

	SVGNAME <- paste("SIMULATION_RESULTS_SETUPNO", SET_NO, "reps_max", reps_max, sep = "_")
	SVGNAME <- paste(SVGNAME, ".svg", sep = "")
	svglite::svglite(file = paste(OUTPUT_ADDRESS, SVGNAME, sep = "/"), 
					width = plot_size_in, 
				   height = plot_size_in) 
	ggplot2::ggplot(data = MIS_PLOT_DATA, 
				 mapping = aes(x = norm, 
							   y = relative_loss, 
							fill = noise)) +
	geom_boxplot() +
	theme_bw() + 
	#theme(legend.position.inside = "top", legend.justification = "right") +
	theme(legend.position.inside = c(0.10, 0.88)) +
	theme(axis.title.x = element_blank()) +
	ylab("Relative Loss") + 
	theme(text = element_text(size = font_scale)) 
	dev.off()
}







