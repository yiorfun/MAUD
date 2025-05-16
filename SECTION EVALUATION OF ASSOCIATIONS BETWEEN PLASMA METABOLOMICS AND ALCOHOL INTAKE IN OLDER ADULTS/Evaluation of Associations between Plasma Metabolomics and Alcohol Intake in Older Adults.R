#################################################
### Real Data Analyses for NMR Biomarker Data ###
#################################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"readr",
	### readr::read_csv(), read data from a CSV file
	"readxl", 
	### readxl::read_excel(), read data from an Excel file
	"xlsx",
	### xlsx::write.xlsx() saves the results in a sheet
	"matlib", 
	### matlib::symMat(), create a symmetric matrix from a vector
	"R.matlab", 
	### R.matlab, convert between R files and MATLAB files
	"matrixcalc",		
	### matrixcalc::vech(), create a vector from a symmetric matrix
	"mmrm",
	### mmrm::mmrm(), perform the mmrm model with various covariance
	"ggplot2",
	### ggplot2::ggplot(), create the forest plots
	"stringr"
	### stringr::str_detect() check a string vector for a specific sub-string
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

BEST_UNBIASED_ESTIMATOR <- function(S, p_vec){
	
	K <- length(p_vec)
	A_temp <- matrix(0, K, K)
	B_temp <- matrix(0, K, K)
	for(k in 1 : K){
		for(kp in 1 : K){
			SUB <- S[(sum(p_vec[0 : (k - 1)]) + 1) : sum(p_vec[0 : k]), 
						  (sum(p_vec[0 : (kp - 1)]) + 1) : sum(p_vec[0 : kp])]
			if(kp == k){
				SUB_ON <- mean(diag(SUB), na.rm = TRUE)
				SUB_OF <- (sum(SUB[upper.tri(SUB)], na.rm = TRUE) + sum(SUB[lower.tri(SUB)], na.rm = TRUE)) / (p_vec[k] ^ 2 - p_vec[k] - sum(is.na(SUB[upper.tri(SUB)])) - sum(is.na(SUB[lower.tri(SUB)])))
				A_temp[k, kp] <- SUB_ON - SUB_OF
				B_temp[k, kp] <- SUB_OF
			} else{
				B_temp[k, kp] <- mean(as.vector(SUB), na.rm = TRUE)
			}	
		}
	}
	A <- A_temp
	B <- (B_temp + t(B_temp)) / 2
	return(list(A = A, B = B))
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
	
		RES_OPTIM <- optim(par = runif(dim_gamma, gamma_ini_lb, gamma_ini_ub), 
							fn = NEG_LOG_LIKELIHOOD_FUN_FOR_GAMMA, 
							gr = NEG_SCORE_FUN_FOR_GAMMA, 
						method = "L-BFGS-B",
					  p_vector = p_vector, 
					         S = S, 
				    		 n = n,			  
			  SPLIT_UB_FOR_SUM = SPLIT_UB_FOR_SUM, 
	    BLOCK_HADAMARD_PRODUCT = BLOCK_HADAMARD_PRODUCT) 
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
	return(gamma_temp)
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
		DIVIDED_SEQUENCES <- c(rep(divided_by, R %/% divided_by), R - sum(rep(divided_by, R %/% divided_by))) 
		### sum = R 
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

SOLVING_GLMM_MODEL <- function(XT_mat, YT_mat){

	RES_GLM_REG <- lm(YT_mat ~ - 1 + XT_mat)
	SSCPE <- crossprod(YT_mat - RES_GLM_REG$fitted.values)
	SigmaHat <- SSCPE / (nrow(YT_mat) - nrow(coef(RES_GLM_REG)))
	### in our setting, GLMM assumes a diagonal covariance matrix
	return(diag(diag(SigmaHat)))
}

#############################
### Loading the data file ###
#############################

INPUT_ADDRESS <- OUTPUT_ADDRESS <- "..."

NMR_temp <- readr::read_csv(paste(INPUT_ADDRESS, "NMR.csv", sep = "/"))
clist_temp <- readr::read_csv(paste(INPUT_ADDRESS, "Clist_SICERS_E.csv", sep = "/"))
covariates_temp <- readRDS(paste(INPUT_ADDRESS, "covariatesv3.Rds", sep = "/"))
biomarkers_source <- readxl::read_excel(paste(INPUT_ADDRESS, "NMR_Biomarkers_249.xlsx", sep = "/"), sheet = "Source")
### The real dataset is available upon request from the authors.
varnames_170 <- as.numeric(colnames(clist_temp))
subject_id <- NMR_temp$eid
covariates <- covariates_temp[, c("Sex", "AgeIns2", "BMI", "AlcoholIntakeFrq", "heavysmk")]
### Sex: 0 female, 1 male 
### https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=31
### AgeIns2: numeric, in year
### https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21003
### BMI: numeric, in Kg/m2
### https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=21001
### AlcoholIntakeFrq: 1 daily, 2 three-four, ..., 6 never [decreasing]
### https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=1558	
### heavysmk: 0 not heavy, 1 heavy
### https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id=1249
covariates$AlcoholIntakeFrq <- 7 - covariates$AlcoholIntakeFrq
### AlcoholIntakeFrq: 1 never, 6 daily [increasing]
biomarkers <- NMR_temp[, - c(1, 251, 252)]
biomarkers_name_raw <- colnames(biomarkers)
biomarkers_name_all <- biomarkers_name_raw[varnames_170]
biomarkers_name_170 <- biomarkers_name_all[1 : 170]
biomarkers_name_rest <- biomarkers_name_all[- (1 : 170)]
### note that the biomarkers have been reordered using the community detection algorithm

###################################
### Removing the missing values ###
###################################

Y_mat_temp <- data.matrix(biomarkers)
### 4140 by 249, with colnames
X_mat_temp <- data.matrix(cbind(intercept = rep(1, 4140), covariates))
### 4140 by 6, first is intercept, with colnames
covariates_name_6 <- colnames(X_mat_temp)
### 6 by 1
row_num_without_NA <- which(rowSums(is.na(X_mat_temp)) == 0)
### rows having NA's = 4140 - 3984 = 156
subject_id_cl <- subject_id[row_num_without_NA] ### 3984 * 1
Y_mat_cl <- Y_mat_temp[row_num_without_NA, ]
### 3984 by 249, with colnames
X_mat_cl <- X_mat_temp[row_num_without_NA, ]
### 3984 by 6, first is intercept, with colnames
YT_mat <- Y_mat_cl[, biomarkers_name_170]
### 3984 by 170, with colnames
YT_rest_mat <- Y_mat_cl[, biomarkers_name_rest]
### 3984 by 79, with colnames

XT_mat <- X_mat_cl
### n = 3984 by p = 6, first is intercept, with colnames
Y_mat <- t(YT_mat)
### R = 170 by n = 3984
X_mat <- t(XT_mat)
### p = 6 by n = 3984, each column = p by 1 xi
### following a definition in VERSION 1.8.0 of paper

p <- ncol(XT_mat)		       ### 6: 1, Sex, AgeIns2, BMI, Ach, heavysmk
n <- nrow(XT_mat)			   ### 3984 	 
r_vec <- c(77, 47, 19, 11, 16) ### vector of features for G = 5
R_mat <- diag(r_vec) 
G <- length(r_vec)   		   ### number of groups, G = 5
R <- sum(r_vec)      		   ### total number of features 170
alpha <- 0.05

inv_mat <- matrix(0, p, p)
for(i in 1 : n){
	inv_mat <- inv_mat + X_mat[, i] %*% t(X_mat[, i])   ### p by 1
}
inv_x_sum_mat <- solve(unname(inv_mat))                 ### p by p

p_adj_method <- "BH"

##############################################################
### Creating the histograms for communities and singletons ###
##############################################################

R_mat_cl <- unname(cor(as.matrix(cbind(YT_mat, YT_rest_mat))))

s_vec <- c(r_vec, 79)
	
for(k in 1 : length(s_vec)){
	for(kp in 1 : length(s_vec)){
		SUB_temp <- R_mat_cl[(sum(s_vec[0 : (k - 1)]) + 1) : sum(s_vec[0 : k]), (sum(s_vec[0 : (kp - 1)]) + 1) : sum(s_vec[0 : kp])]
		if(k == kp){
			
			FILE_NAME_3 <- paste("Boxplot_communities_", k, "_and_", kp, ".pdf", sep = "")
			pdf(paste(OUTPUT_ADDRESS, FILE_NAME_3, sep = "/"), width = 12, height = 21, pointsize = 64)
			boxplot(as.vector(SUB_temp[upper.tri(SUB_temp)]), ylim = c(-1, 1))
			dev.off()
		} else if(k < kp){
			
			FILE_NAME_3 <- paste("Boxplot_communities_", k, "_and_", kp, ".pdf", sep = "")
			pdf(paste(OUTPUT_ADDRESS, FILE_NAME_3, sep = "/"), width = 12, height = 21, pointsize = 64)
			boxplot(as.vector(SUB_temp), ylim = c(-1, 1))
			dev.off()
		}	
	}
}



###################################################################
### Matching our community membership and the source membership ###
###################################################################

biomarkers_source_groups <- rep(0, length(biomarkers_name_all))
for(reps_bio in 1 : length(biomarkers_name_all)){
	for(reps_bio_pm in 1 : dim(biomarkers_source)[1]){
		if(biomarkers_name_all[reps_bio] == biomarkers_source$Description[reps_bio_pm]) {
			biomarkers_source_groups[reps_bio] <- biomarkers_source$Group[reps_bio_pm]
		}
	}
}
biomarkers_name_our_groups <- rep(c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "singletons"), c(77, 47, 19, 11, 16, 79)) ### 249 by 1
biomarkers_name_two_groups <- data.frame(cbind(biomarkers_name_all, biomarkers_name_our_groups, biomarkers_source_groups))
colnames(biomarkers_name_two_groups) <- c("Name", "Our Groups", "Source Groups")

xlsx::write.xlsx(
		  x = biomarkers_name_two_groups,
	   file = paste(OUTPUT_ADDRESS, "Evaluation of Associations between Plasma Metabolomics and Alcohol Intake in Older Adults.xlsx", sep = "/"),
  sheetName = "biomarkers_name_two_groups",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)
### compare group membership

##########################################
### Fitting the GLMM to 249 biomarkers ###
##########################################
Y_all_mat <- as.matrix(rbind(t(YT_mat), t(YT_rest_mat)))
### 249 by 3984
YT_all_mat <- t(Y_all_mat)
### 3984 by 249
R_all <- 249
BetaT_all_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% YT_all_mat 
### p by 249
Beta_all_OLS <- t(BetaT_all_OLS)							   
### 249 by p
E_all_mat <- Y_all_mat - Beta_all_OLS %*% X_mat         
### 249 by 3984
ET_all_mat <- t(E_all_mat) 
### 3984 by 249
S_all <- (n - 1) * cov(ET_all_mat) / n 					   
### 249 * 249
beta_all_OLS <- as.vector(BetaT_all_OLS)           
### by row of Beta_all_OLS

Sigma_all_EST_GLMM <- SOLVING_GLMM_MODEL(XT_mat, YT_all_mat)
COV_beta_all_EST_GLMM <- kronecker(Sigma_all_EST_GLMM, inv_x_sum_mat)

SE_Beta_all_GLMM <- matrix(sqrt(diag(COV_beta_all_EST_GLMM)), nrow = R_all, ncol = p, byrow = TRUE)
T_Beta_all_GLMM <- Beta_all_OLS / SE_Beta_all_GLMM
p_Beta_all_GLMM <- 2 * (1 - pt(abs(T_Beta_all_GLMM), df = n - 1))

colnames(Beta_all_OLS) <- colnames(SE_Beta_all_GLMM) <- colnames(T_Beta_all_GLMM) <- colnames(p_Beta_all_GLMM) <- covariates_name_6
rownames(Beta_all_OLS) <- rownames(SE_Beta_all_GLMM) <- rownames(T_Beta_all_GLMM) <- rownames(p_Beta_all_GLMM) <- c(biomarkers_name_170, biomarkers_name_rest)
### 249 by p

output_all_mat <- data.matrix(cbind(c(biomarkers_name_170, biomarkers_name_rest), Beta_all_OLS, p_Beta_all_GLMM))
colnames(output_all_mat) <- c("name", paste("EST", covariates_name_6, sep = "_"), paste("unadj_p", covariates_name_6, sep = "_"))

xlsx::write.xlsx(
		  x = output_all_mat,
	   file = paste(OUTPUT_ADDRESS, "Evaluation of Associations between Plasma Metabolomics and Alcohol Intake in Older Adults.xlsx", sep = "/"),
  sheetName = "GLMM for all",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

#############################################################
### Fitting the GLMM to 79 biomarkers outside communities ###
#############################################################

Y_rest_mat <- t(YT_rest_mat)
### 79 by 3984
R_rest <- 249 - R 
### 79
BetaT_rest_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% YT_rest_mat 
### p by 79
Beta_rest_OLS <- t(BetaT_rest_OLS)							   
### 79 by p
E_rest_mat <- Y_rest_mat - Beta_rest_OLS %*% X_mat         
### 79 by 3984
ET_rest_mat <- t(E_rest_mat) 
### 3984 by 79
S_rest <- (n - 1) * cov(ET_rest_mat) / n 					   
### 79 * 79
beta_rest_OLS <- as.vector(BetaT_rest_OLS)           
### by row of Beta_rest_OLS

Sigma_rest_EST_GLMM <- SOLVING_GLMM_MODEL(XT_mat, YT_rest_mat)
COV_beta_rest_EST_GLMM <- kronecker(Sigma_rest_EST_GLMM, inv_x_sum_mat)

SE_Beta_rest_GLMM <- matrix(sqrt(diag(COV_beta_rest_EST_GLMM)), nrow = R_rest, ncol = p, byrow = TRUE)
T_Beta_rest_GLMM <- Beta_rest_OLS / SE_Beta_rest_GLMM
p_Beta_rest_GLMM <- 2 * (1 - pt(abs(T_Beta_rest_GLMM), df = n - 1))

colnames(Beta_rest_OLS) <- colnames(SE_Beta_rest_GLMM) <- colnames(T_Beta_rest_GLMM) <- colnames(p_Beta_rest_GLMM) <- covariates_name_6
rownames(Beta_rest_OLS) <- rownames(SE_Beta_rest_GLMM) <- rownames(T_Beta_rest_GLMM) <- rownames(p_Beta_rest_GLMM) <- biomarkers_name_rest
### 79 by p

######################################################################
### Fitting the various model to 170 biomarkers within communities ###
######################################################################

### 1 Calculating the OLS for beta or B ###
	
BetaT_OLS <- solve(t(XT_mat) %*% XT_mat) %*% t(XT_mat) %*% YT_mat 
Beta_OLS <- t(BetaT_OLS)								 		   
### R * p
E_mat <- Y_mat - Beta_OLS %*% X_mat   				           
### R by n
ET_mat <- t(E_mat)   				       				   
### n by R
S <- (n - 1) * cov(ET_mat) / n 					     	   
### R by R
beta_OLS <- as.vector(BetaT_OLS)                     
### by row of B_OLS

### 2 Calculating the SE for beta or B ###
### 2-1 MAUD ###
L_mat <- diag(r_vec)
iters_max <- 100
gamma_ini_lb <- - 1
gamma_ini_ub <- 2 
gamma_EST_MAUD <- SOLVING_SCORE_FUN_FOR_GAMMA_OPTIM(
						   p_vector = r_vec, 
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
Sigma_EST_MAUD <- BLOCK_HADAMARD_PRODUCT(AS_EST_MAUD, BS_EST_MAUD, r_vec)
COV_beta_EST_MAUD <- kronecker(Sigma_EST_MAUD, inv_x_sum_mat)

SE_Beta_MAUD <- matrix(sqrt(diag(COV_beta_EST_MAUD)), nrow = R, ncol = p, byrow = TRUE)
T_Beta_MAUD <- Beta_OLS / SE_Beta_MAUD
p_Beta_MAUD <- 2 * (1 - pt(abs(T_Beta_MAUD), df = n - 1))

colnames(Beta_OLS) <- colnames(SE_Beta_MAUD) <- colnames(T_Beta_MAUD) <- colnames(p_Beta_MAUD) <- covariates_name_6
rownames(Beta_OLS) <- rownames(SE_Beta_MAUD) <- rownames(T_Beta_MAUD) <- rownames(p_Beta_MAUD) <- biomarkers_name_170
### R by p

### 2-2 GLMM ###
Sigma_EST_GLMM <- SOLVING_GLMM_MODEL(XT_mat, YT_mat)
COV_beta_EST_GLMM <- kronecker(Sigma_EST_GLMM, inv_x_sum_mat)

SE_Beta_GLMM <- matrix(sqrt(diag(COV_beta_EST_GLMM)), nrow = R, ncol = p, byrow = TRUE)
T_Beta_GLMM <- Beta_OLS / SE_Beta_GLMM
p_Beta_GLMM <- 2 * (1 - pt(abs(T_Beta_GLMM), df = n - 1))

colnames(SE_Beta_GLMM) <- colnames(T_Beta_GLMM) <- colnames(p_Beta_GLMM) <- covariates_name_6
rownames(SE_Beta_GLMM) <- rownames(T_Beta_GLMM) <- rownames(p_Beta_GLMM) <- biomarkers_name_170

### 3 Generating the results for all biomarkers ###

Beta_all_sep_OLS <- rbind(Beta_OLS, Beta_rest_OLS)
p_Beta_all_GLMM_GLMM <- rbind(p_Beta_GLMM, p_Beta_rest_GLMM)
p_Beta_all_MAUD_GLMM <- rbind(p_Beta_MAUD, p_Beta_rest_GLMM)
SE_Beta_all_MAUD_GLMM <- rbind(SE_Beta_MAUD, SE_Beta_rest_GLMM)


output_GLMM_GLMM_mat <- data.matrix(cbind(c(biomarkers_name_170, biomarkers_name_rest), Beta_all_sep_OLS, p_Beta_all_GLMM_GLMM))
output_MAUD_GLMM_mat <- data.matrix(cbind(c(biomarkers_name_170, biomarkers_name_rest), Beta_all_sep_OLS, p_Beta_all_MAUD_GLMM))

colnames(output_GLMM_GLMM_mat) <- colnames(output_MAUD_GLMM_mat) <- c("name", paste("EST", covariates_name_6, sep = "_"), paste("unadj_p", covariates_name_6, sep = "_"))

xlsx::write.xlsx(
		  x = output_GLMM_GLMM_mat,
	   file = paste(OUTPUT_ADDRESS, "Investigation of the Effect of Alcohol Intake on Plasma Metabolomics.xlsx", sep = "/"),
  sheetName = "GLMM_GLMM",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

xlsx::write.xlsx(
		  x = output_MAUD_GLMM_mat,
	   file = paste(OUTPUT_ADDRESS, "Evaluation of Associations between Plasma Metabolomics and Alcohol Intake in Older Adults", sep = "/"),
  sheetName = "MAUD_GLMM",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

output_decision <- cbind(
c(biomarkers_name_170, biomarkers_name_rest), 
(p.adjust(p_Beta_all_GLMM[, 5], method = p_adj_method) < alpha) * 1, 
(p.adjust(p_Beta_all_GLMM_GLMM[, 5], method = p_adj_method) < alpha) * 1, 
(p.adjust(p_Beta_all_MAUD_GLMM[, 5], method = p_adj_method) < alpha) * 1
)

colnames(output_decision) <- c("name", "GLMM_for_all", "GLMM_GLMM", "MAUD_GLMM")

xlsx::write.xlsx(
		  x = output_decision,
	   file = paste(OUTPUT_ADDRESS, "Evaluation of Associations between Plasma Metabolomics and Alcohol Intake in Older Adults.xlsx", sep = "/"),
  sheetName = "Decision",
  col.names = TRUE,
  row.names = FALSE,
     append = TRUE,
     showNA = TRUE,
   password = NULL
)

#############################
### Creating forest plots ###
#############################

beta_SIG_ACH <- as.numeric(output_decision[, 4])
beta_OLS_ACH <- unname(Beta_all_sep_OLS[, 5])
beta_SE_ACH <- unname(SE_Beta_all_MAUD_GLMM[, 5])
beta_LB_ACH <- beta_OLS_ACH - 1.96 * beta_SE_ACH
beta_UB_ACH <- beta_OLS_ACH + 1.96 * beta_SE_ACH
beta_CL_ACH <- c(ifelse(beta_SIG_ACH[1:170] == 0, "grey", "black"),
				 ifelse(beta_SIG_ACH[-(1:170)] == 0, "grey", "blue"))
beta_HDL_ACH <- stringr::str_detect(c(biomarkers_name_170, biomarkers_name_rest), "HDL")
beta_CL2_ACH <- beta_CL_ACH
for(i in 1 : length(beta_HDL_ACH)){
	if(beta_HDL_ACH[i] & (beta_CL_ACH[i] != "grey") & (i != 57)) beta_CL2_ACH[i] <- "red"
}

FOREST_DF <- data.frame(
	Index = seq(249), ## This provides an order to the data
	Biomarker = c(biomarkers_name_170, biomarkers_name_rest), 
	Alcohol_Intake_Frq = beta_OLS_ACH,
	LL = beta_LB_ACH,
	UL = beta_UB_ACH,
	SIG = beta_SIG_ACH, 
	COL = beta_CL2_ACH
)

plot_within <- ggplot2::ggplot(
	FOREST_DF[1 : R, ], aes(y = Index, x = Alcohol_Intake_Frq)) +
	geom_point(shape = 18, size = 5, colour = beta_CL2_ACH[1 : R]) +  
	geom_errorbarh(aes(xmin = beta_LB_ACH[1 : R], xmax = beta_UB_ACH[1 : R]), height = 0.25, colour = beta_CL2_ACH[1:R]) +
	geom_vline(xintercept = 0, color = "grey", linetype = "longdash", cex = 1, alpha = 0.5) +
	geom_vline(xintercept =-0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	geom_vline(xintercept = 0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	scale_y_continuous(name = "", breaks = 1 : R, labels = FOREST_DF[1:R, ]$Biomarker, trans = "reverse") +
	xlab("Estimates of Coefficients (95% CI)") + 
	ylab(" ") + 
	theme_bw() +
	theme(panel.border = element_blank(),
			panel.background = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			axis.line = element_line(colour = "black"),
			axis.text.y = element_text(size = 12, colour = beta_CL2_ACH[1 : R]),
			axis.text.x.bottom = element_text(size = 12, colour = "black"),
			axis.title.x = element_text(size = 12, colour = "black"))
FILE_NAME_1 <- paste("Forest_Plot_Features_Within_Communities", ".svg", sep = "")
ggsave(paste(OUTPUT_ADDRESS, FILE_NAME_1, sep = "/"), plot = plot_within, device = "svg", width = 20, height = 35, units = "in")

plot_outside <- ggplot2::ggplot(
	FOREST_DF[- (1 : R), ], aes(y = Index, x = Alcohol_Intake_Frq)) +
	geom_point(shape = 18, size = 5, colour = beta_CL2_ACH[- (1 : R)]) +  
	geom_errorbarh(aes(xmin = beta_LB_ACH[- (1 : R)], xmax = beta_UB_ACH[- (1 : R)]), height = 0.25, colour = beta_CL2_ACH[-(1:R)]) +
	geom_vline(xintercept = 0, color = "grey", linetype = "longdash", cex = 1, alpha = 0.5) +
	geom_vline(xintercept =-0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	geom_vline(xintercept = 0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	scale_y_continuous(name = "", breaks = ((R + 1) : 249), labels = FOREST_DF[- (1 : R), ]$Biomarker, trans = "reverse") +
	xlab("Estimates of Coefficients (95% CI)") + 
	ylab(" ") + 
	theme_bw() +
	theme(panel.border = element_blank(),
			panel.background = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			axis.line = element_line(colour = "black"),
			axis.text.y = element_text(size = 12, colour = beta_CL2_ACH[- (1 : R)]),
			axis.text.x.bottom = element_text(size = 12, colour = "black"),
			axis.title.x = element_text(size = 12, colour = "black"))
FILE_NAME_2 <- paste("Forest_Plot_Features_Outside_Communities", ".svg", sep = "")
ggsave(paste(OUTPUT_ADDRESS, FILE_NAME_2, sep = "/"), plot = plot_outside, device = "svg", width = 20, height = 35, units = "in")

left_R <- 1 : (R / 2)
right_R <- (R / 2 + 1) : R

plot_within_Left <- ggplot2::ggplot(
	FOREST_DF[left_R, ], aes(y = Index, x = Alcohol_Intake_Frq)) +
	geom_point(shape = 18, size = 5, colour = beta_CL2_ACH[left_R]) +  
	geom_errorbarh(aes(xmin = beta_LB_ACH[left_R], xmax = beta_UB_ACH[left_R]), height = 0.25, colour = beta_CL2_ACH[left_R]) +
	geom_vline(xintercept = 0, color = "grey", linetype = "longdash", cex = 1, alpha = 0.5) +
	geom_vline(xintercept =-0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	geom_vline(xintercept = 0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	scale_y_continuous(name = "", breaks = left_R, labels = FOREST_DF[left_R, ]$Biomarker, trans = "reverse") +
	xlab("Estimates of Coefficients (95% CI)") + 
	ylab(" ") + 
	theme_bw() +
	theme(panel.border = element_blank(),
			panel.background = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			axis.line = element_line(colour = "black"),
			axis.text.y = element_text(size = 12, colour = beta_CL2_ACH[left_R]),
			axis.text.x.bottom = element_text(size = 12, colour = "black"),
			axis.title.x = element_text(size = 12, colour = "black"))
FILE_NAME_1 <- paste("Forest_Plot_Features_Within_Communities_Left", ".svg", sep = "")
ggsave(paste(OUTPUT_ADDRESS, FILE_NAME_1, sep = "/"), plot = plot_within_Left, device = "svg", width = 20, height = 35, units = "in")

plot_within_Right <- ggplot2::ggplot(
	FOREST_DF[right_R, ], aes(y = Index, x = Alcohol_Intake_Frq)) +
	geom_point(shape = 18, size = 5, colour = beta_CL2_ACH[right_R]) +  
	geom_errorbarh(aes(xmin = beta_LB_ACH[right_R], xmax = beta_UB_ACH[right_R]), height = 0.25, colour = beta_CL2_ACH[right_R]) +
	geom_vline(xintercept = 0, color = "grey", linetype = "longdash", cex = 1, alpha = 0.5) +
	geom_vline(xintercept =-0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	geom_vline(xintercept = 0.1, color = "grey", linetype = "dotted", cex = 1, alpha = 0.5) +
	scale_y_continuous(name = "", breaks = right_R, labels = FOREST_DF[right_R, ]$Biomarker, trans = "reverse") +
	xlab("Estimates of Coefficients (95% CI)") + 
	ylab(" ") + 
	theme_bw() +
	theme(panel.border = element_blank(),
			panel.background = element_blank(),
			panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(), 
			axis.line = element_line(colour = "black"),
			axis.text.y = element_text(size = 12, colour = beta_CL2_ACH[right_R]),
			axis.text.x.bottom = element_text(size = 12, colour = "black"),
			axis.title.x = element_text(size = 12, colour = "black"))
FILE_NAME_1 <- paste("Forest_Plot_Features_Within_Communities_Right", ".svg", sep = "")
ggsave(paste(OUTPUT_ADDRESS, FILE_NAME_1, sep = "/"), plot = plot_within_Right, device = "svg", width = 20, height = 35, units = "in")