######################################
### COMMUNITY DETECTION ALGORITHMS ###
######################################

############################
### Loading the packages ###
############################

REQUIRED_PACKAGES <- c(
	"mvnfast",		
	### mvnfast::rmvn(), fast generate multivariate normal samples
	"matlib", 
	### matlib::symMat(), create a symmetric matrix from a vector
    "ICONS"
	### ICONS::dense() and ICONS::plotMatrix(), extract subnetwork and heatmap plot
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

BLOCK_HADAMARD_PRODUCT <- function(A, B, p_vec){
  K <- length(p_vec)
  COL_temp <- c()
  for(k in 1 : K){
    ROW_temp <- c()
    for(kp in 1 : K){
      if(k == kp){
        ROW_temp <- cbind(ROW_temp, A[k, k] * diag(p_vec[k]) + B[k, k] * matrix(1, p_vec[k], p_vec[k]))
      } else{
        ROW_temp <- cbind(ROW_temp, B[k, kp] * matrix(1, p_vec[k], p_vec[kp]))
      }
    }
    COL_temp <- rbind(COL_temp, ROW_temp)
  }
  M <- COL_temp
  return(M)
}

n <- 50
### sample size
l_vec <- c(30, 30, 40)
### the known partition-size vector
l0 <- 100
### the number of singletons
R <- sum(l_vec) 
S <- R + l0        		

A_0 <- diag(c(0.1, 0.2, 0.3))			
B_0 <- symMat(c(0.9, - 0.430920, -0.353572, 0.8, - 0.333474, 0.7))			
Sigma_0 <- matrix(0, S, S)
Sigma_0[1 : R, 1 : R] <- BLOCK_HADAMARD_PRODUCT(A_0, B_0, l_vec)
Sigma_0[(R + 1) : S, (R + 1) : S] <- diag(l0)

set.seed(20251)
data_mat <- mvrnorm(n = n, mu = rep(0, S), Sigma = Sigma_0)
### n by S
var_name <- paste("Var", seq(S), sep = "_")
colnames(data_mat) <- var_name
set.seed(20252)
name_order <- sample.int(n = S, size = S, replace = FALSE)
data_perm_mat <- data_mat[, name_order]
 
S_mat <- cov(data_mat)
### S by S
S_perm_mat <- cov(data_perm_mat)
### S by S 

### applying k-mean functions to data ###
set.seed(20253)
result_kmeans <- kmeans(x = t(data_perm_mat), centers = 4)
table(result_kmeans$cluster)
data_kmeans <- cbind(
	data_perm_mat[, result_kmeans$cluster == 1],
	data_perm_mat[, result_kmeans$cluster == 2],
	data_perm_mat[, result_kmeans$cluster == 4],
	data_perm_mat[, result_kmeans$cluster == 3])
S_mat_est_kmeans <- cov(data_kmeans)
	
### applying ICON functions to data ### 
set.seed(20254)
result_ICONS <- dense(S_perm_mat) 
S_mat_est_ICONS <- result_ICONS$W_dense
### note that the diagonals are zero
diag(S_mat_est_ICONS) <- diag(S_perm_mat)[result_ICONS$Clist]

OUTPUT_ADDRESS <- "..."
### the default output address

plotMatrix(z = Sigma_0, cex.axis = 1.3, cex.lab = 1.3, save.image = T, filepath = paste(OUTPUT_ADDRESS, "SM_Sigma_0.pdf", sep = ""), format = "pdf")
### heatmap of true covariance matrix based on covariates 1 to S
plotMatrix(z = S_mat, cex.axis = 1.3, cex.lab = 1.3, save.image = T, filepath = paste(OUTPUT_ADDRESS, "SM_S_mat.pdf", sep = ""), format = "pdf")
### heatmap of sample covariance matrix based on covariates 1 to S
plotMatrix(z = S_perm_mat, cex.axis = 1.3, cex.lab = 1.3, save.image = T, filepath = paste(OUTPUT_ADDRESS, "SM_S_perm_mat.pdf", sep = ""), format = "pdf")
### heatmap of sample covariance matrix based permuted covariates
plotMatrix(z = S_mat_est_kmeans, cex.axis = 1.3, cex.lab = 1.3, save.image = T, filepath = paste(OUTPUT_ADDRESS, "SM_S_mat_est_kmeans.pdf", sep = ""), format = "pdf")
### heatmap of sample covariance matrix based permuted covariates
plotMatrix(z = S_mat_est_ICONS, cex.axis = 1.3, cex.lab = 1.3, save.image = T, filepath = paste(OUTPUT_ADDRESS, "SM_S_mat_est_ICONS.pdf", sep = ""), format = "pdf")
### heatmap of sample covariance matrix based permuted covariates

