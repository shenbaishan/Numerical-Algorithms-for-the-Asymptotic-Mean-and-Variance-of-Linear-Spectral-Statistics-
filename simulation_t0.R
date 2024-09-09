library(expm) # for matrix square root
library(mvtnorm) # for generating random variables

generate_sigma <- function(p, theta, eta, scenario) {
  Sigma <- matrix(0, nrow = p, ncol = p)
  
  if (scenario == 1) {
    # Scenario 1
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] <- ifelse(i == j, 1, theta^abs(i - j))
      }
    }
  } else if (scenario == 2) {
    # Scenario 2
    for (i in 1:p) {
      for (j in 1:p) {
        Sigma[i, j] <- ifelse(i == j, 1, eta)
      }
    }
  } else {
    stop("Invalid scenario selected")
  }
  
  return(Sigma)
}

# Function to generate an n x p matrix X with each row uniformly distributed on the unit sphere S^(p-1)
generate_X_on_sphere <- function(n, p) {
  X <- matrix(nrow = n, ncol = p)
  for (i in 1:n) {
    x <- rnorm(p) # Generate a p-dimensional normal vector
    X[i, ] <- x / sqrt(sum(x^2)) # Normalize to lie on the sphere
  }
  return(X)
}

# Function to generate y_j samples
generate_Y <- function(n, p, Gamma, v) {
  X <- generate_X_on_sphere(n, p)
  Y <- matrix(nrow = n, ncol = p)
  for (j in 1:n) {
    # Generate rho_j from the F distribution
    rho_j_squared <- p * rf(1, df1 = p, df2 = v)
    # Then take the square root
    rho_j <- sqrt(rho_j_squared)
    
    # Compute y_j using the spatial-sign transformation
    y_j <- rho_j * Gamma %*% X[j, ]
    
    # Scale y_j to have unit norm
    Y[j, ] <- sqrt(p) * y_j / sqrt(sum(y_j^2))
  }
  return(Y)
}
generate_samples <- function(n, p, theta, eta, scenario) {
  Sigma <- generate_sigma(p, theta, eta, scenario)
  Gamma <- sqrtm(Sigma)
  Y <- generate_Y(n, p, Gamma,v)
  return(Y)
}

system.time({
  # Set the parameters
  n <- 625 # Number of samples
  p <- 500 # Dimension
  theta <- 0.03
  eta <- 0.03 # Parameter for scenario 2
  v <- 5 # Degrees of freedom for the t distribution
  scenario <- 1 # Choose the scenario
  c_n <- p/n
  c_n1 <- p/(n-1)
  alpha <- 0.05 # Significance level
  num_tests <- 2000 # Number of tests
  
  set.seed(123) 
  
  generate_TL <- function(n, p, theta, eta, scenario) {
    Y <- generate_samples(n, p, theta, eta, scenario)
    R_n <- cor(Y) # 计算样本相关性矩阵
    eigenvalues <- eigen(R_n)$values # 计算特征值
    T_L <- sum(log(eigenvalues)) # 计算测试统计量T_L
    return(T_L)
  }
  
  
  # 进行模拟并计算所有T_L的均值
  test_statistics_T_L <- replicate(num_tests, {
    generate_TL(n, p, theta, eta, scenario)
  })
  # 计算均值
  mean_T_L <- mean(test_statistics_T_L)
  #计算方差
  var_T_L <- var(test_statistics_T_L)
  # 标准化均值
 
  
  # 输出均值和标准化后的均值
  cat("Mean of T_L:", mean_T_L, "\n")
  cat("VAR of T_L:",var_T_L , "\n")
  
}) 