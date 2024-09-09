library(minpack.lm)
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


#p <- 2
#theta <- 0.2
#eta <- 2
#scenario <- 1
#Sigma <- generate_sigma(p, theta, eta, scenario)
m_zfun <- function(m, lambda, p, n, x) {
  mr <- m[1]
  mi <- m[2]
  f <- numeric(2)
  f[1] <- -mr / (mr^2 + mi^2)
  f[2] <- mi / (mr^2 + mi^2)
  for (i in 1:p) {
    f[1] <- f[1] + 1/n * (lambda[i] + lambda[i]^2 * mr) / ((1 + lambda[i] * mr)^2 + lambda[i]^2 * mi^2)
    f[2] <- f[2] + 1/n * (-lambda[i]^2 * mi) / ((1 + lambda[i] * mr)^2 + lambda[i]^2 * mi^2)
  }
  f[1] <- f[1] - x
  f[2] <- f[2] - 1/100000
  return(f)
}
fun_optim <- function(m, lambda, p, n, x) {
  f <- m_zfun(m, lambda, p, n, x)
  return(sum(f^2))  # 返回残差的平方和
}


  
  # 初始化参数
  n <- 1000 # Number of samples
  p <- 100 # Dimension
  theta <- 0.05 # Parameter for scenario 1
  
  
  eta <- 0.03 # Parameter for scenario 2
  v <- 5 # Degrees of freedom for the t distribution
  scenario <- 1 # Choose the scenario
  y <- p/n
  y1 <- p/(n-1)
  scenario <- 1
  
  
  
  Sigma <- generate_sigma(p, theta, eta, scenario)
  
  
  Sigma_inv_sqrt <- solve(diag(sqrt(diag(Sigma))))
  R <- Sigma_inv_sqrt %*% Sigma %*% Sigma_inv_sqrt
  lambda <- eigen(R)$values
  I <- diag(1, p)
  
  
  a <- min(lambda) * (1 - sqrt(y))^2 * as.numeric(y < 1)
  b <- max(lambda) * (1 + sqrt(y))^2
  interval <- 0.001
  x <- seq(a, b, by = interval)
  initial_value <- c(0, 0.1)
  M <- matrix(ncol = 2, nrow = 0)  # 存储 m 值
  tau <- 2
  EX <- 0
  
  
  
  optimized_params_matrix <- matrix(0, nrow = length(x), ncol = 2)
  for (j in 1:length(x)) {
    result <- optim(par = initial_value, 
                    fn = fun_optim,  # 使用修改后的函数
                    lambda = lambda, 
                    p = p, 
                    n = n, 
                    x = x[j], 
                    method = "L-BFGS-B", 
                    lower = c(-Inf, 0),  # 确保mi的下界是0
                    upper = c(Inf, Inf), 
                    control = list(maxit = 10000)) # 可以通过control参数调整优化过程，如最大迭代次数
    optimized_params <-result$par
    optimized_params_matrix[j,] <- result$par
  }
  
  
  M <-optimized_params_matrix
  
  
  for (j in 1:nrow(M)) {
    # 从M矩阵中创建复数
    m_complex <- complex(real = M[j, 1], imaginary = M[j, 2])
    R_x <- solve(I + m_complex * R)
    
    # 第一部分
    sum1 <- sum(lambda^2 * m_complex^2 / (1 + lambda * m_complex)^2)
    part1 <-  (1 / x[j]) * Im(log(1 - y / p * sum1))/2
    
    
    part2 <- -1 / x[j] * Im(1 / p * sum1)/2
    
    # 第三部分
    # 计算 R_x R 的乘积
    product_R_x_R <- R_x %*% R
    
    # 提取对角元素
    diagonal_elements <- diag(product_R_x_R)
    
    # 对提取的对角元素平方，然后每个元素乘以 m_complex
    squared_and_multiplied_by_s <- (diagonal_elements^2) * m_complex
    
    # 计算最终的和
    sum3 <- sum(squared_and_multiplied_by_s)
    
    # 对 sum3 取虚部
    part3 <- -2 / n * (1 / x[j]) * Im(sum3)
    
    
    # 第四部分
    R_squared <- R^2  # R 的每个元素的平方
    R_x_squared <- R_x^2  # R_x 的每个元素的平方
    
    # 计算元素级乘积并求和
    sum4 <- sum(R_squared * R_x_squared)
    part4 <- -1 / n * (1 / x[j]) * Im(sum4)
    
    # 将四部分累加到EX
    EX <- EX + part1 + y*(tau - 2)*part2 + part3 + part4
  }
  EX <- (b - a) / (length(x) * pi) * EX
  
  
  Fch <- numeric(length(x))
  Fch <- M[, 2] / y / pi
  
  # 计算偏差
  err <- (x[2] - x[1]) * sum(Fch) - 1
  err <- err / (1 / y / pi * (x[2] - x[1]))
  
  j <- 1
  while (err - M[j, 2] > 0) {
    err <- err - M[j, 2]
    M[j, 2] <- 0
    j <- j + 1
  }
  
  M[j, 2] <- M[j, 2] - err
  
  
  FG <- sum((log(x)) * M[, 2]) / y / pi * (b - a) / length(x)
  
  M_L <- p*FG+EX
  
  cat("Mean of T_L:", M_L, "\n")
  
  
  VAR <- 0
  # 对于VAR1，只考虑s != t的情况
  combinations <- expand.grid(s = 1:(length(x) - 1), t = 1:(length(x) - 1))
  combinations <- subset(combinations, s != t)
  
  #VAR1的计算
  denominators <- with(combinations, (M[s, 1] - M[t, 1])^2 + (M[s, 2] - M[t, 2])^2)
  ln_inputs <- 1 + 4 * M[combinations$s, 2] * M[combinations$t, 2] / denominators
  valid_indices <- which(denominators > 0 & ln_inputs > 0)
  VAR1 <- sum((1 / x[combinations$s[valid_indices]]) * (1 / x[combinations$t[valid_indices]]) * log(ln_inputs[valid_indices]), na.rm = TRUE)
  VAR1 <- 0.5 * VAR1
  
  # 对于VAR2, VAR3, VAR4，考虑所有s和t的组合，包括s=t
  combinations <- expand.grid(s = 1:(length(x) - 1), t = 1:(length(x) - 1))
  n_combinations <- nrow(combinations)
  VAR2 <- VAR3 <- VAR4 <- 0
  
  # 预先计算R^2避免循环中重复计算
  R_squared <- R^2
  
  for (i in 1:n_combinations) {
    s <- combinations$s[i]
    t <- combinations$t[i]
    
    m_complex_s <- complex(real = M[s, 1], imaginary = M[s, 2])
    m_complex_t <- complex(real = M[t, 1], imaginary = M[t, 2])
    m_complex_s_ <- complex(real = M[s, 1], imaginary = -M[s, 2])
    m_complex_t_ <- complex(real = M[t, 1], imaginary = -M[t, 2])
    
    R_x <- solve(I + m_complex_s * R)
    R_y <- solve(I + m_complex_t * R)
    R_x_ <- solve(I + m_complex_s_ * R)
    R_y_ <- solve(I + m_complex_t_ * R)
    
    Rx_diag <- diag(R_x)
    Ry_diag <- diag(R_y)
    Rx_diag_ <- diag(R_x_)
    Ry_diag_ <- diag(R_y_)
    
    Ry_diag_VAR3 <- Ry_diag - 1
    Rx_diag_VAR4 <- Rx_diag - 1
    
    diag_product_VAR2 <- outer(Rx_diag, Ry_diag)
    diag_product_VAR2_ <- outer(Rx_diag_, Ry_diag)
    
    sum_result_VAR2 <- sum(Re(R_squared * diag_product_VAR2))-sum(Re(R_squared * diag_product_VAR2_))
    
    sum_k_VAR3 <- sum(Re(Rx_diag * (Ry_diag_VAR3 / m_complex_t)))-sum(Re(Rx_diag_ * (Ry_diag_VAR3 / m_complex_t)))
    sum_k_VAR4 <- sum(Re((Rx_diag_VAR4 / m_complex_s) * Ry_diag))-sum(Re((Rx_diag_VAR4 / m_complex_s) * Ry_diag_))
    
    VAR2 <- VAR2 + 1 / n * (1 / x[s]) * (1 /  x[t]) * sum_result_VAR2
    VAR3 <- VAR3 + 1 / n * (1 / x[s]) * (1 /  x[t]) * sum_k_VAR3
    VAR4 <- VAR4 + 1 / n * (1 /  x[s]) * (1 / x[t]) * sum_k_VAR4
  }
  
  # 最终VAR的计算
  VAR <- (VAR1 - VAR2 - VAR3 - VAR4) / pi^2 * (b - a)^2 / (length(x)^2 - length(x))
  
  
  

cat("VAR of T_L:",VAR , "\n")
