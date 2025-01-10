
n = 500;
p = 50;
theta_ = 0;
eta = 0.00;
v = 5;
scenario = 2;
c = p/n;
c1 = p / (n - 1);
scenario = 2;
beta_x = 0;
% Generate Sigma matrix using the generate_sigma function
Sigma = generate_sigma(p, theta_, eta, scenario);
Sigma_inv_sqrt = diag(1 ./ sqrt(diag(Sigma)));
R = Sigma_inv_sqrt * Sigma * Sigma_inv_sqrt;
Gamma = sqrtm(Sigma);

% 计算 Sigma 的对角线元素平方根的逆
diagSigmaInvSqrt = diag(1 ./ sqrt(diag(Sigma)));

% 定义单位矩阵 I
I = eye(p);

% 计算矩阵 G
G = diagSigmaInvSqrt * Gamma;

result = compute_contour_integral_part1(R, c) + compute_contour_integral_part2(R, c) + compute_contour_integral_part3(R, c);
disp('函数结果:')
disp(result);