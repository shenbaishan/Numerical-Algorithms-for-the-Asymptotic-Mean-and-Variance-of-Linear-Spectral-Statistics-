function result = compute_contour_integral_part1(R, c)
    lambda = eig(R);
    % 定义 a, b, 和 r
    a = min(lambda) * (1 - sqrt(c))^2 * (c < 1);
    b = max(lambda) * (1 + sqrt(c))^2;
    r = 0.001;
    N = 1000;  % 离散路径的点数
    theta = linspace(0, 2 * pi, N);
    z_path = r * exp(1i * theta);  % 离散化路径 C

    % 定义外层积分的核函数
    function outer_integrand = outer_integrand(theta)
        % z 点
        z = r * exp(1i * theta);
        
        % 内层积分 \int [t * s'(z)]^2 / [s(z) * (1 + t * s(z))^3] dH(t)
        inner_integral = compute_sum(R, z, c);
        
        % 外层积分核函数
        outer_integrand = log(z) * inner_integral  * 1i * r * exp(1i * theta);
    end

    % 计算外层积分（路径积分）
    result = -c / (2 * pi * 1i) * integral(@(theta) outer_integrand(theta), 0, 2 * pi, 'ArrayValued', true);
end
