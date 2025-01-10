function result = compute_sum(R, z, c)
    % 计算公式 \sum_{i=1}^p [lambda_i * s'(z)]^2 / (m(z) * [1 + lambda_i * s(z)]^3)
   
    lambda = eig(R);
    s_prime_z = s_prime(z, c);
    s_z = s_z_(z, c);
    % 逐项计算公式
    result = sum((lambda .* s_prime_z).^2 ./ (s_z .* (1 + lambda .* s_z).^3));
end