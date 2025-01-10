function result = compute_sum_partial_derivative(R, z, c)
    s_prime_z = s_prime(z, c);
    s_z = s_z_(z, c);
    result = 0;
    p = size(R, 1);
    R_z = compute_Rz(R, s_z);
    % 遍历每个 k
    for k = 1:p
        % 单位向量 e_k
        e_k = zeros(p, 1);
        e_k(k) = 1;

        % 计算 e_k^T R(z) R e_k
        term_k = e_k' * R_z * R * e_k;

        % 对 z 的导数部分
        partial_k = s_prime_z * (term_k^2) ...
                    + 2 * s_z * term_k * (e_k' * (compute_Rz_prime(R, s_z,  s_prime_z)) * R * e_k);

        % 累加
        result = result + partial_k;
    end
end