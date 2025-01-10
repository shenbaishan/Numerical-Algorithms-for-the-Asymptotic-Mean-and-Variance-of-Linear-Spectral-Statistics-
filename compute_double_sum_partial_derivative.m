function result = compute_double_sum_partial_derivative(R, z, c)
    s_prime_z = s_prime(z, c);
    s_z = s_z_(z, c);
    R_z = compute_Rz(R, s_z);
    R_z_prime = compute_Rz_prime(R, s_z,  s_prime_z);
    result = 0;
    p = size(R, 1);
    % 遍历所有 k 和 l
    for k = 1:p
        for l = 1:p
            % 提取 r_kl
            r_kl = R(k, l);

            % 定义单位向量 e_k 和 e_l
            e_k = zeros(p, 1); e_k(k) = 1;
            e_l = zeros(p, 1); e_l(l) = 1;

            % 计算 e_k^T R(z) e_l
            term_kl = e_k' * R_z * e_l;

            % 对 z 的导数部分
            partial_kl = 2 * term_kl * (e_k' * R_z_prime * e_l);

            % 累加结果
            result = result + r_kl^2 * partial_kl;
        end
    end
end