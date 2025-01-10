function R_z_prime = compute_Rz_prime(R, s_z, s_z_prime)
    % 输入:
    % R: 矩阵 R
    % s_z: 标量参数 s(z)
    % s_z_prime: 标量参数 s'(z)
    % 输出:
    % R_z_prime: 计算结果 R'(z) = −R(z) ⋅ (s_z_prime ⋅ R) ⋅ R(z)

    % 调用 compute_Rz 计算 R(z)
    R_z = compute_Rz(R, s_z);
    
    % 计算 R'(z)
    R_z_prime = -R_z * (s_z_prime * R) * R_z;
end
