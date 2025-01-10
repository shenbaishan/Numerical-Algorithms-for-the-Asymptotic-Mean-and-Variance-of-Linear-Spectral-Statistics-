

function R_z = compute_Rz(R, s_z)
    % 输入:
    % R: 矩阵 R
    % s_z: 标量参数 s(z)
    % 输出:
    % R_z: 计算结果 R(z) = (I_p + s_z * R)^(-1)

    % 获取矩阵 R 的尺寸
    [p, ~] = size(R);
    
    % 定义单位矩阵 I_p
    I_p = eye(p);
    
    % 计算 R(z)
    R_z = inv(I_p + s_z * R);
end


