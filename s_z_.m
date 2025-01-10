function result = s_z_(z, c)
    % s(z) 的计算公式
    result = (sqrt((z - 1 - c)^2 - 4 * c) + c - z - 1) / (2 * z);
    if imag(result) < 0
    result = -result; % 如果虚部为负，则选择另一个分支
end


