function result = s_prime(z, c)
    % s'(z) 的计算公式
    term1 = - (c - z + sqrt((z - 1 - c)^2 - 4 * c) - 1) / (2 * z^2);
    term2 = - ((2 * c - 2 * z + 2) / (2 * sqrt((z - 1 - c)^2 - 4 * c)) + 1) / (2 * z);
    result = term1 + term2;
end
