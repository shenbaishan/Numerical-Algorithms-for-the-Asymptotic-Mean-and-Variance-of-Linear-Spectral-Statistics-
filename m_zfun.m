function f = m_zfun(m, lambda, p, n, x)
    % 提取实部和虚部
    mr = m(1);
    mi = m(2);
    
    % 初始化结果
    f = zeros(2, 1);
    f(1) = -mr / (mr^2 + mi^2);
    f(2) = mi / (mr^2 + mi^2);
    
    % 向量化计算
    denom = (1 + lambda .* mr).^2 + (lambda .* mi).^2;
    f(1) = f(1) + sum(1/n * (lambda + lambda.^2 * mr) ./ denom);
    f(2) = f(2) + sum(1/n * (-lambda.^2 * mi) ./ denom);
    
    % 调整结果
    f(1) = f(1) - x;
    f(2) = f(2) - 1/100000;
end

function f_value = fun_optim(m, lambda, p, n, x)
    % 计算残差向量
    f = m_zfun(m, lambda, p, n, x);
    % 返回残差平方和
    f_value = sum(f.^2);
end
