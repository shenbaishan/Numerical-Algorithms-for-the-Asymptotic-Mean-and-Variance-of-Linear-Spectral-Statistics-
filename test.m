 R = [4, 2, 1; 2, 5, 3; 1, 3, 6];

    % 参数
    z = 1 + 1i;  % 示例复数参数
    c = 0.5;     % 比例参数

    % 调用函数
    result =  compute_contour_integral_part3(R, c);

    % 显示结果
    disp('compute_double_sum_partial_derivative 函数结果:');
    disp(result);