function Sigma = generate_sigma(p, theta_, eta, scenario)
    % Initialize Sigma matrix
    Sigma = zeros(p, p);

    if scenario == 1
        % Scenario 1
        for i = 1:p
            for j = 1:p
                if i == j
                    Sigma(i, j) = 1;
                else
                    Sigma(i, j) = theta_^abs(i - j);
                end
            end
        end
    elseif scenario == 2
        % Scenario 2
        for i = 1:p
            for j = 1:p
                if i == j
                    Sigma(i, j) = 1;
                else
                    Sigma(i, j) = eta;
                end
            end
        end
    else
        error('Invalid scenario selected');
    end
end


