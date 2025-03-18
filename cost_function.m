% Cost function
% Minimize the power bought from the DNO
function J = cost_function(X, U, e, data, Nc, M, beta_c, beta_d, wt, pv, D)
    J = 0;
    for k = 1:Nc
        for m = 1:M
            J = J + U(k, 6 + 8*(m-1));
        end
    end
end