function c = ineq_constraints(X, U, e, data, Nc, M, beta_c, beta_d, cap, wt, pv, D)
    c = [];

    % The stored energy must be greater than or equal to 0 
    for k = 1:size(X,1)  
        Xk = X(k,:)';
        Uk = U(k,:)';
        c = [c; -Xk; Xk - cap; -Uk];
    end

end