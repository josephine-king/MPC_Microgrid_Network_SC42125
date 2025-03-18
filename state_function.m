% State function
function z = state_function(x, u, Nc, M, beta_c, beta_d, cap, wt, pv, D)
    z = x;
    for m = 1:M
        z(m) = z(m) + beta_c*u(1 + 8*(m-1)) - beta_d*u(2 + 8*(m-1));
    end
end
