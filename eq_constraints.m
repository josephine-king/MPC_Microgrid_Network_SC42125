% Equality constraints
function ceq = eq_constraints(X, U, data, Nc, M, beta_c, beta_d, wt, pv, D)

    ceq = [];
    % Loop over all the time indices in the horizon
    for k = 1:size(X,1)  

        % Get the inputs from the timestep
        Uk = U(k,:);

        % Energy balance
        ubal_lhs = zeros(M,1);
        ubal_rhs = zeros(M,1);
        for m = 1:M
            ubal_lhs(m) = wt(m) + pv(m) - D(m);
            energy_sold = sum(Uk(3 + 8*(m-1) : 5 + 8*(m-1)));
            energy_purchased = sum(Uk(6 + 8*(m-1) : 8 + 8*(m-1)));
            u_char = Uk(1);
            u_dischar = Uk(2);
            ubal_rhs(m) = energy_sold - energy_purchased + beta_c*u_char - beta_d*u_dischar;
        end 
    
        % Energy sold from one grid must match energy bought from other grid
        ums_12_ump_21 = Uk(4)  - Uk(15);
        ums_13_ump_31 = Uk(5)  - Uk(23);
        ums_23_ump_32 = Uk(13) - Uk(24);

        ums_21_ump_12 = Uk(12) - Uk(7);
        ums_31_ump_13 = Uk(20) - Uk(8);
        ums_32_ump_23 = Uk(21) - Uk(16);

        ceq = [ceq; ubal_lhs(1) - ubal_rhs(1); ubal_lhs(2) - ubal_rhs(2); ubal_lhs(3) - ubal_rhs(3)];
        ceq = [ceq; ums_12_ump_21; ums_13_ump_31; ums_23_ump_32; ums_21_ump_12; ums_31_ump_13; ums_32_ump_23];
    end

end