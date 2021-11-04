function Gamma = ScalePow(I_prev, I_new, I_pow)
%     Gamma = I_pow/sum(abs(I));
    error = 100;
    lower = 1/(L1(I_prev + I_new)/I_pow);
    upper = 1/(L1(I_new)/(I_pow - L1(I_prev)));
    while (abs(error)>0.1)
        mid = 0.5*(lower + upper);
        M = I_prev + mid*I_new;
        error = I_pow - L1(M);
        if (error < 0)
            lower = mid;
        else
            upper = mid;
        end
        if (abs(lower - upper) <0.001)
           error = 0;
        end
    end
    Gamma = mid;
end