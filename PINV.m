function I = PINV(Wrench, A, n, I_MaxHigh,I_MaxLow, I_pow)
    Is = pinv(A)*Wrench;
    [Gamma_amp, arg] = ScaleAmp(I_MaxHigh,I_MaxLow,Is);
    Gamma_pow = ScalePow(zeros(length(I_MaxHigh),1), Is, I_pow);
    Gamma = min([Gamma_amp,Gamma_pow,1]);
    if Gamma_pow <=1
%         disp('Power supply saturation')
%         Gamma_pow
    end
    if Gamma_amp ~=1
%         disp('Amp saturation')
%         Gamma_amp
    end
    
    if (Gamma_amp <1) || (Gamma_pow <=1) || (n == 0)
       I = Gamma * Is; 
    else
        I = Is;
    end
end