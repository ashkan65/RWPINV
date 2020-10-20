function I = RWPINV(Wrench, A, n, I_MaxHigh,I_MaxLow, I_pow,I_old, Avail_coils,W)
%     disp('------------ a call -------------')
act = A(:,Avail_coils);
    w = zeros(length(Avail_coils));
    for i = 1:length(Avail_coils)
        w(i,i) = W(Avail_coils(i),Avail_coils(i));
    end
%     I_temp = WPINV3(Wrench,act,w);
    I_temp = WPINV(Wrench,act,w);

    
%     I_temp  = pinv(act) * Wrench;
    I = zeros(8,1);
    for i = 1 : length(Avail_coils)
        I(Avail_coils(i)) = I_temp(i);
    end
    [Gamma_amp,arg] = ScaleAmp(I_MaxHigh - I_old, I_MaxLow - I_old, I);
%      if Gamma_amp ~=1
%         disp('Amp saturation')
%         Gamma_amp;
%     end
    if (Gamma_amp<1)
        for i = 1:length(arg)
            Avail_coils(find(Avail_coils == arg(i))) = [];
        end
    end
    I = Gamma_amp * I;
    Gamma_pow = ScalePow(I_old, I, I_pow); % This does the bysection!
%     if Gamma_pow <=1
%         disp('Power supply saturation')
%         Gamma_pow;
%     end
    if (Gamma_amp ==1) || (Gamma_pow <=1) || (n == 0)
%        Gamma_pow * I_temp;
%        I_current
        if (Gamma_pow <=1)
            I = Gamma_pow * I;
%         else
%             I =  I_temp;
        end
    else
        Wrench =  Wrench - A * I;
        I_old = I_old + I;
        n = n-1;
        I_r = RWPINV(Wrench, A, n, I_MaxHigh, I_MaxLow, I_pow, I_old, Avail_coils,W);
        I = I_r + I;        
    end
end