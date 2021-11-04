function [Gamma, arg] = ScaleAmp(I_MaxHigh,I_MaxLow,Is)
    Gamma = 1;
    arg = [];
%     I_MaxHigh
%     I_MaxLow
%     Is
    for i = 1 :length(Is)
        if (Is(i)>0) && I_MaxHigh(i)/Is(i)<Gamma
            Gamma = I_MaxHigh(i)/Is(i);
            arg = [i];
        end
        if (Is(i)<0) && I_MaxLow(i)/Is(i)<Gamma 
            Gamma = I_MaxLow(i)/Is(i);
            arg = [i];
        end
    end
end
