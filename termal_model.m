function [New_T] = termal_model (Old_T,I,magnet,dtt,Room_temperature)
S = length(Old_T);
room = Room_temperature;
% l = magnet.OR-magnet.IR/magnet.Wire.d;
% R =  magnet.IR : magnet.Wire.d : magnet.OR;
% R = [magnet.IR , R, magnet.OR];

R =  magnet.IR : magnet.Wire.d : 4*magnet.OR;
R = [magnet.IR , R, 4*magnet.OR];

T = room*ones(1,S + 2);
T(2 : (S + 1)) = Old_T;
% alpha = 3.9e-03;
% beta = 35;
% gamma = 7.15e-2;

alpha = 1.9e-03;
beta = 65;
gamma = 7.15e-2;
gamma = 18.15e-2;

% gamma * I^2 %* R(1) * 2* pi * 0.0144 * (1 + alpha * T(1)) 

for i = 2 : S+1  
    New_T(i-1) = dtt * ( ((T(i+1)-T(i))*(R(i)+R(i+1))) +  ((T(i-1)-T(i))*(R(i)+R(i-1))) + gamma * I^2 * R(i) * 2* pi * 0.0144 * (1 + alpha * T(i)) )/ (beta *  (R(i) * 2* pi * 0.0144));
end

New_T = Old_T + New_T;



