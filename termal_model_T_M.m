function [New] = termal_model_T_M (I)
global Old_T m_magnet dtt Room_temperature
for i = 1:8
    Current_T(i,:) = (termal_model (Old_T(i,:),I(i),m_magnet,dtt,Room_temperature));
end
New = max(max(Current_T));
% New = sum(max(Current_T));
end