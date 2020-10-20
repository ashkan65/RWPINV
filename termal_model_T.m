function [Current_T] = termal_model_T (I,Room_temperature)
global Old_T m_magnet dtt
for i = 1:length(I)
    Current_T(i,:) = (termal_model (Old_T(i,:),I(i),m_magnet,dtt, Room_temperature));
end
% New = max(New_T(i));
end