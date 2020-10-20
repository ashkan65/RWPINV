function [e] = Min_cost(I)
    global des act
    out_put = act*I;
%     vecnorm(out_put-des,1)
%     e = vecnorm(out_put-des,1) + 0.0001*vecnorm(I,1)/(0.0001 + 500*vecnorm(out_put-des,1)^2);
    er =out_put-des; 
    e = dot(er,er) + 0.00001*dot(I,I)/(1+dot(er,er)^1);
end