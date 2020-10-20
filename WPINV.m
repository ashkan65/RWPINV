function [I] =WPINV(des,act,W)
% global W
[~,~,V] = svd(act);
IS = pinv(act)*des;
% C = -(V(:,5:8)'*W*V(:,5:8 ))%*V(5:8 ,:)'*W*IS;
% V(:,6:8)';
if length(IS)<6
    disp('NO Answer');
end
C = -pinv(V(:,6:length(IS))'*W*V(: ,6:length(IS)))*V(:,6:length(IS))'*W*IS;
I = IS + V(:,6:length(IS))*C;
end