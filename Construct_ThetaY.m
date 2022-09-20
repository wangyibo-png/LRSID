function [Theta_y] = Construct_ThetaY(A, C, K, s)
%CONSTRUCT_THETAY 此处显示有关此函数的摘要
%   此处显示详细说明
Theta_y = [];
for i = 0:s-2
    Theta_y = [Theta_y;C*A^i*K];
end
end

