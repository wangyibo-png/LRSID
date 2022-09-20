function [Theta_u] = Construct_ThetaU(A, B, C, D, s)
%CONSTRUCT_THETAU 此处显示有关此函数的摘要
%   此处显示详细说明
Theta_u = [D];
for i = 0:s-2
    Theta_u = [Theta_u;C*A^i*B];
end
end