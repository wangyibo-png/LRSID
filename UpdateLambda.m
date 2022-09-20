function [Lambda_new] = UpdateLambda(Lambda, tau, beta, Expre)
%UPDATELAMBDA 此处显示有关此函数的摘要
%   此处显示详细说明
Lambda_new = Lambda - tau*beta*Expre;
end

