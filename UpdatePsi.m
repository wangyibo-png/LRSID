function [Psi] = UpdatePsi(Hs_yhat, Ts_thetau, Us, Ts_thetay, Ys, Lambda, beta, n)
%UPDATEPSI 此处显示有关此函数的摘要
%   此处显示详细说明
Temp = Hs_yhat - Ts_thetau * Us - Ts_thetay * Ys - Lambda/beta;
[U, Sigma, V] = svd(Temp);
Sigma((n+1):end, (n+1):end) = 0;
Psi = U*Sigma*V';
end

