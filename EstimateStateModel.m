function [A, B, C, D, K] = EstimateStateModel(Hs_yhat, Ts_thetau, Ts_thetay, Us, Ys, n, m, p, s, N, U, Y)
%ESTIMATESTATEMODEL 此处显示有关此函数的摘要
%   此处显示详细说明
tempM = Hs_yhat - Ts_thetau*Us - Ts_thetay*Ys;
[U_svd, ~, ~] = svd(tempM);
Os = U_svd(:, 1:n);
C = Os(1:p, :);
A = pinv(Os(1:((s-1)*p), :)) * Os((p+1):end, :);
K = pinv(Os(1:((s-1)*p), :)) * Ts_thetay((p+1):p*s, 1:p);
A = A + K*C;
tempM1 = [];
tempM2 = [];
tempM3 = [];
for i = 1:N
    tempM1 = [tempM1; C*A^(i-1)];
    tempM3 = [tempM3; kron(U((i-1)*m+1:i*m, :)', eye(p))];
    tempM = zeros(p, m*n);
    for j = 1:i-1
        tempM = tempM + kron(U((j-1)*m+1:j*m, :)', C*A^(i-j-1));
    end
    tempM2 = [tempM2; tempM];
end
tempM = [tempM1 tempM2 tempM3];
Theta = pinv(tempM)*Y;
B = reshape(Theta((n+1):(n+n*m)), [n, m]);
D = reshape(Theta((n+n*m+1):end), [p, m]);
% D = Ts_thetau(1:p, 1:m);
% B = pinv(Os(1:(s-1)*p, :)) * Ts_thetau((p+1):end, 1:m);
% B = B + K*D;
end