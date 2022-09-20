function [y_hat, theta_u, theta_y] = UpdateParameter(n, m, p, N, s, Us, Ys, Y, Lambda, Psi, beta, F_yhat, H_yhat, G_yhat, F2_u, H_u, G2_u, F2_y, H_y, G2_y, M)
%UPDATEPARAMETER 此处显示有关此函数的摘要
%   此处显示详细说明
X = Lambda + beta*Psi;
Hs_adj_yhat = hankelize_adj_fft(X, N, F_yhat, H_yhat, G_yhat);
Ts_adj_thetau = toplitz_adj_fft(X, 2*s-1, F2_u, H_u, G2_u*Us);
Ts_adj_thetay = toplitz_adj_fft(X, 2*s-1, F2_y, H_y, G2_y*Ys);
tempTheta = [Y + Hs_adj_yhat; -1*Ts_adj_thetau; -1*Ts_adj_thetay];
tempM = zeros(size(M));
tempM(1:p*N, 1:p*N) = eye(p*N);
tempM = beta*M + tempM;
theta_new = pinv(tempM)*tempTheta;
theta_cell = mat2cell(theta_new, [p*N p*s*m p*(s-1)*p], 1);
y_hat = cell2mat(theta_cell(1));
theta_u = cell2mat(theta_cell(2));
theta_y = cell2mat(theta_cell(3));
end