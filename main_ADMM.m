function [Y_hat, A_hat, B_hat, C_hat, D_hat, K_hat] = main_ADMM(U, Y, m, n, p, N, s, A, B, C, D, K, beta, rho, beta_max, tau, epsilon, max_loop)
%MAIN_ADMM ADMM算法的主流程
%   U为输入数据,U=[u(1);u(2);...;u(N)];Y为输出数据,Y=[y(1);y(2);...;y(N)];
%   s为hankel矩阵行块数;
Y_hat = Y;
A_new = A - K*C;
B_new = B - K*D;
Theta_u = Construct_ThetaU(A_new, B_new, C, D, s);
theta_u = vec(Theta_u);
Theta_y = Construct_ThetaY(A_new, C, K, s);
theta_y = vec(Theta_y);
Us = hankelize(U, m, N, s);
Ys = hankelize(Y, p, N, s);
[~, F2_u, H_u, ~, G2_u] = matrix_fft(p, m, s, 2*s-1, s);
[~, F2_y, H_y, ~, G2_y] = matrix_fft(p, p, s, 2*s-1, s-1);
[F_yhat, ~, H_yhat, G_yhat, ~] = matrix_fft(p, 1, s, N, N);
M11 = matrix_HrHr(N, F_yhat, H_yhat, G_yhat);
M12 = matrix_HrTr(N, 2*s-1, F_yhat, F2_u, H_yhat, H_u, G_yhat, G2_u*Us);
M13 = matrix_HrTr(N, 2*s-1, F_yhat, F2_y, H_yhat, H_y, G_yhat, G2_y*Ys);
M22 = matrix_TrTr(2*s-1, F2_u, F2_u, H_u, H_u, G2_u*Us, G2_u*Us);
M23 = matrix_TrTr(2*s-1, F2_u, F2_y, H_u, H_y, G2_u*Us, G2_y*Ys);
M33 = matrix_TrTr(2*s-1, F2_y, F2_y, H_y, H_y, G2_y*Ys, G2_y*Ys);
M = [M11 -M12 -M13; -M12' M22 M23; -M13' M23' M33];
Hs_yhat = hankelize_fft(Y_hat, p, 1, N, H_yhat, F_yhat, G_yhat);
Ts_thetau = toplitz_fft(theta_u, p, m, 2*s-1, F2_u, H_u, G2_u);
Ts_thetay = toplitz_fft(theta_y, p, p, 2*s-1, F2_y, H_y, G2_y);
Lambda = zeros(p*s, N-s+1);
theta = [Y_hat; theta_u; theta_y];
for count = 1:max_loop
    Psi = UpdatePsi(Hs_yhat, Ts_thetau, Us, Ts_thetay, Ys, Lambda, beta, n);
    [Y_hat, theta_u, theta_y] = UpdateParameter(n, m, p, N, s, Us, Ys, Y, Lambda, ...
        Psi, beta, F_yhat, H_yhat, G_yhat, F2_u, H_u, G2_u, F2_y, H_y, G2_y, M);
    Hs_yhat = hankelize_fft(Y_hat, p, 1, N, H_yhat, F_yhat, G_yhat);
    Ts_thetau = toplitz_fft(theta_u, p, m, 2*s-1, F2_u, H_u, G2_u);
    Ts_thetay = toplitz_fft(theta_y, p, p, 2*s-1, F2_y, H_y, G2_y);
    Lambda = UpdateLambda(Lambda, tau, beta, Hs_yhat - Ts_thetau*Us - Ts_thetay*Ys - Psi);
    beta = min(beta_max, rho*beta);
    theta_new = [Y_hat; theta_u; theta_y];
    thresh = norm(theta_new-theta, 2)/norm(theta, 2);
    if(thresh < epsilon)
        break;
    else
        theta = theta_new;
    end
end
[A_hat, B_hat, C_hat, D_hat, K_hat] = EstimateStateModel(Hs_yhat, Ts_thetau, Ts_thetay, Us, Ys, n, m, p, s, N, U, Y);
end

