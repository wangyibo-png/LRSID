function [ratio2] = SearchRatio(A, B, C, D, u, w_origin, v_origin, n, p, N, ratio1, ratio2_min, ratio2_max, SNR_t)
%SEARCHRATIO 此处显示有关此函数的摘要
%   此处显示详细说明
Q = eye(n)*ratio1;
epsilon = 0.01;
while(true)
    ratio2 = (ratio2_min+ratio2_max)/2;
    R = eye(p)*ratio2;
    w = sqrt(Q)*w_origin;
    v = sqrt(R)*v_origin;
    [~, ~, ~, SNR] = KalmanF(A, B, C, D, Q, R, u, w, v, n, p, N);
    if(abs(SNR-SNR_t)<= epsilon)
        break;
    elseif(SNR>SNR_t)
        ratio2_min = ratio2;
    else
        ratio2_max = ratio2;
    end
end
end

