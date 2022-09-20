function [F1, F2, H, G1, G2] = matrix_fft(n, m, r, N, s)
%MATRIX_FFT 返回与fft相关的矩阵
%   n*m为基础矩阵维度,r为hankel矩阵的block行数,N为总可用数据量
F = fft(eye(N));
H = F(:, (N-r+1):N);
G = F(:, (N-r+1):-1:1);
F1 = kron(eye(m), kron(F, eye(n)));
F2 = kron(eye(m), kron(F(:, (end-s+1):end), eye(n)));
H = kron(ones(m, 1), kron(H, eye(n)));
G = kron(G, ones(n, 1));
index1 = kron((1:N-r+1), eye(m));
index1(index1==0) = N-r+2;
index2 = kron((N-r+1:-1:1), eye(m));
index2(index2==0) = N-r+2;
Gc1 = mat2cell([G zeros(N*n, 1)], N*n, ones(1, N-r+2));
G1 = cell2mat(Gc1(index1));
Gc2 = mat2cell([G zeros(N*n, 1)], N*n, ones(1, N-r+2));
G2 = cell2mat(Gc2(index2));
end