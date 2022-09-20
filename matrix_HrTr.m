function [M] = matrix_HrTr(N1, N2, F1, F2, H1, H2, G1, G2)
%MATRIX_HRTR 产生辅助矩阵M
%   <Hr(y), Tr(x)> = y'*M*x; Hr*(Tr(x)) = Mx
M = real(1/(N1*N2)*F1'*((H1*H2').*((G2*G1').'))*F2);
end