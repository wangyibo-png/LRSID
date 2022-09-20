function [M] = matrix_TrTr(N, F1, F2, H1, H2, G1, G2)
%MATRIX_TRTR 产生辅助矩阵M
%   x'*M*x = <Tr(x), Tr(x)>; Tr*(Tr(x)) = Mx
M = real(1/N^2*F1'*((H1*H2').*((G2*G1').'))*F2);
end