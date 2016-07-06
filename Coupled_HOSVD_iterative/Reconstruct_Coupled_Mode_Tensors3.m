% Function to take multple coupled tensors and their singular factors and
% reproduce their coupled Mode by using existing Factors

function [Refactored_CMode]=Reconstruct_Coupled_Mode_Tensors3(Tensor_A,Tensor_B,Tensor_C,Singular_Factors_A,Singular_Factors_B,Singular_Factors_C,Coupled_Mode_A,Coupled_Mode_B,Coupled_Mode_C)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N' having 'O' modes coupled
%                         with Tensor_B and Tensor_C
% Tensor_B              : Tensor 'B' of Mode 'M' having 'O' modes coupled
%                         with Tensor_A and Tensor_C
% Tensor_C              : Tensor 'P' of Mode 'P' having 'O' modes coupled
%                         with Tensor_A and Tensor_B
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% Singular_Factors_B    : Singular Factors of Tensor 'B' in cell array
%                         format
% Singular_Factors_C    : Singular Factors of Tensor 'C' in cell array
%                         format
% Coupled_Mode_A        : Index of coupled mode of Tensor 'A' with other
%                         tensors
% Coupled_Mode_B        : Index of coupled mode of Tensor 'B' with other
%                         tensors
% Coupled_Mode_C        : Index of coupled mode of Tensor 'C' with other
%                         tensors
%
%
% Output
% Refactored_CMode      : Reconstructed Coupled factor of the Tensor
%
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 05/07/2016

%

[Ten_A]=Tensor_Multiply_Specific_Factors(Tensor_A,Coupled_Mode_A,Singular_Factors_A);
[Ten_B]=Tensor_Multiply_Specific_Factors(Tensor_B,Coupled_Mode_B,Singular_Factors_B);
[Ten_C]=Tensor_Multiply_Specific_Factors(Tensor_C,Coupled_Mode_C,Singular_Factors_C);

Matricize_A=tenmat(Ten_A,Coupled_Mode_A);
Matricize_B=tenmat(Ten_B,Coupled_Mode_B);
Matricize_C=tenmat(Ten_C,Coupled_Mode_C);

[Refactored_CMode,~,~]=svd([Matricize_A.data Matricize_B.data Matricize_C.data]);

end