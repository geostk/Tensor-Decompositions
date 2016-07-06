% Main Script to execute coupled HOSVD for given two tensors and their
% coupled modes and find the core tensors and their singular factors
% 
% 
% Tensor_A          : First Tensor of dimension 'N' coupled with 'P'
%                     Modes of second  and third Tensor
% Tensor_B          : Second Tensor of dimension 'M' coupled with 'P'
%                     Modes of first and third Tensor
% Tensor_C          : Third Tensor of dimension 'K' coupled with 'P'
%                     dimensions of first and Second Tensor
% Coupled_Modes_ABC : Array of order 3*P, each column specifies the mode
%                     in Tensor A coupled with which mode in Tensor B and
%                     Tensor C (each tensors coupled modes are row entries)
% option            : Takes an arugment '1' or '2' and returns left or right 
%                     singular vecctors. 1-Left or 2-Right
% Rank              : Rank of the core tensor as desired
% 
% Author            : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update       : 22/06/2016

% % 



% clear all
clc
option=1;
Rank=2; 
Coupled_Modes_ABC=[1;1;1];

% % Decomposing Tensors with and without Coupling 
Tensor_A=Tensor_Activity;
Tensor_B=Tensor_Click;
Tensor_C=Tensor_Impression;


[Coupled_Core_Tensor_A,Coupled_Singular_Factors_A,Coupled_Core_Tensor_B,Coupled_Singular_Factors_B,Coupled_Core_Tensor_C,Coupled_Singular_Factors_C]=Coupled_HOSVD(Tensor_A,Tensor_B,Tensor_C,Coupled_Modes_ABC,option,Rank);

[Core_Tensor_A,Singular_Factors_A]=Decompose_Tensore_HOSVD(Tensor_A,option,Rank);

[Core_Tensor_B,Singular_Factors_B]=Decompose_Tensore_HOSVD(Tensor_B,option,Rank);

[Core_Tensor_C,Singular_Factors_C]=Decompose_Tensore_HOSVD(Tensor_C,option,Rank);

% % Reconstructing Tensors

[Reconstructed_Coupled_Tensor_A]=Reconstruct_Tensor(Coupled_Core_Tensor_A,Coupled_Singular_Factors_A);

[Reconstructed_Coupled_Tensor_B]=Reconstruct_Tensor(Coupled_Core_Tensor_B,Coupled_Singular_Factors_B);

[Reconstructed_Coupled_Tensor_C]=Reconstruct_Tensor(Coupled_Core_Tensor_C,Coupled_Singular_Factors_C);

[Reconstructed_Tensor_A]=Reconstruct_Tensor(Core_Tensor_A,Singular_Factors_A);

[Reconstructed_Tensor_B]=Reconstruct_Tensor(Core_Tensor_B,Singular_Factors_B);

[Reconstructed_Tensor_C]=Reconstruct_Tensor(Core_Tensor_C,Singular_Factors_C);


% %  ERROR Computation

Error_Coupled_Tensor_A = norm(Tensor_A)-norm(Reconstructed_Coupled_Tensor_A)
Error_Tensor_A = norm(Tensor_A)-norm(Reconstructed_Tensor_A)

Error_Coupled_Tensor_B = norm(Tensor_B)-norm(Reconstructed_Coupled_Tensor_B)
Error_Tensor_B = norm(Tensor_B)-norm(Reconstructed_Tensor_B)

Error_Coupled_Tensor_C = norm(Tensor_C)-norm(Reconstructed_Coupled_Tensor_C)
Error_Tensor_C = norm(Tensor_C)-norm(Reconstructed_Tensor_C)






