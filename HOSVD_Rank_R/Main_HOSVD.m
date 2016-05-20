% Main Script to execute coupled HOSVD for given two tensors and their
% coupled modes and find the core tensors and their singular factors
% 
% 
% Tensor_A          : First Tensor of dimension 'N' coupled with 'P'
%                     dimensions of second Tensor
% Rank_A            : Number of Modes which will form Tensor A
% Tensor_B          : Seconf Tensor of dimension 'M' coupled with 'P'
%                     dimensions of first Tensor
% Rank_B            : Number of Modes which will form Tensor B
% Coupled_Modes_AB  : Array of order 2*P, each column specifies the mode
%                     in Tensor A coupled with which mode in Tensor B
% option            : Takes an arugment '1' or '2' and returns left or right 
%                     singular vecctors. 1-Left or 2-Right
% Rank              : Rank of the core tensor as desired
% 
% Author            : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update       : 17/05/2016

% % 



clear all
clc
% randseed(7);
rng('default') 


Dimension_A=[4,3,5];
Rank_A=4;
Dimension_B=[3,4,5];
Rank_B=4;
Coupled_Modes_AB=[1 2;2 1];
MU_AB=[3 1];
Sigma_AB=[.9 .4];

option=1;
Rank=3; %should be less than min(#dimensions over all tensors)

[Tensor_A,Tensor_B]= Generate_Coupled_Tensors(Dimension_A,Rank_A,Dimension_B,Rank_B,Coupled_Modes_AB,MU_AB,Sigma_AB);

[Core_Tensor_A,Singular_Factors_A,Core_Tensor_B,Singular_Factors_B]=Coupled_HOSVD(Tensor_A,Tensor_B,Coupled_Modes_AB,option,Rank);