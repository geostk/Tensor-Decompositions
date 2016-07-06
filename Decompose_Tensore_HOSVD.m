% Function to find HOSVD of a single complete tensor

function [Core_Tensor_A,Singular_Factors_A]=Decompose_Tensore_HOSVD(Tensor_A,option,Rank)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N'
% option                : Takes an arugment '1' or '2' and returns left or right 
%                         singular vecctors. 1-Left or 2-Right
% Rank                  : Rank reduction in SVD
% 
% 
% Output
% Core_Tensor_A         : Core Tensor of Tensor 'A' obtained using coupled
%                         HOSVD
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% 
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 21/06/2016

% % 

iter_A=ndims(Tensor_A);
Singular_Factors_A=cell(iter_A,1);

for i=1:iter_A
    Singular_Factors_A{i,1}=find_svd_mode(Tensor_A,i,option,Rank);
end

% Computing the core tensors

[Core_Tensor_A]=find_Core_Tensor(Tensor_A,Singular_Factors_A);

end
        
        
        

    