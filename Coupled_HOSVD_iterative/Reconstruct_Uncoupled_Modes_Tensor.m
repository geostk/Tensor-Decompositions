
% Function to reconstruct uncoupled Modes of the Tensor

function [S_Factors]=Reconstruct_Uncoupled_Modes_Tensor(Tensor_A,Singular_Factors_A,Rank,Coupled_Modes_A,S_Factors_A)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N' 
% Rank                  : Rank reduction in SVD
% Coupled_Modes_A       : index to the coupled modes of Tensor_A
% S_FactorsA            : Singular Factors of the Coupled Modes
%                         reconstrcuted Earlier
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% 
% Output
% S_FactorsA            : Singular Factors of the UnCoupled Modes
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 05/07/2016

%

iter_A=ndims(Tensor_A);
S_Factors=S_Factors_A;
for i=1:iter_A
        Ten_A=Tensor_A;
        if(~any(i==Coupled_Modes_A))
            
            [Ten_A]=Tensor_Multiply_Specific_Factors(Tensor_A,i,Singular_Factors_A);
        
        end
        
        U=tenmat(Ten_A,i);
        [u,~,~]=svd(U.data);
        S_Factors{i,1}=u(:,1:Rank);
        
end

end