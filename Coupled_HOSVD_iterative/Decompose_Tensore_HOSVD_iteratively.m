% Function to find HOSVD of a single complete tensor

function [Core_Tensor_A,Singular_Factors_A]=Decompose_Tensore_HOSVD_iteratively(Tensor_A,Rank,Error_Threshold,Max_iterations)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N'
% Rank                  : Rank reduction in SVD
% Error_Threshold       : Allowable error tolerance limit
% Max_iterations        : Maximum allowable iterations    
% 
% 
% Output
% Core_Tensor_A         : Core Tensor of Tensor 'A' obtained using
%                         iterative HOSVD
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% 
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 04/07/2016

% % 



% Decomposing tensor HOSVD first time 

iter_A=ndims(Tensor_A);
Singular_Factors_A=cell(iter_A,1);

for i=1:iter_A
    Singular_Factors_A{i,1}=find_svd_mode(Tensor_A,i,Rank);
end

% Computing the core tensors

[Core_Tensor_A]=find_Core_Tensor(Tensor_A,Singular_Factors_A);

% Reconstrcut Tensor

[Reconstructed_Tensor_A]=Reconstruct_Tensor(Core_Tensor_A,Singular_Factors_A);

Error=norm(Tensor_A)-norm(Reconstructed_Tensor_A);
iter_number=1;
Singular_Factors=Singular_Factors_A;

while((Error > Error_Threshold) && (iter_number <= Max_iterations))

    S_Factors=cell(iter_A,1);
    for i=1:iter_A
        Ten_A=Tensor_A;
        
        for j=1:iter_A
            
            if(j~=i)
                
                Ten_A=ttm(Ten_A,Singular_Factors{j,1}',j);
            end
        end
        
        U=tenmat(Ten_A,i);
        [u,s,v]=svd(U.data);
        S_Factors{i,1}=u(:,1:Rank);
        
    end
    [Core_Tr]=find_Core_Tensor(Tensor_A,S_Factors);
    [RT]=Reconstruct_Tensor(Core_Tr,S_Factors);
    Error=norm(Tensor_A)-norm(RT)
    Singular_Factors=S_Factors;
    
    
 iter_number=iter_number+1;   
 Core_Tensor_A=Core_Tr;
 Singular_Factors_A=S_Factors;
end        


end

