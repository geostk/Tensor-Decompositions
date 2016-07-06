% Function to find core tensor (Tucker decomposition)

function [Core_Tensor]=find_Core_Tensor(Tensor_A,Singular_Factors_A)

% Inputs
% Tensor_A             : Tensor 'A' of Mode 'n' 
% Singular_Factors_A   : Left Singular values obtained by using svd on each
%                        of tensor (cell array or N*1 : in order with modes)
% 
% Outputs
% Core_Tensor          : Core Tensor of Tensor A (computed according to
%                        HOSVD)
% 
% Author               : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update          : 17/05/2016

% % 


iter=ndims(Tensor_A);

Kproduct=Tensor_A;
for i=1:iter
    Kproduct= ttm(Kproduct,Singular_Factors_A{i,1}',i);
end

Core_Tensor=Kproduct;

end