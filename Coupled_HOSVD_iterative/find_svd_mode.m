% Function to find left/right singular values of a matrix

function [Singular_factor]=find_svd_mode(Tensor_A,Mode,Rank)

% Inputs
% Tensor_A          : Tensor A of mode N
% Mode              : Specifying the matricization Mode of tensor A
% Rank              : Rank reduction in SVD
% 
% Outputs
% Singular_factor   : Singular factor of Tensor_A matricized in Mode 'I'
% 
% Author            : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update       : 17/05/2016

% % 

Matricze_Tensor=tenmat(Tensor_A,Mode);

[U,S,V]=svd(Matricze_Tensor.data);

Singular_factor=U(:,1:Rank);
    
end
    