% Function to find left/right singular values of a matrix

function [Singular_factor]=find_svd_mode(Tensor_A,Mode,option,Rank)

% Inputs
% Tensor_A          : Tensor A of mode N
% Mode              : Specifying the matricization Mode of tensor A
% option            : Takes an arugment '1' or '2' and returns left or right 
%                     singular vecctors. 1-Left or 2-Right
% Rank              : Rank reduction in SVD
% 
% Outputs
% Singular_factor   : Singular factor of Tensor_A matricized in Mode 'I'
% 
% Author            : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update       : 17/05/2016

% % 

Matricxe_Tensor=tenmat(Tensor_A,Mode);

[U,S,V]=svd(Matricxe_Tensor.data);


if(option==1)
    Singular_factor=U(:,1:Rank);
    
else
    Singular_factor=V(:,1:Rank);
    
end

end
    