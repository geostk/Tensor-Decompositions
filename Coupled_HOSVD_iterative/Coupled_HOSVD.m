% Function to take three tensors of mode 'N' and 'M' 
% coupled in 'P' number of modes
% Coupled HOSVD Decomposition

function [Core_Tensor_A,Singular_Factors_A,Core_Tensor_B,Singular_Factors_B,Core_Tensor_C,Singular_Factors_C]=Coupled_HOSVD(Tensor_A,Tensor_B,Tensor_C,Coupled_Modes_ABC,option,Rank)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N' having 'P' modes coupled with
%                         Tensor B and Tensor C
% Tensor_B              : Tensor 'B' of Mode 'M' having 'P' modes coupled with
%                         Tensor A and Tensor C
% Tensor_C              : Tensor 'C' of Mode 'M' having 'P' modes coupled with
%                         Tensor A and Tensor B
% Coupled_Modes_ABC     : Array of order 3*P, each column specifies the mode
%                         in Tensor A coupled with which mode in Tensor B
%                         and Tensor C similar for other columns
% option                : Takes an arugment '1' or '2' and returns left or right 
%                         singular vecctors. 1-Left or 2-Right
% Rank                  : Rank reduction in SVD
% 
% 
% Output
% Core_Tensor_A         : Core Tensor of Tensor 'A' obtained using coupled
%                         HOSVD
% Core_Tensor_B         : Core Tensor of Tensor 'B' obtained using coupled
%                         HOSVD
% Core_Tensor_C         : Core Tensor of Tensor 'C' obtained using coupled
%                         HOSVD
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% Singular_Factors_B    : Singular Factors of Tensor 'B' in cell array
%                         format
% Singular_Factors_C    : Singular Factors of Tensor 'C' in cell array
%                         format
% 
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 22/06/2016

% % 

iter_A=ndims(Tensor_A);
Singular_Factors_A=cell(iter_A,1);

for i=1:iter_A
    Singular_Factors_A{i,1}=find_svd_mode(Tensor_A,i,option,Rank);
end

iter_B=ndims(Tensor_B);
Singular_Factors_B=cell(iter_B,1);

for i=1:iter_B
    Singular_Factors_B{i,1}=find_svd_mode(Tensor_B,i,option,Rank);
end


iter_C=ndims(Tensor_C);
Singular_Factors_C=cell(iter_C,1);

for i=1:iter_C
    Singular_Factors_C{i,1}=find_svd_mode(Tensor_C,i,option,Rank);
end

    
Number_coupled_modes=size(Coupled_Modes_ABC,2);

% Computing the coupled Singular vectors

for i=1:Number_coupled_modes
    Coupled_mode_A=Coupled_Modes_ABC(1,i);
    Coupled_mode_B=Coupled_Modes_ABC(2,i);
    Coupled_mode_C=Coupled_Modes_ABC(3,i);
    
    Matricize_A=tenmat(Tensor_A,Coupled_mode_A);
    Matricize_B=tenmat(Tensor_B,Coupled_mode_B);
    Matricize_C=tenmat(Tensor_C,Coupled_mode_C);
    
    [U,S,V]=svd([Matricize_A.data Matricize_B.data Matricize_C.data]);
    
    if (option==1)
        Singular_Factors_A{Coupled_mode_A,1}=U(:,1:Rank);
        Singular_Factors_B{Coupled_mode_B,1}=U(:,1:Rank);
        Singular_Factors_C{Coupled_mode_C,1}=U(:,1:Rank);
        
    else
        Singular_Factors_A{Coupled_mode_A,1}=V(:,1:Rank);
        Singular_Factors_B{Coupled_mode_B,1}=V(:,1:Rank);
        Singular_Factors_C{Coupled_mode_C,1}=V(:,1:Rank);
    end
    
end

% Computing the core tensors

[Core_Tensor_A]=find_Core_Tensor(Tensor_A,Singular_Factors_A);
[Core_Tensor_B]=find_Core_Tensor(Tensor_B,Singular_Factors_B);
[Core_Tensor_C]=find_Core_Tensor(Tensor_C,Singular_Factors_C);

end
        
        
        

    