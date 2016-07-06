% Function to decompose tensor in coupled modes
function [Singular_Factors_A,Singular_Factors_B,Singular_Factors_C]=Decompose_Tensor_Coupled_Modes_HOSVD3(Tensor_A,Tensor_B,Tensor_C,Rank,Coupled_Modes_ABC,S_Factors_A,S_Factors_B,S_Factors_C)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N' having 'O' modes coupled
%                         with Tensor_B and Tensor_C
% Tensor_B              : Tensor 'B' of Mode 'M' having 'O' modes coupled
%                         with Tensor_A and Tensor_C
% Tensor_C              : Tensor 'P' of Mode 'P' having 'O' modes coupled
%                         with Tensor_A and Tensor_B
% Rank                  : Rank reduction in SVD
% Coupled_Modes_ABC     : Array of order 3*P, each column specifies the mode
%                         in Tensor A coupled with which mode in Tensor B
%                         and Tensor C similar for other columns
% S_Factors_A           : Singular Factors of Tensor 'A' in cell array
%                         format before decomposing in coupled modes
% S_Factors_B           : Singular Factors of Tensor 'B' in cell array
%                         format before decomposing in coupled modes
% S_Factors_C           : Singular Factors of Tensor 'C' in cell array
%                         format before decomposing in coupled modes
%
%
% Output
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format after decomposing in coupled modes
% Singular_Factors_B    : Singular Factors of Tensor 'B' in cell array
%                         format after decomposing in coupled modes
% Singular_Factors_C    : Singular Factors of Tensor 'C' in cell array
%                         format after decomposing in coupled modes
%
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 05/07/2016

% %
Number_coupled_modes=size(Coupled_Modes_ABC,2);

Singular_Factors_A=S_Factors_A;
Singular_Factors_B=S_Factors_B;
Singular_Factors_C=S_Factors_C;


% Computing the coupled Singular vectors

for i=1:Number_coupled_modes
    Coupled_mode_A=Coupled_Modes_ABC(1,i);
    Coupled_mode_B=Coupled_Modes_ABC(2,i);
    Coupled_mode_C=Coupled_Modes_ABC(3,i);
    
    Matricize_A=tenmat(Tensor_A,Coupled_mode_A);
    Matricize_B=tenmat(Tensor_B,Coupled_mode_B);
    Matricize_C=tenmat(Tensor_C,Coupled_mode_C);
    
    [U,~,~]=svd([Matricize_A.data Matricize_B.data Matricize_C.data]);
    
    Singular_Factors_A{Coupled_mode_A,1}=U(:,1:Rank);
    Singular_Factors_B{Coupled_mode_B,1}=U(:,1:Rank);
    Singular_Factors_C{Coupled_mode_C,1}=U(:,1:Rank);
    
    
end


end