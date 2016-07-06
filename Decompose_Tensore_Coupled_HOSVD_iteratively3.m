%% Function to find HOSVD of a single complete tensor

function [Core_Tensor_A,Singular_Factors_A,Core_Tensor_B,Singular_Factors_B,Core_Tensor_C,Singular_Factors_C]=Decompose_Tensore_Coupled_HOSVD_iteratively3(Tensor_A,Tensor_B,Tensor_C,Rank,Error_Threshold,Max_iterations,Coupled_Modes_ABC)

%Input
% Tensor_A              : Tensor 'A' of Mode 'N' having 'O' modes coupled
%                         with Tensor_B and Tensor_C
% Tensor_B              : Tensor 'B' of Mode 'M' having 'O' modes coupled
%                         with Tensor_A and Tensor_C
% Tensor_C              : Tensor 'P' of Mode 'P' having 'O' modes coupled
%                         with Tensor_A and Tensor_B
% Rank                  : Rank reduction in SVD
% Error_Threshold       : Allowable error tolerance limit
% Max_iterations        : Maximum allowable iterations
% Coupled_Modes_ABC     : Array of order 3*P, each column specifies the mode
%                         in Tensor A coupled with which mode in Tensor B
%                         and Tensor C similar for other columns
%
%
% Output
% Core_Tensor_A         : Core Tensor of Tensor 'A' obtained using coupled
%                         iterative HOSVD
% Singular_Factors_A    : Singular Factors of Tensor 'A' in cell array
%                         format
% Core_Tensor_B         : Core Tensor of Tensor 'B' obtained using coupled
%                         iterative HOSVD
% Singular_Factors_B    : Singular Factors of Tensor 'B' in cell array
%                         format
% Core_Tensor_C         : Core Tensor of Tensor 'C' obtained using coupled
%                         iterative HOSVD
% Singular_Factors_C    : Singular Factors of Tensor 'C' in cell array
%                         format
%
% Author                : Sunny Verma (sunnyverma.iitd@gmail.com)
% Last_Update           : 05/07/2016

%%

% Decomposing tensor HOSVD first time

iter_A=ndims(Tensor_A);
[S_Factors_A]=Decompose_Tensore_HOSVD_Cou_UnCo(Tensor_A,Rank,Coupled_Modes_ABC(1,:));

% Tensor B

iter_B=ndims(Tensor_B);
[S_Factors_B]=Decompose_Tensore_HOSVD_Cou_UnCo(Tensor_B,Rank,Coupled_Modes_ABC(2,:));

% Tensor C

iter_C=ndims(Tensor_C);
[S_Factors_C]=Decompose_Tensore_HOSVD_Cou_UnCo(Tensor_C,Rank,Coupled_Modes_ABC(3,:));

% Decompose coupled modes

Number_coupled_modes=size(Coupled_Modes_ABC,2);
[Singular_Factors_A,Singular_Factors_B,Singular_Factors_C]=Decompose_Tensor_Coupled_Modes_HOSVD3(Tensor_A,Tensor_B,Tensor_C,Rank,Coupled_Modes_ABC,S_Factors_A,S_Factors_B,S_Factors_C);

% Computing the core tensors

[Core_Tensor_A]=find_Core_Tensor(Tensor_A,Singular_Factors_A);
[Core_Tensor_B]=find_Core_Tensor(Tensor_B,Singular_Factors_B);
[Core_Tensor_C]=find_Core_Tensor(Tensor_C,Singular_Factors_C);

% Reconstrcut Tensor

[Reconstructed_Tensor_A]=Reconstruct_Tensor(Core_Tensor_A,Singular_Factors_A);
[Reconstructed_Tensor_B]=Reconstruct_Tensor(Core_Tensor_B,Singular_Factors_B);
[Reconstructed_Tensor_C]=Reconstruct_Tensor(Core_Tensor_C,Singular_Factors_C);

%%

% initiallizing things before moving to iterative error reduction

Error_A=norm(Tensor_A)-norm(Reconstructed_Tensor_A)
Error_B=norm(Tensor_B)-norm(Reconstructed_Tensor_B)
Error_C=norm(Tensor_C)-norm(Reconstructed_Tensor_C)

Error= Error_A + Error_B + Error_C;
iter_number=1;
ISingular_FA=Singular_Factors_A;
ISingular_FB=Singular_Factors_B;
ISingular_FC=Singular_Factors_C;

%%

while((Error > Error_Threshold) && (iter_number <= Max_iterations))
    
    S_FactorsA=cell(iter_A,1);
    S_FactorsB=cell(iter_B,1);
    S_FactorsC=cell(iter_C,1);
     
    % Computing the Factors of Coupled Modes
    
    for i=1:Number_coupled_modes
        
        Coupled_Mode_A=Coupled_Modes_ABC(1,i);
        Coupled_Mode_B=Coupled_Modes_ABC(2,i);
        Coupled_Mode_C=Coupled_Modes_ABC(3,i);
        
        [Refactored_CMode]=Reconstruct_Coupled_Mode_Tensors3(Tensor_A,Tensor_B,Tensor_C,Singular_Factors_A,Singular_Factors_B,Singular_Factors_C,Coupled_Mode_A,Coupled_Mode_B,Coupled_Mode_C);
        
        S_Factors_A{Coupled_Mode_A,1}=Refactored_CMode(:,1:Rank);
        S_Factors_B{Coupled_Mode_B,1}=Refactored_CMode(:,1:Rank);
        S_Factors_C{Coupled_Mode_C,1}=Refactored_CMode(:,1:Rank);
        
        
    end
    
    [S_FactorsA]=Reconstruct_Uncoupled_Modes_Tensor(Tensor_A,Singular_Factors_A,Rank,Coupled_Modes_ABC(1,:),S_Factors_A);
    [S_FactorsB]=Reconstruct_Uncoupled_Modes_Tensor(Tensor_B,Singular_Factors_B,Rank,Coupled_Modes_ABC(2,:),S_Factors_B);
    [S_FactorsC]=Reconstruct_Uncoupled_Modes_Tensor(Tensor_C,Singular_Factors_C,Rank,Coupled_Modes_ABC(3,:),S_Factors_C);
    
%     Reconstructing Tensors based on iterative Calculations

    [Core_TrA]=find_Core_Tensor(Tensor_A,S_FactorsA);
    [RTA]=Reconstruct_Tensor(Core_TrA,S_FactorsA);
    
    [Core_TrB]=find_Core_Tensor(Tensor_A,S_FactorsA);
    [RTB]=Reconstruct_Tensor(Core_TrB,S_FactorsA);
    
    [Core_TrC]=find_Core_Tensor(Tensor_A,S_FactorsA);
    [RTC]=Reconstruct_Tensor(Core_TrC,S_FactorsA);
    
%     New Error
    
    ErrorA=norm(Tensor_A)-norm(RTA)
    ErrorB=norm(Tensor_B)-norm(RTB)
    ErrorC=norm(Tensor_C)-norm(RTC)
    
    Error=ErrorA + ErrorB + ErrorC
      
    
    iter_number=iter_number+1;
    Core_Tensor_A=Core_TrA;
    Core_Tensor_B=Core_TrB;
    Core_Tensor_C=Core_TrC;
    Singular_Factors_A=S_FactorsA;
    Singular_Factors_B=S_FactorsB;
    Singular_Factors_C=S_FactorsC;
end


end

