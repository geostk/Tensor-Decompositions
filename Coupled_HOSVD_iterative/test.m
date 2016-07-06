X=tensor(rand(4,5,3));
Rank=3;


Tensor_A=X;
% Decomposing tensor HOSVD first time 

iter_A=ndims(Tensor_A);
Singular_Factors_A=cell(iter_A,1);

for i=1:iter_A
    Singular_Factors_A{i,1}=find_svd_mode(Tensor_A,i,Rank);
end

% Computing the core tensors

[Core_Tensor_A]=find_Core_Tensor(Tensor_A,Singular_Factors_A);

% Reconstrcut Tensor

[RTensor_A]=Reconstruct_Tensor(Core_Tensor_A,Singular_Factors_A);

Error=norm(Tensor_A)-norm(RTensor_A)

        
a=ttm(Tensor_A,Singular_Factors_A{2,1}',2);
b=ttm(a,Singular_Factors_A{3,1}',3);
U1=tenmat(b,1);

a=ttm(Tensor_A,Singular_Factors_A{1,1}',1);
b=ttm(a,Singular_Factors_A{3,1}',3);
U2=tenmat(b,2);

a=ttm(Tensor_A,Singular_Factors_A{2,1}',2);
b=ttm(a,Singular_Factors_A{1,1}',1);
U3=tenmat(b,3);

S_factors=cell(3,1);

[u,s,v]=svd(U1.data);
S_factors{1,1}=u(:,1:3);

[u,s,v]=svd(U2.data);
S_factors{2,1}=u(:,1:3);

[u,s,v]=svd(U3.data);
S_factors{3,1}=u(:,1:3);

[Core_Tensor]=find_Core_Tensor(Tensor_A,S_factors);

% Reconstrcut Tensor

[RT]=Reconstruct_Tensor(Core_Tensor,S_factors);

Error=norm(Tensor_A)-norm(RT)

% SF=Singular_Factors_A;

Singular_Factors_A=S_factors;
  
B=tucker_als(X,3);

norm(X)-norm(full(B))
