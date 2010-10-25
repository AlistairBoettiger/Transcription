
clear all;

n=2;

% Specific Model Constants
k = sym('k','real');

% Parallel submatrices
G{1} = [[ - k, k, 0];
        [k, -2*k, k];
        [0, k, -k]];
G{2} = [[ - k, k];
        [k, -k,]];    
    
% G{2} = [[ - k, k, 0];
%         [k, -2*k, k];
%         [0, k, -k]];

R = cell(1,n); L=R; Lambda = cell(1,n); lam = cell(1,n); 
for i=1:n
    [R{i},Lambda{i}] = eig(G{i}); % columns of R are the evecs, diagnol of Lambda are the evals  
    [L{i},Lambda{i}] = eig(G{i}'); % columns of R are the evecs, diagnol of Lambda are the evals  
    lam{i} = diag(Lambda{i});
end

% outer products for evals and evecs of big matrix
Lambdas = kron(lam{1:end}); % eigenvalues of big matrix 
Rs = kron(R{1:end}); % columns of this are the Right eigenvectors of the big matrix 
Ls = kron(L{1:end});

% s / x is indexing the elements of eigenvectors 
ind = Lambdas~=0;
nind = logical(1-ind); 

% M(derivative, s or f)

 M(1,1) = sum(Rs(1,nind).*Ls(end,nind)); % zero eigenvalue terms
 M(1,2) = sum(Rs(end,nind).*Ls(end,nind));

 % for first moments
 M(2,1) = sum(1./Lambdas(ind)'.*Rs(1,ind).*Ls(end,ind)); % start --> finish
 M(2,2) = sum(1./Lambdas(ind)'.*Rs(end,ind).*Ls(end,ind)); % finish --> finish 
 
 % for second moments
 M(3,1) = sum(2./Lambdas(ind)'.^2.*Rs(1,ind).*Ls(end,ind));
 M(3,2) = sum(2./Lambdas(ind)'.^2.*Rs(end,ind).*Ls(end,ind));

% moments of big matrices
 m1 = (M(2,1)*M(1,2) - M(1,1)*M(2,2) )/M(1,2)^2;
 m2 = (-M(3,1)*M(1,2)^2 + M(1,1)*M(1,2)*M(3,2) + 2*M(1,2)*M(2,1)*M(2,2) - 2*M(1,1)*M(2,2)^2)/M(1,2)^3; 
 
 




% limit rho --> infty
% E^x[exp(-lambda tau)] = (lambda*I - G)^{-1}_{x,f}/(lambda*I - G)^{-1}_{f,f}  
% (lambda I - G)^{-1}_{xy} = sum_{i_1..i_k}[(1-lambda prod_k[lambda_{i_k}])^{-1}*
% prod_k[r^k_{i_k}(x_k) l^k_{i_k}(y_k)]]


