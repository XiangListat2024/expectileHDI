%% Real Data Analysis -- gene data
% Author : Xiang Li
% E-mail : 12135003@zju.edu.cn
% Release: 1.0
% Release data: 2024/05/26

%% Load the data
X = csvread('x_real.csv',1); 
Y = csvread('y_real.csv',1);
n = size(X,1);
p = size(X,2);
X = [ones(n,1) X];  % add 1's as the 1st column of X corresponding to the intercept term

%% Five-fold cross-validation (Lasso-Lasso)
fold_index = zeros(5,floor(n/5));
fold_index(1,:) = linspace(1,floor(n/5) ,floor(n/5));
fold_index(2,:) = linspace(floor(n/5) + 1, 2*floor(n/5), floor(n/5));
fold_index(3,:) = linspace(2* floor(n/5) + 1, 3*floor(n/5), floor(n/5));
fold_index(4,:) = linspace(3* floor(n/5) + 1, 4*floor(n/5), floor(n/5));
fold_index(5,:) = linspace(4* floor(n/5) + 1, 5*floor(n/5), floor(n/5));  

% Define a series of expectile-levels 
tau_vec = [0.1, 0.3, 0.5, 0.7, 0.9 ];
lambda_vec = [];
for i = 1: length(tau_vec)
    lambda_vec(i,:) = [cv_gene_la(X,Y,fold_index,tau_vec(i)) , tau_vec(i)];
end

%% Compute initial Lasso estimators
%def a record matrix.
result_lasso_mat = zeros(length(tau_vec), p+1 );

for i = 1: length(tau_vec)
 result_lasso_mat(i,:) = only_lasso(Y,X,lambda_vec(i,1),lambda_vec(i,2));
 result_lasso_mat(i,:) = result_lasso_mat(i,:) .* (abs(result_lasso_mat(i,:)) > 1e-03);
end 


%% Compute de-biased estimators for each single variable and the significant genes.

alpha = 0.05;

de_result_lasso_mat_01 = zeros(p+1,6);
de_result_lasso_mat_03 = zeros(p+1,6);
de_result_lasso_mat_05 = zeros(p+1,6);
de_result_lasso_mat_07 = zeros(p+1,6);
de_result_lasso_mat_09 = zeros(p+1,6);


for i = 1: p+1
de_result_lasso_mat_01(i,:) = debiased_lasso_spa(Y,X, result_lasso_mat(1,:) ,tau_vec(1),lambda_vec(1,1),i,alpha);
de_result_lasso_mat_03(i,:) = debiased_lasso_spa(Y,X, result_lasso_mat(2,:) ,tau_vec(2),lambda_vec(2,1),i,alpha);
de_result_lasso_mat_05(i,:) = debiased_lasso_spa(Y,X, result_lasso_mat(3,:) ,tau_vec(3),lambda_vec(3,1),i,alpha);
de_result_lasso_mat_07(i,:) = debiased_lasso_spa(Y,X, result_lasso_mat(4,:) ,tau_vec(4),lambda_vec(4,1),i,alpha);
de_result_lasso_mat_09(i,:) = debiased_lasso_spa(Y,X, result_lasso_mat(5,:) ,tau_vec(5),lambda_vec(5,1),i,alpha);
end

% Mark the significant genes at given alpha with different expectile levels.
ind_01  = find( de_result_lasso_mat_01(:,4) >1e-4);
ind_03  = find( de_result_lasso_mat_03(:,4) >1e-4);
ind_05  = find( de_result_lasso_mat_05(:,4) >1e-4);
ind_07  = find( de_result_lasso_mat_07(:,4) >1e-4);
ind_09  = find( de_result_lasso_mat_09(:,4) >1e-4);


