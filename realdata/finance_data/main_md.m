%% Real Data Analysis -- finance macro data
% Author : Xiang Li
% E-mail : 12135003@zju.edu.cn
% Release: 1.0
% Release data: 2024/05/26


%% Import data.
%read from the second-row and the first-column.
data = csvread('clean_2023-10-monthly.csv',1,0);

% Here the response variable is S&P500 col-ind = 74.
Y = data(:,74);
X = data;
X(:,74) = [];

%generate lag-terms from initial data and form the covariates used in the
%regression
    X_lag1 = X;
    Y(1) = [];  
    X_lag1(end,:) = [];
    X_lag2 = X_lag1;
    Y(1) = [];  
    X_lag1(1,:) = [];
    X_lag2(end,:) = [];
    X_lag3 = X_lag2;
    Y(1) = [];  
    X_lag1(1,:) = [];
    X_lag2(1,:) = [];
    X_lag3 (end,:) = [];
[n,~] = size(Y); % Add an intercept item.
X_combine = [ones(n,1),X_lag1,X_lag2,X_lag3];     

X = X_combine;
[n,p] = size(X);
p = p-1;

%% Five-fold cross-validation part
fold_index = zeros(5,floor(n/5));
fold_index(1,:) = linspace(1,floor(n/5) ,floor(n/5));
fold_index(2,:) = linspace(floor(n/5) + 1, 2*floor(n/5), floor(n/5));
fold_index(3,:) = linspace(2* floor(n/5) + 1, 3*floor(n/5), floor(n/5));
fold_index(4,:) = linspace(3* floor(n/5) + 1, 4*floor(n/5), floor(n/5));
fold_index(5,:) = linspace(4* floor(n/5) + 1, 5*floor(n/5), floor(n/5));

%define tau level (expectile)

tau_vec = [ 0.1, 0.3 , 0.5, 0.7, 0.9 ];
lambda_vec = [];

% Chosen by lasso 
for i = 1: length(tau_vec)
    lambda_vec(i,:) = [cv_md(X,Y,fold_index,tau_vec(i)) , tau_vec(i)];
end


%% Compute initial Lasso estimates
%def a record matrix.
result_lasso_mat = zeros(length(tau_vec), p+1 );

for i = 1: length(tau_vec)
 result_lasso_mat(i,:) = only_lasso(Y,X,lambda_vec(i,1),lambda_vec(i,2));
 result_lasso_mat(i,:) = result_lasso_mat(i,:) .* (abs(result_lasso_mat(i,:)) > 1e-03);
end 

%% Compute de-biased estimators for each single variable.

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



%% Test for group variables
% find the index of the critical-varibale corresponding with M1SL,BOGMBASE [ 65,68 ];
% set the specific index sets of the critical-variable, respectively.
test_g_65 = [ 65, 65+p/3 ,65+2*p/3 ];

test_g_68 = [ 68, 68+p/3 ,68+2*p/3 ];
 

for i = 1: length(lambda_vec(:,1))
    
    gp65_result_la(i,:) = debiased_lasso_groupG(Y,X,result_lasso_mat(i,:),tau_vec(i),lambda_vec(i,1),test_g_65,alpha);
  
    gp68_result_la(i,:) = debiased_lasso_groupG(Y,X,result_lasso_mat(i,:),tau_vec(i),lambda_vec(i,1),test_g_68,alpha);
  
end
