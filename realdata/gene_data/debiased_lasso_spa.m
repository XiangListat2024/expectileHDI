%% This function is core-part for the node-wise procedure with Lasso method 
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26
% Y--respond variable
% Z--covariates
% coeff--initial estimator of the coefficients
% tau --expectile level
% lambda_j -- tuning parameter for the node-procedure
% kk -- the index of the variable(single) in the test problem
% alpha -- the significant level
function [output] = debiased_lasso_spa(Y,Z,coeff,tau,lambda_j,kk,alpha)
% Obtain the corresponding size from matrix itself.
  [n,p] = size(Z);
% Regernate the design.
  Q = Z;
% Calculate the error.
  epsilon = Y - Q * coeff';
% Calculate the weight matrix.  
  w_i = abs(tau - (epsilon < 0 )).^(0.5);
  W = diag(w_i);
% Weight the design by the empirical weight matrix W 
  Q_w = W * Q;
% Define the H_matrix and D_matrix (by the dia-value).
  H = ones(1,p);
  %d_i = zeros(p,1);
% Use cvx to solve specific nodewise-lasso problem.
  %for i = 1 : p
      Q_1 = Q_w(:,kk);
      Q_2 = Q_w;
      Q_2(:,kk) = [];
      % from cvx calculations
      coeff_nodewise = nodewise_lasso(Q_1,Q_2,lambda_j);
      H(1:kk-1) = -coeff_nodewise(1:kk-1);
      H(kk+1:p) = -coeff_nodewise(kk:p-1);
      d_i  = coeff_nodewise(p);
  %end
 D =  d_i;
 % Generate the pesduo inverse of 
 Theta = D^(-1) * H;

 % Generate the de-biased lasso 

 output_debiased = (coeff(kk) + Theta * Q_w' * W * epsilon/n);
% Generate the estimator of the variance  
 estimated_variance = ( Theta * Q_w' *  (W * diag(epsilon))) * (Theta * Q_w' *  (W * diag(epsilon)))' /n;

 upper_bound = output_debiased + norminv(1-alpha/2) * sqrt(  estimated_variance /n) ;
 lower_bound = output_debiased - norminv(1-alpha/2) * sqrt(  estimated_variance /n) ;

 student_test_stat =  sqrt(n)* output_debiased / sqrt(estimated_variance);
 judege_1 = (0 > upper_bound) + (0< lower_bound);
 ci_length = upper_bound - lower_bound;
 judege = abs( student_test_stat) > norminv(1-alpha/2);
 p_value = 2*(1-normcdf(abs(student_test_stat)));
 output = [output_debiased,estimated_variance,student_test_stat,judege,p_value,ci_length ];
end