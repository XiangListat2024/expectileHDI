%% This function is core-part for the node-wise procedure with SCAD method (Group test)
% eg, H_0 : \beta_1 = \beta_2 = \beta_3 =0
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26
% Y--respond variable
% Z--covariates
% coeff--initial estimator of the coefficients
% tau --expectile level
% lambda_j -- tuning parameter for the node-procedure
% G -- the index set of the variables(group) in the test problem
% alpha -- the significant level

function [output] = debiased_scad_groupG(Y,Z,coeff,tau,lambda_j,G,alpha)
%   Nodewise_lasso for test problem H_0: \beta_{G} = 0
%   Obtain the corresponding size from matrix itself.
  [n,p] = size(Z);
  [~,k_G] = size(G);
  coeff_nodewise = zeros(k_G,p);
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
 
  Theta = zeros(k_G,p);
  debiased = zeros(k_G,1)';
  %d_i = zeros(p,1);
% Use cvx to solve specific nodewise-lasso problem.
  for i = 1 : k_G
       H = ones(1,p);
      Q_1 = Q_w(:,G(i));
      Q_2 = Q_w;
      Q_2(:,G(i)) = [];
      % from cvx calculations
      coeff_nodewise(i,:) = nodewise_lla_spa(Q_1,Q_2,0.95*lambda_j,3.7);
      H(1:G(i)-1) = -coeff_nodewise(i,1:G(i)-1);
      H(G(i)+1:p) = -coeff_nodewise(i,G(i):p-1);
      d_i  = coeff_nodewise(i,p);
  %end
      D =  d_i;
 % Generate the pesduo inverse of 
     Theta(i,:) = D^(-1) * H;
  
 % Generate the de-biased lasso 

 debiased(i) = (coeff(G(i)) + Theta(i,:) * Q_w' * W * epsilon/n );
  
% Generate the estimator of the variance 
  end
hat_Omega =  ( Theta * Q_w' *  (W * diag(epsilon))) * (Theta * Q_w' *  (W * diag(epsilon)))' /n;
output_debiased = n * (debiased/hat_Omega) * debiased' ;
judege = output_debiased > chi2inv(1-alpha,k_G);
p_value =  (1-chi2cdf(output_debiased,k_G));
output = [output_debiased,judege,p_value ];
end