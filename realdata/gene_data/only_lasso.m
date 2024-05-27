%% This function is aimed to calculate an estimator under the regularized expectile framework with Lasso regularizer.
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26
function [output] = only_lasso(Q_1,Q_2,lambda,tau)
% Q_1 response ; Q_2 covariates; tau-expectile level
[n,p] = size(Q_2);
tic;
  cvx_begin 
    cvx_quiet(true)
    variables gamma_1(p,1)
    minimize  sum(tau*square_pos(Q_1 -Q_2*gamma_1) + (1-tau)*square_pos( -(Q_1 -Q_2*gamma_1))) /n + 2*lambda*norm(gamma_1,1)
  cvx_end  
toc;

output  = [gamma_1'];
end