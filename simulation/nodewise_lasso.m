%% This function is used for calculating the precision matrix
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26
function [output] = nodewise_lasso(Q_1,Q_2,lambda)
[n,p] = size(Q_2);
tic;
  cvx_begin
    cvx_quiet(true)
    variables gamma_1(p)
    minimize sum(square(Q_1 -Q_2*gamma_1)) /n + 2*lambda*norm(gamma_1,1)
  cvx_end  
toc;

tau_squre = sum(square(Q_1 -Q_2*gamma_1)) /n + lambda*norm(gamma_1,1);
output  = [gamma_1',  tau_squre ];
end