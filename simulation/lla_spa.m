%% This function is aimed to calculate an estimator under the regularized expectile framework 
% with SCAD regularizer.
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26

function [output] = lla_spa(Y,Z,lambda,a,tau)
%Y--respond Z--covariates, lambda-tuning parameter,  a-- crucial parameter for SCAD, tau--expectile level 
lla = tic;
[n,p] = size(Z);

gamma_1 = zeros(p,1);
omega = SCAD_deriv(gamma_1,lambda,a);  
k = 10 ;% maximum iteration steps
gamma_matrix = zeros(k,p);
output = zeros(1,p);
scad = tic;
for count = 1:k   
    cvx_begin
    cvx_quiet(true)
    variables  gamma_1(p,1)
    minimize sum((1 - tau) * square_pos(-(Y-Z*gamma_1)) + tau * square_pos(Y -Z*gamma_1))/n + 2* norm(omega.*gamma_1,1); 
    cvx_end  
    gamma_matrix(count,:) = gamma_1';
    omega = SCAD_deriv(gamma_1,lambda,a);
end
toc(scad)
 if norm(output,1)<1e-4
    output = [gamma_matrix(end,:)];
 end
toc(lla)
end

