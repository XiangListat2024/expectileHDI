%% This function is main-part for the test problem with SCAD-SCAD method 
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26

%% Part one .  Initialize the model to be tested
cvx_solver mosek  % most efficient solver  
n = 300; % n--sample size ,
p = 40; % p--dimension of covariate 
errtype =1; % 1 -- N(0,1) , 2--t_4
modtype =2; % 1 -- homo case , 2--hete case
tau = 0.1;  % expectile level
alpha = 0.05; % significant level
xi =0.75;   % parameter for the covariance matrix of the design under Toeplitz case
t = 0.7;    % parameter for the hete case in the model 
a = 3.7;    % parameter for the SCAD 
%initialize the true coefficients
 gamma_tr = zeros(p,1);
 gamma_tr(15,1)  =1;
 gamma_tr(6,1)   =1;
 gamma_tr(12,1)  =1;
 gamma_tr(20,1)  =1;
 gamma_tr (2,1) = 0/sqrt(n);
 gamma_tr (1,1) = 1/sqrt(n);

REP = 1000;  % the repititions 

lasso_result = zeros(REP,p);  % record the lasso_estimator
result = zeros(REP,6);        % record related results for de-baising version.
 

q = 20231111;   %random seed setting
dat  = Samples_generation_despa(n,errtype,p,tau,gamma_tr,q,modtype,xi,t);
%dat = Samples_generation_despa_scalefree(n,errtype,p,tau,gamma_tr,q,modtype,t);

Y = dat(:,1);
error = dat(:,2);
Z = dat(:,3:end); 


%% Part 2---Five-fold cross-validation
fold_index = zeros(5,60); % Take n =300 for example , split the sample
fold_index(1,:) = linspace(1,60,60);
fold_index(2,:) = linspace(61,120,60);
fold_index(3,:) = linspace(121,180,60);
fold_index(4,:) = linspace(181,240,60);
fold_index(5,:) = linspace(241,300,60);  

nlambda = 50;  % set the grid for the choice of tuning parameter lambda
X = Z;
lambda1 = linspace(0.01,1.26,nlambda)*sqrt(log(p)/n);
rss2 = zeros(1,nlambda);
for j = 1:nlambda
    for k = 1:5
       [j,k]
       x = removerows(Z,'ind',fold_index(k,:));
       y = removerows(Y,'ind',fold_index(k,:));
       beta2 = lla_spa(y,x,lambda1(j),a,tau)' ;
       beta2 = beta2 .* (abs(beta2) > 1e-03);
       sz = sum (abs(beta2) > 1e-03);
       rss2(j) = rss2(j) + sum( tau*square_pos(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta2) +(1-tau)*square_pos(-(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta2)) )/length(fold_index(1,:)) ;
    end
 end
rss2 = rss2/5;
rss2
[~,minind] = min(rss2);
lambda =lambda1(minind);


%% Part 3---De-biasing procedure

for i =1 :REP
    q = i+1101;
   %dat = Samples_generation_despa(n,errtype,p,tau,gamma_tr,q,modtype,xi,t);
    dat = Samples_generation_despa_scalefree(n,errtype,p,tau,gamma_tr,q,modtype,t);

    Y = dat(:,1);
    error = dat(:,2);
    Z = dat(:,3:end); 
    lasso_result(i,:) =  lla_spa(Y,Z,lambda,a,tau);
    lasso_result(i,:)  = lasso_result(i,:) .*(abs(lasso_result(i,:)) >1e-4);
    
i
    %Here we take single variable for example.
    r_2esult(i,:) =  debiased_lasso_spa_scad(Y,Z, lasso_result(i,:) ,tau,lambda,2,alpha);
 
end

%%Part 4. save the results.

%save('12-8-homo-300400-scalefreedirac4-full-e1-6-025-01-scsc.mat')
 

