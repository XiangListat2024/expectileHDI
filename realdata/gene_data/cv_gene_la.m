%% Cross validation with the Lasso regularizer.
% X-- covariates Y-- response
% fold_index -- cv index
% tau -- expectile level
% Author : Xiang Li
% E-mail : 12135003@zju.edu.cn
function [lambda_use] = cv_gene_la(X,Y,fold_index,tau)
nlambda = 100;
[n,p] = size(X);
lambda1 = exp(linspace(log(0.001),log(50),nlambda));
%rss1 = zeros(1,nlambda);
rss1 = zeros(1,nlambda);
for j = 1:nlambda
    for k = 1:5
       [j,k]
       x = removerows(X,'ind',fold_index(k,:));
       y = removerows(Y,'ind',fold_index(k,:));
       beta1 = only_lasso(y, x, lambda1(j),tau)';
       beta1 = beta1 .* (abs(beta1) > 1e-03);
       rss1(j) = rss1(j) +  sum( tau*square_pos(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta1) +(1-tau)*square_pos(-(Y(fold_index(k,:))-X(fold_index(k,:),:)*beta1)) )/length(fold_index(1,:)) ;
    end
end
rss1 = rss1/5;
%rss1
[~,minind] = min(rss1);
lambda_use = lambda1(minind) ;

end