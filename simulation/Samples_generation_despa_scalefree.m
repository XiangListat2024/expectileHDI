%% This function is aimed to creat the samples with scale-free graph covariance matrix
%  arguments:
%  n-- the sample size
%  p-- the dimension
%  tau-- the expectile level
%  errtype-- errtype ==1 for the N(0,1) error,  errtype ==2 for the t_4 error
%  gamma_tr-- true coefficients for the model
%  k-- random seed
%  modtype-- ==1 for the homo model case, ==2 for the hete model case
%  t--  a parameter for the hete case in the model 
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26

function [dat] = Samples_generation_despa_scalefree(n,errtype,p,tau,gamma_tr,k,modtype,t)
rng(k);
mu = zeros(1,p);
sigama = zeros(p,p);

for a = 1:p
for b = 1:p
% For different choices of the matrix. The folllowing condition needs to be changed.
sigama(a,b) = 0.3 *( abs(a-b)>1e-4  );
end
end
eig1 = eig(sigama);
for a = 1:p
for b = 1:p
D(a,b) = (a == b) *( 1 + 2*(a>199)  );
end
end
gene1 = D*( sigama + ( 0.2+ abs(min(eig1)) )*eye(p))*D;
sigama = inv(gene1);


Z = mvnrnd(mu,sigama,n);
Z_new = Z;

for a = 1:n
    Z_new(a,1) = normcdf(Z(a,1));
end
    efun_m =[-0.861592112923972,-1.15470052347519;-0.549155820917804,-0.707106775371358;-0.337119881442050,-0.426824202004354;-0.161657506250776,-0.203079914354021;0,0;0.161657506250776,0.203079914354021;0.337119881442050,0.426824202004354;0.549155820917804,0.707106775371358;0.861592112923972,1.15470052347519];
if (modtype ==1)
    if (errtype==1)

    erro0 = normrnd(0,1,n,1) - efun_m(10*tau,1);

    end
    if(errtype==2)

    erro0 = trnd(4,n,1) - efun_m(10*tau,2);

    end
    if (errtype ==3)
        qq =unifrnd(0,1,n,1);
        erro0  = (qq>0.1).*(normrnd(0,1,n,1)- efun_m(10*tau,1) )+(qq<=0.1).*(trnd(4,n,1) - efun_m(10*tau,2) );
    end

    Y =  Z*gamma_tr +erro0;
end

if(modtype ==2)
    if (errtype==1)
    erro0 = normrnd(0,1,n,1)  ;
    end
    if(errtype==2)
    erro0 = trnd(4,n,1) ;
    end
    if (errtype ==3)
        qq =unifrnd(0,1,n,1);
        erro0  = (qq>0.1).*(normrnd(0,1,n,1))+(qq<=0.1).*(trnd(4,n,1));
    end

Y =   Z*gamma_tr + t*Z_new(:,1).*erro0 ;
end
dat =[Y,erro0,Z];
end 