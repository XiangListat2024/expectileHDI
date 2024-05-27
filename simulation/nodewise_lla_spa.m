function [output] = nodewise_lla_spa(Y,Z,lambda,a)
%LLA  进行scad  
lla = tic;
[n,p] = size(Z);
gamma_1 = zeros(p,1);
omega = SCAD_deriv(gamma_1,lambda,a);  
k =  10 ; % set maximum iterations
gamma_matrix = zeros(k,p);
output = zeros(1,p+1);
scad = tic;
for count = 1:k   
    cvx_begin
    cvx_quiet(true)
    variables  gamma_1(p,1)
    minimize sum( square(  Y-Z*gamma_1  ) )/n + 2* norm(omega.*gamma_1,1); 
    cvx_end  
    gamma_matrix(count,:) = gamma_1';
    omega = SCAD_deriv(gamma_1,lambda,a);
    if count >1 && norm(gamma_matrix(count,:) - gamma_matrix(count-1,:),1) <1e-3
      output = [gamma_matrix(count,:),  sum(square( Y-Z*gamma_matrix(count,:)')) /n+ norm(SCAD_deriv(gamma_matrix(count,:),lambda,a)*gamma_matrix(count,:)',1) ];
        break
    end
end
toc(scad)
 if norm(output,1)<1e-4
    gamma_1 = [gamma_matrix(end,:)];
    tau_squre_new = sum(square( Y-Z*gamma_1')) /n + SCAD_deriv(gamma_1,lambda,a)*gamma_1';

    output = [ gamma_1,tau_squre_new];
 end
toc(lla)
end