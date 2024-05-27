%% This function is aimed to calculate the first derivative of the SCAD penalty
%  arguments:
%  x-- the covaraites need to be regularized
%  lambda-- the tuning parameter for the regularized framework
%  a-- another tuning parameter for the SCAD 
% Author: Xiang Li
% E-mail: 12135003@zju.edu.cn
% Release: 1.0
% Release date: 2024/05/26

function SD = SCAD_deriv(x, lambda , a )
  x = abs(x);
  SD1 = (a*lambda >= x) .* (a * lambda - x) / ((a - 1) * lambda);
  SD = (x <= lambda) * lambda + (x > lambda) .* SD1 * lambda;
end
