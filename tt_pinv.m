function [M,invS,invN,N] = tt_pinv(tt,l)
% Keisuke Fujii

d=tt.d; % dimensionality of tensor
r=tt.r; % rank of tensor
cr=tt.core; % tensor core
ps=tt.ps; % core position in tt.core
n=tt.n;

% compute tensor core
core_l =  cr(ps(l):ps(l+1)-1); 
core_l = reshape(core_l,[r(l)*n(l),r(l+1)]);
[~,S,~] = svd(core_l,'econ');
invS = diag(1./diag(S)); 

% compute M and invN
M = tt_tensor;
M.d = l; M.r = r(1:l+1); M.n = n(1:l);
M.core = cr(1:ps(l+1)); M.ps = ps(1:l+1);

N = tt_tensor;
N.d = d-l; N.r = r(l+1:end); N.n = n(l+1:end);
N.core = cr(ps(l+1):end);
N.ps = ps(l+1:end) - ps(l+1) + 1;
 
tmp = pinv(full(N)); 
invN = tt_matrix(tmp) ;

end
