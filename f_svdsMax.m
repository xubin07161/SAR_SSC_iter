function [U,S,V] = f_svdsMax(A)
%   See also SVD, EIGS.

[m,n] = size(A);
B = [sparse(m,m) A; A' sparse(n,n)];
k=1;

if nnz(A) == 0
    U = eye(m,k);
    S = zeros(k,k);
    V = eye(n,k);
    return
end

bk = 1;
if isreal(A)
    bsigma = 'LA';
else
    bsigma = 'LR';
end
boptions.tol = 5e-6 / sqrt(2);
boptions.disp = 0;


[W,d] = eigs(B,bk,bsigma,boptions);
% [W,d] = f_eigsMax(B);

U = sqrt(2) * W(1:m);
S = d;
V = sqrt(2) * W(m+(1:n));

end


