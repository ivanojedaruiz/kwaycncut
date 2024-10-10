function [x,f] = uzawa(Ar,A,B,c,tol,maxit); 

addpath('/Users/y_l39/COMPUTE/ifem-master/solver'); 
option.solvermaxit = 100;
option.preconditioner = 'V';

r = 100000000000000; 
nf = size(A,1); 
ng = size(B,1);
[row_ind,col_ind,valb] = find(B);
w = sparse(zeros(nf,1));
w(col_ind(:)) = c(row_ind(:));

%f = -A*w;
f = r*B'*B*w;
g = sparse(zeros(ng,1));
u = zeros(nf,1); 
p = zeros(ng,1); 
err = 1;
itnum = 0; 
AugA = A + r*B'*B;

while and(err > tol, itnum < maxit) 
    itnum = itnum + 1;
    f2 = f - B'*p; 
    ubackslash = AugA\f; 
    % uminres = minres(AugA,f,1/r);
    % uamg = amg(AugA,f,option,[]);

    u2=ubackslash;
    p2 = p + r*B*u2;
    res1 = f - A*u2 - B'*p2; 
    res2 = -B*u2;
    res = [res1;res2]; 
    err = norm(res)/norm(f);
    u = u2; 
    p = p2; 
end 
stopERR = err; 

fprintf('Augmented Lagrangian Uzawa Method is applied:\n');
fprintf('Iteration: %2.0u, Stop-error = %8.2e\n',itnum, stopERR);

%x = u + w;



x = ubackslash;
% n=size(x,1);
% norm(ubackslash-uminres,2)/n
% norm(ubackslash-uamg,2)/n
% norm(uminres-uamg,2)/n
% pause 
% 
% max(ubackslash)
% max(uamg)
% max(uminres)
% pause
