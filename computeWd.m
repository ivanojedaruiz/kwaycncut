function [W,d] = computeWd(Img)
% Compute semilarity matrix
W = ICgraph(Img);
W = sparsifyc(W,1e-6);
% check for matrix symmetry
if max(max(abs(W-W'))) > 1e-10 %voir (-12) 
    disp(max(max(abs(W-W'))));
    error('W not symmetric');
end
d = sum(abs(W),2);
