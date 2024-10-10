function W = computeW(imageX,dataW,emag,ephase)
% W = computeW(imageX,dataW,emag,ephase)
% Timothee Cour, Stella Yu, Jianbo Shi, 2004.
 [p,q] = size(imageX);
% s=dataW.sample_rate
% r=dataW.sampleRadius
% pause
%tic
%[w_i,w_j] = global_neigh(p,q,dataW.sampleRadius,dataW.sample_rate); %cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate);

%[w_i(1:20),w_j(1:20)]
%toc
% tic
[w_i,w_j] = cimgnbmap([p,q],dataW.sampleRadius,dataW.sample_rate);
% toc

%[m,id]=max(i-w_i)  %Make sure sample_rate=1 before using this lines
%sum(j-w_j)

%W = affinityic_param2(emag,ephase,w_i,w_j,max(emag(:)) * dataW.edgeVariance);
W = affinityic(emag,ephase,w_i,w_j,max(emag(:)) * dataW.edgeVariance);
%W = max(W,W');
W = W/max(W(:));
