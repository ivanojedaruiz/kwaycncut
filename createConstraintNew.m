function [B,c] = createConstraintNew(pts,nr,nc)
%
%   purpose:
%      create the matrix B for the constraints Bv=c
%   Input:
%      pts1: points selected and labeled for cluster 1
%      pts2: points selected and labeled for cluster 2
%      nr:   number of rows of image I
%      nc:   number of columns of image I
%   Output:
%      C: the constraints matrix
%
n = nr*nc;

for i=1:size(pts,2)
    BP{i} = round(pts{i});
    x{i} = BP{i}(:,1);
    y{i} = BP{i}(:,2);
    m{i} = size(BP{i},1);
end
%% homogenous

% if (boolMustLinkOnly == 0)
%     m = m1+m2-1;
%     C = zeros(m1+m2-1,n);
%     idx1 = (x1(1)-1)*nr+y1(1);
%     for i = 2:m1
%         idx = (x1(i)-1)*nr+y1(i);
%         C(i-1, idx1)=1;
%         C(i-1, idx)= -1;
%     end
%     for i = m1+1:m1+m2
%         idx = (x2(i-m1)-1)*nr+y2(i-m1);
%         C(i-1, idx1)=1;
%         C(i-1, idx) =1;
%     end
% else
%     m = m1-1;
%     C = zeros(m1-1,n);
%     idx1 = (x1(1)-1)*nr+y1(1);
%     for i = 2:m1
%         idx = (x1(i)-1)*nr+y1(i);
%         C(i-1, idx1)=1;
%         C(i-1, idx)= -1;
%     end
% end
%% inhomogenous

% if (boolMustLinkOnly == 0)
%% Initialize
    mtotal=0;
    for i = 1:size(pts,2)
        m{i} = size(pts{i},1);
        mtotal = mtotal + m{i};
    end
    c = zeros(mtotal,size(pts,2));
    B = zeros(mtotal,n);

%% Generate values 
for i = 1:m{1}
        idx = (x{1}(i)-1)*nr+y{1}(i);
        B(i, idx)=1;
        c(i,1) = 1;%1/sqrt(n);
end

start = 0;
for k = 2:size(pts,2)
    start = start + m{k-1};
    for i = start + 1: start + m{k}
        idx = (x{k}(i-start)-1)*nr+y{k}(i-start);
        B(i, idx)=1;
        c(i,k) = 1;%1/sqrt(n);
    end
end 

% for i = m{k}+1:mtotal
%         i-m{k};
%         idx = (x{k+1}(i-m{k})-1)*nr+y{k+1}(i-m{k});
%         B(i, idx)=1;
%         c(i,2) = 1;%1/sqrt(n);
% end

% else
%     m = m1;
%     B = zeros(m,n);
%     c = zeros(m,1);
%     for i = 1:m
%         idx = (x1(i)-1)*nr+y1(i);
%         B(i, idx)=1;
%         c(i) = 1/sqrt(n);
%     end
% end


% B = zeros(m1+m2,n);
% c = zeros(m1+m2,1);
% for i = 1:m1
%     j = (x1(i)-1)*nr+y1(i);
%     B(i,j) = 1;
%     c(i) = 1/sqrt(n);
% end
% for i = m1+1:m1+m2;
%     j = (x2(i-m1)-1)*nr+y2(i-m1);
%     B(i,j) = 1;
%     c(i) = -1/sqrt(n);
% end