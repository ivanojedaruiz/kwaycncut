% K-way CNCut CODE
% This code does not use a ground truth
% By Ivan Ojeda-Ruiz

%%%%%%%%%%%%%%%%%%%%
clear all
close all

%% Load Image
imgName = 'elephant';
Imgo = imread([imgName,'.jpg']);

%% Preprocess Image
Imgr = Imgo;%imresize(Imgo,[NaN 100]);
Img=im2gray(Imgr);
[nr,nc]=size(Img);
[nro,nco]=size(Imgo);
n=nr*nc;
Img=rescale(double(Img));

%% Initialize matrix lexicographical locator
u.loc = reshape(linspace(1,n,n),nr,nc);

%% Regular Intervening Contours
[W,D] = computeWd(Img);

%% Choose Points
% choice = 0 uses saved points
% choice = 1 allows the user to click on points
choice = 0;
%filename = ['newpts_',imgName];
if choice ==0
    filename = [imgName,'pts'];
    load(filename)
    %     for k=1:size(pts,2)
    %         newpts{k}(:,2)=pts{k}(:,1)%.*nc./nco;
    %         newpts{k}(:,1)=pts{k}(:,2)%.*nr./nro;
    %         pts{k}(:,2)=newpts{k}(:,2);
    %         pts{k}(:,1)=newpts{k}(:,1);
    %     end
    sets = size(pts,2);
else if choice ==1
        pts = [];
        sets =size(gtseg,3)+3;
        figure();clf;imagesc(Img);truesize([700,700]);
        colormap(gray);
        axis off;set(gca,'position',[0 0 1 1],'units','normalized');
        pts = manuallyChoosePoints(imgName,sets);
        %save([imgName,'pts'],'pts')
end
end

%% Use boundary Constraints and decrease brigthness
% [newImg,pts] = BoundConst(Img,pts,10);
% sets = size(pts,2);
% [W,D] = computeWd(newImg);

save('newpts','pts')

%% Display Points Chosen
figure();clf;
h = imagesc(Imgr);
colormap(gray)
hold on;
% set(gca,'position',[0 0 1 1],'units','normalized');
for i =1:size(pts,2)
    h = plot(pts{i}(:,1),pts{i}(:,2),'o','MarkerSize',5);
    set(h,'linewidth',5);
end
hold off;

%addpts = 1
%while addpts = 1
%% Compute the solution
clear x y
disp('computing B c');
[B, c] = createConstraintNew(pts,nr,nc);
B = sparse(B);
disp('computing L');
d = sum(abs(W),2);
L = spdiags(d,0,size(W,2),size(W,2)) - sparse(W);
disp('computing P');
tp = tic;
Dinvsqrt = 1./(sqrt(d));
D = spdiags(1./(sqrt(d)),0,size(W,2),size(W,2));
A = spmtimesd(L,Dinvsqrt,Dinvsqrt);
tol_ppm =10^(-5);
k_ppm =10^7;
for i =1:sets
    [x(:,i),f(:,i)] = uzawa(A,A,B,c(:,i),tol_ppm,k_ppm);
end





clear h
%% Display results
thresh=0.5;
x(:,4)=x(:,4)*0.7;
[Imgt1,Imgt2,Label]=postProcesskway(x,Img,Img,pts,thresh);
bw = edge(Label,0.01);
figure();clf;
showmask(double(im2gray(Img)),imdilate(bw,ones(2,2)));
hold on;
for i = 1:sets
    h = scatter(pts{i}(:,1),pts{i}(:,2),100,"MarkerFaceColor",rand(1,3))
    %set(h,'linewidth',2*nc/150);
end
hold off;
axis off;
%disp('would you like to add points: 0 = no, 1 = yes')
%if addpts = 1

%end



%% Plot Segments individually

  
figure()
imagesc(Label)
axis off
colorbar

figure()
%h = tiledlayout(1,4, 'TileSpacing', 'compact', 'Padding', 'none');
for i=1:size(x,2)
    subplot(3,4,i)
    imagesc(reshape(x(:,i),nr,nc))
    colorbar
    %truesize
    axis off
end
sgtitle('CNCUT Fuzzy Values')