%% Algebraic Multigrid Method
% Coarsening algorithm as proposed in 
% Achi Brandt (2002)
% Code written by Iván Ojeda-Ruiz

clear all; close all;

%% Berkeley Dataset
imgName = 'cowboy'; 
Imgo = imread([imgName,'.jpg']);

%% Simple Example
%Imgo =[0,0,0,0,0;0,1,1,1,0;0,1,1,1,0;0,0,1,0,0;0,0,0,0,0]

%% Preprocess Image
Imgr = imresize(Imgo,[NaN 100]);
Img=im2gray(Imgr);
[nr,nc]=size(Img);
n=nr*nc;
Img=rescale(double(Img)); %Change values to [0,1]

%% Display grayscale image
imagesc(Img)
colormap('gray')

%% Affinity Matrix 
[W{1},D] = computeWd(Img);

%% Initialize location mapping u and nodes C
u.list=linspace(1,n,n);
u.loc = reshape(linspace(1,n,n),size(Img,1),size(Img,2));
C{1}=u.list;

%% STEP I: Obtain Coarsening, block pix (input: A, u.list,alphat)
level = 10;
alphat=.1; 
alpha = 2;
Phi{1} = ones(1,n);
parents{1}=0;
parents{n}=0;
modW{1}=W{1};
Athresh = modW{1};
cImg=reshape(Img,n,1);
Cthresh = C{1};
Vol = ones(1,n);
for i =1:level
    clear Cloc
    %% Find Coarse Nodes
    [C{i+1},Cloc,Pyramid{i},parents] = blockpixpyramidrand(modW{i},Cthresh,C{i},alphat,parents,n,cImg,Vol); 
    % Cloc is the location of C{i+1} nodes in C{i}.
    % This is needed in order to find the rows of the C{i+1} nodes in W{i}
    % and thus being able to generate P{i}.
    
    %% Compute Coarse Affinity (W) and prolongation (P)
    [P{i},W{i+1},Vol] = blockweight(modW{i},C{i+1},Cloc,Vol);
    
    %% Compute Saliency
    tic
    sigma=i+1;
    Phi{i+1} = zeros(1,size(P{i},2));
    Gamma{i} = zeros(1,size(P{i},2));
    for k=1:size(P{i},2)
        Phi{i+1}(k)=sum(Phi{i}'.*P{i}(:,k));
        Gamma{i}(k)=(sum(W{i}(k,:))/(Phi{i+1}(k)^alpha))*(2^sigma);
    end
    
    % Here we use 0.0001*(average of Gamma) so that
    % the saliency exclusion does not affect early stages
    % of the coarsening process.
    clear Athresh Cthresh
    nGamma{i} = Gamma{i}/max(Gamma{i});
    threshpts{i} = find(nGamma{i}<.0001*mean(nGamma{i}));
    Cthresh = C{i+1}(threshpts{i});
    
    %% Compute G and modified A (input: Img, W, mu)
    mu=1;
    [G{i},modW{i+1}]=modifyW(Img,P,C{i+1},W{i+1},mu);
    disp(['Coarsening level ', num2str(i), ' completed'])
    
    Athresh = modW{i+1}(threshpts{i},threshpts{i});
end

%% Use all vectors as separate segments
nbSegments = size(C{level+1},2)
Seg = eye(nbSegments);
%% Apply prolongation to obtain fine solution
SegStruct{size(P,2)+1} = Seg;
[~,SegStructNum{size(P,2)+1}] = max(SegStruct{size(P,2)+1},[],2);
unitSeg = eye(size(Seg,1),size(Seg,1));
for i = size(P,2):-1:1
    SegStruct{i} = P{i}*SegStruct{i+1};
    unitSeg = P{i}*unitSeg;
    [~,SegStructNum{i}] = max(SegStruct{i},[],2);
end
Seg = SegStruct{1};
figure()
for i=1:level+1
    ptsSegStruct{i} = segtopts(C{i},SegStructNum{i},nbSegments,u);
    subplot(5,5,i)
    plotsetseg(Imgr,ptsSegStruct{i})
end

start =level-3;
for i =start:size(C,2)-1
    for j =1:size(C{i+1},2)
            [pts{i-start+1}(j,2),pts{i-start+1}(j,1)]=find(C{i+1}(j)==u.loc);
    end
end

%% Perform Ncut with constraint
%% Compute the Laplacian
for k=1:size(pts,2)
clear c d L Dinvsqrt D X newpts
    for j=1:size(pts{k},1)
        newpts{j} = pts{k}(j,:);
    end
    disp('computing B c');
    [B{k}, c] = createConstraintNew(newpts,nr,nc);
    B{k} = sparse(B{k});
    disp('computing L');
    d = sum(abs(W{1}),2);
    L = spdiags(d,0,size(W{1},2),size(W{1},2)) - sparse(W{1});
    disp('computing P');
    tp = tic;
    Dinvsqrt = 1./sqrt(d);
    D = spdiags(1./sqrt(d),0,size(W{1},2),size(W{1},2));
    A = spmtimesd(L,Dinvsqrt,Dinvsqrt);
    tol_ppm =10^(-5);
    k_ppm =10^7;
    tic
    for j=1:size(c,2)
        [x{k}(:,j),f{k}(:,j)] = uzawa(A,A,B{k},c(:,j),tol_ppm,k_ppm); 
    end 
    cputimedecoupled=toc
end

    clear bw
    thresh=.29;
    figure()
    h = tiledlayout(3,size(pts,2), 'TileSpacing', 'compact', 'Padding', 'none');
    for i=start+1:level+1
        nexttile
        plotsetseg(Imgr,ptsSegStruct{i})
    end
    for k=1:size(pts,2)
        clear X
        nsets = size(pts{k},1);
        X=x{k};
        nexttile
        [Imgt1,Imgt2,Label]=postProcesskway(X,Img,Imgr,pts,thresh);
        bw{k} = edge(Label,0.1); 
        imagesc(Label)
    end
    for k=1:size(pts,2)
        nexttile
        h = showmask(double(im2gray(Imgr)),imdilate(bw{k},ones(2,2)));
        axis off
    end
    
    % thresh=0.5;
    % for k =1:size(x,2)
    %     clear sets pts
    %     pts = ptsSegStruct{end-size(x,2)+k};
    %     sets = size(pts,2);
    %     for i=1:size(x{k},2)
    %         y{k} = x{k}
    %     end
    %     [Imgt1,Imgt2,Label]=postProcesskway(y{k},Img,Img,pts,thresh);
    %     figure()
    %     imagesc(Label)
    %     axis off
    %     truesize
    %     bw = edge(Label,0.01);
    %     figure();clf;
    %     showmask(double(im2gray(Imgr)),imdilate(bw,ones(2,2)),[0,1]);
    %     axis off;
    %     truesize
    %     colorSet = {'b','g','y','c','m'};
    %     hold on;
    %     for i = 1:sets
    %         scatter(pts{i}(:,1),pts{i}(:,2),nc/3,'o',colorSet{i},'filled');
    %     end
    %     hold off;
    % end
    % sgtitle('CNCUT for coarse level constraints')

    %% Display Fuzzy Values 
    figure()
    h = tiledlayout(3,5, 'TileSpacing', 'compact', 'Padding', 'none');
    for i=1:size(x{2},2)
        nexttile
        imagesc(reshape(x{2}(:,i),nr,nc))%,[0,1])
        colorbar
        axis off
    end
    sgtitle('CNCUT Fuzzy Values')
    
    figure()
    h = tiledlayout(1,nbSegments, 'TileSpacing', 'compact', 'Padding', 'none');
    for i=1:nbSegments
        nexttile
        imagesc(reshape(SegStruct{1}(:,i),nr,nc))
        colorbar
        axis off
    end
    sgtitle('Prolongation (Achi) Fuzzy Values')
    
