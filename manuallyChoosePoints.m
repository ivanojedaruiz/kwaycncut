function pts = manuallyChoosePoints(imgName,numLabelSets)
bb = [0, 0.4470, 0.7410];
oo = [0.8500 0.3250 0.0980];
yy = [0.9290 0.6940 0.1250];
pp = [0.4940 0.1840 0.5560];
gg = [0.4660 0.6740 0.1880];
pp = [0.3010 0.7450 0.9330];
rr = [0.6350 0.0780 0.1840];
colorSet = {bb,'b','g','y','c','m','r','k',oo,yy,pp,gg,pp,rr'};
     hold on;
    for i = 1:numLabelSets
        [x, y] = ginput;
        [m,n]=size(x);
        if m>0
            pts{i} = [x,y];
            scatter(pts{i}(:,1),pts{i}(:,2),300,"MarkerFaceColor",rand(1,3));
        end
    end
    hold off;
end