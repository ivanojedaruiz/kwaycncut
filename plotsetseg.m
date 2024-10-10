function [] = plotsetseg(Imgo,pts)

h = imagesc(Imgo);
colormap(gray)
hold on;

colorSet = {'co','go','yo','co','mo','b*','g*','y*','c*','m*',...
    'b+','g+','y+','c+','m+','b.','g.','y.','c.','m.','bx','gx','yx','cx','mx'};
% set(gca,'position',[0 0 1 1],'units','normalized');
colorSet
for i =1:size(pts,2)
    h = plot(pts{i}(:,1),pts{i}(:,2),colorSet{1},'MarkerSize',5);
    set(h,'linewidth',5);
end
hold off;
%truesize([2000 1000])
axis off;