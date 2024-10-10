function pts = segtopts(nodes,seg,nbSegments,u)
coarseset = [];
for i = 1:nbSegments
    coarseset{i} = nodes(find(seg==i));
end

%% Use coarse nodes from pyramid sets
k=1;
for i =1:size(coarseset,2)
    for j =1:size(coarseset{i},2)
        [pts{k}(j,2),pts{k}(j,1)]=find(coarseset{i}(j)==u.loc);
    end
    k=k+1;
end