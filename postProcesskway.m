function [Imgt1,Imgt2,Label] = postProcesskway(X,Img,Imgo,pts,thresh)

[nr,nc,~] = size(Img);
Seg = zeros(nr*nc,1); 
%% Exclusion condition
for i =1:nr*nc%size(C,2)
    [~,id]=max(X(i,:)); %rescale
    Seg(i)=id;
end 

%% Use kmeans
% size(X)
%  Seg = kmeans(X,size(X,2));
% unique(Seg)
% pause 
%% Segmentation by Threshold
% for i=1:nr*nc
%     id = find(X(i,:)>thresh);
%     if size(id,2)==1
%         Seg(i)=id;
%     end
%     if size(id,2)>1
%         Seg(i)=id(1)+size(X,2);
%     end
% end
       
Label = reshape(Seg,nr,nc);

 % imagesc(Label);
 % pause
% axis off;

Labelb = Label(1,1);
Imgt1 = Imgo;
for i = 1:nr
    for j = 1:nc
        if(Label(i,j)==Labelb)
            Imgt1(i,j,:) = 0;
        end
    end
end

Imgt2 = Imgo;
for i = 1:nr
    for j = 1:nc
        if(Label(i,j)~=Labelb)
            Imgt2(i,j,:) = 0;
        end
    end
end

%% Plot result of two partitions
% figure();
% colormap(gray);
% subplot(1,2,1);imagesc(Imgt1);title('First Partition');%colormap(gray);clf;
% subplot(1,2,2);imagesc(Imgt2);title('Second Partition');%colormap(gray);clf;
% truesize([400,400])
% pause 

% axis off;
%set(gca,'position',[0 0 1 1],'units','normalized');