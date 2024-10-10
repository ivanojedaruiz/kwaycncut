function [C,Cloc,Pyramid,parents]=blockpixpyramidrand(A,oldC,Call,alphat,parents,totalpixels,Img,Vol)
    
    n = size(A,2)
    %sort randomly
    %[sd,r]=sort(rand(1,size(A,2)));
    
    %sort by coupling
    %[sd,r]=sort(full(sum(A,1)));
    
    %sort by intensity
%     cImg = reshape(Img,n,1);
%     size(Call)
%     size(oldC)
%     [sd,r]=sort(full(cImg(oldC)));  
    
    %sort by volume
    [sd,r]=sort(full(Vol),'descend');
    
    C = Call(r(1));
    Cloc = r(1);
    Pyramid{1} = 1;
    %Pyramid{totalpixels}= 0;
    ASum = alphat*sum(A,1);
    inCthresh = [];
    for i = 2:n
        inCthresh = [inCthresh, size(find(Call(r(i))==oldC),2)];
        if max(A(r(i),Cloc)) < ASum(r(i)) || inCthresh(i-1)==1     
        %if max(A(r(i),Cloc)) < ASum(r(i)) && inCthresh(i-1)==1
            Cloc = [Cloc,r(i)];
            C = [C,Call(r(i))];
            Pyramid{r(i)} = Call(r(i));
            parents{Call(r(i))} = 0;
        else
            Pyramid{r(i)} = Call(find(A(r(i),:)./max(A(r(i),:))>.2)); %only for childs
            parents{Call(r(i))} = Pyramid{r(i)};
        end
        if mod(i,1000)==0
            disp(['Block Pixel iteration: ',num2str(i)])
        end
    end
    Cloc = sort(Cloc);
    %max(Cloc)
    C = sort(C);
%     max(C)
%     Call(Cloc)==C
%     pause
end