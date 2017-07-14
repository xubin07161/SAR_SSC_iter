function  [Index] = BlockMatch(X, patNum,Neighbor_arr,Num_arr, SelfIndex_arr)
Index   =  zeros(patNum,length(Num_arr));
Dist=c_CalDist(X,Neighbor_arr,Num_arr, SelfIndex_arr);
Dist_max=max(Dist(:));
Dist(Dist==-1)=Dist_max*2;
idx = f_FindMinK(Dist,patNum);
for  i  =  1 : length(Num_arr)
    Index(:,i)=Neighbor_arr(idx(:,i),i);
end
end

function [idx] = f_FindMinK(X,k)
[M,N]=size(X);
idx=zeros(k,N);
X_max=max(X(:));
for j=1 : k-1
    [~,idx_new]=min(X,[],1);
    X(idx_new+[0:M:M*(N-1)])=X_max;
    idx(j,:)=idx_new;
end
[~,idx_new]=min(X,[],1);
idx(k,:)=idx_new;
end