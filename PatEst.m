function [ EPat W ] = PatEst( NL_mat, Self_arr, Sigma_arr, CurPat, Par, Dict )

EPat = zeros(size(CurPat));
W    = zeros(size(CurPat));
Data=zeros(Par.patsize*Par.patsize,length(Self_arr)*Par.patnum);
M_Data=zeros(Par.patsize*Par.patsize,length(Self_arr)*Par.patnum);
for i = 1 : length(Self_arr)
    Temp = CurPat(:, NL_mat(1:Par.patnum,i));
    M_Temp  =   repmat(mean(Temp,2),1,Par.patnum);
    Temp    =   Temp-M_Temp;
    Data(:,Par.patnum*i-Par.patnum+1:Par.patnum*i) = Temp;%./Sigma_arr(Self_arr(i));
    M_Data(:,Par.patnum*i-Par.patnum+1:Par.patnum*i) = M_Temp;
end
% errT=0.42*mean(Sigma_arr);
errT=0.68*mean(Sigma_arr);
param.L=floor(Par.patsize*Par.patsize/2);
param.eps=errT^2*Par.patsize*Par.patsize;
ind_groups=int32([0:Par.patnum:(length(Self_arr)-1)*Par.patnum]);
Coefs=mexSOMP(Data,Dict,ind_groups,param);
for kk=1 : full(max(sum(isnan(Coefs))))
    Coefs(isnan(Coefs))=0;
end
E_Data=Dict*Coefs+M_Data;
for i = 1 : length(Self_arr)
    EPat(:,NL_mat(1:Par.patnum,i))  = EPat(:,NL_mat(1:Par.patnum,i))+...
        E_Data(:,Par.patnum*i-Par.patnum+1:Par.patnum*i);%.*Sigma_arr(Self_arr(i));
    W(:,NL_mat(1:Par.patnum,i))     = W(:,NL_mat(1:Par.patnum,i))+...
        ones(Par.patsize*Par.patsize,size(NL_mat(1:Par.patnum,i),1));
end
end
