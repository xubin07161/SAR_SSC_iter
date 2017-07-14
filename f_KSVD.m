function [Dictionary,output] = f_KSVD(...
    NL_mat, Self_arr, Sigma_arr, CurPat, Par,... % 
    param)

Dictionary=param.initialDictionary;

%normalize the dictionary.
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.
totalErr = zeros(1,param.numIteration);

% the K-SVD algorithm starts here.
if Par.patnum*length(Self_arr)>param.maxBlocks
    groupNum=ceil(param.maxBlocks/Par.patnum);
    tmpIdx=randperm(length(Self_arr));
    NL_mat=NL_mat(:,tmpIdx(1:groupNum));
    Self_arr=Self_arr(:,tmpIdx(1:groupNum));
end

paraMex.L=floor(Par.patsize*Par.patsize/2);
% paraMex.eps=(param.errorGoal)^2*size(Dictionary,1);
ind_groups=int32([0]);
CoefMatrix=[];
Data=zeros(Par.patsize*Par.patsize,length(Self_arr)*Par.patnum);
for  i =  1 : length(Self_arr)                                 % For each keypatch group
    Data(:,Par.patnum*i-Par.patnum+1:Par.patnum*i) = CurPat(:, NL_mat(1:Par.patnum,i));
end

Sigma_mean=mean(Sigma_arr);
ind_groups=int32([0:Par.patnum:(length(Self_arr)-1)*Par.patnum]);

for iterNum = 1:param.numIteration
    CoefMatrix=zeros(param.K,length(Self_arr)*Par.patnum);
    paraMex.eps=(Sigma_mean*Par.patsize)^2*0.5;
    CoefMatrix=mexSOMP(Data,Dictionary,ind_groups,paraMex);
%     for  i      =  1 : length(Self_arr)                                 
%         paraMex.eps=(Sigma_arr(Self_arr(i))*Par.patsize)^2*0.5;
%         CoefMatrix(:,Par.patnum*i-Par.patnum+1:Par.patnum*i)=...
%             mexSOMP(Data(:,Par.patnum*i-Par.patnum+1:Par.patnum*i),Dictionary,ind_groups,paraMex);
%     end
    for kk=1 : full(max(sum(isnan(CoefMatrix))))
        CoefMatrix(isnan(CoefMatrix))=0;
    end
    
%     replacedVectorCounter = 0;
	rPerm = randperm(size(Dictionary,2));
    for j = rPerm
        [Dictionary(:,j),CoefMatrix] = I_findBetDictElem(Data,...
            Dictionary,j,CoefMatrix);
    end

    if (iterNum>1 & param.displayProgress)
        disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',...
            num2str(length(find(CoefMatrix))/size(Data,2))]);        
    end    
    Dictionary = I_clearDictionary(Dictionary,CoefMatrix(1:end,:),Data);
end

output.CoefMatrix = CoefMatrix;
% Dictionary = Dictionary;
%==========================================================================
%  findBetterDictionaryElement

function [betDictElem,CoefMatrix] = I_findBetDictElem(Data,Dictionary,j,CoefMatrix)
relevantDataIndices = find(CoefMatrix(j,:)); 
if (length(relevantDataIndices)<1) 
%     [~,i]=max(sum((Data-Dictionary*CoefMatrix).^2));
    i=1;
    betDictElem = Data(:,i);
    betDictElem = sign(betDictElem(1))*betDictElem./sqrt(betDictElem'*betDictElem);
    CoefMatrix(j,:) = 0;
%     NewVectorAdded = 1;
    return;
end

tmpCoefMatrix = CoefMatrix(:,relevantDataIndices); 
tmpCoefMatrix(j,:) = 0;
errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix); 

% [betDictElem,singularValue,betaVector] = svds(errors,1);
[betDictElem,singularValue,betaVector] = f_svdsMax(errors);
CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';

%==========================================================================
%  I_clearDictionary
function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)
T2 = 0.99;
T1 = 3;
K=size(Dictionary,2);
Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms
tmp=logical(diag(ones(1,K)));
G=Dictionary'*Dictionary; G(tmp)=0;%G = G-diag(diag(G));
for jj=1:K
    if max(G(jj,:))>T2 || length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,
        [~,pos]=max(Er);
        Er(pos(1))=0;
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary;
        G(tmp)=0;%G = G-diag(diag(G));
    end;
end;

