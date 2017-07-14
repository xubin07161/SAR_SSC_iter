function [ Dict ] = f_LearnDict( NL_mat, Self_arr, Sigma_arr, CurPat, Par )

Pn=ceil(sqrt(Par.K));
DicOri = Generate_DCT_Matrix(Par.patsize,Pn);
DicOri = DicOri./repmat(sqrt(sum(DicOri.^2)),[size(DicOri,1),1]);
% Par.learn=1;
Dict=DicOri;

if Par.learn==0,  Dict=DicOri; end
if Par.learn==1
    par2.K = Par.K;
    par2.numIteration = 5;
    par2.errorFlag = 1;
    errT=Sigma_arr(1);
    par2.errorGoal = errT^2*Par.patsize*Par.patsize*0.6;
    par2.preserveDCAtom = 0;
    par2.initialDictionary = DicOri(:,1:par2.K);
    par2.InitializationMethod =  'GivenMatrix';
    par2.maxBlocks=5000;
    par2.displayProgress = 1;    
    [Dict,~] = f_KSVD(NL_mat, Self_arr, Sigma_arr, CurPat, Par, par2);    
end

end
