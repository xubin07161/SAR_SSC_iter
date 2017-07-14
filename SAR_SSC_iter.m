function [img_est] = SAR_SSC_iter(img, varargin)
%% parameters
Par.isOri = 0;                              % ��ԭͼ��

for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'ENL'))
        Par.ENL = varargin{argI+1};       % �����ı�׼��
    end
    if (strcmp(varargin{argI}, 'ori'))
        Par.img_ori = varargin{argI+1};     % ԭͼ��������psnr��ssim
        Par.isOri = 1;                      % ��ԭͼ��
    end
    if (strcmp(varargin{argI}, 'patnum'))
        Par.patnum = varargin{argI+1};       
    end
    if (strcmp(varargin{argI}, 'patsize'))
        Par.patsize = varargin{argI+1};       
    end
    if (strcmp(varargin{argI}, 'lamada'))
        Par.lamada = varargin{argI+1};       
    end
    if (strcmp(varargin{argI}, 'SW'))
        Par.SearchWin = varargin{argI+1};       
    end
    if (strcmp(varargin{argI}, 'delta'))
        Par.delta = varargin{argI+1};       
    end
end

Par.SearchWin =   40;                       % nonlocal��������Ϊ2*SearchWin+1
Par.delta     =   0.03;                    % ÿ�ε������������
Par.Innerloop =   1;                        % ���ڼ��ٿ�ƥ��Ĵ���
Par.learn = 1;                              % �Ƿ�ѧϰ�ֵ�

if Par.ENL<1.5
    Par.patnum = 15;
else
    Par.patnum = 10;
end

if Par.ENL<=1.5
    Par.patsize = 9;
elseif Par.ENL<=3
    Par.patsize = 8;
elseif Par.ENL<=6
    Par.patsize = 7;
else
    Par.patsize = 6;
end

Par.Iter          =   6;
Par.lamada        =   0.54;
if prod(size(img))<=300^2
    Par.step = 2;                           % �ο��鲽��
elseif prod(size(img))<=600^2
    Par.step = 3;
else
    Par.step = 4;
end

Par.K=max(256,4*Par.patsize*Par.patsize);
% Par.K=4*Par.patsize*Par.patsize;
%% �����任
Par.sigma=sqrt(psi(1,Par.ENL));

img(isnan(img))=0; % consider the NAN case
img=abs(img); % consider the nonnegative case
img(img==0)=10^(-8); % consider the zero case

img=log(img)-(psi(0,Par.ENL)-log(Par.ENL));
%% Ѱ��ÿһ���ο����������index
img_est = img;
[Height, Width] = size(img);
% ����ÿһ���ο����������±�(Neighbor_arr)������(Num_arr)��ÿһ���ο�����±�(Self_arr)
[Neighbor_arr Num_arr Self_arr] = NeighborIndex(img, Par);        
NL_mat = zeros(Par.patnum,length(Num_arr));          % ���nonlocalѡȡ��Ŀ��±�

%%
for iter = 1 : Par.Iter
    img_est = img_est + Par.delta*(img - img_est);
%     [CurPat Sigma_arr] = Im2Patch( img_est, img, Par );
    CurPat = im2colstep(img_est,[Par.patsize Par.patsize]);
    N_CurPat = im2colstep(img,[Par.patsize Par.patsize]);
    Sigma_arr = Par.lamada*sqrt(abs(repmat(Par.sigma^2,1,size(CurPat,2))-...
        mean((N_CurPat-CurPat).^2)));
    
    if(iter==1)
        Sigma_arr = Par.sigma * ones(size(Sigma_arr));                      
    end
    tmp = mean(Sigma_arr)
    
    if (mod(iter-1,Par.Innerloop)==0)
        NL_mat = BlockMatch(CurPat, Par.patnum, Neighbor_arr, Num_arr, Self_arr);
    end
    if (mod(iter-1,2)==0)||(iter==Par.Iter)
        Dict = f_LearnDict( NL_mat, Self_arr, Sigma_arr, CurPat, Par );
    end
    [EPat W] = PatEst( NL_mat, Self_arr, Sigma_arr, CurPat, Par, Dict );
    
%     img_est =  Patch2Im( EPat, W, Par.patsize, Height, Width );
    img_est = col2imstep(EPat,[Height Width],[Par.patsize Par.patsize]);
    W_Img = col2imstep(W,[Height Width],[Par.patsize Par.patsize]);
    img_est  =  img_est./(W_Img+eps);
    
    if Par.isOri==1
        fprintf( 'Iter = %2.3f, PSNR = %2.3f, SSIM=%2.4f \n', iter,...
            f_PSNR(sqrt(Par.img_ori), sqrt(exp(img_est))),...
            ssim(sqrt(Par.img_ori)*255, sqrt(exp(img_est))*255));
    end
%     figure,imshow(sqrt(exp(img_est)));
end
img_est=exp(img_est);
end