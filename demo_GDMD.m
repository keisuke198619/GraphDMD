% demo_GDMD.m
% Keisuke Fujii
 
clear; close all

% You need to download Tensor Train (TT) Toolbox :
% https://github.com/oseledets/TT-Toolbox
% Please insert the path of the folder below:
path_TT = '..\TT-Toolbox';
addpath(genpath(path_TT));

% Generate Training/Test Data
params.noise_level = 0.01;
params.lt = 20; % time length [frame]
params.T = 20; % time length [s]
nTask = 10 ; % number of tasks
Asize = 64 ; % adjacent matrix size
lam(1) = 0.99 ; % eigenvalue 1
lam(2) = 0.9 ; % eigenvalue 2
eps_ttdmd = 1e-1 ; % tolerance of TT-decomposition and basic DMD

% base adjacency matrix A1 and A2
sigma = Asize*1 ;
Sigma = [sigma 0; 0 sigma];
x1 = Asize-1:-1:0; x2 = 0:Asize-1;
[X1,X2] = meshgrid(x1,x2);
A0 = 200*mvnpdf([X1(:) X2(:)],[0 0],Sigma);
A0 = reshape(A0,length(x2),length(x1));
A1 = A0+A0' ;
A2 = rot90(A0,1)+rot90(A0,3);
for a = 1:Asize
    A1(a,a) = 0 ;
    A2(a,a) = 0 ;
end

% create sequences of adjacent matrices
t = linspace(0,params.T,params.lt);
dt = t(2) - t(1);
rs = randn(size(A1,1),size(A1,2),nTask,length(t)) ; % independent of task
for n = 1:nTask
    for tt = 1:length(t)
        tmp = A1*(lam(1).^t(tt)) + A2*(lam(2).^t(tt)) + rs(:,:,n,tt)*params.noise_level ;
        x_t{n}(:,:,tt) = tmp ; % for tensor DMD
        x_tm{n}(:,tt) = tmp(:) ; % for basic DMD
    end
end

%% DMD and tensor DMD

r = 2 ;
% perform basic DMD-------------------
for t = 1:nTask
    X = x_tm{t}(:,1:end-1) ;
    Y = x_tm{t}(:,2:end) ;
    [sLambda{t},sPhi{t}] = DMD(X,Y,eps_ttdmd,dt) ; 
    slam(t,1:r) = sLambda{t}(1:r);
    for rr = 1:r
        mPhiSM(:,rr,t) = sPhi{t}(:,rr) ;
        mPhiS(:,:,rr,t) = reshape(mPhiSM(:,rr,t),Asize,Asize);
    end
end
disp('DMD is done')
mPhiS = mean(abs(mPhiS),4) ;

% perform tensor DMD------------------
for t = 1:nTask
    ttX = tt_tensor(x_t{t}(:,:,1:end-1),eps_ttdmd) ;
    ttY = tt_tensor(x_t{t}(:,:,2:end  ),eps_ttdmd) ;
    [Phi{t},Lambda{t}] = tt_dmd(ttX,ttY,dt) ; 
    % plot(real(Lambda),imag(Lambda),'ko') ;
    lamT(t,1:r) = diag(Lambda{t}(1:r,1:r));
    mPhiTT(:,:,1,t) = Phi{t}(:,:,1) ;
    mPhiTT(:,:,2,t) = Phi{t}(:,:,2) ;
    disp(['TDMD ',num2str(t),' is done'])
end
mPhiT = mean(abs(mPhiTT),4) ;

% Figure Phi
range0 = [0 0.05] ;
range1 = [0 0.01] ;
range2 = [0 0.01] ;
figure(1)
subplot(2,3,1); imshow(A1,range0);
subplot(2,3,2); imshow(mPhiS(:,:,1),range1);
subplot(2,3,3); imshow(mPhiT(:,:,1),range2);
subplot(2,3,4); imshow(A2,range0);
subplot(2,3,5); imshow(mPhiS(:,:,2),range1);
subplot(2,3,6); imshow(mPhiT(:,:,2),range2);

% Table w
Tab(1,:) = lam ;
Tab(2,:) = mean(slam,1) ;
Tab(3,:) = mean(lamT,1) ;