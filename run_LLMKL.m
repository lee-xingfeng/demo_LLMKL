% For Knowledge-Based Systems.
clc;
warning off

ds = {'AR_840n_768d_120c_uni'};
iData = 1;
data_dir = fullfile(pwd);
dataset = ds{iData}
data_file = fullfile([dataset, '.mat']);
kernel_file = fullfile([dataset, '_allkernel.mat']);
load(data_file)
load(kernel_file)

X=X';
X = NormalizeFea(X,0);

lambda1=1e-2;
lambda2=1e3;
lambda3=3;
lambda4=1e2; 



mu = 1.2; 
param.lambda1=lambda1;
param.lambda2=lambda2;
param.lambda3=lambda3;
param.lambda4=lambda4;
param.eta = 20;
param.mu = mu;
fprintf('%.5f %.4f %.4f %.4f',lambda1,lambda2,lambda3,lambda4);
[result,Z,E]= LLMKL(K,X,y,param);
result
imagesc(Z);colormap(jet);axis off;