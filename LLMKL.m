function [result,Z,E]= LLMKL(H,X,y,param)
nCluster=length(unique(y));

nn=length(y);
m = size(H,3);

AA = zeros(1,m);
BB = ones(1,m);
CC = 1;
lb =zeros(1,m);

K = zeros(nn);
g = ones(1,m)/m;
for j = 1:m
    K = K + g(j)*H(:,:,j);
end

E = zeros(nn);
Z=eye(nn);
epsilon = 1e-5; maxIter = 1e3;  max_rho = 1e10;
lambda1 = param.lambda1;
lambda2 = param.lambda2;
lambda3 = param.lambda3;
lambda4 = param.lambda4;
mu = param.mu;
eta = param.eta;
Y1=zeros(nn);

distX = L2_distance_1(X,X);
DD = lambda3*distX;

BTB=K;
obj=[];
for iter = 1:maxIter
    %Update Z
     Z= (BTB+eps*eye(nn))\(BTB-DD);
     Z = Z - diag(diag(Z));
    for ii = 1:size(Z,2)
        idx = 1:size(Z,2);
        idx(ii) = [];
        Z(ii,idx) = EProjSimplex_new(Z(ii,idx));
    end
    
    %Update H
    HI=zeros(nn);
    for j = 1:m
        HI = HI + g(j)*H(:,:,j);
    end
    Kold=K;
    K=((mu+lambda2)*eye(nn))\((mu*(BTB+E)-Y1)+lambda2*HI);
    K(find(K<0))=0;
    
    %Update B
    G=K-E+Y1/mu;
    GWave=G-(eye(nn)-2*Z'+Z*Z')/(2*mu);
    B = solveB(GWave,mu/lambda1);
    BTB=B'*B;
    
    % Update E
    tmp = K - BTB + Y1/mu;
    E = max(0,tmp - lambda4/mu)+min(0,tmp + lambda4/mu);
    
    %update g
    for a = 1:m
        for b = 1:m
            M(a,b) = trace( H(:,:,a)*H(:,:,b));
        end
        AA(a)= 2*lambda2*trace(K*H(:,:,a));
    end
    options = optimoptions('quadprog','Algorithm','interior-point-convex','Display','none');
    g = quadprog( 2*M*lambda2,-AA,[],[],BB,CC,lb,[],[],options);
    
    leq = K - BTB - E;
    Y1 = Y1 + mu*leq;
    mu = min(max_rho,mu*eta);
    %obj(iter)=norm(K-Kold,'fro');  
    %obj(iter)=0.5*trace((eye(nn)-2*Z'+Z*Z')*BTB)+lambda1*rank(B)+0.5*lambda2*norm(K-HI,'fro')+lambda3*trace(DD'*Z);

    if((iter>5)&(norm(K-Kold,'fro') < norm(Kold,'fro') * epsilon))
        break
    end
end

Z = Z - diag(diag(Z));
Z = abs(Z);
Z=(Z+Z')/2;
for ij=1:5 % 5 -> 20
    CKSym = BuildAdjacency(thrC(Z,0.7));
    grps = SpectralClustering(CKSym,nCluster);
    grps = bestMap(y,grps);
    res(ij,:) =  ClusteringMeasure(grps,y);
end
result(1,1)=mean(res(:,1));
result(2,1)=mean(res(:,2));
result(3,1)=mean(res(:,3));
result(1,2)=std(res(:,1));
result(2,2)=std(res(:,2));
result(3,2)=std(res(:,3));
end
