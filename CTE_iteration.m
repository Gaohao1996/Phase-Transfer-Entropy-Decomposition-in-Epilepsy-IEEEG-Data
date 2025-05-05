function te=CTE_iteration(x,p,d,i,j)
% x data S (samples) x n (variables) dati 
% i driver j target
% p order of the model
[~, n]=size(x);
ind=setdiff(1:n,[i j]);
%% X,Y and Z variables
X = x(:,i);
Y  = x(:,j);
Z = x(:,ind);

% %% Select best model order and delay for CTE calculation based on AIC
% GC_X = copnorm(X);
% GC_Y = copnorm(Y);
% GC_Z = copnorm(Z);
% [best_L, best_delay, ~, ~, ~] = optimize_te_aic_delay(GC_X, GC_Y, GC_Z, p, d);

te= conditional_TE_delay(X,Y,Z,p,d);


% NS_p=NS-p;
% y=zeros(NS_p,1);
% Y=zeros(NS_p,p);
% X=zeros(NS_p,p);
% Z=zeros(NS_p,p*(n-2));
% y=x(p+1:NS,j);
% for h=1:p
% Y(:,h)=x(h:NS_p+h-1,j);
% X(:,h)=x(h:NS_p+h-1,i);
%     for m=1:n-2
%         Z(:,(n-2)*(h-1)+m)=x(h:NS_p+h-1,ind(m));
%     end
% end
% [~,~,R] = regress(y,[Y Z]);v=var(R);
% [~,~,R] = regress(y,[Y X Z]);v1=var(R);
% te=.5*log(v/v1);