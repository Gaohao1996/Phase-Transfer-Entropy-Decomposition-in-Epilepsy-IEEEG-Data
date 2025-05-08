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
te= conditional_TE_delay(X,Y,Z,p,d);

