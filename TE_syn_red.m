function [drivers_red,drivers_syn,g_red,g_syn]=TE_syn_red(x,i,j,p,d)
%transfer entropy i -> j 
%p=order of the model
%d=delay between the variables
xxx=size(x);
drivers_red=[i j];
ind=setdiff(1:xxx(2),[i j]); %the varying set of condition variables
m=length(ind);
% GC_X = copnorm(x(:,i));
% GC_Y = copnorm(x(:,j));
% [best_L, best_delay, ~, ~, ~] = optimize_te_aic_delay(GC_X, GC_Y, [], orders, delays); % best order and 
TE= CTE_iteration(x(:,[i j]),p,d,i,j); %
g=mean(TE);
g_red(1)=g;g_syn(1)=g;

%find the subset  Zm (drivers_red) that minimize the conditional TE: TE(x→y|Zm)
icont=1;

%save the best orders and delays for red and syn iterations, repectively (for surrogate TE generation)
% best_order_dealys_red = ones(2,m);
% best_order_dealys_syn = ones(2,m);
while m > 0
    icont=icont+1;
%     drivers_red_update = [drivers_red ind(1)];
%     condition_variables_dt = setdiff(drivers_red_update,[i j]);
%     GC_Z = copnorm(x(:,condition_variables_dt));
%     [best_L_red, best_delay_red, ~, ~, ~] = optimize_te_aic_delay(GC_X, GC_Y, GC_Z, window_lengths, delays);
%     best_order_dealys_red(1,m) = 
%     best_order_dealys_red(2,m) = 
    for h=1:m
        dt=[drivers_red ind(h)];
        TE= CTE_iteration(x(:,dt),p,d,i,j);
        TE_mean=mean(TE);
        deltas(h)=TE_mean;
    end
    [g_red(icont), ii]=min(deltas);
    drivers_red=[drivers_red ind(ii)];
    ind=setdiff(ind,ind(ii));m=length(ind);
    clear deltas;
end
%%find the subset  ZM (drivers_syn) that maximize the conditional TE: TE(x→y|ZM)
drivers_syn=[i j];
ind=setdiff(1:xxx(2),[i j]);
n=length(drivers_syn);
m=length(ind);
icont=1;
while m > 0
    icont=icont+1;
    for h=1:m
        dt=[drivers_syn ind(h)];
        TE=CTE_iteration(x(:,dt),p,d,i,j);
        TE_mean=mean(TE);
        deltas(h)=TE_mean;
    end
    [g_syn(icont), ii]=max(deltas);
    drivers_syn=[drivers_syn ind(ii)];
    ind=setdiff(ind,ind(ii));m=length(ind);
    clear deltas;
end
