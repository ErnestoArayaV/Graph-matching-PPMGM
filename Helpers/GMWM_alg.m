%Greedy algorithm GMWM, maximum version
%Inputs: 
%M           <------ matrix of 'costs' 
%param       <------ negative parameter to 'emulate' the erase step

function X=GMWM_alg(M,param)
    n=length(M);
    X=zeros(n,n);
    for i=1:n
        [~,I] = max(M(:));
        [I_row, I_col] = ind2sub(size(M),I);
        X(I_row, I_col)=1;
        M(I_row,:)=param*ones(1,n);
        M(:,I_col)=param*ones(n,1);
    end 
        