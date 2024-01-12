%Matching by the projected power method
%Inputs: 
%A,B <----- Square matrices to be matched.
%in_point<----- Initialization matrix (should be a permutation matrix).
%maxiter<----- maximum PPM iterations
%Initialization: 

function P = matching_proj_it_in(A, B,maxiter,in_point) %[P,inter,is_subset] = matching_proj_it_in(A, B,maxiter,in_point,P_gt)
%
%n = size(A, 1);
%Alig= zeros(maxiter,1);

X=A*in_point'*B;
%indaux=find(all((P_gt-in_point)==0,2));
for i=1:maxiter
    P=GMWM_alg(X',-2000);
    %Alig(i)=(1/n)*sum(sum(X.*(grounT)));
    X=A*P'*B;
end

P=GMWM_alg(X',-2000);
%[M1r, M1c] = linear_sum_assignment(-X');
%P = full(sparse(M1r, M1c, 1, n, n));
%P=P';

%X = reshape(X, [n n]);

%P=X';