function P_init=initial_perm(n,initial_ball,P_0)
    P_init=P_0;                             %latent (g.t) permutation
   %P_rnd = eye(n);                        %if we assume g.t=Id
    p1=randperm(n,initial_ball);           %choose 'initial_ball' indices at random
    p2=p1(randpermfull(length(p1)));       %generate a derrangement on those indices
    P_init(:,p2) = P_init(:,p1);