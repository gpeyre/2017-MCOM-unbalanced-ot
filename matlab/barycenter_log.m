function [nu,gamma] = barycenter_log(mu,c,epsilon,w,options)

% barycenter_log - stabilized sinkhorn for barycenters over log domain with acceleration
%
%   [nu,gamma] = barycenter_log(mu,c,epsilon,w,options);
%
%   mu are marginals.
%   c is cost
%   epsilon is regularization
%
%   options.niter is the number of iterations.
%   options.tau is an avering step.
%       - tau=0 is usual sinkhorn
%       - tau<0 produces extrapolation and can usually accelerate.
%
%   options.rho controls the amount of mass variation. Large value of rho
%   impose strong constraint on mass conservation. rho=Inf (default)
%   corresponds to the usual OT (balanced setting).
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter', 1000);
tau  = getoptions(options, 'tau', -.5);
verb = getoptions(options, 'verb', 1);
rho = getoptions(options, 'rho', Inf);

tau_v = getoptions(options, 'tau_v', -.5);

% number of input marginals
L = length(mu);

lambda = rho/(rho+epsilon);
if rho==Inf
    lambda=1;
end

N = size(mu{1},1);
H1 = ones(N,1);
H2 = ones(N,1);

ave = @(tau, u,u1)tau*u+(1-tau)*u1;

lse = @(A)log(sum(exp(A),2));
M = @(u,v)(-c+u*H2'+H1*v')/epsilon;

% kullback divergence
H = @(p)-sum( p(:).*(log(p(:)+1e-20)-1) );
KL  = @(h,p)sum( h(:).*log( h(:)./p(:) ) - h(:)+p(:) );
KLd = @(u,p)sum( p(:).*(exp(-u(:))-1) );
dotp = @(x,y)sum(x(:).*y(:));

err = [];
for k=1:L
    u{k} = zeros(N,1);
    v{k} = zeros(N,1);
end
disp_rate = getoptions(options,'disp_rate', 1);

for i=1:niter
    if verb==1
        progressbar(i,niter);
    end
    % update u
    for k=1:L
        u{k} = ave(tau, u{k}, ...
            lambda*epsilon*log(mu{k}) - lambda*epsilon*lse( M(u{k},v{k}) ) + lambda*u{k} );
        if nargin>1
            gamma{k} = exp(M(u{k},v{k}));
        end
    end
    % update nu
    Lnu = mu{1}*0;
    for k=1:L
        LSEv{k} = lse( M(u{k},v{k})' );
        Lnu = Lnu + w(k)*( LSEv{k}  ); % + v{k}/epsilon
    end
    % update v
    for k=1:L
        v{k} = ave(tau_v, v{k}, ...
            epsilon*Lnu - epsilon*LSEv{k} + v{k} );
        if nargin>1
            gamma{k} = exp(M(u{k},v{k}));
        end
    end
    if mod(i,disp_rate)==1
        nu = exp(Lnu);
        clf; hold on;
        plot(mu{1}, 'r');
        plot(mu{2}, 'b');
        plot(nu, 'k');
        axis tight; drawnow;
    end
end


end