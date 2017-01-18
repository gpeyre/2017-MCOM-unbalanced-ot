%%
% Test for Sinkhorn between two empirical distributions.

addpath('toolbox/');

% helper functions
plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

rep = 'results/sinkhorn/';
[~,~] = mkdir(rep);

%%
% Load input points.

% number of points.
N = 200;

% test name
name = 'randompointsimplex';
name = 'pointclouds-2d';

switch name
    case 'pointclouds-2d'
        rand('state', 123);
        % samples in a square
        x = rand(2,N)-.5;
        % samples in an annulus
        theta = 2*pi*rand(1,N);
        r = .8 + .2*rand(1,N);
        y = [cos(theta).*r; sin(theta).*r];  
        % uniform distributions
        mu = ones(N,1)/N;
        nu = ones(N,1)/N;
        % cost: squared L2 norm
        x2 = sum(x.^2,1); y2 = sum(y.^2,1);
        c = repmat(y2,N,1)+repmat(x2.',1,N)-2*x.'*y;
        
    case 'randompointsimplex'
        addpath('toolbox_shortest_paths/');
        mu=randomPointinDSimplex(N);
        nu=randomPointinDSimplex(N);
        % random cost matrix
        % Generate Random graph
        sparsity=.05;
        M = sprand(N,N,sparsity);
        M = M-diag(diag(M));
        M = (M+M')/2; 
        c = floyd_warshall_all_sp(M);   
        
    otherwise
        error('Unknown');        
end

if exist('x');
    % display clouds
    clf; hold on;
    plotp(x, 'b');
    plotp(y, 'r');
    axis('off'); axis('equal');
end

%%
% Parameters.

epsilon = .5*(.1)^2;
options.niter = 500;

%%
% Test Sinkhorn for different values of extrapolation. tau.

tau_list = [-.85 -.8 -.5 -.3 0 .3 ];
err = []; lgd = {};
for k=1:length(tau_list)
    options.tau = tau_list(k);
    [u,v,gamma,Wprimal,Wdual,err(:,k)] = sinkhorn_log(mu,nu,c,epsilon,options);
    lgd{k} = ['\tau=' num2str(tau_list(k))];
end

%%
% Display error decays. 

E = log10(err/err(1));
clf; hold on;
for k=1:length(tau_list)
    t = (k-1)/(length(tau_list)-1);
    plot(E(:,k), 'color', [t 0 1-t]);
end`
legend(lgd);
saveas(gcf, [rep name '-conv-eps' num2str(round(epsilon*10000)) '.png'], 'png');
