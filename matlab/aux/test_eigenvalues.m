%%
% test for eigenvalues of Sinkhorn matrix.

addpath('toolbox/');

% helper functions
plotp = @(x,col)plot(x(1,:)', x(2,:)', 'o', 'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', col, 'LineWidth', 2);

rep = 'results/';
[~,~] = mkdir(rep);

% number of points.
N = 200;

% test name
name = 'pointclouds-2d';
name = 'gaussians';


switch name
    case 'gaussians'
        x = (0:N-1)/N;
        y = x;
        gaussian = @(m,s)exp( -(x'-m).^2/(2*s^2) );
        normalize = @(x)x/sum(x(:));
        vmin = .05;
        mu = normalize(vmin + gaussian(.2,.05));
        nu = normalize(vmin + gaussian(.7,.1));
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
end


% cost: squared L2 norm
x2 = sum(x.^2,1); y2 = sum(y.^2,1);
c = repmat(y2,N,1)+repmat(x2.',1,N)-2*x.'*y;

% display clouds
if strcmp(name, 'pointclouds-2d')
    % display clouds
    clf; hold on;
    plotp(x, 'b');
    plotp(y, 'r');
    axis('off'); axis('equal');
else
    plot([mu nu]);
end

epsilon = (.015)^2;
options.tau = 0;
options.niter = 8000;

[u,v,gamma,Wprimal,Wdual,err] = sinkhorn_log(mu,nu,c,epsilon,options);

clf;
plot(log10(err)); axis tight;

H = diag(1./mu)*gamma*diag(1./nu)*gamma';

[U,S] = eig(H); S = diag(S);
clf; plot(real(S), '.-'); axis tight;



