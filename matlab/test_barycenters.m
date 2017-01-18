%%
% test for sinkhorn barycenters.

addpath('toolbox/');

rep = 'results/barycenters/';
[~,~] = mkdir(rep);

% number of points
N = 150; 

%%
% Helpers.

SetAR = @(ar)set(gca, 'PlotBoxAspectRatio', [1 ar 1]);
myplot = @(x,y,c,a)area(x, y, 'FaceColor', c, 'EdgeColor', c, 'FaceAlpha', a, 'EdgeAlpha', a);
x = (0:N-1)'/N;
normalize = @(x)x/sum(x(:));
gauss = @(m,s)exp(-(x-m).^2/(2*s^2));

%%
% Load two Gaussians. 


sigma = .03;
mu = {...
    1.7*gauss(.1,sigma) + 0.6*gauss(.7,sigma), ...
    gauss(.25,sigma)  + gauss(.9,sigma), ...
   };
vmin = 1e-3;
for k=1:2
    mu{k} =  normalize( vmin + mu{k} );
end

%%
% Display inputs.

clf; hold on;
myplot(x,mu{1}, 'r',1);
myplot(x,mu{2}, 'b',1);
axis tight; box on;

%% 
% Cost matrix.

[Y,X] = meshgrid(x,x);
c = (X-Y).^2;

%% 
% Parameters.

epsilon = (.02).^2;
options.null = 0;
% Change this to switch between balanced and unbalanced
if not(isfield(options, 'rho'))
    options.rho = 0.3; % unbalanced
    options.rho = Inf; % balanced
end
options.niter = 4000;
 
%%
% compute barycenters.

K = 10; % number of computed barycenters
nu = {};
for k=1:K
    t = (k-1)/(K-1);
    w = [1-t t];
    options.tau = 0;
    options.tau_v = 0;
    options.disp_rate = 50;
    [nu{k},gamma] = barycenter_log(mu,c,epsilon,w,options);
end

%% 
% Rendering

s = 0; F = [];
for k=[1:K]
    s = s+1;
    t = (k-1)/(K-1);
    clf;
    hold on;
    h = myplot(x, nu{1}, 'b', .25); alpha(h,.25);
    h = myplot(x, nu{end}, 'r', .25); alpha(h,.25);
    h = myplot(x, nu{k}, [t 0 1-t], .8); alpha(h,.8);
    axis([0 1 0 1.05*max([nu{1}(:);nu{end}(:)])]);
    set(gca, 'XTick', [], 'YTick', []); 
    box on; SetAR(1/4);
    drawnow;    
    if k==1
     %   set(gca,'nextplot','replacechildren','visible','off')
    end
    f = getframe; 
    F(:,:,:,s) = f.cdata;
end

%% 
% Save as .gif file.

% find colormap         
A = permute(F, [1 2 4 3]);
A = reshape(A, [size(A,1) size(A,2)*size(A,3) 3]);
[~,map] = rgb2ind(uint8(A),254,'nodither');
map(end+1,:) = 0; map(end+1,:) = 1; 
% convert
im = [];
for s=1:size(F,4);
    im(:,:,1,s) = rgb2ind(uint8(F(:,:,:,s)),map,'nodither');
end
% save
imwrite(im+1,map,[rep 'interp-rho' num2str(options.rho) '.gif'], ...
        'DelayTime',0,'LoopCount',inf);
