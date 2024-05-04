function [h,xd] = plot_ds_model_3D_2Dslice_in_3D_plot(fig, ds, attractor, limits, dimensions, varargin)

quality='low';

if nargin > 4
    quality = varargin{1};
end

if strcmpi(quality,'high')
    nx=400;
    ny=400;
elseif strcmpi(quality,'medium')
    nx=200;
    ny=200;
else
    nx=50;
    ny=50;
end

axlim           = limits;
ax_x            = linspace(axlim(1),axlim(2),nx); % computing the mesh points along each axis
ax_y            = linspace(axlim(3),axlim(4),ny); % computing the mesh points along each axis
ax_z            = linspace(axlim(5),axlim(6),ny);
[d1_tmp, d2_tmp,d3_tmp]= meshgrid(ax_x,ax_y,ax_z);  % meshing the input domain

% d3_tmp = zeros(size(d1_tmp));
d1_vec = d1_tmp(:);
d2_vec = d2_tmp(:);
d3_vec = d3_tmp(:);

dimensions

x = repmat(attractor,[1 length(d1_vec)]);

x(1,:) = d1_vec';
x(2,:) = d2_vec';
x(3,:) = d3_vec';

xd = feval(ds, x);
fprintf('done\n');

if dimensions==[1 2]
h = streamslice(d1_tmp, d2_tmp,d3_tmp,reshape(xd(dimensions(1),:),ny,nx,nx), reshape(xd(dimensions(2),:),ny,nx,nx),reshape(xd(3,:),ny,nx,nx),[],[],0);
elseif dimensions==[1 3]
h = streamslice(d1_tmp, d2_tmp,d3_tmp,reshape(xd(dimensions(1),:),ny,nx,nx),reshape(xd(2,:),ny,nx,nx), reshape(xd(dimensions(2),:),ny,nx,nx),[],0,[]);
else
h = streamslice(d1_tmp, d2_tmp,d3_tmp,reshape(xd(1,:),ny,nx,nx),reshape(xd(dimensions(1),:),ny,nx,nx), reshape(xd(dimensions(2),:),ny,nx,nx),0,[],[]);
end


set(h,'LineWidth', 0.5)

end