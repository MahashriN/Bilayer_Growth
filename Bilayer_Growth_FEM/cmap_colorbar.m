function cmap_colorbar()
lvls = [ -0.1, 0, 0.1];
% load('blueyellow.mat')
load('brownyellow.mat')
% find the indices where the level difference is "big" (0.1):
lvl_idx = find(diff(lvls) > 0.001); % (closer to 0.1 than to 0.05)
n_big = numel(lvl_idx);
% create a colormap with an extra color for each big level difference:
% cmap = BY; 
cmap = twocolors;
% binary code map color
% change indices in lvls to indices in the colormap:
cmap_idx = lvl_idx;
for ii = 1:n_big
    cmap_idx(ii) = lvl_idx(ii) + nnz(lvl_idx <= lvl_idx(ii));
end
% duplicate the colors in the colormap at the big level indices:
cmap(cmap_idx,:) = cmap(cmap_idx-1,:);
% apply the colormap, create the colorbar, and set the ticks to lvls:
colormap(cmap);
cb = colorbar();
cb.Ticks = lvls;
colorbar('off')
end
