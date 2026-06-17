clear all; close all; clc;
tic

%% Model selection : set exactly one kinetics and one growth law equal to 1.

Schnakenberg = 1;
Brusselator  = 0;

linear       = 1;
exponen      = 0;
logistic     = 0; flipped = 0; %saturation

%% Model and simulation parameters

showanimation = 1;
makegif       = 1;
drawperframe  = 10000;
save_data     = 0;

% Coupling and growth parameters
eta = 2;
gamma1 = 0.05;
gamma2 = 0.05;

% Logistic saturation setting
Linf1 = 0;
Linf2 = 0;

%% Output directory

if Schnakenberg
    kineticsName = 'Schnakenberg';
elseif Brusselator
    kineticsName = 'Brusselator';
end

if linear
    growthName = 'Linear';
elseif exponen
    growthName = 'Exponential';
elseif logistic
    growthName = 'Logistic';
end

folder = fullfile(pwd,'output',kineticsName,growthName);

if ~exist(folder,'dir')
    mkdir(folder);
end

%% Run information

gr1=round(gamma1,5); gr1=num2str(gr1);
gr2=round(gamma2,5); gr2=num2str(gr2);

etaStr=num2str(round(eta,4));

runID = datestr(datetime('now'),'yyyymmdd_HHMMSS');

prefix = sprintf('%s_%s_eta_%s_g1_%s_g2_%s',...
    runID,...
    growthName,...
    etaStr,...
    gr1,gr2);

runPrefix = fullfile(folder,prefix);

diary([runPrefix,'.txt']);

%% Backup current script

fprintf('saving to %s\n',folder);

origFile = [mfilename("fullpath") '.m'];

newbackup = fullfile(folder, sprintf('%s_script_backup.txt',prefix));

copyfile(origFile, newbackup);

fprintf('Backup created: %s\n', newbackup);

%% Kinetics parameters

if Schnakenberg
    % 
    a=0.05; b=1.4;
    D=[1 50 40 800];
    %
    f=@(u1,v1,u2,v2,eta) a - u1 + u1.^2.*v1 + eta*(u2-u1);
    g=@(u1,v1,u2,v2,eta) b - u1.^2.*v1 + eta*(v2-v1);
    %
    u1_0=a+b; v1_0=b/(a+b)^2;
    u2_0=a+b; v2_0=b/(a+b)^2;
elseif Brusselator
    %
    a=3; b=4.9;
    D = [1 10 15 115];
    %
    f = @(u1,v1,u2,v2,eta) a - (1+b)*u1 + u1.^2.*v1 + eta*(u2-u1);
    g = @(u1,v1,u2,v2,eta) b*u1 - u1.^2.*v1 + eta*(v2-v1);
    %
    u1_0=a; v1_0=b/a;
    u2_0=a; v2_0=b/a;
end
Du1=D(1); Dv1=D(2); Du2=D(3); Dv2=D(4); 

%% Growth functions

if linear
    L = @(t,gamma) 1+gamma*t;
elseif exponen
    L = @(t,gamma) exp(gamma*t);
elseif logistic
    L = @(t,gamma,Linf) Linf ./ (1 + (Linf-1).*exp(-gamma.*t));
end

if linear
    h = @(t,gamma) gamma/L(t,gamma);
elseif exponen
    h = @(t,gamma) gamma;
elseif logistic
    h = @(t,gamma,Linf) gamma .* (Linf - L(t,gamma,Linf)) ./ Linf;
end

%% Simulation summary

fprintf('\n=====================================\n');
fprintf(' Bilayer Growth RD Simulation\n');
fprintf('=====================================\n');
fprintf(' Kinetics : %s\n', kineticsName);
fprintf(' Growth   : %s\n', growthName);
fprintf(' eta      : %s\n', etaStr);
fprintf(' gamma1   : %s\n', gr1);
fprintf(' gamma2   : %s\n', gr2);
fprintf(' u1_0      : %.4f\n', u1_0);
fprintf(' v1_0      : %.4f\n', v1_0);
fprintf(' u2_0      : %.4f\n', u2_0);
fprintf(' v2_0      : %.4f\n', v2_0);
fprintf('=====================================\n\n');

%% Load Mesh Data: p-nodes; t-elements;
L1=[0,1];
L2=[0,1];
Ne = 400; 

% Node coordinates
p1 = linspace(L1(1), L1(2), Ne + 1)'; 
p2 = linspace(L2(1), L2(2), Ne + 1)';
% Element connectivity
t = [1:Ne; 2:Ne+1]'; 
t2=t; t1=t;

n2=length(p2);
n1=length(p1);

%% Assemble Matrices
ord=3;
[S_1,M_1]=AssembleGlobalMatrices1D(p1,t1,ord);
[S_2,M_2]=AssembleGlobalMatrices1D(p2,t2,ord);

%% time discretization
T=1000;
dt=0.01;
nt=T/dt+1;
stopti=nt;

%% Preallocate storage
max_steps = nt;
u1_store = zeros(max_steps,n1);
v1_store = zeros(max_steps,n1);
u2_store = zeros(max_steps,n2);
v2_store = zeros(max_steps,n2);
tt = zeros(max_steps,1); 

%% Initial condition
rng(0); %change random seed
Perturbations1 = 0.001*randn(n1,1);
Perturbations2 = 0.001*randn(n2,1);
u1=u1_0 + Perturbations1;
v1=v1_0 + Perturbations1;
u2=u2_0 + Perturbations2;
v2=v2_0 + Perturbations2;

%% set up figure
giffile = [runPrefix,'_pattern','.gif'];
fig = figure('Color','w');

subplot(2,2,1)
hold on
u1_fig=plot(p1, u1, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$u_s$','Interpreter','latex');
figtitle1=title('t=0');
set(gca,'FontSize',12);
hold off  

subplot(2,2,2)
hold on
u2_fig=plot(p2, u2, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$u_c$','Interpreter','latex');
figtitle2=title('t=0');
set(gca,'FontSize',12);
hold off  

subplot(2,2,3)
hold on
v1_fig=plot(p1, v1, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$v_s$','Interpreter','latex');
set(gca,'FontSize',12);
hold off  

subplot(2,2,4)
hold on
v2_fig=plot(p2, v2, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$v_c$','Interpreter','latex');
set(gca,'FontSize',12);
hold off  

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

    if logistic
        Lcur1  = L(t,gamma1,Linf1);
        sigma1 = h(t,gamma1,Linf1);
        Lcur2  = L(t,gamma2,Linf2);
        sigma2 = h(t,gamma2,Linf2);
    else
        Lcur1  = L(t,gamma1);
        sigma1 = h(t,gamma1);
        Lcur2  = L(t,gamma2);
        sigma2 = h(t,gamma2);
    end

    scale1 = (1 ./ Lcur1).^2;
    scale2 = (1 ./ Lcur2).^2;

    % physical coordinates for each layer (for plotting)
    x_phys1 = p1 * Lcur1;
    x_phys2 = p2 * Lcur2;

    % Right-hand side assembly
    ord=2;
    Fu1 = ReactKineInt1D(p1,t1,p2,t2,u1,v1,u2,v2,eta,f,ord);
    Gv1 = ReactKineInt1D(p1,t1,p2,t2,u1,v1,u2,v2,eta,g,ord);
    Fu2 = ReactKineInt1D(p2,t2,p1,t1,u2,v2,u1,v1,eta,f,ord);
    Gv2 = ReactKineInt1D(p2,t2,p1,t1,u2,v2,u1,v1,eta,g,ord);

    u1_new = (M_1 + dt/2*Du1*scale1*S_1)\ ...
             ((M_1 - dt/2*Du1*scale1*S_1)*u1 + dt*(Fu1 - sigma1*u1));

    v1_new = (M_1 + dt/2*Dv1*scale1*S_1)\ ...
             ((M_1 - dt/2*Dv1*scale1*S_1)*v1 + dt*(Gv1 - sigma1*v1));

    u2_new = (M_2 + dt/2*Du2*scale2*S_2)\ ...
             ((M_2 - dt/2*Du2*scale2*S_2)*u2 + dt*(Fu2 - sigma2*u2));

    v2_new = (M_2 + dt/2*Dv2*scale2*S_2)\ ...
             ((M_2 - dt/2*Dv2*scale2*S_2)*v2 + dt*(Gv2 - sigma2*v2));
    
    % Update and store
    u1 = u1_new; v1 = v1_new;
    u2 = u2_new; v2 = v2_new;

    u1_store(ti,:) = u1;
    v1_store(ti,:) = v1;
    u2_store(ti,:) = u2;
    v2_store(ti,:) = v2;

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        if showanimation
            u1_fig.XData = x_phys1;
            u1_fig.YData = u1;
            v1_fig.XData = x_phys1;
            v1_fig.YData = v1;
            u2_fig.XData = x_phys2;
            u2_fig.YData = u2;
            v2_fig.XData = x_phys2;
            v2_fig.YData = v2;
            figtitle1.String = ['t = ', num2str(t, '%.4f')];
            figtitle2.String = ['t = ', num2str(t, '%.4f')];
            drawnow;
        end
        if makegif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ti == 1 
                imwrite(imind, cm, giffile, 'gif', 'Loopcount', inf);
            else
                imwrite(imind, cm, giffile, 'gif', 'WriteMode', ...
                    'append', 'DelayTime', 0);
            end
        end
    end

    if ti == stopti
        break;
    end

end

%% saving final pattern
saveas(fig,[runPrefix,'_final.png']);
saveas(fig,[runPrefix,'_final.fig']);

%% Save raw solution data and metadata

if save_data
    save([runPrefix,'_u1_store.mat'], ...
        'u1_store', ...
        '-v7.3');

    save([runPrefix,'_v1_store.mat'], ...
        'v1_store', ...
        '-v7.3');

    save([runPrefix,'_u2_store.mat'], ...
        'u2_store', ...
        '-v7.3');

    save([runPrefix,'_v2_store.mat'], ...
        'v2_store', ...
        '-v7.3');

    save([runPrefix,'_time.mat'], ...
        'tt', ...
        '-v7.3');

    save([runPrefix,'_metadata.mat'], ...
        ...
        'p1','p2', ...
        'L1','L2', ...
        'Ne', ...
        ...
        'gamma1','gamma2', ...
        'Schnakenberg','Brusselator',...
        'linear','exponen','logistic',...
        ...
        'a','b','D','eta', ...
        ...
        'dt','T', ...
        ...
        '-v7.3');
end

%% Simulation figure settings

t_all = tt(1:ti);
if logistic
    L_all1 = L(t_all,gamma1,Linf1);
    L_all2 = L(t_all,gamma2,Linf2);
else    
    L_all1 = L(t_all,gamma1);
    L_all2 = L(t_all,gamma2);
end

% Fixed plotting grid

Nplot = Ne;

L1max = max(L_all1(:));
L2max = max(L_all2(:));
Xgrid1 = linspace(0,L1max,Nplot);
Xgrid2 = linspace(0,L2max,Nplot);

% Allocate interpolated arrays

U_rescaled_u1 = nan(ti,Nplot);
U_rescaled_v1 = nan(ti,Nplot);

U_rescaled_u2 = nan(ti,Nplot);
U_rescaled_v2 = nan(ti,Nplot);

% Interpolate evolving domains onto fixed grid

for k = 1:ti

    tcur = t_all(k);

    % Current physical domains
    if logistic
        Lcur1 = L(tcur,gamma1,Linf1);
        Lcur2 = L(tcur,gamma2,Linf2);
    else
        Lcur1 = L(tcur,gamma1);
        Lcur2 = L(tcur,gamma2);
    end

    x_phys1 = p1 * Lcur1;
    x_phys2 = p2 * Lcur2;

    % Current solution slices
    u1row = u1_store(k,:);
    v1row = v1_store(k,:);

    u2row = u2_store(k,:);
    v2row = v2_store(k,:);

    % Linear interpolation
    U_rescaled_u1(k,:) = interp1( ...
        x_phys1,...
        u1row,...
        Xgrid1,...
        'linear',...
        NaN);

    U_rescaled_v1(k,:) = interp1( ...
        x_phys1,...
        v1row,...
        Xgrid1,...
        'linear',...
        NaN);

    U_rescaled_u2(k,:) = interp1( ...
        x_phys2,...
        u2row,...
        Xgrid2,...
        'linear',...
        NaN);

    U_rescaled_v2(k,:) = interp1( ...
        x_phys2,...
        v2row,...
        Xgrid2,...
        'linear',...
        NaN);

end

%% Figure font settings

fs_axis  = 16;
fs_ticks = 16;
fs_cb    = 12;

%% Surface layer u_1

figure('Color','w')

imagesc(Xgrid1,t_all,U_rescaled_u1)

axis xy

xt = xticks;
yt = yticks;
xticklabels(string(xt/100))
yticklabels(string(yt/1000))

xlabel('$x_1\,(10^2)$','Interpreter','latex','FontSize',fs_axis);
ylabel('$t\,(10^3)$','Interpreter','latex','FontSize',fs_axis);

% title('$u_1$','Interpreter','latex')

xlim([0,max(L_all1)])

set(gca,...
    'FontSize',fs_ticks,...
    'LineWidth',1.2);

set(gca,'Color','w')

set(findobj(gca,'Type','image'),...
    'AlphaData',~isnan(U_rescaled_u1))

validU1 = U_rescaled_u1(~isnan(U_rescaled_u1));

clim([ ...
    prctile(validU1,2), ...
    prctile(validU1,98)])

colormap(parula)

cb = colorbar;
cb.FontSize = fs_cb;

hold on

plot(L1(1)+L_all1,...
     t_all,...
     'r-',...
     'LineWidth',2)

plot(L1(1)+zeros(size(t_all)),...
     t_all,...
     'r-',...
     'LineWidth',1)

hold off

drawnow

% Export
exportgraphics(gcf,...
    [runPrefix,'_growthu1.png'],...
    'Resolution',600);

exportgraphics(gcf,...
    [runPrefix,'_growthu1.pdf'],...
    'ContentType','image');

%% Surface layer v_1

figure('Color','w')

imagesc(Xgrid1,t_all,U_rescaled_v1)

axis xy

xt = xticks;
yt = yticks;
xticklabels(string(xt/100))
yticklabels(string(yt/1000))

xlabel('$x_1\,(10^2)$','Interpreter','latex','FontSize',fs_axis);
ylabel('$t\,(10^3)$','Interpreter','latex','FontSize',fs_axis);

% title('$v_1$','Interpreter','latex')

xlim([0,max(L_all1)])

set(gca,...
    'FontSize',fs_ticks,...
    'LineWidth',1.2);

set(gca,'Color','w')

set(findobj(gca,'Type','image'),...
    'AlphaData',~isnan(U_rescaled_v1))

validV1 = U_rescaled_v1(~isnan(U_rescaled_v1));

clim([ ...
    prctile(validV1,2), ...
    prctile(validV1,98)])

colormap(parula)

cb = colorbar;
cb.FontSize = fs_cb;

hold on

plot(L1(1)+L_all1,...
     t_all,...
     'r-',...
     'LineWidth',2)

plot(L1(1)+zeros(size(t_all)),...
     t_all,...
     'r-',...
     'LineWidth',1)

hold off

drawnow

% Export
exportgraphics(gcf,...
    [runPrefix,'_growthv1.png'],...
    'Resolution',600);

exportgraphics(gcf,...
    [runPrefix,'_growthv1.pdf'],...
    'ContentType','image');

%% Core layer u_2

figure('Color','w')

imagesc(Xgrid2,t_all,U_rescaled_u2)

axis xy

xt = xticks;
yt = yticks;
xticklabels(string(xt/100))
yticklabels(string(yt/1000))

xlabel('$x_2\,(10^2)$','Interpreter','latex','FontSize',fs_axis);
ylabel('$t\,(10^3)$','Interpreter','latex','FontSize',fs_axis);

% title('$u_2$','Interpreter','latex')

xlim([0,max(L_all2)])

set(gca,...
    'FontSize',fs_ticks,...
    'LineWidth',1.2);

set(gca,'Color','w')

set(findobj(gca,'Type','image'),...
    'AlphaData',~isnan(U_rescaled_u2))

validU2 = U_rescaled_u2(~isnan(U_rescaled_u2));

clim([ ...
    prctile(validU2,2), ...
    prctile(validU2,98)])

colormap(parula)

cb = colorbar;
cb.FontSize = fs_cb;

hold on

plot(L2(1)+L_all2,...
     t_all,...
     'r-',...
     'LineWidth',2)

plot(L2(1)+zeros(size(t_all)),...
     t_all,...
     'r-',...
     'LineWidth',1)

hold off

drawnow

% Export
exportgraphics(gcf,...
    [runPrefix,'_growthu2.png'],...
    'Resolution',600);

exportgraphics(gcf,...
    [runPrefix,'_growthu2.pdf'],...
    'ContentType','image');

%% Core layer v_2

figure('Color','w')

imagesc(Xgrid2,t_all,U_rescaled_v2)

axis xy

xt = xticks;
yt = yticks;
xticklabels(string(xt/100))
yticklabels(string(yt/1000))

xlabel('$x_2\,(10^2)$','Interpreter','latex','FontSize',fs_axis);
ylabel('$t\,(10^3)$','Interpreter','latex','FontSize',fs_axis);

% title('$v_2$','Interpreter','latex')

xlim([0,max(L_all2)])

set(gca,...
    'FontSize',fs_ticks,...
    'LineWidth',1.2);

set(gca,'Color','w')

set(findobj(gca,'Type','image'),...
    'AlphaData',~isnan(U_rescaled_v2))

validV2 = U_rescaled_v2(~isnan(U_rescaled_v2));

clim([ ...
    prctile(validV2,2), ...
    prctile(validV2,98)])

colormap(parula)

cb = colorbar;
cb.FontSize = fs_cb;

hold on

plot(L2(1)+L_all2,...
     t_all,...
     'r-',...
     'LineWidth',2)

plot(L2(1)+zeros(size(t_all)),...
     t_all,...
     'r-',...
     'LineWidth',1)

hold off

drawnow

% Export
exportgraphics(gcf,...
    [runPrefix,'_growthv2.png'],...
    'Resolution',600);

exportgraphics(gcf,...
    [runPrefix,'_growthv2.pdf'],...
    'ContentType','image');

%% Optional: plot growth histories for verification
figure('Color','w');
plot(L1(1) + L_all1, t_all, 'b-', 'LineWidth', 2); hold on;
plot(L2(1) + L_all2, t_all, 'r--', 'LineWidth', 2);
ylabel('Time t'); xlabel('Right boundary x');
legend('Surface right boundary','Core right boundary');
title('Layer boundaries vs time (differential growth)');
set(gca,'FontSize',12);

saveas(gcf,[runPrefix,'_rates.png']);

%%
fprintf('\nDone!\n');
toc
diary OFF
