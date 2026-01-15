clear all; close all; clc;
tic

linear=01;
exponen=0;
logistic=0;

Linf=1000; % logistic-saturation

%% pre-setting
showanimation=1;
makegif=1;
drawperframe=500;
tol=1e-9;

% coupling strenth
eta = 2;

% Different growth rates for surface and core layers
gammaS = 0.05;
gammaC = 0.1;

if linear
T=15000;
elseif exponen
    if gammaS<gammaC
        gamma=gammaS;
    else
        gamma=gammaC;
    end
T=round(log(1500)/gamma,0);
elseif logistic
T = 2500;
end

%% saving folder and name
folder='D:\MATLAB Output\Bilayer_Growth_FEM\';
time = datestr(datetime('now'),'yyyymmdd_HHMMSS');
prefix = [folder,time];
diary([prefix,'.txt']);
fprintf('saving to %s\n',folder);

%%
FileNameAndLocation = mfilename('fullpath');
[filepath,name,~] = fileparts(FileNameAndLocation);
ext='.m';
origFile = fullfile(filepath,[name ext]);
newbackup = fullfile(folder, sprintf('%s_%s.txt',time,name));
copyfile(origFile, newbackup);
fprintf('Backup created: %s\n', newbackup);

%% Load Mesh Data: p-nodes; t-elements;
LS=[0,1];
LC=[0,1];
Ne = 400; 

% Node coordinates
pS = linspace(LS(1), LS(2), Ne + 1)'; 
pC = linspace(LC(1), LC(2), Ne + 1)';
% Element connectivity
t = [1:Ne; 2:Ne+1]'; 
tC=t; tS=t;

nC=length(pC);
nS=length(pS);

%% Parameters
a=0.05; b=1.4;
D=[1 50 40 800];
DuS=D(1); DvS=D(2); DuC=D(3); DvC=D(4); 

%% Define the Source Term + Coupling Term
f=@(u1,v1,u2,v2,eta) a - u1 + u1.^2.*v1 + eta*(u2-u1);
g=@(u1,v1,u2,v2,eta) b - u1.^2.*v1 + eta*(v2-v1);

%% Initial Steady states
uS0=a+b; vS0=b/(a+b)^2;
uC0=a+b; vC0=b/(a+b)^2;

fprintf(['\nInitial steady state for\n' ...
    'surface layer is: uS0=%.4f, vS0=%.4f\n'...
    'core layer is: uC0=%.4f, vC0=%.4f\n'], ...
    uS0,vS0,uC0,vC0);

%% Growthrates
if linear
L=@(t,gamma) 1+gamma*t;
elseif exponen
L=@(t,gamma) exp(gamma*t);
elseif logistic
L = @(t,gamma) Linf ./ (1 + (Linf-1).*exp(-gamma.*t));
end

if linear
h=@(t,gamma) gamma/L(t,gamma);
elseif exponen
h=@(t,gamma) gamma;
elseif logistic
h = @(t,gamma) gamma .* (Linf - L(t,gamma)) ./ Linf;
end

%% Assemble Matrices
ord=3;
[S_S,M_S]=AssembleGlobalMatrices1D(pS,tS,ord);
[S_C,M_C]=AssembleGlobalMatrices1D(pC,tC,ord);

%% time discretization
% T=1000;
dt=0.01;
nt=T/dt+1;
stopti=nt;

%% Preallocate storage
max_steps = nt;
uS_store = zeros(max_steps,nS);
vS_store = zeros(max_steps,nS);
uC_store = zeros(max_steps,nC);
vC_store = zeros(max_steps,nC);
tt = zeros(max_steps,1); 

%%
grS=round(gammaS,4);
grC=round(gammaC,4);
fprintf('\nFor growth_rateS=%.4f, growth_rateC=%.4f,\n',grS,grC);
gr=['_',num2str(grS),'_',num2str(grC)];

for e=1:length(eta)
eval=round(eta(e),4);
fprintf('\nFor eta=%.4f,\n',eval);
eval=['_',num2str(eval)];

pattern_detected = false;
pattern_end=30.05;
stopti=nt;

%% Initial condition
rng(0); %change random seed
Perturbations2 = 0.01*randn(nS,1);
Perturbations1 = 0.01*randn(nC,1);
uS=uS0 + Perturbations2;
vS=vS0 + Perturbations2;
uC=uC0 + Perturbations1;
vC=vC0 + Perturbations1;

%% set up figure
giffile = [prefix,eval,gr,'_pattern','.gif'];
fig = figure('Color','w');%,'WindowState', 'maximized');

subplot(2,2,1)
hold on
uS_fig=plot(pS, uS, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$u_s$','Interpreter','latex');
figtitleS=title('t=0');
set(gca,'FontSize',12);
hold off  

subplot(2,2,2)
hold on
uC_fig=plot(pC, uC, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$u_c$','Interpreter','latex');
figtitleC=title('t=0');
set(gca,'FontSize',12);
hold off  

subplot(2,2,3)
hold on
vS_fig=plot(pS, vS, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$v_s$','Interpreter','latex');
set(gca,'FontSize',12);
hold off  

subplot(2,2,4)
hold on
vC_fig=plot(pC, vC, 'b.-', 'LineWidth', 2, 'MarkerSize', 15);
xlabel('$x$','Interpreter','latex');
ylabel('$v_c$','Interpreter','latex');
set(gca,'FontSize',12);
hold off  

%% Simulation iteration
for ti = 1:nt
    t = dt * (ti - 1);
    tt(ti) = t;

    LcurS = L(t,gammaS);
    sigmaS = h(t,gammaS);
    scaleS = (1 ./ LcurS).^2;

    LcurC = L(t,gammaC);
    sigmaC = h(t,gammaC);
    scaleC = (1 ./ LcurC).^2;

    % physical coordinates for each layer (for plotting)
    x_physS = pS * LcurS;
    x_physC = pC * LcurC;

    % --- Right-hand side assembly ---
    ord=2;
    FuS = ReactKineInt1D(pS,tS,pC,tC,uS,vS,uC,vC,eta(e),f,ord);
    GvS = ReactKineInt1D(pS,tS,pC,tC,uS,vS,uC,vC,eta(e),g,ord);
    FuC = ReactKineInt1D(pC,tC,pS,tS,uC,vC,uS,vS,eta(e),f,ord);
    GvC = ReactKineInt1D(pC,tC,pS,tS,uC,vC,uS,vS,eta(e),g,ord);

    uS_new = (M_S + dt/2*DuS*scaleS*S_S)\ ...
             ((M_S - dt/2*DuS*scaleS*S_S)*uS + dt*(FuS - sigmaS*uS));

    vS_new = (M_S + dt/2*DvS*scaleS*S_S)\ ...
             ((M_S - dt/2*DvS*scaleS*S_S)*vS + dt*(GvS - sigmaS*vS));

    uC_new = (M_C + dt/2*DuC*scaleC*S_C)\ ...
             ((M_C - dt/2*DuC*scaleC*S_C)*uC + dt*(FuC - sigmaC*uC));

    vC_new = (M_C + dt/2*DvC*scaleC*S_C)\ ...
             ((M_C - dt/2*DvC*scaleC*S_C)*vC + dt*(GvC - sigmaC*vC));
    
    % Pattern detection using norms
    erruS = norm(uS_new - uS)/norm(uS);
    errvS = norm(vS_new - vS)/norm(vS);
    erruC = norm(uC_new - uC)/norm(uC);
    errvC = norm(vC_new - vC)/norm(vC);
    maxErr = max([erruS, errvS, erruC, errvC]);

    % Update and store
    uS = uS_new; vS = vS_new;
    uC = uC_new; vC = vC_new;

    uS_store(ti,:) = uS;
    vS_store(ti,:) = vS;
    uC_store(ti,:) = uC;
    vC_store(ti,:) = vC;

    %% Plot and gif
    if mod(ti, drawperframe) == 1
        if showanimation
            uS_fig.XData = x_physS;
            uS_fig.YData = uS;
            vS_fig.XData = x_physS;
            vS_fig.YData = vS;
            uC_fig.XData = x_physC;
            uC_fig.YData = uC;
            vC_fig.XData = x_physC;
            vC_fig.YData = vC;
            figtitleS.String = ['t = ', num2str(t, '%.4f')];
            figtitleC.String = ['t = ', num2str(t, '%.4f')];
            drawnow;
        end
        if makegif
            frame = getframe(fig);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            if ti == 1
                imwrite(imind, cm, giffile, 'gif', 'Loopcount', inf);
            else
                imwrite(imind, cm, giffile, 'gif', 'WriteMode', 'append', 'DelayTime', 0);
            end
        end
    end

    if ~pattern_detected && maxErr < tol
        fprintf('\nPattern fully developed at t = %.5f\n', t);
        pattern_detected = true;
        stopti = ti + round(pattern_end / dt);
    end

    if ti == stopti
        break;
    end

end
%% saving final pattern
saveas(fig,[prefix,eval,gr,'_final.png']);
saveas(fig,[prefix,eval,gr,'_final.fig']);

%% Space-time graphs
t_all = tt(1:ti);
L_allS = L(t_all,gammaS);
L_allC = L(t_all,gammaC);
Lmax = max([L_allS(:); L_allC(:)]); % overall max length required

xmin_plot = min(LS(1), LC(1));
Xgrid = linspace(xmin_plot, xmin_plot + Lmax, Ne); 

U_rescaled_uS = nan(ti, numel(Xgrid));
U_rescaled_vS = nan(ti, numel(Xgrid));
U_rescaled_uC = nan(ti, numel(Xgrid));
U_rescaled_vC = nan(ti, numel(Xgrid));

for k = 1:ti

    tcur = t_all(k);
    LcurS = L(tcur,gammaS);
    LcurC = L(tcur,gammaC);
    x_physS = pS * LcurS;
    x_physC = pC * LcurC;

    uSrow = uS_store(k, :);
    vSrow = vS_store(k, :);
    uCrow = uC_store(k, :);
    vCrow = vC_store(k, :);

    U_rescaled_uS(k, :) = interp1(x_physS, uSrow, Xgrid, 'linear', NaN);
    U_rescaled_vS(k, :) = interp1(x_physS, vSrow, Xgrid, 'linear', NaN);
    U_rescaled_uC(k, :) = interp1(x_physC, uCrow, Xgrid, 'linear', NaN);
    U_rescaled_vC(k, :) = interp1(x_physC, vCrow, Xgrid, 'linear', NaN);
end

%% Create the space-time figure

% ---------- Surface layer: u_S ----------
figure('Color','w')
imagesc(Xgrid, t_all, U_rescaled_uS);
axis xy;
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
title('$u_1$','Interpreter','latex');
set(gca,'Color','w')         
set(findobj(gca,'Type','image'),'AlphaData',~isnan(U_rescaled_uS))
validUS = U_rescaled_uS(~isnan(U_rescaled_uS));
clim([prctile(validUS, 2), prctile(validUS, 98)])
colormap(parula); colorbar;
hold on;
plot(LS(1) + L_allS, t_all, 'r-', 'LineWidth', 2); 
plot(LS(1) + zeros(size(t_all)), t_all, 'r-', 'LineWidth', 1);
hold off;
set(gca,'FontSize',14);
saveas(gcf,[prefix,eval,gr,'_growthu1.png']);
saveas(gcf,[prefix,eval,gr,'_growthu1.eps']);
if exponen || logistic
saveas(gcf,[prefix,eval,gr,'_growthu1.fig']);
end

%%
save([prefix,eval,gr,'_growthu1_data.mat'], ...
     'U_rescaled_uS', ...
     'Xgrid', ...
     't_all', ...
     'L_allS', ...
     'LS', ...
     '-v7.3');

%%
% % ---------- Surface layer: v_S ----------
figure('Color','w')
imagesc(Xgrid, t_all, U_rescaled_vS);
axis xy;
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
title('$v_1$','Interpreter','latex');
set(gca,'Color','w')        
set(findobj(gca,'Type','image'),'AlphaData',~isnan(U_rescaled_vS))
validVS = U_rescaled_vS(~isnan(U_rescaled_vS));
clim([prctile(validVS,2), prctile(validVS,98)])
colormap(parula); colorbar;
hold on;
plot(LS(1) + L_allS, t_all, 'r-', 'LineWidth', 2);
plot(LS(1) + zeros(size(t_all)), t_all, 'r-', 'LineWidth', 1);
hold off;
set(gca,'FontSize',12);
saveas(gcf,[prefix,eval,gr,'_growthv1.png']);
saveas(gcf,[prefix,eval,gr,'_growthv1.eps']);
if exponen || logistic
saveas(gcf,[prefix,eval,gr,'_growthv1.fig']);
end

%%
save([prefix,eval,gr,'_growthv1_data.mat'], ...
     'U_rescaled_vS', ...
     '-v7.3');

%%
% ---------- Core layer: u_C ----------
figure('Color','w')
imagesc(Xgrid, t_all, U_rescaled_uC);
axis xy;
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
title('$u_2$','Interpreter','latex');
set(gca,'Color','w')         
set(findobj(gca,'Type','image'),'AlphaData',~isnan(U_rescaled_uC))
validUC = U_rescaled_uC(~isnan(U_rescaled_uC));
clim([prctile(validUC,2), prctile(validUC,98)])
colormap(parula); colorbar;
hold on;
plot(LC(1) + L_allC, t_all, 'r-', 'LineWidth', 2); 
plot(LC(1) + zeros(size(t_all)), t_all, 'r-', 'LineWidth', 1);
hold off;
set(gca,'FontSize',14);
saveas(gcf,[prefix,eval,gr,'_growthu2.png']);
saveas(gcf,[prefix,eval,gr,'_growthu2.eps']);
if exponen || logistic
saveas(gcf,[prefix,eval,gr,'_growthu2.fig']);
end

%%
save([prefix,eval,gr,'_growthu2_data.mat'], ...
     'U_rescaled_uC', ...
     'Xgrid', ...
     't_all', ...
     'L_allC', ...
     'LC', ...
     '-v7.3');

%%
% % ---------- Core layer: v_C ----------
figure('Color','w')
imagesc(Xgrid, t_all, U_rescaled_vC);
axis xy;
xlabel('$x$','Interpreter','latex');
ylabel('$t$','Interpreter','latex');
title('$v_2$','Interpreter','latex');
set(gca,'Color','w')        
set(findobj(gca,'Type','image'),'AlphaData',~isnan(U_rescaled_vC))
validVC = U_rescaled_vC(~isnan(U_rescaled_vC));
clim([prctile(validVC,2), prctile(validVC,98)])
colormap(parula); colorbar;
hold on;
plot(LC(1) + L_allC, t_all, 'r-', 'LineWidth', 2);
plot(LC(1) + zeros(size(t_all)), t_all, 'r-', 'LineWidth', 1);
hold off;
set(gca,'FontSize',12);
saveas(gcf,[prefix,eval,gr,'_growthv2.png']);
saveas(gcf,[prefix,eval,gr,'_growthv2.eps']);
if exponen || logistic
saveas(gcf,[prefix,eval,gr,'_growthv2.fig']);
end

%%
save([prefix,eval,gr,'_growthv2_data.mat'], ...
     'U_rescaled_vC', ...
     '-v7.3');

%% Optional: plot growth histories for verification
figure('Color','w');
plot(t_all, LS(1) + L_allS, 'b-', 'LineWidth', 2); hold on;
plot(t_all, LC(1) + L_allC, 'r--', 'LineWidth', 2);
xlabel('Time t'); ylabel('Right boundary x');
legend('Surface right boundary','Core right boundary');
title('Layer boundaries vs time (differential growth)');
set(gca,'FontSize',12);

saveas(gcf,[prefix,eval,gr,'_rates.png']);

end

%%
fprintf('\nDone!\n');
toc
diary OFF
