clear all; close all; clc;
tic

linear=01;
exponen=0;
logistic=0; 

Linf=1000;

%% pre-setting

etaa = 2;

gammaS = 0.2;
gammaC = 0.1;

t_snaps = linspace(0, T, 6); %[0 1000:1000:5000 10000 15000]; %
Kmax    = 500; 

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
folder='D:MATLAB Output\Bilayer_Growth_FEM\';
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

%% Parameters
a=0.05; b=1.4;
D=[1 50 40 800];

dt=0.001;
nt=T/dt+1;

%% Define the Source Term + Coupling Term
f=@(u1,v1,u2,v2,eta) a - u1 + u1.^2.*v1 + eta*(u2-u1);
g=@(u1,v1,u2,v2,eta) b - u1.^2.*v1 + eta*(v2-v1);

%% Steady states at t=0
uS0=a+b; vS0=b/(a+b)^2;
uC0=a+b; vC0=b/(a+b)^2;

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

%% TDHSS
for e=1:length(etaa)
    for gj=1:length(gammaC)
        gi=gj;
grS=round(gammaS(gi),4);
grC=round(gammaC(gj),4);
fprintf('\nFor growth_rateS=%.4f, growth_rateC=%.4f,\n',grS,grC);
gr=['_',num2str(grS),'_',num2str(grC)];

eval=round(etaa(e),4);
fprintf('\nFor eta=%.4f,\n',eval);
eval=['_',num2str(eval)];

%% Integrate spatially homogeneous base ODE
tspan = [0, T];
y0 = [uS0; vS0; uC0; vC0];
baseODE = @(t,y) [
    f(y(1),y(2),y(3),y(4),etaa(e)) - h(t,gammaS(gi))*y(1);
    g(y(1),y(2),y(3),y(4),etaa(e)) - h(t,gammaS(gi))*y(2);
    f(y(3),y(4),y(1),y(2),etaa(e)) - h(t,gammaC(gj))*y(3);
    g(y(3),y(4),y(1),y(2),etaa(e)) - h(t,gammaC(gj))*y(4);
];
[t_base, y_base] = ode45(baseODE, tspan, y0);

fprintf('\nHomogeneous base ODE integrated for 0 <= t <= %.2f.\n',T);

uS_t = y_base(:,1);  vS_t = y_base(:,2);
uC_t = y_base(:,3);  vC_t = y_base(:,4);

%%
figure('Units','inches','Position',[1 1 6.5 5]); 

% ----------- U-PANEL -----------
subplot(2,1,1)
plot(t_base, uS_t, 'b',  'LineWidth',1.6); hold on
plot(t_base, uC_t, 'm--','LineWidth',1.6);

xlim([0 500]);
ylabel('Concentration','FontSize',12);
title('Time-dependent homogeneous base states','FontSize',14);
set(gca,'FontSize',11,'LineWidth',1);

lgd = legend('$u_S$','$u_C$');
set(lgd,'Interpreter','latex','FontSize',12,'Location','northeast');

% inset
axes('Position',[0.50 0.70 0.25 0.20]); % inset properly placed
box on
plot(t_base, uS_t,'b','LineWidth',1); hold on
plot(t_base, uC_t,'m--','LineWidth',1);
xlim([100 150]);
set(gca,'FontSize',9,'LineWidth',1);

% ----------- V-PANEL -----------
subplot(2,1,2)
plot(t_base, vS_t, 'r',  'LineWidth',1.6); hold on
plot(t_base, vC_t, 'k--','LineWidth',1.6);

xlim([0 500]);
xlabel('Time','FontSize',12);
ylabel('Concentration','FontSize',12);
set(gca,'FontSize',11,'LineWidth',1);

lgd = legend('$v_S$','$v_C$');
set(lgd,'Interpreter','latex','FontSize',12,'Location','northeast');

% inset
axes('Position',[0.50 0.22 0.25 0.20]);
box on
plot(t_base, vS_t,'r','LineWidth',1); hold on
plot(t_base, vC_t,'k--','LineWidth',1);
xlim([100 150]);
set(gca,'FontSize',9,'LineWidth',1);

saveas(gcf,[prefix,eval,gr,'_TDHSS.png']);
% saveas(gcf,[prefix,eval,gr,'_TDHSS.fig']);
    end
end

%% ==== Frozen-time dispersion at selected snapshots ====
% Kmax     = 60;
eta_grid = linspace(0,10,200);      
k_grid   = 0:Kmax;
[kmesh, emesh] = meshgrid(k_grid, eta_grid);

% choose snapshots
% t_snaps = linspace(0, T, 6);

for ts = 1:numel(t_snaps)
    t = t_snaps(ts);
    yb = interp1(t_base, y_base, t);
    uS = yb(1); vS = yb(2); uC = yb(3); vC = yb(4);

    % Jacobians of Schnakenberg (no dependence on eta)
    fu = @(u,v) -1 + 2*u.*v;  fv = @(u,v) u.^2;
    gu = @(u,v) -2*u.*v;      gv = @(u,v) -u.^2;

    J_S = [fu(uS,vS)  fv(uS);  gu(uS,vS)  gv(uS)];
    J_C = [fu(uC,vC)  fv(uC);  gu(uC,vC)  gv(uC)];

    % growth + diffusion scalings
    L1 = L(t, gammaS);   L2 = L(t, gammaC);
    h1 = h(t, gammaS);   h2 = h(t, gammaC);

    DM = diag([D(1)/L1^2, D(2)/L1^2, D(3)/L2^2, D(4)/L2^2]);
    HM = diag([h1, h1, h2, h2]);
    JM = [J_S, zeros(2); zeros(2), J_C];

    Mr = zeros(size(kmesh));  % max Re(lambda)

    for ii = 1:numel(kmesh)
        k  = kmesh(ii);
        et = emesh(ii);
        rho = (k*pi)^2;               % Neumann mode index eigenvalue
        CM  = et*[-eye(2) eye(2); eye(2) -eye(2)];

        Mk  = -rho*DM - HM + JM + CM;
        ev  = eig(Mk);
        Mr(ii) = max(real(ev));
    end

    figure('Color','w');
    surf(k_grid, eta_grid, Mr, 'linestyle','none'); 
    view(2); 
    shading interp;
    xlabel('mode index k'); ylabel('\eta'); 
    title(sprintf('max Re(\\lambda), t=%.1f', t));
    clim([-0.1 0.1]);
    cmap_colorbar()
    % colormap('jet')
    colorbar; 
    set(gca,'FontSize',14);
    saveas(gcf, sprintf('%s%s_disp_t%.1f.png', prefix, gr, t));
    saveas(gcf, sprintf('%s%s_disp_t%.1f.fig', prefix, gr, t));
end

%% ==== Frozen-time dispersion for a single eta value ====

% etaa = 2;          % <---- set desired eta value here
% Kmax = 60;               % highest mode index
k_grid = 0:Kmax;
% t_snaps = linspace(0, T, 6);

fprintf('\nComputing frozen-time dispersion: eta = %.2f\n', etaa);
figure('Color','w');
for ts = 1:numel(t_snaps)

    t = t_snaps(ts);
    yb = interp1(t_base, y_base, t);
    uS = yb(1); vS = yb(2); uC = yb(3); vC = yb(4);

    % Jacobians for Schnakenberg
    fu = @(u,v) -1 + 2*u.*v;  fv = @(u,v) u.^2;
    gu = @(u,v) -2*u.*v;      gv = @(u,v) -u.^2;

    J_S = [fu(uS,vS)  fv(uS);  gu(uS,vS)  gv(uS)];
    J_C = [fu(uC,vC)  fv(uC);  gu(uC,vC)  gv(uC)];

    % growth + diffusion scaling
    L1 = L(t, gammaS);
    L2 = L(t, gammaC);
    h1 = h(t, gammaS);
    h2 = h(t, gammaC);

    DM = diag([D(1)/L1^2, D(2)/L1^2, D(3)/L2^2, D(4)/L2^2]);
    HM = diag([h1, h1, h2, h2]);
    JM = [J_S, zeros(2); zeros(2), J_C];

    Mr = zeros(size(k_grid)); % max Re(lambda(k))

    for i = 1:length(k_grid)
        k = k_grid(i);
        rho = (k*pi)^2;  % Neumann eigenvalue (no domain length scaling)
        CM = etaa * [-eye(2) eye(2); eye(2) -eye(2)];

        Mk = -rho*DM - HM + JM + CM;
        ev = eig(Mk);
        Mr(i) = max(real(ev));
    end

    % ---- Plot ----
    plot(k_grid, Mr,'LineWidth', 1.8,'DisplayName', sprintf('t = %.1f', t)); 
    hold on;
    yline(0,'k--','linewidth',2,'HandleVisibility', 'off');  % stability threshold
    
    xlabel('mode index k','FontSize',12);
    ylabel('max Re(\lambda)','FontSize',12);
    title(sprintf('Frozen-time dispersion: t = %.1f, \\eta = %.1f', t, etaa), 'FontSize',14);
    grid on; set(gca,'FontSize',12);
    xlim([0,Kmax]);
    % ylim([min(Mr)-0.5,1])
end
legend('show', 'Location', 'northeast');
set(legend,'Interpreter','latex','FontSize',11);

% ---- Save figure ----
saveas(gcf, sprintf('%s%s%s_disp_t%.1f.png', prefix, eval, gr, t));
saveas(gcf, sprintf('%s%s%s_disp_t%.1f.fig', prefix, eval, gr, t));    

%%
toc
diary off
