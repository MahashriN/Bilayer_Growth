clear all; close all; clc;
tic

%% Model selection : set exactly one kinetics and one growth law equal to 1.

Schnakenberg = 1;
Brusselator  = 0;

linear       = 1;
exponen      = 0;
logistic     = 0; 

k_eta        = 0;
single_eta   = 1;
TDHSS_graph  = 0; % optional

%% Model and simulation parameters

% Coupling and growth parameters
eta = 2;
gamma1 = 0.05;
gamma2 = 0.05;

% Logistic saturation setting
Linf1 = 0;
Linf2 = 0;

T = 15000;

Kmax = 200;

t_snaps = linspace(0, T, 4);

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

folder = fullfile(pwd,'output',kineticsName,[growthName '_SA']);

if ~exist(folder,'dir')
    mkdir(folder);
end

%% Run information

gr1=num2str(round(gamma1,5));
gr2=num2str(round(gamma2,5));

etaStr=num2str(round(eta,4));

runID = datestr(datetime('now'),'yyyymmdd_HHMMSS');

prefix = sprintf('%s_%s_eta_%s_g1_%s_g2_%s',...
    runID,...
    growthName,...
    etaStr,...
    gr1,gr2);

if logistic
    prefix = sprintf('%s_Linf1_%s_Linf2_%s',...
        prefix,num2str(Linf1),num2str(Linf2));
end

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
    %
    fu = @(u,v) -1 + 2*u.*v;  fv = @(u,v) u.^2;
    gu = @(u,v) -2*u.*v;      gv = @(u,v) -u.^2;
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
    %
    fu = @(u,v) -(b+1) + 2*u.*v; fv = @(u,v) u.^2;
    gu = @(u,v) b - 2*u.*v;      gv = @(u,v) -u.^2;
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
fprintf(' Frozen-Time Stability Analysis\n');
fprintf('=====================================\n');
fprintf(' Kinetics : %s\n', kineticsName);
fprintf(' Growth   : %s\n', growthName);
fprintf(' eta      : %s\n', etaStr);
fprintf(' gamma1   : %s\n', gr1);
fprintf(' gamma2   : %s\n', gr2);
fprintf('=====================================\n\n');

%% TDHSS: Integrate spatially homogeneous base ODE
tspan = [0, T];
y0 = [u1_0; v1_0; u2_0; v2_0];
if logistic
    hfun1=@(t) h(t,gamma1,Linf1);
    hfun2=@(t) h(t,gamma2,Linf2);
else
    hfun1=@(t) h(t,gamma1);
    hfun2=@(t) h(t,gamma2);
end
baseODE = @(t,y) [
    f(y(1),y(2),y(3),y(4),eta) - hfun1(t)*y(1);
    g(y(1),y(2),y(3),y(4),eta) - hfun1(t)*y(2);
    f(y(3),y(4),y(1),y(2),eta) - hfun2(t)*y(3);
    g(y(3),y(4),y(1),y(2),eta) - hfun2(t)*y(4);
];
[t_base, y_base] = ode45(baseODE, tspan, y0);

fprintf('\nHomogeneous base ODE integrated for 0 <= t <= %.2f.\n',T);

u1_t = y_base(:,1);  v1_t = y_base(:,2);
u2_t = y_base(:,3);  v2_t = y_base(:,4);

%% Optional: homogeneous base-state visualization
if TDHSS_graph
figure('Units','inches','Position',[1 1 6.5 5]); 

% ----------- U-PANEL -----------
subplot(2,1,1)
plot(t_base, u1_t, 'b',  'LineWidth',1.6); hold on
plot(t_base, u2_t, 'm--','LineWidth',1.6);

xlim([0 500]);
ylabel('Concentration','FontSize',12,'Interpreter','latex');
set(gca,'FontSize',11,'LineWidth',1);

lgd = legend('$u_1$','$u_2$');
set(lgd,'Interpreter','latex','FontSize',12,'Location','northeast');

% inset
axes('Position',[0.50 0.70 0.25 0.20]); % inset properly placed
box on
plot(t_base, u1_t,'b','LineWidth',1); hold on
plot(t_base, u2_t,'m--','LineWidth',1);
xlim([100 150]);
set(gca,'FontSize',9,'LineWidth',1);

% ----------- V-PANEL -----------
subplot(2,1,2)
plot(t_base, v1_t, 'r',  'LineWidth',1.6); hold on
plot(t_base, v2_t, 'k--','LineWidth',1.6);

xlim([0 500]);
xlabel('$t$','FontSize',12,'Interpreter','latex');
ylabel('Concentration','FontSize',12,'Interpreter','latex');
set(gca,'FontSize',11,'LineWidth',1);

lgd = legend('$v_1$','$v_2$');
set(lgd,'Interpreter','latex','FontSize',12,'Location','northeast');

% inset
axes('Position',[0.50 0.22 0.25 0.20]);
box on
plot(t_base, v1_t,'r','LineWidth',1); hold on
plot(t_base, v2_t,'k--','LineWidth',1);
xlim([100 150]);
set(gca,'FontSize',9,'LineWidth',1);

exportgraphics(gcf,...
    [runPrefix,'_TDHSS.png'],...
    'Resolution',600);
end

%% Frozen-time dispersion at selected snapshots 

if k_eta
eta_grid = linspace(0,10,400);
k_grid   = linspace(0,Kmax,400);
[kmesh, emesh] = meshgrid(k_grid, eta_grid);

for ts = 1:numel(t_snaps)
    t = t_snaps(ts);
    yb = interp1(t_base, y_base, t);
    u1_base = yb(1); v1_base = yb(2); u2_base = yb(3); v2_base = yb(4);

    J_1 = [fu(u1_base,v1_base)  fv(u1_base);  gu(u1_base,v1_base)  gv(u1_base)];
    J_2 = [fu(u2_base,v2_base)  fv(u2_base);  gu(u2_base,v2_base)  gv(u2_base)];

    % growth + diffusion scalings
    if logistic
        L1 = L(t, gamma1,Linf1);   L2 = L(t, gamma2,Linf2);
        h1 = h(t, gamma1,Linf1);   h2 = h(t, gamma2,Linf2);
    else
        L1 = L(t, gamma1);   L2 = L(t, gamma2);
        h1 = h(t, gamma1);   h2 = h(t, gamma2);
    end
    
    DM = diag([D(1)/L1^2, D(2)/L1^2, D(3)/L2^2, D(4)/L2^2]);
    HM = diag([h1, h1, h2, h2]);
    JM = [J_1, zeros(2); zeros(2), J_2];

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
    shading flat;
    xlabel('$k$','Interpreter','latex'); 
    ylabel('$\eta$','Interpreter','latex'); 
    clim([-0.1 0.1]);
    cmap_colorbar()
    cb = colorbar;
    cb.FontSize = 12;
    set(gca,'FontSize',16);
    
    % Export
    exportgraphics(gcf,...
        sprintf('%s_t_%.1f.png',runPrefix,t),...
        'Resolution',600);
    savefig(gcf,...
        sprintf('%s_t_%.1f.fig',runPrefix,t));
end
end

%% Frozen-time dispersion for a single eta value 
if single_eta
k_grid = linspace(0,Kmax,400);

fprintf('\nComputing frozen-time dispersion: eta = %.2f\n', eta);
figure('Color','w');
for ts = 1:numel(t_snaps)

    t = t_snaps(ts);
    yb = interp1(t_base, y_base, t);
    u1_base = yb(1); v1_base = yb(2); u2_base = yb(3); v2_base = yb(4);

    J_1 = [fu(u1_base,v1_base)  fv(u1_base);  gu(u1_base,v1_base)  gv(u1_base)];
    J_2 = [fu(u2_base,v2_base)  fv(u2_base);  gu(u2_base,v2_base)  gv(u2_base)];

    % growth + diffusion scaling
    if logistic
    L1 = L(t, gamma1,Linf1);   L2 = L(t, gamma2,Linf2);
    h1 = h(t, gamma1,Linf1);   h2 = h(t, gamma2,Linf2);
    else
    L1 = L(t, gamma1);   L2 = L(t, gamma2);
    h1 = h(t, gamma1);   h2 = h(t, gamma2);
    end
    
    DM = diag([D(1)/L1^2, D(2)/L1^2, D(3)/L2^2, D(4)/L2^2]);
    HM = diag([h1, h1, h2, h2]);
    JM = [J_1, zeros(2); zeros(2), J_2];

    Mr = zeros(size(k_grid)); % max Re(lambda(k))

    for i = 1:length(k_grid)
        k = k_grid(i);
        rho = (k*pi)^2;  % Neumann eigenvalue (no domain length scaling)
        CM = eta * [-eye(2) eye(2); eye(2) -eye(2)];

        Mk = -rho*DM - HM + JM + CM;
        ev = eig(Mk);
        Mr(i) = max(real(ev));
    end
    % Diagnostic: dominant unstable modes
    [maxMr,x_pos]=findpeaks(Mr);

    fprintf('\nSnapshot t = %.2f\n',t);
    fprintf('\npeak growthrate    k\n');
    disp([maxMr(:) k_grid(x_pos(:))']);

    %
    plot(k_grid, Mr,'LineWidth', 3,'DisplayName', sprintf('$t = %d$', t)); 
    hold on;
    yline(0,'k--','linewidth',2,'HandleVisibility', 'off');  % stability threshold
    xlabel('$k$','FontSize',16,'Interpreter','latex');
    ylabel('$\max\Re(\lambda)$','FontSize',16,'Interpreter','latex');
    grid on; set(gca,'FontSize',16);
    xlim([0,Kmax]);
end
legend('show', 'Location', 'northeast');
set(legend,'Interpreter','latex','FontSize',14);
ylim([-1.5,1])

% Export
exportgraphics(gcf,...
    sprintf('%s_t_%.1f.png',runPrefix,t),...
    'Resolution',600);

exportgraphics(gcf,...
    sprintf('%s_t_%.1f.pdf',runPrefix,t),...
    'ContentType','vector');

savefig(gcf,...
    sprintf('%s_t_%.1f.fig',runPrefix,t));
end

%%
toc
diary off
