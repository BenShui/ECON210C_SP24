% ECON 210C Homework 1
% Author: Jiahui Shui

%% Parameters
gam = 1;
phi = 1;
chi = 1;
beta = 0.99;
rho = 0.99;
sigma_m = 0.01;
nu_arr = [0.25, 0.5, 0.999, 2, 4];

% time
T = 500;

%% Calibration
nu_len = length(nu_arr);
theta_arr = zeros(1, nu_len);
for i = 1:nu_len
    nu = nu_arr(i);
    f = @(theta) (1-beta)^(-1/nu)*(theta/(1-theta))^(1/nu)*((1-theta)/chi*((1-theta+theta*(1-beta)^(-(1-nu)/nu)*(theta/(1-theta))^((1-nu)/nu)))^((nu-gam)/(1-nu)))^(1/(phi+gam))-1;
    theta_arr(i) = fzero(f, 0.01);
end


%% Setup Matrices
I = sparse(eye(T));
O = sparse(zeros(T,T));
Ip1 = sparse(diag(ones(1, T-1), 1));
Im1 = sparse(diag(ones(1, T-1), -1));

%% Generating m
m = zeros(T,1);
m(1) = 1;
for j=2:T
    m(j) = rho*m(j-1);
end

%% Iteration
for i=1:nu_len
    nu = nu_arr(i);
    theta = theta_arr(i);
    % Market Clearing Condition Block
    Phigmc = I;
    Phigmy = -I;
    Phigmq = O; Phigmx = O; Phigmw = O;
    
    Phieulc = nu*I-nu*Ip1;
    Phieulx = (gam-nu)*I-(gam-nu)*Ip1;
    Phieulq = I;
    Phieuly = O;
    Phieulw = O;
 
    dHdY = [Phigmc, Phigmq, Phigmx, Phigmy, Phigmw; Phieulc, Phieulq, Phieulx, Phieuly, Phieulw];
    % There's extra term dHdU here
    Phigmn = O;
    Phieuln = O;
    Phigmp = O;
    Phieulp = I-Ip1;
    dHdU = [Phigmn, Phigmp; Phieuln, Phieulp];

    % Firm Block
    Phiyn = I;
    Phiyp = O;
    Phiyep = O;
    Phiwp = I;
    Phiwn = O;
    Phiwep = O;
    dYFdU = [Phiyn, Phiyp;Phiwn, Phiwp];
    dYFdZ = [Phiyep; Phiwep];

    % Households Block
    omega = (theta^(1/nu)*(1-beta)^(1-1/nu))/((1-theta)^(1/nu)+theta^(1/nu)*(1-beta)^(1-1/nu));
    eta = beta/(nu*(1-beta));
    Phicn = -eta*phi/(gam*eta-(gam-nu)*omega*eta)*I;
    Phicp = (gam-nu)*omega*eta/(gam*eta-(gam-nu)*omega*eta)*I;
    Phiqn = -phi/(gam*eta-(gam-nu)*omega*eta)*I;
    Phiqp = gam/(gam*eta-(gam-nu)*omega*eta)*I;
    Phicm = -(gam-nu)*omega*eta/(gam*eta-(gam-nu)*omega*eta)*I;
    Phiqm = -gam/(gam*eta-(gam-nu)*omega*eta)*I;
    Phixn = (1-omega)*Phicn;
    Phixp = (1-omega)*Phicp - omega*I;
    Phixm = (1-omega)*Phicm + omega*I;
    dYHdU = [Phicn, Phicp;Phiqn, Phiqp; Phixn, Phixp];
    dYHdZ = [Phicm; Phiqm; Phixm];

    % Stack dYdU and dYdZ
    dYdU = [dYHdU; dYFdU];
    dYdZ = [dYHdZ; dYFdZ];

    % Compute dHdU and dHdZ
    dHdU = dHdY * dYdU + dHdU;
    dHdZ = dHdY * dYdZ;

    % Jacobian
    dUdZ = -inv(dHdU)*dHdZ;
    dYdZ = dYdU*dUdZ+dYdZ;
    dXdZ = [dUdZ;dYdZ];

    % Compute IRF
    X = dXdZ * m;
    B = reshape(X, [T, 7]);
    
    % Unpack B
    n = B(: ,1);
    p = B(:, 2);
    c = B(:, 3);
    q = B(:, 4);
    x = B(:, 5);
    y = B(:, 6);
    w = B(:, 7);
    
    % Plot c
    figure(1);
    plot(1:length(c), c*100, '-', 'Linewidth', 1.5);
    hold on

    figure(2);
    plot(1:length(p), p*100, '-', 'Linewidth', 1.5);
    hold on

    figure(3);
    plot(1:length(q), q*100, '-', 'Linewidth', 1.5);
    hold on
end

%% Check Sum
sum(abs(y-n)<1e-5)/T == 1
sum(abs(w-p)<1e-5)/T == 1
sum(abs(gam*c+phi*n+(gam-nu)*omega*eta*q)<1e-5)/T==1
sum(abs(m-p-c+eta*q)<1e-5)/T == 1
sum(abs(y-c)<1e-5)/T==1
sum(abs(q+Phieulp*p+Phieulc*c+Phieulx*x)<1e-5)/T==1
sum(abs(x-(1-omega)*c-omega*(m-p))<1e-5)/T==1

%% Finish Plot
figure(1)
hold off
legend(strcat('\nu = ',num2str(nu_arr')))
title('Consumption IRF')
xlabel('Time')
ylabel('Consumption Change (%)')
grid on

figure(2)
hold off
legend(strcat('\nu = ',num2str(nu_arr')))
title('Price IRF')
xlabel('Time')
ylabel('Price Change (%)')
grid on

figure(3)
hold off
legend(strcat('\nu = ',num2str(nu_arr')))
title('Nominal Interest Rate IRF')
xlabel('Time')
ylabel('Nominal Interest Rate Change (%)')
grid on