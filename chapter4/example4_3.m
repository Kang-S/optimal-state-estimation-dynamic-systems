% This program shows how a linear oscillator is used to examine 
% the minimum variance properties of the Consider Kalman Filter 
% on dynamic systems.

% Optimal Estimation of Dynamic Systems (2nd ed.) by Crassidis and Junkins
% Example 4.3

% Written by Drew P. Woodbury 2/11

% Preamble
clear; clc; close all;

set(0,'DefaultAxesFontName','Cambria')
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultTextInterpreter','latex')

% Initialization

% Time
t0 = 0;
tstep = 0.1;
tf = 10;
time = t0:tstep:tf;

N = length(time);       % Define size of vectors

omega = 1;              % Oscillator's natural frequency (rad/s)
nTests = 1000;          % Number of Monte Carlo runs

% Initial truth and estimate values
x0true = 2.3;           % State 1
v0true = 0.2;           % State 2
btrue = 0.04;           % Parameter

x0est = 3;
v0est = 0;
best = 0.08;

% Define measurement noise and variance
sigmax = 0.1;
sigmav = 0.1;
R = diag([sigmax^2 sigmav^2]);
invR = inv(R);

% Generate minimum variance weight matrices
W1 = diag(1/sigmax^2*ones(N,1));
W2 = diag([1/sigmax^2*ones(N,1);1/sigmav^2*ones(N,1)]);
W3 = W2;

% Generate initial covariance matrices
% 1-sigma values is error between initial truth and estimated values
Pxx0 = diag([(x0est - x0true)^2, (v0est - v0true)^2]);
invPxx0 = inv(Pxx0);
Pxp0 = [0;0];
Ppp0 = (best - btrue)^2;
invPpp0 = inv(Ppp0);

% Preallocate memory to increase speed
xTLS1err = zeros(2,nTests);
xTLS2err = zeros(2,nTests);
xCLSerr = zeros(2,nTests);

% Main Loop
for k = 1:nTests
    % Generate Measurements
    
    xtrue = zeros(2,N);
    ymeas = zeros(2,N);

    Hx1 = [1 0];
    Hx1Big = zeros(N,2);
    Hx2Big = zeros(2*N,2);

    for i = 1:N

        Phi = [cos(omega*time(i)), 1/omega*sin(omega*time(i));...
            -omega*sin(omega*time(i)), cos(omega*time(i))];

        xtrue(:,i) = Phi*[x0true; v0true];

        ytrue = xtrue(:,i) + [0;btrue];

        ymeas(:,i) = ytrue + diag([sigmax; sigmav])*randn(2,1);
        
        % Fill batch measurement sensitivity matrices
        Hx1Big(i,:) = Hx1*Phi;

        Hx2Big(i,:) = Hx1Big(i,:);
        Hx2Big(N+i,:) = [0 1]*Phi;

    end

    Hx3Big = Hx2Big;
    Hp3Big = [zeros(N,1); ones(N,1)];

    % Scenario 1 - Position Measurements Only, Traditional Least Squares

    PxxTLS1 = inv(Hx1Big'*W1*Hx1Big + invPxx0);
    x0estTLS1 = PxxTLS1*(Hx1Big'*W1*ymeas(1,:)' + invPxx0*[x0est;v0est]);

    xTLS1err(:,k) = x0estTLS1 - [x0true;v0true];

    % Scenario 2 - Position and Velocity Measurements, Traditional Least Squares

    PxxTLS2 = inv(Hx2Big'*W2*Hx2Big + invPxx0);
    
    x0estTLS2 = PxxTLS2*(Hx2Big'*W2*[ymeas(1,:) ymeas(2,:)]' + invPxx0*[x0est;v0est]);

    xTLS2err(:,k) = x0estTLS2 - [x0true;v0true];

    % Scenario 3 - Position and Velocity Measurements, Consider Least Squares
    
    Mxx = inv(Pxx0 - Pxp0*invPpp0*Pxp0');
    Mpp = inv(Ppp0 - Pxp0'*invPxx0*Pxp0);
    Mxp = -Mxx*Pxp0*invPpp0;

    PxxCLS = inv(Mxx + Hx3Big'*W3*Hx3Big - (Hx3Big'*W3*Hp3Big + Mxp)*inv(Hp3Big'*W3*Hp3Big + Mpp)*(Hp3Big'*W3*Hx3Big + Mxp'));
    PxpCLS = -PxxCLS*(Hx3Big'*W3*Hp3Big + Mxp)*inv(Hp3Big'*W3*Hp3Big + Mpp);

    x0estCLS = (PxxCLS*Hx3Big' + PxpCLS*Hp3Big')*W3*[ymeas(1,:) ymeas(2,:)]' + (PxxCLS*Mxx + PxpCLS*Mxp')*[x0est;v0est] + (PxxCLS*Mxp + PxpCLS*Mpp)*best;

    xCLSerr(:,k) = x0estCLS - [x0true;v0true];
    
end

% Means and Standard Deviations

% Position
muTLS1x = mean(xTLS1err(1,:));
stdTLS1x = std(xTLS1err(1,:));
muTLS2x = mean(xTLS2err(1,:));
stdTLS2x = std(xTLS2err(1,:));
muCLSx = mean(xCLSerr(1,:));
stdCLSx = std(xCLSerr(1,:));

positionTable = {'Position (x1)', 'Scenario 1', 'Scenario 2', 'Scenario 3';...
    'Mean', muTLS1x, muTLS2x, muCLSx;...
    'Stan. Dev.', stdTLS1x, stdTLS2x, stdCLSx};

% Velocity
muTLS1v = mean(xTLS1err(2,:));
stdTLS1v = std(xTLS1err(2,:));
muTLS2v = mean(xTLS2err(2,:));
stdTLS2v = std(xTLS2err(2,:));
muCLSv = mean(xCLSerr(2,:));
stdCLSv = std(xCLSerr(2,:));

velocityTable = {'Velocity (x2)', 'Scenario 1', 'Scenario 2', 'Scenario 3';...
    'Mean', muTLS1v, muTLS2v, muCLSv;...
    'Stan. Dev.', stdTLS1v, stdTLS2v, stdCLSv};

disp(positionTable)
disp(velocityTable)

% Plots

figure
hist(xTLS1err(1,:))
hold on; plot(3*sqrt(PxxTLS1(1,1))*ones(1,2),ylim,'r-',-3*sqrt(PxxTLS1(1,1))*ones(1,2),ylim,'r-'); hold off;
xlabel('Postion Estimate Error'); ylabel('Number of Occurrences')
legend('Position Error','3$\sigma$ Covariance Bound')
title('Traditional Least Squares ($x_1$) - Measuring Position Only')

figure
hist(xTLS1err(2,:))
hold on; plot(3*sqrt(PxxTLS1(2,2))*ones(1,2),ylim,'r-',-3*sqrt(PxxTLS1(2,2))*ones(1,2),ylim,'r-'); hold off;
xlabel('Velocity Estimate Error'); ylabel('Number of Occurrences')
legend('Velocity Error','3$\sigma$ Covariance Bound')
title('Traditional Least Squares ($x_2$) - Measuring Position Only')

figure
hist(xTLS2err(1,:))
hold on; plot(3*sqrt(PxxTLS2(1,1))*ones(1,2),ylim,'r-',-3*sqrt(PxxTLS2(1,1))*ones(1,2),ylim,'r-'); hold off;
xlabel('Position Estimate Error'); ylabel('Number of Occurrences')
legend('Position Error','3$\sigma$ Covariance Bound')
title('Traditional Least Squares ($x_1$) - Measuring Position and Velocity')

figure
hist(xTLS2err(2,:))
hold on; plot(3*sqrt(PxxTLS2(2,2))*ones(1,2),ylim,'r-',-3*sqrt(PxxTLS2(2,2))*ones(1,2),ylim,'r-'); hold off;
xlabel('Velocity Estimate Error'); ylabel('Number of Occurrences')
legend('Velocity Error','3$\sigma$ Covariance Bound')
title('Traditional Least Squares ($x_2$) - Measuring Position and Velocity')

figure
hist(xCLSerr(1,:))
hold on; plot(3*sqrt(PxxCLS(1,1))*ones(1,2),ylim,'r-',-3*sqrt(PxxCLS(1,1))*ones(1,2),ylim,'r-'); hold off;
xlabel('Position Estimate Error'); ylabel('Number of Occurrences')
legend('Position Error','3$\sigma$ Covariance Bound')
title('Consider Least Squares ($x_1$) - Measuring Position and Velocity')

figure
hist(xCLSerr(2,:))
hold on; plot(3*sqrt(PxxCLS(2,2))*ones(1,2),ylim,'r-',-3*sqrt(PxxCLS(2,2))*ones(1,2),ylim,'r-'); hold off;
xlabel('Velocity Estimate Error'); ylabel('Number of Occurrences')
legend('Velocity Error','3$\sigma$ Covariance Bound')
title('Consider Least Squares ($x_2$) - Measuring Position and Velocity')