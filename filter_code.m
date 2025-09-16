%% EKF for differential-drive robot with XY position measurements
clear; clc;

% ------------ user inputs ------------
file = '\\nas01.itap.purdue.edu\puhome\My Documents\AAE - ECL Lab\Lab 1\Data\vr_008_vl_010_t20.csv';   % your CSV
L    = 0.25;                      % wheelbase [m]  <-- set correctly
vr   = 0.08;                      % right wheel speed [m/s]
vl   = 0.10;                      % left  wheel speed [m/s]

sigma_meas_xy = 0.005;            % meas. std dev in X,Y [m]
sigma_vr      = 0.01;             % wheel-speed std dev [m/s]
sigma_vl      = 0.01;             % wheel-speed std dev [m/s]
P0 = diag([0.01^2, 0.01^2, (15*pi/180)^2]);   % initial covariance
% -------------------------------------

% Load and pre-process
T = readtable(file, 'VariableNamingRule','preserve');
t  = T.("time[s]");
xm = T.("wr_x[mm]")*1e-3;   % m
ym = T.("wr_y[mm]")*1e-3;   % m
N  = numel(t);

dt = [diff(t); 0];
good = dt > 0;
if any(~good(1:end-1)), dt(~good) = median(dt(good)); end
dt(end) = dt(end-1);

% State: x = [X; Y; theta]
xhat = zeros(3,N);
P    = zeros(3,3,N);

% Initialize from first measurement(s)
xhat(1,1) = xm(1);
xhat(2,1) = ym(1);
if N >= 2
    xhat(3,1) = atan2(ym(2)-ym(1), xm(2)-xm(1));
else
    xhat(3,1) = 0;
end
P(:,:,1) = P0;

% Const matrices
H = [1 0 0;
     0 1 0];
R = (sigma_meas_xy^2)*eye(2);
M = diag([sigma_vr^2, sigma_vl^2]);

for k = 2:N
    dtk = dt(k-1);
    th  = xhat(3,k-1);

    % Controls
    V     = 0.5*(vr + vl);
    omega = (vr - vl)/L;

    % ---- Predict (Euler discretization) ----
    xpred = xhat(:,k-1);
    xpred(1) = xpred(1) + V*cos(th)*dtk;
    xpred(2) = xpred(2) + V*sin(th)*dtk;
    xpred(3) = xpred(3) + omega*dtk;
    xpred(3) = atan2(sin(xpred(3)), cos(xpred(3)));   % wrap to [-pi,pi]

    Fk = [1, 0, -V*sin(th)*dtk;
          0, 1,  V*cos(th)*dtk;
          0, 0,  1];

    Gk = [0.5*cos(th)*dtk, 0.5*cos(th)*dtk;
          0.5*sin(th)*dtk, 0.5*sin(th)*dtk;
          (1/L)*dtk,      -(1/L)*dtk];

    Qk = Gk*M*Gk';
    %Qk = eye(3);
    Ppred = Fk*P(:,:,k-1)*Fk' + Qk;

    % ---- Update with position measurement ----
    zk     = [xm(k); ym(k)];
    ytilde = zk - H*xpred;
    S      = H*Ppred*H' + R;
    K      = (Ppred*H')/S;            % same as Ppred*H'*inv(S)

    xupd   = xpred + K*ytilde;
    Pupd   = (eye(3) - K*H)*Ppred;

    xhat(:,k) = xupd;
    P(:,:,k)  = Pupd;
end

% Save & plot
results = table(t, xhat(1,:)', xhat(2,:)', xhat(3,:)', ...
    'VariableNames', {'time_s','x_hat_m','y_hat_m','theta_hat_rad'});
writetable(results, 'ekf_results.csv');   % -> ekf_results.csv

figure; hold on;
plot(xm, ym, 'LineWidth',1.5, 'DisplayName','Measured (X,Y)');
plot(xhat(1,:), xhat(2,:), 'LineWidth',1.5, 'DisplayName','EKF estimate');
axis equal; grid on;
xlabel('X [m]'); ylabel('Y [m]');
title('Measured vs EKF-estimated trajectory');
legend('Location','best');
