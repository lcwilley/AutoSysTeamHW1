clear
clc

% rng(0); % Fix the random number generator for debugging

dt = 0.1;
Tf = 20;
time = 0:dt:Tf;
NN = length(time);
% Field of operation is 20m by 20m, presumably 10 m from 0 in each
% direction.

% Robot initial conditions
x0 = -5; %m
y0 = -3; %m
th0 = pi/2; %rad, 90 degrees

% Create linear and angular velocity commands
vc = 1 + 0.5*cos(2*pi*0.2*time);
wc = -0.2 + 2*cos(2*pi*0.6*time);
uu = [vc; wc];

% Create noisy experienced velocities
alph = [0.1, 0.01, 0.01, 0.1];
vv = vc + sqrt(alph(1)*vc.^2+alph(2)*wc.^2).*randn(1,NN);
ww = wc + sqrt(alph(3)*vc.^2+alph(4)*wc.^2).*randn(1,NN);

% Initialize truth data vectors, as well as observations
xx = [x0, zeros(1,NN-1)];
yy = [y0, zeros(1,NN-1)];
th = [th0, zeros(1,NN-1)];

% Fill in truth data based on experienced velocities
for tt = 2:NN
    xx(tt) = xx(tt-1) - vv(tt)/ww(tt)*sin(th(tt-1)) + vv(tt)/ww(tt)*sin(th(tt-1)+ww(tt)*dt);
    yy(tt) = yy(tt-1) + vv(tt)/ww(tt)*cos(th(tt-1)) - vv(tt)/ww(tt)*cos(th(tt-1)+ww(tt)*dt);
    th(tt) = th(tt-1) + ww(tt)*dt;
end

% For Grading (includes xx,yy,th, vv, om, tt)
% load('hw2_soln_data.mat')
XX = [xx;yy;th];
% ww = om;

% Problem-defined range and bearing noise
sig_r = 0.1; %m
sig_th = 0.05; %rad
sig_rth = [sig_r; sig_th];

% Define landmark positions
NL = 3; % Number of landmarks
Lm = zeros(2,NL);
Lm(:,1) = [6;4];
Lm(:,2) = [-7;8];
Lm(:,3) = [6;-4];
% Lm = [Lm [-9;-9]];

% Initialize observation data
zz = zeros(2,NL,NN);
z_true = zeros(2,NL,NN);
% keyboard
for tt = 1:NN
    lm_sub_xy = Lm - XX(1:2,tt);
    z_true(:,:,tt) = [sqrt(sum(lm_sub_xy.^2, 1));
            atan2(lm_sub_xy(2,:),lm_sub_xy(1,:)) - XX(3,tt)];
    %
    z_true(2,:,tt) = rad_wrap_pi(z_true(2,:,tt));

    noise = sig_rth .* randn(2,NL);
    %
    zz(:,:,tt) = zz(:,:,tt) + noise;
    %

end

% Initialize mu and sigma, then the extended Kalman filter
mu0 = [x0;y0;th0];
sig0 = [1 0 0; 0 1 0; 0 0 0.1];
ekf = EKF(mu0,sig0,dt,sig_r,sig_th,alph,Lm,NN);

% Run the states through the EKF
for tt = 2:NN
    ekf.update(tt,uu(:,tt),zz(:,:,tt));
end


%% Plotting
[mu_h, sig_h, mu_bar_h, sig_bar_h, gains] = ekf.get_estimates();

labels = {'X Position', 'Y Position', 'Heading', 'X Error', 'Y Error',...
    'Heading Error'};

% Plot the robot animation
robotPlot = mobileRobotVis(XX(:,1),Lm,dt,Tf);
for tt = 2:NN
    robotPlot.updatePlot(XX(:,tt));
    pause(0.001)
end

% Plot true and estimated states versus time
figure(2)
for ii = 1:3
    subplot(3,1,ii); hold on;
    plot(time,XX(ii,:),'color',[0 0.6588 0.8039]);
    plot(time,mu_h(ii,:),'color',[0 0.3098 0.4196],'LineStyle','--');
    xlabel('Time'); ylabel(labels(ii))
    legend('Actual','Measured');
end

figure(3)
% Plot error and covariance over time
for ii = 1:3
    subplot(3,1,ii); hold on;
    plot(time,XX(ii,:)-mu_h(ii,:),'color',[0 0.6588 0.8039]);
    plot(time,2*sqrt(sig_h(ii,:)),'color',[0 0.3098 0.4196])
    plot(time,-2*sqrt(sig_h(ii,:)),'color',[0 0.3098 0.4196])
    xlabel('Time'); ylabel(labels(3+ii))
    legend('Error','95% Confidence');
end

var_labels = {'X ', 'Y ', 'Theta '};
meas_labels = {'Range Gain', 'Bearing Gain'};

% Plot gains over time
figure(4)
for ii = 1:3
    for jj = 1:2
        subplot(3,2,2*(ii-1)+jj)
        plot(time,reshape(gains(ii,jj,:),[1,NN]));
        xlabel('Time'); ylabel(strcat(var_labels(ii),meas_labels(jj)));
    end
end

% Plot actual range/bearing with measurments
% figure(3)
% for ii = 1:3
%     subplot(2,3,2*ii-1)
%     plot(time,reshape(z_true(1,ii,:),[1,NN]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(zz(1,ii,:),[1,NN]),'LineWidth',0.5);
%     subplot(2,3,2*ii)
%     plot(time,reshape(z_true(2,ii,:),[1,NN]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(zz(2,ii,:),[1,NN]),'LineWidth',0.5);
% end

% Create a visualization of the robot in the horizontal plane
%   Why doesn't the velocity motion model work for straight line motion?
%   Because
% Simulate the range and bearing measurements
% Implement the EKF algorithm, Plot the following:
%   True and estimated states versus time
%   Error and 95% uncertainty versus time
%	Gains versus time
% Change input velocity, landmark locations, sensor noise, and control
%  noise
% Increase and decrease the number of landmarks
% EXTRA: Modify the simulation and EKF to handle curved and straight-line
%  motions
% The discrete equations are simple without turning. With the curved
%  velocities, it has omega in the denomenator, and you have to look at
%  conditional statements to make sure that you deal with weird conditions
