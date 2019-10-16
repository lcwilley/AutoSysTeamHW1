clear
clc

% rng(0); % Fix the random number generator for debugging

dt = 0.1;
Tf = 20;
time = 0:dt:Tf;
N = length(time);
% Field of operation is 20m by 20m, presumably 10 m from 0 in each
% direction.

% Robot initial conditions
x0 = -5; %m
y0 = -3; %m
th0 = pi/2; %rad, 90 degrees

% Create linear and angular velocity commands
vc = 1 + 0.5*cos(2*pi*0.2*time);
wc = -0.2 + 2*cos(2*pi*0.6*time);
u = [vc; wc];

% Create noisy experienced velocities
alph = [0.1, 0.01, 0.01, 0.1];
v = vc + sqrt(alph(1)*vc.^2+alph(2)*wc.^2).*randn(1,N);
w = wc + sqrt(alph(3)*vc.^2+alph(4)*wc.^2).*randn(1,N);

% Initialize truth data vectors, as well as observations
x = [x0, zeros(1,N-1)];
y = [y0, zeros(1,N-1)];
th = [th0, zeros(1,N-1)];

% Fill in truth data based on experienced velocities
for t = 2:N
    x(t) = x(t-1) - v(t)/w(t)*sin(th(t-1)) + v(t)/w(t)*sin(th(t-1)+w(t)*dt);
    y(t) = y(t-1) + v(t)/w(t)*cos(th(t-1)) - v(t)/w(t)*cos(th(t-1)+w(t)*dt);
    th(t) = th(t-1) + w(t)*dt;
end

% For Grading (includes x,y,th, v, om, t)
load('hw2_soln_data.mat')
X = [x;y;th];
w = om;

% Problem-defined range and bearing noise
sig_r = 0.1; %m
sig_th = 0.05; %rad

% Define landmark positions
NL = 3; % Number of landmarks
L = zeros(2,NL);
L(:,1) = [6;4];
L(:,2) = [-7;8];
L(:,3) = [6;-4];
% L = [L [-9;-9]];

% Initialize observation data
z = zeros(2,NL,N);
z_true = zeros(2,NL,N);
for t = 1:N
    z_true(:,:,t) = [sqrt((L(1,:)-x(t)).^2+(L(2,:)-y(t)).^2);
        atan2(L(2,:)-y(t),L(1,:)-x(t))-th(t)];
    z(:,:,t) = [sqrt((L(1,:)-x(t)).^2+(L(2,:)-y(t)).^2)+sig_r*randn(1,NL);
        atan2(L(2,:)-y(t),L(1,:)-x(t))-th(t)+sig_th*randn(1,NL)];
end

% Initialize mu and sigma, then the extended Kalman filter
mu0 = [x0;y0;th0];
sig0 = [1 0 0; 0 1 0; 0 0 0.1];
ekf = EKF(mu0,sig0,dt,sig_r,sig_th,alph,L,N);

% Run the states through the EKF
for t = 2:N
    ekf.update(t,u(:,t),z(:,:,t));
end


%% Plotting
[mu_h, sig_h, mu_bar_h, sig_bar_h, gains] = ekf.get_estimates();

labels = {'X Position', 'Y Position', 'Heading', 'X Error', 'Y Error',...
    'Heading Error'};

% Plot the robot animation
robotPlot = mobileRobotVis(X(:,1),L,dt,Tf);
for t = 2:N
    robotPlot.updatePlot(X(:,t));
    pause(0.01)
end

% Plot true and estimated states versus time
figure(2)
for i = 1:3
    subplot(3,1,i); hold on;
    plot(time,X(i,:),'color',[0 0.6588 0.8039]);
    plot(time,mu_h(i,:),'color',[0 0.3098 0.4196],'LineStyle','--');
    xlabel('Time'); ylabel(labels(i))
    legend('Actual','Measured');
end

figure(3)
% Plot error and covariance over time
for i = 1:3
    subplot(3,1,i); hold on;
    plot(time,X(i,:)-mu_h(i,:),'color',[0 0.6588 0.8039]);
    plot(time,2*sqrt(sig_h(i,:)),'color',[0 0.3098 0.4196])
    plot(time,-2*sqrt(sig_h(i,:)),'color',[0 0.3098 0.4196])
    xlabel('Time'); ylabel(labels(3+i))
    legend('Error','95% Confidence');
end

var_labels = {'X ', 'Y ', 'Theta '};
meas_labels = {'Range Gain', 'Bearing Gain'};

% Plot gains over time
figure(4)
for i = 1:3
    for j = 1:2
        subplot(3,2,2*(i-1)+j)
        plot(time,reshape(gains(i,j,:),[1,N]));
        xlabel('Time'); ylabel(strcat(var_labels(i),meas_labels(j)));
    end
end

% Plot actual range/bearing with measurments
% figure(3)
% for i = 1:3
%     subplot(2,3,2*i-1)
%     plot(time,reshape(z_true(1,i,:),[1,N]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(z(1,i,:),[1,N]),'LineWidth',0.5);
%     subplot(2,3,2*i)
%     plot(time,reshape(z_true(2,i,:),[1,N]),'LineWidth',2);
%     hold on;
%     plot(time,reshape(z(2,i,:),[1,N]),'LineWidth',0.5);
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



