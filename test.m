%% Clear workspace
close all
clear;clc

addpath('fcns')
addpath('models')

%% Find equations of motion
[h, x, u, h_num, params] = findeoms();

%% Define waypoints and starting position
% initial state in simulation frame (X-forward, Z-down)
x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
% target state in Crazyflie frame (X-forward, Z-up)
XE = [[0; 0; 0.5],...
      [1; 0; 1],...
      [1; 1; 0.5],...
      [0; 1; 1],...
      [0; 0; 0.5]];
XE = [XE;zeros(9,size(XE,2))];

xTol = 0.05; % hover position tolerance
tMax = 8; % max simulation time

% define waypoint hover times; also chooses the flight video to compare
tHover = [1,1,1,1,1]; physDataName = '04A03';
% tHover = [1,0,0,0,1]; physDataName = '03B02';

%% Initialize data structures
% depending on whether the PID controller is reset after each waypoint
% we initialize the data structure a little differently
resetPID_atWaypoint = false;
if resetPID_atWaypoint
    data = struct('t',[],'x',[],'u',[],'e',[],...
                  'int',[],'m',[],'xd',[],'xn',[]);
else
    data = [];
end
% timestamps when a new waypoint is set
% i.e. the previous waypoint has been satisfied
switch_times = [];

%% Iterate through waypoints
t0 = 0;
for i = 1:size(XE,2)
    xe = XE(:,i); % update waypoint
    tStab = tHover(i); % update waypoint hover time
    if i ~= 1
        x0 = data.x(:,end); % update initial state
        t0 = data.t(end); % update initial time
    end
    % simulate
    [data_wp] = simulate(h_num,params,xe,x0,xTol,t0,tStab,tMax,[]);
    % store data
    if resetPID_atWaypoint
        data.t = [data.t(:,1:end-1), data_wp.t];
        data.x = [data.x(:,1:end-1), data_wp.x];
        data.u = [data.u(:,1:end-1), data_wp.u];
        data.e = [data.e(:,1:end-1), data_wp.e];
        data.int = [data.int(:,1:end-1), data_wp.int];
        data.m = [data.m(:,1:end-1), data_wp.m];
        data.xd = [data.xd(:,1:end-1), data_wp.xd];
        data.xn = [data.xn(:,1:end-1), data_wp.xn];
    else
        data = simulate(h_num,params,xe,x0,xTol,t0,tStab,tMax,data);
    end
    % store waypoint starting timestamp
    switch_times(i) = t0;
end
% append simulation end time
switch_times(i+1) = max(data.t);

%% Show video

% Parse data from simulation
t_s = data.t; % time, s
o_s = data.x(1:3, :); % position vector, m
hy_s = data.x(4, :); % yaw, rad
hp_s = data.x(5, :); % pitch, rad
hr_s = data.x(6, :); % roll, rad
v = data.x(7:9, :); % velocity vector, m/s
w = data.x(10:12, :); % angular velocity vector, rad/s

% could give the name of a file to save a movie
% of the visualization to, like 'test'
% the file format extension is automatically added
videoPath = 'videos/';
videoName = [];
if isempty(videoName)
    moviefile = [];
else
    moviefile = [videoPath videoName '.mp4'];
end

% scale o to make quadcopter look smaller
O_s = 6*o_s;

dirpath = ['C:\Users\tigre\Desktop\Directory\UIUC' ...
           '\ae483\crazyflie_python_controls\lab3_data_mat'];
dat = open([dirpath '\dat_' physDataName '.mat']);

dat.switches(end+1) = max(dat.t);
switches = dat.switches;

o_r = [dat.px; dat.py; dat.pz] .* [1; -1; -1];
O_r = 6*o_r;
O_r(3,:) = O_r(3,:);
t_r = dat.t;
hy_r = zeros(size(dat.px));
hp_r = hy_r;
hr_r = hy_r;

T = {t_s,t_r};
O = {O_s,O_r};
HY = {hy_s,hy_r};
HP = {hp_s,hp_r};
HR = {hr_s,hr_r};
C = {'c','r'};

close(figure(1))
figure(1);
visualize(T, O, HY, HP, HR, C, moviefile); 
hold off

%% Show path

close(figure(2));
figure(2); hold on; box on; grid on;
oo_s = o_s .* [1;-1;-1];
plot3(oo_s(1,:),oo_s(2,:), oo_s(3,:),'k',...
      'HandleVisibility','Off'); % plot quadrotor trajectory in Z-up frame
plot3(oo_s(1,1:75:end), ...
      oo_s(2,1:75:end), ...
      oo_s(3,1:75:end),'b.','markersize',12,...
      'DisplayName','Simulation'); % plot points along trajectory equally spaced in time

oo_r = o_r .* [1;-1;-1];
plot3(oo_r(1,:),oo_r(2,:), oo_r(3,:),'k',...
      'HandleVisibility','Off'); % plot quadrotor trajectory in Z-up frame
plot3(oo_r(1,1:end), ...
      oo_r(2,1:end), ...
      oo_r(3,1:end),'r.','markersize',12,...
      'DisplayName','Experiment'); % plot points along trajectory equally spaced in time
xlabel('X, m');
ylabel('Y, m');
zlabel('Z, m');
latexify(16,14,15)
expand
axis equal
view(3)
setgrid
legend('Location','Best')
hold off

%% Show time-of-flights
close(figure(3))
figure(3)
hold on
alpha = 0.55;
alm = 1 - alpha;
simColors = [alm, alm, 1.0; ...
             alm, 1.0, alm];
irlColors = [1.0, alm, alm; ...
             1.0, 1.0, alm];
for i = 1:length(switch_times)-1
    x = switch_times(i); y = 0;
    w = switch_times(i+1) - switch_times(i); h = 1;
    c = simColors(mod(i+1,2)+1,:);
    rectangle('Position',[x,y,w,h],'FaceColor',c);
    x = switches(i); y = 1;
    w = switches(i+1) - switches(i); h = 1;
    c = irlColors(mod(i+1,2)+1,:);
    rectangle('Position',[x,y,w,h],'FaceColor',c);
end
annotation('Textbox', [0, 0.37, 0.1, 0.25], ...
           'String', 'Sim', ...
           'HorizontalAlignment', 'right', ...
           'Margin', 0, 'EdgeColor', 'none', ...
           'Interpreter', 'latex')
annotation('Textbox', [0, 0.64, 0.1, 0.25], ...
           'String', 'Real', ...
           'HorizontalAlignment', 'right', ...
           'Margin', 0, 'EdgeColor', 'none', ...
           'Interpreter', 'latex')
xlabel('Time, s')
set(gca,'TickDir','out');
set(gca,'XMinorTick','on')
set(gca,'YTickLabel',[]);
xlim([0,floor(max(switch_times(end),switches(end)))+1])
hold off
latexify(16,4)
expand(0.1,0.02,0,0.02)

%% Show Axis Comparisons
close(figure(4))
figure(4)

simStyle = 'b-';
expStyle = 'r-.';

subplot(3,1,1)
hold on
plot(data.t,data.x(1,:),simStyle,'LineWidth',1.5)
plot(dat.t,dat.px,expStyle,'LineWidth',1.5)
ylabel('X, m')
hold off
setgrid
legend('Simulation','Experiment','Location','Best')

subplot(3,1,2)
hold on
plot(data.t,-data.x(2,:),simStyle,'LineWidth',1.5)
plot(dat.t,dat.py,expStyle,'LineWidth',1.5)
ylabel('Y, m')
hold off
setgrid

subplot(3,1,3)
hold on
plot(data.t,-data.x(3,:),simStyle,'LineWidth',1.5)
plot(dat.t,dat.pz,expStyle,'LineWidth',1.5)
ylabel('Z, m')
xlabel('Time, s')
hold off
setgrid

latexify(24,14,15)