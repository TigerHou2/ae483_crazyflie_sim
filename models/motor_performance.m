%% motor_performance.m
%
% Performs polynomial fitting to find the relation
%   between PWM, RPM, and thurst.

close all
clear;clc

addpath('..\fcns')

thrust = [  0 1.6 4.8 7.9 10.9 13.9 17.3 21 ...
            24.4 28.6 32.8 37.3 41.7 46 51.9 57.9  ];
PWM = [ 0 6.25 12.5 18.75 25 31.25 37.5 43.25 ...
        50 56.25 62.5 68.75 75 81.25 87.5 93.75  ] / 100 * 65535;
RPM = [ 0 4485 7570 9374 10885 12277 13522 14691 ...
        15924 17174 18179 19397 20539 21692 22598 23882  ];

PWM_vect = linspace(0,65535,1000);
    
[p,S] = polyfit(PWM,RPM,3);
Rsq = 1 - (S.normr/norm(RPM - mean(RPM)))^2;
disp(['R^2 = ' num2str(Rsq)])

figure(1)
plot(PWM,RPM,'ro','MarkerSize',6,'LineWidth',1.5)
hold on
plot(PWM_vect,polyval(p,PWM_vect),'b-.','Linewidth',1.25)
hold off
xlabel('PWM')
ylabel('RPM')
latexify(16,9,16)
legend('Test Data', 'Polyfit', 'Location', 'Best')
setgrid