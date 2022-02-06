clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')
global C
global Em T
global BoundaryX BoundaryY

global Vtherm
global nElectrons
global t_mn



C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                    % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s²


Em = 0.26 * C.m_0;                    % Mass of the Electron

T = 300;                            % Semiconductor temperature
            
t_mn = 0.2e-12;                     % Mean time between collisions
                    
TimeSteps = 300;                   % Number of time steps

nElectrons = 1000;                   % Number of electrons

dt = 1e-14;                         % Time Step

Vtherm = sqrt(2 * C.kb * T/Em);

 
Standard_deviation = Vtherm/(sqrt(2));
delta_t = 7.5 * 10 ^ -15;
array_1 = zeros(1, 1000);
tmp_arary_1 = (1:1:1000);
BoundaryX = 200 * 10 ^ -9;
BoundaryY = 100 * 10 ^ -9;
tmn = 0.2 * 10 ^ -12;
Px = rand(1, nElectrons) .* BoundaryX;
Py = rand(1, nElectrons) .* BoundaryY;
Theta = rand(1,nElectrons)*2*pi;
Vx = randn(1,nElectrons) .*Standard_deviation;
Vy = randn(1,nElectrons) .*Standard_deviation;
BoundaryX = 200 * 10 ^ -9;
BoundaryY = 100 * 10 ^ -9;

inbox1 = true;
while inbox1 == true
    inbox = ((Px <= (1.30 * BoundaryX/2) & (Px >= (0.90 * BoundaryX/2))) & ((Py < (BoundaryY/3)) | Py >= (2*BoundaryY/3)));
    if (sum(inbox) > 0)
        Px(inbox) = rand(1, sum(inbox)) .* BoundaryX;
        Py(inbox) = rand(1, sum(inbox)) .* BoundaryY;
    else
        inbox1 = false;
    end
    
end

for i = 1:1000
    
    x_old = Px;
    Px = x_old + delta_t*Vx;
    y_old = Py;
    Py = y_old + delta_t*Vy;
    j1 = (Py >= BoundaryY);
    j2 = (Py < 0);
    Vy(j1) = -Vy(j1);
    Vy(j2) = -Vy(j2);
    j3 = ((Px) >= BoundaryX);
    j4 = ((Px) < 0); %
    x_old(j3) = x_old(j3) - BoundaryX;
    Px(j3) = Px(j3) - BoundaryX;
    x_old(j4) = x_old(j4) + BoundaryX;
    Px(j4) = Px(j4) + BoundaryX;
    
    % if(  inbox==0&& middle==0);
%     inbox = ((Px <= (1.40 * BoundaryX/2) & (Px >= (0.80 * BoundaryX/2))) & ((Py < (BoundaryY/3)) | Py >= (2*BoundaryY/3)));
%     middle = ((x_old<= (1.40 * BoundaryX/2) & (x_old>= (0.80 * BoundaryX/2))) & ((y_old > (BoundaryY/3)) & y_old <= (2*BoundaryY/3)));
%     
%     plot (Px,Py,'b.');
%     pause(0.05)
%     hold on
%     
% end 

    inbox = ((Px <= (1.40 * BoundaryX/2) & (Px >= (0.80 * BoundaryX/2))) & ((Py < (BoundaryY/3)) | Py >= (2*BoundaryY/3)));
    middle = ((x_old<= (1.40 * BoundaryX/2) & (x_old>= (0.80 * BoundaryX/2))) & ((y_old > (BoundaryY/3)) & y_old <= (2*BoundaryY/3)));
    
    plot (Px,Py,'b.');
    pause(0.05)
    hold on
   
    
    %line =0;
    
    % if (line-1 == 0)
%         
%         
%         % Plotting the boxes
%         line([0.80*BoundaryX/2 0.80*BoundaryX/2], [BoundaryY 2*BoundaryY/3]);
%     line([1.40*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY 2*BoundaryY/3]);
%     line([0.80*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY BoundaryY]);
%     line([0.80*BoundaryX/2 1.40*BoundaryX/2], [2*BoundaryY/3 2*BoundaryY/3]);
%         
%      
%     elseif (line-1 < 1000)
%         line([0.80*BoundaryX/2 0.80*BoundaryX/2], [0 BoundaryY/3]);
%     line([1.40*BoundaryX/2 1.40*BoundaryX/2], [0 BoundaryY/3]);
%     line([0.80*BoundaryX/2 1.40*BoundaryX/2], [0 0]);
%     line([0.80*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY/3 BoundaryY/3]);
%     else
%         
%         hold off
%     end

    line([0.80*BoundaryX/2 0.80*BoundaryX/2], [BoundaryY 2*BoundaryY/3]);
    line([1.40*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY 2*BoundaryY/3]);
    line([0.80*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY BoundaryY]);
    line([0.80*BoundaryX/2 1.40*BoundaryX/2], [2*BoundaryY/3 2*BoundaryY/3]);
    
    line([0.80*BoundaryX/2 0.80*BoundaryX/2], [0 BoundaryY/3]);
    line([1.40*BoundaryX/2 1.40*BoundaryX/2], [0 BoundaryY/3]);
    line([0.80*BoundaryX/2 1.40*BoundaryX/2], [0 0]);
    line([0.80*BoundaryX/2 1.40*BoundaryX/2], [BoundaryY/3 BoundaryY/3]);
    
    
    vrms = sqrt ((Vx.^2)+(Vy.^2));
    Vx(inbox&(~middle)) = -Vx(inbox&(~middle));
    Vy(inbox&middle) = -Vy(inbox&middle);
    meanfreetime = (Vtherm * (C.q_0))/200;
    meanfreepath = mean(vrms) * meanfreetime;
    vrms = sqrt ((Vx.^2)+(Vy.^2));
    show_m = (sqrt (2) * (mean(vrms)^2) * Em )/ (C.kb);
    array_1 (1,i)=  show_m;
    
    [x_mesh, y_mesh] = meshgrid(0:(BoundaryX/10):BoundaryX, 0:(BoundaryY/10):BoundaryY);
    electron_mat = zeros(11, 11);
    temperture_mat = zeros(11, 11);
    numelec_t = 0;
    t_vel = 0;
    
    for j = 1:10
        efx_min = x_mesh(1, j);
        efx_max = x_mesh(1,j+1);
        for k = 1:10
            efy_min = y_mesh(k, 1);
            efy_max = y_mesh(k+1, 1);
            for m = 1:nElectrons
                if((Px(m) > efx_min) && (Px(m) < efx_max) && ((Py(m) > efy_min) && Py(m) < efy_max))
                    numelec_t = numelec_t + 1;
                    %numelec_t=numelec_t++;
                    electron_mat(j, k) = electron_mat(j, j) + 1;
                    t_vel = t_vel + sqrt((Vx(m) .^ 2) + (Vx(m) .^ 2));
                    temperture_mat(j, k) = ((sqrt(2)*(t_vel/numelec_t) ^ 2) * Em) / (C.kb);
                end
            end
            t_vel = 0;
            numelec_t = 0;
        end
    end
end
 
fprintf('The Mean Free Time is = %f\n', meanfreetime);
fprintf('The Mean Free Path is = %f\n', meanfreepath);
figure (1)
title(['Average temperature value = ' num2str(show_m)]);
xlabel('particle on x axis');
ylabel('particle on y axis');
figure (2)
plot(tmp_arary_1, array_1);
title('Temperature across Time');
ylabel('Temperature');
xlabel('Time');
hold on
figure(3); histogram(vrms, 15);
title('Histogram of Thermal Velocities');
xlabel('Xlable');
ylabel('Ylable');
figure(4); surf(electron_mat);
title('Density Mapping');
xlabel('Xlable');
ylabel('Ylable');
figure(5); surf(temperture_mat);
title('Temperature Mapping');
xlabel('Xlable');
ylabel('Ylable');
 

