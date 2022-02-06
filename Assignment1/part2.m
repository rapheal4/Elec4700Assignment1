
clearvars
clearvars -GLOBAL
close all
set(0,'DefaultFigureWindowStyle', 'docked')


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
BoundaryX = 200e-9;                    % X boundary
BoundaryY = 100e-9;                    % Y boundary
T = 300;                            % Semiconductor temperature
            
t_mn = 0.2e-12;                     % Mean time between collisions
                    
TimeSteps = 1000;                   % Number of time steps

nElectrons = 1000;                   % Number of electrons

dt = 1e-14;                         % Time Step



PPx(1: nElectrons) = rand(nElectrons, 1) * BoundaryX;
PPy(1: nElectrons) = rand(nElectrons, 1) * BoundaryY;

Vtherm = sqrt(2 * C.kb * T/Em);

Vx(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
Vy(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));

myColors = ['r' 'b' 'g' 'y' 'm' ];
myColorTyp = 1;

Pscat = 1 - exp(-(dt/t_mn));

a = randi(nElectrons,5,1);

TAvgp = 300;

for i=2:TimeSteps
    
   if(Pscat > rand())
      Vx(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
      Vy(1: nElectrons) = randn(nElectrons, 1) * (Vtherm/sqrt(2));
   end
    
   Px(1: nElectrons) = PPx(1: nElectrons) + (Vx .* dt);
   Py(1: nElectrons) = PPy(1: nElectrons) + (Vy .* dt);

   Vy((Py>BoundaryY) | (Py<0)) = -Vy((Py>BoundaryY) | (Py<0));
   
   for j=1:5
       subplot(2,1,1);
       plot([PPx(a(j)) Px(a(j))], [PPy(a(j)) Py(a(j))],myColors(j));
       xlim([0 BoundaryX]);
       ylim([0 BoundaryY]);
   end
   pause(0.1)
   hold on
   title('2-D plot of particle trajectories');
   
   VxAbs = abs(Vx);
   VyAbs = abs(Vy);
   
   TAvg = (mean((VxAbs.^2)+ (VyAbs.^2)) * Em)/(2 * C.kb);
 
   
   subplot(2,1,2);
   plot([i-1 i],[TAvgp TAvg],'r');
   xlim([0 TimeSteps]);
   ylim([0 800]);
   pause(0.1)
   hold on
   title('Average Temperature');
   
   Px(Px>BoundaryX) = Px(Px>BoundaryX)-BoundaryX;
   Px(Px<0) = BoundaryX;
   
   PPx = Px;
   PPy = Py;
   

%    
   
end


vtherm = sqrt(2 * C.kb * T/Em);
v_change = dt/vtherm;
t_mn = 0.2e-12;
 nElectrons =1000; 
X = rand(2,nElectrons);
Y = rand(1,nElectrons);
positionX(1,:) = X(1,:)*BoundaryX;
positionY(1,:) = Y(1,:)*BoundaryY;
 
left = positionX > 0.9e-7;
right = positionX < 1.1e-7;
check_X = left & right;
top = positionY > 0.66e-7;
bottom = positionY < 0.5e-7;
top = top & check_X;
bottom = bottom & check_X;
BOX = top | bottom;
 
while(sum(BOX) > 0)
    
    x1 = rand(1,sum(BOX));
    y1 = rand(1,sum(BOX));
    positionX(BOX) = x1*BoundaryX;
    positionY(BOX) = y1*BoundaryY;
    
    left = positionX > 0.9e-7;
    right = positionX < 1.1e-7;
    check_X = left & right;
    top = positionY > 0.66e-7;
    bottom = positionY < 0.5e-7;
    top = top & check_X;
    bottom = bottom & check_X;
    BOX = top | bottom;
end

angle(1,:) = X(2,:)*2*pi;
sigma = sqrt(2 * C.kb * T/Em)/4;
max_boltz_disttribution = makedist('Normal',vtherm,sigma);
velocity = random(max_boltz_disttribution,1,nElectrons);
figure(6)
hist(velocity)
title('Particle Velocity Histogram')
X_velocity = v_change*velocity(1,:).*cos(angle(1,:));
Y_velocity = v_change*velocity(1,:).*sin(angle(1,:));
PSCAT = 1 - exp(-v_change/t_mn);
mfp_vec = zeros(1,nElectrons);
