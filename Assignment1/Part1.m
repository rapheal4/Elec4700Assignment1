global C
global Em T
global BoundaryX BoundaryY
global PPx Px PPy Py Vx Vy
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

%set colour
myColors = ['r' 'b' 'g' 'y' 'm' ];
myColorTyp = 1;



a = randi(nElectrons,5,1);
%a = horzcat (nElectrons,5,1);

TAvgp1 = 300;
for j=2:TimeSteps
    
   Px(1: nElectrons) = PPx(1: nElectrons) + (Vx .* dt);
   Py(1: nElectrons) = PPy(1: nElectrons) + (Vy .* dt);

   Vy((Py>=BoundaryY) | (Py<0)) = -Vy((Py>BoundaryY) | (Py<0));
   
   for k=1:5
       subplot(2,1,1);
       plot([PPx(a(k)) Px(a(k))], [PPy(a(k)) Py(a(k))],myColors(k));
       %plot(repmat([PPx(a(j)) Px(a(j))], [PPy(a(j)) Py(a(j))],myColors(j))
      % x_pos = rand(1, nElectrons) .* BoundaryX;
       %y_pos = rand(1, nElectrons) .* BoundaryY; 
       xlim([0 BoundaryX]);
       ylim([0 BoundaryY]);
   end
   pause(0.1)
   hold on
   title('The 2-D plot of particle trajectories');
   
   VxAbs = abs(Vx);
   VyAbs = abs(Vy);
   
   TT1=(VxAbs.^2);
   TT2=(VyAbs.^2);
   TemperatureAvg = (mean(TT1+ TT2) * Em)/(2 * C.kb);
 
   subplot(2,1,2);
   plot([j-1 j],[TAvgp1 TemperatureAvg],'r');
   xlim([0 TimeSteps]);
   ylim([0 600]);
   pause(0.1)
   hold on
   title('The Average Temperature');
   
   Px(Px>BoundaryX) = Px(Px>BoundaryX)-BoundaryX;
   Px(Px<0) = BoundaryX;
   
   PPx = Px;
   PPy = Py;
   
   TAvgp1 = TemperatureAvg;
   
end