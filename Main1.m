clear all
close all 

dx=0.0025*4; % input('Please enter delta x: ')
T_Heat=100; % input('Please input the boundary heating conditions: ');
T_Cool=0; % input('Please input the boundary cooling conditions: ');
tol = 0.01; % input('Please enter the temperature convergence criterion: ');
xH = 0.3; % X Position of Centre of Hole
yH = 0.3; % Y Position of Centre of Hole
tol = 0.01;
optsteps = 4;

%Create dx and dy steps 
dxstep = linspace(dx*optsteps, dx, optsteps);
dystep = dxstep; 

phiint = []; 
Xinit = [];
Yinit = []; 

for i = 1:optsteps;
    [phiconv, X, Y, xc, yc, errsave] = Heat_Transfer_Problem_Updated_4(dxstep(i),dystep(i), tol, xH,yH, T_Heat, T_Cool, phiint, Xinit, Yinit, optsteps, i); 
    phiint = phiconv;
    Xinit = X; 
    Yinit = Y;
    i;
end 


%% Plot Result

%Surface Contour Plot
P = figure;
surf(X,Y,flipud(phiconv)); %flipud flips the matrix across the horizontal: needs to be done due to the way the matrix is defined.
xlabel('X Position / m');
ylabel('Y Position / m');
title('Temperature Distribution / °C');

%Surface Contour Plot
S = figure;
surf(X,Y,flipud(phiconv));
view(2); %See from above
% colorbar
axis equal
xlabel('X Position / m');
ylabel('Y Position / m');
title('Temperature Distribution / °C (Top View');
zmax = max(max(phiconv));
line(xc, yc, zmax*ones(size(xc,2),1), 'Color', 'k');

%Plot Error Convergence 
C = figure;
plot(1:size(errsave,2), errsave, '-b');
title('Change in Temperature vs. Iteration Number'); 
xlabel('Iteration Number');
ylabel('Error Change');