function [phiconv, X, Y, xc, yc, errsave] = HeatTransferProblem(dx,dy,tol,xH,yH,T_Heat, T_Cool, phiold, Xinit, Yinit, optsteps, currentloop)
%% Define Variables
length = 0.8;
width = 0.4;

nx=ceil(length/dx);
ny=ceil(width/dy);
 
if currentloop~=optsteps    
    if ceil(nx)~=nx || ceil(ny)~=ny
       disp('The dx or dy values provided are not exactly divisible by the lengths or width of the plate.');
    end
end

%Plate Dimensions
length = 0.8;
width = 0.4;
a=0.1; % Horizontal length of plate being T_Heated or cooled
b = 0.2; % Vertical length of plate being heated or cooled
R = 0.05; %Radius of Hole

%% Define Mesh 

%Calculate Location of Mesh Grid Points
[X,Y] = meshgrid(linspace(0, length, nx) , linspace(0, width, ny)); 

%Define phiold by interpolating between old values and new values
if ~isempty(Xinit)
    phiold = griddata(Xinit, Yinit, phiold, X, Y);
end

%Calculate Number of Nodes to represent Boundary Conditions 
apos=ceil((a/length)*nx); % Number of Mesh Nodes for Horizontal Length of of plate being Heated or Cooled
bpos=ceil((b/width)*ny); % Number of Mesh Nodes for Vertical Length of of plate being Heated or Cooled

%Find Nodes Representing the Circle
deg = pi/180 * (0:0.1:360); %Discretize Circle Angles
xc = xH + R*cos(deg);
yc = yH + R*sin(deg);

Mat_c = flipud(inpolygon(X, Y, xc, yc)); % Uses xc and yc to create a polygon to represent a circle, and evaluates whether the coordinates X,Y are in the polygon. Returns a matrix of 0s or 1s based on if the coordinate is in the polygon.[Row_c, Col_c] = find(Mat_c);

% Find the coordinates on the edge of the circle
Mat_c_o=zeros(ny,nx);
for i=1:nx;
    for j=1:ny;
        if Mat_c(j,i)==1;
            if Mat_c(j,i+1)==0 || Mat_c(j,i-1)==0 || Mat_c(j+1,i)==0 || Mat_c(j-1,i)==0;
               Mat_c_o(j,i) = 1;
            end
        end
    end
end

No_Coord_c_o = sum(sum(Mat_c_o)); %Number of Coordinates on the edge of the circle
No_Coord_c = sum(sum(Mat_c)); %Number of Coordinates in Circle

xrange = 2:ny-1;
yrange = 2:nx-1;

%% Define Storage Vectors and Variables

%Only for the first Optimization Loop 
if isempty(phiold)
    phiold = zeros(ny, nx);
end

err=tol+1;
errsave = zeros(1,10000); %Create vector for saving Delta T, set length so vector length doesn't change every loop
count = 0; %Counter for Iteration Number

%% Set Initial Boundary Conditions
% Initial boundary conditions for rectangular plate
%Dirichlet Boundary Conditons
phiold((1:bpos),1) = T_Heat; %Left Side
phiold(1,(1:apos)) = T_Heat; %Top Side
phiold((ny-bpos+1):ny,nx) = T_Cool; %Right Side
phiold(ny,nx-apos+1:nx) = T_Cool; %Bottom Side

phiold = SetBoundaryConditions(phiold, bpos, apos, ny, nx, Mat_c_o, No_Coord_c_o, Mat_c); %Set Insulation and Circle Boundary Conditions Just in Case
phinew = phiold; %Done because the Dirichelet Boundary Conditions only need to be set once

%% Set Up Loop
while err >tol
    
    count = count + 1; %Increment Loop Counter    
    
    %For every position other than the boundary, calculate phinew 
    phinew(xrange, yrange) = 0.25*(phiold(xrange+1, yrange) + phiold(xrange-1, yrange) + phiold(xrange, yrange+1) + phiold(xrange, yrange-1));
  
    % Set Boundary Conditions
    phinew = SetBoundaryConditions(phinew, bpos, apos,ny,nx,Mat_c_o,No_Coord_c_o,Mat_c);

    %Calculate error by finding difference between phinew and phiold
    err = sum(sum(abs(phinew(2:end-1, 2:end-1)-phiold(2:end-1, 2:end-1))));
    errsave(count) = err; %Save the error
    
    %Swap old and new for next iteration after calculating error
    phiold = phinew;
   
end

%Reduce Errsave Vector to Used Values
errsave = errsave(1:count);

phiconv = phinew; %Save final value of phinew to Output
end

%% Nested Function: Boundary Conditions

function phinew = SetBoundaryConditions(phinew, bpos, apos, ny, nx, Mat_c_o, No_Coord_c_o, Mat_c)

%Insulated Boundary Conditions
    phinew((bpos+1:ny),1) = phinew((bpos+1:ny),2); %Left Side
    phinew(1,(apos+1:nx))= phinew(2,(apos+1:nx)); %Top Side
    phinew(1:ny-bpos,nx) = phinew(1:ny-bpos,nx-1);   %Right Side
    phinew(ny,1:nx-apos) = phinew(ny-1,1:nx-apos);%Bottom Side
    
%Boundary conditions for the circle
    Avetemp=(sum(sum((phinew).*(Mat_c_o))))/No_Coord_c_o; %Find the average temperature of the edge values by summing and using matrix operations
    phinew=(phinew.*(~Mat_c))+(Avetemp*Mat_c); %Use matrix operations to input this average temperature value into phinew
end