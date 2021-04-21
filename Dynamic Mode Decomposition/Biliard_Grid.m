
%Biliard Grid
%Lama Hamadeh

close all;
clear all;

% % Notify me that it's beginning the simulation
%fprintf('\n Simulating...\n\n');
% %-----------------------------

%Construct the flow map
%The flow map bascially takes the initial pair (s0,p0) and evolve it so
%it hits the boundary/edge where a reflection/bounce happens there 
%with a specific direction which produces a new pair (s1,p1).
%Once the first bounces is set, constructing the iterations of the 
%subsequent bounces would be straightforward. 

%Step 1: Construct the xy flow map mesh

%boundary length 
L = 1; 
%number of points on S axis
ns = 25; 
%boundary variable
s = linspace(0,L,ns); 
%construct coordinates meshgrid
[X,Y] = meshgrid(s,s); 
% mesh needs X,Y and Z so create z
Z = zeros(size(X));
%Visualise the grid
figure;
mesh(X,Y,Z,'Marker','o','EdgeColor',"k") %or surf
axis equal tight
view(2)
set(gca,'ytick',[])
xlabel('$S$','Interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
set(gca,'FontSize',16)
hold on
% %-----------------------------

%iteration (this is "Step 5" after constructing the initial bounce)
%for k = 0:10 %iterations
%for s = 0:L/ns:L %starting points
%for p = -1:1 %starting directions
%Step 2: Compute both line equations in form of Ax + By + C = 0.

%Since the slope-intercept line formula does not hold in case of a vertical
%line, i.e., x2-x1 = 0 and hence the gradient m=y2-y1/x2-x1 = NaN, then the
%more generic approach that should be capable of defining every line is to
%write the line equation in its standard form as: Ax + By + C = 0, where 
%here a vertical line would simply mean B = 0.

for iter = 1 : 5

%Compute first line equation | L1: Ax + By + C = 0 
%first line coordinates
x1L1 = 0.7;  y1L1 = 0;
x2L1 = 1  ;  y2L1 = 0.5;

%Setting the determinant to be zero to get the line's standrad formula
A = y1L1 - y2L1;
B = x2L1 - x1L1;
C = x1L1*y2L1 - x2L1*y1L1; %or C=-A*x1-B*y1
%L1 standard line equation
L1 = (-C/B)+(-A/B)*s; % -A/B is the gradient and -C/B is the intercept
%plotting on the grid
plot(s,L1,'b','LineWidth',2)
axis([min(s) max(s) min(s) max(s)])
%find the angle between L1 and the norm N
L1v = [x2L1,y2L1] - [x1L1,y1L1]; %compute the x and y projections of L1
%defining the vertical norm
vnorm = x1L1;
Nvx = [0 vnorm];  %the verticle normal vector 
%plot the verticle norm line
y=get(gca,'ylim');
hold on
plot([vnorm vnorm],y,'g','LineWidth',2)
%dot product between L1v and Nvx to get the angle
theta1 = acos((L1v(1,1) * Nvx(1,1)+L1v(1,2) * Nvx(1,2))/...
    (sqrt((L1v(1,1))^2+(L1v(1,2))^2)*sqrt((Nvx(1,1))^2+(Nvx(1,2))^2))); 
%define the angular variable
p1 = sin(theta1);
%show results on the screen
fprintf('The angle between the first line and the vertical norm: theta1=%g, and p1=%g\n',...
    theta1*(180/pi),p1');

%constructing the first bounce, i.e., evolving L1 so it hits the 
%"next" boundary and reflects to produce the first bounce    

%Compute second line equation | L2: Dx + Ey + F = 0 
%second line coordinates
x1L2 = x2L1;    y1L2 = y2L1;

if x1L2 == 1
    x2L2 = randi(length(s));     x2L2 = s(x2L2); %randomly generated
    y2L2 = 1;
    if x2L2-x1L2 == 0 %avoiding zero gradient
        x2L2 = randi(length(s));     x2L2 = s(x2L2);
    end
end

if y1L2 == 1
    y2L2 = randi(length(s));     y2L2 = s(y2L2); %randomly generated
    x2L2 = 0;
    if y2L2-y1L2 == 0 %avoiding undefined gradient 
        y2L2 = randi(length(s));     y2L2 = s(y2L2);
    end
end

if x1L2 == 0
    x2L2 = randi(length(s));     x2L2 = s(x2L2); %randomly generated
    y2L2 = 0;
    if x2L2-x1L2 == 0 %avoiding zero gradient
        x2L2 = randi(length(s));     x2L2 = s(x2L2);
    end
end

if y1L2 == 0
    y2L2 = randi(length(s));     y2L2 = s(y2L2); %randomly generated
    x2L2 = 1;
    if y2L2-y1L2 == 0 %avoiding undefined gradient
        y2L2 = randi(length(s));     y2L2 = s(y2L2);
    end
end

%Setting the determinant to be zero to get the line's standrad formula
D = y1L2 - y2L2; 
E = x2L2 - x1L2;
F = x1L2*y2L2 - x2L2*y1L2; %or F=-D*x2-E*y2
%L2 stamdard line equation
L2 = (-F/E)+(-D/E)*s; % -D/E is the gradient and -F/E is the intercept
%plotting on the grid
plot(s,L2,'r','LineWidth',2)
axis([0 1 0 1])
%find the angle between L2 and the norm N
L2v = [x2L2,y2L2] - [x1L2,y1L2]; %compute the x and y projections of L2

%if the gradient positive, compute the angle with the vertical norm
if -D/E > 0
    %defining the vertical norm
    vnorm = x1L2;
    Nvx = [0 vnorm];  %the verticle normal vector 
    %plot the verticle norm line
    y=get(gca,'ylim');
    hold on
    plot([vnorm vnorm],y,'g','LineWidth',2)
    %dot product between L1v and Nvx to get the angle
    theta2 = acos((L2v(1,1) * Nvx(1,1) + L2v(1,2) * Nvx(1,2))/...
        (sqrt((L2v(1,1))^2+(L2v(1,2))^2)*sqrt((Nvx(1,1))^2+(Nvx(1,2))^2))); 
    %define the angular variable
    p2 = sin(theta2);
    %show results on the screen
    fprintf('The angle between the first line and the vertical norm: theta2=%g, and p2=%g\n',...
            theta2*(180/pi),p2');
end

%if the gradient negative, compute the angle with the horizontal norm
if -D/E < 0
    %defining the horizontal norm
    hnorm = y1L2;
    Nhy = [hnorm 0]; %the normal vector
    %plot the horizontal norm line
    x=get(gca,'xlim');
    hold on
    plot(x,[hnorm hnorm],'g','LineWidth',2)
    %dot product between L2 and N to get the angle
    theta2 = acos((L2v(1,1) * Nhy(1,1) + L2v(1,2) * Nhy(1,2))/...
        (sqrt((L2v(1,1))^2+(L2v(1,2))^2)*sqrt((Nhy(1,1))^2+(Nhy(1,2))^2))); 
    theta2 = pi-theta2;
    %define the angular variable
    p2 = sin(theta2);
    %show results on the screen
    fprintf('The angle between the second line and the horizontal norm: theta2=%g, and p2=%g\n',...
            theta2*(180/pi),p2');
end

pause(2)

end


% %-------

%Step 3: Before finding the intersection point coordinate, check whether 
%the lines are parallel or not by ensuring if determinant is zero 
%lines are parallel.

Coeff_m = [A B ; D E]; %matrix of coefficients (A in Ax=b)
determinant = det(Coeff_m);

if determinant ==0
    fprintf('Both lines are parallel and so no unique point where they meet exists!');

    % %-------

    %Step 4: Find the coordinates of the intersection point (solving Ax=b)

    else
        Coeff_v = [-C ; -F]; %vector of coefficients (b in Ax=b)
        inter_coor = Coeff_m\Coeff_v; %coordinates vector: x and y (x in Ax=b)
        fprintf('Both lines meet at the point (x=%g , y=%g)\n',...
            inter_coor(1,1),inter_coor(2,1)); %show the intersection point coordinates
end


% %-----------------------------

% %Visualise space-phase coordinates

% %Defining space variables
% L = 4; %boundary length 
% 
% ns = 25; %number of points on S axis
% np = 25; %number of points on P axis
% 
% s = linspace(0,L,ns); %boundary variable
% p = linspace(-1,1,np); % direction/angle variable
% 
% [S,P] = meshgrid(s,p); %construct coordinates meshgrid
% % mesh needs X,Y and Z so create z
% Z = zeros(size(S));
% 
% %Visualise the grid
% figure;
% mesh(S,P,Z,'Marker','o','MarkerFaceColor','k','EdgeColor',"k")
% axis equal tight
% view(2)
% xlabel('$S$','Interpreter','latex')
% ylabel('$P$','Interpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'FontSize',16)
% hold on
% %-----------------------------




