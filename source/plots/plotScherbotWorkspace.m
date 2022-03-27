%%%%%%%%%%%%%%%%%%%%%%%%%%%%% parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1000; % point to draw circles
wsColor = [0.8 0.8 0.8]; % grey color for workspace
invSingWidth = 5;
L = model.L;
Rleft = L(2)+L(3);
Rright = L(4)+L(5);
Cleft = [0;0];
Cright = [L(1);0];
Rinleft = L(3)-L(2);
Rinright = L(5)-L(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% circles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parametric circle
theta = linspace(0,2*pi,N);
% left inner circle
Xcileft = Cleft(1)+Rinleft*cos(theta);
Ycileft = Cleft(2)+Rinleft*sin(theta);
% right inner circle
Xciright = Cright(1)+Rinleft*cos(theta);
Yciright = Cright(2)+Rinleft*sin(theta);

% intersection points of outer circles
R = norm(Cleft-Cright); % distance between circles
temp = 0.5*[Cleft(1)+Cright(1); Cleft(2)+Cright(2)] + ...
              (Rleft^2-Rright^2)/(2*R^2)*[Cright(1)-Cleft(1); Cright(2)-Cleft(2)];
sol1 = temp+0.5*sqrt(2*(Rright^2+Rleft^2)/R^2- ...
       (Rleft^2-Rright^2)^2/R^4-1)*[Cright(2)-Cleft(2); Cleft(1)-Cright(1)];
sol2 = temp-0.5*sqrt(2*(Rright^2+Rleft^2)/R^2- ...
       (Rleft^2-Rright^2)^2/R^4-1)*[Cright(2)-Cleft(2); Cleft(1)-Cright(1)];

% find parametric angle
theta_s = atan(sol2(2)/sol2(1));
% left angle
theta = linspace(-theta_s,theta_s,N);
% left outer circle
Xcleft = Cleft(1)+Rleft*cos(theta);
Ycleft = Cleft(2)+Rleft*sin(theta);
% right cangle
theta = linspace(pi-theta_s,pi+theta_s,N);
Xcright = Cright(1)+Rright*cos(theta);
Ycright = Cright(2)+Rright*sin(theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% drawing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% draw and fill outer circle
plot(ax,Xcleft,Ycleft,'b','LineWidth',invSingWidth)
htemp = fill(ax,Xcleft,Ycleft, wsColor);
set(htemp,'EdgeColor','none')
plot(ax,Xcright,Ycright,'b','LineWidth',invSingWidth)
htemp = fill(ax,Xcright,Ycright, wsColor);
set(htemp,'EdgeColor','none')
% draw and fill
plot(ax,Xcileft,Ycileft,'b','LineWidth',invSingWidth)
plot(ax,Xciright,Yciright,'b','LineWidth',invSingWidth)
fill(ax,Xcileft,Ycileft, 'w')
fill(ax,Xciright,Yciright, 'w')
