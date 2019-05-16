% B-spline curve estimation (control points estimation) from data about points on curve.
% Author : rudraa (Implemented using the bspline library by Levente Hunyadi)
% https://www.mathworks.com/matlabcentral/fileexchange/27374-b-splines
%
clear all; 
% spline order
k = 3; %Order in IGa notation is one less than the notation here. So k=3 is a 2nd order curve.
numKnots=40;
t = [0 0 linspace(0,1,numKnots) 1 1 ]; %k repetitions of 0 and 1 at the ends of the knot vector
%
R=20.0;
r=1*R;
%
%populate M (2XN points on the curve... where N is some large number that is sufficient for estimating the control points)
%tube
%theta=linspace(4*R,r,4*numKnots); 

%base line
%theta=linspace(17*R,R+r,ceil(0.5*numKnots)); 
%M=[theta; zeros(size(theta))];

%base arc 
theta=linspace(pi/2,0,ceil(2*numKnots)); theta=theta(2:end);
M=[R+r-r*cos(theta);r-r*sin(theta)];
%M=[M(1,:) R+r-r*cos(theta);M(2,:) 0.5*(r-r*sin(theta))];

%line
theta=linspace(r,10*R,6*numKnots); 
%theta=linspace(0,10*R,2*numKnots); 
M=[M(1,:) R*ones(size(theta));M(2,:) theta];
%M=[R*ones(size(theta));theta];

%cap 
%theta=linspace(0, 0.8*pi/2,ceil(numKnots)); theta=theta(2:end);
%M=[M(1,:) R*cos(theta); M(2,:) 10*R+R*sin(theta)];

%shift
M(2,:)=M(2,:)-min(M(2,:));

%line
%theta=linspace(R+r,1.5*R,ceil(numKnots)); theta=theta(2:end);
%M=[M(1,:) theta; M(2,:) R*zeros(size(theta))];

D = bspline_estimate(k,t,M);
C = bspline_deboor(k,t,D);

% plot control points and spline
figure;
hold all;
plot(M(1,:), M(2,:), 'k');
plot(D(1,:), D(2,:), 'rx');
plot(C(1,:), C(2,:), 'c');
legend('data', 'estimated control points', 'estimated curve', ...
    'Location', 'Best');
hold off;
axis equal;
%save to file
order=k-1;
knots=t;
controlPoints=D;
save(['tubeFullWithBaseFor1R4Helix' num2str(numKnots) '.mat'],'order','knots','controlPoints','-v6') %version 6 format needed to read into python using scipy.io.loadmat
