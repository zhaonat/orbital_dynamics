% RUNGE KUTTA 1D Example EOM

%1D example with a projectile launched with some initial velocity 
% in x, y coordinate space
% verifies 1D Runge Kutta implementation is correcty

close all
clear


syms x y v t %position_x, position_y, initial_velocity, time
t0 = 0;

tf = 2.0;
n = 20; 
v0 = 5;
h = (tf-t0)/n;
m = 1; g= 9.8;
%% specification of f(x,y,t) and g(x,y,t)
f = symfun(-v^2, [v t])
g = symfun(-(g - 1/2*v^2) , [v t]);

v_t = zeros(n,1);
v_t(1) = v0;
result = [v0, t0];

for i = 1:n-1
    
    k1 = double(h*f(v0, t0));
    k2 = double(h*f(v0 + k1/2, t0+h/2));
    k3 = h*double(f(v0 + k2/2, t0+h/2));
    k4 = h*double(f(v0 + k3, t0+h/2));
    
    k = (k1/6 + k2/3 + k3/3 + k4/6);
    v0 = v0+k;
    v_t(i+1) = v_t(i) + k;
    result = [result; v_t(i+1), t0+h];
    t0 = t0+h;
    
end

%% visualize solution
figure;
plot(result(:,2), result(:,1), '-o')
hold on;
%plot(result(:,2), 1/(1/5+(result(:,2))), 'x');

%% Symbolic Solution analytic solution
syms a y(t)
cond = y(0) == 5;
eqn = diff(y,t) == -y^2;
soln = dsolve(eqn, cond);
fplot(soln)
xlim([0 max(result(2,:))])
ylim([0 5])