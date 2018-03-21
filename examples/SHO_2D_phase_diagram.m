%% Runge Kutta Benchmark

px0 = 1; py0 =1;
x0 = 1; y0 = 1;
tf = 1000; t0 = 0; N = 4000;
[pxAns, pyAns, xAns, yAns, tAns, l,n] = RK4HHSimpleFunc(px0, py0, x0, y0, tf, t0, N); 

figure;
plot(pyAns, yAns);
xlabel('y')
ylabel('py')
title('Phase Map for SHO in y Coordinate')
figure;
plot(pxAns, xAns);
xlabel('x')
ylabel('px')
title('phase Map for SHO in x Coordinate')