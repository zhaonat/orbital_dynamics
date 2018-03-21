%% Runge Kutta Function

function [pxAns, pyAns, xAns, yAns, tAns, l, n] = RK4HHFunc(px0, py0, x0, y0, tf, t0, N)
    h = (tf-t0)/N;
    %% specification of f(x,y,t) and g(x,y,t)
    f = @(px, py, x, y, t) (-1*(x + 2*x*y));
    %f = @(px, py, x, y, t) -x; %dpx/dt
    %g = symfun(px, [px py x y t]); %dx/dt
    g = @(px, py, x, y, t) px;
    a =  @(px, py, x, y, t) (-1*(y + x^2 - y^2));
    %a = @(px, py, x,y, t) -y; %dpy/dt
    %b = symfun(py, [px py x y t]); %dy/dt
    b = @(px, py, x, y, t) py;

    y_t = zeros(N,1);
    x_t = zeros(N,1);
    px_t = zeros(N,1);
    py_t = zeros(N,1);

    %% Setup Initial Conditions
%     x0 = .1; px0 = 0;
%     y0 = .1; py0 = 0;

    y_t(1) = y0;
    x_t(1) = x0;
    py_t(1) = py0;
    px_t(1) = px0;

    %% Calculate Energy Associated with the Initial Position

    %% Need Four Sets of Coefficients?

    result = [px0, x0, py0, y0, t0];

    for i = 1:N-1    
        %% updates for px, x 
        k1 = h*double(f(px0, py0, x0, y0, t0)); %dpx/dt
        l1 = h*double(g(px0, py0, x0, y0, t0)); %dx/dt
        m1 = h*double(a(px0, py0, x0, y0, t0)); %dpy/dt
        n1 = h*double(b(px0, py0, x0, y0, t0));  %dy/dt

        k2 = h*double(f(px0 + k1/2, py0 + m1/2, x0 + l1/2, y0 + n1/2, t0 + h/2));
        l2 = h*double(g(px0 + k1/2, py0 + m1/2, x0 + l1/2, y0 + n1/2, t0 + h/2));
        m2 = h*double(a(px0 + k1/2, py0 + m1/2, x0 + l1/2, y0 + n1/2, t0 + h/2));
        n2 = h*double(b(px0 + k1/2, py0 + m1/2, x0 + l1/2, y0 + n1/2, t0 + h/2));

        k3 = h*double(f(px0 + k2/2, py0 + m2/2, x0 + l2/2, y0 + n2/2, t0 + h/2));
        l3 = h*double(g(px0 + k2/2, py0 + m2/2, x0 + l2/2, y0 + n2/2, t0 + h/2));
        m3 = h*double(a(px0 + k2/2, py0 + m2/2, x0 + l2/2, y0 + n2/2, t0 + h/2));
        n3 = h*double(b(px0 + k2/2, py0 + m2/2, x0 + l2/2, y0 + n2/2, t0 + h/2));

        k4 = h*double(f(px0 + k3, py0 + m3, x0 + l3, y0 + n3, t0 + h/2));
        l4 = h*double(g(px0 + k3, py0 + m3, x0 + l3, y0 + n3, t0 + h/2));
        m4 = h*double(a(px0 + k3, py0 + m3, x0 + l3, y0 + n3, t0 + h/2));
        n4 = h*double(b(px0 + k3, py0 + m3, x0 + l3, y0 + n3, t0 + h/2));


        k = (k1/6 + k2/3 + k3/3 + k4/6);
        l = (l1/6 + l2/3 + l3/3 + l4/6);
        m = (m1/6 + m2/3 + m3/3 + m4/6);
        n = (n1/6 + n2/3 + n3/3 + n4/6);

        %% update the positions and momentum
        x0 = x0 + l; y0 = y0 + n; t0 = t0+h;
        px0 = px0 + k; py0 = py0 + m;

        %% Update All Results Data Structurs
        x_t(i+1) = x_t(i) + l;
        y_t(i+1) = y_t(i) + n;
        px_t(i+1) = px_t(i) + k;
        py_t(i+1) = py_t(i) + m;
        arr = [px_t(i+1) x_t(i+1) py_t(i+1) y_t(i+1) t0+h];
        result = [result; arr];

    end
    pxAns = result(:,1);
    tAns = result(:,5);
    xAns = result(:,2);
    pyAns = result(:,3);
    yAns = result(:,4);


end