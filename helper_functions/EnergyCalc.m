%% Energy Calculator

function Energy = EnergyCalc(px, py, x, y)
    KE = 0.5*(px^2 + py^2);
    PE = 0.5*(x^2 + y^2 + 2*x^2*y - (2/3)*y^3);
    Energy = KE+PE;
end