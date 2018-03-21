%% Henon Hailes, analysis of stability with the time step deltaT = h;
clear 
close all

E = [ 1/12];
    
for energy = E
    figure;

     for maxStep = 10000:2000:50000
        Ans = [0 0];

       j = figure;
       for xinit = 0:0
         for yinit = -0.25:0.05:0.25
    
           py0 = 0;
           x0 = xinit
           y0 = yinit
           px0 = (2*(energy - py0^2/2 - 0.5*(x0^2 + y0^2 + 2*x0^2*y0 - (2/3)*y0^3)))^.5;
           E_initial = EnergyCalc(px0, py0, x0, y0);
           t0 =0; tf = 10000; 
           n = maxStep;
           h = (tf-t0)/n;
           [pxAns, pyAns, xAns, yAns, tAns, dx, dy] = RK4HHFunc(px0, py0, x0, y0, tf, t0, n);


           for k = 1:length(yAns)
              if(abs(xAns(k)) < 0.015)
                 Ans = [Ans; yAns(k), pyAns(k)]; 
              end
           end
           plot(Ans(:,1), Ans(:,2), '.','color',rand(1,3));
           title(strcat('phase map with time step (h) = ', num2str(h),' E= ', num2str(energy)));
           xlabel('y coordinate')
           ylabel('py (y-momentum)')
           hold on;
       
         end   
         end
         
       
         saveas(j, strcat('phase map with time step (h) = ', num2str(h),' E= ', num2str(energy),'.png'))
         j2 = figure;
         plot(xAns, yAns, '-.','color', rand(1,3));
         title(strcat('trajectory (x vs y) with time step (h) = ', num2str(h),' E= ', num2str(energy)))
         xlabel( 'x coordinate')
         ylabel('y coordinate')
         saveas(j2, strcat('trajectory with time step (h) = ', ...
              num2str(h),' E= ', num2str(energy),'.png'))
    end
%     figure;
%     s =plot(Ans(:,1), Ans(:,2), '.','color',rand(1,3));
%   
%     dim = [.2 .5 .3 .3];
%     title(strcat('Phase Map at Energy = ', num2str(energy)));
%     str = strcat('px0=',num2str(px0),' x0=',num2str(x0),' y0=', num2str(y0));
%     figure;
%     scatter3(xAns, yAns, pyAns);
    
end