%% Runge Kutta total phase diagram creation
clear
close all

E = [1/12 1/8];
    
for energy = E
    fig1 = figure;
       Ans = [0 0];
%     for xinit = -0.1:0.02:0.1
%     for yinit = -.4:0.01:.4
     for xinit = 0:0
     for yinit = -0.25:0.05:0.25
       AnsTemp = [0 0];
       py0 = 0;
       x0 = xinit
       y0 = yinit
       px0 = (2*(energy - py0^2/2 - 0.5*(x0^2 + y0^2 + 2*x0^2*y0 - (2/3)*y0^3)))^.5;
       E_initial = EnergyCalc(px0, py0, x0, y0);
       t0 =0; tf = 10000; 
       n = 100000;
       [pxAns, pyAns, xAns, yAns, tAns, dx, dy] = RK4HHFunc(px0, py0, x0, y0, tf, t0, n);

      
       for k = 1:length(yAns)
          if(abs(xAns(k)) < 0.01)
             Ans = [Ans; yAns(k), pyAns(k)]; 
             AnsTemp = [AnsTemp; yAns(k), pyAns(k)];
          end
       end
       plot(AnsTemp(:,1), AnsTemp(:,2),'.');
       hold on;
    end
    end
    fig2 = figure;
    s = plot(Ans(:,1), Ans(:,2), '.','color',rand(1,3));
 
    dim = [.2 .5 .3 .3];
    title(strcat('Phase Map at Energy = ', num2str(energy)));
    str = strcat('px0=',num2str(px0),' x0=',num2str(x0),' y0=', num2str(y0));
    annotation('textbox',dim,'String',str,'FitBoxToText','on');
    
    fig3 = figure;
    scatter3(xAns, yAns, pyAns);
    saveas(fig1, strcat('E=',num2str(energy),' phaseDiagramDifferentInitialPoints.png'));
    saveas(fig2, strcat('E=',num2str(energy),'phaseDiagramTotal.png'));
    saveas(fig3, strcat('E=',num2str(energy),'pyAsfunctionofxandy.png'));
    
    
    
end