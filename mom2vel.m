%% mom2vel
%  Calculate the velocity on cell surfaces from the momentum at cell centers
function [fluxU,fluxV]=mom2vel(rhoU,rhoV,rho)
    global Ifim Ifi Ila Ilap Jfim Jfi Jla Jlap;

    % Calculation
    fluxU = 0.5/rho * (rhoU(Ifim:Ila,:) + rhoU(Ifi:Ilap,:));
    fluxV = 0.5/rho * (rhoV(:,Jfim:Jla) + rhoV(:,Jfi:Jlap));
end
