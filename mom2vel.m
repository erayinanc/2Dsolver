%% mom2vel
%  Calculate the velocity on cell surfaces from the momentum at cell centers
function [fluxU,fluxV]=mom2vel(rhoU,rhoV,rho)
    % Initialisation
    [Ima,Jma]=size(rhoU);

    % Calculation
    fluxU = 0.5/rho * (rhoU(1:Ima-1,:) + rhoU(2:Ima,:));
    fluxV = 0.5/rho * (rhoV(:,1:Jma-1) + rhoV(:,2:Jma));
end