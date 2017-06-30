%% Caculate the Diffusive fluxes
%  fluxDifX, fluxDifY over the cell surfaces as a function of
%  a conserved scalar field Phi, the diffusivity D and a 
%  grid-spacing of deltaX

function [fluxDifX,fluxDifY]=calcFluxDif(Phi,dx,D)
    % Initialisation
    [Ima,Jma]=size(Phi);
    % Calculation
    fluxDifX = D/dx.*(Phi(2:Ima,:) - Phi(1:Ima-1,:));
    fluxDifY = D/dx.*(Phi(:,2:Jma) - Phi(:,1:Jma-1));
end