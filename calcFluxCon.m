%% Caculate the convective fluxes
%  fluxConX, fluxConY over the cell surfaces as a function of
%  a conserved scalar field Phi from the fluxes U and V 
%  with different schemes, DDS will never work!

function [fluxConX,fluxConY]=calcFluxCon(Phi,fluxU,fluxV,dt,deltaX,scheme)
    % Initialisation
    [Ima,Jma]=size(Phi);

    % Calculation
    switch(scheme)
        case 0 % CDS
            fluxConX = 0.5 .* fluxU(1:Ima-1,:).*(Phi(1:Ima-1,:)+Phi(2:Ima,:));
            fluxConY = 0.5 .* fluxV(:,1:Jma-1).*(Phi(:,1:Jma-1)+Phi(:,2:Jma));
        case 1 % UDS
            fluxConX = max(0,fluxU(1:Ima-1,:)).*(Phi(1:Ima-1,:)) + ...
                       min(0,fluxU(1:Ima-1,:)).*(Phi(2:Ima,:));
            fluxConY = max(0,fluxV(:,1:Jma-1)).*(Phi(:,1:Jma-1)) + ...
                       min(0,fluxV(:,1:Jma-1)).*(Phi(:,2:Jma));
        case 2 % DDS
            fluxConX = fluxU(1:Ima-1,:).*(Phi(2:Ima,:));
            fluxConY = fluxV(:,1:Jma-1).*(Phi(:,2:Jma));
        case 3 % UDCDS         
            % Calculate weigths
            ww(1:Ima-1,:) = 0.5.*(fluxU(1:Ima-1,:)*dt/deltaX + 1);
            we(2:Ima,:)   = 1.0-ww(1:Ima-1,:);
            ws(:,1:Jma-1) = 0.5.*(fluxV(:,1:Jma-1)*dt/deltaX + 1);
            wn(:,2:Jma)   = 1.0-ws(:,1:Jma-1);

            % Calculation
            fluxConX = fluxU(1:Ima-1,:)...
                                .* ( ww(1:Ima-1,:).*Phi(1:Ima-1,:)...
                                   + we(2:Ima,:)  .*Phi(2:Ima,:) );
            fluxConY = fluxV(:,1:Jma-1)...
                                .* ( ws(:,1:Jma-1).*Phi(:,1:Jma-1)...
                                   + wn(:,2:Jma)  .*Phi(:,2:Jma) );
    end
end