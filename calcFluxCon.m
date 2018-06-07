%% Caculate the convective fluxes
%  fluxConX, fluxConY over the cell surfaces as a function of
%  a conserved scalar field Phi from the fluxes U and V 
%  with different schemes, DDS will never work!
function [fluxConX,fluxConY]=calcFluxCon(Phi,fluxU,fluxV,dt,scheme)
    global dx Ifim Ifi Ila Ilap Jfim Jfi Jla Jlap;

    % Calculation
    switch(scheme)
        case 0 % CDS
            fluxConX = 0.5 .* fluxU(Ifim:Ila,:).*(Phi(Ifim:Ila,:)+Phi(Ifi:Ilap,:));
            fluxConY = 0.5 .* fluxV(:,Jfim:Jla).*(Phi(:,Jfim:Jla)+Phi(:,Jfi:Jlap));
        case 1 % UDS
            fluxConX = max(0,fluxU(Ifim:Ila,:)).*(Phi(Ifim:Ila,:)) + ...
                       min(0,fluxU(Ifim:Ila,:)).*(Phi(Ifi:Ilap,:));
            fluxConY = max(0,fluxV(:,Jfim:Jla)).*(Phi(:,Jfim:Jla)) + ...
                       min(0,fluxV(:,Jfim:Jla)).*(Phi(:,Jfi:Jlap));
        case 2 % UDCDS         
            % Calculate weigths
            ww(Ifim:Ila,:) = 0.5.*(fluxU(Ifim:Ila,:)*dt/dx + 1);
            we(Ifi:Ilap,:) = 1.0-ww(Ifim:Ila,:);
            ws(:,Jfim:Jla) = 0.5.*(fluxV(:,Jfim:Jla)*dt/dx + 1);
            wn(:,Jfi:Jlap) = 1.0-ws(:,Jfim:Jla);

            % Calculation
            fluxConX = fluxU(Ifim:Ila,:)...
                                .* ( ww(Ifim:Ila,:).*Phi(Ifim:Ila,:)...
                                   + we(Ifi:Ilap,:).*Phi(Ifi:Ilap,:) );
            fluxConY = fluxV(:,Jfim:Jla)...
                                .* ( ws(:,Jfim:Jla).*Phi(:,Jfim:Jla)...
                                   + wn(:,Jfi:Jlap).*Phi(:,Jfi:Jlap) );
    end
end
