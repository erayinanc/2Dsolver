%% initialisation of fields
% different cases are defined here for scalars
function [ RhoU, RhoV, phi] = initialise( u, v, JJet, vprofile, pprofile )
    global Ima Jma nG Jmid ;

    % Boundary conditions
    dom = zeros(Ima+2*nG,Jma+2*nG);
    % velocities
    switch(vprofile)
        case 0 % linear flow
            RhoU = dom+u;RhoV = dom+v;
        case 1 % vertical flow
            RhoU = dom+u;
            for i=1:1:Ima+2*nG;
                RhoU(i,:) = u.*(0.5.*sin((i-1)/2/pi)+0.5);
            end
            RhoV = dom+v;
        case 2 % jet flow
            RhoU = dom.*0;
            RhoV = dom+v;
            Jmidp=Jmid+1;
            
            RhoU(1:2,:) = 0.0*u;
            RhoU(1:2,Jmid-floor(JJet/2)-1:Jmidp+ceil(JJet/2)+1) = 0.5*u;
            RhoU(1:2,Jmid-floor(JJet/2):Jmidp+ceil(JJet/2)) = 1.0*u;
        case 3 % jet flow-filled
            RhoU = dom.*0;
            RhoV = dom+v;
            Jmidp=Jmid+1;
            
            RhoU(1:2,:) = 0.0*u;
            RhoU(1:end,Jmid-floor(JJet/2)-1:Jmidp+ceil(JJet/2)+1) = 0.5*u;
            RhoU(1:end,Jmid-floor(JJet/2):Jmidp+ceil(JJet/2)) = 1.0*u;
        case 4 % vortical flow
            RhoU = dom.*0;
            RhoV = dom.*0;
            a = linspace(-u,u,Ima+2*nG);
            b = linspace(-v,v,Jma+2*nG);

            for i = 1:Ima+2*nG;
                RhoV(i,:) = a(i);
            end;
            RhoV = RhoV./max(max(RhoV));

            for j = 1:Jma+2*nG;
                RhoU(:,j) = -b(j);
            end;
            RhoU = RhoU./max(max(RhoU));
    end

    % scalars
    phi = dom;            % transported function
    switch(pprofile)
        case 0       
            phi = RhoU./max(max(RhoU)); 
        case 1     
            phi(round(Ima/2)-5:round(Ima/2)+5,round(Jma/2)-5:round(Jma/2)+5) = 1;  
    end
end

