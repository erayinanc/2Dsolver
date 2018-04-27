%% 2Dsolver 
% 2D solver with primitive pressure correction
% by Eray Inanc (eray.inanc@uni-due.de)
% 4 convection schemes are included
% The divergence field for pressure correction is computed either ... 
% cell centers or cell faces
% Def: CFL 0.2, D = 2e-5, dx=2e-3, epsMax=1e-3 is working for jets
clear;clc;clf;

%% parameters run
u           = 1;            % axial velocity
v           = 0;            % radial velocity
rjet        = 10;           % jet annulus radius
CFL         = 0.2;          % timestep condition
tmax        = 10.0;         % max val of t
D           = 2e-5;         % viscosity
rho         = 1.2;          % density

%% parameters net 
Ima         = 100;          % axial points
Jma         = 100;          % radial points
dx          = 2e-4;         % spacing between points
nG          = 1;            % ghost cells (depending on the scheme)

%% case 
% 0:vertical / 1:sin / 2:jet / 3:jet-filled / 4:vortical
vprofile    = 3;            % Initial momentum          
pprofile    = 0;            % Initial passive scalar 0:left / 1:middle
iniTurb     = true;         % Random turbulence at inlet

%% parameters scheme
% 0:cds / 1:uds / 2:dds / 3:udcds 
schemePsi   = 3;            % Passive scalars
schemeVel   = 3;            % Momentum
scheme      = [schemePsi,schemeVel,schemeVel];

%% pressure corrections options
divP        = 0;            % compute div from 0:fluxes / 1:cell-centers
nItMax      = 100000;       % maximum Poisson iterations
epsMax      = 1.0e-4;       % stop limit for Possion interations
omega       = 0.5;          % relaxation factor

%% boundary conditions
% 0: dummy / 1: Dirichlet (fixed) / 2: Neumann (0 grad) / 3: Periodic (use in pairs)
momBC       = [2,2,2,0];    % Momentum BCs / E / W / N / S
psiBC       = [3,3,2,0];    % Scalar BCs     E / W / N / S
fixedValMom = [0,0,0,0];    % fixed value for Momentum at BC if momBC=1
fixedVal    = [0,0,0,0];    % fixed value for Scalar at BC if psiBC=1

%% post processing options
postProc    = true;         % see fields on the run
postProcEnd = true;         % see fields at the end
postOut     = 100;          % output every postOut step
figPos      = [0 500 1200 400]; % figure position

%% initialise and create class of Phi elements
[rhoPhi.rhoU, rhoPhi.rhoV, rhoPhi.rhoPsi] = ...
    initialise( u, v, rjet, Ima, Jma, nG, vprofile, pprofile );

% save ini velocity profile
iniU = rhoPhi.rhoU;

% allocate predicted quantity matrices
arrayTab = fieldnames(rhoPhi);
for ii = 1: length(fieldnames(rhoPhi))
    str = [arrayTab{ii}];
    rhoPhiP.(str) = rhoPhi.(str);
end

%% naming convention
% Ifim - Ifi - Ifip - .......- Ilam - Ila - Ilap (Ifim/Ilap ghost cells)
Ifim=1; Ifi=nG+1; Ifip=Ifi+1; Ilam= Ima; Ila=Ima+nG; Ilap=Ila+nG;
Jfim=1; Jfi=nG+1; Jfip=Jfi+1; Jlam= Jma; Jla=Jma+nG; Jlap=Jla+nG;

%% main program
% stop button to end the loop
h = uicontrol('Style', 'PushButton', 'String', 'Stop', ...
              'Callback', 'delete(gcbo)');

% initial time width
dt = CFL*dx/max(max(max(rhoPhi.rhoU,rhoPhi.rhoV))); 

nsteps = round(tmax/dt,1);
t=0;nt=0;dt=0;
for n=1:round(tmax/dt)
    % velocities at cell faces
    [fluxU,fluxV] = mom2vel(rhoPhi.rhoU,rhoPhi.rhoV,rho);

    % density weighted scalar
    for ii = 1:length(fieldnames(rhoPhi))
        str = [arrayTab{ii}];
        
        % convective fluxes
        [fluxConX,fluxConY] = calcFluxCon(rhoPhi.(str),fluxU,fluxV,dt,dx,scheme(4-ii));
        % diffusive fluxes
        [fluxDifX,fluxDifY] = calcFluxDif(rhoPhi.(str),dx,D);
        % apply fluxes
        rhoPhiP.(str)(Ifi:Ila,Jfi:Jla) = rhoPhi.(str)(Ifi:Ila,Jfi:Jla) ... 
            - rho*dt/dx*( fluxConX(Ifi:Ila,Jfi:Jla)-fluxConX(Ifim:Ilam,Jfi:Jla)...
                        + fluxConY(Ifi:Ila,Jfi:Jla)-fluxConY(Ifi:Ila,Jfim:Jlam)) ...
            + rho*dt/dx*( fluxDifX(Ifi:Ila,Jfi:Jla)-fluxDifX(Ifim:Ilam,Jfi:Jla)...
                        + fluxDifY(Ifi:Ila,Jfi:Jla)-fluxDifY(Ifi:Ila,Jfim:Jlam));
    end
    
    % Small primitive perturbations at inlet to trigger pseudo turbulence
    if iniTurb; rhoPhiP.rhoV(Ifim,:) = 0.5*(rand(1,Jma+2*nG)-0.5).*iniU(1,:); end;

    % update t, nt and dt
    t   = t + dt;
    nt  = nt + 1;
    dt  = CFL*dx/max(max(max(rhoPhiP.rhoU,rhoPhiP.rhoV))); 
          
    % Calculate divergence inside domain
    switch(divP)
        case 0 
            [fluxUP,fluxVP] = mom2vel(rhoPhiP.rhoU,rhoPhiP.rhoV,rho);
            [div] = calcDivergenceFlux(fluxUP,fluxVP);
        case 1 
            [div] = calcDivergence(rhoPhiP.rhoU,rhoPhiP.rhoV);
    end

    % Poisson solver
    [p,nIt,eps] = PoissonSolver(div,rhoPhiP,rho,dt,dx,omega,nItMax,epsMax); 

    % apply pseudo pressure to predicted mom
    rhoPhi.rhoU(Ifi:Ila,Jfi:Jla) = rhoPhiP.rhoU(Ifi:Ila,Jfi:Jla) - ...
        dt/dx/2*(p(Ifip:Ilap,Jfi:Jla)-p(Ifim:Ilam,Jfi:Jla));
    rhoPhi.rhoV(Ifi:Ila,Jfi:Jla) = rhoPhiP.rhoV(Ifi:Ila,Jfi:Jla) - ...
        dt/dx/2*(p(Ifi:Ila,Jfip:Jlap)-p(Ifi:Ila,Jfim:Jlam));
    
    % apply predicted fields
    rhoPhi.rhoPsi = rhoPhiP.rhoPsi;
            
    % apply boundary conditions
    % 1:E / 2:W / 3:N / 4:S
    for bc= 1:4; 
        switch(momBC(bc)); 
            case 1 % Dirichlet (fixed)
                for ii = 1:2;
                    str = [arrayTab{ii}];
                    if bc==1; rhoPhi.(str)(:,Jfim) = fixedValMom(bc); end;
                    if bc==2; rhoPhi.(str)(:,Jlap) = fixedValMom(bc); end;
                    if bc==3; rhoPhi.(str)(Ilap,:) = fixedValMom(bc); end;
                    if bc==4; rhoPhi.(str)(Ifim,:) = fixedValMom(bc); end;
                end
            case 2 % Neumann (0 grad)
                for ii = 1:2;
                    str = [arrayTab{ii}];
                    if bc==1; rhoPhi.(str)(:,Jfim) = rhoPhi.(str)(:,Jfi); end;
                    if bc==2; rhoPhi.(str)(:,Jlap) = rhoPhi.(str)(:,Jla); end;
                    if bc==3; rhoPhi.(str)(Ilap,:) = rhoPhi.(str)(Ila,:); end;
                    if bc==4; rhoPhi.(str)(Ifim,:) = rhoPhi.(str)(Ifi,:); end;
                end
            case 3 % Periodic (use in pairs)
                for ii = 1:2;
                    str = [arrayTab{ii}];
                    if bc==1; rhoPhi.(str)(:,Jfim) = rhoPhi.(str)(:,Jla); end;
                    if bc==2; rhoPhi.(str)(:,Jlap) = rhoPhi.(str)(:,Ifi); end;
                    if bc==3; rhoPhi.(str)(Ilap,:) = rhoPhi.(str)(Ifi,:); end;
                    if bc==4; rhoPhi.(str)(Ifim,:) = rhoPhi.(str)(Ila,:); end;
                end
        end
        switch(psiBC(bc)); 
            case 1 % Dirichlet (fixed)
                    if bc==1; rhoPhi.rhoPsi(:,Jfim) = fixedVal(bc); end;
                    if bc==2; rhoPhi.rhoPsi(:,Jlap) = fixedVal(bc); end;
                    if bc==3; rhoPhi.rhoPsi(Ilap,:) = fixedVal(bc); end;               
                    if bc==4; rhoPhi.rhoPsi(Ifim,:) = fixedVal(bc); end;
            case 2 % Neumann (0 grad)
                    if bc==1; rhoPhi.rhoPsi(:,Jfim) = rhoPhi.(str)(:,Jfi); end;
                    if bc==2; rhoPhi.rhoPsi(:,Jlap) = rhoPhi.(str)(:,Jla); end;
                    if bc==3; rhoPhi.rhoPsi(Ilap,:) = rhoPhi.(str)(Ila,:); end;               
                    if bc==4; rhoPhi.rhoPsi(Ifim,:) = rhoPhi.(str)(Ifi,:); end;
            case 3 % Periodic (use in pairs)
                    if bc==1; rhoPhi.rhoPsi(:,Jfim) = rhoPhi.(str)(:,Jla); end;
                    if bc==2; rhoPhi.rhoPsi(:,Jlap) = rhoPhi.(str)(:,Ifi); end;
                    if bc==3; rhoPhi.rhoPsi(Ilap,:) = rhoPhi.(str)(Ifi,:); end;               
                    if bc==4; rhoPhi.rhoPsi(Ifim,:) = rhoPhi.(str)(Ila,:); end;
        end
    end    

    % Cut backflow 
    rhoPhi.rhoU(Ilap,:) = max(0,rhoPhi.rhoU(Ila,:));
    
    % Post-Procs
    if postProc && (mod(nt,postOut)==0 || nt == 1)
        xx = linspace(-round(Jma/2),round(Jma/2),Jma+2*nG)*dx;
        yy = linspace(0,Ima+2*nG,Ima+2*nG)*dx;
    
        % Transported scalar
        ax1 = subplot(1,3,1); 
        contourf(ax1,xx,yy,rhoPhi.rhoPsi,'LineStyle','none','LevelStep',.1);
        title({['\itnt = ',num2str(nt)]});
        xlabel({'mm','Psi'});
        colorbar;
        colormap(ax1,'hot')
        ax1.CLim = [0,1];
        axis equal;

        % Axial velocity
        ax2 = subplot(1,3,2);
        contourf(ax2,xx,yy,rhoPhi.rhoU,'LineStyle','none','LevelStep',.1);
        title({['\itt = ',num2str(t),'ms']});
        xlabel({'mm','U'});
        colorbar;
        colormap(ax2,'jet')
        ax2.CLim = [min(min(rhoPhi.rhoU)),max(max(rhoPhi.rhoU))];
        axis equal;
        
        % Div. field
        ax3 = subplot(1,3,3);
        contourf(ax3,xx,yy,div,'LineStyle','none','LevelStep',.1); colorbar;
        xlabel({'mm','div U_{pred}'});
        colorbar; 
        axis equal;

        % Update plots
        shg;

        % Plot position
        bx = gcf; bx.Color = [1 1 1]; bx.Resize = 'off'; bx.ToolBar = 'none'; bx.MenuBar = 'none';
        bx.Position = figPos;
    end

    % check the stop button
    drawnow
    if ~ishandle(h)
        break;
    end
    
    % printout some stuff
    display(['nt: ',num2str(nt),' | dt: ',num2str(dt),' | nIt: ',...
        num2str(nIt), ' | eps: ',num2str(eps)]);
end;

%% Post-Procs end
if postProcEnd
    xx = linspace(-round(Jma/2),round(Jma/2),Jma+2*nG)*dx;
    yy = linspace(0,Ima+2*nG,Ima+2*nG)*dx;
   
    % Transported scalar
    ax1 = subplot(1,3,1); %#ok<UNRCH>
    contourf(ax1,xx,yy,rhoPhi.rhoPsi,'LineStyle','none','LevelStep',.1);
    title({['\itnt = ',num2str(nt)]});
    ylabel('mm');
    xlabel({'mm','Psi'});
    colorbar;
    colormap(ax1,'hot')
    ax1.CLim = [0,1];
    axis equal;
   
    % Axial velocity
    ax2 = subplot(1,3,2);
    contourf(ax2,xx,yy,rhoPhi.rhoU,'LineStyle','none','LevelStep',.1);
    title({['\itt = ',num2str(t),'ms']});
    xlabel({'mm','U'});
    colorbar;
    colormap(ax2,'jet')
    ax2.CLim = [min(min(rhoPhi.rhoU)),max(max(rhoPhi.rhoU))];
    axis equal;
    
    % Div. field
    ax3 = subplot(1,3,3);
    contourf(ax3,xx,yy,div,'LineStyle','none','LevelStep',.1); colorbar;
    xlabel({'mm','div U_{pred}'});
    colorbar; 
    axis equal;
    
    shg;
end
% END %
