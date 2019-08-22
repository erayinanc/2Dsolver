%% 2Dsolver 
% 2D solver with primitive pressure correction
% by Eray for lecture Numerics and Flow Sim.
% CFL 0.2, D = 2e-5, dx=2e-3, epsMax=1e-3 is all right for jets
clear;clc;clf;

%% parameters run
u           = 1;            % axial velocity
v           = 0;            % radial velocity
rjet        = 20;            % jet annulus radius
CFL         = 0.2;          % timestep condition
tmax        = 10.0;         % max val of t
D           = 2e-5;         % viscosity
rho         = 1.2;          % density

%% parameters net 
global Ima Jma dx nG;
Ima         = 100;           % axial points
Jma         = 80;            % radial points
dx          = 2e-4;         % spacing between points
nG          = 1;            % ghost cells (depending on the scheme)

%% case 
% 0:vertical / 1:sin / 2:jet / 3:jet-filled / 4:vortical
vprofile    = 3;            % Initial momentum          
pprofile    = 0;            % Initial passive scalar 0:left / 1:middle
iniTurb     = false;         % Random turbulence at inlet

%% parameters scheme
% 0:cds / 1:uds / 2:udcds 
schemePsi   = 1;            % Passive scalars
schemeVel   = 0;            % Momentum
scheme      = [schemePsi,schemeVel,schemeVel];

%% pressure corrections options
divP        = 1;            % compute div from 0:fluxes / 1:cell-centers
nItMax      = 100000;       % maximum Poisson iterations
epsMax      = 1.0e-4;       % stop limit for Possion interations
omega       = 1.0;          % relaxation factor

%% boundary conditions
% 0: dummy / 1: Dirichlet (fixed) / 2: Neumann (0 grad) /  
% 3: Neumann (0 grad) clipped backflow / 4: Periodic (use in pairs)
momBC       = [0,3,4,4];    % Momentum BCs / E / W / N / S
psiBC       = [0,2,4,4];    % Scalar BCs   / E / W / N / S
fixedValMom = [0,0,0,0];    % fixed value for Momentum at BC if momBC=1
fixedVal    = [0,0,0,0];    % fixed value for Scalar at BC if psiBC=1
%% post processing options
postProc    = true;         % see fields on the run
postProcEnd = true;         % see fields at the end
postOut     = 100;          % output every postOut step
figPos      = [0 1000 2400 800]; % figure position

%% initialise and create class of Phi elements
[rhoPhi.rhoU,rhoPhi.rhoV,rhoPhi.rhoPsi] = initialise(u,v,rjet,vprofile,pprofile);

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
global Ifim Ifi Ifip Ilam Ila Ilap Imid Jfim Jfi Jfip Jlam Jla Jlap Jmid;
Ifim=1; Ifi=nG+1; Ifip=Ifi+1; Ilam= Ima; Ila=Ima+nG; Ilap=Ila+nG; Imid=ceil(Ilap/2);
Jfim=1; Jfi=nG+1; Jfip=Jfi+1; Jlam= Jma; Jla=Jma+nG; Jlap=Jla+nG; Jmid=ceil(Jlap/2);
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
        [fluxConX,fluxConY] = calcFluxCon(rhoPhi.(str),fluxU,fluxV,dt,scheme(4-ii));
        % diffusive fluxes
        [fluxDifX,fluxDifY] = calcFluxDif(rhoPhi.(str),D);
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
    t = t + dt;
    nt = nt + 1;
    dt = CFL*dx/max(max(max(rhoPhiP.rhoU,rhoPhiP.rhoV))); 
          
    % Calculate divergence inside domain
    switch(divP)
        case 0 
            [fluxUP,fluxVP] = mom2vel(rhoPhiP.rhoU,rhoPhiP.rhoV,rho);
            [div]=calcDivergenceFlux(fluxUP,fluxVP,rho);
        case 1 
            [div]=calcDivergence(rhoPhiP.rhoU,rhoPhiP.rhoV);
    end

    % Poisson solver
    [p,nIt,eps]=PoissonSolver(div,rho,dt,omega,nItMax,epsMax); 

    % apply pseudo pressure to predicted mom
    rhoPhi.rhoU(Ifi:Ila,Jfi:Jla) = rhoPhiP.rhoU(Ifi:Ila,Jfi:Jla) - ...
        dt/dx/2*(p(Ifip:Ilap,Jfi:Jla)-p(Ifim:Ilam,Jfi:Jla));
    rhoPhi.rhoV(Ifi:Ila,Jfi:Jla) = rhoPhiP.rhoV(Ifi:Ila,Jfi:Jla) - ...
        dt/dx/2*(p(Ifi:Ila,Jfip:Jlap)-p(Ifi:Ila,Jfim:Jlam));
    
    % apply predicted fields
    rhoPhi.rhoPsi = rhoPhiP.rhoPsi;
            
    % apply boundary conditions
    [rhoPhi]=boundary(rhoPhi,arrayTab,momBC,psiBC,fixedVal,fixedValMom);
    
    % Cut backflow 
    rhoPhi.rhoU(Ilap,:) = max(0,rhoPhi.rhoU(Ila,:));
    
    % Post-Procs
    if postProc && (mod(nt,postOut)==0 || nt == 1)
        output(rhoPhi,div,t,nt,figPos); 
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

% Post-Procs end
if postProcEnd
    output(rhoPhi,div,t,nt,figPos);
end

%%END
