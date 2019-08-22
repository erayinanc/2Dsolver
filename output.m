%% Output
%  Output properties!

function []=output(rhoPhi,div,t,nt,figPos) 
    global Ima Jma dx nG Ifim Ifi Ifip Ilam Ila Ilap Jfim Jfi Jfip Jlam Jla Jlap;

    xx = linspace(-round(Jma/2),round(Jma/2),Jma+2*nG)*dx;
    yy = linspace(0,Ima+2*nG,Ima+2*nG)*dx;

    ax1 = subplot(1,3,1); %#ok<*UNRCH>
    contourf(ax1,xx,yy,rhoPhi.rhoPsi,'LineStyle','none','LevelStep',.1);
    title({['\itnt = ',num2str(nt)]});
    xlabel({'mm','Psi'});
    colorbar;
    colormap(ax1,'hot')
    ax1.CLim = [0,1];
    axis equal;

    ax2 = subplot(1,3,2);
    contourf(ax2,xx,yy,rhoPhi.rhoU,'LineStyle','none','LevelStep',.1);
    title({['\itt = ',num2str(t),'ms']});
    xlabel({'mm','U'});
    colorbar;
    colormap(ax2,'jet')
    ax2.CLim = [min(min(rhoPhi.rhoU)),max(max(rhoPhi.rhoU))];
    axis equal;

    ax3 = subplot(1,3,3);
    contourf(ax3,xx,yy,div,'LineStyle','none','LevelStep',.1); colorbar;
    xlabel({'mm','div U_{pred}'});
    colorbar; 
    axis equal;

    shg;

    bx = gcf; bx.Color = [1 1 1]; bx.Resize = 'off'; bx.ToolBar = 'none'; bx.MenuBar = 'none';
    bx.Position = figPos;
end   





























