function [fig]=zoom_Function(I,ptx,lx,pty,ly,title)
    R = insertShape(I,"Rectangle",[ptx pty lx ly],LineWidth=2,Color="red");
    fig=figure('Name',title);
    colormap gray;
    axis off image; 
    zoom(0);
    imshow(R);
    hold on
    axes('position',[.10 .10 .35 .35])
    box on % put box around new pair of axes
    x = I(pty:pty+ly,ptx:ptx+lx); % range box
    imshow(x)
    axis tight
end
