%my colormap, purple to cyan.
function [colors]=plasma(m)

    if nargin < 1
       f = get(groot,'CurrentFigure');
       if isempty(f)
          m = size(get(groot,'DefaultFigureColormap'),1);
       else
          m = size(f.Colormap,1);
       end
    end

    Pu=[ 0.94   0.98    0.13]';
    Gr=[0.8000    0.2784    0.4706]';
    Cy=[0.0510    0.0314    0.5294]';

    
    x=linspace(0,1,m);
    colors=interp1([0 .5 1],[Pu, Gr,Cy]',x);
end
