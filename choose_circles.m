function [c r] = choose_circles(cin,rin)
% circle maker: 
% present window with default axis.
% first click in window picks first circle center; circle now shown 
%   with chosen center, changes radius with mouse; next mouse click
%   sets radius.
% subsequent clicks make subsequent circles. (should check for overlapping
%   circles)
% pressing enter stops the process (like ginput)
% pressing + zooms axis in (limit to circles extent?)
% pressing - zooms axis out
% if given some circles, plot those and return with any extra circles
%   added
%
%
% This file is part of MCSC.
% 
% MCSC is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% MCSC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with MCSC.  If not, see <http://www.gnu.org/licenses/>.

% Copyright Everett Kropf, 2013


%% clicky figure w/default axis
cfh = figure;
init_figpos(cfh)
cah = axes;
set(cah,'dataaspectratio',[1 1 1])
axis([-5 5 -5 5])
grid on


%% gui data
gdata.center_set = 0;
gdata.cah = cah;
gdata.circh = [];
gdata.eit = exp((0:200)'*2i*pi/200);
gdata.axlim = 5;

if nargin
  gdata.c = cin;
  gdata.r = rin;
  gdata.n = length(rin);
  for j = 1:gdata.n
    circ = cin(j) + rin(j)*gdata.eit;
    line(real(cin(j)),imag(cin(j)),'color','k','marker','x')
    line(real(circ),imag(circ),'color','k','linestyle','-')
  end
else
  gdata.c = [];
  gdata.r = [];
  gdata.n = 0;
end

guidata(cfh,gdata);


%% set...go!
set(cfh,'pointer','fullcrosshair')
set(cfh,'windowbuttondownfcn',@do_click)
set(cfh,'windowkeypressfcn',@do_keypress)

uiwait(cfh)


%% cleanup
if ishandle(cfh)
  set(cfh,'pointer',arrow)
  set(cfh,'currentaxes',cah)
  grid off
  set(cfh,'windowbuttondownfcn','')
  set(cfh,'windowkeypressfcn','')
end

if nargout
  if ishandle(cfh)
    gdata = guidata(cfh);
    c = gdata.c(:);
    r = gdata.r(:);
  else
    c = [];
    r = [];
  end
end

end


%---------------------------------------------------------------
function do_click(src,evt) %#ok<INUSD>
  gdata = guidata(src);
  pt = get(gdata.cah,'currentpoint');
  x = pt(1,1); y = pt(1,2);
  
  if ~gdata.center_set
    gdata.center_set = 1;
    set(src,'windowbuttonmotionfcn',@do_move)
    
    line(x,y,'marker','x','color','k')
    gdata.n = gdata.n + 1;
    gdata.c(gdata.n) = complex(x,y);
  else
    gdata.center_set = 0;
    set(src,'windowbuttonmotionfcn','')
    
    n = gdata.n;
    gdata.r(n) = abs(gdata.c(n) - complex(x,y));
    circ = gdata.c(n) + gdata.r(n)*gdata.eit;
    line(real(circ),imag(circ),'color','k','linestyle','-')
    
    if ~isempty(gdata.circh)
      delete(gdata.circh)
      gdata.circh = [];
    end
  end

  guidata(src,gdata)
end

%---------------------------------------------------------------
function do_move(src,evt) %#ok<INUSD>
  gdata = guidata(src);
  pt = get(gdata.cah,'currentpoint');
  x = pt(1,1); y = pt(1,2);
  
  if ~isempty(gdata.circh)
    delete(gdata.circh)
  end
  
  r = abs(gdata.c(gdata.n) - complex(x,y));
  circ = gdata.c(gdata.n) + r*gdata.eit;
  gdata.circh = line(real(circ),imag(circ),'color','r','linestyle','-');
  
  guidata(src,gdata)
end

%---------------------------------------------------------------
function do_keypress(src,evt)
  gdata = guidata(src);
  set_axis = 0;
  
  switch evt.Key
    case 'return'
      if isempty(evt.Modifier) && ~gdata.center_set
        uiresume(src)
        disp('All done.')
      end
    case 'subtract'
      lim = gdata.axlim + 1;
      set_axis = 1;
    case 'add'      
      if gdata.axlim > 1
        lim = gdata.axlim - 1;
        set_axis = 1;
      end
  end
  
  if set_axis
    axis(gdata.cah,[-lim lim -lim lim])
    gdata.axlim = lim;
    guidata(src,gdata)
  end
end

%---------------------------------------------------------------
function init_figpos(fh)
  set(fh,'units','normalized')
  pos = get(fh,'position');
  set(fh,'position',[0.05 0.2 pos(3) 0.65])
  set(fh,'units','pixels')
  pos = get(fh,'position');
  set(fh,'position',[pos(1:2) 1.2*pos(4) pos(4)])
end
