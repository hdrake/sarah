
ifig = 2;
figure(ifig), clf
ipr = 'n';
fsize = 15;

cmap = redblue(1024);
%cmap(1,:) = 1; % set lower limit of colormap to white

runid = input(' enter the data file name (ext = .dat): ','s');

cdir = pwd;

deta = 0.005;
vc0 = [0.001 0.001]';
veta = [deta:deta:3];


fid=fopen([runid,'.dat'], 'rb','l');  % 'l' is for little endian
                                       % binary format on the alpha
% read all the initial information
hr1    = fread(fid, 1, 'int32');
  nout   = fread(fid, 1, 'int32');
  tend   = fread(fid, 1, 'float64');
hr1    = fread(fid, 1, 'int32');
hr1    = fread(fid, 1, 'int32');
  mx     = fread(fid, 1, 'int32');
  my     = fread(fid, 1, 'int32');
  meqn   = fread(fid, 1, 'int32');
  mbc    = fread(fid, 1, 'int32');
hr1    = fread(fid, 1, 'int32');
hr1    = fread(fid, 1, 'int32');
  dt     = fread(fid, 1, 'float64');
  dx     = fread(fid, 1, 'float64');
  dy     = fread(fid, 1, 'float64');
  gamma  = fread(fid, 1, 'float64');
  cf     = fread(fid, 1, 'float64');
  cb     = fread(fid, 1, 'float64');
hr1    = fread(fid, 1, 'int32');
hr1    = fread(fid, 1, 'int32');
  ixbc    = fread(fid, 1, 'int32');
  iybc     = fread(fid, 1, 'int32');
hr1    = fread(fid, 1, 'int32');
hr1 = fread(fid, 1, 'int32');
  x   = fread(fid, [mx,1], 'float64');
hr1 = fread(fid, 1, 'int32');
hr1 = fread(fid, 1, 'int32');
  y   = fread(fid, [my,1], 'float64');
hr1 = fread(fid, 1, 'int32');
hr1 = fread(fid, 1, 'int32');
  b   = fread(fid, [mx,my], 'float64');
hr1 = fread(fid, 1, 'int32');
hr1 = fread(fid, 1, 'int32');
  cbxy = fread(fid, [mx,my], 'float64');
hr1 = fread(fid, 1, 'int32');

x0 = x(1)  -dx/2;
xf = x(end)+dx/2;
y0 = y(1)  -dy/2;
yf = y(end)+dy/2;

for i = 1:nout+1
 disp('start read')
% ** read the h, Qu, Qv fields *****************
   hr1 = fread(fid, 1, 'int32');
    t   = fread(fid, 1, 'float64');
   hr1 = fread(fid, 1, 'int32');
   hr1 = fread(fid, 1, 'int32');
    h   = fread(fid, [mx,my], 'float64');
   hr1 = fread(fid, 1, 'int32');
   hr1 = fread(fid, 1, 'int32');
    qu  = fread(fid, [mx,my], 'float64');
   hr1 = fread(fid, 1, 'int32');
   hr1 = fread(fid, 1, 'int32');
    qv  = fread(fid, [mx,my], 'float64');
   hr1 = fread(fid, 1, 'int32');
% ******************************************
 disp('end read')

   time(i) = t;
   
   u = qu./h;
   v = qv./h;
   u(  find(h < 0.001) ) = 0;
   v(  find(h < 0.001) ) = 0;
   qu( find(h < 0.001) ) = 0;
   qv( find(h < 0.001) ) = 0;
   h(  find(h < 0.001) ) = 0;
   eta  = h+b;

  figure(1)
  subplot(1,2,1)
  hold on;
  %contourf(x,y,h',20)
  contourf(x,y,qu',10)
  contour(x,y,b',"k")
  colorbar()
  %caxis([0,0.7])
  caxis([-1.,1.])
  colormap(cmap)
  hold off;
  
  subplot(1,2,2)
  hold on;
  plot(y,b(400,:),"k");
  if ~mod(i-1,5)
      plot(y,h(400,:)+b(400,:));
  end
  set(gca, 'XDir','reverse')
  
  saveas(gcf,"figures/"+string(i)+".png");
  %pause
  
end

status = fclose(fid);




