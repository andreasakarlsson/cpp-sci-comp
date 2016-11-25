#! /usr/bin/octave -qf
fname = input('', 's');  ## read from stdin
fid = fopen(strcat(fname,'.bin'),'r');
c = fread(fid,'double');
x = c(1:length(c)/2);
y = c(length(c)/2+1:end);

f = figure;
set(f, 'visible', 'off')
plot(x,y,'.');
my_axis = [-11 6 -1 4];
axis(my_axis)
axis equal
set(gcf, 'PaperPosition', [-1 0 14 4]);
set(gcf, 'PaperSize', [12 4]);
print(strcat(fname, '.pdf'), '-dpdf');
exit;
