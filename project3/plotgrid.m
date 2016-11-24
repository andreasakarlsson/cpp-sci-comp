#! /usr/bin/octave -qf
fid = fopen('outfile.bin','r'); c = fread(fid,'double');
x = c(1:length(c)/2);  y = c(length(c)/2+1:end);

f = figure;
set(f, "visible", "off")
plot(x,y,'.');
my_axis = [-11 6 -1 4];
axis(my_axis)
axis equal
set(gcf, 'PaperPosition', [0 0 16 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [17 5]); %Set the paper to have width 5 and height 5.
print("gridpoints.pdf", "-dpdf");
