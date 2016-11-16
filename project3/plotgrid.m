#! /usr/bin/octave -qf
fid = fopen('outfile.bin','r'); c = fread(fid,'double');
x = c(1:length(c)/2);  y = c(length(c)/2+1:end);

f = figure;
set(f, "visible", "off")
plot(x,y,'o');
print("gridpoints.png", "-dpng");
