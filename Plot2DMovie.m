clear all
close all

axes = readmatrix('xy.csv');
Files = dir('Data');
num_files = length(Files);
cd Data
data{:} = zeros(num_files-2);
for i = 3:num_files
   file = Files(i).name;
   data{i-2} = readmatrix(file);
end
cd ..
%% Video
v = VideoWriter('Spikes.mp4', 'MPEG-4');
v.FrameRate = 10;  % arbitrary
open(v)
f=figure;
pause(0.2) % let plot wake up
for i=1:num_files - 2
    surf(axes(:,1),axes(:,2),data{i});
    im = getframe(f);
    writeVideo(v,im)
end
close(v)
