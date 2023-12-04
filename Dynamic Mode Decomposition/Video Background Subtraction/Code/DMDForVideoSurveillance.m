clear all,
close all,

% %first video
% [fgmc, bgmc] = dmdBgSub('monte_carlo_low.mp4', 100, 0.2);
% n = size(fgmc, 3); %number of snapshots/frames
% fgmcvid = VideoWriter('monte_carlo_FG'); %create the foreground video
% bgmcvid = VideoWriter('monte_carlo_BG'); %create the background video
% open(fgmcvid); %should be start with when using writeVideo
% open(bgmcvid);
% for i=1:n
%     writeVideo(fgmcvid, fgmc(:, :, i));
%     writeVideo(bgmcvid, bgmc(:, :, i));
% end
% close(fgmcvid); %should be end with when using writeVideo
% close(bgmcvid);
%-----------------
%second video
[fgski, bgski] = dmdBgSub('ski_drop_low.mp4', 50, 0.2);
n = size(fgski, 3); %number of snapshots/frames
fgskivid = VideoWriter('ski_drop_FG'); %create the foreground video
bgskivid = VideoWriter('ski_drop_BG'); %create the background video
open(fgskivid); %should be start with when using writeVideo
open(bgskivid);
for i=1:n
    writeVideo(fgskivid, fgski(:, :, i));
    writeVideo(bgskivid, bgski(:, :, i));
end
close(fgskivid);
close(bgskivid); %should be end with when using writeVideo#
%-----------------