function fPauseVideo(videoPauseTime,thisVideo)

% pause video
%----------------------
% pauseVideo = 1;
for dummyFrames = 1:round(videoPauseTime*thisVideo.FrameRate)
    %frame = getframe;
    %writeVideo(thisVideo,frame);
    writeVideo(thisVideo,getframe);
end

end
