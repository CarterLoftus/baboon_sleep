function [VEDBApoint]= calc_vedba(AccInMs2, MovAvgWindow)
if nargin==1
    MovAvgWindow=12; %points for the moving avg. if 3 point per sec than 15 is 5 sec
end
halfWindowSize=ceil(MovAvgWindow/2);
VEDBApoint=nan(size(AccInMs2));%here i store odba value for each point
%% loop on the point in the measurment excluding the edges to allow MovAvg calc
for PointCnt=(halfWindowSize+1):(length(AccInMs2)-halfWindowSize-1)
    Moveavg=nanmean(AccInMs2((PointCnt-halfWindowSize):(PointCnt+halfWindowSize-1)));
    %difference (the Dynamic component)
    Diff=AccInMs2(PointCnt)-Moveavg;
    VEDBApoint(PointCnt)=(Diff^2)^0.5;%per point ODBA
end
