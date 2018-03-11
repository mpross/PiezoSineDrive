close all;
data=load('FA_2018_69c.dat');
rawTime=(data(:,1)-data(2,1))*24*3600;
rawBRS=data(:,2)-data(1,2);

sampF=1/(rawTime(4)-rawTime(3));

[ABRS, F] = asd2(rawBRS,1/sampF, 1, 1, @hann);

w0=F(find(ABRS==max(ABRS(find(F>1e-3)))));
Q=1e5;
I=0.59;
M=4.5;
g=9.81;

driveF=[0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045];
periods=15;
lastTime=0;
amp=[];
for j=1:8
    
    if(j>1)
        startTime=lastTime;
        endTime=sampF/driveF(j)*periods+lastTime;
    else
        startTime=10;
        endTime=sampF/driveF(j)*periods;
    end    
    lastTime=endTime;
    time=rawTime(startTime:endTime);
    BRS=7e-6/4*rawBRS(startTime:endTime);

    [ABRS, F] = asd2(BRS,1/sampF, 1, 1, @hann);
    
    amp=[amp;max(ABRS(find(and(F<5e-3,F>0.9e-3))))/sqrt((endTime-startTime)/sampF)];  
    
    hold on
    figure(1)
    plot(time,BRS);
    
    figure(3)
    loglog(F,abs(ABRS));
    
end

startTime=5.7e4*sampF;
endTime=length(rawTime);

time=rawTime(startTime:endTime);
BRS=7e-6/4*rawBRS(startTime:endTime);

[ABRS, F] = asd2(BRS,1/sampF, 1, 1, @hann);

norm=max(ABRS(find(F>5e-2)))/sqrt((endTime-startTime)/sampF);

amp=[amp;norm]./norm;

driveF=[driveF 0.2];

figure(2)
plot(rawTime,rawBRS);

fun=@(x,w)(w.^2-M*g*x(1)/I)./sqrt((w.^2-w0^2).^2+w0^4/Q^2);
x0=3e-7;
options = optimset('Display','iter','TolX', 1e-30, 'TolFun', 1e-10, 'MaxFunEvals', 4000, 'MaxIter', 4000);
d=lsqcurvefit(fun,x0,driveF,amp',0,[],options)
w=1e-4:1e-4:1;

d=3.6e-7

figure(4)
loglog(driveF,amp','.',w,abs((w.^2-M*g*d/I)./sqrt((w.^2-w0^2).^2+w0^4/Q^2)));
ylim([1e-1,1e1]);
grid on