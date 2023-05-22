function [J,v,a,n]=bilatFilt(x,y,t,nr,ss,sr,di)
%
%[J, v, a, n] = bilatFilt(x, y, t, nr, ss, sr, di)
%Applies a bilateral filter to a time series of gaze coordinates
%
%INPUTS
%x:     Time series of gaze coordinates, horizontal axis (deg visual angle)
%y:     Time series of gaze coordinates, vertical axis (deg visual angle)
%t:     Time stamp of each sample in x and y (seconds)
%nr:    Sliding window size (number of samples, 3 or 4 works well with my
%       60Hz data (~50-70ms)
%ss:    Standard deviation of Gaussian along the time dimension (number
%       of samples, can play around with this, but 3 works pretty well
%       my 60Hz data (~50ms)
%sr:    Standard deviation of Gaussian along the spatial dimension
%       (can play around with this, 3 works pretty well with my data)
%di:    Display figures? 0 = no, 1 = yes
%
%OUTPUTS
%J:     Structure with filtered time series
%v:     Structure with unfiltered x, y, and absolute velocity (deg/s)
%a:     Structure with unfiltered x, y, and absolute acceleration (deg/s^2)
%n:     Structure with various counting variables


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Unfiltered velocity (used for visual comparison)
v.x=diff(x,1)./diff(t,1);
v.y=diff(y,1)./diff(t,1);
v.r=sqrt(v.x.^2+v.y.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Unfiltered acceleration (used for visual comparison)
a.x=diff(x,2)./diff(t,2);
a.y=diff(y,2)./diff(t,2);
a.r=sqrt(a.x.^2+a.y.^2);
a.r(a.r>1e4)=1e4;                                                           
%Acceleration is impressively noisy, this mitigates crazy values, mainly for visualization purposes to keep the axes reasonable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Initialize variables and pre-allocate filtered output
s.s=ss;
s.r=sr;
n.r=nr;
n.s=length(x);
J.x=nan(1,n.s);
J.y=nan(1,n.s);
J.r=nan(1,n.s);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Filter within a sliding window

%tlb testing adjustments for bilateral filter
% trying out sliding window that avoids invalid samples marked by nan

for i=n.r+1:n.s-n.r
    
    p=i+(-n.r:n.r); % p - current idxs to sample
    
    if isnan(x(i)) % don't do anything with an invalid sample
        J.x(i) = nan;
        J.y(i) = nan;
    else
        f=exp(-(p-i).^2/(2*s.s^2));
        g=exp(-(x(p)-x(i)).^2/(2*s.r^2)); % calculate the difference between all samples in the window and the current sample
        J.x(i)=nansum(f.*g.*x(p))/nansum(f.*g); % add smoothing over nans
        g=exp(-(y(p)-y(i)).^2/(2*s.r^2));
        J.y(i)=nansum(f.*g.*y(p))/nansum(f.*g);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Filtered velocity
J.vx=diff(J.x,1)./diff(t,1);
J.vy=diff(J.y,1)./diff(t,1);
J.v=sqrt(J.vx.^2+J.vy.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Filtered acceleration
J.ax=diff(J.x,2)./diff(t,2);
J.ay=diff(J.y,2)./diff(t,2);
J.a=sqrt(J.ax.^2+J.ay.^2);
J.a(J.a>1e4)=1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Display figures if desired
if di == 1
    
    figure;
    hold on;
    plot(J.x,'r','LineWidth',2);
    plot(x,'k');
    
    figure;
    hold on;
    plot(J.y,'r','LineWidth',2);
    plot(y,'k');
    
    figure;
    hold on;
    plot(J.v,'r','LineWidth',2);
    plot(v.r,'k');
    
    figure;
    hold on;
    plot(J.a,'r','LineWidth',2);
    plot(a.r,'k');
    
end

