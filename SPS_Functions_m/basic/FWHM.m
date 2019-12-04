function width = FWHM(x,y)
% calculates the Full-Width at Half-Maximum (FWHM) of the waveform y(x)
% width = FWHM(x,y)
% 'width' is in units of 'x' 


cutoff=0.5; %this is the normalized height at which the width is calculated.

y = y / max(y);
N = length(y);
if y(1) < cutoff                  % find index of center (max or min) of pulse
    [~,centerindex]=max(y);
else
    [~,centerindex]=min(y);
end

%find first crossing
ii = 2;
while sign(y(ii)-cutoff) == sign(y(ii-1)-cutoff)
    ii = ii+1;
end

interp = (cutoff-y(ii-1)) / (y(ii)-y(ii-1));
t1 = x(ii-1) + interp*(x(ii)-x(ii-1));

%start search for next crossing at center
ii = centerindex+1;
while ((sign(y(ii)-cutoff) == sign(y(ii-1)-cutoff)) && (ii <= N-1))
    ii = ii+1;
end

if ii ~= N
    interp = (cutoff-y(ii-1)) / (y(ii)-y(ii-1));
    t2= x(ii-1) + interp*(x(ii)-x(ii-1));
    width = t2- t1;
else
    width = NaN;
end

end


