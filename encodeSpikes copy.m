function [spike_times, min, max] = encodeSpikes(stimulus, min, max)
nRF = 7; % neurons per variable
tdelta = 5; % msecs
Centers = zeros(size(stimulus,1),nRF);

if(isempty(min))
    min = stimulus;
    max = stimulus;
else
    max(max<stimulus) = stimulus(max<stimulus);
    min(stimulus<min) = stimulus(stimulus<min);
end

for i=0:nRF-1
    Centers(:,i+1) = min + (max-min).*((2*(i+1)-3)/2)/(nRF-2);
end
beta = 1;
widths = 1/beta * (max-min)./(nRF-2);
GRFs = getGRFs(stimulus,Centers,widths);

if(~isempty(min))
    spike_times = round(getSpikeTimes(GRFs,tdelta,0.2))
    plot(spike_times(2,:),spike_times(1,:),'o');
end
end

function spike_times = getSpikeTimes(G,t,thresh)
m = t/(thresh-1);
b = t - m*thresh;
s = m.*G + b;
st = reshape(s',1,size(s,1)*size(s,2));
n = find(st(:)' <= t);
spike_times = [n;st(n)];
end

function G = getGRFs(x,C,w)
if(sum(sum(w))==0)
    G = zeros(size(C));
else
    G = bsxfun(@minus,x,C);
    G = bsxfun(@rdivide,G,w);
    G = exp(-G.^2);
end
end