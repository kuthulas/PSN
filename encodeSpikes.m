function X = encodeSpikes(stimulus)
clc
nRF = 7;
tdelta = 5;
Centers = zeros(size(stimulus,1),nRF);

max = 20*ones(size(stimulus));
min = zeros(size(stimulus));

for i=0:nRF-1
    Centers(:,i+1) = min + (max-min).*((2*(i+1)-3)/2)/(nRF-2);
end
beta = 1;
width = 1/beta * (max-min)./(nRF-2);
GRFs = getGRFs(stimulus,Centers,width);
spike_times = round(getSpikeTimes(GRFs,tdelta,0.2));

% Decoding experiment
sr = (tdelta+rand).*ones(1,size(stimulus,1)*nRF);
sr(1,spike_times(1,:)) = spike_times(2,:);
sr = reshape(sr,nRF,size(stimulus,1))';
mr = tdelta/(0.2-1);
br = tdelta - mr*0.2;
Gr = (sr - br)./mr;
Gr = sqrt(-log(Gr));
Grn = -Gr;
Gr = bsxfun(@times,Gr,width);
Grn = bsxfun(@times,Grn,width);
x = Gr + Centers;
xn = Grn + Centers;
recon = Inf.*ones(1,size(stimulus,1)*nRF);
recon(1,spike_times(1,:)) = 1;
recon = reshape(recon,nRF,size(stimulus,1))';
Z = [x.*recon,xn.*recon];
X = zeros(size(stimulus));
for i=1:size(stimulus,1)
    Zi = Z(i,:);
    Zi(abs(Zi)==Inf)=[];
    X(i) = mode(Zi);
end
X = abs(stimulus-X);
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