function [ptout, Nframe] = cutframe(pt, fragsize, overlap)

noover  = 1 - overlap;
Nlength = size(pt,1);
Move    = fragsize * noover;
Nframe  = ceil( (Nlength-fragsize*overlap) / Move );
data    = [pt; zeros(Move*Nframe+fragsize*overlap-Nlength,1)];
ptout   = zeros(Nframe, fragsize);

for n = 1:Nframe
    ptout(n,:) = data(1+(n-1)*Move:fragsize+(n-1)*Move,1)';
end