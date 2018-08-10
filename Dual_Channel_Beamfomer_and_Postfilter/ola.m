function yout2 = ola(yin,fragsize,overlap)

noover  = 1 - overlap;
Move    = fragsize * noover;
Nframe  = size(yin,1);
NFFT    = size(yin,2);
yout    = zeros(1, Nframe*Move+NFFT);

for n = 1:Nframe
    framestart  = 1 + (n-1)*Move;
    frameend    = fragsize + (n-1)*Move;
    yout(framestart:frameend) = yout(framestart:frameend) + yin(n,1:fragsize);
end

yout2 = yout(1:Nframe*Move+fragsize*overlap);