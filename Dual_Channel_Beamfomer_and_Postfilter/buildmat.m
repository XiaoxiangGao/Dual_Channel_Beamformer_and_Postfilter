function datamatout = buildmat(Y, nrow)

if mod(nrow,2) == 0
    error('The length of the smoothing window in frequency domain must be odd!');
end
ncol = length(Y);
offset = (nrow-1)/2;
datamat = zeros(nrow,ncol+offset*2);
for n = 1:nrow
    datamat(n,1+(nrow-n):end-(n-1)) = Y;
end
datamatout = datamat(:,offset+1:end-offset);
