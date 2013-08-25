function W = avgw_read(fn)

if nargin < 1
  error('Specify filename');
end

f = fopen(fn, 'r');
W.np = fread(f, 1, 'int');
W.dr = fread(f, 1, 'single');
W.nfun = fread(f, 1, 'int');
W.n = fread(f, W.nfun, 'int');
W.i0 = W.n;
W.i0(1) = 1;
W.i0(2:end) = cumsum(W.n(1:end-1))+1;
W.in = W.i0 + W.n - 1;
m = sum(W.n);
W.s = fread(f, m, 'single');
W.nz = fread(f, W.nfun, 'int');
W.Icut = fread(f, W.nfun, 'single');
fclose(f);
