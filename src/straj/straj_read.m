function H = dhist_read(fn)

if nargin < 1
  error('Specify filename');
end

f = fopen(fn, 'r');
H.dr = fread(f, 1, 'double');
H.np = fread(f, 1, 'int');
H.n = fread(f, 1, 'int');
H.nfun = fread(f, 1, 'int');
H.lm = fread(f, H.nfun, 'int');
H.ld = fread(f, H.nfun, 'int');
H.hst = fread(f, H.n, 'uint16');
fclose(f);
