function quant = quantize(Vector, CB)

CB     = squeeze(CB);
Vector = squeeze(Vector);
CBsize = size(CB, 2);
max    = 0;
index  = 0;

for i = 1 : CBsize
    value = norm(CB(:, i)' * Vector, 'fro');
    if value >= max
        max = value;
        index = i;
    end
end
quant = CB(:, index);

end
