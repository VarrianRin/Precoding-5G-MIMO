function [CB, N1, N2, O1, O2] = codebook(M)

% codebook builder

if M ~= 16
    error('!!!UNKNOWN NUMBER OF ANTENNAS!!!');
end

N1 = 8;
N2 = 2;
O1 = 4;
O2 = 4;

Vh = zeros(N1, O1 * N1);
for i11 = 0 : N1*O1 - 1
    Vh(:, i11 + 1) = exp(1i*2*pi*(0 : N1-1)*i11 / (O1*N1)).';
end

Vv = zeros(N2, O2 * N2);
for i12 = 0 : N2*O2 - 1
    Vv(:, i12 + 1) = exp(1i*2*pi*(0 : N2-1)*i12 / (O2*N2)).';
end

CB = kron(Vv, Vh);

end

