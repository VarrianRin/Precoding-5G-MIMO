function cb = lloyd(vec, B, graph, d, df, rank, Nrank)

% -- init --
vec     = squeeze(vec);
[M, I]  = size(vec);        % M - size of vector, I - volume of statistic
no_iter = 90;               % number of iterations
N       = 2^B;              % number of vectors in codebook

cb      = zeros(M, N);

page = struct('anchor', zeros(M, 1), 'vector', zeros(M, I), 'number', 0, 'loss', zeros(I, no_iter),...   % number - number of vectors connected to anchor
              'mean', 0, 'corr', zeros(no_iter, 1));                                                     % mean - mean loss
repmat(page, N, 1);

for n = 1 : N
    page(n).loss = zeros(I, no_iter);
    page(n).mean = zeros(no_iter, 1);
    page(n).corr = zeros(no_iter, 1);
end

init = randperm(I);

for n = 1 : N
    page(n).anchor = vec(:, init(n));
end


% -- iterations --
for iter = 1 : no_iter
    
    iter
    for n = 1 : N
        page(n).vector = zeros(M, I);
        page(n).number = 0;
    end

    % -- step 1 (trying to find the best anchor for vectors) --
    for i = 1 : I

        max = 0;
        ind = 0;

        for n = 1 : N

            Frob = norm(page(n).anchor' * vec(:, i), 'fro');
            if Frob >= max
                ind = n;
                max = Frob;
            end
        end

        page(ind).vector(:, page(ind).number + 1) = vec(:, i);
        page(ind).number = page(ind).number + 1;
    end

    % -- step 2 (trying to find the best anchor in a set of vectors) --
    for n = 1 : N

        last = page(n).number;
        if last > 0

            sum = zeros(M, M);
            for i = 1 : last
                sum = sum + page(n).vector(:, i);
            end

            [U, ~, ~]       = svd(sum);
            page(n).anchor  = U(:, 1);

            sum = 0;
            for i = 1 : last
                prod = norm(page(n).anchor' * page(n).vector(:, i), 'fro');
                page(n).loss(i, iter) = 1 - prod;
                sum  = sum + prod;
            end

            page(n).corr(iter) = 1/last * sum;
            page(n).mean(iter) = mean(page(n).loss(1 : last, iter));
        end
    end
end

% --- codebook write ---
for n = 1 : N
    cb(:, n) = page(n).anchor;
end

if exist('graph', 'var') && graph == 1
    lloyd_graph(page, no_iter, N, d, df, rank, Nrank);
end

end



function me_loss = lloyd_graph(page, no_iter, N, d, df, rank, Nrank)

me_loss = zeros(no_iter, 1);   % mean loss

% --- mean calculation ---
for n = 1 : N
    me_loss(:) = me_loss(:) + page(n).mean;
end
me_loss = me_loss / N;

name   = ['functions_stats\stats\mean_loss_' num2str(rank) '(' num2str(Nrank) ')_N=' num2str(N) '_d=' num2str(d) '_df=' num2str(df) '.mat'];
save(name, 'me_loss');

end
















