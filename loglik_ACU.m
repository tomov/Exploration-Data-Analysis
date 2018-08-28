function loglik = loglik_ACU(x, data)

latents = ACU(data, x);

N = length(data.block);

loglik = 0;
for n = 1:N
    if ~data.timeout(n)
        loglik = loglik + latents.loglik(n, data.choice(n));
    end
end
