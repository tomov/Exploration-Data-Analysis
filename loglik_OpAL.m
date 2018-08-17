function loglik = loglik_OpAL(x, data)

latents = OpAL(data, x);

N = length(data.block);

loglik = 0;
for n = 1:N
    if ~data.timeout(n)
        loglik = loglik + latents.loglik(n, data.choice(n));
    end
end
