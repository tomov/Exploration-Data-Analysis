function latents = AU(data, x)
    
    % Actor with uncertainty (AU) model
    %
    % USAGE: latents = AU(data, x)
    %
    % INPUTS:
    %   data - single subject data
    %   x (optional) - parameters
  
    if nargin < 2
        G_0R = 0.1;
        N_0R = 0.1;
        G_0S = 0.1;
        N_0S = 0.1;
        alpha = 0.1;
        beta = 0.1;
        a = 2;
        b = 2;
    else
        G_0R = x(1);
        N_0R = x(2);
        G_0S = x(3);
        N_0S = x(4);
        alpha = x(5);
        beta = x(6);
        a = x(7);
        b = x(8);
    end

    N = length(data.block);
    
    for n = 1:N
        
        % initialization at the start of each block
        if n == 1 || data.block(n)~=data.block(n-1)
            switch data.cond(n)
                case 1 % RS
                    G = [G_0R G_0S];
                    N = [N_0R N_0S];
                case 2 % SR
                    G = [G_0S G_0R];
                    N = [N_0S N_0R];
                case 3 % RR 
                    G = [G_0R G_0R];
                    N = [N_0R N_0R];
                case 4 % SS
                    G = [G_0S G_0S];
                    N = [N_0S N_0S];
            end
        end

        loglik = a * G - b * N - logsumexp(a * G - b * N); % avoid numeric underflow

        c = data.choice(n);
        r = data.reward(n);

        % store latents
        latents.G(n,:) = G;
        latents.N(n,:) = N;
        latents.loglik(n,:) = loglik;
        
        % update
        Q = G(c) - N(c);
        G(c) = G(c) + alpha * rect_pos(r - Q) - beta * G(c);
        N(c) = N(c) + alpha * rect_neg(r - Q) - beta * N(c);
    end
end


function x = rect_pos(x)
    if x <= 0
        x = 0;
    end
end

function x = rect_neg(x)
    if x >= 0
        x = 0;
    end
end
