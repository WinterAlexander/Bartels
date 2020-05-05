function x = fcg(A, b, x)
% fcg Unfinished function for computing filtered conjugate gradients

    r = b - A * x;
    % filter rf
    p = r;
    rsold = r' * r;

    for i = 1:length(b)
        Ap = A * p;
        % filter Ap
        alpha = rsold / (p' * Ap);
        x = x + alpha * p;
        r = r - alpha * Ap;
        rsnew = r' * r;
        if sqrt(rsnew) < 1e-10
              break;
        end
        p = r + (rsnew / rsold) * p;
        % filter p
        rsold = rsnew;
    end
end