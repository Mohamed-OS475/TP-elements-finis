function [tilde_AA, tilde_LL] = elimine(AA, LL, Refneu)

tilde_AA = AA;
tilde_LL = LL;
n = size(LL,1);
for i=1:n
    if (Refneu(i) == 1)
        tilde_LL(i) = 0;
        tilde_AA(i,:) = 0;
        tilde_AA(:,i) = 0;
        tilde_AA(i, i) = 1;
    end
end
end
