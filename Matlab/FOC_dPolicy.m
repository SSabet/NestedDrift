function d = FOC_dPolicy(Va,MUc,a, chi0, chi1)
    d = min(Va./MUc - 1 + chi0,0).*a/chi1 +  max(Va./MUc - 1 - chi0,0).*a/chi1;     
end