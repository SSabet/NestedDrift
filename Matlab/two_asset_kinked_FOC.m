function d = two_asset_kinked_FOC(V_a,MUc,a, chi0, chi1)
d = min(V_a./MUc - 1 + chi0,0).*a/chi1 +  max(V_a./MUc - 1 - chi0,0).*a/chi1;