function adj_p = f_Bonferroni(p)
adj_p = p*length(p(:));
adj_p(adj_p>1) = 1;
end
