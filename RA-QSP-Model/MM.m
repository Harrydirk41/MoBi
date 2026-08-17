function [ out ] = MM( X, Vm, Km, n )

out = Vm * X.^n ./(Km.^n + X.^n);

end

