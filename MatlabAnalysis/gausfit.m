function y = gausfit(params,x)

y = params(1) * exp(-1*(x - params(2)).^2 ./ (2*params(3)^2));

return;