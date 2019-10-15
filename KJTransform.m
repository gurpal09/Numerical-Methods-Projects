function [x,y] = KJTransform(c, Xi, Eta)

x = Xi .* (1 + c^2 ./ (Eta.^2 + Xi.^2));
y = Eta .* (1 - c^2 ./ (Eta.^2 + Xi.^2));

end
