function [b,a] = compute(b1,a1,b2,a2)

b = conv(b1,a2) - conv(b2,a1);
a = 2*conv(a2,a1);

end