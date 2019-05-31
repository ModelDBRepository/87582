function y=filtras250(f,Q,x)

a=-pi*f/Q;
b=sqrt((2*pi*f)^2-a*a);
t=1:251;
h=(1/b)*exp(a*t).*sin(b*t);

y=conv(h,x);
