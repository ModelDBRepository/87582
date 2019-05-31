function y=histereze(x,delta_x)


if (x<0.5)
  if delta_x<0
    y=-x+1./(1+exp(-delta_x)*(1-x)./x);
  else
    y=0.25*delta_x;
  end
else
  if delta_x<0
    y=0.25*delta_x;
  else
    y=-x+1./(1+exp(-delta_x)*(1-x)./x);
  end
end
