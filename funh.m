function sigmoid_out=funh(x)
%sigmoid_out=1./(1+exp(-x));
if x<0
  sigmoid_out=0;
else
  if x>1
    sigmoid_out=1;
    else
      sigmoid_out=x;
  end
end
