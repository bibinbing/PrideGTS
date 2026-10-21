function scatnxm(a,b);
% Scatter plots of the columns of a against b
%
% Code by Stephen Roberts
  

  if (length(a) ~= length(b))
    error('a and b must be the same length!');
  end;
  
  if (size(a,1)<size(a,2))
    a = a';
  end;
  if (size(b,1)<size(b,2))
    b = b';
  end;
  
  N = min(size(a)); M = min(size(b));
  if (N>M)
    c=b; b=a; a=c;
    c=M; M=N; N=c;
  end;
  for n=1:N,
    for m=1:M,
      subplot(N,M,(n-1)*M+m),plot(a(:,n),b(:,m),'.');
      axis('off');
    end;
  end;
  
  return;