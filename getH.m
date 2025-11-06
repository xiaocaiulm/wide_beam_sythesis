function G = getH(y,m,n)
    G = zeros(2*m,2*n);
    G(1:m,1:n)=real(y);
    G(m+1:2*m, n+1:2*n) = -real(y);
    G(1:m,n+1:2*n) = imag(y);
    G(m+1:2*m,1:n) = imag(y);
end