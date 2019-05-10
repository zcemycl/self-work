function Spiral(n)
s = flipud(spiral(n)); 
t = isprime(s);
im2bw(t);
imshow(t);
axis('equal')
axis('off')


end 