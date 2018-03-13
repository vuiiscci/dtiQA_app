function bing = binghamDist(e1, e2, e3, k1, k2)
%
% Gets a Bingham distribution object.
%
% Args:
%  e1, e2, e3 - eigenvectors of distribution, sorted by magnitude, largest first (remembering k1 < k2 < k3 < 0)  
%  
%  k1, k2 - concentration parameters, k1 < k2 < 0.
%
% If k1 == k2, distribution is Watson with concentration |k| and mean axis e3
% If k2 == 0, distribution is Watson with concentration k1 and mean axis e1
%
% example:
%
%  bing = binghamDist([1 0 0], [0 1 0], [0 0 1], -100, -50);
%  vec = bing.nextVector();
%  matrix = [vec.x, vec.y, vec.z];
%
%  aVec = numerics.Vector3D(0.5, 0.866025, 0);
%  p = bing.pdf(aVec)
%  

evecs = javaArray('numerics.Vector3D', 3);

evecs(1) = numerics.Vector3D(e1(1), e1(2), e1(3));
evecs(2) = numerics.Vector3D(e2(1), e2(2), e2(3));
evecs(3) = numerics.Vector3D(e3(1), e3(2), e3(3));

bing = javaMethod('getBinghamDistribution', 'numerics.BinghamDistribution', evecs, k1, k2, numerics.MTRandom(839462));

