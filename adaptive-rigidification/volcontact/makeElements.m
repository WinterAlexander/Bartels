function [el] = makeElements( p, t )
% MAKEELEMENTS Creates elements from triangles.
%   makeElements( p, t ) returns an array of structures given 2D
%   points p (stored as rows), and triangles t defined by indices
%   (1 indexed and stored as rows).

el = repmat(struct('t',1), size(t,1), 1 );
for index = 1:size(t,1)
    i = t(index,1);
    j = t(index,2);
    k = t(index,3);
    Dm = [ p(i,:)' - p(k,:)', p(j,:)' - p(k,:)' ];
    el(index).t = t(index,:);
    el(index).Bm = inv( Dm );
    area = 1/2 * det( Dm );
    el(index).area = area;
end