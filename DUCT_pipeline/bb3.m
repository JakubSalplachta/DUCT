function bb = bb3(mask)
% 3D bounding box
%   Computing the bounding box of non-zero voxels in the input 3D volume
[a,b,c]=size(mask);

%Bounding box
minZ = c;
maxZ = 1;
for z = 1:c
    s = sum(sum(mask(:,:,z)));
    if s > 0
        minZ = z;
        break;
    end
end
for z = c:-1:minZ
    s = sum(sum(mask(:,:,z)));
    if s > 0
        maxZ = z;
        break;
    end
end
minY = a;
maxY = 1;
for y = 1:a
    s = sum(sum(mask(y,:,:)));
    if s > 0
        minY = y;
        break;
    end
end
for y = a:-1:minY
    s = sum(sum(mask(y,:,:)));
    if s > 0
        maxY = y;
        break;
    end
end
minX = b;
maxX = 1;
for x = 1:b
    s = sum(sum(mask(:,x,:)));
    if s > 0
        minX = x;
        break;
    end
end
for x = b:-1:minX
    s = sum(sum(mask(:,x,:)));
    if s > 0
        maxX = x;
        break;
    end
end

bb = [minY, maxY, minX, maxX, minZ, maxZ];

end

