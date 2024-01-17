function grids = makegrids(bmin, bmax, I, bgrid_type, amin, amax, J, agrid_type,z)

switch bgrid_type
    case 'Linear'
        grids.b = linspace(bmin,bmax,I)';
    case 'NL_symmetric'
        grids.b = bmin+(1-real(exp((0:I-1)'*1i*pi/(I-1))))*(bmax-bmin)/2;
    case 'NL_lefts'
        grids.b = bmin+(1-real(exp((0:I-1)'*1i*pi/2/(I-1))))*(bmax-bmin);
    case 'NL_rights'
        grids.b = bmin+flip(real(exp((0:I-1)'*1i*pi/2/(I-1))))*(bmax-bmin);
    otherwise
        grids.b = linspace(bmin,bmax,I)';
end

switch agrid_type
    case 'Linear'
        grids.a = linspace(amin,amax,J);
    case 'NL_symmetric'
        grids.a = amin+(1-real(exp((0:J-1)*1i*pi/(J-1))))*(amax-amin)/2;
    case 'NL_lefts'
        grids.a = amin+(1-real(exp((0:J-1)*1i*pi/2/(J-1))))*(amax-amin);
    case 'NL_rights'
        grids.a = amin+flip(real(exp((0:J-1)*1i*pi/2/(J-1))))*(amax-amin);
    otherwise
        grids.a = linspace(amin,amax,J);
end

[grids.aaa, grids.bbb, grids.zzz] = meshgrid(grids.a,grids.b,z);

end