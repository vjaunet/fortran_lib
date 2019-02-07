%
%
%--------------------------------------------------------------

%read whole plane PIV data
fid=fopen('sys1_polar.bin','r','ieee-le');
piv.typeofgrid = char(fread(fid,1,'uchar'));

%read the header
piv.nr         =fread(fid,1,'int');
piv.ntheta     =fread(fid,1,'int');
piv.ncomponents=fread(fid,1,'int');
piv.nsamples   =fread(fid,1,'int');
piv.dx         =fread(fid,1,'float');
piv.dy         =fread(fid,1,'float');
piv.x0         =fread(fid,1,'float');
piv.y0         =fread(fid,1,'float');
piv.pix_step   =fread(fid,1,'int');
piv.fs         =fread(fid,1,'float');
piv.z_pos      =fread(fid,1,'float');

%read stagnation conditions
piv.ncgen      =fread(fid,1,'int');
if (piv.ncgen > 0)
    piv.cgen = fread(fid,piv.ncgen,'float');
end

%read comments
piv.comments =char(fread(fid,500,'char'));


%read the velocity data
piv.u =single(zeros(piv.nr,piv.ntheta,piv.ncomponents, ...
                             piv.nsamples));
for is=1:piv.nsamples
    for ic=1:piv.ncomponents
        piv.u(:,:,ic,is) =fread(fid,[piv.nr piv.ntheta],'float');
    end
end
fclose(fid);

%build the spatial coordinates
piv.r=single(zeros(piv.nr,1));
piv.theta=single(zeros(piv.ntheta,1));
for ir=1:piv.nr
    for it=1:piv.ntheta
        piv.r(ir) =  (ir-1) * piv.dx * cos((it-1)*piv.dy);
        piv.theta(it) = (ir-1) * piv.dx * cos((it-1)*piv.dy);
    end
end

%plot a PIV field
pcolor(piv.u(:,:,1,1));
