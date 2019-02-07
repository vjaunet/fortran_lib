%
%
%    Reading screech data
%     - the first 6 sensors are on an azimuthal array
%     - the seventh is a farfield microphone in the nozzle plane
%
%
%------------------------------------------------------------------

clear all;

%read whole plane PDATA data
fid=fopen('mach110_screech_pa.bin','r','ieee-le');
pdata.version = char(fread(fid,4,'uchar'));
%if (pdata.version ~= "v0.1") break;

%read the header
pdata.nsensors   =fread(fid,1,'int');
pdata.nsamples   =fread(fid,1,'int');
pdata.fs         =fread(fid,1,'float'); %smapling freq
pdata.Pa         =fread(fid,1,'float'); %atmospheric pressure

%read comments
pdata.comments =char(fread(fid,500,'char'));


%read the pressure and senors localisation data
pdata.p =single(zeros(pdata.nsensors,pdata.nsamples));
pdata.x =single(zeros(pdata.nsensors,2));

pdata.x = fread(fid,[pdata.nsensors 2],'float');
pdata.p = fread(fid,[pdata.nsamples pdata.nsensors],'float');
fclose(fid);
