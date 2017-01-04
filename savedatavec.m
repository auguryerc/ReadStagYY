function [DAT] = savedatavec(DATAX,DATAY,thetad,phid)

    [THD PHD] = meshgrid(thetad,phid);
    ntheta    = length(thetad);
    nphi      = length(phid);

    DAT  = zeros(ntheta*nphi,3);
    for ith = 1:ntheta
        for iph = 1:nphi                   
            ilg         = (ith-1)*nphi + iph;   
            DAT(ilg,1)  = PHD(iph,ith);
            DAT(ilg,2)  = THD(iph,ith);
            DAT(ilg,3)  = DATAX(ith,iph);
            DAT(ilg,4)  = DATAY(ith,iph);
        end
    end