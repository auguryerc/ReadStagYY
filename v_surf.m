clear all
close all

%--------- General information------------------ 

FileFormat = 'n';

number_start = 1;
incr         = 1;

%tw=number_start:number_end;             %imagefile

netalaw   = 1;         %number of viscosity laws used
    
dlim    = 0;  %depth of transition
part    = 0;  % 0=bottom, 1=surface
   
for etalaw=1:netalaw
if (etalaw==1)          
number_end =343;
directory='/Users/coltice/Desktop/Work/Projets/STAGVIZ';
refname ='Longrun';
end

%loop on steps of pertubation
for step_number=number_start:incr:number_end 
        
  step_number
 
   %% if excluded cells : vol and redflag
    if dlim~=0
        Type='velocity'; 
        [Vol,~,~,~,~,~,~,RedFlag,~,~,~,~,~,~,z,vs,v,~,~,~]=...
            ReadStag3D_LeaV(directory,refname,step_number,Type);
        zlim = 1-dlim ; %value of z coordinate of last/first excluded cells
        if (step_number==number_start)
            RedFlagdeep=RedFlag;
            if(part==1) % 1=upper mantle
                iz=1;
                while z(iz)<=(zlim)  % for layered visco width of transition zone=2/48
                    RedFlagdeep(:,:,iz)=1;
                   iz=iz+1;
                end
            else           % 0=lower mantle
                iz=size(z);
                while z(iz)>=(zlim)
                    RedFlagdeep(:,:,iz)=1;
                    iz=iz-1;
                end
            end
        end
    else
        Type='velocity'; 
        [Vol,~,~,~,~,~,~,RedFlag,~,~,~,~,~,~,~,vs,v,~,~,~]=...
            ReadStag3D_LeaV(directory,refname,step_number,Type);
        if (step_number==number_start)
        RedFlagdeep=RedFlag;
        end
    end
 
    nredflag=size(v,1);
    
    Vol=Vol(~logical(RedFlagdeep));
    Vtot=sum(Vol);
        
    %% velocity stuff  
    v_lin=v(~logical(RedFlagdeep)).*(Vol);
    vpart(step_number-number_start+1)=(sum(v_lin)/Vtot)^(1/2);
    vmax(step_number-number_start+1)=max(v(~logical(RedFlagdeep)))^(1/2);
    vsurf(step_number-number_start+1)=vs;
        
    %dataref=load(sprintf('%s%s_time.dat',directory,refname)); % reading vrms in time.dat file for comparison 
    %vrms(step_number)=dataref(step_number,9);
                
end
  
    %vmoy(etalaw)=mean(vrms);
    vsurfmoy(etalaw)=mean(vsurf);
    vpartmoy(etalaw)=mean(vpart);

  vsurf
   clear vsurf
   %clear vrms
   clear dataref
   
vsurfmoy,%vmoy

end
  
%% saving data
%vsurfmoy,vmoy,vpartmoy
%Mforsave=[vsurfmoy,vmoy]
%save('vsurf.txt','Mforsave','-ascii');


