#!/bin/bash



################################################################
#   SCRIPT TO COMPILE SOURCE CODE OF THE KINEMATIC INVERSION   #
#                AND ALL ITS RELATED CODES                     #
################################################################



# ===================================================
# DEFINE AND RETRIEVE PARAMETERS
# ===================================================
# Parameters:
dt=`more input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $4}'`
ntp=`more input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $3}'`
fr1=`more input_dynamic_ellipse.txt | sed '30q;d' | awk '{print $1}'`
fr2=`more input_dynamic_ellipse.txt | sed '30q;d' | awk '{print $2}'`
units=`more input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $5}'`
nsp=5000
ncp=16
nrp=8



# ===================================================
# EDIT CRITICAL SCRIPTS
# ===================================================
# Edit "dimension.inc" file:
cd Source/fd3d_subs/
sed -i "4s/.*/	parameter   (dt=${dt}, nsp=${nsp},ncp=${ncp},nrp=${nrp},ntp=${ntp},/" dimension.inc 
if  [ ${units} -eq 1 ]; then
    # If units are in displacement, use order 4:
    sed -i "5s/.*/     *  nrtp=nrp*ntp, fr2=${fr2},fr1=${fr1}, bp_order=4)/" dimension.inc

elif [ ${units} -eq 2 ]; then
    # If units are in velocity, use order 2:
    sed -i "5s/.*/     *  nrtp=nrp*ntp, fr2=${fr2},fr1=${fr1}, bp_order=2)/" dimension.inc
fi


# Edit filters:
if  [ ${units} -eq 1 ]; then
    # If units are in displacement, use order 4:
    sed -i "383s/.*/         call XAPIIR(taz,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,dt,1)/" convm_dynsource.f
    sed -i "384s/.*/         call XAPIIR(tay,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,dt,1)/" convm_dynsource.f
    sed -i "385s/.*/         call XAPIIR(tax,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,dt,1)/" convm_dynsource.f

elif [ ${units} -eq 2 ]; then
    # If units are in velocity, use order 2:
    sed -i "383s/.*/         call XAPIIR(taz,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)/" convm_dynsource.f
    sed -i "384s/.*/         call XAPIIR(tay,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)/" convm_dynsource.f
    sed -i "385s/.*/         call XAPIIR(tax,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,dt,1)/" convm_dynsource.f
fi


# Edit integrations according to the units:
if [ ${units} -eq 1 ]; then
    sed -i "315s/.*/          do ir=1,nr/" convm_dynsource.f
    sed -i "316s/.*/               ux(jf,ir)=ux(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "317s/.*/               uy(jf,ir)=uy(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "318s/.*/               uz(jf,ir)=uz(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "319s/.*/               if(abs(ux(jf,ir)).gt.amax)amax=abs(ux(jf,ir))/" convm_dynsource.f
    sed -i "320s/.*/          enddo/" convm_dynsource.f
elif [ ${units} -eq 2 ]; then
    sed -i "315s/.*/c          do ir=1,nr/" convm_dynsource.f
    sed -i "316s/.*/c               ux(jf,ir)=ux(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "317s/.*/c               uy(jf,ir)=uy(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "318s/.*/c               uz(jf,ir)=uz(jf,ir)\/(ai * omega)/" convm_dynsource.f
    sed -i "319s/.*/c               if(abs(ux(jf,ir)).gt.amax)amax=abs(ux(jf,ir))/" convm_dynsource.f
    sed -i "320s/.*/c          enddo/" convm_dynsource.f
else
    echo "ERROR: You must specify 1 or 2 for units variable"
    echo "       when you are defining displacement or velocity"
    exit
fi



# ===================================================
# COMPILE STAGE
# ===================================================
# Compile source code of the dynamic inversion:
echo "---------------------------------------------"
echo "Compiling source code of dynamic inversion:"
echo "---------------------------------------------"
cd ../
make clean
make all
make fd3d_direct
make fd3d_na_mpi

# Copying main executables to primary folder:
cp fd3d_direct ../
cp fd3d_na_mpi ../
