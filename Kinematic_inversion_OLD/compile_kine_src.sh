#!/bin/bash



################################################################
#   SCRIPT TO COMPILE SOURCE CODE OF THE KINEMATIC INVERSION   #
#                AND ALL ITS RELATED CODES                     #
################################################################



# ===================================================
# DEFINE AND RETRIEVE PARAMETERS
# ===================================================
# Main parameters:
nsp=1000
ncp=40
nrp=40
ntp=`more input_kine_ellipse.txt | sed '7q;d' | awk '{print $3}'`
dt=`more input_kine_ellipse.txt | sed '7q;d' | awk '{print $4}'`
units=`more input_kine_ellipse.txt | sed '7q;d' | awk '{print $5}'`



# ===================================================
# EDIT CRITICAL SCRIPTS
# ===================================================
# Edit "kine_param.inc" file:
cd src-covm_inkm/
sed -i "11s/.*/	parameter  (dt=${dt},nsp=${nsp},ncp=${ncp},nrp=${nrp},ntp=${ntp})/" kine_param.inc

# Edit filters and physical units:
cd kine_subs/
if [ ${units} -eq 1 ]; then
    # If units are in displacement, use order 4:
    sed -i "127s/.*/      icc=1/" convm_new_src.f
    sed -i "285s/.*/         call XAPIIR(taz,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,pas,1)/" convm_new_src.f
    sed -i "286s/.*/         call XAPIIR(tay,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,pas,1)/" convm_new_src.f
    sed -i "287s/.*/         call XAPIIR(tax,ntp,'BU',0.0,0.0,4,'BP',fr1,fr2,pas,1)/" convm_new_src.f

elif [ ${units} -eq 2 ]; then
    # If units are in velocity, use order 2:
    sed -i "127s/.*/      icc=2/" convm_new_src.f
    sed -i "285s/.*/         call XAPIIR(taz,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)/" convm_new_src.f
    sed -i "286s/.*/         call XAPIIR(tay,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)/" convm_new_src.f
    sed -i "287s/.*/         call XAPIIR(tax,ntp,'BU',0.0,0.0,2,'BP',fr1,fr2,pas,1)/" convm_new_src.f
fi



# ===================================================
# COMPILE STAGE
# ===================================================
# Compile source code of the kinematic inversion:
echo "---------------------------------------------"
echo "Compiling source code of kinematic inversion:"
echo "---------------------------------------------"
cd ../
make clean
make all

# Compile supplementary codes:
echo ""
echo "------------------------------"
echo "Compiling supplementary codes:"
echo "------------------------------"
cd ../Stations/
ifx stations.f -o stations
cd ../Faille/
ifx faille.f -o faille
cd ../AXITRA_F90/
make
cd ../
