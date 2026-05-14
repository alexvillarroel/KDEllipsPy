#!/bin/bash


##########################################################
#    SCRIPT TO AUTOMATE THE DYNAMIC INVERSION PROCESS    #
#               USING A ELLIPTICAL PATCH                 #
##########################################################

#
# ===========================
# AUTHOR: Carlos Herrera R.
# ===========================
#

# Name of the event (event processing folder):
ev_name=Evento



# ==============================================
# CHANGING PARAMETERS OF PLOTTING SCRIPTS:
# ==============================================
# Removing previous results:
if [ -f ${ev_name}/momentos0.dat ]; then
   rm ${ev_name}/momentos*
fi

# Getting parameters:
ncores=`more input_dynamic_ellipse.txt | sed '36q;d' | awk '{print $5}'`



# ==============================================
# EDITING "station_in" FILE AND RUN "stations":
# ==============================================
# Editing file:
source_pos=`more input_dynamic_ellipse.txt | sed '12q;d' | awk '{print $0}'`
station_pos=`more input_dynamic_ellipse.txt | grep "dst" | awk '{print $2, $3, $4, $5}'`
nst=`more input_dynamic_ellipse.txt | sed '57q;d' | awk '{print $1}'`

printf "%s\n" "${source_pos}" "${nst}" "${station_pos}"| column -t > station_in
mv station_in Stations/

# Run stations:
echo "-----------------------"
echo "Run stations executable"
echo "-----------------------"
cd Stations/
./stations < station_in
cp station ../${ev_name}/
cp station ../AXITRA_F90/



# ==============================================
# EDITING "faille.in" FILE AND RUN "faille":
# ==============================================
cd ../Faille/
# Editing file:
flength=`more ../input_dynamic_ellipse.txt | sed '23q;d' | awk '{print $0}'`
fhypo=`more ../input_dynamic_ellipse.txt | sed '24q;d' | awk '{print $0}'`
depth_m=`more ../input_dynamic_ellipse.txt | sed '12q;d' | awk '{print 1000*$3}'`
num_faultss=`more ../input_dynamic_ellipse.txt | sed '25q;d' | awk '{print $0}'`
strike=`more ../input_dynamic_ellipse.txt | sed '17q;d' | awk '{print $1}'`
dip=`more ../input_dynamic_ellipse.txt | sed '17q;d' | awk '{print $2}'`
rake=`more ../input_dynamic_ellipse.txt | sed '17q;d' | awk '{print $3}'`

sed -i "1s/.*/${flength}/" faille.in
sed -i "2s/.*/${fhypo}/" faille.in
sed -i "3s/.*/${depth_m}/" faille.in
sed -i "4s/.*/${num_faultss}/" faille.in
sed -i "5s/.*/${strike}/" faille.in
sed -i "6s/.*/${dip}/" faille.in
sed -i "7s/.*/${rake}/" faille.in

# Run faille:
echo "-----------------------"
echo "Run faille executable"
echo "-----------------------"
./faille
cp axi.hist ../${ev_name}/
cp source ../${ev_name}/
cp axi.hist ../AXITRA_F90/
cp source ../AXITRA_F90/




# ==============================================
# EDITING "axi.data" FILE AND RUN "axitra.out":
# ==============================================
cd ../AXITRA_F90/

# Editing file:
nfreq=`more ../input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $3/4}'`
f1=`more ../input_dynamic_ellipse.txt | sed '30q;d' | awk '{print $1}'`
f2=`more ../input_dynamic_ellipse.txt | sed '30q;d' | awk '{print $2}'`
tl=`more ../input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $2-$1}'`
ns=`more ../input_dynamic_ellipse.txt | sed '25q;d' | awk '{print $1*$2}'`

sed -i "3s/.*/ NFREQ   =          ${nfreq},/" axi.data
sed -i "4s/.*/ Fr1     =          ${f1},/" axi.data
sed -i "5s/.*/ Fr2     =          ${f2},/" axi.data
sed -i "6s/.*/ TL      =          ${tl},/" axi.data
sed -i "8s/.*/ NR      =          ${nst},/" axi.data
sed -i "9s/.*/ NS      =          ${ns},/" axi.data

# Run Axitra excecutable:
echo "-------------------------"
echo "Run axitra.out executable"
echo "-------------------------"
./axitra.out

# Moving useful files to the earthquake folder:
mv axi.res ../${ev_name}/
cp axi.head ../${ev_name}/



# ==============================================
# EDITING "fd3d.in", "na.in" and "fd3d_param":
# ==============================================
# Editing fd3d.in file:
cd ../${ev_name}/
npts=`more ../input_dynamic_ellipse.txt | sed '7q;d' | awk '{print $3}'`

sed -i "5s/.*/${npts}/" fd3d.in
sed -i "6s/.*/${nst}/" fd3d.in


# Editing na.in:
alg_type=`more ../input_dynamic_ellipse.txt | sed '36q;d' | awk '{print $1}'`
num_iter=`more ../input_dynamic_ellipse.txt | sed '36q;d' | awk '{print $2}'`
samp_size=`more ../input_dynamic_ellipse.txt | sed '36q;d' | awk '{print $3}'`
cells=`more ../input_dynamic_ellipse.txt | sed '36q;d' | awk '{print $4}'`

sed -i "4s/.* :/${alg_type}        :/" na.in
sed -i "5s/.* :/${num_iter}        :/" na.in
sed -i "6s/.* :/${samp_size}        :/" na.in
sed -i "7s/.* :/${samp_size}        :/" na.in
sed -i "8s/.* :/${cells}        :/" na.in


# Editing fd3d_param file with the intervals to search in each inverted parameter:
npar=`more ../input_dynamic_ellipse.txt | sed '42q;d' | awk '{print $0}'`
p1=`more ../input_dynamic_ellipse.txt | sed '43q;d' | awk '{print $0}'`
p2=`more ../input_dynamic_ellipse.txt | sed '44q;d' | awk '{print $0}'`
p3=`more ../input_dynamic_ellipse.txt | sed '45q;d' | awk '{print $0}'`
p4=`more ../input_dynamic_ellipse.txt | sed '46q;d' | awk '{print $0}'`
p5=`more ../input_dynamic_ellipse.txt | sed '47q;d' | awk '{print $0}'`
p6=`more ../input_dynamic_ellipse.txt | sed '48q;d' | awk '{print $0}'`
p7=`more ../input_dynamic_ellipse.txt | sed '49q;d' | awk '{print $0}'`
p8=`more ../input_dynamic_ellipse.txt | sed '50q;d' | awk '{print $0}'`
p9=`more ../input_dynamic_ellipse.txt | sed '51q;d' | awk '{print $0}'`
p10=`more ../input_dynamic_ellipse.txt | sed '52q;d' | awk '{print $0}'`

sed -i "1s/.*/   ${npar}/" fd3d_param
sed -i "2s/.*/   ${p1}/" fd3d_param
sed -i "3s/.*/   ${p2}/" fd3d_param
sed -i "4s/.*/   ${p3}/" fd3d_param
sed -i "5s/.*/   ${p4}/" fd3d_param
sed -i "6s/.*/   ${p5}/" fd3d_param
sed -i "7s/.*/   ${p6}/" fd3d_param
sed -i "8s/.*/   ${p7}/" fd3d_param
sed -i "9s/.*/   ${p8}/" fd3d_param
sed -i "10s/.*/   ${p9}/" fd3d_param
sed -i "11s/.*/   ${p10}/" fd3d_param



# ==============================================
# EXECUTING INVERSION AND FORWARD METHOD:
# ==============================================
# Inverse method (editing the numbers of cores):
sed -i "12s/.*/cores=${ncores}/" fd3d_na.run

# Run inversion:
echo "----------------------------------"
echo "Run inversion with ${ncores} cores"
echo "----------------------------------"
./fd3d_na.run


# Run forward method:
echo "-----------------------"
echo "Run Direct Method"
echo "-----------------------"
./fd3d_direct.run

# Show a reminder about the rotation of components in the calculation of misfit:
rot=`grep "r_obs" ../Source/fd3d_subs/calcmisfit.f -r | head -n 1 | awk '{print $1}'`
if [ -z ${rot} ]; then
   echo "=============================================================="
   echo "REMINDER: You are modeling traces in the N-E horizontal system"
   echo "=============================================================="
else
   echo "=============================================================="
   echo "REMINDER: You are modeling traces in the R-T horizontal system"
   echo "=============================================================="
fi


