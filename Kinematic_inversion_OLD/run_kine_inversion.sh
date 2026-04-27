#!/bin/bash


##########################################################
#   SCRIPT TO AUTOMATE THE KINEMATIC INVERSION PROCESS   #
#               USING AN ELLIPTICAL PATCH                #
##########################################################



# Copying pre-processed data to fine_files folder:
cp DATA/real_disp_* Event/kine_files/



# ==============================================
# EDITING "station_in" FILE AND RUN "stations":
# ==============================================
# Editing file:
source_pos=`more input_kine_ellipse.txt | sed '12q;d' | awk '{print $0}'`
nst=`more input_kine_ellipse.txt | sed '55q;d' | awk '{print $0}'`
station_pos=`more input_kine_ellipse.txt | grep "kst" | awk '{print $2, $3, $4, $5}'`

printf "%s\n" "${source_pos}" "${nst}" "${station_pos}"| column -t > station_in
mv station_in Stations/

# Run stations:
echo "-----------------------"
echo "Run stations executable"
echo "-----------------------"
cd Stations/
./stations < station_in
cp station ../Event/
cp station ../AXITRA_F90/


# Update azimuth values for plots of rotated components:
echo "---------------------------------"
echo "Calculating azimuth and distances"
echo "---------------------------------"
cd ../Scripts/
python calc_azi_times.py
cd ../



# ==============================================
# EDITING "faille.in" FILE AND RUN "faille":
# ==============================================
cd Faille/
# Editing file:
flength=`more ../input_kine_ellipse.txt | sed '23q;d' | awk '{print $0}'`
fhypo=`more ../input_kine_ellipse.txt | sed '24q;d' | awk '{print $0}'`
depth_m=`more ../input_kine_ellipse.txt | sed '12q;d' | awk '{print 1000*$3}'`
num_faultss=`more ../input_kine_ellipse.txt | sed '25q;d' | awk '{print $0}'`
strike=`more ../input_kine_ellipse.txt | sed '17q;d' | awk '{print $1}'`
dip=`more ../input_kine_ellipse.txt | sed '17q;d' | awk '{print $2}'`
rake=`more ../input_kine_ellipse.txt | sed '17q;d' | awk '{print $3}'`

sed -i "1s/.*/${flength}/" faille.in
sed -i "2s/.*/${fhypo}/" faille.in
sed -i "3s/.*/${depth_m}/" faille.in
sed -i "4s/.*/${num_faultss}/" faille.in
sed -i "5s/.*/${strike}/" faille.in
sed -i "6s/.*/${dip}/" faille.in
sed -i "7s/.*/${rake}/" faille.in

# Run faille:
echo "---------------------"
echo "Run faille executable"
echo "---------------------"
./faille
cp axi.hist ../Event/
cp source.dat ../Event/
cp axi.hist ../AXITRA_F90/
cp source.dat ../AXITRA_F90/



# ==============================================
# EDITING "kine.in" FILE:
# ==============================================
cd ../Event/kine_files/
# Editing file:
nsource=`more ../../input_kine_ellipse.txt | sed '25q;d' | awk '{print $1*$2}'`
f_param=`echo ${fhypo} ${num_faultss}`
n_ellipse=`more ../../input_kine_ellipse.txt | sed '31q;d' | awk '{print $1}'`
slip_ini=`more ../../input_kine_ellipse.txt | sed '31q;d' | awk '{print $2}'`
slip_shape=`more ../../input_kine_ellipse.txt | sed '31q;d' | awk '{print $3}'`
obs_param=`more ../../input_kine_ellipse.txt | sed '7q;d' | awk '{print $1, $2, $3}'`

sed -i "6s/.*/${nsource}/" kine.in
sed -i "7s/.*/${f_param}/" kine.in
sed -i "8s/.*/${n_ellipse}/" kine.in
sed -i "9s/.*/${slip_ini}/" kine.in
sed -i "10s/.*/${slip_shape}/" kine.in
sed -i "11s/.*/${nst}/" kine.in
sed -i "13s/.*/${obs_param}/" kine.in
sed -i "15s/.*/${obs_param}/" kine.in
sed -i "17s/.*/${obs_param}/" kine.in



# ==============================================
# EDITING "axi.data" FILE:
# ==============================================
cd ../../AXITRA_F90/
# Editing file:
nfreq=`more ../input_kine_ellipse.txt | sed '7q;d' | awk '{print $3/4}'`
tl=`more ../input_kine_ellipse.txt | sed '7q;d' | awk '{print $2-$1}'`  
f1=`more ../input_kine_ellipse.txt | sed '31q;d' | awk '{print $4}'`
f2=`more ../input_kine_ellipse.txt | sed '31q;d' | awk '{print $5}'`

sed -i "3s/.*/ nfreq=${nfreq},/" axi.data
sed -i "4s/.*/ fr1=${f1},/" axi.data
sed -i "5s/.*/ fr2=${f2},/" axi.data
sed -i "6s/.*/ tl=${tl},/" axi.data
sed -i "8s/.*/ nr=${nst},/" axi.data
sed -i "9s/.*/ ns=${nsource},/" axi.data



# ==============================================
# RUN "axitra.out" and edit "axi.head" FILE:
# ==============================================
# Run excecutable:
echo "-------------------------"
echo "Run axitra.out executable"
echo "-------------------------"
./axitra.out

# Editing axi.head in main folder:
xl=`more axi.data | grep "xl" | sed 's/\=/ /g' | awk '{print $2}'`
ikmax=`more axi.data | grep "ikmax" | sed 's/\=/ /g' | awk '{print $2}'`
t0=`more ../input_kine_ellipse.txt | sed '31q;d' | awk '{print $6}'`
nc=`more axi.data | grep "nc" | sed 's/\=/ /g' | awk '{print $2}'`

# Using axi.head2 format:
sed -i "2s/.*/ NC      =          ${nc}/" axi.head
sed -i "3s/.*/ NFREQ   =          ${nfreq},/" axi.head
sed -i "4s/.*/ TL      =          ${tl},/" axi.head
sed -i "5i\ AW      =   2.00000000000000     ," axi.head
sed -i "6s/.*/ NR      =          ${nst},/" axi.head
sed -i "7s/.*/ NS      =          ${nsource},/" axi.head
sed -i "8s/.*/ XL      =         ${xl}/" axi.head
sed -i "9s/.*/ T0      =         ${t0},/" axi.head
sed -i "10s/.*/ FR1     =          ${f1},/" axi.head
sed -i "11s/.*/ FR2     =          ${f2},/" axi.head
sed -i "12s/.*/ IKMAX   =         ${ikmax}/" axi.head
sed -i "13s/.*/ UCONV   =  1.000000000000000E-003,/" axi.head

# Moving output files to earthquake folder:
mv axi.res ../Event/
cp axi.head ../Event/



# ==============================================
# Edit "kine_param" AND "na.in" FILES:
# ==============================================
cd ../Event/kine_files/
# Editing kine_param:
npar=`more ../../input_kine_ellipse.txt | sed '37q;d' | awk '{print $0}'`
p1=`more ../../input_kine_ellipse.txt | sed '38q;d' | awk '{print $0}'`
p2=`more ../../input_kine_ellipse.txt | sed '39q;d' | awk '{print $0}'`
p3=`more ../../input_kine_ellipse.txt | sed '40q;d' | awk '{print $0}'`
p4=`more ../../input_kine_ellipse.txt | sed '41q;d' | awk '{print $0}'`
p5=`more ../../input_kine_ellipse.txt | sed '42q;d' | awk '{print $0}'`
p6=`more ../../input_kine_ellipse.txt | sed '43q;d' | awk '{print $0}'`
p7=`more ../../input_kine_ellipse.txt | sed '44q;d' | awk '{print $0}'`

sed -i "1s/.*/ ${npar}/" kine_param
sed -i "2s/.*/   ${p1}/" kine_param
sed -i "3s/.*/   ${p2}/" kine_param
sed -i "4s/.*/   ${p3}/" kine_param
sed -i "5s/.*/   ${p4}/" kine_param
sed -i "6s/.*/   ${p5}/" kine_param
sed -i "7s/.*/   ${p6}/" kine_param
sed -i "8s/.*/   ${p7}/" kine_param

# Editing na.in:
cd ../
alg_type=`more ../input_kine_ellipse.txt | sed '50q;d' | awk '{print $1}'`
iter=`more ../input_kine_ellipse.txt | sed '50q;d' | awk '{print $2}'`
ss1=`more ../input_kine_ellipse.txt | sed '50q;d' | awk '{print $3}'`
ss_other=`more ../input_kine_ellipse.txt | sed '50q;d' | awk '{print $4}'`
n_cell=`more ../input_kine_ellipse.txt | sed '50q;d' | awk '{print $5}'`

sed -i "4s/.* :/${alg_type}        :/" na.in
sed -i "5s/.* :/${iter}        :/" na.in
sed -i "6s/.* :/${ss1}        :/" na.in
sed -i "7s/.* :/${ss_other}        :/" na.in
sed -i "8s/.* :/${n_cell}        :/" na.in



# ==============================================
# EXECUTING INVERSION AND FORWARD METHOD:
# ==============================================
# Inverse method:
echo "-----------------------"
echo "Run Inverse Method"
echo "-----------------------"
../bin/kine_na

# Direct method:
echo "-----------------------"
echo "Run Direct Method"
echo "-----------------------"
../bin/kine_direct

echo "PROCESSING FINISHED !!"
