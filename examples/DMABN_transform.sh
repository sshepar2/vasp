#####################################################################################
# bash script used to call on transform.py multiple times with other commands       #
#                                                                                   #
# Execute in terminal prompt as:                                                    #
#                                                                                   #
# $ ./DMABN_rot_stretch.sh                                                          #
#                                                                                   #
# Uses the POSCAR file in this directory.                                           #
#                                                                                   #
# 1. Use rotate functionality to rotate the dimethylamino end of a DMABN            #
#    molecule by 90 degrees out of the molecular plane.                             #
#                                                                                   #
# 2. Use the stretch functionality to make the c axis larger such that              #
#    all parts of the newly rotated DMABN molecule remain 10 Ang. apart             #
#    from their images.                                                             #
#                                                                                   #
# 3. Use the translate functionality to center the entire molecule in the           #
#    new unit cell. This is purely for illustrative purposes and not for            #
#    any practical purposes.                                                        #
#                                                                                   #
#####################################################################################

################
# 1. Rotating  # - two calls of transform.py since the atoms are separated in two groups in the coordinate list
################

# rotating the carbon atoms in the two methyl groups

printf '1\n1\n[90,"D",[0.40822,0.50000,0.49887],[0.46319,0.50001,0.49945]]\n[9,10]\n' | python ../transform.py

# the two arrays in the above command are the direct coordinates of the nitrogen atom in the dimethylamino 
# group and the carbon atom in the benzene ring it connects to, respectively, defining the rotation vector.

# saving the original POSCAR which is now saved as POSCAR_old

cat POSCAR_old > POSCAR_original

# rotating the 3 hydrogens in each methyl group (6 hydrogens) by the same rotation as above

printf '1\n1\n[90,"D",[0.40822,0.50000,0.49887],[0.46319,0.50001,0.49945]]\n[16,21]\n' | python ../transform.py

# saving new POSCAR after this step which is now saved as POSCAR_old

cat POSCAR_old > POSCAR_step1

#################
# 2. Stretching #
#################

# stretching c lattice vector to 14 Ang.

printf '2\n1\n[[],[],[14]]\n' | python ../transform.py

# saving POSCAR from this step

cat POSCAR > POSCAR_step2

##################
# 3. Translating #
##################

# translating whole molecule by 2 Ang. along c vector

printf '0\n1\n["C",0,0,2]\n[]\n' | python ../transform.py

#  saving POSCAR from this step and renaming POSCAR_original back to POSCAR

cat POSCAR > POSCAR_step3
cat POSCAR_original > POSCAR

rm POSCAR_original
rm POSCAR_old

# do not need the z.py file

rm z.py

# will now have the following files:
   # POSCAR - which should be the original POSCAR prior to running this bash script
   # POSCAR_step1 - POSCAR after performing rotation of the dimethylamino group
   # POSCAR_step2 - POSCAR after performing stretch on POSCAR_step1
   # POSCAR_step3 - POSCAR after performing translate on POSCAR_step2
