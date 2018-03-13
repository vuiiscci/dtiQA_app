#!/bin/bash -x

SCRIPTDIR=${0%/*}

# Set this to something else in order to use a specific Camino
BINPATH=$SCRIPTDIR/../../bin

PATH=${BINPATH}:${PATH}

TEST_COUNTER=1

if [ $# -gt 0 ]; then
    TEST_COUNTER=$1
fi

# Check where we are

if [[ $0 != "./tests.sh" ]]; then
    echo "Test script must be run from directory containing the script"
    exit 1
fi	

D2T="double2txt 6"
F2T="float2txt 6"
	

# Get rid of output from previous tests
rm -rf ${SCRIPTDIR}/crossing/results ${SCRIPTDIR}/torus/results ${SCRIPTDIR}/bedpost/results

mkdir ${SCRIPTDIR}/crossing/results ${SCRIPTDIR}/torus/results ${SCRIPTDIR}/bedpost/results

# NOTE: the file torus_06.Bdouble.gz has data dimensions (80, 32, 33). For some reason that now escapes me,
# I made a bunch of ROIs (80,32,12) - reading only part of the available data (which track lets you do). I'm in the
# process of replacing these.
#
# Some ROIs have dims (80,32,33). You have to be careful to specify whether you want all the tensor data to be read,
# by defining the input space correctly with -header etc. What you can't do is let the input space be defined
# as (80,32,12), then read in an (80,32,33) bgmask or something. 
# 
# I think I have got rid of all such instances at this point, but watch out for any I've missed.




# Track in single ROI and get connection probability maps 
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 or 1
track -tracker euler -stepsize 1 -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr -iterations 5000 -curvethresh 80 -curveinterval 2 | procstreamlines -iterations 5000 -seedfile torus/rois/torus_single_roi.hdr -outputcp -outputroot torus/results/torus_single_roi_cp_ -gzip

# Test removed 
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 1

# feed these to cpstats
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 2
for operation in max min median mean; do
imagestats -stat ${operation} -images torus/results/torus_single_roi_cp_1_* -outputroot torus/results/torus_single_roi_cp_output_${operation} -gzip
done

# Compare results
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 3
for operation in max min median mean; do
    imagessd torus/results/torus_single_roi_cp_output_${operation}.hdr torus/expectedResults/torus_single_roi_cp_output_${operation}.hdr
done

# Output tracts as from two ROIs in one pass
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 4
track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -randomseed 188 -seedfile torus/rois/torus_multi_roi.hdr | procstreamlines -seedfile torus/rois/torus_multi_roi.hdr  -outputroot torus/results/torus_multi_roi_streamlines_

diff --brief -N torus/results/torus_multi_roi_streamlines_1.Bfloat torus/expectedResults/torus_multi_roi_streamlines_1.Bfloat
diff --brief -N torus/results/torus_multi_roi_streamlines_2.Bfloat torus/expectedResults/torus_multi_roi_streamlines_2.Bfloat


# Split by ROI
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 5

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -randomseed 188 -seedfile torus/rois/torus_multi_roi.hdr -outputroot torus/results/torus_multi_roi_streamlines_

cat torus/results/torus_multi_roi_streamlines_1.Bfloat torus/results/torus_multi_roi_streamlines_2.Bfloat | procstreamlines -header torus/rois/torus_multi_roi.hdr > torus/results/torus_multi_roi_combined_streamlines.Bfloat

diff -N torus/results/torus_multi_roi_combined_streamlines.Bfloat torus/expectedResults/torus_multi_roi_combined_streamlines.Bfloat


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 6
cat torus/results/torus_multi_roi_streamlines_2.Bfloat | procstreamlines -seedfile torus/rois/torus_multi_roi.hdr -outputroot torus/results/torus_multi_roi_streamlines_regionindex_ -regionindex 2

# test similarity - use F2T
$F2T < torus/results/torus_multi_roi_streamlines_2.Bfloat > torus/results/torus_multi_roi_streamlines_2.txt
$F2T < torus/results/torus_multi_roi_streamlines_regionindex_2.Bfloat > torus/results/torus_multi_roi_streamlines_regionindex_2.txt

diff -N torus/expectedResults/torus_multi_roi_streamlines_regionindex_2.txt torus/expectedResults/torus_multi_roi_streamlines_regionindex_2.txt


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 7

# Single ROI streamline tracking with waypoint
track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_waypoints.hdr  | procstreamlines -waypointfile torus/rois/torus_waypoints.hdr  > torus/results/torus_single_roi_waypoint.Bfloat

diff -N --brief torus/results/torus_single_roi_waypoint.Bfloat torus/expectedResults/torus_single_roi_waypoint.Bfloat


# Single ROI streamline tracking with exclusion (discard)
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 8

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr | procstreamlines -exclusionfile torus/rois/torus_exclusion.hdr  > torus/results/torus_single_roi_discard_exclude.Bfloat

diff -N --brief torus/results/torus_single_roi_discard_exclude.Bfloat torus/expectedResults/torus_single_roi_discard_exclude.Bfloat


# Single ROI streamline tracking with exclusion (truncate)
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 9

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr | procstreamlines -exclusionfile torus/rois/torus_exclusion.hdr -truncateinexclusion  > torus/results/torus_single_roi_trunc_exclude.Bfloat

# CP maps with waypoints
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 10

track -tracker euler -stepsize 0.5 -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr -iterations 100 | procstreamlines -iterations 100 -seedfile torus/rois/torus_single_roi.hdr -outputsc -outputroot torus/results/torus_single_roi_waypoint_sc_ -waypointfile torus/rois/torus_waypoints.hdr -gzip

imagestats -stat max -images torus/results/torus_single_roi_waypoint_sc_1_* -outputroot torus/results/torus_single_roi_waypoint_sc_max_1 -gzip

imagessd torus/results/torus_single_roi_waypoint_sc_max_1.img.gz torus/expectedResults/torus_single_roi_waypoint_sc_max_1.img.gz


# Get target CP maps from two ROIs
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 11

track -tracker euler -stepsize 1 -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -seedfile torus/rois/torus_multi_roi.hdr -iterations 500 -randomseed 188 | procstreamlines -iterations 500 -seedfile torus/rois/torus_multi_roi.hdr -outputcp -outputroot torus/results/torus_multi_roi_target_cp_ -targetfile torus/rois/torus_targets.hdr -gzip  

for i in 1 2 ; do for j in 1 2; do 
imagessd torus/results/torus_multi_roi_target_cp_${i}_${j}_1.img.gz torus/expectedResults/torus_multi_roi_target_cp_${i}_${j}_1.img.gz
done
done

# Test voxel to physical transformation with NIFTI

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 12

# Track with 1mm seed points, then output results in 2mm space
nii2dt -inputfile torus/torus_06_qformdt.nii.gz | track -inputmodel dt -header torus/torus_06_qformdt.nii.gz -seedfile torus/rois/torus_06_qform_1mm_seeds.nii.gz | procstreamlines -seedfile torus/rois/torus_06_qform_seed.nii.gz -outputacm -outputroot torus/results/torus_qform_1mm_seed_2mm_

nii2dt -inputfile torus/torus_06_qformdt.nii.gz | track -tracker euler -interpolator linear -stepsize 0.5 -inputmodel dt -header torus/torus_06_qformdt.nii.gz -seedfile torus/rois/torus_06_qform_1mm_seeds.nii.gz | procstreamlines -waypointfile torus/rois/torus_06_qform_1mm_waypoints.nii.gz -outputacm -outputroot torus/results/torus_qform_1mm_

imagessd torus/results/torus_qform_1mm_acm_sc.nii.gz  torus/expectedResults/torus_qform_1mm_acm_sc.nii.gz   
imagessd torus/results/torus_qform_1mm_seed_2mm_acm_sc.nii.gz  torus/expectedResults/torus_qform_1mm_seed_2mm_acm_sc.nii.gz  


# Get CBS from the single ROI
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 13

track -tracker euler -stepsize 1 -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -seedfile torus/rois/torus_single_roi.hdr -iterations 200 -randomseed 188 | procstreamlines -iterations 200 -seedfile torus/rois/torus_single_roi.hdr -outputcbs -outputsc -outputroot torus/results/torus_single_roi_cbs_ -targetfile torus/rois/torus_targets.hdr -gzip

imagessd torus/results/torus_single_roi_cbs_labels_1_1.img.gz torus/expectedResults/torus_single_roi_cbs_labels_1_1.img.gz 
imagessd torus/results/torus_single_roi_cbs_labelsc_1_1.img.gz torus/expectedResults/torus_single_roi_cbs_labelsc_1_1.img.gz

# CP with two PDs in seed
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 14
track -inputmodel multitensor -inputfile crossing/crossing_05_0505_90.multitensor.Bdouble.gz -seedfile crossing/rois/crossing_multi_roi.hdr | procstreamlines -seedfile crossing/rois/crossing_multi_roi.hdr -outputcp -outputroot crossing/results/crossing_multi_roi_cp_ -iterations 1 -gzip

imagessd crossing/results/crossing_multi_roi_cp_1_1_1.img.gz crossing/expectedResults/crossing_multi_roi_cp_1_1_1.img.gz 
imagessd crossing/results/crossing_multi_roi_cp_1_2_1.img.gz crossing/expectedResults/crossing_multi_roi_cp_1_2_1.img.gz 

imagessd crossing/results/crossing_multi_roi_cp_2_1_1.img.gz crossing/expectedResults/crossing_multi_roi_cp_2_1_1.img.gz 
imagessd crossing/results/crossing_multi_roi_cp_2_2_1.img.gz crossing/expectedResults/crossing_multi_roi_cp_2_2_1.img.gz 

imagessd crossing/results/crossing_multi_roi_cp_2_1_2.img.gz crossing/expectedResults/crossing_multi_roi_cp_2_1_2.img.gz 
imagessd crossing/results/crossing_multi_roi_cp_2_2_2.img.gz crossing/expectedResults/crossing_multi_roi_cp_2_2_2.img.gz 



# Target CP with two PDs in seed
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 15
track -inputmodel multitensor -inputfile crossing/crossing_05_0505_90.multitensor.Bdouble.gz -seedfile crossing/rois/crossing_multi_roi.hdr | procstreamlines -seedfile crossing/rois/crossing_multi_roi.hdr -outputsc -outputroot crossing/results/crossing_multi_roi_targetsc_ -iterations 1 -targetfile crossing/rois/crossing_targets.hdr -gzip

imagessd crossing/results/crossing_multi_roi_targetsc_1_1_1.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_1_1_1.img.gz 
imagessd crossing/results/crossing_multi_roi_targetsc_1_2_1.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_1_2_1.img.gz 

imagessd crossing/results/crossing_multi_roi_targetsc_2_1_1.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_2_1_1.img.gz 
imagessd crossing/results/crossing_multi_roi_targetsc_2_2_1.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_2_2_1.img.gz 

imagessd crossing/results/crossing_multi_roi_targetsc_2_1_2.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_2_1_2.img.gz 
imagessd crossing/results/crossing_multi_roi_targetsc_2_2_2.img.gz  crossing/expectedResults/crossing_multi_roi_targetsc_2_2_2.img.gz 


# CBS with two PDs in seed
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 16
track -inputmodel multitensor -inputfile crossing/crossing_05_0505_90.multitensor.Bdouble.gz -seedfile crossing/rois/crossing_multi_roi.hdr | procstreamlines -seedfile crossing/rois/crossing_multi_roi.hdr -outputcbs -outputroot crossing/results/crossing_multi_roi_ -iterations 1 -targetfile crossing/rois/crossing_targets.hdr -gzip 

imagessd crossing/results/crossing_multi_roi_labels_1_1.img.gz crossing/expectedResults/crossing_multi_roi_labels_1_1.img.gz     
imagessd crossing/results/crossing_multi_roi_labelsc_1_1.img.gz crossing/expectedResults/crossing_multi_roi_labelsc_1_1.img.gz  

imagessd crossing/results/crossing_multi_roi_labels_1_1.img.gz crossing/expectedResults/crossing_multi_roi_labels_1_1.img.gz     
imagessd crossing/results/crossing_multi_roi_labelsc_1_1.img.gz crossing/expectedResults/crossing_multi_roi_labelsc_1_1.img.gz   

imagessd crossing/results/crossing_multi_roi_labels_2_1.img.gz crossing/expectedResults/crossing_multi_roi_labels_2_1.img.gz     
imagessd crossing/results/crossing_multi_roi_labelsc_2_1.img.gz crossing/expectedResults/crossing_multi_roi_labelsc_2_1.img.gz  

imagessd crossing/results/crossing_multi_roi_labels_2_2.img.gz crossing/expectedResults/crossing_multi_roi_labels_2_2.img.gz     
imagessd crossing/results/crossing_multi_roi_labelsc_2_2.img.gz crossing/expectedResults/crossing_multi_roi_labelsc_2_2.img.gz       

# OOGL binary format - removed
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 17


# targetprobs2txt
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 18

targetprobs2txt -inputroot torus/expectedResults/torus_multi_roi_target_cp_ -seedfile torus/rois/torus_multi_roi.hdr -targetfile torus/rois/torus_targets.hdr -regionindex 1 > torus/results/torus_multi_roi_target_cp_r1.txt

targetprobs2txt -inputroot torus/expectedResults/torus_multi_roi_target_cp_ -seedfile torus/rois/torus_multi_roi.hdr -targetfile torus/rois/torus_targets.hdr -regionindex 2 > torus/results/torus_multi_roi_target_cp_r2.txt

diff -N torus/results/torus_multi_roi_target_cp_r1.txt torus/expectedResults/torus_multi_roi_target_cp_r1.txt
diff -N torus/results/torus_multi_roi_target_cp_r2.txt torus/expectedResults/torus_multi_roi_target_cp_r2.txt

# for crossing
targetprobs2txt -inputroot crossing/expectedResults/crossing_multi_roi_targetsc_ -seedfile crossing/rois/crossing_multi_roi.hdr -targetfile crossing/rois/crossing_targets.hdr -pd 2 -regionindex 2 > crossing/results/crossing_multi_roi_target2txt.txt

diff -N crossing/results/crossing_multi_roi_target2txt.txt crossing/expectedResults/crossing_multi_roi_target2txt.txt



# Test -seedindex - deprecated
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 19


# Test bootstrap tracking

# Generate bootstrap data

cat torus/torus_06.Bdouble.gz | gunzip -c | datasynth -inputmodel dt -snr 16 -schemefile schemefiles/M1_N30_b1000.scheme -inputdatatype double -seed 273 -outputfile torus/torus_bootstrap_1.Bfloat.gz -gzip

cat torus/torus_06.Bdouble.gz | gunzip -c | datasynth -inputmodel dt -snr 16 -schemefile schemefiles/M1_N30_b1000.scheme -inputdatatype double -seed 274 -outputfile torus/torus_bootstrap_2.Bfloat.gz -gzip

cat crossing/crossing_05_0505_90.multitensor.Bdouble.gz | gunzip -c | datasynth -inputmodel multitensor -snr 16 -schemefile schemefiles/M1_N30_b1000.scheme -gzip -inputdatatype double -seed 273 -outputfile crossing/crossing_bootstrap_1.Bfloat.gz

cat crossing/crossing_05_0505_90.multitensor.Bdouble.gz | gunzip -c | datasynth -inputmodel multitensor -snr 16 -schemefile schemefiles/M1_N30_b1000.scheme -gzip -inputdatatype double -seed 274 -outputfile crossing/crossing_bootstrap_2.Bfloat.gz


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 20

track -tracker euler -interpolator prob_nn -stepsize 1 -inputmodel repbs_dt -randomseed 39395 -inputdatatype float -bsdatafiles torus/torus_bootstrap_1.Bfloat.gz torus/torus_bootstrap_2.Bfloat.gz -schemefile schemefiles/M1_N30_b1000.scheme -inversion 1 -iterations 1000 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 1000 -outputsc -outputroot torus/results/torus_bootstrap_sc_1_ -gzip

imagessd torus/results/torus_bootstrap_sc_1_1_1_1.img.gz torus/expectedResults/torus_bootstrap_sc_1_1_1_1.img.gz


# Bootstrap with crossing
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 21

track -inputmodel repbs_multitensor -randomseed 39395 -inputdatatype float -tracker euler -stepsize 1 -bsdatafiles crossing/crossing_bootstrap_1.Bfloat.gz crossing/crossing_bootstrap_2.Bfloat.gz -schemefile schemefiles/M1_N30_b1000.scheme -inversion 31 -voxclassmap crossing/crossing_vc.Bint -iterations 1000 -seedfile crossing/rois/crossing_seedpoint_0_16_7.hdr | procstreamlines -seedfile crossing/rois/crossing_seedpoint_0_16_7.hdr -iterations 1000 -outputsc -outputroot crossing/results/crossing_bootstrap_sc_1_ -gzip

imagessd crossing/results/crossing_bootstrap_sc_1_1_1_1.img.gz crossing/expectedResults/crossing_bootstrap_sc_1_1_1_1.img.gz



echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 22

# Single ROI streamline tracking with waypoint, turn on loop checking
track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_waypoints.hdr  | procstreamlines -truncateloops -waypointfile torus/rois/torus_waypoints.hdr > torus/results/torus_single_roi_waypoint_looptrunc.Bfloat

# There are no loops, so should get the same output as without -truncateloops
diff --brief -N torus/results/torus_single_roi_waypoint_looptrunc.Bfloat torus/expectedResults/torus_single_roi_waypoint.Bfloat



# VTK streamlines
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 23

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr | vtkstreamlines -seedfile torus/rois/torus_single_roi.hdr > torus/results/torus_single_roi.vtk

diff -N torus/results/torus_single_roi.vtk torus/expectedResults/torus_single_roi.vtk



# Wild bootstrap

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 24

track -inputmodel wildbs_dt -randomseed 39395 -inputdatatype float -inputfile torus/torus_bootstrap_1.Bfloat.gz -schemefile schemefiles/M1_N30_b1000.scheme -iterations 1000 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 1000 -outputsc -outputroot torus/results/torus_wildbs_sc_1_ -gzip 

imagessd torus/results/torus_wildbs_sc_1_1_1_1.img.gz torus/expectedResults/torus_wildbs_sc_1_1_1_1.img.gz

# Use bgmask, anisfile
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 25

track -inputmodel wildbs_dt -randomseed 39395 -inputdatatype float -inputfile torus/torus_bootstrap_1.Bfloat.gz -schemefile schemefiles/M1_N30_b1000.scheme -iterations 1000 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -header torus/rois/torus_halfmask.hdr -bgmask torus/rois/torus_halfmask.hdr | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 1000 -outputsc -outputroot torus/results/torus_wildbs_mask_sc_1_ -gzip

track -inputmodel wildbs_dt -randomseed 39395 -inputdatatype float -inputfile torus/torus_bootstrap_1.Bfloat.gz -schemefile schemefiles/M1_N30_b1000.scheme -iterations 1000 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -header torus/rois/torus_halfmask.hdr -anisthresh 1 -anisfile torus/rois/torus_halfmask.hdr | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 1000 -outputsc -outputroot torus/results/torus_wildbs_anismask_sc_1_ -gzip

imagessd torus/results/torus_wildbs_mask_sc_1_1_1_1.img.gz torus/expectedResults/torus_wildbs_mask_sc_1_1_1_1.img.gz

imagessd torus/results/torus_wildbs_anismask_sc_1_1_1_1.img.gz torus/expectedResults/torus_wildbs_anismask_sc_1_1_1_1.img.gz


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 26

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -seed 22382 -seedfile torus/rois/torus_halfmask.hdr | procstreamlines -outputroot torus/results/torus_ -iterations 1 -seedfile torus/rois/torus_halfmask.hdr -outputacm -outputcp -gzip

imagessd torus/results/torus_acm_cp.img.gz torus/expectedResults/torus_acm_cp.img.gz



echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 27

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -seedfile torus/rois/torus_multi_roi.hdr | vtkstreamlines -seedfile torus/rois/torus_multi_roi.hdr > torus/results/torus_multi_roi.vtk

diff -N torus/results/torus_multi_roi.vtk torus/expectedResults/torus_multi_roi.vtk


# Tract statistics

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 28


dtfit torus/torus_bootstrap_1.Bfloat.gz schemefiles/M1_N30_b1000.scheme > torus/torus_bs_1.inv1.Bdouble

fa < torus/torus_bs_1.inv1.Bdouble > torus/torus_bs_1.fa.img

analyzeheader -datatype double -datadims 80 32 33 -voxeldims 2 2 2 > torus/torus_bs_1.fa.hdr

rm -f torus/results/tractStats.Bdouble

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble -header torus/torus_bs_1.fa.hdr -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  > torus/results/single_streamline.Bfloat

for i in mean var min max median meanvar sum; do 
    
     cat torus/results/single_streamline.Bfloat | tractstats -tractstat ${i} -scalarfile torus/torus_bs_1.fa.hdr  >> torus/results/tractStats.Bdouble

done

cat torus/results/single_streamline.Bfloat | tractstats -tractstat length -datadims 100 100 100 -voxeldims 2 2 2 >> torus/results/tractStats.Bdouble


# diff doesn't work too well on small binary files

cat torus/results/tractStats.Bdouble | $D2T > torus/results/tractStats.txt

diff -N -b torus/results/tractStats.txt torus/expectedResults/tractStats.txt


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 29

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble  -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  | tractstatimage -scalarfile torus/torus_bs_1.fa.hdr -countintersect -tractstat mean -imagestat mean -outputroot torus/results/torus_faMean_all -gzip

imagessd torus/results/torus_faMean_all.img.gz torus/expectedResults/torus_faMean_all.img.gz


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 30

# Test deterministic interpolation

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble -seedfile torus/rois/torus_seedpoint_17_0_7.hdr -tracker euler -interpolator linear -stepsize 0.2 > torus/results/single_streamline_vec_interp.Bfloat

cat torus/results/single_streamline_vec_interp.Bfloat | $F2T > torus/results/single_streamline_vec_interp.txt
cat torus/expectedResults/single_streamline_vec_interp.Bfloat | $F2T > torus/results/single_streamline_vec_interp_expected.txt

diff -N --brief torus/results/single_streamline_vec_interp.txt torus/expectedResults/single_streamline_vec_interp_expected.txt



echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 31


# Test NC
track  -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr -tracker euler -interpolator prob_nn -stepsize 0.2 | procstreamlines  -seedfile torus/rois/torus_seedpoint_17_0_7.hdr -outputcp -outputroot torus/results/torus_single_seed_cp_ncinterp_ -gzip

imagessd torus/results/torus_single_seed_cp_ncinterp_1_1_1.img.gz torus/expectedResults/torus_single_seed_cp_ncinterp_1_1_1.img.gz


# Bayesian tracking on bootstrap torus data

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 32

cat torus/torus_06.Bdouble.gz | gunzip -c | datasynth -inputmodel dt -snr 16 -schemefile schemefiles/bmx7.scheme -inputdatatype double -seed 273 -outputfile torus/torus_bmx7_snr16.Bfloat.gz -gzip

track -randomseed 188 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -inputmodel bayesdirac_dt -inputfile torus/torus_bmx7_snr16.Bfloat.gz -schemefile schemefiles/bmx7.scheme -iterations 50 -curvepriorg 1.0 | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 50 -outputsc -outputroot torus/results/torus_single_seed_bayescyldt_ -gzip 

track -randomseed 188 -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -inputmodel bayesdirac -inputfile torus/torus_bmx7_snr16.Bfloat.gz -schemefile schemefiles/bmx7.scheme -iterations 50 -curvepriork 1.0 | procstreamlines -seedfile torus/rois/torus_seedpoint_17_0_7.hdr  -iterations 50 -outputsc -outputroot torus/results/torus_single_seed_bayesballstick_ -gzip 

imagessd torus/results/torus_single_seed_bayescyldt_1_1_1.img.gz torus/expectedResults/torus_single_seed_bayescyldt_1_1_1.img.gz

imagessd torus/results/torus_single_seed_bayesballstick_1_1_1.img.gz torus/expectedResults/torus_single_seed_bayesballstick_1_1_1.img.gz


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 33


# Check -curvethresh, -ipthresh equivalent

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -seedfile torus/rois/torus_seedpoint_17_0_7.hdr -ipthresh 0.996 > torus/results/ipThresh.Bfloat

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -seedfile torus/rois/torus_seedpoint_17_0_7.hdr -curvethresh 5.1264 > torus/results/curveThresh.Bfloat

diff -N torus/results/ipThresh.Bfloat torus/results/curveThresh.Bfloat


# test end points 

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 34

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz -seedfile torus/rois/torus_endpoints_seed.hdr | procstreamlines -endpointfile torus/rois/torus_endpoints.hdr | vtkstreamlines -header torus/rois/torus_endpoints.hdr > torus/results/endpoints.vtk

diff -N torus/results/endpoints.vtk torus/expectedResults/endpoints.vtk


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 35

# test nifti seeds / waypoints - no longer needed, this gets tested later


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 36

# test MHA seeds / waypoints - still allowed but deprecated, transform to physical space too poorly defined


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 37

# TEND

track -inputmodel multitensor -inputfile crossing/crossing_05_0505_90.snr20.Bdouble.gz -tracker euler -stepsize 1 -interpolator tend -seedfile crossing/rois/crossing_multi_roi.hdr -tendg 0.6 -gzip > crossing/results/crossing_multi_roi_tend.Bfloat.gz

diff -N crossing/results/crossing_multi_roi_tend.Bfloat.gz crossing/expectedResults/crossing_multi_roi_tend.Bfloat.gz


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 38

# tract stats, output raw values
cat torus/results/single_streamline.Bfloat | tractstats -tractstat none -scalarfile torus/torus_bs_1.fa.hdr | $D2T > torus/results/tractFA.txt

diff -N torus/results/tractFA.txt torus/expectedResults/tractFA.txt


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 39

# Tend with probabilistic nearest neighbour

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble -randomseed 188 -seedfile torus/rois/torus_single_roi.hdr -interpolator tend_prob_nn -tracker euler -stepsize 1 -iterations 100 | procstreamlines -header torus/rois/torus_single_roi.hdr -outputacm -outputroot torus/results/torus_tend_nc_

imagessd torus/results/torus_tend_nc_acm_sc.hdr torus/expectedResults/torus_tend_nc_acm_sc.hdr

# Test fix of bug with non-contiguous labels 

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 40

for i in 1 2 3; do 
track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -randomseed 188 -seedfile torus/rois/torus_multi_spacedlabels_roi.hdr -regionindex $i > torus/results/torus_spaced_labels_r${i}.Bfloat

diff -N torus/results/torus_spaced_labels_r${i}.Bfloat torus/expectedResults/torus_spaced_labels_r${i}.Bfloat

done


# ACM with targets
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 41

track -inputmodel dt -inputfile torus/torus_06.Bdouble.gz  -randomseed 188 -seedfile torus/rois/torus_waypoints.hdr | procstreamlines -outputacm -targetfile torus/rois/torus_targets.hdr -outputroot torus/results/torus_acm_targets_ -gzip

imagessd torus/results/torus_acm_targets_acm_target_sc.img.gz torus/expectedResults/torus_acm_targets_acm_target_sc.img.gz

# vtkstreamlines, no scalars
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 42

vtkstreamlines -inputfile torus/expectedResults/torus_spaced_labels_r64.Bfloat -outputfile torus/results/torus_spaced_labels_r64.vtk

diff -N torus/expectedResults/torus_spaced_labels_r64.vtk torus/results/torus_spaced_labels_r64.vtk


# conmat
echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 43

track -inputfile torus/torus_06.Bdouble.gz -inputmodel dt -seedfile torus/rois/torus_conmat_targets.nii.gz -anisthresh 0.25 | conmat -targetfile torus/rois/torus_conmat_targets.nii.gz -outputroot torus/results/torus_conmat_ -targetnamefile torus/rois/conmat_names.txt -scalarfile torus/fa.nii.gz -tractstat max

diff -N torus/expectedResults/torus_conmat_sc.csv torus/results/torus_conmat_sc.csv

diff -N torus/expectedResults/torus_conmat_ts.csv torus/results/torus_conmat_ts.csv


# cbsmat

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 44

track -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 188 -seedfile torus/rois/torus_cbs_seeds.nii.gz -iterations 100 | cbsmat -targetfile torus/rois/torus_cbs_targets.nii.gz -seedfile torus/rois/torus_cbs_seeds.nii.gz -outputroot torus/results/torus_cbstest_

diff -N torus/expectedResults/torus_cbstest_sc.csv torus/results/torus_cbstest_sc.csv

for i in sc labels; do
  imagessd torus/results/torus_cbstest_cbs_${i}.nii.gz torus/expectedResults/torus_cbstest_cbs_${i}.nii.gz
done


# Test some new tracking features 

echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 45

# User can set how regularly we check curvature

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble -randomseed 1088 -seedfile torus/rois/torus_06_qform_seed.nii.gz -tracker fact -curvethresh 30 -curveinterval 30 | vtkstreamlines -outputfile torus/results/torus_curve_30.vtk

diff -N torus/results/torus_curve_30.vtk torus/expectedResults/torus_curve_30.vtk


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 46

# RK4 

track -inputmodel dt -inputfile torus/torus_bs_1.inv1.Bdouble -randomseed 1088 -seedfile torus/rois/torus_06_qform_seed.nii.gz -tracker rk4 -stepsize 1 -interpolator linear | vtkstreamlines -outputfile torus/results/torus_rk4.vtk

diff -N torus/results/torus_rk4.vtk torus/expectedResults/torus_rk4.vtk


track -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 1880 -iterations 100 -seedfile torus/rois/torus_single_seed.nii.gz -tracker euler -stepsize 0.5 -interpolator nn | vtkstreamlines -outputfile torus/results/torus_euler_prob.vtk

diff -N torus/results/torus_euler_prob.vtk torus/expectedResults/torus_euler_prob.vtk

track -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 1880 -iterations 100 -seedfile torus/rois/torus_single_seed.nii.gz -tracker euler -stepsize 0.5 -interpolator prob_nn | vtkstreamlines -outputfile torus/results/torus_euler_prob_nn.vtk

diff -N torus/results/torus_euler_prob_nn.vtk torus/expectedResults/torus_euler_prob_nn.vtk

track -inputmodel pico -inputfile torus/torus_06.bingham.Bdouble.gz -randomseed 1880 -iterations 100 -seedfile torus/rois/torus_single_seed.nii.gz -tracker rk4 -stepsize 0.5 -interpolator linear | vtkstreamlines -outputfile torus/results/torus_rk4_prob.vtk


diff -N torus/results/torus_rk4_prob.vtk torus/expectedResults/torus_rk4_prob.vtk


echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 47

# DWI linear interpolation

track -randomseed 1880 -seedfile torus/rois/torus_cbs_seeds.nii.gz -inputmodel dwi_dt -inputfile torus/torus_bmx7_snr16.Bfloat.gz -schemefile schemefiles/bmx7.scheme -tracker euler -interpolator dwi_linear -stepsize 0.5 -model ldt_wtd | vtkstreamlines -outputfile torus/results/torus_dwi_interp.vtk

diff -N torus/results/torus_dwi_interp.vtk torus/expectedResults/torus_dwi_interp.vtk

cat crossing/crossing_05_0707_90.multitensor.Bdouble.gz | gunzip -c | datasynth -inputmodel multitensor -snr 20 -schemefile schemefiles/M6_N54_b1500.scheme -outputfile crossing/crossing_05_0707_90_snr20.Bfloat.gz -gzip

track -randomseed 1880 -seedfile crossing/rois/crossing_multi_rois.nii.gz -inputmodel dwi_multitensor -inputfile  crossing/crossing_05_0707_90_snr20.Bfloat.gz -schemefile schemefiles/M6_N54_b1500.scheme -tracker euler -interpolator dwi_linear -stepsize 1 -model pospos ldt_wtd -voxclassmap crossing/crossing_vc.Bint | vtkstreamlines -outputfile crossing/results/crossing_dwi_interp.vtk

diff -N crossing/results/crossing_dwi_interp.vtk crossing/expectedResults/crossing_dwi_interp.vtk

track -randomseed 1880 -seedfile crossing/rois/crossing_multi_rois.nii.gz -inputmodel dwi_multitensor -inputfile  crossing/crossing_05_0707_90_snr20.Bfloat.gz -schemefile schemefiles/M6_N54_b1500.scheme -tracker euler -interpolator tend -stepsize 1 -model pospos ldt_wtd -voxclassmap crossing/crossing_vc.Bint | vtkstreamlines -outputfile crossing/results/crossing_tend_interp.vtk

diff -N crossing/results/crossing_tend_interp.vtk crossing/expectedResults/crossing_tend_interp.vtk



echo "TRACKING TEST $TEST_COUNTER"; TEST_COUNTER=$(($TEST_COUNTER + 1)) # $1 + 48

# BedpostX

track -randomseed 1880 -seedfile bedpost/rois/seed.nii.gz -inputmodel bedpostx_dyad -bedpostxdir bedpost/data -seedfile bedpost/rois/seed.nii.gz -interpolator prob_nn -tracker euler -stepsize 0.25 -outputfile bedpost/results/bedpost_dyad.Bfloat

procstreamlines -inputfile bedpost/results/bedpost_dyad.Bfloat -seedfile bedpost/rois/seed.nii.gz -outputacm -outputroot bedpost/results/dyad_

imagessd bedpost/results/dyad_acm_sc.nii.gz bedpost/expectedResults/dyad_acm_sc.nii.gz 


track -randomseed 1880 -seedfile bedpost/rois/seed.nii.gz -inputmodel bedpostx -bedpostxdir bedpost/data -seedfile bedpost/rois/seed.nii.gz -interpolator prob_nn -tracker euler -stepsize 0.25 -iterations 100 -outputfile bedpost/results/bedpostx.Bfloat

procstreamlines -inputfile bedpost/results/bedpostx.Bfloat -seedfile bedpost/rois/seed.nii.gz -outputacm -outputroot bedpost/results/bedpostx_

imagessd bedpost/results/bedpostx_acm_sc.nii.gz bedpost/expectedResults/bedpostx_acm_sc.nii.gz 
