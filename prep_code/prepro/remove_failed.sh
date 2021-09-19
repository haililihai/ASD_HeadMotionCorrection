#!/bin/bash
rundir=`pwd`
WD=`pwd`/../..

#for site in usm um nyu pitt2 stanford ohsu kki
for site in TR2_200 TR2_240 TR3_120 TR3_240
do
    for sublist in `cat ${WD}/${site}/sub.csv`
    do
        sub=$(echo ${sublist}|cut -d ',' -f 1)
        if [ -d "${WD}/${site}/rest/${sub}/01prepro" ]; then
            rm -rfv "${WD}/${site}/rest/${sub}/01prepro"
        fi
        if [ -d "${WD}/${site}/rest/${sub}/02prepro" ]; then
            rm -rfv "${WD}/${site}/rest/${sub}/02prepro"
        fi
        if [ -f "${WD}/${site}/rest/${sub}/01/rest.nii" ] || [ -d "${WD}/${site}/rest/${sub}/01/prepro" ]
        then
            echo "=== ${site}/${sub}/01 ==="
            rm -rvf ${WD}/${site}/rest/${sub}/01/prepro
            fslchfiletype NIFTI_GZ ${WD}/${site}/rest/${sub}/01/rest.nii
            cd ${WD}/${site}/anat/${sub}/01
            fslchfiletype NIFTI_GZ anat.nii
            rm -rvf `ls ./*|grep -v "^./anat.nii.gz$"`
            cd ${rundir}
        fi
        if [ -f "${WD}/${site}/rest/${sub}/02/rest.nii" ] || [ -d "${WD}/${site}/rest/${sub}/02/prepro" ]
        then
            echo "=== ${site}/${sub}/02 ==="
            rm -rvf ${WD}/${site}/rest/${sub}/02/prepro
            fslchfiletype NIFTI_GZ ${WD}/${site}/rest/${sub}/02/rest.nii
            cd ${WD}/${site}/anat/${sub}/02
            fslchfiletype NIFTI_GZ anat.nii
            rm -rvf `ls ./*|grep -v "^./anat.nii.gz$"`
            cd ${rundir}
        fi
    done
done