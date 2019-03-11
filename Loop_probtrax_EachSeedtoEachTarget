#!/bin/bash

#supply a seeds.txt and a targets.txt file located in each subject's folder in /YOUR_BEDPOST_DIR/CoRegistration/ that gives #path to a set of seeds and targets in patient diffusion space.

homepath=/YOUR_HOME(ALLsubjects)_DIRECTORY/
cd $homepath

transform3=/YOUR_RIGID_TRANSFORM_IITFA_to_MNI.mat
transform4=/YOUR_AFFINE_TRANSFORM_IITFA_to_MNI.mat
transform5=/YOUR_NONLINEAR_WARP_IITFA_to_MNI.nii.gz

for d in $homepath/* #loop through all control scans
do
	cd $homepath/$d/YOUR_BEDPOST_DIR 
		#make post txform target directories
		mkdir probtrax
		mkdir -p probtrax/specific_targets
		
		#assign relevant txform files
		reference=/YOUR_REF_BRAIN.nii.gz
		transform1=/YOURSUBJ_to_IIT_affine_transform.mat
		transform2=/YOURSUBJ_TO_IIT_nonlinear_warp.nii.gz
		
		#do probtrackx2
		while read p
		do
		  #get the name of the seed we're using for the output directory
      seedname=`echo $p`
      #Get total number of paths emanating from seed. 5000*voxel#
		  max=`fslstats $p -V | awk -F ' ' '{print $1}'`
		  max_paths=`echo "$((5000 *$max))"`
		
      #save SEED max_paths into appropriate file
		  echo $max_paths >> /DATA_OUTPUT_FOLDER/max_paths_SEED.txt

        while read a  
			   do
			    #get target name
			    targetname=`echo $a`
          #get the total number of paths that emanate from target, when used as a seed.
			    maxTarget=`fslstats $a -V | awk -F ' ' '{print $1}'`
			    max_paths_Targ=`echo "$((5000 *$maxTarget))"`
			    #save TARGET max_paths into appropriate file
			    echo $max_paths_Targ >> /DATA_OUTPUT_FOLDER/max_paths_Targ.txt    
           
           ##Waypoint calculation, using callosal exclusion, stop mask, A:B calculation
            probtrackx2 -s bedpost_dir.bedpostX/merged -m bedpost_dir.bedpostX/nodif_brain_mask.nii.gz -x $p -- 
            dir="$seedname"_"$targetname"_termmask --waypoints=$a --stop=$a --avoid=CoRegistration/CORPUSCALLOSUM.nii.gz -l -c 
            0.2 -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --opd -V 1
			      
            #save waytotal number into variable, move resulting probtrax folder into appropriate directory
			      waytotalX=$(head -n 1 *$targetname*/waytotal)
            
            #make a normalized paths file for brining to MNI space
		        fslmaths *$targetname*/fdt_paths.nii.gz -div $max_paths fdt_paths_normalized.nii.gz
		        mv *normalized* *$targetname*
            #move folder into subject specific subdirectory
			      mv *$targetname* probtrax/specific_targets
		
		       ##Waypoint calculation, using callosal exclusion, stop mask, B:A calculation, store waytotal as variable, move 
            resulting probtrax folder into appropriate directory	
            probtrackx2 -s bedpost_dir.bedpostX/merged -m bedpost_dir.bedpostX/nodif_brain_mask.nii.gz -x $a -- 
            dir="$targetname"_"$seedname"_termmask --waypoints=$p --stop=$p --
            avoid=CoRegistration/MNIvolumes_WarpedTo_SubjDiff/WM_Callosal_Body_MNI_Juelich_thresh40_DiffSpace.nii.gz -l -c 0.2 
            -S 2000 --steplength=0.5 -P 5000 --fibthresh=0.01 --distthresh=0.0 --sampvox=0.0 --opd -V 1
            
            #save waytotal number into variable, move resulting probtrax folder into appropriate directory
			      waytotalB=$(head -n 1 *$targetname*/waytotal)
			      mv *$targetname* probtrax/specific_targets
	

		        #warp seed to target paths into MNI space
		        $ANTSPATH/antsApplyTransforms -d 3 -i probtrax/specific_targets/"$seedname"_"$targetname"_termmask/fdt_paths_normalized.nii.gz -r $reference -o       
            /YOUR_SUBJ_MNI_OUTPUTDIR/"$seedname"_"$targetname"_MNIpaths.nii.gz -t $transform5 -t       
            $transform4 -t $transform3 -t $transform2 -t $transform1 -v
		

		
            #store waytotal values into signle textfile
			echo $waytotalX >> /DATA_OUTPUT_FOLDER/waytotals.txt
			echo $waytotalB >> $homepath/sbsnider/response_processing/waytotalsBtoA.txt
			  done <CoRegistration/targets.txt
			
		  done < CoRegistration/seeds.txt


      cd $homepath


done















