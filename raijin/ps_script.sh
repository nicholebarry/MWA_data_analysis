#!/bin/bash

######################################################################################
# Top level script to integrate healpix cubes and run power spectrum code.
#
# A file path to the fhd directory is needed.
# 
# A file path to a text file listing observation ids OR preintegrated subcubes is
# needed.
# 
# If a text file of observation ids to be used in integration is specified, the obs 
# ids are assumed to be seperated by newlines.
#
# If a text file of preintegrated subcubes is specified, the format should be
# the name of the save file seperated by newlines.  "even_cube.sav" and "odd_cube.sav"
# is not necessary to include, as both will be used anyways.  The subcubes are
# assumed to be in <fhd_directory>/Healpix/. If elsewhere in the FHD directory, the 
# name of the subcubes must specify this in the text file as Other_than_Healpix/<name>.
#
# Set -ps to 1 to skip integration and make cubes only.
# 
# NOTE: print statements must be turned off in idl_startup file (e.g. healpix check)
######################################################################################

#Parse flags for inputs
while getopts ":d:f:p:w:n:m:o:l:" option
do
   case $option in
        d) FHDdir="$OPTARG";;			#file path to fhd directory with cubes
        f) integrate_list="$OPTARG";;		#txt file of obs ids or subcubes or a single obsid
        p) priority=$OPTARG;;           	#priority level for grid engine qsub
        w) wallclock_time=$OPTARG;;     	#Time for execution in grid engine
        n) ncpus=$OPTARG;;             	#Number of slots for grid engine
        m) mem=$OPTARG;;                	#Memory per core for grid engine
	o) ps_only=$OPTARG;;			#Flag for skipping integration to make PS only
        l) legacy=$OPTARG;;                     #Use legacy directory structure. hacky solution to a silly problem.
        \?) echo "Unknown option: Accepted flags are -d (file path to fhd directory with cubes), -f (obs list or subcube path or single obsid), "
	    echo "-p (priority for grid engine), -w (wallclock time), -n (number of slots), -m (memory allocation), "
	    echo "-o (make ps only without integration), and -l (legacy flag for old file structure),"
	    echo "-h (hold int/ps script on a running job id), and -t (apply a tukey window filter during ps),"
            exit 1;;
        :) echo "Missing option argument for input flag"
           exit 1;;
   esac
done

#Manual shift to the next flag
shift $(($OPTIND - 1))

#Throw error if no file path to FHD directory
if [ -z ${FHDdir} ]
then
   echo "Need to specify a file path to a FHD directory with cubes: Example /nfs/complicated_path/fhd_mine/"
   exit 1
fi

#Throw error if file path does not exist
if [ ! -d "$FHDdir" ]
then
   echo "Argument after flag -d is not a real directory. Argument should be the file path to the location of cubes to integrate."
   exit 1
fi

#Remove extraneous / on FHD directory if present
if [[ $FHDdir == */ ]]; then FHDdir=${FHDdir%?}; fi

#Error if integrate_list is not set
if [ -z ${integrate_list} ]
then
    echo "Need to specify obs list file path or preintegrated subcubes list file path with option -f"
    exit 1
fi

#Warning if integrate list filename does not exist
if [ ! -e "$integrate_list" ]
then
    echo "Integrate list is either not a file or the file does not exist!"
    echo "Assuming the integrate list is a single observation id."

    if [ -z ${ps_only} ]
    then
        echo "ps_only flag must be set if integrate list is a single observation id. Set -o 1 if desired function"
        exit 1
    fi 
    version=$integrate_list  #Currently assuming that the integrate list is a single obsid
else
    version=$(basename $integrate_list) # get filename
    version="${version%.*}" # strip extension
fi

#Default priority if not set.
if [ -z ${priority} ]; then priority=0; fi

#Set typical wallclock_time for standard PS with obs ids if not set.
if [ -z ${wallclock_time} ]; then wallclock_time=10:00:00; fi

#Set typical slots needed for standard PS with obs ids if not set.
if [ -z ${ncpus} ]; then ncpus=1; fi

#Set typical memory needed for standard PS with obs ids if not set.
if [ -z ${mem} ]; then mem=25G; fi

#Set default to do integration
if [ -z ${ps_only} ]; then ps_only=0; fi

#Set default to use current file structure
if [ -z ${legacy} ]; then legacy=0; fi

### NOTE this only works if idlstartup doesn't have any print statements (e.g. healpix check)
PSpath=$(idl -e 'print,rootdir("eppsilon")')

#Versions made during integrate list logic check above
echo Version is $version

if [ ! -e "$integrate_list" ]
then
    first_line=$integrate_list
else
    first_line=$(head -n 1 $integrate_list)
fi

first_line_len=$(echo ${#first_line})

rm -f ${FHDdir}/Healpix/${version}_int_chunk*.txt # remove any old chunk files lying around

exit_flag=0
n_pol=2 #default
#Check that cubes or integrated cubes are present, print and error if they are not
if [ "$ps_only" -ne "1" ]; then 	#only if we're integrating
while read line
do
   if [ "$first_line_len" == 10 ]; then
      if [ "$legacy" -ne "1" ]; then
	  if ! ls $FHDdir/Healpix/$line*cube*.sav &> /dev/null; then
              echo Missing cube for obs $line
	      if [ -z "$hold" ]; then
		  exit_flag=1
	      fi
	  fi
      pol_files=$(($(ls $FHDdir/Healpix/$line*cube*.sav | wc -l) / 2))
      if [[ "$pol_files" -eq 2 ]]; then 
          n_pol=2; else
              if [[ "$pol_files" -eq 4 ]]; then
                  n_pol=4; else 
                      echo Non-standard amount of cubes per obs $line
                      exit_flag=1
              fi
      fi 
      else
	  if ! ls $FHDdir/$line*cube*.sav &> /dev/null; then
	      echo Missing cube for obs $line
	      exit_flag=1
	  fi
      fi
   else
      if [[ "$first_line" != */* ]]; then
	 check=$FHDdir/Healpix/$line*.sav
      else
	 check=$FHDdir/$line*.sav
      fi
      if ! ls $check &> /dev/null; then
	    echo Missing save file for $line
	    exit_flag=1
      fi
   fi
done < $integrate_list
fi

if [ "$exit_flag" -eq 1 ]; then exit 1; fi

if [ "$first_line_len" == 10 ]; then
    
    # read in obs ids 100 at a time and divide into chunks to integrate in parallel mode
    obs=0   

    while read line
    do
        ((chunk=obs/100+1))		#integer division results in chunks labeled 0 (first 100), 1 (second 100), etc
        echo $line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt	#put that obs id into the right txt file
        ((obs++))			#increment obs for the next run through
    done < $integrate_list
    nchunk=$chunk 			#number of chunks we ended up with

else

    if [[ "$first_line" != */* ]]; then
   
        chunk=0 
        while read line
        do
            echo $line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt        #put that obs id into the right txt file
        done < $integrate_list
        nchunk=$chunk                       #number of chunks we ended up with
    
    else

        chunk=0 
        while read line
        do
            echo $line >> ${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt        #put that obs id into the right txt file
        done < $integrate_list
        nchunk=$chunk                       #number of chunks we ended up with

    fi

fi

if [[ "$n_pol" -eq 2 ]]; then pol_list=( XX YY );fi
if [[ "$n_pol" -eq 4 ]]; then pol_list=( XX YY XY YX );fi

unset idlist
if [ "$ps_only" -ne "1" ]; then   
    if [ "$nchunk" -gt "1" ]; then

        # set up files for master integration
        sub_cubes_list=${FHDdir}/Healpix/${version}_sub_cubes.txt
        rm $sub_cubes_list # remove any old lists

        # launch separate chunks
        for chunk in $(seq 1 $nchunk); do
	    chunk_obs_list=${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt
	    outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
	    errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	    for evenodd in even odd; do
		for pol in ${pol_list[@]}; do

    #Make the qsub command given the input parameters.
    echo "#!/bin/bash" > ./chunk_int_run_command.sh
    echo "#PBS -r y" >> ./chunk_int_run_command.sh
    echo "#PBS -N int" >> ./chunk_int_run_command.sh
    echo "#PBS -l mem=$mem" >> ./chunk_int_run_command.sh
    echo "#PBS -l walltime=${wallclock_time}" >> ./chunk_int_run_command.sh
    echo "#PBS -l ncpus=$ncpus" >> ./chunk_int_run_command.sh
    echo "#PBS -l software=idl" >> ./chunk_int_run_command.sh
    echo "#PBS -l jobfs=0MB" >> ./chunk_int_run_command.sh
    echo "#PBS -l wd" >> ./chunk_int_run_command.sh
    echo "#PBS -q express" >> ./chunk_int_run_command.sh
    echo "#PBS -P ru04" >> ./chunk_int_run_command.sh
    echo "#PBS -e ${outfile}" >> ./chunk_int_run_command.sh
    echo "#PBS -o ${errfile}" >> ./chunk_int_run_command.sh

    echo "module load idl/8.4"  >> ./chunk_int_run_command.sh

    #Create a name for the obs txt file based off of inputs
    evenoddpol_file_paths=${FHDdir}/Healpix/${version}_int_chunk${chunk}_${evenodd}${pol}_list.txt
    #clear old file paths
    echo "rm $evenoddpol_file_paths"  >> ./chunk_int_run_command.sh

    #***Fill the obs text file with the obsids to integrate
    echo "nobs=0"  >> ./chunk_int_run_command.sh
    echo "while read line"  >> ./chunk_int_run_command.sh
    echo "do"  >> ./chunk_int_run_command.sh
        echo "evenoddpol_file=${FHDdir}/Healpix/${line}_${evenodd}_cube${pol}.sav"  >> ./chunk_int_run_command.sh
        echo "echo \$evenoddpol_file >> $evenoddpol_file_paths"  >> ./chunk_int_run_command.sh
        echo "((nobs++))"  >> ./chunk_int_run_command.sh
    echo "done < \"$chunk_obs_list\""  >> ./chunk_int_run_command.sh
    #***

    #***If the integration has been split up into chunks, name the save file specifically off of that.
    echo "if [ \"$chunk\" -gt \"0\" ]; then"  >> ./chunk_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_int_chunk${chunk}_${evenodd}_cube${pol}.sav"  >> ./chunk_int_run_command.sh
    echo "else"  >> ./chunk_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_${evenodd}_cube${pol}.sav"  >> ./chunk_int_run_command.sh
    echo "fi"  >> ./chunk_int_run_command.sh
    #***

    #***Run the integration IDL script
    echo "idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncpus -e integrate_healpix_cubes -args \"$evenoddpol_file_paths\" \"\$save_file_evenoddpol\" &"  >> ./chunk_int_run_command.sh
    #***

    job="$(qsub ./chunk_int_run_command.sh)"
                    if [ "$chunk" == 0 ]; then
                    if [ "$evenodd" == "even" ]; then
                        if [ "$pol" == "XX" ]; then
                            id_list=${job}
                            else id_list="${id_list} -e ${job}"
                        fi
                        else id_list="${id_list} -e ${job}"
                    fi
                    else id_list="${id_list} -e ${job}"
                    fi
                    #else id_list="${id_list} -e ${job}"
          
 
		done
	    done
	    echo Combined_obs_${version}_int_chunk${chunk} >> $sub_cubes_list # trick it into finding our sub cubes
        done

    output="$(qstat | grep -e ${id_list} | wc -l)"
    while [ "$output" -ne "0" ]; do
        sleep 60
        output="$(qstat | grep -e ${id_list} | wc -l)"
    done



        # master integrator
        chunk=0
        outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
        errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	for evenodd in even odd; do
	    for pol in ${pol_list[@]}; do
                #message=$(sbatch --dependency=afterok:$idlist --mem=$mem -t $wallclock_time -n $ncores --export=file_path_cubes=$FHDdir,obs_list_path=$sub_cubes_list,version=$version,chunk=$chunk,ncores=$ncores,legacy=$legacy,evenodd=$evenodd,pol=$pol -e $errfile -o $outfile ${PSpath}../pipeline_scripts/bash_scripts/ozstar/integrate_slurm_job.sh)
        	#message=($message)
		#if [[ "$evenodd" = "even" ]] && [[ "$pol" = "XX" ]]; then idlist_master=${message[3]}; else idlist_master=${idlist_master}:${message[3]}; fi

    #Make the qsub command given the input parameters.
    echo "#!/bin/bash" > ./master_int_run_command.sh
    echo "#PBS -r y" >> ./master_int_run_command.sh
    echo "#PBS -N int" >> ./master_int_run_command.sh
    echo "#PBS -l mem=$mem" >> ./master_int_run_command.sh
    echo "#PBS -l walltime=${wallclock_time}" >> ./master_int_run_command.sh
    echo "#PBS -l ncpus=$ncpus" >> ./master_int_run_command.sh
    echo "#PBS -l software=idl" >> ./master_int_run_command.sh
    echo "#PBS -l jobfs=0MB" >> ./master_int_run_command.sh
    echo "#PBS -l wd" >> ./master_int_run_command.sh
    echo "#PBS -q express" >> ./master_int_run_command.sh
    echo "#PBS -P ru04" >> ./master_int_run_command.sh
    echo "#PBS -e ${outfile}" >> ./master_int_run_command.sh
    echo "#PBS -o ${errfile}" >> ./master_int_run_command.sh

    echo "module load idl/8.4"  >> ./master_int_run_command.sh

    #Create a name for the obs txt file based off of inputs
    evenoddpol_file_paths=${FHDdir}/Healpix/${version}_int_chunk${chunk}_${evenodd}${pol}_list.txt
    #clear old file paths
    echo "rm $evenoddpol_file_paths"  >> ./master_int_run_command.sh

    #***Fill the obs text file with the obsids to integrate
    echo "nobs=0"  >> ./master_int_run_command.sh
    echo "while read line"  >> ./master_int_run_command.sh
    echo "do"  >> ./master_int_run_command.sh
        echo "evenoddpol_file=${FHDdir}/Healpix/\${line}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
        echo "echo \$evenoddpol_file >> $evenoddpol_file_paths"  >> ./master_int_run_command.sh
        echo "((nobs++))"  >> ./master_int_run_command.sh
    echo "done < $sub_cubes_list"  >> ./master_int_run_command.sh
    #***

    #***If the integration has been split up into chunks, name the save file specifically off of that.
    echo "if [ \"$chunk\" -gt \"0\" ]; then"  >> ./master_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_int_chunk${chunk}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
    echo "else"  >> ./master_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
    echo "fi"  >> ./master_int_run_command.sh
    #***

    #***Run the integration IDL script
    echo "idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncpus -e integrate_healpix_cubes -args \"$evenoddpol_file_paths\" \"\$save_file_evenoddpol\" &"  >> ./master_int_run_command.sh
    #***

    job="$(qsub ./master_int_run_command.sh)"
                    if [ "$evenodd" == "even" ]; then
                        if [ "$pol" == "XX" ]; then
                            id_list=${job}
                            else id_list="${id_list} -e ${job}"
                        fi
                        else id_list="${id_list} -e ${job}"
                    fi


	    done
	done
    output="$(qstat | grep -e ${id_list} | wc -l)" 
    while [ "$output" -ne "0" ]; do
        sleep 60
        output="$(qstat | grep -e ${id_list} | wc -l)"
    done

    else

        # Just one integrator
        mv ${FHDdir}/Healpix/${version}_int_chunk1.txt ${FHDdir}/Healpix/${version}_int_chunk0.txt
        chunk=0
        chunk_obs_list=${FHDdir}/Healpix/${version}_int_chunk${chunk}.txt
        outfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_out.log
        errfile=${FHDdir}/Healpix/${version}_int_chunk${chunk}_err.log
	for evenodd in even odd; do
	    for pol in ${pol_list[@]}; do
                #message=$(sbatch ${hold_str} --mem=$mem -t ${wallclock_time} -n ${ncores} --export=file_path_cubes=$FHDdir,obs_list_path=$chunk_obs_list,version=$version,chunk=$chunk,ncores=$ncores,legacy=$legacy,evenodd=$evenodd,pol=$pol -e $errfile -o $outfile ${PSpath}../pipeline_scripts/bash_scripts/ozstar/integrate_slurm_job.sh)
       		#message=($message)
		#if [[ "$evenodd" = "even" ]] && [[ "$pol" = "XX" ]]; then idlist_int=${message[3]}; else idlist_int=${idlist_int}:${message[3]}; fi

    #Make the qsub command given the input parameters.
    echo "#!/bin/bash" > ./master_int_run_command.sh
    echo "#PBS -r y" >> ./master_int_run_command.sh
    echo "#PBS -N int" >> ./master_int_run_command.sh
    echo "#PBS -l mem=$mem" >> ./master_int_run_command.sh
    echo "#PBS -l walltime=${wallclock_time}" >> ./master_int_run_command.sh
    echo "#PBS -l ncpus=$ncpus" >> ./master_int_run_command.sh
    echo "#PBS -l software=idl" >> ./master_int_run_command.sh
    echo "#PBS -l jobfs=0MB" >> ./master_int_run_command.sh
    echo "#PBS -l wd" >> ./master_int_run_command.sh
    echo "#PBS -q express" >> ./master_int_run_command.sh
    echo "#PBS -P ru04" >> ./master_int_run_command.sh
    echo "#PBS -e ${outfile}" >> ./master_int_run_command.sh
    echo "#PBS -o ${errfile}" >> ./master_int_run_command.sh

    echo "module load idl/8.4"  >> ./master_int_run_command.sh

    #Create a name for the obs txt file based off of inputs
    evenoddpol_file_paths=${FHDdir}/Healpix/${version}_int_chunk${chunk}_${evenodd}${pol}_list.txt
    #clear old file paths
    echo "rm $evenoddpol_file_paths"  >> ./master_int_run_command.sh
        
    #***Fill the obs text file with the obsids to integrate
    echo "nobs=0"  >> ./master_int_run_command.sh
    echo "while read line"  >> ./master_int_run_command.sh
    echo "do"  >> ./master_int_run_command.sh
        echo "evenoddpol_file=${FHDdir}/Healpix/\${line}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
        echo "echo \$evenoddpol_file >> $evenoddpol_file_paths"  >> ./master_int_run_command.sh
        echo "((nobs++))"  >> ./master_int_run_command.sh
    echo "done < \"$chunk_obs_list\""  >> ./master_int_run_command.sh
    #***

    #***If the integration has been split up into chunks, name the save file specifically off of that.
    echo "if [ \"$chunk\" -gt \"0\" ]; then"  >> ./master_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_int_chunk${chunk}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
    echo "else"  >> ./master_int_run_command.sh
        echo "save_file_evenoddpol=$FHDdir/Healpix/Combined_obs_${version}_${evenodd}_cube${pol}.sav"  >> ./master_int_run_command.sh
    echo "fi"  >> ./master_int_run_command.sh
    #***

    #***Run the integration IDL script
    echo "idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncpus -e integrate_healpix_cubes -args \"$evenoddpol_file_paths\" \"\$save_file_evenoddpol\""  >> ./master_int_run_command.sh
    #***

    job="$(qsub ./master_int_run_command.sh)"
                    if [ "$evenodd" == "even" ]; then
                        if [ "$pol" == "XX" ]; then
                            id_list=${job}
                            else id_list="${id_list} -e ${job}"
                        fi
                        else id_list="${id_list} -e ${job}"
                    fi

	    done
	done
    output="$(qstat | grep -e ${id_list} | wc -l)"
    while [ "$output" -ne "0" ]; do
        sleep 60
        output="$(qstat | grep -e ${id_list} | wc -l)"
    done


    fi
else
    echo "Running only ps code" # Just PS if flag has been set
fi
outfile=${FHDdir}/ps/logs/${version}_ps_out
errfile=${FHDdir}/ps/logs/${version}_ps_err


mkdir -p ${FHDdir}/ps/logs
mkdir -p ${FHDdir}/ps/data/uvf_cubes
mkdir -p ${FHDdir}/ps/data/kspace_cubes
mkdir -p ${FHDdir}/ps/data/2d_binning
mkdir -p ${FHDdir}/ps/data/1d_binning

mkdir -p ${FHDdir}/ps/plots/slices
mkdir -p ${FHDdir}/ps/plots/2d_binning
mkdir -p ${FHDdir}/ps/plots/1d_binning

if [[ "$n_pol" -eq 2 ]]; then pol_list=( xx yy );fi
if [[ "$n_pol" -eq 4 ]]; then pol_list=( xx yy xy yx );fi

evenodd_list=( even odd )
cube_type_list=( weights dirty model res )


for cube_type in ${cube_type_list[@]}; do
for pol in ${pol_list[@]}; do
    echo $pol
    for evenodd in ${evenodd_list[@]}; do
        echo $evenodd
        evenodd_initial="$(echo $evenodd | head -c 1)"
        echo $cube_type
	    if [ "$cube_type" == "dirty" ]; then 
                   if [ "$evenodd" == "even" ]; then
                        if [ "$pol" == "xx" ]; then            
echo "in wait loop"
    output="$(qstat | grep -e ${id_list} | wc -l)"
    while [ "$output" -ne "0" ]; do
        sleep 60
        output="$(qstat | grep -e ${id_list} | wc -l)"
    done
                        fi
                   fi
            fi


            #message=$(sbatch ${hold_str_cubes} --mem=$mem -t ${wallclock_time} -n ${ncores} --export=file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,ncores=$ncores,cube_type=$cube_type,pol=$pol,evenodd=$evenodd,image_filter_name=$tukey_filter -e ${errfile}_${pol}_${evenodd}_${cube_type}.log -o ${outfile}_${pol}_${evenodd}_${cube_type}.log -J PS_${pol}${evenodd_initial}_${cube_type} ${PSpath}${pipe_path}PS_list_slurm_job.sh)
            #message=($message)

    #Make the qsub command given the input parameters.
    echo "#!/bin/bash" > ./ps_run_command.sh
    echo "#PBS -r y" >> ./ps_run_command.sh
    echo "#PBS -N PS_${pol}${evenodd_initial}_${cube_type}" >> ./ps_run_command.sh
    echo "#PBS -l mem=$mem" >> ./ps_run_command.sh
    echo "#PBS -l walltime=${wallclock_time}" >> ./ps_run_command.sh
    echo "#PBS -l ncpus=$ncpus" >> ./ps_run_command.sh
    echo "#PBS -l software=idl" >> ./ps_run_command.sh
    echo "#PBS -l jobfs=1MB" >> ./ps_run_command.sh
    echo "#PBS -l wd" >> ./ps_run_command.sh
    echo "#PBS -q express" >> ./ps_run_command.sh
    echo "#PBS -P ru04" >> ./ps_run_command.sh
    echo "#PBS -e ${errfile}_${pol}_${evenodd}_${cube_type}.log" >> ./ps_run_command.sh
    echo "#PBS -o ${outfile}_${pol}_${evenodd}_${cube_type}.log" >> ./ps_run_command.sh

    echo "module load idl/8.4"  >> ./ps_run_command.sh

#***Get the obsid file paths
echo "nobs=0" >> ./ps_run_command.sh
echo "while read line" >> ./ps_run_command.sh
echo "do" >> ./ps_run_command.sh
  echo "((nobs++))" >> ./ps_run_command.sh
echo "done < \"$integrate_list\"" >> ./ps_run_command.sh

echo "input_file=${FHDdir}/"  >> ./ps_run_command.sh

echo "arg_string=\"\${input_file} ${version} ${cube_type} ${pol} ${evenodd}\""   >> ./ps_run_command.sh

echo "PSpath=\$(idl -e 'print,rootdir(\"eppsilon\")')"   >> ./ps_run_command.sh
echo "pipe_path=\"../pipeline_scripts/eppsilon_IDL_wrappers/\""   >> ./ps_run_command.sh

echo "idl -IDL_DEVICE ps -IDL_CPU_TPOOL_NTHREADS $ncpus -e raijin_ps_job -args \$arg_string"  >> ./ps_run_command.sh


    job="$(qsub ./ps_run_command.sh)"
                if [ "$cube_type" == "weights" ]; then
                    if [ "$evenodd" == "even" ]; then
                        if [ "$pol" == "xx" ]; then
                            id_list=${job}
                            else id_list="${id_list} -e ${job}"
                        fi
                        else id_list="${id_list} -e ${job}"
                    fi
                 else id_list="${id_list} -e ${job}"
                 fi

        done
    done
done

#final plots
#if [[ -n ${tukey_filter} ]]; then plot_walltime=1:00:00; else plot_walltime=00:20:00; fi
#sbatch --dependency=afterok:$id_list --mem=$mem -t ${wallclock_time} -n ${ncores} --export=file_path_cubes=$FHDdir,obs_list_path=$integrate_list,version=$version,ncores=$ncores,image_filter_name=$tukey_filter -e ${errfile}_plots.log -o ${outfile}_plots.log -J PS_plots ${PSpath}${pipe_path}PS_list_slurm_job.sh
