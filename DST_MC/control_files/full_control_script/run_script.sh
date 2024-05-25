#!/bin/bash

#todo : the script should be copied to the corresponding control folder

# note : full directory

control_folder_directory=`pwd`
running_argument=$1
topic_focus=$2
data_type=$3
file_folder_directory=$4
number_of_jobs=$5
special_name=$6

dir_to_runXY_stability="/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC/runXY_stability"
dir_to_eta_hist_merge="/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC/eta_hist_merge"
dir_to_ana_map="/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC/ana_map_folder"
dir_to_DST_MC="/sphenix/user/ChengWei/INTT/INTT_dNdeta_repo/DST_MC"
used_map="ana_map_v1.h"

control_dir_to_data_type=${control_folder_directory}/for_${data_type}

condor_monitor_func() {
    N_CW_condor_job=100;
    time while [[ $N_CW_condor_job -ne 1 ]]; 
    do 
        N_CW_condor_job=`condor_q | grep ecie9969 | wc -l`

        printf '=%.0s' $(seq $N_CW_condor_job); 
        echo $N_CW_condor_job;
        echo `condor_q | grep "Total for ecie9969"`
        sleep 2; 
    done;
}

if [[ "$running_argument" != "run_confirm" && \
      "$running_argument" != "remove_sub_sh" && \
      "$running_argument" != "clear_control_file" && \
      "$running_argument" != "clear_ana_file_folder" && \
      "$running_argument" != "help" && \
      "$running_argument" != "list" && \
      "$running_argument" != "condor_monitor" && \
      "$running_argument" != "create_new" && \
      "$running_argument" != "condor_submit" && \
      "$running_argument" != "get_merged_result" && \
      "$running_argument" != "full_clear" ]]; then
    echo "Incorrect running_argument, exiting the operation"
    exit 1
fi

# note : ===============================================================================================================================================================================
if [[ "${running_argument}" == "help" ]]; then 
    echo
    echo "please do: sh run_script.sh [running_argument] [topic_focus] [data_type] [file_folder_directory] [number_of_jobs/special_name]"
    echo
    exit 1
fi

# note : ===============================================================================================================================================================================
if [[ "${running_argument}" == "list" ]]; then
    echo '====================================================================================================================================';
    echo '| running_argument (5): create_new                                                                                                 |';
    echo '| running_argument (4): condor_submit / run_confirm / get_merged_result / full_clear / clear_ana_file_folder                       |';
    echo '| running_argument (3): clear_control_file                                                                                         |';
    echo '| running_argument (1): help / list / condor_monitor / remove_sub_sh                                                               |';
    echo '|                                                                                                                                  |';
    echo '| 2. topic_focus: avg_vtxXY / evt_vtxZ / evt_vtxXY / evt_tracklet / evt_tracklet_inner_phi_rotate                                  |';
    echo '|                                                                                                                                  |';
    echo '| 3. data_type: data / MC                                                                                                          |';
    echo '|                                                                                                                                  |';
    echo '| 4. File_directory                                                                                                                |';
    echo '|                                                                                                                                  |';
    echo '| 5. number_of_jobs                                                                                                                |';
    echo '|                                                                                                                                  |';
    echo '| 6. The_sub_folder_name_for_special                                                                                               |';
    echo '|                                      1.             2.            3.                4.                  5.             (6.)      |';
    echo '| please do: sh run_script.sh [running_argument] [topic_focus] [data_type] [file_folder_directory] [number_of_jobs] [special_name] |';
    echo '====================================================================================================================================';
    exit 1
fi

# note : ===============================================================================================================================================================================
if [[ "${running_argument}" == "condor_monitor" ]]; then 
    condor_monitor_func
    exit 1
fi

# note : ===============================================================================================================================================================================
if [[ "${running_argument}" == "remove_sub_sh" ]]; then 
    rm for_MC/*/run_root_sub_*.sh
    rm for_data/*/run_root_sub_*.sh

    exit 1
fi

# note : ===============================================================================================================================================================================
if [[ ${#running_argument} -ne 0 ]]; then
    echo the running_argument: ${running_argument}
fi

if [[ "$data_type" != "data" && "$data_type" != "MC" ]]; then
    echo incorrect data_type, exit the operation
    exit 1
fi

if [[ "$topic_focus" != "avg_vtxXY" && \
      "$topic_focus" != "evt_vtxZ" && \
      "$topic_focus" != "evt_vtxXY" && \
      "$topic_focus" != "evt_tracklet" && \
      "$topic_focus" != "evt_tracklet_inner_phi_rotate" ]]; then
    echo "incorrect topic_focus, please try [avg_vtxXY], [evt_vtxZ], [evt_vtxXY], [evt_tracklet], [evt_tracklet_inner_phi_rotate]" 
    echo exit the program
    exit 1
fi

if [[ ${#control_folder_directory} -ne 0 ]]; then
    echo the control_folder_directory: ${control_folder_directory}
fi

if [[ ${#file_folder_directory} -ne 0 ]]; then
    echo the file_folder_directory: ${file_folder_directory}
fi

if [ "$data_type" == "MC" ]; then 
    data_type_id=1;
else 
    data_type_id=0;
fi

# note : ===============================================================================================================================================================================
file_directory_in_used_map=`grep "${data_type}_input_directory" ${dir_to_ana_map}/${used_map} | awk '{print $4}' | tr -d '";'`
print_message1="| the one in the ${used_map}:"
print_message2="| the one passed in the $0:"
max_length=$((${#print_message1} > ${#print_message2} ? ${#print_message1} : ${#print_message2}))
echo "================================================================================================================================|"
echo "| double confirmation of the input file directory                                                                               |"
printf "%-${max_length}s %-${max_length}s\n" "$print_message1" "${file_directory_in_used_map}"
printf "%-${max_length}s %-${max_length}s\n" "$print_message2" "${file_folder_directory}"
echo "================================================================================================================================|"
if [[ "${file_directory_in_used_map}" != "${file_folder_directory}" ]]; then
    echo '!!!!!'the defined file directory in the ${used_map} is different from what you passed to the $0, 'please take a look!'
    exit 1
fi

# note : ===============================================================================================================================================================================
N_found_special_tag=$(( `grep "special_tag" ${dir_to_ana_map}/${used_map} | wc -l` + `grep "special_tag" ${dir_to_DST_MC}/*.h | wc -l` + `grep "special_tag" ${control_dir_to_data_type}/${topic_focus}/*_mother.C | wc -l` ))

# note : if we have the sixth argument for the speical tests
if [ ! -z "$6" ]; then
    echo '!!!!! it seems to be a special run, with the name:' $special_name
    read -r -p "Shall we move on? " response
    if [[ "$response" != "yes" ]]; then
        echo "Operation aborted."
        exit 1
    fi
    
    # note : check the name of the special run
    if [[ `echo $special_name | grep ${topic_focus} | wc -l` -ne 1 ]]; then
        echo the $special_name seems not to be the branch of the ${topic_focus}, please modify the way how you name it.
        exit 1
    fi

    # note : check the speical tag
    if [[ ${N_found_special_tag} -eq 0 ]]; then 
        echo
        echo '!!!! we do not see any special_tag in the following files, pleaes check'
        echo ${dir_to_ana_map}/${used_map}
        echo ${dir_to_DST_MC}/*.h | sed -e "s/ /\n/g"
        echo ${control_dir_to_data_type}/${topic_focus}/*_mother.C | sed -e "s/ /\n/g"
        echo 
        exit 1
    else 
        echo found ${N_found_special_tag} in the followings, please confirm
        grep "special_tag" ${dir_to_ana_map}/${used_map}
        grep "special_tag" ${dir_to_DST_MC}/*.h
        grep "special_tag" ${control_dir_to_data_type}/${topic_focus}/*_mother.C
        echo 
        read -r -p "Are they looked fine ? " response
        if [[ "$response" != "yes" ]]; then
            echo "Operation aborted."
            exit 1
        fi
    fi 

    dir_to_file_sub_folder=${file_folder_directory}/${special_name}
else 
    if [[ ${N_found_special_tag} -ne 0 ]]; then
        echo found ${N_found_special_tag} in the followings, 'but we are in the normal mode!!! please confirm'
        grep "special_tag" ${dir_to_ana_map}/${used_map}
        grep "special_tag" ${dir_to_DST_MC}/*.h
        grep "special_tag" ${control_dir_to_data_type}/${topic_focus}/*_mother.C
        exit 1
    fi

    dir_to_file_sub_folder=${file_folder_directory}/${topic_focus}
fi

# note : ===============================================================================================================================================================================
if [[ "${running_argument}" == "create_new" ]]; then

    if [[ ${#control_folder_directory} == 0 || ${#file_folder_directory} == 0 ]]; then
        echo No control file directory in the input or the file_folder_directory, exit the operation
        exit 1
    fi

    if [ -z "$5" ]; then 
        echo "No number of cores in the argument!"
        exit 1
    fi

    echo "there seems to be a new ${data_type} file! Creating the folders for it under the file_directory"
    mkdir ${control_dir_to_data_type}/${topic_focus}/condor_outputs
    mkdir $dir_to_file_sub_folder
    mkdir $dir_to_file_sub_folder/CW_log
    
    if [[ "$topic_focus" == "evt_vtxZ" || "$topic_focus" == "evt_vtxXY" || "$topic_focus" == "evt_tracklet" || "$topic_focus" == "evt_tracklet_inner_phi_rotate" ]]; then
        mkdir $dir_to_file_sub_folder/complete_file
    fi

    cd ${control_dir_to_data_type}/${topic_focus}
    if [ -f "run_condor_job_mother_copied.job" ]; then
        rm run_condor_job_mother_copied.job
    fi
    cp run_condor_job_mother.job run_condor_job_mother_copied.job

    if [ -f "run_root_mother_copied.sh" ]; then
        rm run_root_mother_copied.sh
    fi
    cp run_root_mother.sh run_root_mother_copied.sh

    if [ -f "run_condor_copied.sh" ]; then
        rm run_condor_copied.sh
    fi
    cp run_condor.sh run_condor_copied.sh 

    sed -i "s|control_file_to_be_updated|${control_dir_to_data_type}|g" run_condor_job_mother_copied.job
    sed -i "s|condor_output_to_be_updated|${file_folder_directory}|g" run_condor_job_mother_copied.job
    
    sed -i "s|control_file_to_be_updated|${control_dir_to_data_type}|g" run_root_mother_copied.sh
    
    sed -i "s/number_of_jobs_to_be_updated/${number_of_jobs}/g" run_condor_copied.sh
    
    cd ${control_folder_directory}

    exit 1

# note : ===============================================================================================================================================================================
elif [[ "${running_argument}" == "condor_submit" ]]; then
    echo "the script is going to submit the condor jobs"

    cd ${control_dir_to_data_type}/${topic_focus}
    echo '-----> check number of jobs!: ~~'"`grep \"i<=\" run_condor_copied.sh`"'~~'  
    read -r -p "The number of jobs has been printed out above, are you sure about that? Type 'yes' to proceed: " response
    if [[ "$response" != "yes" ]]; then
        echo "Operation aborted."
        exit 1
    fi

    if [[ `ls -1 $dir_to_file_sub_folder/CW_log | wc -l` -ne 0 ]]; then 
        rm $dir_to_file_sub_folder/CW_log/condor*
    fi

    if [[ `ls -1 ${control_dir_to_data_type}/${topic_focus}/run_root_sub_*.sh | wc -l` -ne 0 ]]; then 
        rm ${control_dir_to_data_type}/${topic_focus}/run_root_sub_*.sh
    fi  

    if [[ "$topic_focus" == "avg_vtxXY" ]]; then
        if [[ `ls -d -1 $dir_to_file_sub_folder/runXY_* | wc -l` -ne 0 ]]; then
            rm -r $dir_to_file_sub_folder/runXY_*
        fi
    fi

    if [[ "$topic_focus" == "evt_vtxZ" ]]; then
        if [[ `ls -d -1 $dir_to_file_sub_folder/complete_file/evtZ_* | wc -l` -ne 0 ]]; then
            rm -r $dir_to_file_sub_folder/complete_file/evtZ_*
        fi
    fi

    if [[ "$topic_focus" == "evt_vtxXY" ]]; then
        if [[ `ls -d -1 $dir_to_file_sub_folder/complete_file/evtXY_* | wc -l` -ne 0 ]]; then
            rm -r $dir_to_file_sub_folder/complete_file/evtXY_*
        fi
    fi

    if [[ "$topic_focus" == "evt_tracklet" || "$topic_focus" == "evt_tracklet_inner_phi_rotate" ]]; then
        if [[ `ls -d -1 $dir_to_file_sub_folder/complete_file/evtTracklet_* | wc -l` -ne 0 ]]; then
            rm -r $dir_to_file_sub_folder/complete_file/evtTracklet_*
        fi
    fi

    sh run_condor_copied.sh
    cd ${control_folder_directory}

    condor_monitor_func

    exit 1

# note : ===============================================================================================================================================================================
elif [[ "${running_argument}" == "get_merged_result" ]]; then 
    source /opt/sphenix/core/bin/sphenix_setup.sh -n ana.410

    if [[ "$topic_focus" == "avg_vtxXY" ]]; then
        echo running the merged file ana for $topic_focus
        cd ${dir_to_runXY_stability}
        ./run_avg_vtxXY_merge_result $data_type_id $dir_to_file_sub_folder
        cd ${control_folder_directory}  
    fi

    if [[ "$topic_focus" == "evt_vtxZ" ]]; then
        echo running the merged file ana for $topic_focus
        root -l -b -q ${control_dir_to_data_type}/${topic_focus}/MakeEvtZPlots.cpp
    fi

    if [[ "$topic_focus" == "evt_vtxXY" ]]; then
        echo running the merged file ana for $topic_focus
        mkdir $dir_to_file_sub_folder/complete_file/file_merged_folder
        cd $dir_to_file_sub_folder/complete_file
        ls evtXY_*/evt_XY_tree.root > file_list.txt
        cd ${control_dir_to_data_type}/${topic_focus}
        root -l -b -q make_evtXY_plot.cpp
        cd ${control_folder_directory}

        exit 1
    fi    

    if [[ "$topic_focus" == "evt_tracklet" || "$topic_focus" == "evt_tracklet_inner_phi_rotate" ]]; then
        echo running the merged file ana for $topic_focus
        cd $dir_to_file_sub_folder/complete_file
        
        if [[ ! -d "merged_file_folder" ]]; then mkdir merged_file_folder; fi
        
        if [[ ! -f "file_list_full.txt" ]]; then
            ls $dir_to_file_sub_folder/complete_file/evtTracklet_*/INTT_final_hist_info.root > file_list_full.txt
        fi

        if [[ "$topic_focus" == "evt_tracklet" ]]; then
            if [[ ! -f "file_list_full_shuf.txt" ]]; then shuf "file_list_full.txt" > "file_list_full_shuf.txt"; fi

            if [[ ! -f "file_list_half_corr.txt" || ! -f "file_list_half_test.txt" ]]; then
                Nhalf_evtTracklet=$(( `grep ".root" file_list_full_shuf.txt | wc -l` / 2 ))
                head -n $Nhalf_evtTracklet "file_list_full_shuf.txt" > "file_list_half_corr.txt"
                tail -n $(( `grep ".root" file_list_full_shuf.txt | wc -l` - Nhalf_evtTracklet )) "file_list_full_shuf.txt" > "file_list_half_test.txt"
            fi
        else # note : for the "evt_tracklet_inner_phi_rotate"
            echo for the evt_tracklet_inner_phi_rotate, we copied the file from evt_tracklet
            cp ${file_folder_directory}/evt_tracklet/complete_file/file_list_half_corr.txt $dir_to_file_sub_folder/complete_file
            cp ${file_folder_directory}/evt_tracklet/complete_file/file_list_half_test.txt $dir_to_file_sub_folder/complete_file
            
            sed -i "s/evt_tracklet/${topic_focus}/g" file_list_half_corr.txt
            sed -i "s/evt_tracklet/${topic_focus}/g" file_list_half_test.txt

        fi        

        echo N_lines of file_list_full.txt: `grep ".root" file_list_full.txt | wc -l`
        if [[ -f "file_list_full_shuf.txt" ]]; then echo N_lines of file_list_full_shuf.txt: `grep ".root" file_list_full_shuf.txt | wc -l`; fi
        echo N_lines of file_list_half_corr.txt: `grep ".root" file_list_half_corr.txt | wc -l`
        echo N_lines of file_list_half_test.txt: `grep ".root" file_list_half_test.txt | wc -l`

        sort file_list_half_corr.txt > sorted_corr.txt
        sort file_list_half_test.txt > sorted_test.txt

        echo "================================================================================================================================"
        echo "======Check identical lines====================================================================================================="
        echo "================================================================================================================================"
        comm -12 sorted_corr.txt sorted_test.txt
        echo "================================================================================================================================"
        echo "======Do you got something above ?=============================================================================================="
        echo "================================================================================================================================"

        if [[ `comm -12 sorted_corr.txt sorted_test.txt | wc -l` -ne 0 ]]; then
            echo '!!!! something wrong in split the files, please take a look. Abort the program now'
            exit 1
        fi

        cd merged_file_folder
        if [[ -f "merged_hist_half_corr.root" && -f "merged_hist_half_test.root" && -f "merged_hist_full.root" ]]; then 
            root -l -b -q ${control_dir_to_data_type}/${topic_focus}/merged_tracklet_ana.C
        else 
            echo there is no final merged files found, so submit the job in the background
            cd ${dir_to_eta_hist_merge}
            # echo '!!!! you are in the directory:' `pwd`
            # echo '!!!! for the current safety, please run the followings manually'
            # echo "./run_merge $dir_to_file_sub_folder/complete_file file_list_full.txt merged_hist_full.root"
            # echo "./run_merge $dir_to_file_sub_folder/complete_file file_list_half_corr.txt merged_hist_half_corr.root"
            # echo "./run_merge $dir_to_file_sub_folder/complete_file file_list_half_test.txt merged_hist_half_test.root"

            if [[ -f "$dir_to_file_sub_folder/complete_file/recod_merged_hist_full.txt" ]]; then rm $dir_to_file_sub_folder/complete_file/recod_merged_hist_full.txt; fi
            if [[ -f "$dir_to_file_sub_folder/complete_file/recod_merged_hist_half_corr.txt" ]]; then rm $dir_to_file_sub_folder/complete_file/recod_merged_hist_half_corr.txt; fi
            if [[ -f "$dir_to_file_sub_folder/complete_file/recod_merged_hist_half_test.txt" ]]; then rm $dir_to_file_sub_folder/complete_file/recod_merged_hist_half_test.txt; fi

            ./run_merge $dir_to_file_sub_folder/complete_file file_list_full.txt merged_hist_full.root > $dir_to_file_sub_folder/complete_file/recod_merged_hist_full.txt 2>&1 &
            ./run_merge $dir_to_file_sub_folder/complete_file file_list_half_corr.txt merged_hist_half_corr.root > $dir_to_file_sub_folder/complete_file/recod_merged_hist_half_corr.txt 2>&1 &
            ./run_merge $dir_to_file_sub_folder/complete_file file_list_half_test.txt merged_hist_half_test.root > $dir_to_file_sub_folder/complete_file/recod_merged_hist_half_test.txt 2>&1 &

            N_merge_finish_tag=100;
            time while [[ $N_merge_finish_tag -ne 3 ]]; 
            do 
                N_merge_finish_tag=`grep "the merged file is closed" $dir_to_file_sub_folder/complete_file/recod_*.txt | wc -l`

                echo N run_merge jobs running: $(( `ps aux | grep run_merge | wc -l` - 1 ))
                printf '=%.0s' $(seq $N_merge_finish_tag); 
                echo $N_merge_finish_tag;
                sleep 2; 
            done;

            if [[ $N_merge_finish_tag == 3 ]]; then 
                pkill run_merge
                root -l -b -q ${control_dir_to_data_type}/${topic_focus}/merged_tracklet_ana.C
            fi
        fi

        exit 1
        # cd ${control_dir_to_data_type}/${topic_focus}
        # root -l -b -q ${control_dir_to_data_type}/${topic_focus}/MakeEvtZPlots.cpp
    fi

# note : ===============================================================================================================================================================================
elif [[ "${running_argument}" == "full_clear" ]]; then 
    if [[ ${#file_folder_directory} -eq 0 ]]; then
        echo no file directory input abort the program
        exit 1
    fi

    echo deleting the following items
    echo ${control_dir_to_data_type}/${topic_focus}/condor_outputs
    echo ${control_dir_to_data_type}/${topic_focus}/run_condor_job_mother_copied.job
    echo ${control_dir_to_data_type}/${topic_focus}/run_root_mother_copied.sh
    echo ${control_dir_to_data_type}/${topic_focus}/run_condor_copied.sh
    echo $dir_to_file_sub_folder

    read -r -p "Do you want to continue? Type 'yes' to proceed: " response
    if [[ "$response" != "yes" ]]; then
        echo "Operation aborted."
        exit 1
    fi

    rm -r ${control_dir_to_data_type}/${topic_focus}/condor_outputs
    rm ${control_dir_to_data_type}/${topic_focus}/run_condor_job_mother_copied.job
    rm ${control_dir_to_data_type}/${topic_focus}/run_root_mother_copied.sh
    rm ${control_dir_to_data_type}/${topic_focus}/run_condor_copied.sh
    rm -r $dir_to_file_sub_folder

    exit 1

# note : ===============================================================================================================================================================================
elif [[ "${running_argument}" == "clear_control_file" ]]; then 

    echo deleting the following items
    echo ${control_dir_to_data_type}/${topic_focus}/condor_outputs
    echo ${control_dir_to_data_type}/${topic_focus}/run_condor_job_mother_copied.job
    echo ${control_dir_to_data_type}/${topic_focus}/run_root_mother_copied.sh
    echo ${control_dir_to_data_type}/${topic_focus}/run_condor_copied.sh

    read -r -p "Do you want to continue? Type 'yes' to proceed: " response
    if [[ "$response" != "yes" ]]; then
        echo "Operation aborted."
        exit 1
    fi

    rm -r ${control_dir_to_data_type}/${topic_focus}/condor_outputs
    rm ${control_dir_to_data_type}/${topic_focus}/run_condor_job_mother_copied.job
    rm ${control_dir_to_data_type}/${topic_focus}/run_root_mother_copied.sh
    rm ${control_dir_to_data_type}/${topic_focus}/run_condor_copied.sh

    exit 1

# note : ===============================================================================================================================================================================
elif [[ "${running_argument}" == "clear_ana_file_folder" ]]; then 
    if [[ ${#file_folder_directory} -eq 0 ]]; then
        echo no file directory input abort the program
        exit 1
    fi

    echo deleting the following items
    echo $dir_to_file_sub_folder

    read -r -p "Do you want to continue? Type 'yes' to proceed: " response
    if [[ "$response" != "yes" ]]; then
        echo "Operation aborted."
        exit 1
    fi

    rm -r $dir_to_file_sub_folder

    exit 1

elif [[ "${running_argument}" == "run_confirm" ]]; then 
    echo N condor_out with finished: `grep "finished" $dir_to_file_sub_folder/CW_log/condor_*.out | wc -l`
    exit 1
else 
    echo '!!!!! Something wrong, you better check !!!!!'
fi

echo "This is the end of the script"
