#! /usr/global/bin/python

# Pipeline Wrapper Version (last revision's modified date)
last_modified = 'Last modified by Eric Earl, 6/25/2015'


# Add system paths for python modules while script is running.
import sys,os,numpy,argparse,csv,re
from scipy.io import loadmat

# Program Description
prog_descrip = """%(prog)s: use -h option for usage."""

def main(argv=sys.argv):
    pipeline_script_dir = os.path.dirname(argv[0])
    
    # Setting up the argument parsing for the script
    arg_parser = argparse.ArgumentParser(description=prog_descrip)
    
    arg_parser.add_argument('-a','--analysis-group-folder', nargs=1, action='store', required=True, type=os.path.abspath,
                            help=('Path to the analyses_v2 group folder containing subjects in subject list text file.'),
                            dest='group_folder'
                           )
    arg_parser.add_argument('-l','--subject-list', nargs=1, action='store', required=True, type=os.path.abspath,
                            help=('The subject list text file (which can also be used for the graph analysis GUI).'),
                            dest='subject_list'
                           )
    arg_parser.add_argument('-o','--output-csv', nargs=1, action='store', required=True, type=os.path.abspath,
                            help=('Path to output movement report CSV file.'),
                            dest='output_csv'
                           )
    arg_parser.add_argument('-d','--demo-file', nargs=1, action='store', required=False, type=os.path.abspath, default='',
                            help=('Optional demographics file.'),
                            dest='demo_file'
                           )
#    arg_parser.add_argument('-c','--demo-columns', nargs='*', action='store', required=False, type=int, default=[],
#                            help=('Optional demographics file columns to read.'),
#                            dest='demo_columns'
#                           )
    arg_parser.add_argument('-v','-version', action='version', version='%(prog)s: ' + last_modified,
                            help=("Return script's last modified date.")
                           )
    
    args = arg_parser.parse_args()
    
    analyses_v2_group_folder = args.group_folder[0] # example: /group_shares/FAIR_ADHD/CYA/analyses_v2/ADHD-HumanYouth-OHSU/FAIRPRE10_TR2pt5_RAD50pt0_SKIPSEC5pt0-FCON8-SEEDCORREL_1
    GAGUI_subject_list_file = args.subject_list[0] # A Graph Analysis GUI subject list text file containing a list
    output_csv = args.output_csv[0] # The path to and name of the output CSV file including the .csv extension
    
    if args.demo_file:
        demo_file_path = args.demo_file[0]
    
#    if args.demo_columns:
#        keep_columns = args.demo_columns
    
    subject_list = []
    subject_list_file = open(GAGUI_subject_list_file)
    for line in subject_list_file:
        subject_list.append(line.strip())
    subject_list_file.close()
    
    subject_dict = {}
    for subject in subject_list:
        try:
            prepost_exist_flag = True
            prepost_mat = loadmat(os.path.join(analyses_v2_group_folder,subject,'motion','FD.mat'))
        except:
            print 'No prepost MAT file found here: ' + os.path.join(analyses_v2_group_folder,subject,'motion','FD.mat')
            prepost_exist_flag = False
            
        try:
            power2014_exist_flag = True
            power2014_mat = loadmat(os.path.join(analyses_v2_group_folder,subject,'motion','power_2014_motion.mat'))
        except:
            print 'No power2014 MAT file found here: ' + os.path.join(analyses_v2_group_folder,subject,'motion','power_2014_motion.mat')
            power2014_exist_flag = False
            
        out_dict = {'prepost':[],'prepost_exist':prepost_exist_flag,'power2014':[],'power2014_exist':power2014_exist_flag}
        
        for fd_index in range(51):
            if prepost_exist_flag:
                prepost_dict = {}
                prepost_dict['frame_removal'] = prepost_mat['FD_data'][0,fd_index].__getitem__('frame_removal')[0,0]
                prepost_dict['format_string'] = str(prepost_mat['FD_data'][0,fd_index].__getitem__('format_string')[0,0])
                prepost_dict['skip'] = int(prepost_mat['FD_data'][0,fd_index].__getitem__('skip')[0,0])
                prepost_dict['total_frame_count'] = int(prepost_mat['FD_data'][0,fd_index].__getitem__('total_frame_count')[0,0])
                prepost_dict['remaining_frame_count'] = int(prepost_mat['FD_data'][0,fd_index].__getitem__('remaining_frame_count')[0,0])
                prepost_dict['epi_TR'] = float(prepost_mat['FD_data'][0,fd_index].__getitem__('epi_TR')[0,0])
                prepost_dict['FD_threshold'] = float(prepost_mat['FD_data'][0,fd_index].__getitem__('FD_threshold')[0,0])
                prepost_dict['remaining_seconds'] = float(prepost_mat['FD_data'][0,fd_index].__getitem__('remaining_seconds')[0,0])
                try:
                    prepost_dict['remaining_frame_mean_FD'] = float(prepost_mat['FD_data'][0,fd_index].__getitem__('remaining_frame_mean_FD')[0,0])
                except:
                    print 'No "remaining_frame_mean_FD" in ' + subject + ' prepost MAT file.'
                
                out_dict['prepost'].append(prepost_dict)
            
            if power2014_exist_flag:
                single_fd_row = []
                for dvar_index in range(51):
                    power2014_dict = {}
                    power2014_dict['frame_removal'] = power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('frame_removal')[0,0]
                    power2014_dict['format_string'] = str(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('format_string')[0,0])
                    power2014_dict['skip'] = int(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('skip')[0,0])
                    power2014_dict['total_frame_count'] = int(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('total_frame_count')[0,0])
                    power2014_dict['remaining_frame_count'] = int(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('remaining_frame_count')[0,0])
                    power2014_dict['epi_TR'] = float(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('epi_TR')[0,0])
                    power2014_dict['FD_threshold'] = float(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('FD_threshold')[0,0])
                    power2014_dict['remaining_seconds'] = float(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('remaining_seconds')[0,0])
                    try:
                        power2014_dict['remaining_frame_mean_FD'] = float(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('remaining_frame_mean_FD')[0,0])
                        power2014_dict['remaining_frame_mean_DVAR'] = float(power2014_mat['motion_data'][fd_index,dvar_index,0].__getitem__('remaining_frame_mean_DVAR')[0,0])
                    except:
                        print 'No "remaining_frame_mean_FD" or "remaining_frame_mean_DVAR" in ' + subject + ' power2014 MAT file.'
                    
                    single_fd_row.append(power2014_dict)
                    
                out_dict['power2014'].append(single_fd_row)
        
        subject_dict[subject] = out_dict
    
    if args.demo_file:
        demographic_reader = csv.DictReader(open(demo_file_path))
        demographics = [demographic for demographic in demographic_reader]
        
        for subject in subject_dict.keys():
            print subject
            sub,date_string = subject.split('+')
            visit = date_string[:8]
            sub_number = sub.replace('-','')
            try:
                new_sub_number = sub_number[:re.search('[A-Z]|[a-z]',somestr).start()]
            except:
                new_sub_number = sub_number
            
            csv_row_index = [i for i,row in enumerate(demographics) if new_sub_number in row['Mergeid']]
            
            try:
                subject_dict[subject]['demo0'] = demographics[csv_row_index[0]]['Y1_Adhd_Status']
                #subject_dict[subject]['demo1'] = demographics[csv_row_index[0]]['W1_C_Dteamdcsn_Overallsb']
                subject_dict[subject]['demo1'] = demographics[csv_row_index[0]]['Y1_C_Dteamdcsn_Overallsb']
                subject_dict[subject]['demo2'] = demographics[csv_row_index[0]]['Dob']
                subject_dict[subject]['demo3'] = demographics[csv_row_index[0]]['Sex']
            except:
                print 'No demo information for ' + new_sub_number
    
    csvfile = open(output_csv,'w')
    csvfile.write('Subject,Visit,Y1_Adhd_Status,Y1_C_Dteamdcsn_Overallsb,Dob,Sex,FD 0.2 % Frames,FD 0.2 Minutes,FD 0.2 mean FD,FD 0.3 % Frames,FD 0.3 Minutes,FD 0.3 mean FD,Power (FD 0.2/DVAR 20) % Frames,Power (FD 0.2/DVAR 20) Minutes,Power (FD 0.2/DVAR 20) mean FD,Power (FD 0.2/DVAR 20) mean DVAR,Power (FD 0.3/DVAR 20) % Frames,Power (FD 0.3/DVAR 20) Minutes,Power (FD 0.3/DVAR 20) mean FD,Power (FD 0.3/DVAR 20) mean DVAR\n')
#variable changed to Y1_Dteamdcsn_Overallsb    
    for subject in subject_dict.keys():
        sub,date_string = subject.split('+')
        visit = date_string[:8]
        
        try:
            demo0 = subject_dict[subject]['demo0']
            demo1 = subject_dict[subject]['demo1']
            demo2 = subject_dict[subject]['demo2']
            demo3 = subject_dict[subject]['demo3']
        except:
            demo0 = ''
            demo1 = ''
            demo2 = ''
            demo3 = ''
        
        if subject_dict[subject]['prepost_exist']:
            prepost_fdpt2_percent_frames = float( 100.0 * subject_dict[subject]['prepost'][20]['remaining_frame_count'] / subject_dict[subject]['prepost'][20]['total_frame_count'] )
            prepost_fdpt3_percent_frames = float( 100.0 * subject_dict[subject]['prepost'][30]['remaining_frame_count'] / subject_dict[subject]['prepost'][30]['total_frame_count'] )
            prepost_fdpt2_minutes = float( subject_dict[subject]['prepost'][20]['remaining_seconds'] / 60.0 )
            prepost_fdpt3_minutes = float( subject_dict[subject]['prepost'][30]['remaining_seconds'] / 60.0 )
            try:
                prepost_fdpt2_fdmean = float( subject_dict[subject]['prepost'][20]['remaining_frame_mean_FD'] )
                prepost_fdpt3_fdmean = float( subject_dict[subject]['prepost'][30]['remaining_frame_mean_FD'] )
            except:
                print 'No "remaining_frame_mean_FD" numbers available.'
                prepost_fdpt2_fdmean = ''
                prepost_fdpt3_fdmean = ''
        else:
            prepost_fdpt2_percent_frames = ''
            prepost_fdpt3_percent_frames = ''
            prepost_fdpt2_minutes = ''
            prepost_fdpt3_minutes = ''
            prepost_fdpt2_fdmean = ''
            prepost_fdpt3_fdmean = ''
        
        if subject_dict[subject]['power2014_exist']:
            power2014_fdpt2_dvar20_percent_frames = float( 100.0 * subject_dict[subject]['power2014'][20][20]['remaining_frame_count'] / subject_dict[subject]['power2014'][20][20]['total_frame_count'] )
            power2014_fdpt3_dvar20_percent_frames = float( 100.0 * subject_dict[subject]['power2014'][30][20]['remaining_frame_count'] / subject_dict[subject]['power2014'][30][20]['total_frame_count'] )
            power2014_fdpt2_dvar20_minutes = float( subject_dict[subject]['power2014'][20][20]['remaining_seconds'] / 60.0 )
            power2014_fdpt3_dvar20_minutes = float( subject_dict[subject]['power2014'][30][20]['remaining_seconds'] / 60.0 )
            try:
                power2014_fdpt2_dvar20_fdmean = float( subject_dict[subject]['power2014'][20][20]['remaining_frame_mean_FD'] )
                power2014_fdpt3_dvar20_fdmean = float( subject_dict[subject]['power2014'][30][20]['remaining_frame_mean_FD'] )
                power2014_fdpt2_dvar20_dvarmean = float( subject_dict[subject]['power2014'][20][20]['remaining_frame_mean_DVAR'] )
                power2014_fdpt3_dvar20_dvarmean = float( subject_dict[subject]['power2014'][30][20]['remaining_frame_mean_DVAR'] )
            except:
                print 'No "remaining_frame_mean_FD" or "remaining_frame_mean_DVAR" numbers available.'
                power2014_fdpt2_dvar20_fdmean = ''
                power2014_fdpt3_dvar20_fdmean = ''
                power2014_fdpt2_dvar20_dvarmean = ''
                power2014_fdpt3_dvar20_dvarmean = ''
        else:
            power2014_fdpt2_dvar20_percent_frames = ''
            power2014_fdpt3_dvar20_percent_frames = ''
            power2014_fdpt2_dvar20_minutes = ''
            power2014_fdpt3_dvar20_minutes = ''
            power2014_fdpt2_dvar20_fdmean = ''
            power2014_fdpt3_dvar20_fdmean = ''
            power2014_fdpt2_dvar20_dvarmean = ''
            power2014_fdpt3_dvar20_dvarmean = ''
    
        csvfile.write(','.join([ sub , 
                                 visit , 
                                 demo0 , 
                                 demo1 , 
                                 demo2 , 
                                 demo3 , 
                                 str(prepost_fdpt2_percent_frames) , 
                                 str(prepost_fdpt2_minutes) , 
                                 str(prepost_fdpt2_fdmean) , 
                                 str(prepost_fdpt3_percent_frames) , 
                                 str(prepost_fdpt3_minutes) , 
                                 str(prepost_fdpt3_fdmean) , 
                                 str(power2014_fdpt2_dvar20_percent_frames) , 
                                 str(power2014_fdpt2_dvar20_minutes) , 
                                 str(power2014_fdpt2_dvar20_fdmean) , 
                                 str(power2014_fdpt2_dvar20_dvarmean) , 
                                 str(power2014_fdpt3_dvar20_percent_frames) , 
                                 str(power2014_fdpt3_dvar20_minutes) , 
                                 str(power2014_fdpt3_dvar20_fdmean) , 
                                 str(power2014_fdpt3_dvar20_dvarmean)
                                 ]) + '\n')
        
    csvfile.close()

# Exit when complete
if __name__ == "__main__":
    sys.exit(main())
