# from email.policy import default
from optparse import OptionParser
from src.parsing.input_file_parser import InputParser as InputParser
from src.all_fe_aly_cls.overall_fe_aly_cls import singleside_ANA_ABFE_OVERALL as singleside_ANA_ABFE_OVERALL
import os




class optParser():  
    def __init__(self, fakeArgs):
        parser = OptionParser()
        parser.add_option('-i', '--input', dest='input', help="The file name of input, which recording the analysis settings. Default: 'input.txt'", default='input.txt')
        if fakeArgs:
            self.option, self.args = parser.parse_args(fakeArgs)
        else:
            self.option, self.args = parser.parse_args()

opts = optParser('')
#####################################################################################################
# keep for jupyter notebook test
# fakeArgs = '-i input.txt'
# opts = optParser(fakeArgs.strip().split())
#####################################################################################################
input_parser = InputParser(opts.option.input)
basic_settings = input_parser.get_Basic_settings()

# Basic settings
simulation_pack = str(basic_settings['simulation_software'])
file_str = str(basic_settings['file_directory']) # For openmm, it is the directory that stores the state_s*.csv or state_g*_s*.csv files.
file_prefix = str(basic_settings['file_prefix'])
file_suffix = str(basic_settings['file_suffix'])
fraction = float(basic_settings['fraction'])
output_csv_filename = str(basic_settings['output_csv_filename'])
energy_unit = str(basic_settings['energy_unit'])
std_mode = str(basic_settings['std_mode'])
cal_wins = str(basic_settings['calculation_windows'])
ifplot_du = basic_settings['plot_du']
ifsubsample = basic_settings['subsample']

fe_cls = singleside_ANA_ABFE_OVERALL(file_str, simulation_pack, file_prefix, file_suffix, [0,1,2,3], '|', ifsubsample)

if cal_wins == 'all':
    lambda_info_list = list(fe_cls.com_u_nk_pd.index.names)
    lambda_info_list.pop(0)
    simulation_win_lst = [i for i,j in fe_cls.com_u_nk_pd.groupby(lambda_info_list, sort=False)]
    force_field_cal_win_lst = list(fe_cls.com_u_nk_pd.columns)
    cal_win_lst = simulation_win_lst
    cal_win_lst = [cal_win_lst, ]
else:
    with open(cal_wins, 'r') as f:
        f_content = f.readlines()
        cal_win_lst = [eval(line.strip()) for line in f_content]
        import sys 
        if not isinstance(cal_win_lst, list):
            print('Error: The content of the win_lst file is not a list!')
            sys.exit()
print(cal_win_lst[0])

def jobs(fe_cal_cls, store_path, cal_win, energy_unit, std_mode, fraction, ifplot_du):
    fe_cal_cls.cal_fe(store_path, output_csv_filename, cal_win, energy_unit, 'BAR', std_mode, ifplot_du, fraction)

if len(cal_win_lst) == 1:
    jobs(fe_cls, 'fe_cal_out', cal_win_lst[0], energy_unit, std_mode, fraction, ifplot_du)
else:
    for idx_ in range(0, len(cal_win_lst)):
        win_lst = cal_win_lst[idx_]
        print(win_lst)
        jobs(fe_cls, f'fe_cal_out_win_lst_{idx_}', win_lst, energy_unit, std_mode, fraction, ifplot_du)
