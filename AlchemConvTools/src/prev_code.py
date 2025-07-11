import pandas as pd
# import numpy as np
from glob import glob
import os
import shutil

def gen_toget_csvs_dict(csvs_base_path, ifall=True, ratio=0.5, forward=True, csv_names=None):
    toget_csvs_dict = {}
    def get_range(total_num, ratio, forward):
        if forward:
            start = 0
            end = int(total_num * ratio)
        else:
            start = int(total_num * (1 - ratio))
            end = total_num
        return start, end
    if ifall:
        for file in os.listdir(csvs_base_path):
            if file.startswith('state') and file.endswith('.csv'):
                one_csv = One_csv_data(os.path.join(csvs_base_path, file))
                num_rows = len(one_csv.df)
                toget_csvs_dict[file] = get_range(num_rows, ratio, forward)
    elif ifall == False and csv_names is not None:
        for file in csv_names:
            if file.startswith('state') and file.endswith('.csv'):
                one_csv = One_csv_data(os.path.join(csvs_base_path, file))
                num_rows = len(one_csv.df)
                toget_csvs_dict[file] = get_range(num_rows, ratio, forward)
    return toget_csvs_dict

def recover_meta_data(csvs_base_path):
    if os.path.exists(os.path.join(csvs_base_path, 'metadata')):
        for file in os.listdir(os.path.join(csvs_base_path, 'metadata')):
            # print(file)
            if file.startswith('state') and file.endswith('.csv'):
                src_file_path = os.path.join(csvs_base_path, 'metadata', file)
                dst_file_path = os.path.join(csvs_base_path, file)
                # print('Copy {src_file_path} to {dst_file_path}')
                shutil.copy(src_file_path, dst_file_path)
    else:
        os.makedirs(os.path.join(csvs_base_path, 'metadata'))
        for file in os.listdir(csvs_base_path):
            # print(file)
            if file.startswith('state') and file.endswith('.csv'):
                src_file_path = os.path.join(csvs_base_path, file)
                dst_file_path = os.path.join(csvs_base_path, 'metadata', file)
                # print('Copy {src_file_path} to {dst_file_path}')
                shutil.copy(src_file_path, dst_file_path)  

class One_mol_data():
    def __init__(self, root_path, mol_name, toget_dict, ratio=0.5, forward=True):
        '''
        toget_dict: dict, like: {'complex': ['restraints', 'electrostatics', 'sterics_group1'], 'ligand': ['electrostatics', 'sterics_group1'] }
        '''
        # print('init____-')
        self.root_path = root_path # /nfs/export3_25T/Bygroup_FEP_data/BRD4_bygroup/1st_batch/
        self.mol_name = mol_name # 3mxf
        self.toget_dict = toget_dict
        for side, themo_ends in self.toget_dict.items():
            for themo_end in themo_ends:
                # print(f'aaaaaa : {themo_end}')
                csvs_base_path = os.path.join(self.root_path, self.mol_name, 'openmm_run', side, themo_end, 'sample_csv_data', )
                recover_meta_data(csvs_base_path) # generate metadata
                toget_csvs_dict = gen_toget_csvs_dict(csvs_base_path, ifall=True, ratio=ratio, forward=forward)
                # print('aaaaaaaaaaaaaaaa')
                print(toget_csvs_dict)
                onethemoend_obj = One_themo_end_data(csvs_base_path, toget_csvs_dict)
                onethemoend_obj.regen_toget_csvs()

class One_themo_end_data():
    def __init__(self, csvs_base_path, toget_csvs_dict=False):
        self.current_path = os.getcwd()
        self.csvs_base_path = csvs_base_path
        if toget_csvs_dict:
            self.toget_csvs_dict = toget_csvs_dict


    def regen_toget_csvs(self, ):
        for file, sample_range in self.toget_csvs_dict.items():
            os.chdir(self.csvs_base_path)
            onecsv_obj = One_csv_data(file, sample_range)
            onecsv_obj.generate_csv()
            os.chdir(self.current_path)

class One_csv_data():
    def __init__(self, csv_path, sample_range=None):
        # print('-------')
        self.csv_path = csv_path
        self.df = pd.read_csv(self.csv_path, delimiter='|', header=[0,])
        # print(self.df)
        # print(self.csv_path)
        if sample_range is not None:
            self.df = self.df.iloc[sample_range[0]:sample_range[1], :]
        
        # print(self.df)
        self.df = self.generate_time_column(self.df)

    def generate_csv(self):
        self.df.to_csv(self.csv_path, index=False, sep='|')

    def generate_time_column(self,df):
        num_rows = len(df)  # 获取 DataFrame 的行数
        start = 0.2  # 起始值
        step = 0.2  # 步长
        new_time_column = [start + i * step for i in range(num_rows)]
        new_time_column = [round(i, 2) for i in new_time_column]

        # print(new_time_column)
        df['times(ps)'] = new_time_column  # 替换 'time(ps)' 列
        # print(df)
        return df


def main():
    root_path = '/nfs/export3_25T/Bygroup_FEP_data/CDK2_bygroup/3rd_batch'
    mol_name = '3ddq'
    # toget_dict = {'complex': ['electrostatics', 'sterics_group1', 'sterics_group4', 'sterics_group8', 'sterics_group9','sterics_group10', 'sterics_group11', 'sterics_group13','sterics_group14', 'sterics_group15', 'sterics_group20', 'sterics_group23', 'sterics_group25', ], 'ligand': ['electrostatics', 'sterics_group5', 'sterics_group7', 'sterics_group10', 'sterics_group12', 'sterics_group25', ]}
    toget_dict = {'complex': ['electrostatics','sterics_group3', 'sterics_group5', 'sterics_group7', 'sterics_group9', 'sterics_group10', 'sterics_group16', 'sterics_group17', 'sterics_group18', 'sterics_group20', 'sterics_group21', 'sterics_group22', 'sterics_group23'],  } #forward
    toget_dict = {'complex': ['restraints', 'sterics_group4', 'sterics_group6', 'sterics_group8', 'sterics_group13', 'sterics_group14', ], } #backward 
    # toget_dict = {'ligand': ['electrostatics', 'sterics_group1', ]}
    # toget_dict = {'ligand': ['electrostatics', 'sterics_group1', 'sterics_group5', 'sterics_group7', 'sterics_group8', 'sterics_group10', 'sterics_group12', 'sterics_group13', 'sterics_group20', 'sterics_group25'], }
    # toget_dict = {'complex': ['electrostatics',], }
    # toget_dict = {'complex': ['electrostatics', 'sterics_group1', 'sterics_group2', 'sterics_group3', 'sterics_group4', 'sterics_group5', 'sterics_group6', 'sterics_group7', 'sterics_group8',  'sterics_group10', 'sterics_group11', 'sterics_group12', 'sterics_group13', 'sterics_group14', 'sterics_group15', 'sterics_group16', 'sterics_group17', 'sterics_group18', 'sterics_group19', 'sterics_group20', 'sterics_group21'], }
    ratio = 0.5
    forward = True
    one_mol_obj = One_mol_data(root_path, mol_name, toget_dict, ratio, forward)

if __name__ == '__main__':
    main()
