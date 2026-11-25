#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
from dataclasses import dataclass
from typing import Tuple, List, Dict
import parmed as pmd
import logging
import numpy as np
import contextlib
from concurrent.futures import ProcessPoolExecutor
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class BoxParameters:
    """存储盒子参数的数据类"""
    box_size: str  # 盒子三个维度的尺寸，空格分隔的字符串
    side_x: int    # x方向的边长
    wat_num: int   # 水分子数量

@contextlib.contextmanager
def working_directory(path):
    """一个用于临时更改工作目录的上下文管理器
    
    Args:
        path: 要切换到的目录路径
    """
    prev_cwd = os.getcwd()
    os.makedirs(path, exist_ok=True)  # 如果目录不存在则创建
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)

class APBSConfig:
    """APBS配置类"""
    def __init__(self, num_parallel: int = 5):
        self.num_parallel = num_parallel
        self.solvent_dielectric = {
            'het': 97.0,
            'hom': 1.0
        }
        # 定义不同计算场景的配置
        self.calc_configs = {
            'complex_protein': {
                'pqr_file': 'prot_only.pqr',
                'output_prefix': 'protein_RIP'
            },
            'complex_ligand': {
                'pqr_file': 'lig_in_prot.pqr',
                'output_prefix': 'ligand_RIP'
            },
            'ligand_only': {
                'pqr_file': 'lig_in_prot.pqr',
                'output_prefix': 'ligand_RIP'
            }
        }

@dataclass
class PhysicalConstants:
    """物理常数数据类"""
    # 介电常数
    eps_het: float = 97.0    # TIP3P水的相对介电常数
    eps_hom: float = 1.0     # 真空的相对介电常数
    
    # 温度相关
    temperature: float = 300  # 温度 (K)
    kB: float = 0.0083144621  # 玻尔兹曼常数，单位kJ/(mol K)
    
    # 静电相关
    coulomb_factor: float = 138.93545585  # (4*pi*eps_0)**-1，单位(kJ nm)/(e**2 mol)
    xi_CB: float = 2.38008   # Born自相关项，见论文方程(21)
    ELS: float = -2.837297   # 立方晶格和（维格纳）积分常数
    
    # 溶剂相关
    ys: float = 0.00764      # TIP3P水的四极矩迹（相对于其单个范德华相互作用位点），单位e nm²
    
    @property
    def kBT(self) -> float:
        """热能 (kJ/mol)"""
        return self.kB * self.temperature
    
    @property
    def e0(self) -> float:
        """真空介电常数"""
        return 1/(self.coulomb_factor*4*np.pi)

class ChargeCorrection:
    def __init__(self, complex_crd, complex_top, ligand_crd, ligand_top, ligand_residue_name='MOL', chg_mol: int = -1, chg_com: int = 3,
                 apbs_config: APBSConfig = None, phys_const: PhysicalConstants = None):
        """初始化电荷校正计算
        
        Args:
            complex_crd: 复合物的坐标文件
            complex_top: 复合物的拓扑文件
            ligand_crd: 配体的坐标文件
            ligand_top: 配体的拓扑文件
            ligand_residue_name: 配体的残基名称
            chg_mol: 分子的电荷
            chg_com: 复合物的电荷
            apbs_config: APBS配置对象
            phys_const: 物理常数对象
        """
        self.complex_crd = complex_crd
        self.complex_top = complex_top
        self.ligand_crd = ligand_crd
        self.ligand_top = ligand_top
        self.ligand_residue_name = ligand_residue_name
        self.chg_mol = chg_mol
        self.chg_com = chg_com
        self.chg_pro = chg_com - chg_mol
        logger.info("Starting charge correction workflow")            
        # 1. 准备PDB文件
        logger.info("Step 1: Preparing PDB files")
        self.prepare_pdb_files()

        # Get box parameters of complex
        self.complex_box_parms = self.get_cubic_box_parms("complex.pdb")
        # Get box parameters of ligand
        self.ligand_box_parms = self.get_cubic_box_parms("ligand.pdb")
        # 2. 准备PQR文件
        logger.info("Step 2: Preparing PQR files")
        self.prepare_pqr_files()
        self.apbs_config = apbs_config or APBSConfig()
        self.phys_const = phys_const or PhysicalConstants()

    def check_cubic_box(self, structure: pmd.Structure) -> bool:
        """
        Check if the box is cubic
        
        Args:
            structure: ParmEd Structure object
            
        Returns:
            bool: True if box is cubic, False otherwise
            
        Raises:
            ValueError: With detailed information if box is not cubic
            
        Note:
            Box is considered cubic if:
            - All angles are within 0.1° of 90°
            - All sides are within 1.0 Å of each other
        """
        # Check if box information exists
        if structure.box is None or not isinstance(structure.box, np.ndarray):
            raise ValueError("No valid box information found in structure")
            
        # Get box dimensions and angles
        box = structure.box
        if len(box) != 6:
            raise ValueError("Invalid box format: expected 6 values (a,b,c,alpha,beta,gamma)")
        
        a, b, c = box[0:3]
        alpha, beta, gamma = box[3:6]
        
        # Check if angles are 90° (within 0.1° tolerance)
        angles = [alpha, beta, gamma]
        if not all(abs(angle - 90.0) < 0.1 for angle in angles):
            raise ValueError(
                f"Box angles must be 90°. Current angles: "
                f"alpha={alpha:.2f}°, beta={beta:.2f}°, gamma={gamma:.2f}°"
            )
        
        # Check if sides are equal (within 1.0 Å tolerance)
        TOLERANCE = 1.0  # Å
        if not (abs(a - b) < TOLERANCE and abs(b - c) < TOLERANCE):
            raise ValueError(
                f"Box must be approximately cubic (sides within {TOLERANCE}Å). Current sides: "
                f"a={a:.3f}, b={b:.3f}, c={c:.3f} Å"
            )
        
        return True

    def prepare_pdb_files(self):
        """基于parmed 准备PDB文件"""
        # 从复合物的拓扑文件中提取配体的残基名称
        logger.info("Preparing PDB files of complex")
        complex_struct = pmd.load_file(self.complex_top, self.complex_crd)
        self.check_cubic_box(complex_struct)
        complex_struct.save("complex.pdb", overwrite=True)

        logger.info("Preparing PDB files of ligand")
        ligand_struct = pmd.load_file(self.ligand_top, self.ligand_crd)
        self.check_cubic_box(ligand_struct)
        ligand_struct.save("ligand.pdb", overwrite=True)

    def prepare_pqr_files(self):
        """从AMBER输出文件准备PQR文件"""
        def call_ambpdb(prmtop_file, crd_file, pqr_file):
            """调用ambpdb生成PQR文件
            
            Args:
                prmtop_file: 拓扑文件路径
                crd_file: 坐标文件路径
                pqr_file: 输出的PQR文件路径
            """
            try:
                with open(crd_file, 'r') as crd, open(pqr_file, 'w') as pqr:
                    process = subprocess.Popen(
                        ["ambpdb", "-p", str(prmtop_file), "-pqr", "-aatm", "-bres", "-ext"],
                        stdin=crd,
                        stdout=pqr,
                        stderr=subprocess.PIPE
                    )
                    _, stderr = process.communicate()
                    
                    if process.returncode != 0:
                        raise RuntimeError(f"ambpdb failed with error: {stderr.decode()}")
                        
                logger.info(f"Successfully generated {pqr_file}")
                
            except Exception as e:
                logger.error(f"Error running ambpdb: {str(e)}")
                raise
        
        try:
            # 调用ambpdb生成com.pqr
            logger.info("Generating PQR files using ambpdb")
            call_ambpdb(self.complex_top, self.complex_crd, "com.pqr")
            
            # 提取蛋白质和配体原子并生成所有需要的pqr文件
            logger.info("Extracting PQR files from com.pqr")
            with open("com.pqr", 'r') as f:
                lines = f.readlines()
                
            # 提取ATOM行到p.pqr
            with open("p.pqr", 'w') as f:
                for line in lines:
                    if line.startswith('ATOM'):
                        f.write(line)
                        
            # 提取MOL行到l.pqr，并将HETATM改为ATOM
            with open("l.pqr", 'w') as f:
                for line in lines:
                    if self.ligand_residue_name in line:
                        modified_line = line.replace('HETATM', 'ATOM  ')
                        f.write(modified_line)
            
            # 创建零电荷版本
            for infile, outfile in [('l.pqr', 'l0.pqr'), ('p.pqr', 'p0.pqr')]:
                logger.debug(f"Processing {infile} -> {outfile}")
                with open(infile) as f, open(outfile, 'w') as out:
                    for i, line in enumerate(f, 1):
                        try:
                            # 跳过非坐标行
                            if line.startswith(('TER', 'END', 'REMARK')):
                                out.write(line)
                                continue
                            
                            # 解析PQR格式的固定宽度字段
                            # 保持原始格式，只修改电荷为0
                            fields = [
                                line[:30],           # 保持atom/residue信息不变
                                line[30:38],         # x coordinate
                                line[38:46],         # y coordinate
                                line[46:54],         # z coordinate
                                "0.0000",            # 将电荷设为0
                                line[62:70],         # 保持radius不变
                                line[70:]            # 保持其余部分不变
                            ]
                            
                            # 使用固定宽度格式写入新行
                            new_line = (f"{fields[0]}"     # atom/residue info
                                    f"{fields[1]}"       # x
                                    f"{fields[2]}"       # y
                                    f"{fields[3]}"       # z
                                    f"{float(fields[4]):8.4f}"  # charge (0.0000)
                                    f"{fields[5]}"       # radius
                                    f"{fields[6]}")      # rest of line
                            
                            out.write(new_line)
                            
                        except (ValueError, IndexError) as e:
                            if line.strip():  # 只报告非空行的错误
                                logger.error(f"Error processing line {i} in {infile}:")
                                logger.error(f"Line content: {line.rstrip()}")
                                logger.error(f"Error: {str(e)}")
                                raise ValueError(f"Failed to process {infile} at line {i}: {str(e)}")
                                
                logger.info(f"Successfully processed {infile} -> {outfile}")
                
            # 处理complex目录
            with working_directory('complex'):
                logger.info("Processing complex directory")
                # 复制所需文件
                for file in ['l.pqr', 'l0.pqr', 'p0.pqr', 'p.pqr']:
                    if os.path.exists(f"../{file}"):
                        with open(f"../{file}", 'r') as src, open(file, 'w') as dst:
                            dst.write(src.read())
                    
                # 创建所需的组合文件
                with open("lig_only.pqr", 'w') as f:
                    with open("l.pqr", 'r') as src:
                        f.write(src.read())
                        
                with open("lig_in_prot.pqr", 'w') as f:
                    with open("p0.pqr", 'r') as src1:
                        f.write(src1.read())
                    with open("l.pqr", 'r') as src2:
                        f.write(src2.read())

                with open("prot_only.pqr", 'w') as f:
                    with open("p.pqr", 'r') as src1:
                        f.write(src1.read())
                    with open("l0.pqr", 'r') as src2:
                        f.write(src2.read())

            # 处理ligand目录
            with working_directory('ligand'):
                logger.info("Processing ligand directory")
                # 复制所需文件
                for file in ['l.pqr', 'l0.pqr']:
                    if os.path.exists(f"../{file}"):
                        with open(f"../{file}", 'r') as src, open(file, 'w') as dst:
                            dst.write(src.read())
                    
                # 创建所需的文件
                for src, dst in [
                    ('l.pqr', 'lig_only.pqr'),
                    ('l.pqr', 'lig_in_prot.pqr'),
                    ('l0.pqr', 'prot_only.pqr')
                ]:
                    with open(src, 'r') as src_f, open(dst, 'w') as dst_f:
                        dst_f.write(src_f.read())

            logger.info("PQR file processing completed")
            
        except Exception as e:
            logger.error(f"Error during PQR file processing: {str(e)}")
            raise


    def get_cubic_box_parms(self, pdb_file: str) -> BoxParameters:
        """计算立方体盒子的参数
        
        Args:
            pdb_file: PDB文件路径
            
        Returns:
            BoxParameters对象，包含盒子尺寸、x边长和水分子数量
        """
        box_sizes = []
        
        # 读取PDB文件
        with open(pdb_file, 'r') as f:
            lines = f.readlines()
            
        # 处理每个维度（x, y, z对应第6-8列）
        for i in range(6, 9):
            coords = []
            for line in lines:
                if line.startswith(('ATOM', 'HETATM')):
                    try:
                        coord = float(line.split()[i-1])
                        coords.append(coord)
                    except (IndexError, ValueError):
                        continue
            
            if not coords:
                raise ValueError(f"No valid coordinates found for dimension {i-5}")
                
            max_x = max(coords)
            min_x = min(coords)
            side = max_x - min_x
            box_sizes.append(side)
            
        # 计算水分子数量
        wat_num = sum(1 for line in lines if 'WAT' in line) // 3
        
        # 构建返回值
        box_size = ' '.join(map(str, box_sizes))
        side_x = int(box_sizes[0])
        
        return BoxParameters(box_size=box_size, side_x=side_x, wat_num=wat_num)

    @staticmethod
    def run_single_apbs(args: Tuple[str, str]) -> None:
        """
        Run APBS calculation with given arguments

        Args:
            args: Tuple of (working_dir, apbs_in_file)
        """
        working_dir, apbs_in_file = args
        with working_directory(working_dir):
            logger.info(f"Running APBS in {working_dir} with {apbs_in_file}")
            subprocess.run(["apbs", apbs_in_file], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    def prepare_apbs_input(self, template_file: str, output_file: str, 
                          box_size: str, dielectric: float, calc_type: str,
                          pqr_file: str, output_prefix: str) -> None:
        """准备APBS输入文件
        
        Args:
            template_file: 模板文件路径
            output_file: 输出文件路径
            box_size: 盒子尺寸
            dielectric: 介电常数
            calc_type: 计算类型 ('het' 或 'hom')
            pqr_file: 要计算的PQR文件
            output_prefix: 输出文件前缀
        """
        with open(template_file, 'r') as f:
            content = f.read()
        
        # 替换参数
        content = content.replace('_BOX_SIZE_', box_size)
        content = content.replace('_SOLVENT_DIELECTRIC_', str(dielectric))
        content = content.replace('_PQR_FILE_', pqr_file)
        content = content.replace('_OUTPUT_FILE_', f"{output_prefix}_{calc_type}")
        
        with open(output_file, 'w') as f:
            f.write(content)

    def run_apbs_calculations(self, parallel: bool = True) -> None:
        """运行APBS计算
        
        Args:
            parallel: 是否使用并行计算
        """
        try:
            # 获取脚本所在目录的路径
            script_dir = os.path.dirname(os.path.abspath(__file__))
            template_path = os.path.join(script_dir, 'apbs_tmp.in')
            
            if not os.path.exists(template_path):
                raise FileNotFoundError(f"APBS template file not found at {template_path}")
            
            # 准备APBS计算任务列表
            apbs_tasks = []
            
            # 处理complex目录
            with working_directory('complex'):
                # 复制模板文件
                shutil.copy(template_path, 'apbs_tmp.in')
                
                # 生成complex中protein和ligand的计算文件
                for calc_scenario in ['complex_protein', 'complex_ligand']:
                    config = self.apbs_config.calc_configs[calc_scenario]
                    for calc_type in ['het', 'hom']:
                        # complex中的ligand只需要het计算
                        if calc_scenario == 'complex_protein' and calc_type == 'hom':
                            continue
                            
                        output_file = f'apbs_{calc_scenario}_{calc_type}.in'
                        dielectric = self.apbs_config.solvent_dielectric[calc_type]
                        self.prepare_apbs_input(
                            'apbs_tmp.in', 
                            output_file,
                            self.complex_box_parms.box_size,
                            dielectric,
                            calc_type,
                            config['pqr_file'],
                            config['output_prefix']
                        )
                        apbs_tasks.append(('complex', output_file))
                        
            # 处理ligand目录
            with working_directory('ligand'):
                # 复制模板文件
                shutil.copy(template_path, 'apbs_tmp.in')
                
                # 生成ligand的计算文件
                config = self.apbs_config.calc_configs['ligand_only']
                for calc_type in ['het', 'hom']:
                    output_file = f'apbs_{calc_type}.in'
                    dielectric = self.apbs_config.solvent_dielectric[calc_type]
                    self.prepare_apbs_input(
                        'apbs_tmp.in', 
                        output_file,
                        self.ligand_box_parms.box_size,
                        dielectric,
                        calc_type,
                        config['pqr_file'],
                        config['output_prefix']
                    )
                    apbs_tasks.append(('ligand', output_file))

            logger.info(f"Prepared {len(apbs_tasks)} APBS calculations")
            
            # 运行APBS计算
            if parallel:
                # 并行执行
                with ProcessPoolExecutor(max_workers=self.apbs_config.num_parallel) as executor:
                    list(executor.map(self.run_single_apbs, apbs_tasks))
            else:
                # 串行执行
                for task in apbs_tasks:
                    self.run_single_apbs(task)
                    
            logger.info("All APBS calculations completed")
            
        except Exception as e:
            logger.error(f"Error during APBS calculations: {str(e)}")
            raise

    def run_calculations(self):
        """运行所有计算（APBS和RIP）的接口方法"""
        # 这里将添加调用APBS和其他Python脚本的接口
        pass

    def calculate_RIP(self, dx_file: str, net_charge: float, eps: float, calc_type: str = '') -> Dict[str, float]:
        """计算残余积分电势 (RIP)
        
        Args:
            dx_file: APBS输出的DX文件路径
            net_charge: 系统的净电荷
            eps: 相对介电常数 (het: 97.0 或 hom: 1.0)
            calc_type: 计算类型标识（如'protein_het', 'ligand_het'等）
            
        Returns:
            Dict[str, float]: 包含计算结果的字典
        """
        logger.info(f"Calculating RIP for {calc_type} (dx_file: {dx_file})")
        
        # 读取DX文件
        with open(dx_file, 'r') as f:
            lines = f.readlines()
            
        # 提取网格间距
        dx = float(lines[6].split()[1])  # x方向间距
        dy = float(lines[7].split()[2])  # y方向间距
        dz = float(lines[8].split()[3])  # z方向间距
        
        # 计算体积元素 (Å³ -> nm³)
        dV = (dx * dy * dz) / 1000.0
        
        # 提取电势数据
        data_lines = filter(lambda x: x[0] in '0123456789-.', lines)
        
        # 计算总电势和网格点数
        tot_pot = 0.0
        numpts = 0
        for line in data_lines:
            values = [float(x) for x in line.split()]
            tot_pot += sum(values)
            numpts += len(values)
            
        # 计算总体积
        V = dV * numpts
        
        # 计算盒子积分 B_X (kJ nm³/(mol e))
        B_X = tot_pot * self.phys_const.kBT * V / numpts
        
        # 计算参考电势 B_QX
        B_QX = (net_charge * self.phys_const.xi_CB * self.phys_const.coulomb_factor 
                / eps * (V ** (2.0 / 3.0)))
        
        # 计算RIP
        RIP = B_X - B_QX
        
        # 创建结果字典
        results = {
            'gridpoints': numpts,
            'volume': V,
            'B_X': B_X,
            'B_QX': B_QX,
            'RIP_kj': RIP,
            'RIP_kcal': RIP / 4.184
        }
        
        # 保存计算结果到输出文件
        out_file = f"{calc_type}_RIP.out"
        with open(out_file, 'w') as f:
            f.write(f"RIP Calculation Results for {calc_type}\n")
            f.write("-" * 50 + "\n")
            f.write(f"Input DX file: {dx_file}\n")
            f.write(f"Net charge: {net_charge}\n")
            f.write(f"Dielectric constant: {eps}\n")
            f.write("-" * 50 + "\n")
            f.write(f"Grid points: {numpts}\n")
            f.write(f"Volume: {V:.4f} nm³\n")
            f.write(f"B_X: {B_X:.4f} kJ nm³/(mol e)\n")
            f.write(f"B_QX: {B_QX:.4f} kJ nm³/(mol e)\n")
            f.write(f"RIP: {RIP:.4f} kJ nm³/(mol e)\n")
            f.write(f"RIP: {RIP/4.184:.4f} kcal nm³/(mol e)\n")
        
        logger.info(f"RIP calculation completed for {calc_type}")
        logger.info(f"Results saved to {out_file}")
        
        return results

    def cal_complex_rip(self) -> Dict[str, Dict]:
        """计算复合物体系的RIP值，使用多进程并行计算
        
        Returns:
            Dict[str, Dict]: 包含蛋白质和配体RIP计算结果的字典
        """
        logger.info("Starting complex RIP calculations")
        results = {}
        
        # 定义要计算的任务
        tasks = [
            ('protein_het', 'protein_RIP_het.dx', self.chg_pro, self.phys_const.eps_het),
            ('ligand_het', 'ligand_RIP_het.dx', self.chg_mol, self.phys_const.eps_het),
            ('ligand_hom', 'ligand_RIP_hom.dx', self.chg_mol, self.phys_const.eps_hom)
        ]
        
        # 在complex目录下执行计算
        with working_directory('complex'):
            # 使用进程池并行计算
            with ProcessPoolExecutor(max_workers=3) as executor:
                # 提交所有任务
                future_to_task = {
                    executor.submit(
                        self.calculate_RIP, 
                        dx_file=task[1],
                        net_charge=task[2],
                        eps=task[3],
                        calc_type=task[0]
                    ): task[0] for task in tasks
                }
                
                # 收集结果
                for future in future_to_task:
                    task_name = future_to_task[future]
                    try:
                        results[task_name] = future.result()
                        logger.info(f"Complex {task_name} RIP calculation completed")
                    except Exception as e:
                        logger.error(f"Error in complex {task_name} RIP calculation: {str(e)}")
                        raise
        
        logger.info("Complex RIP calculations completed")
        return results

    def cal_ligand_rip(self) -> Dict[str, Dict]:
        """计算配体体系的RIP值，使用多进程并行计算
        
        Returns:
            Dict[str, Dict]: 包含配体RIP计算结果的字典
        """
        logger.info("Starting ligand RIP calculations")
        results = {}
        
        # 定义要计算的任务
        tasks = [
            ('ligand_het', 'ligand_RIP_het.dx', self.chg_mol, self.phys_const.eps_het),
            ('ligand_hom', 'ligand_RIP_hom.dx', self.chg_mol, self.phys_const.eps_hom)
        ]
        
        # 在ligand目录下执行计算
        with working_directory('ligand'):
            # 使用进程池并行计算
            with ProcessPoolExecutor(max_workers=2) as executor:
                # 提交所有任务
                future_to_task = {
                    executor.submit(
                        self.calculate_RIP, 
                        dx_file=task[1],
                        net_charge=task[2],
                        eps=task[3],
                        calc_type=task[0]
                    ): task[0] for task in tasks
                }
                
                # 收集结果
                for future in future_to_task:
                    task_name = future_to_task[future]
                    try:
                        results[task_name] = future.result()
                        logger.info(f"Ligand {task_name} RIP calculation completed")
                    except Exception as e:
                        logger.error(f"Error in ligand {task_name} RIP calculation: {str(e)}")
                        raise
        
        logger.info("Ligand RIP calculations completed")
        return results

    def calculate_binding_correction(self) -> Dict[str, Dict[str, float]]:
        """计算结合自由能的完整校正
        
        Returns:
            Dict: 包含复合物侧和配体侧校正项的字典
        """
        results = {
            'complex': {},
            'ligand': {},
            'binding': {}
        }
        
        # 计算复合物侧的校正
        self.complex_rips = self.cal_complex_rip()
        results['complex'] = self.one_side_charge_correction(
            self.complex_rips,
            self.complex_box_parms.side_x,
            self.complex_box_parms.wat_num,
            True
        )
        
        # 计算配体侧的校正
        self.ligand_rips = self.cal_ligand_rip()
        results['ligand'] = self.one_side_charge_correction(
            self.ligand_rips,
            self.ligand_box_parms.side_x,
            self.ligand_box_parms.wat_num,
            False
        )
        
        # 计算结合自由能校正（配体侧减去复合物侧）
        for correction_type in [
            'ddG_RIP',
            'ddG_NET',
            'ddG_USV',
            'ddG_EMP',
            'ddG_DSC',
            'ddG_total',
        ]:
            results['binding'][correction_type] = (
                results['ligand'][correction_type] - 
                results['complex'][correction_type]
            )
        
        return results

    def one_side_charge_correction(self, rip_results: Dict[str, Dict], box_length: float, wat_num: int, is_complex: bool = True) -> Dict[str, float]:
        """计算单侧（复合物或配体）的校正项
        
        Args:
            rip_results: RIP计算结果字典，包含'protein_het', 'ligand_het', 'ligand_hom'等键
            box_length: 模拟盒子的边长（单位：Å）
            wat_num: 水分子数量
            is_complex: 是否是复合物侧的计算（默认True）
            
        Returns:
            Dict[str, float]: 包含各校正项的字典（单位：kJ/mol）
        """
        # 提取RIP值
        IP = rip_results.get('protein_het', {}).get('RIP_kj', 0)  # 如果是配体侧，这个值为0
        IL = rip_results['ligand_het']['RIP_kj']
        IL_hom = rip_results['ligand_hom']['RIP_kj']
        
        # 根据计算侧设置电荷
        protein_charge = self.chg_pro if is_complex else 0
        #protein_charge=0
        # 转换盒子长度到nm
        L = box_length / 10.0  # Å -> nm
        
        # 计算配体半径
        # RL = ((IL - IL_hom)/(self.phys_const.coulomb_factor/2 * 4*np.pi/3 * 
        #              (1-1/self.phys_const.eps_het) * abs(self.chg_mol)))**(1/2)
        RL = (abs(IL - IL_hom)/(self.phys_const.coulomb_factor/2 * 4*np.pi/3 * 
                     (1-1/self.phys_const.eps_het) * abs(self.chg_mol)))**(1/2)
        
        
        # 计算各校正项
        q_total = protein_charge + self.chg_mol

        ddg_rip = ((IP + IL) * q_total - IP * protein_charge) / L**3
        delta_q2 = (protein_charge**2) - (q_total**2)
        # Rocklin charge correction split:
        #   ΔG_NET(L) = - ξ_LS/(8π ε0) * ΔQ^2 / L
        #   ΔG_USV(L) = + ξ_LS/(8π ε0) * (1 - 1/ε_S) * ΔQ^2 / L
        # Here, 1/(8π ε0) == coulomb_factor/2 in our unit system.
        prefactor = self.phys_const.ELS * self.phys_const.coulomb_factor / 2

        ddg_net = -prefactor * delta_q2 / L
        ddg_usv = prefactor * (1 - 1 / self.phys_const.eps_het) * delta_q2 / L

        ddg_emp = (
            -self.phys_const.coulomb_factor
            / 2
            * (16 * np.pi**2 / 45)
            * (1 - 1 / self.phys_const.eps_het)
            * delta_q2
            * RL**5
            / L**6
        )
        delta_q = -self.chg_mol
        ddg_dsc = -((self.phys_const.ys * wat_num * delta_q) / (6 * self.phys_const.e0 * L**3))

        corrections = {
            'ddG_RIP': ddg_rip,
            'ddG_NET': ddg_net,
            'ddG_USV': ddg_usv,
            'ddG_EMP': ddg_emp,
            'ddG_DSC': ddg_dsc,
        }
         
        print(-self.phys_const.coulomb_factor/2)
        print(protein_charge)
        print(self.chg_mol)
        print(L**6)
        print( RL**5)

        # 计算总校正能
        corrections['ddG_total'] = (
            corrections['ddG_RIP']
            + corrections['ddG_NET']
            + corrections['ddG_USV']
            + corrections['ddG_EMP']
            + corrections['ddG_DSC']
        )
        
        return corrections

    def save_correction_results(self, results: Dict[str, Dict[str, float]], output_file: str = "charge_corrections.csv") -> None:
        """将校正结果保存为CSV文件
        
        Args:
            results: calculate_binding_correction返回的结果字典
            output_file: 输出CSV文件的路径，默认为'charge_corrections.csv'
        """
        try:
            import pandas as pd
            
            # 创建数据列表
            data = []
            
            # 遍历每个系统（complex, ligand, binding）
            for system, corrections in results.items():
                # 遍历该系统的所有校正项
                for correction_type, value in corrections.items():
                    data.append({
                        'System': system,
                        'Correction_Type': correction_type,
                        'Value_kJ_mol': value,
                        'Value_kcal_mol': value / 4.184
                    })
            
            # 创建DataFrame
            df = pd.DataFrame(data)
            
            # 重新排序行，使binding结果在最后
            system_order = {'complex': 0, 'ligand': 1, 'binding': 2}
            df['System_Order'] = df['System'].map(system_order)
            df = df.sort_values(['System_Order', 'Correction_Type'])
            df = df.drop('System_Order', axis=1)
            
            # 保存到CSV文件
            df.to_csv(output_file, index=False, float_format='%.4f')
            logger.info(f"Results saved to {output_file}")
            
            # 打印结果摘要
            print("\nCharge Correction Results Summary:")
            print("-" * 60)
            for system in ['complex', 'ligand', 'binding']:
                print(f"\n{system.capitalize()} System:")
                system_data = df[df['System'] == system]
                for _, row in system_data.iterrows():
                    print(f"{row['Correction_Type']:12s}: {row['Value_kJ_mol']:8.4f} kJ/mol "
                          f"({row['Value_kcal_mol']:8.4f} kcal/mol)")
            
        except ImportError:
            logger.warning("pandas not installed. Saving results in simple text format.")
            with open(output_file, 'w') as f:
                f.write("System,Correction_Type,Value_kJ_mol,Value_kcal_mol\n")
                for system, corrections in results.items():
                    for correction_type, value in corrections.items():
                        f.write(f"{system},{correction_type},{value:.4f},{value/4.184:.4f}\n")
            logger.info(f"Results saved to {output_file}")

    def save_details_results(self, results: Dict[str, Dict[str, float]], output_file: str = "charge_corrections_details.csv") -> None:
        """将详细的校正结果保存为CSV文件
        
        Args:
            results: calculate_binding_correction返回的结果字典
            output_file: 输出CSV文件的路径，默认为'charge_corrections_details.csv'
        """
        try:
            import pandas as pd
            
            # 创建数据列表
            data = []
            
            # 遍历每个系统（complex和ligand）
            for side in ['complex', 'ligand']:
                # 获取RIP值
                if side == 'complex':
                    rip_results = self.complex_rips
                    box_length = float(self.complex_box_parms.side_x)/10
                    wat_num = self.complex_box_parms.wat_num
                    q_prot = self.chg_pro
                else:
                    rip_results = self.ligand_rips
                    box_length = float(self.ligand_box_parms.side_x)/10
                    wat_num = self.ligand_box_parms.wat_num
                    q_prot = 0
                
                # 提取RIP值
                i_l_het = rip_results['ligand_het']['RIP_kj']
                i_l_hom = rip_results['ligand_hom']['RIP_kj']
                i_p_het = rip_results.get('protein_het', {}).get('RIP_kj', 0)
                
                # 计算配体半径
                r_l = (abs(i_l_het - i_l_hom)/(self.phys_const.coulomb_factor/2 * 4*np.pi/3 * 
                       (1-1/self.phys_const.eps_het) * abs(self.chg_mol)))**(1/2)
                
                # 获取校正项
                corrections = results[side]
                
                # 添加到数据列表
                data.append({
                    'Side': side,
                    'Ns': wat_num,
                    'L(nm)': box_length,
                    'Q_prot': q_prot,
                    'Q_liga': self.chg_mol,
                    'I_L_het(kJ nm3/mol/e)': i_l_het,
                    'I_L_hom(kJ nm3/mol/e)': i_l_hom,
                    'I_P_het(kJ nm3/mol/e)': i_p_het,
                    'R_L(nm)': r_l,
                    'dG_net(kJ/mol)': corrections['ddG_NET'],
                    'dG_usv(kJ/mol)': corrections['ddG_USV'],
                    'dG_rip(kJ/mol)': corrections['ddG_RIP'],
                    'dG_emp(kJ/mol)': corrections['ddG_EMP'],
                    'dG_dsc(kJ/mol)': corrections['ddG_DSC'],
                    'dG_total(kJ/mol)': corrections['ddG_total']
                })
            
            # 创建DataFrame
            df = pd.DataFrame(data)
            
            # 保存到CSV文件
            df.to_csv(output_file, index=False, float_format='%.4f')
            logger.info(f"Detailed results saved to {output_file}")
            
            # 打印结果摘要
            print("\nDetailed Charge Correction Results:")
            print("-" * 100)
            print(df.to_string(float_format=lambda x: '{:.4f}'.format(x)))
            print("-" * 100)
            
        except ImportError:
            logger.warning("pandas not installed. Saving results in simple text format.")
            with open(output_file, 'w') as f:
                # 写入表头
                headers = ['Side', 'Ns', 'L(nm)', 'Q_prot', 'Q_liga', 'I_L_het(kJ nm3/mol/e)', 'I_L_hom(kJ nm3/mol/e)', 
                          'I_P_het(kJ nm3/mol/e)', 'R_L(nm)', 'dG_net(kJ/mol)', 'dG_usv(kJ/mol)', 'dG_rip(kJ/mol)', 'dG_emp(kJ/mol)', 'dG_dsc(kJ/mol)', 'dG_total(kJ/mol)']
                f.write(','.join(headers) + '\n')
                
                # 写入数据
                for row in data:
                    values = [str(row[h]) for h in headers]
                    f.write(','.join(values) + '\n')
            logger.info(f"Detailed results saved to {output_file}")

    def workflow(self) -> Dict[str, Dict[str, float]]:
        """执行完整的电荷校正计算工作流
        
        工作流程：
        1. 准备PDB文件
        2. 准备PQR文件
        3. 运行APBS计算
        4. 计算RIP值
        5. 计算校正项
        6. 保存结果
        
        Returns:
            Dict[str, Dict[str, float]]: 包含所有校正项的计算结果
        """
        try:
            
            # 3. 运行APBS计算
            logger.info("Step 3: Running APBS calculations")
            self.run_apbs_calculations(parallel=True)
            
            # 4&5. 计算RIP值和校正项
            logger.info("Step 4: Calculating binding free energy corrections")
            results = self.calculate_binding_correction()
            
            # 6. 保存结果
            logger.info("Step 5: Saving results")
            self.save_correction_results(results)
            self.save_details_results(results)
            
            logger.info("Workflow completed successfully")
            return results
            
        except Exception as e:
            logger.error(f"Error in workflow: {str(e)}")
            raise

def main():
    """主函数"""
    import argparse
    
    # 创建命令行参数解析器
    parser = argparse.ArgumentParser(
        description='Calculate charge corrections for binding free energy',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s com.prmcrd com.prmtop lig.prmcrd lig.prmtop -1 3
  %(prog)s com.prmcrd com.prmtop lig.prmcrd lig.prmtop -1 3 -r LIG -o results.csv
        """
    )
    
    # 必需的位置参数
    parser.add_argument('complex_crd', help='Complex coordinate file')
    parser.add_argument('complex_top', help='Complex topology file')
    parser.add_argument('ligand_crd', help='Ligand coordinate file')
    parser.add_argument('ligand_top', help='Ligand topology file')
    parser.add_argument('chg_mol', type=int, help='Ligand charge')
    parser.add_argument('chg_com', type=int, help='Complex charge')
    
    # 可选参数
    parser.add_argument('-r', '--residue', default='MOL', 
                       help='Ligand residue name (default: MOL)')
    parser.add_argument('-o', '--output', default='charge_corrections.csv',
                       help='Output file name (default: charge_corrections.csv)')
    
    # 解析命令行参数
    args = parser.parse_args()
    
    try:
        # 创建ChargeCorrection实例
        charge_corr = ChargeCorrection(
            complex_crd=args.complex_crd,
            complex_top=args.complex_top,
            ligand_crd=args.ligand_crd,
            ligand_top=args.ligand_top,
            ligand_residue_name=args.residue,
            chg_mol=args.chg_mol,
            chg_com=args.chg_com
        )
        
        # 运行工作流
        results = charge_corr.workflow()
        
        # 保存结果到指定文件
        charge_corr.save_correction_results(results, args.output)
        
        logger.info("Calculation completed successfully")
        
    except Exception as e:
        logger.error(f"Error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()
