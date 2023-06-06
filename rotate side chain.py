# python==3.9  biopython==1.81  numpy==1.21.2
# http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html
import sys
import math
from Bio.PDB import PDBParser, PDBIO, rotaxis, Vector
import numpy as np
import warnings
from Bio import BiopythonWarning
import argparse

# 创建 ArgumentParser 对象
parser = argparse.ArgumentParser(description='这是一个对单个残基旋转侧链二面角的简单脚本，旋转后可能会发生原子碰撞等问题，请仔细使用。\n'
                                             '运行环境python==3.9 biopython==1.81 numpy==1.21.2，其他环境未经过测试。\n'
                                             '二面角定义来自于http://www.mlb.co.jp/linux/science/garlic/doc/commands/dihedrals.html。\n'
                                             '注意！此脚本可能无法正确处理残基PRO与LEU，如要处理，请仔细查看结果。\n'
                                             '可在‘文件输入与参数控制区’调整输入文件名以及需要旋转的二面角以及二面角采样间隔，也可通过 python.py -i -ri -ci -c1 -c2 -c3- c4- -c5 ，后一种可覆盖前一种的参数'
                                             '注意：采样间隔是度（°）而非弧度（rad）,且都是整数。\n')
# 添加命令行选项和参数
parser.add_argument('-i', '--input', type=str, help='input pdb file.')
parser.add_argument('-ri', '--residue_id', type=int, help='id of the rotated residue.')
parser.add_argument('-ci', '--chain_id', type=str, help='chain of the rotated residue.')
parser.add_argument('-c1', '--chi1', type=int, help='Whether to rotate chi1 dihedral angle 1 rotate 0 do not rotate.')
parser.add_argument('-c2', '--chi2', type=int, help='Whether to rotate chi2 dihedral angle 1 rotate 0 do not rotate.')
parser.add_argument('-c3', '--chi3', type=int, help='Whether to rotate chi3 dihedral angle 1 rotate 0 do not rotate.')
parser.add_argument('-c4', '--chi4', type=int, help='Whether to rotate chi4 dihedral angle 1 rotate 0 do not rotate.')
parser.add_argument('-c5', '--chi5', type=int, help='Whether to rotate chi5 dihedral angle 1 rotate 0 do not rotate.')
parser.add_argument('-itv1', '--interval_degree_chi1', type=int, help='Chi1 Dihedral Sampling Interval.')
parser.add_argument('-itv2', '--interval_degree_chi2', type=int, help='Chi2 Dihedral Sampling Interval.')
parser.add_argument('-itv3', '--interval_degree_chi3', type=int, help='Chi3 Dihedral Sampling Interval.')
parser.add_argument('-itv4', '--interval_degree_chi4', type=int, help='Chi4 Dihedral Sampling Interval.')
parser.add_argument('-itv5', '--interval_degree_chi5', type=int, help='Chi5 Dihedral Sampling Interval.')
# 解析命令行选项和参数
args = parser.parse_args()

########################## 文件输入与参数控制区 ################################
pdb_file = 'test.pdb' #输入文件
residue_id = 405                  #残基号
chain_id = 'A'                    #残基所在的链
chi1 = 1                          #是否旋转chi1二面角 1旋转 0不旋转
chi2 = 0                          #是否旋转chi2二面角 1旋转 0不旋转
chi3 = 0                          #是否旋转chi3二面角 1旋转 0不旋转
chi4 = 0                          #是否旋转chi4二面角 1旋转 0不旋转
chi5 = 0                          #是否旋转chi5二面角 1旋转 0不旋转
interval_degree_chi1 = 50         #chi1二面角采样间隔
interval_degree_chi2 = 50         #chi2二面角采样间隔
interval_degree_chi3 = 50         #chi3二面角采样间隔
interval_degree_chi4 = 50         #chi4二面角采样间隔
interval_degree_chi5 = 50         #chi5二面角采样间隔
###############################################################
if args.input is not None  :
    if args.residue_id is None or args.chain_id is None or args.interval_degree_chi5 is None or args.interval_degree_chi4 is None or args.interval_degree_chi3 is None or args.interval_degree_chi2 is None or args.chi1 is None or args.chi5 is None or args.chi4 is None or args.chi3 is None :
        print("参数不全。")
        sys.exit(1)
    pdb_file = args.input  # 输入文件
    residue_id = args.residue_id  # 残基号
    chain_id = args.chain_id  # 残基所在的链
    chi1 = args.chi1  # 是否旋转chi1二面角 1旋转 0不旋转
    chi2 = args.chi2  # 是否旋转chi2二面角 1旋转 0不旋转
    chi3 = args.chi3  # 是否旋转chi3二面角 1旋转 0不旋转
    chi4 = args.chi4  # 是否旋转chi4二面角 1旋转 0不旋转
    chi5 = args.chi5  # 是否旋转chi5二面角 1旋转 0不旋转
    interval_degree_chi1 = args.interval_degree_chi1  # chi1二面角采样间隔
    interval_degree_chi2 = args.interval_degree_chi2  # chi2二面角采样间隔
    interval_degree_chi3 = args.interval_degree_chi3  # chi3二面角采样间隔
    interval_degree_chi4 = args.interval_degree_chi4  # chi4二面角采样间隔
    interval_degree_chi5 = args.interval_degree_chi5  # chi5二面角采样间隔
    print(chi1)

#################################################################

def calc_dihedral(p, q, r, s):
    pq = q - p
    qr = r - q
    rs = s - r

    a = pq.normalized()
    b = qr.normalized()
    c = rs.normalized()
    a1 = np.array([a[0], a[1], a[2]])
    b1 = np.array([b[0], b[1], b[2]])
    c1 = np.array([c[0], c[1], c[2]])

    # 法向量
    n1 = np.cross(a1, b1)
    n2 = np.cross(b1, c1)
    # 侧向量
    m1 = np.cross(n1, b1)

    x = n1.dot(n2)
    y = m1.dot(n2)

    return math.atan2(y, x)

def read_pdb_structure(pdb_file):
    # 禁用 BiopythonWarning 警告信息
    warnings.simplefilter('ignore', BiopythonWarning)
    # 解析PDB文件
    parser = PDBParser()
    structure = parser.get_structure('input_structure', pdb_file)

    # 去掉氢原子
    for a in structure.get_atoms():
        if a.element == "H":
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            if atom.element == "H":
                                residue.detach_child(atom.get_id())
    return structure

def rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi_angle_degrees, no_rotation_atom, rotation_axis_atom):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    n1 = residue[rotation_axis_atom[0]].get_vector()
    n2 = residue[rotation_axis_atom[1]].get_vector()

    # 计算要旋转的二面角（弧度）
    rotate_chi_angle_rad = math.radians(rotate_chi_angle_degrees)
    rotation_angle = rotate_chi_angle_rad
    # 旋转轴是定义为从前一个chi原子(n2)指向当前chi原子(n1)的归一化向量。
    rotation_axis = (n2 - n1).normalized()

    # 对残基中的每个原子进行循环，并使用rotaxis函数应用旋转。rotaxis函数接受一个角度和一个旋转轴，
    # 返回一个旋转矩阵，用于围绕该轴旋转向量。使用dot乘积将旋转矩阵应用于表示原子坐标相对于前一个chi原子的向量(atom.get_vector() - ca)，
    # 然后将结果坐标平移回前一个chi原子的新位置(ca + ...)。
    # if语句检查原子是否为骨架原子(N、CA、C、O)或用于定义旋转轴的侧链原子(CB)。
    # 如果原子不是这些原子之一，则意味着它是需要随侧链一起旋转的侧链原子。然后使用set_coord函数将原子的坐标更新为旋转后的位置
    for atom in residue:
        if atom.get_name() not in no_rotation_atom:
            # 转为np矩阵
            a = atom.get_vector() - n1
            b = np.array([a[0], a[1], a[2]])
            c = rotaxis(rotation_angle, rotation_axis).dot(b)
            d = Vector(c)
            atom.set_coord(n1 + d)
    return structure

def SaveStructure(structure, pdb_file, chain_id, rotate_chi_angle_degrees):
    # 保存旋转后的结构
    io = PDBIO()
    io.set_structure(structure)
    output_file = f"{pdb_file[0:-4]}_rotated_{chain_id}_{residue_id}_{rotate_chi_angle_degrees}.pdb"
    io.save(output_file)
    print(f"旋转{rotate_chi_angle_degrees}°后的结构已保存到 {output_file}")

def Chi1(structure, residue_id, chain_id):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    # 检查残基是否有侧链
    if residue.resname == 'GLY' or residue.resname == 'ALA':
        print("残基为GLY或者ALA，侧链无法旋转。")
        sys.exit(1)
    chi1_rotation_axis_atom = ['CA', 'CB']
    chi1_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB']

    return chi1_no_rotation_atom,chi1_rotation_axis_atom

def Chi2(structure, residue_id, chain_id):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    # 检查残基是否有chi2
    if residue.resname in ['ALA', 'CYS', 'GLY', 'SER', 'THR', 'VAL']:
        print(f"残基为{residue.resname}，无chi2二面角,只有ARG,ASN,ASP,GLN,GLU,HIS,HIE,LEU,LYS,MET,PHE,PRO,TRP,TYR有chi2二面角")
        sys.exit(1)

    if residue.resname == 'ILE':
        chi2_rotation_axis_atom = ['CB', 'CG1']
    else:
        chi2_rotation_axis_atom = ['CB', 'CG']

    if residue.resname in ['ARG','ASN','ASP','GLN','GLU','HIS','LEU','LYS','MET','PHE','TRP','TYR','PRO']:
        chi2_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG']
    if residue.resname in ['ILE']:
        chi2_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG1','CG2']

    return chi2_no_rotation_atom, chi2_rotation_axis_atom

def Chi3(structure, residue_id, chain_id):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    # 检查残基是否有chi3
    if residue.resname not in ['ARG','GLN','GLU','LYS','MET']:
        print(f"残基是{residue.resname}，无chi3二面角，只有ARG,GLN,GLU,LYS,MET有chi3二面角")
        sys.exit(1)

    if residue.resname == 'MET':
        chi3_rotation_axis_atom = ['CG', 'SD']
    else:
        chi3_rotation_axis_atom = ['CG', 'CD']

    if residue.resname in ['ARG','GLN','GLU','LYS']:
        chi3_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG','CD']
    if residue.resname in ['MET']:
        chi3_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG','SD']

    return chi3_no_rotation_atom, chi3_rotation_axis_atom

def Chi4(structure, residue_id, chain_id):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    # 检查残基是否有chi4
    if residue.resname not in ['ARG','LYS']:
        print(f"残基是{residue.resname}，无chi4二面角，只有ARG,LYS有chi4二面角")
        sys.exit(1)

    if residue.resname == 'ARG':
        chi4_rotation_axis_atom = ['CD', 'NE']
    if residue.resname == 'LYS':
        chi4_rotation_axis_atom = ['CD', 'CE']

    if residue.resname == 'ARG':
        chi4_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG','CD','NE']
    if residue.resname == 'LYS':
        chi4_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB','CG','CD', 'CE']

    return chi4_no_rotation_atom, chi4_rotation_axis_atom

def Chi5(structure, residue_id, chain_id):
    # 获取指定的残基
    chain = structure[0][chain_id]
    residue = chain[residue_id]

    # 检查残基是否有chi5
    if residue.resname != 'ARG':
        print(f"残基是{residue.resname}，无chi5二面角,只有ARG有chi5二面角。")
        sys.exit(1)

    chi5_rotation_axis_atom = ['NE','CZ']

    chi5_no_rotation_atom = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE','CZ']

    return chi5_no_rotation_atom, chi5_rotation_axis_atom




if chi1 == 0:
    interval_degree_chi1=361
if chi2 == 0:
    interval_degree_chi2=361
if chi3 == 0:
    interval_degree_chi3=361
if chi4 == 0:
    interval_degree_chi4=361
if chi5 == 0:
    interval_degree_chi5=361

chis = []
if chi1 == 1:
    chis.append('chi1')
if chi2 == 1:
    chis.append('chi2')
if chi3 == 1:
    chis.append('chi3')
if chi4 == 1:
    chis.append('chi4')
if chi5 == 1:
    chis.append('chi5')
structure = read_pdb_structure(pdb_file)
chain = structure[0][chain_id]
residue = chain[residue_id]

for rotate_chi1_angle_degrees in range(0, 361, interval_degree_chi1):
    if 'chi1' in chis:
        chi1_no_rotation_atom, chi1_rotation_axis_atom = Chi1(structure, residue_id, chain_id)
        structure = rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi1_angle_degrees, chi1_no_rotation_atom, chi1_rotation_axis_atom)
    for rotate_chi2_angle_degrees in range(0, 361, interval_degree_chi2):
        if 'chi2' in chis:
            chi2_no_rotation_atom, chi2_rotation_axis_atom = Chi2(structure, residue_id, chain_id)
            structure = rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi2_angle_degrees, chi2_no_rotation_atom, chi2_rotation_axis_atom)
        for rotate_chi3_angle_degrees in range(0, 361, interval_degree_chi3):
            if 'chi3' in chis:
                chi3_no_rotation_atom, chi3_rotation_axis_atom = Chi3(structure, residue_id, chain_id)
                structure = rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi3_angle_degrees, chi3_no_rotation_atom, chi3_rotation_axis_atom)
            for rotate_chi4_angle_degrees in range(0, 361, interval_degree_chi4):
                if 'chi4' in chis:
                    chi4_no_rotation_atom, chi4_rotation_axis_atom = Chi4(structure, residue_id, chain_id)
                    structure = rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi4_angle_degrees, chi4_no_rotation_atom, chi4_rotation_axis_atom)
                for rotate_chi5_angle_degrees in range(0, 361, interval_degree_chi5):
                    if 'chi5' in chis:
                        chi5_no_rotation_atom, chi5_rotation_axis_atom = Chi5(structure, residue_id, chain_id)
                        structure = rotate_sidechain_chi(structure, residue_id, chain_id, rotate_chi5_angle_degrees, chi5_no_rotation_atom, chi5_rotation_axis_atom)
                    SaveStructure(structure, pdb_file, chain_id, f"chi1_{rotate_chi1_angle_degrees}_chi2_{rotate_chi2_angle_degrees}_chi3_{rotate_chi3_angle_degrees}_chi4_{rotate_chi4_angle_degrees}_chi5_{rotate_chi5_angle_degrees}")
