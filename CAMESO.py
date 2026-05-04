import csv
from collections import defaultdict
from tqdm import tqdm
import pandas as pd
import os
import sys
from concurrent.futures import ThreadPoolExecutor
import multiprocessing
class Catanormcalculator:
    def __init__(self):
        # ========== 1. 分子量定义（阳离子形式）==========
        self.molecular_weights = {
            # 普通氧化物（直接使用）
            'SiO2': 60.08, 'TiO2': 79.90, 'FeO': 71.85, 'MnO': 70.94,
            'MgO': 40.30, 'CaO': 56.08, 'BaO': 153.33, 'SrO': 103.62,
            'NiO': 74.69, 'ZrO2': 123.22, 'SnO2': 150.71, 'CO2': 44.01,
            'H2O': 18.02, 'F': 19.00, 'Cl': 35.45, 'S': 32.07,
            
            # 阳离子形式（已预转换）
            'PO2.5': 70.97,   # P⁵⁺（P2O5→2PO2.5）
            'CrO1.5': 76.00,  # Cr³⁺（Cr2O3→2CrO1.5）
            'AlO1.5': 50.98,  # Al³⁺（Al2O3→2AlO1.5）
            'FeO1.5': 79.85,  # Fe³⁺（Fe2O3→2FeO1.5）
            'NaO0.5': 30.99,  # Na⁺（Na2O→2NaO0.5）
            'KO0.5': 47.10    # K⁺（K2O→2KO0.5）
        }

        # 氧化物到阳离子的映射
        self.oxide_to_cation = {
            'P2O5': 'PO2.5',
            'Cr2O3': 'CrO1.5',
            'Al2O3': 'AlO1.5',
            'Fe2O3': 'FeO1.5',
            'Na2O': 'NaO0.5',
            'K2O': 'KO0.5'
        }
        # ========== 2. 数据存储初始化 ==========
        self.oxide_wt_percent = defaultdict(float)  # 存储转换后的组分（如PO2.5）
        self.cation_proportions = defaultdict(float)  # 阳离子毫摩尔量
        self.minerals = defaultdict(float)  # 生成的矿物量
        self.desilication_steps = []  # 脱硅化步骤记录
        self.negative_quartz = 0.0  # 负石英量

    # ========== 3. 映射转换 ==========
    def convert_oxides_to_cations(self):
        converted = defaultdict(float)

        for oxide, wt in self.oxide_wt_percent.items():
            if oxide in self.oxide_to_cation:
                cation = self.oxide_to_cation[oxide]

                
                # 保持质量一致性（质量换算）
                converted[cation] += wt 
            else:
                converted[oxide] += wt  # 直接保留普通氧化物

        self.oxide_wt_percent = converted   

    # ========== 4. 阳离子比例计算 ==========
    def calculate_cation_proportions(self):
        """根据氧化物重量百分数，计算初始阳离子毫摩尔量"""
        self.cation_proportions = defaultdict(float)
        for component, wt_percent in self.oxide_wt_percent.items():
            if wt_percent > 0:
                self.cation_proportions[component] = (wt_percent / self.molecular_weights[component]) * 1000



    def convert_cations_to_percentage(self):
        """将阳离子数转换为阳离子百分比"""
        total = sum(self.cation_proportions.values())
        if total == 0:
            return

        for component, value in self.cation_proportions.items():
            self.cation_proportions[component] = (value / total) * 100

    def merge_cation_equivalents(self):
        """将某些氧化物合并为等效阳离子"""
        # Fe²⁺ = FeO + MnO + NiO
        if 'MnO' in self.cation_proportions:
            self.cation_proportions['FeO'] += self.cation_proportions['MnO']
            del self.cation_proportions['MnO']
        if 'NiO' in self.cation_proportions:
            self.cation_proportions['FeO'] += self.cation_proportions['NiO']
            del self.cation_proportions['NiO']

        # Ca²⁺ = CaO + BaO + SrO
        if 'BaO' in self.cation_proportions:
            self.cation_proportions['CaO'] += self.cation_proportions['BaO']
            del self.cation_proportions['BaO']
        if 'SrO' in self.cation_proportions:
            self.cation_proportions['CaO'] += self.cation_proportions['SrO']
            del self.cation_proportions['SrO']

    # ========== 5. 矿物计算步骤 ==========
    def calculate_calcite_and_cassiterite(self):
        """方解石(Cc)和锡石(Ct)"""
        if 'CO2' in self.cation_proportions:
            self.minerals['Cc'] = 2 * self.cation_proportions['CO2']
            self.cation_proportions['CaO'] -= self.cation_proportions['CO2']
        
        if 'SnO2' in self.cation_proportions:
            self.minerals['Ct'] = self.cation_proportions['SnO2']
            del self.cation_proportions['SnO2']

    def calculate_accessory_minerals(self):
        """副矿物计算"""
        # 磷灰石(Ap) - 使用PO2.5
        if 'PO2.5' in self.cation_proportions:
            po = self.cation_proportions['PO2.5']
            f = self.cation_proportions.get('F', 0)
            
            if po <= 3 * f:
                self.minerals['Ap'] = 3 * po
                self.cation_proportions['CaO'] -= 1.667 * po
                self.cation_proportions['F'] -= 0.333 * po
            else:
                self.minerals['Ap'] = (2.667 * po) + f
                self.cation_proportions['CaO'] -= 1.667 * po
                self.cation_proportions['F'] = 0
        
        # 萤石(Fr)
        if 'F' in self.cation_proportions and self.cation_proportions['F'] > 0:
            self.minerals['Fr'] = 1.5 * self.cation_proportions['F']
            self.cation_proportions['CaO'] -= 0.5 * self.cation_proportions['F']
            self.cation_proportions['F'] = 0
        
        # 石盐(Hl)
        if 'Cl' in self.cation_proportions and 'NaO0.5' in self.cation_proportions:
            self.minerals['Hl'] = 2 * self.cation_proportions['Cl']
            self.cation_proportions['NaO0.5'] -= self.cation_proportions['Cl']
            self.cation_proportions['Cl'] = 0
        
        # 黄铁矿(Pr)
        if 'S' in self.cation_proportions:
            self.minerals['Pr'] = 1.5 * self.cation_proportions['S']
            self.cation_proportions['FeO'] -= 0.5 * self.cation_proportions['S']
            self.cation_proportions['S'] = 0
        
        # 铬铁矿(Cm)
        if 'CrO1.5' in self.cation_proportions:
            self.minerals['Cm'] = 1.5 * self.cation_proportions['CrO1.5']
            self.cation_proportions['FeO'] -= 0.5 * self.cation_proportions['CrO1.5']
            self.cation_proportions['CrO1.5'] = 0
        
        # 锆石(Z)
        if 'ZrO2' in self.cation_proportions:
            self.minerals['Z'] = 2 * self.cation_proportions['ZrO2']
            self.cation_proportions['SiO2'] -= self.cation_proportions['ZrO2']
            self.cation_proportions['ZrO2'] = 0
        
        # 钛铁矿(Il)
    def calculate_ilmenite(self):
        """计算钛铁矿（Il）"""
        if 'TiO2' in self.cation_proportions and 'FeO' in self.cation_proportions:
            tio2 = self.cation_proportions['TiO2']
            feo = self.cation_proportions['FeO']
            ilmenite = 2 * min(feo, tio2)
            self.minerals['Il'] = ilmenite
            if feo >= tio2:
                self.cation_proportions['FeO'] -= tio2
                self.cation_proportions['TiO2'] = 0
            else:
                self.cation_proportions['TiO2'] -= feo
                self.cation_proportions['FeO'] = 0
        
    def calculate_potassium_feldspar(self):
            """钾长石(Or)和钾硅酸盐(Ks)计算"""
            if 'KO0.5' in self.cation_proportions and 'AlO1.5' in self.cation_proportions:
                ko = self.cation_proportions['KO0.5']
                al = self.cation_proportions['AlO1.5']
                
                # 1. 钾长石(Or) = 5 × min(KO0.5, AlO1.5)
                or_amount = 5 * min(ko, al)
                if or_amount > 0:
                    self.minerals['Or'] = or_amount
                    consumed = min(ko, al)
                    self.cation_proportions['AlO1.5'] -= consumed
                    self.cation_proportions['KO0.5'] -= consumed
                    self.cation_proportions['SiO2'] -= 3 * consumed
                
                # 2. 剩余KO0.5转为钾硅酸盐(Ks)
                if 'KO0.5' in self.cation_proportions and self.cation_proportions['KO0.5'] > 0:
                    ko_remaining = self.cation_proportions['KO0.5']
                    self.minerals['Ks'] = 1.5 * ko_remaining
                    self.cation_proportions['SiO2'] -= 0.5 * ko_remaining
                    self.cation_proportions['KO0.5'] = 0

    def calculate_sodium_feldspar(self):
        """钠长石(Ab)计算"""
        if 'NaO0.5' in self.cation_proportions and 'AlO1.5' in self.cation_proportions:
            nao = self.cation_proportions['NaO0.5']
            al = self.cation_proportions['AlO1.5']
            
            ab_amount = 5 * min(nao, al)
            if ab_amount > 0:
                self.minerals['Ab'] = ab_amount
                consumed = min(nao, al)
                self.cation_proportions['AlO1.5'] -= consumed
                self.cation_proportions['NaO0.5'] -= consumed
                self.cation_proportions['SiO2'] -= 3 * consumed

    def calculate_sodium_silicate(self):
        """钠硅酸盐(Ns)计算"""
        if 'NaO0.5' in self.cation_proportions and self.cation_proportions['NaO0.5'] > 0:
            nao = self.cation_proportions['NaO0.5']
            self.minerals['Ns'] = 1.5 * nao
            self.cation_proportions['SiO2'] -= 0.5 * nao
            self.cation_proportions['NaO0.5'] = 0

    def calculate_anorthite(self):
        """钙长石(An)计算"""
        if 'CaO' in self.cation_proportions and 'AlO1.5' in self.cation_proportions:
            cao = self.cation_proportions['CaO']
            al = self.cation_proportions['AlO1.5']
            
            if cao <= 0.5 * al:
                an_amount = 5 * cao
                self.minerals['An'] = an_amount
                self.cation_proportions['AlO1.5'] -= 2 * cao
                self.cation_proportions['SiO2'] -= 2 * cao
                self.cation_proportions['CaO'] = 0
            else:
                an_amount = 2.5 * al
                self.minerals['An'] = an_amount
                self.cation_proportions['CaO'] -= 0.5 * al
                self.cation_proportions['SiO2'] -= al
                self.cation_proportions['AlO1.5'] = 0

    def calculate_titanite_and_rutile(self):
       #榍石(Tn)和金红石(Ru)"""
        if 'TiO2' in self.cation_proportions and 'CaO' in self.cation_proportions:
            tio2 = self.cation_proportions['TiO2']
            cao = self.cation_proportions['CaO']
            
            if tio2 <= cao:
                self.minerals['Tn'] = 3 * tio2
                self.cation_proportions['CaO'] -= tio2
                self.cation_proportions['SiO2'] -= tio2
                self.cation_proportions['TiO2'] = 0
            else:
                self.minerals['Tn'] = 3 * cao
                self.cation_proportions['TiO2'] -= cao
                self.cation_proportions['SiO2'] -= cao
                self.cation_proportions['CaO'] = 0
        
        # 金红石(Ru) = 剩余TiO2
        if 'TiO2' in self.cation_proportions and self.cation_proportions['TiO2'] > 0:
            self.minerals['Ru'] = self.cation_proportions['TiO2']
            self.cation_proportions['TiO2'] = 0

    def calculate_corundum(self):
        """刚玉(Cor)"""
        if 'AlO1.5' in self.cation_proportions and self.cation_proportions['AlO1.5'] > 0:
            self.minerals['Cor'] = self.cation_proportions['AlO1.5']
            self.cation_proportions['AlO1.5'] = 0

    def calculate_iron_oxides(self):
        #铁氧化物（磁铁矿Mt和赤铁矿Hm）"""
        if 'FeO1.5' in self.cation_proportions and 'FeO' in self.cation_proportions:
            fe3 = self.cation_proportions['FeO1.5']  # Fe³⁺
            fe2 = self.cation_proportions['FeO']     # Fe²⁺
            
            if fe3 <= 2 * fe2:
                self.minerals['Mt'] = 1.5 * fe3
                self.cation_proportions['FeO'] -= 0.5 * fe3
                self.cation_proportions['FeO1.5'] = 0
            else:
                self.minerals['Mt'] = 3 * fe2
                self.cation_proportions['FeO1.5'] -= 2 * fe2
                self.cation_proportions['FeO'] = 0
        
        # 赤铁矿(Hm) = 剩余FeO1.5
        if 'FeO1.5' in self.cation_proportions and self.cation_proportions['FeO1.5'] > 0:
            self.minerals['Hm'] = self.cation_proportions['FeO1.5']
            self.cation_proportions['FeO1.5'] = 0

    def calculate_wollastonite(self):
        """硅灰石(Wo) = 所有剩余CaO，按化学式 CaSiO3"""
        if 'CaO' in self.cation_proportions and self.cation_proportions['CaO'] > 0:
            cao = self.cation_proportions['CaO']
            self.minerals['Wo'] = 2 * cao  # 化学计量系数：1 CaO → 2 Wo
            self.cation_proportions['SiO2'] -= cao    # 消耗等量SiO2
            self.cation_proportions['CaO'] = 0        # CaO全部耗尽

    def calculate_pyroxenes(self):
        """辉石类矿物：顽火辉石(En)处理MgO，铁辉石(Fs)处理FeO"""
        en = 0.0
        fs = 0.0

        if 'MgO' in self.cation_proportions and self.cation_proportions['MgO'] > 0:
            mgo = self.cation_proportions['MgO']
            en = 2 * mgo
            self.cation_proportions['SiO2'] -= mgo
            self.cation_proportions['MgO'] = 0

        if 'FeO' in self.cation_proportions and self.cation_proportions['FeO'] > 0:
            feo = self.cation_proportions['FeO']
            fs = 2 * feo
            self.cation_proportions['SiO2'] -= feo
            self.cation_proportions['FeO'] = 0

        if en + fs > 0:
            self.minerals['Hy'] = en + fs
            x_en = en / (en + fs) * 100
            self.x_en = x_en

            
 

    def calculate_diopside(self):
        """透辉石(Di)"""
        if 'Hy' in self.minerals and 'Wo' in self.minerals:
            hy = self.minerals['Hy']
            wo = self.minerals['Wo']
            
            if hy < wo:
                self.minerals['Di'] = 2 * hy
                self.minerals['Wo'] -= hy
                del self.minerals['Hy']
            else:
                self.minerals['Di'] = 2 * wo
                self.minerals['Hy'] -= wo
                del self.minerals['Wo']

    def calculate_quartz(self):
        """石英(Q)和脱硅化处理"""
        if 'SiO2' in self.cation_proportions:
            sio2 = self.cation_proportions['SiO2']
            if sio2 > 0:
                self.minerals['Q'] = sio2
            else:
                self.negative_quartz = abs(sio2)
                self.handle_desilication()

    def handle_desilication(self):
        """处理负石英（按优先级替代矿物）"""
        # 脱硅化优先级顺序
        steps = [
            ('Hy', self.desilicate_hypersthene),
            ('Tn', self.desilicate_titanite),
            ('Ab', self.desilicate_albite),
            ('Or', self.desilicate_orthoclase),
            ('Lc', self.desilicate_leucite),
            ('Wo', self.desilicate_wollastonite),
            ('Di', self.desilicate_diopside)
        ]
        
        for mineral, func in steps:
            if self.negative_quartz <= 0:
                break
            if mineral in self.minerals and self.minerals[mineral] > 0:
                func()

    # ========== 6. 脱硅化具体步骤 ==========
    def desilicate_hypersthene(self):
        """替代紫苏辉石(Hy → Ol)"""
        hy = self.minerals['Hy']
        d = self.negative_quartz
        
        if d <= hy / 4:
            ol = 3 * d
            self.minerals['Ol'] = ol
            self.minerals['Hy'] = hy - 4 * d
            self.negative_quartz = 0
            
        else:
            ol = (3/4) * hy
            self.minerals['Ol'] = ol
            self.negative_quartz = d - hy / 4
            del self.minerals['Hy']
           

    def desilicate_titanite(self):
        """替代榍石(Tn → Pf)"""
        tn = self.minerals['Tn']
        d = self.negative_quartz
        
        if d <= tn / 3:
            pf = 2 * d
            self.minerals['Pf'] = pf
            self.minerals['Tn'] = tn - 3 * d
            self.negative_quartz = 0
            
        else:
            pf = (2/3) * tn
            self.minerals['Pf'] = pf
            self.negative_quartz = d - tn / 3
            del self.minerals['Tn']
           

    # ...（其他脱硅化方法，保持原逻辑）
    def desilicate_albite(self):
        """替代钠长石 (Ab → Ne)"""
        if 'Ab' in self.minerals:
            ab = self.minerals['Ab']
            d = self.negative_quartz

            if ab >= 2.5 * d:
                # Ab 足够完全反应
                ne = 1.5 * d
                self.minerals['Ne'] = ne
                self.minerals['Ab'] = ab - 2.5 * d
                self.negative_quartz = 0
                
            else:
                # Ab 不够，部分反应
                ne = 0.6 * ab
                self.minerals['Ne'] = ne
                self.negative_quartz = d - 0.4 * ab
                del self.minerals['Ab']
                


    def desilicate_orthoclase(self):
        """脱硅替代：钾长石(Or) → 白榴石(Lc)，不累加、不记录日志"""

        or_val = self.minerals.get('Or', 0)
        d = self.negative_quartz

        if or_val <= 0 or d <= 0:
            return

        if d <= or_val / 5:
            # Or 足够
            lc = 4 * d
            self.minerals['Lc'] = lc
            self.minerals['Or'] = or_val - 5 * d
            self.negative_quartz = 0

            if self.minerals['Or'] <= 1e-6:
                del self.minerals['Or']

        else:
            # Or 不足，部分脱硅
            lc = 0.8 * or_val
            self.minerals['Lc'] = lc
            self.negative_quartz = d - 0.2 * or_val
            del self.minerals['Or']


    def desilicate_leucite(self):
        """替代白榴石(Lc → Kp)"""
        if 'Lc' in self.minerals:
            lc = self.minerals['Lc']
            d = self.negative_quartz
            
            if d <= lc / 4:
                # D 小于或等于 Lc/4，生成钾霞石 Kp = 3 * D
                kp = 3 * d
                self.minerals['Kp'] = kp
                self.minerals['Lc'] = lc - 4 * d
                self.negative_quartz = 0
            else:
                # D 大于 Lc/4，生成钾霞石 Kp = 0.75 * Lc
                kp = 0.75 * lc
                self.minerals['Kp'] = kp
                self.negative_quartz = d - 0.25 * lc
                del self.minerals['Lc']

    def desilicate_wollastonite(self):
        """替代硅灰石(Wo → Cs)"""
        if 'Wo' in self.minerals:
            wo = self.minerals['Wo']
            d = self.negative_quartz
            
            if d <= wo / 4:
                cs = 3 * d
                self.minerals['Cs'] = cs
                self.minerals['Wo'] = wo - 4 * d
                self.negative_quartz = 0
               
            else:
                cs = 0.75 * wo
                self.minerals['Cs'] = cs
                self.negative_quartz = d - 0.25 * wo
                del self.minerals['Wo']
                

    def desilicate_diopside(self):
        """替代透辉石(Di → Cs + Ol)，累加已生成的 Cs 和 Ol"""

        di = self.minerals.get('Di', 0)
        d = self.negative_quartz

        if di <= 0 or d <= 0:
            return

        if d <= di / 4:
            # 分支 1：Di 足够
            delta = 1.5 * d
            self.minerals['Cs'] = self.minerals.get('Cs', 0) + delta
            self.minerals['Ol'] = self.minerals.get('Ol', 0) + delta
            self.minerals['Di'] = di - 4 * d
            self.negative_quartz = 0

            if self.minerals['Di'] <= 1e-6:
                del self.minerals['Di']

        else:
            # 分支 2：Di 不足
            delta = 0.375 * di
            self.minerals['Cs'] = self.minerals.get('Cs', 0) + delta
            self.minerals['Ol'] = self.minerals.get('Ol', 0) + delta
            self.negative_quartz = d - 0.25 * di
            del self.minerals['Di']

    def split_pyroxenes_from_hy(self):
        """如果存在紫苏辉石(Hy)，则拆分为顽火辉石(En)和铁辉石(Fs)，并删除Hy。"""
        if 'Hy' not in self.minerals:
            return  # 没有紫苏辉石，无需处理

        hy_total = self.minerals['Hy']
        
        # 获取形成Hy时的en和fs比例（可用已有变量或重新计算）
        # 由于原始en/(en+fs)*100 = x_en，可以反推出比例
        # 假设你保存过en和fs比例，可以从其他属性中拿
        # 如果没有保存过比例，这里必须有一个默认值或另设传参机制

        # 示例：这里假设你在形成Hy时保存了x_en值：
        x_en = getattr(self, 'x_en', None)
        if x_en is None:
            return  # 没有x_en，无法拆分

        en_ratio = x_en / 100
        fs_ratio = 1 - en_ratio

        en_amount = hy_total * en_ratio
        fs_amount = hy_total * fs_ratio

        # 累加到En和Fs中
        self.minerals['En'] = self.minerals.get('En', 0) + en_amount
        self.minerals['Fs'] = self.minerals.get('Fs', 0) + fs_amount

        # 删除Hy
        del self.minerals['Hy']
               

    # ========== 7. 主计算流程 ==========
    def calculate_all_minerals(self):
        """执行完整的矿物计算流程"""
        self.convert_oxides_to_cations()
        self.calculate_cation_proportions()
        self.convert_cations_to_percentage()
        self.merge_cation_equivalents()       # 1. 计算阳离子比例
        self.calculate_calcite_and_cassiterite()  # 2. 方解石/锡石
        self.calculate_accessory_minerals()       # 3. 副矿物
        self.calculate_ilmenite()                 # 钛铁矿
        self.calculate_potassium_feldspar()
        self.calculate_sodium_feldspar()
        self.calculate_sodium_silicate()
        self. calculate_anorthite()               # 4. 长石类
        self.calculate_titanite_and_rutile()      # 5. 榍石/金红石
        self.calculate_corundum()                # 6. 刚玉
        self.calculate_iron_oxides()             # 7. 铁氧化物
        self.calculate_wollastonite()            # 8. 硅灰石
        self.calculate_pyroxenes()               # 9. 辉石类
        self.calculate_diopside()                # 10. 透辉石
        self.calculate_quartz()                  # 11. 石英/脱硅化
        self.split_pyroxenes_from_hy()           #拆分紫苏辉石

    






class MesonormCalculator(Catanormcalculator):
    def __init__(self, is_mafic=True):
        #:param is_mafic: True表示中性岩(走路径A)，False表示基性岩(走路径B)"""
        super().__init__()
        self.is_mafic = is_mafic
        self.remaining_Al = 0.0  # Al'
        self.remaining_Ca = 0.0  # Ca'
        self.minerals['Act'] = 0.0  # 阳起石
        self.minerals['Ed'] = 0.0   # 浅闪石
        self.minerals['Ri'] = 0.0   # 钠闪石
        self.minerals['Sp'] = 0.0   # 尖晶石

    def calculate_all_minerals(self):
        """Mesonorm特有计算流程"""
        self.convert_oxides_to_cations()
        # 1. 公共初始步骤（与Catanorm相同）
        self.calculate_cation_proportions()
        self.convert_cations_to_percentage()
        self.merge_cation_equivalents()
        # 2. 相同计算的矿物（直接继承）
        self.calculate_calcite_and_cassiterite()
        self.calculate_accessory_minerals()
        
        # 3. Mesonorm特有榍石计算
        self.calculate_titanite_mesonorm()  
        
        # 4. 相同的长石计算
        self.calculate_potassium_feldspar()
        self.calculate_sodium_feldspar()
        
        # 5. 钠闪石计算（新增）
        self.calculate_riebeckite()
        
        # 6. 相同的钠硅酸盐计算
        self.calculate_sodium_silicate()
        
        # 7. 相同的铁氧化物计算
        self.calculate_iron_oxides()
        
        # 8. 合并FeO和MgO并记录Al'和Ca'
        self.prepare_for_amphibole()
        
        # 9. 相同的钙长石和刚玉计算
        self.calculate_anorthite()
        self.calculate_corundum()
        
        # 10. 钾长石转黑云母（新增）
        self.convert_orthoclase_to_biotite()
        
        # 11. 路径选择
        if self.is_mafic:
            self.calculate_path_a()  # 中性岩路径
        else:
            self.calculate_path_b()  # 基性岩路径
        
        # 12. 硅饱和检查与脱硅化
        self.check_silica_saturation()
        
       

    # ========== Mesonorm特有方法 ==========
    def calculate_titanite_mesonorm(self):
        """Mesonorm榍石计算（直接使用TiO₂）"""
        tio2 = self.cation_proportions.get('TiO2', 0)
        if tio2 > 0:
            self.minerals['Tn'] = 3 * tio2
            self.cation_proportions['CaO'] -= tio2
            self.cation_proportions['SiO2'] -= tio2
            self.cation_proportions['TiO2'] = 0

    def calculate_riebeckite(self):
        """钠闪石 (Ri) 计算：基于 NaO0.5 和 FeO1.5，且需满足 NaO0.5 > AlO1.5"""
        na = self.cation_proportions.get('NaO0.5', 0)
        fe3 = self.cation_proportions.get('FeO1.5', 0)
        al = self.cation_proportions.get('AlO1.5', 0)  # 新增：获取 Al 含量

        # 条件：Na 和 Fe³⁺ 均需存在，且 Na > Al
        if na > 0 and fe3 > 0 and na > al:
            if na <= fe3:
                ri = 7.5 * na
                self.minerals['Ri'] = ri
                # 消耗组分
                self.cation_proportions['FeO1.5'] -= na
                self.cation_proportions['NaO0.5'] = 0
                self.cation_proportions['FeO'] -= 1.5 * na
                self.cation_proportions['SiO2'] -= 4 * na
            else:
                ri = 7.5 * fe3
                self.minerals['Ri'] = ri
                # 消耗组分
                self.cation_proportions['NaO0.5'] -= fe3
                self.cation_proportions['FeO1.5'] = 0
                self.cation_proportions['FeO'] -= 1.5 * fe3
                self.cation_proportions['SiO2'] -= 4 * fe3

    def prepare_for_amphibole(self):
       #准备角闪石计算（合并FeO+MgO，记录Al'和Ca'）"""
        # 合并FeO和MgO
        if 'MgO' in self.cation_proportions:
            self.cation_proportions['FeO'] += self.cation_proportions['MgO']
            del self.cation_proportions['MgO']
        self.remaining_Al = self.cation_proportions.get('AlO1.5', 0)
        self.remaining_Ca = self.cation_proportions.get('CaO', 0)

    def convert_orthoclase_to_biotite(self):
        """钾长石(Or) 转黑云母(Bi) - 按比例分支规则处理"""

        or_amount = self.minerals.get('Or', 0)
        mgfe = self.cation_proportions.get('FeO', 0)  # MgO已合并进FeO

        if or_amount > 0 and mgfe > 0:
            threshold = 0.6 * or_amount

            if mgfe <= threshold:
                bi = 2.667 * mgfe
                self.minerals['Bi'] = bi
                self.minerals['Or'] -= 1.667 * mgfe
                self.cation_proportions['FeO'] = 0
            else:
                bi = 1.6 * or_amount
                self.minerals['Bi'] = bi
                self.cation_proportions['FeO'] -= 0.6 * or_amount
                self.minerals['Or'] = 0

    def calculate_corundum(self):
        """刚玉(Cor)计算（mesonorm特有逻辑：AlO1.5 > CaO时才触发）"""
        al = self.cation_proportions.get('AlO1.5', 0)
        ca = self.cation_proportions.get('CaO', 0)
        
        if al > 0 and al > ca:  # 新增条件
            self.minerals['Cor'] = al
            self.cation_proportions['AlO1.5'] = 0

    def calculate_path_a(self):
        """中性岩路径计算（阳起石 → 紫苏辉石 → 硅灰石）"""

        feo = self.cation_proportions.get('FeO', 0)  # Mg+Fe 总量（已合并）
        cao = self.cation_proportions.get('CaO', 0)

        # === 1. 生成阳起石（Act） ===
        if feo > 0 and cao > 0:
            if feo <= 2.5 * cao:
                act = 3 * feo
                self.minerals['Act'] = act
                self.cation_proportions['CaO'] -= 0.4 * feo
                self.cation_proportions['SiO2'] -= 1.6 * feo
                self.cation_proportions['FeO'] = 0
            else:
                act = 7.5 * cao
                self.minerals['Act'] = act
                self.cation_proportions['FeO'] -= 2.5 * cao
                self.cation_proportions['SiO2'] -= 4 * cao
                self.cation_proportions['CaO'] = 0

        # === 2. 剩余 FeO 转为紫苏辉石（Hy）===
        feo = self.cation_proportions.get('FeO', 0)
        if feo > 0:
            self.minerals['Hy'] = 2 * feo
            self.cation_proportions['SiO2'] -= feo
            self.cation_proportions['FeO'] = 0

        # === 3. 剩余 CaO 转为硅灰石（Wo）===
        cao = self.cation_proportions.get('CaO', 0)
        if cao > 0:
            self.minerals['Wo'] = 2 * cao
            self.cation_proportions['SiO2'] -= cao
            self.cation_proportions['CaO'] = 0

    def calculate_path_b(self):
        """基性岩路径计算（复用Catanorm方法）"""
        # 直接调用父类方法
        self.calculate_wollastonite()
        self.calculate_pyroxenes()
        self.calculate_diopside()

    def check_silica_saturation(self):
        """硅饱和检查（复用父类逻辑）"""
        if 'SiO2' in self.cation_proportions:
            sio2 = self.cation_proportions['SiO2']
            if sio2 > 0:
                self.minerals['Q'] = sio2
            else:
                self.negative_quartz = abs(sio2)
                self.handle_desilication_mesonorm()

    def handle_desilication_mesonorm(self):
        """Mesonorm特有脱硅化流程"""
        while self.negative_quartz > 1e-6:  # 浮点精度控制
            original_negative = self.negative_quartz

            # 1. 替换阳起石（Act → Ed）
            if 'Act' in self.minerals:
                self.replace_actinolite()
                if self.negative_quartz <= 1e-6:
                    break

            # 2. 替换紫苏辉石（Hy → Ol）
            if 'Hy' in self.minerals:
                self.desilicate_hypersthene()
                if self.negative_quartz <= 1e-6:
                    break

            # 3. 替换橄榄石+刚玉 → 尖晶石
            if 'Ol' in self.minerals or 'Cor' in self.minerals:
                self.replace_olivine_corundum()
                if self.negative_quartz <= 1e-6:
                    break

            # 4. 替换钠长石（Ab → Ne）
            if 'Ab' in self.minerals:
                self.desilicate_albite()
                if self.negative_quartz <= 1e-6:
                    break

             # 5. 替换钾长石（Or → Lc）
            if 'Or' in self.minerals:
                self.desilicate_orthoclase()
                if self.negative_quartz <= 1e-6:
                    break

            # 6. 替换白榴石（Lc → Kp）
            if 'Lc' in self.minerals:
                self.desilicate_leucite()
                if self.negative_quartz <= 1e-6:
                    break

            # 7. 替换硅灰石（Wo → Cs）
            if 'Wo' in self.minerals:
                self.desilicate_wollastonite()
                if self.negative_quartz <= 1e-6:
                    break

            # 8. 替换透辉石（Di → Cs）
            if 'Di' in self.minerals:
                self.desilicate_diopside()
                if self.negative_quartz <= 1e-6:
                    break


           # 若没有任何变化，则中止，防死循环
            if abs(original_negative - self.negative_quartz) < 1e-6:
                break


    def replace_actinolite(self):
        """阳起石（Act）+ 钠长石（Ab） → 浅闪石（Ed），用于脱硅，部分资源可参与"""

        act = self.minerals.get('Act', 0)
        ab = self.minerals.get('Ab', 0)
        d = self.negative_quartz

        if act <= 0 or ab <= 0 or d <= 0:
            return  # 条件不足，不做任何反应

        # 最大可替换的石英量受限于资源（按完全反应的比例）
        max_q_from_act = act / 3.75
        max_q_from_ab = ab / 1.25
        true_q = min(d, max_q_from_act, max_q_from_ab)

        if true_q <= 0:
            return

        # 正常反应，生成 Ed
        ed_amount = 4 * true_q
        act_used = 3.75 * true_q
        ab_used = 1.25 * true_q

        self.minerals['Ed'] += ed_amount
        self.minerals['Act'] -= act_used
        self.minerals['Ab'] -= ab_used
        self.negative_quartz -= true_q  # Q减少

        # 清除为0的矿物
        if self.minerals['Act'] <= 1e-6:
            del self.minerals['Act']
        if self.minerals['Ab'] <= 1e-6:
            del self.minerals['Ab']


    def replace_olivine_corundum(self):
        """步骤3：橄榄石 + 刚玉 → 尖晶石（Sp），按比例消耗"""

        ol = self.minerals.get('Ol', 0)
        cor = self.minerals.get('Cor', 0)
        d = self.negative_quartz

        if ol <= 0 or cor <= 0 or d <= 0:
            return

        # 可生成的石英量受限于资源和当前负 SiO₂
        q = min(d, ol / 3, cor / 4)

        if q <= 0:
            return

        # 生成 Sp
        self.minerals['Sp'] = self.minerals.get('Sp', 0) + 6 * q

        # 消耗原料
        self.minerals['Ol'] -= 3 * q
        self.minerals['Cor'] -= 4 * q
        self.negative_quartz -= q

        # 清理为 0 的矿物
        if self.minerals['Ol'] <= 1e-6:
            del self.minerals['Ol']
        if self.minerals['Cor'] <= 1e-6:
            del self.minerals['Cor']


        


def load_config(config_file="config.txt"):
    config_path = resource_path(config_file)
    config = {}
    with open(config_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.strip() and not line.strip().startswith("#"):
                key, value = line.strip().split("=")
                config[key.strip()] = value.strip()
    return config

def resource_path(relative_path):
    if getattr(sys, 'frozen', False):
        # PyInstaller 打包后运行的路径
        base_path = os.path.dirname(sys.executable)
    else:
        # 脚本运行时的路径
        base_path = os.path.dirname(os.path.abspath(__file__))
    full_path = os.path.join(base_path, relative_path)
    print(f"Resource path: {full_path}")  # 输出资源路径，用于调试
    return full_path

config = load_config(resource_path("config.txt"))

def process_sample(row_dict, calculator_class, kwargs):
    calc = calculator_class(**kwargs)
    calc.oxide_wt_percent = defaultdict(float)

    for oxide, value in row_dict.items():
        if oxide.lower() == "sample":
            continue
        if pd.notna(value):
            calc.oxide_wt_percent[oxide] = value
    calc.calculate_all_minerals()

    preferred_order = ["Q", "Or", "Ab", "An", "Bi", "Act","Ed","Ri", "En","Fs" "Wo", "Di", "Mt"]
    result = {}

    # 添加矿物数据，保留三位小数
    for key in preferred_order:
        if key in calc.minerals:
            value = calc.minerals[key]
            print(f"Key: {key}, Value: {value}, Type: {type(value)}")  # 添加调试输出
            if isinstance(value, (int, float)):  # 如果是数字类型
                result[key] = round(value, 3)
            else:
                result[key] = value  # 对非数字类型值进行处理（例如字典），这里直接保留

    for key, value in calc.minerals.items():
        if key not in preferred_order:
            print(f"Key: {key}, Value: {value}, Type: {type(value)}")  # 添加调试输出
            if isinstance(value, (int, float)):  # 如果是数字类型
                result[key] = round(value, 3)
            else:
                result[key] = value  # 对非数字类型值进行处理（例如字典），这里直接保留

    return result




def process_csv(input_file, output_file, calculator_class, **kwargs):
    df = pd.read_csv(resource_path(input_file))
    print(f"读取了 {len(df)} 条样本数据")

    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count()) as executor:
        futures = [executor.submit(process_sample, row.to_dict(), calculator_class, kwargs) for _, row in df.iterrows()]
        mineral_results = [future.result() for future in futures]
        mineral_df = pd.DataFrame(mineral_results)

    # 将矿物数据拼接到氧化物数据右侧
    final_df = pd.concat([df.reset_index(drop=True), mineral_df.reset_index(drop=True)], axis=1)

    final_df.to_csv(resource_path(output_file), index=False)



if __name__ == "__main__":
    config = load_config()

    input_file = config["input_file"]
    output_file = config["output_file"]
    method = config.get("method", "1")
    path = config.get("path", "1")

    if method == "1":
        process_csv(input_file, output_file, Catanormcalculator)
    elif method == "2":
        is_mafic = True if path == "1" else False
        process_csv(input_file, output_file, MesonormCalculator, is_mafic=is_mafic)
    else:
        print("无效的方法，请检查 config.txt")
