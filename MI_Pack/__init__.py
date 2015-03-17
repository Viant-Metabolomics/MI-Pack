#!/usr/local/bin/python

from SqliteManager import *
import os, sys, time
import copy
from copy import deepcopy
from collections import OrderedDict
from time import strftime

import math
import operator

from MI_Pack import EFScpp

try:
    import pp
except ImportError:
    _has_pp = False
else:
    _has_pp = True

try:
    import MySQLdb
except ImportError:
    _has_MySQLdb = False
else:
    _has_MySQLdb = True


class Atoms:
    def __init__(self):

        self.lib = OrderedDict([("C", 12.000000),
            ("(13C)", 13.003355),
            ("H", 1.007825),
            ("N", 14.003074),
            ("O", 15.994915),
            ("P", 30.973763),
            ("S", 31.972072),
            ("(34S)", 33.967868),
            ("[NaCl]", 22.989770+34.968853),
            ("[Na(37Cl)]", 22.989770+36.965903)])


        #C	4	12.000  50
        #H	1	1.0078250319	100
        #O	2	15.9949146223	10
        #N	3	14.0030740074	5
        #P	5	30.973763	5
        #S	2	31.972071	5

        #, "F":18.998403, "B":11.009305, "K":38.963708, "Na":22.989770,
        #"Cl":34.968853, "Si":27.976928, "Se":79.916521, "Mg":23.985045, "Br":78.9183,
        #"I":126.904477, "Ca":39.9625, "Ba":137.905236, "Fe": 55.934939, "R":0.000000}


        #self.lib = OrderedDict([("C", 12.000000),
        #    ("H", 1.0078250319),
        #    ("N", 14.0030740074),
        #    ("O", 15.9949146223),
        #    ("P", 30.973763),
        #    ("S", 31.972071)])
    
        self.limits = {
            "C":range(0, 35),
            "(13C)":range(0,2),
            "H":range(0, 73),
            "N":range(0, 16),
            "O":range(0, 20),
            "P":range(0, 8),
            "S":range(0, 9),
            "(34S)":range(0,2),
            "[NaCl]": range(0,2),
            "[Na(37Cl)]": range(0,2)}

        #self.limits = {
        #    "C":range(0, 101),
        #    "H":range(0, 201),
        #    "N":range(0, 11),
        #    "O":range(0, 51),
        #    "P":range(0, 11),
        #    "S":range(0, 11)}
        
        self.atom_bonds = {'C':4, 'H':1, 'N':3, 'O':2, 'P':3, 'S':2}

    def calculate_limits(self, mass):
        for atom in self.limits:
            min_atoms = min([int(self.limits[atom][len(self.limits[atom])-1]), int(mass/self.lib[atom]) + 1])
            if  min_atoms <= 0:
                self.limits[atom] = [0]
            else:
                self.limits[atom] = range(self.limits[atom][0], min_atoms + 1)

        return self.limits
        
    def alter_mass(self, name, mass):
        if type(mass) != float:
            print "mass: wrong format"
        else:
            if name in self.lib.keys():
                self.lib[name] = mass
            else:
                print "Name ion not in library"

    def alter_name(self, name, name_new):
        if name in self.lib.keys():
            mass = self.lib[name]
            del self.lib[name]
            self.lib[name_new] = mass
        else:
            print "Name ion not in library"

    def alter_limits(self, name, min_max):
        if name in self.limits.keys():
            self.limits[name] = range(min_max[0], min_max[1]+1)
        else:
            print "Name ion not in library"

    def add(self, name, mass):
        self.lib[name] = mass
        
    def add_limit(self, name, min, max):
        self.limits[name] = range(min, max+1)

    def delete(self, name):
        if name == "*":
            for name in self.lib.keys():
                del self.lib[name]
                try:
                    del self.limits[name]
                except:
                    pass
        else:
            del self.lib[name]     
            try:
                del self.limits[name]
            except:
                pass

    def __str__(self):
        print "Atoms/Molecules in library:"
        for key in self.lib:
            print "%s\t%s" % (key, self.lib[key])
        print
        print "Atom limits:"
        for key in self.limits:
            print "%s\t%s-%s" % (key, min(self.limits[key]), max(self.limits[key]))
        print
        return ""

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_atoms_length = self.lib.keys()
        self.sort_atoms_length.sort(cmp=bylength)
        return self.sort_atoms_length

class Ions:    
    def __init__(self, ion_mode):
        self.ion_mode = ion_mode.upper()
        self.e = 0.0005486
        
        if self.ion_mode == "POS":
            self.lib = {
                "[M+H]+":1.007825-self.e, 
                "[M+Na]+":22.989770-self.e, 
                "[M+K]+":38.963708-self.e,
                "[M+(41K)]+":40.961825-self.e}#,
                #"[M+(6Li)]+":6.015123-self.e,
                #"[M+Li]+":7.016005-self.e}
        
        elif self.ion_mode == "NEG":
            self.lib = {
                "[M-H]-": -(1.007825-self.e),
                "[M+Cl]-": 34.968853+self.e,
                "[M+(37Cl)]-": 36.965903+self.e,
                "[M+Na-2H]-":(22.989770-(2*1.007825))+self.e,
                "[M+K-2H]-":(38.963708-(2*1.007825))+self.e,
                "[M+Hac-H]-": 59.0138536}
                #"[M-2H]-": -2*(1.007825-self.e),
            
    def alter_mass(self, name, mass):
        if type(mass) != float:
            print "mass: wrong format"
        else:
            if name in self.lib.keys():
                self.lib[name] = mass
            else:
                print "Name ion not in library"

    def alter_name(self, name, name_new):
        if name in self.lib.keys():
            mass = self.lib[name]
            del self.lib[name]
            self.lib[name_new] = mass
        else:
            print "Name ion not in library"

    def add(self, name, mass):
        self.lib[name] = mass
            
    def delete(self, name):
        if name == "*":
            for name in self.lib.keys():
                del self.lib[name]
        else:
            del self.lib[name]
    
    def __str__(self):
        print "Ions in library:"
        for key in self.lib:
            print "%s\t%s" % (key, self.lib[key])
        print
        return ""

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_ions_length = self.lib.keys()
        self.sort_ions_length.sort(cmp=bylength)
        return self.sort_ions_length

class Isotopes:
    def __init__(self, ion_mode):
        self.ion_mode = ion_mode.upper()
        
        self.lib = [["C", "(13C)", 1.003355, 100.0, 1.1],
            ["S", "(34S)", 1.995796, 100.0, 4.21]] 
            
        if self.ion_mode == "POS":
            self.lib.append(["K", "(41K)", 1.998117, 100.0, 6.73])
            self.lib.append(["(6Li)", "Li", 1.000882, 7.42, 100])
            
        elif self.ion_mode == "NEG":
            self.lib.append(["Cl", "(37Cl)", 1.997050, 100.0, 24.23])    

    def add(self, record):
        if len(record) == 5 and record not in self.lib and type(record[0]) == str and type(record[1]) == str and type(record[2]) == float and type(record[3]) == float and type(record[4]) == float:
            self.lib.append(record)
        else:
            print "Check format record to add"

    def delete(self, record):
        if record == "*":
            self.lib = []
        else:
            if record in self.lib:
                self.lib.remove(record)
            else:
                print "Entry not in library"

    def __str__(self):
        print "Isotopes in library:"
        for entry in self.lib:
            print "\t".join(map(str, entry))
        print
        return ""

class Formula:
    
    def __init__(self, formula, atoms, ions, isotopes = []):
        
        self.atoms = atoms
        self.ions = ions
        self.isotopes = isotopes
        self.formula = formula
        
        self.count_atoms = {}
        self.count_ions = {}

        self.exact_mass = 0.0
        self.mass = 0.0

        self.order_atoms = []

        self.formula_print_it = ""

    def count(self, debug=0):

        import re
        
        self.temp_formula = self.formula
        
        for atom in self.atoms.lib:
            self.count_atoms[atom] = 0
        for ion in self.ions.lib:
            self.count_ions[ion] = 0

        l_ions = self.ions.sort_length()
        l_ions.reverse()

        for ion in l_ions:
            if ion in self.temp_formula:
                self.temp_formula =  self.temp_formula.replace(ion,"")
                self.count_ions[ion] += 1

        l_atoms = self.atoms.sort_length()
        l_atoms.reverse()

        for atom in l_atoms:
            str_pattern = []
            if "(" in atom or "[" in atom:

                num = len(atom)
                atom_search = atom.replace("(","\D").replace(")","\D").replace("[","\D").replace("]","\D")
                
                str_pattern.append("%s\d*" % (atom_search))                
                pattern_char = re.compile(r'(%s)' % ("|".join(map(str,str_pattern))))
                result = pattern_char.findall(self.temp_formula)
                for item in result:
                    if atom in item:
                        self.temp_formula = self.temp_formula.replace(item, "")
                        count = item.replace(atom, "")
                        if len(count) == 0:
                            self.count_atoms[atom] += 1
                        else:
                            self.count_atoms[atom] += int(count)

            else:
                num = len(atom)
                str_pattern.append("%s\d*" % (atom))                
                pattern_char = re.compile(r'(%s)' % ("|".join(map(str,str_pattern))))
                result = pattern_char.findall(self.temp_formula)
                for item in result:
                    if atom in item:
                        self.temp_formula = self.temp_formula.replace(item, "")
                        count = item.replace(atom, "")
                        if len(count) == 0:
                            self.count_atoms[atom] += 1
                        else:
                            self.count_atoms[atom] += int(count)
        if len(self.temp_formula) > 0:
            if debug == 1:
                print "Wrong Format Before:", self.formula, "After:", self.temp_formula
            for atom in self.atoms.lib:
                self.count_atoms[atom] = 0
            for ion in self.ions.lib:
                self.count_ions[ion] = 0
            return False

        return self.count_atoms, self.count_ions

    def parent(self, isotopes = []):

        if isotopes != []:
            self.isotopes = isotopes
        
        self.count()
        
        if self.isotopes != []:
            for ip in self.isotopes.lib:
                if ip[1] in self.count_atoms:
                    self.count_atoms[ip[0]] += self.count_atoms[ip[1]]
                    self.count_atoms[ip[1]] = 0

        set_to_zero = []
        total_atoms = 0
        for a in self.count_atoms:
            if "[" and "]" in a and self.count_atoms[a] != 0:
                set_to_zero.append(a)
            else:
                total_atoms += self.count_atoms[a]
            if total_atoms > 0 and len(set_to_zero) > 0:
                for a in set_to_zero:
                    self.count_atoms[a] = 0

        return self.print_it(self.count_atoms)
        
    def calc_mass(self):

        self.count()
        
        self.exact_mass = 0.0
        self.mass = 0.0
        
        for atom in self.count_atoms:
            if self.count_atoms[atom] != 0:
                self.exact_mass += (self.count_atoms[atom] * (self.atoms.lib[atom]))
                
        self.mass += self.exact_mass
        for ion in self.count_ions:
            if self.count_ions[ion] != 0:
                self.mass += (self.count_ions[ion] * (self.ions.lib[ion]))

        return self.exact_mass, self.mass
    
    def HNOPS(self):
        
        self.parent() # update count_atoms 
        
        self.HNOPS_rule = {"HC":0, "NOPSC":0}
        
        if (self.count_atoms['C']) == 0 or self.count_atoms['H'] == 0:
            self.HNOPS_rule["HC"]= 0
        elif (self.count_atoms['C']) == 0 and self.count_atoms['H'] == 0:
            self.HNOPS_rule["HC"]= 0
        elif (self.count_atoms['C']) != 0 and self.count_atoms['H'] != 0:
            if float(self.count_atoms['H'])/float((self.count_atoms['C'])) > 0 and float(self.count_atoms['H']/(self.count_atoms['C'])) < 6:
                self.HNOPS_rule["HC"]= 1
            if float(self.count_atoms['H'])/float((self.count_atoms['C'])) >= 6:
                self.HNOPS_rule["HC"]= 0

        #NOPS
        NOPS_check = []
        for atom in ['N','O','P','S']:
            try:
                check = float(float(self.count_atoms[atom]))/float((self.count_atoms['C']))
                NOPS_check.append(check)
            except ZeroDivisionError:
                NOPS_check.append(float(0))

        if NOPS_check[0] >= float(0) and NOPS_check[0] <= float(4) and NOPS_check[1] >= float(0) and NOPS_check[1] <= float(3) and NOPS_check[2] >= float(0) and NOPS_check[2] <= float(2) and NOPS_check[3] >= float(0) and NOPS_check[3] <= float(3):
            self.HNOPS_rule["NOPSC"] = 1
        if NOPS_check[0] > float(4) or NOPS_check[1] > float(3) or NOPS_check[2] > float(2) or NOPS_check[3] > float(3):
            self.HNOPS_rule["NOPSC"] = 0
        return self.HNOPS_rule
        
    def lewis_senior(self):

        self.parent() # update count_atoms 
        
        self.lewis_senior_rule = {"lewis":0, "senior":0}
        #atoms_sum = sum(self.count_atoms.values())

        atoms_sum = sum(self.count_atoms.values())

        lewis_sum = 0

        for atom in self.atoms.atom_bonds:
            if atom == 'S':
                lewis_sum += self.atoms.atom_bonds[atom] * (self.count_atoms[atom])

            elif atom == 'C':
                lewis_sum += self.atoms.atom_bonds[atom] * (self.count_atoms[atom])
            
            else:
                lewis_sum += self.atoms.atom_bonds[atom] * self.count_atoms[atom]
        
        if lewis_sum%2 == 0:
            self.lewis_senior_rule["lewis"] = 1
        if lewis_sum%2 != 0:
            self.lewis_senior_rule["lewis"] = 0
        if lewis_sum >= ((atoms_sum-1)*2):
            self.lewis_senior_rule["senior"] = 1
        if lewis_sum < ((atoms_sum-1)*2):
            self.lewis_senior_rule["senior"] = 0
        self.count()
        atoms_sum = sum(self.count_atoms.values())
        return self.lewis_senior_rule

    def print_it(self, atom_count={}, ion_count={}):

        if len(atom_count) > 0:
            self.count_atoms = atom_count
        if len(ion_count) > 0:
            self.count_ions = ion_count
        if atom_count == {} and ion_count == {}:
            self.count()
        
        self.formula_print_it = ""
        self.order_atoms = ["C", "(13C)", "H", "N", "O", "P", "S", "(34S)"]
        for atom in self.atoms.sort_length():
            if atom not in self.order_atoms:
                self.order_atoms.append(atom)

        for atom in self.order_atoms:
            try:
                if abs(self.count_atoms[atom]) > 1:
                    self.formula_print_it += atom + str(self.count_atoms[atom])
                if abs(self.count_atoms[atom]) == 1:
                    self.formula_print_it += atom
            except:
                pass

        for ion in self.ions.sort_length():
            try:
                if self.count_ions[ion] > 1:
                    self.formula_print_it += ion + str(self.count_ions[ion])
                if self.count_ions[ion] == 1:
                    self.formula_print_it += ion
            except:
                pass
        if self.formula_print_it == "":
            return False
        else:
            return self.formula_print_it

    def diff(self, inp_formula):

        self.count()

        self.inp_formula = Formula(inp_formula, self.atoms, self.ions, self.isotopes)
        self.inp_formula.count()
        
        self.formulae_diff = ""
        
        for item in self.count_atoms.keys():
             self.inp_formula.count_atoms[item] -= self.count_atoms[item]        

        for atom in self.atoms.sort_length():
            if self.inp_formula.count_atoms[atom] == 1:
                self.formulae_diff += "-" + str(atom)

            if self.inp_formula.count_atoms[atom] == -1:
                self.formulae_diff += str(atom)
                
            if self.inp_formula.count_atoms[atom] < -1:
                self.formulae_diff += str(atom) + str(abs(self.inp_formula.count_atoms[atom]))

            if self.inp_formula.count_atoms[atom] > 1:
                self.formulae_diff += "-" + str(atom) + str(self.inp_formula.count_atoms[atom])

        return self.formulae_diff



class Calculate:
    
    def __init__(self):
        self.mass = 0.0
        self.ppm = 1.0
    
    def tolerance_peak(self, mass, ppm):
        self.mass = mass
        self.ppm = ppm
        self.min_tol = mass - (mass * 0.000001 * ppm)
        self.max_tol = mass + (mass * 0.000001 * ppm)
        return self.min_tol, self.max_tol

    def abs_error(self, mass, ppm):
        self.mass = mass
        self.ppm = ppm
        self.error = (mass * 0.000001 * ppm)
        return self.error
        
    def ppm_error(self, mass, theo_mass):
        return float(theo_mass-mass) / (theo_mass * 0.000001)

    def load_surface(self, fname_surface):

        try:
            import rpy
        except ImportError:
            print "Please, install R (http://www.r-project.org/) and (Rpy: http://rpy.sourceforge.net/rpy.html)"
            return

        try:
            rpy.r("library(GenKern)")
        except:
            rpy.r.install_packages("GenKern")
            rpy.r("library(GenKern)")
        rpy.r('load(file = "%s")' % (fname_surface))

    def mass_error_mz_pair_diff(self, avr_mrange, mz_pair_diff, name_surface):
        
        try:
            import rpy
        except ImportError:
            print "Please, install R (http://www.r-project.org/) and (Rpy: http://rpy.sourceforge.net/rpy.html)"
            return
        
        set_default_mode(BASIC_CONVERSION)
        rpy.r("mrange = %s" % (avr_mrange))
        rpy.r("mz_pair_diff = %s" % (mz_pair_diff))
        x = rpy.r("index_mrange = nearest(%s$x, mrange, outside=TRUE)" % (name_surface))
        y = rpy.r("index_mz_pair_diff = nearest(%s$y, mz_pair_diff, outside=TRUE)" % (name_surface))
        out = rpy.r("%s$z[index_mrange, index_mz_pair_diff]" % (name_surface))
        set_default_mode(NO_CONVERSION)
        return out     

    def mean(self, values):
        return sum(values) / len(values)

    def stdev(self, values):
        n = len(values)
        mean_value = self.mean(values)
        for value in values:
            std = std + (a - mean_value)**2
        std = math.sqrt(std / float(n-1))
        return round(mean_value,1), round(std,1)


def IntensityMatrix(file_matrix, file_mz = None):
    
    intensities = OrderedDict()

    if file_mz != None:
        file_mz = file_mz
        f = open(file_mz)
        for line in f.readlines():
            intensities[line.replace("\n","")] = [] # SIM-stitch output

    f = open(file_matrix)
    if file_mz == None:
        header = f.readline()        
    lines = f.readlines()
    
    if file_mz != None: #read SIM-stitch output
        for line in lines:
            line = line.replace("\n","").split("\t")
            i = 0
            for value in line:
                if value != "":
                    intensities[intensities.keys()[i]].append(float(value))
                    i += 1
    else:
        for line in lines:
            temp = []
            line = line.replace("\n","").split("\t")
            for value in line[1:]:
                temp.append(float(value))
            intensities[line[0]] = temp
    return intensities


class PeakList:
    
    def __init__(self, file_peaklist):

        self.file_peaklist = file_peaklist
        self.intensities = OrderedDict()
        self.mzs = OrderedDict()
        
        f = open(self.file_peaklist)
        self.header = f.readline()

        for line in f.readlines():
            line = line.replace("\n","").split("\t")
            self.intensities[float(line[0])] = round(float(line[1]),2)
            self.mzs[round(float(line[1]),2)] = float(line[0])

    def list_mzs(self):
        return self.intensities.keys()

    def list_intensities(self):
        return self.mzs.keys()

    def intensity(self, mz_inp):
        return self.intensities[mz_inp]

    def mz(self, intensity_inp):
        return self.mzs[intensity_inp]

    """
    def mz_pair_diffs(self, n = 2, peaks=[]):
        if len(peaks) == 0: 
            peaks = self.intensities.keys()
        if n == 0:
            yield []
        else:
            for i in xrange(len(peaks)):
                for cc in self.mz_pair_diffs(n-1, peaks[i+0:]):
                    if len([peaks[i]]+cc) == 2:

                        if peaks[i] < cc[0]:
                            self.intensity_a = self.intensity(peaks[i])
                            self.intensity_b = self.intensity(cc[0])
                            self.mean_mz_ab = (peaks[i] + cc[0])/2

                            yield {"mz_a":peaks[i], "intensity_a":self.intensity_a, "mz_b":cc[0], 
                               "intensity_b":self.intensity_b, "mz_pair_diff":cc[0] - peaks[i], "mean_mz_ab":self.mean_mz_ab}
                            
                        elif peaks[i] > cc[0] or peaks[i] == cc[0]:
                            self.intensity_a = self.intensity(peaks[i])
                            self.intensity_b = self.intensity(cc[0])
                            self.mean_mz_ab = (peaks[i] + cc[0])/2
                            
                            yield {"mz_a":peaks[i], "intensity_a":self.intensity_a, "mz_b":cc[0], 
                               "intensity_b":self.intensity_b, "mz_pair_diff":peaks[i] - cc[0], "mean_mz_ab":self.mean_mz_ab}

                    else:
                        yield [peaks[i]]+cc

    """
    
    def mz_pair_diffs(self, ions = {}, peaks=[]):
        if len(peaks) == 0: 
            peaks = self.intensities.keys()

        for i in xrange(len(peaks)):
            for j in xrange(len(peaks)):
                if peaks[i] <= peaks[j]:

                    if len(ions.keys()) > 0:
                        
                        for k in ions:
                            for l in ions:
                                #print k, l
                                #print peaks[i], peaks[i]-ions[k]
                                #print peaks[j], peaks[j]-ions[l]
                                #print 
                                yield {"mz_a":peaks[i], "intensity_a":self.intensity(peaks[i]), "ion_a":k, "mz_b":peaks[j], 
                                "intensity_b":self.intensity(peaks[j]), "ion_b":l, "mz_pair_diff_neutral":(peaks[i]-ions[k]) - (peaks[j]-ions[l]), "mz_pair_diff":(peaks[i]-peaks[j]), "mean_mz_ab":(peaks[i] + peaks[j])/2}
                                #raw_input()
                    else:
                        yield {"mz_a":peaks[i], "intensity_a":self.intensity(peaks[i]), "ion_a":None, "mz_b":peaks[j], 
                        "intensity_b":self.intensity(peaks[j]), "ion_b":None, "mz_pair_diff_neutral":(peaks[i]-peaks[j]), "mz_pair_diff":(peaks[i]-peaks[j]), "mean_mz_ab":(peaks[i] + peaks[j])/2}
                        #raw_input()
                        
class Identification:

    def __init__(self, fn_peak_list, project_name, ion_mode, ppm, atoms, ions, isotopes, fn_MIDB = None):

        self.fn_MIDB = fn_MIDB
        
        if os.path.exists(fn_peak_list) == 1 and os.path.isfile(fn_peak_list) == 1:
            self.fn_peak_list = fn_peak_list
            self.name_peak_list = os.path.basename(fn_peak_list).split(".")[0]
            self.path = os.path.dirname(self.fn_peak_list)
        else:
            print "Filename %s does not exist" % (fn_peak_list)
            sys.exit()
           
        self.peak_list = PeakList(self.fn_peak_list)
        self.project_name = project_name

        self.db_file_name = os.path.join("%s" % (project_name))
        
        self.ion_mode = ion_mode.upper()
        self.ppm = ppm
        
        self.ions = ions
        self.atoms = atoms
        self.isotopes = isotopes
        
        self.calc = Calculate()
        self.num_peaks = len(self.peak_list.intensities)

        self.transformation_type = 0

        self.subset_MIDB_select = OrderedDict()
        self.subset_MIDB_exclude = OrderedDict()
        
        self.subset_MIDB_KEGG_pathways = []
        self.subset_MIDB_KEGG_organisms = []
        
        self.subset_MIDB_table = ""
        
        self.total_peaks = 0
        self.peaks_done = 0
        
        self.database = SQLiteManager(self.db_file_name)
        if fn_MIDB != None and fn_MIDB != "":
            self.database.attach(self.ion_mode, self.db_file_name, fn_MIDB)
        
        self.output_list = [] #output
        #self._summary("", add=False)

    def MIDB(self, type_subset, inp_subset):
        
        self.types_subset = ["MIDB", "select", "exclude"]

        if type_subset not in self.types_subset:
            print "Type subset does not exist. Please, use one of the following types: 'MIDB', 'select' or 'exclude'"
            sys.exit()

        if type(inp_subset) != dict:
            print "Dictionary required"
            sys.exit()
        
        elif type(inp_subset) == dict:
            self.inp_subset = inp_subset
            dbs_ordered = inp_subset.keys()
            dbs_ordered.sort()            
            
        #print type_subset
        if type_subset == "MIDB":
            
            for db in dbs_ordered:
                
                temp_cpds = []
                
                if inp_subset[db] != "*":
                    cpds = self.database.select("class_compound", ["DISTINCT id"], "WHERE class in ('%s')" % (", ".join(map(str,inp_subset[db]))))

                    if db == "KEGG_COMPOUND":
                        
                        self.subset_MIDB_KEGG_organisms = inp_subset[db]

                    for cpd in cpds:
                        temp_cpds.append(str(cpd[0]))
                    self.subset_MIDB_select[db] = temp_cpds
                    
                elif inp_subset[db] == "*":
                    self.subset_MIDB_select[db] = "*"

        if type_subset == "select" or type_subset == "exclude":  
            for db in dbs_ordered:
                temp_cpds = []
                if len(inp_subset[db]) != 0:
                    for id in inp_subset[db]:
                        if db == "KEGG_COMPOUND":
                            if len(id) == 3:
                                classes = self.database.select("class_compound", ["distinct class"])
                                for class_name in classes:
                                    if id == class_name[0]:
                                        cpds = self.database.select("class_compound", ["distinct id"], "where class = '%s'" % (id))
                                        for cpd in cpds:
                                            if str(cpd[0]) not in temp_cpds:
                                                temp_cpds.append(str(cpd[0]))
                                        if id not in self.subset_MIDB_KEGG_organisms and type_subset == "select": # map()
                                            self.subset_MIDB_KEGG_organisms.append(id) # map()
                                          
                            elif len(id) == 7 or len(id) == 8:
                                KEGG_maps = self.database.select("map_compound", ["distinct mapid"])
                                for record in KEGG_maps:
                                    
                                    if id == str(record["mapid"]):  
                                        cpds = self.database.select("map_compound", ["distinct KEGG_COMPOUND"], "where mapid = '%s'" % (id))
                                        for cpd in cpds:
                                            if str(cpd["KEGG_COMPOUND"]) not in temp_cpds:
                                                temp_cpds.append(str(cpd["KEGG_COMPOUND"]))
                                            if str(record["mapid"]) not in self.subset_MIDB_KEGG_pathways  and type_subset == "select":
                                                self.subset_MIDB_KEGG_pathways.append(str(record["mapid"]))
                            else:
                                if "C" in id and id not in temp_cpds and len(id) == 6:
                                    temp_cpds.append(id)
                                    
                        else: # other databases
                            temp_cpds.append(id)

                    if type_subset == "select":
                        self.subset_MIDB_select[db] = temp_cpds                    
                    elif type_subset == "exclude":
                        self.subset_MIDB_exclude[db] = temp_cpds
        
        
    def SPS(self):
        
        print
        print "Single-peak search - %s" % (self.fn_peak_list)
        print "------------------------------------------------------------------"
        
        self.total_peaks = len(self.peak_list.intensities)
        self.peaks_done, perc = 0, 0.0

        #########################################################################
        # SQL TEMPLATE SEARCH
        #########################################################################

        if len(self.subset_MIDB_select.keys()) == 0:
            print "No database selected"
            sys.exit()
        else:   
            for db in self.subset_MIDB_select:
                for db2 in self.subset_MIDB_exclude:
                    if db == db2:
                        temp_ids = deepcopy(self.subset_MIDB_exclude[db])
                        for id in temp_ids:
                            if id in self.subset_MIDB_select[db]:
                                del self.subset_MIDB_select[db][self.subset_MIDB_select[db].index(id)]
                                del self.subset_MIDB_exclude[db][self.subset_MIDB_exclude[db].index(id)]            

        sql_str = ""
        sql_str_dbs = ""
        
        
        for db in self.subset_MIDB_select:
            sql_str_dbs += "%s text DEFAULT NULL," % (db)
            if len(self.subset_MIDB_select[db]) != 0 and type(self.subset_MIDB_select[db]) == list:
                sql_str += " AND %s in ('%s')\n" % (db, "', '".join(map(str,self.subset_MIDB_select[db])))
        
        for db in self.subset_MIDB_exclude:
            if len(self.subset_MIDB_exclude[db]) != 0 and type(self.subset_MIDB_exclude[db]) == list:
                sql_str += " AND %s not in ('%s')\n" % (db, "', '".join(map(str,self.subset_MIDB_exclude[db])))

        sql_str +=  " AND (%s IS NOT NULL)" % (" IS NOT NULL OR ".join(map(str,self.subset_MIDB_select.keys())))
        
        # check ions
        default_ions = Ions(self.ion_mode)
        if self.ions != None and type(self.ions) == type(default_ions):
            if self.ions.lib != default_ions.lib:
                pass
                #sql_str += " AND ion in ('%s')" % ("', '".join(map(str,self.ions.lib.keys())))
        # check ions

        self.database.execute("DROP TABLE IF EXISTS SPS_%s_%s;" % (self.ion_mode, self.name_peak_list))
        self.database.execute("""
            CREATE TABLE SPS_%s_%s (
            mz decimal(12,7) DEFAULT NULL,
            intensity decimal(18,8) DEFAULT NULL,
            formula text DEFAULT NULL,
            ion char(20) DEFAULT NULL,
            name text DEFAULT NULL,
            exact_mass decimal(12,7) DEFAULT NULL,
            mass decimal(12,7) DEFAULT NULL,
            ppm_error decimal(12,2) DEFAULT NULL,
            %s
            primary key  (%s, ion, mz, intensity, name)
            );""" % (self.ion_mode, self.name_peak_list, sql_str_dbs, ", ".join(map(str,self.subset_MIDB_select))))

        for mz in copy.deepcopy(self.peak_list.intensities):
            for ion in self.ions.lib:
                #print ion
                self.calc.tolerance_peak(mz, self.ppm)
                columns_select = ["ion", "formula", "name", "mz", "intensity", "mass", "ppm_error", "exact_mass"]

                values_to_insert = ['"%s"' % (ion), "formula", "name", str(mz), str(self.peak_list.intensities[mz]), "exact_mass+%s" % (self.ions.lib[ion]), "(%s-(exact_mass+%s))/((exact_mass+%s)*0.000001)" % (mz, self.ions.lib[ion], self.ions.lib[ion]), "exact_mass"]

                values_to_insert.extend(self.subset_MIDB_select.keys()) # ADD DB NAMES
                columns_select.extend(self.subset_MIDB_select.keys()) # ADD DB NAMES

                #print "WHERE exact_mass >= %s AND exact_mass <= %s%s" % (self.calc.min_tol-self.ions.lib[ion], self.calc.max_tol-self.ions.lib[ion], sql_str)
                #raw_input()

                self.database.insert_select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list),
                columns_select, values_to_insert,
                "compounds",
                "WHERE exact_mass >= %s AND exact_mass <= %s%s" % (self.calc.min_tol-self.ions.lib[ion], self.calc.max_tol-self.ions.lib[ion], sql_str))

            self.peaks_done += 1
            if float(self.peaks_done)/self.total_peaks >= perc:
                str_perc = str(perc * 100.0) + "%"
                perc += 0.1
                print str_perc,
        print
        print "------------------------------------------------------------------"
        print "Finished"
        print
        self.database.save()
        self.database.index("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), "SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["mz"])
        self.database.save()


    def TM(self, transformation_type = 1, transformations = False, fname_surface=""):

        self.transformation_type = transformation_type
        
        if len(fname_surface) > 0 and type(fname_surface) == str:
            Calculate.load_surface(fname_surface)
        
        sql_str = "" # Testing
        
        self.subset_MIDB_table = ""
        if "KEGG_COMPOUND" in self.subset_MIDB_select:
            if type(self.subset_MIDB_select["KEGG_COMPOUND"]) == list and len(self.subset_MIDB_select["KEGG_COMPOUND"]) != 0:

                self.database.drop("SUBSET_rpairs")
                self.database.create_select("SUBSET_rpairs", ["*"], "rpairs",
                "where KEGG_COMPOUND_a in ('%s') and KEGG_COMPOUND_b in ('%s')%s" % ("', '".join(map(str,self.subset_MIDB_select["KEGG_COMPOUND"])), "', '".join(map(str,self.subset_MIDB_select["KEGG_COMPOUND"])), sql_str), 1) #1 = temporary table

                if transformations == True:
                    self.database.create_table("SUBSET_unqiue_transformations", self.ion_mode)
                    self.database.index("SUBSET_mass_pair_diff", "SUBSET_unique_transformations", ["mass_pair_diff", "type"])
                
                self.database.index("SUBSET_KEGG_COMPOUND_a_cr", "SUBSET_rpairs", ["KEGG_COMPOUND_a"])
                self.database.index("SUBSET_KEGG_COMPOUND_b_cr", "SUBSET_rpairs", ["KEGG_COMPOUND_b"])
                self.database.index("SUBSET_type_cr", "SUBSET_rpairs", ["type"])
                                
                self.subset_MIDB_table = "SUBSET_" ### IMPORTANT


        if transformations == True:

            print
            print "Peak-pair differences to transformations - %s" % (self.fn_peak_list)
            print "------------------------------------------------------------------"

            self.database.create_table("TM", self.ion_mode, self.name_peak_list)

            self.total_peaks = len(self.peak_list.intensities)
            niters = ((self.total_peaks * self.total_peaks * len(self.ions.lib) * len(self.ions.lib))/2) + (2 * len(self.ions.lib) * self.total_peaks)

            iter_done, perc = 0, 0.0

            columns_insert = ["mz_a", "mz_b", "intensity_a", "intensity_b", "ion_a", "ion_b",
                "transformation", "mz_pair_diff", "mz_pair_diff_neutral", "mass_pair_diff", "ppm_error", "type"]

            for item in self.peak_list.mz_pair_diffs(self.ions.lib):

                if item["mz_pair_diff_neutral"] == 0.0:
                    mass_error_mz_pair_diff = 0.0
       
                elif len(fname_surface) > 0 and type(fname_surface) == str:
                    mass_error_mz_pair_diff = self.calc.mass_error_mz_pair_diff(item["mean_mz_ab"], item["mz_pair_diff"], "surface")
                else:
                    mass_error_mz_pair_diff = self.calc.abs_error(item["mz_a"], self.ppm) + self.calc.abs_error(item["mz_b"], self.ppm)
                
                min_tol_mz_pair_diff = item["mz_pair_diff_neutral"] - mass_error_mz_pair_diff
                max_tol_mz_pair_diff = item["mz_pair_diff_neutral"] + mass_error_mz_pair_diff

                if item["mz_a"] <= item["mz_b"]:

                    columns_select_1 = [item["mz_a"], item["mz_b"], item["intensity_a"], item["intensity_b"],
                    "'%s'" % (item["ion_a"]), "'%s'" % (item["ion_b"]), "transformation", item["mz_b"] - item["mz_a"], item["mz_pair_diff_neutral"], "mass_pair_diff",
                    "(ABS(mass_pair_diff)-%s)/abs(mass_pair_diff*0.000001)" % (abs(item["mz_pair_diff_neutral"])), "type"]

                    where_1 = """WHERE
                    (mass_pair_diff >= %s AND mass_pair_diff <= %s)
                    AND type <= %s%s""" % (abs(min_tol_mz_pair_diff), abs(max_tol_mz_pair_diff), self.transformation_type, sql_str)

                    columns_select_2 = [item["mz_b"], item["mz_a"], item["intensity_b"], item["intensity_a"],
                    "'%s'" % (item["ion_b"]), "'%s'" % (item["ion_a"]), "transformation", item["mz_a"] - item["mz_b"], item["mz_pair_diff_neutral"], "mass_pair_diff",
                    "(ABS(mass_pair_diff)-%s)/abs(mass_pair_diff*0.000001)" % (abs(item["mz_pair_diff_neutral"])), "type"]

                    where_2 = """WHERE
                    (mass_pair_diff <= -%s AND mass_pair_diff >= -%s)
                    AND type <= %s%s""" % (abs(min_tol_mz_pair_diff), abs(max_tol_mz_pair_diff), self.transformation_type, sql_str)                 

                    self.database.insert_select("transformations_%s_%s" % (self.ion_mode, self.name_peak_list), columns_insert, columns_select_1, "%sunique_transformations" % (self.subset_MIDB_table), where_1)
                    self.database.insert_select("transformations_%s_%s" % (self.ion_mode, self.name_peak_list), columns_insert, columns_select_2, "%sunique_transformations" % (self.subset_MIDB_table), where_2)

                elif item["mz_a"] >= item["mz_b"]:
                
                    columns_select_1 = [item["mz_a"], item["mz_b"], item["intensity_a"], item["intensity_b"],
                    "'%s'" % (item["ion_a"]), "'%s'" % (item["ion_b"]), "transformation", item["mz_b"] - item["mz_a"], "mass_pair_diff",
                    "(ABS(mass_pair_diff)-%s)/abs(mass_pair_diff*0.000001)" % (abs(item["mz_pair_diff_neutral"])), "type"]

                    where_1 = """WHERE
                    (mass_pair_diff >= -%s AND mass_pair_diff <= -%s)
                    AND type <= %s%s""" % (abs(max_tol_mz_pair_diff), abs(min_tol_mz_pair_diff), self.transformation_type, sql_str)

                    columns_select_2 = [item["mz_b"], item["mz_a"], item["intensity_b"], item["intensity_a"],
                    "'%s'" % (item["ion_b"]), "'%s'" % (item["ion_a"]), "transformation", item["mz_a"] - item["mz_b"], "mass_pair_diff",
                    "(ABS(mass_pair_diff)-%s)/abs(mass_pair_diff*0.000001)" % (abs(item["mz_pair_diff_neutral"])), "type"]

                    where_2 = """WHERE
                    (mass_pair_diff >= %s AND mass_pair_diff <= %s)
                    AND type <= %s%s""" % (abs(min_tol_mz_pair_diff), abs(max_tol_mz_pair_diff), self.transformation_type, sql_str)
                    
                    self.database.insert_select("transformations_%s_%s" % (self.ion_mode, self.name_peak_list), columns_insert, columns_select_1, "%sunique_transformations_%s" % (self.subset_MIDB_table, self.ion_mode), where_1)
                    self.database.insert_select("transformations_%s_%s" % (self.ion_mode, self.name_peak_list), columns_insert, columns_select_2, "%sunique_transformations_%s" % (self.subset_MIDB_table, self.ion_mode), where_2)
                    self.database.save()
            
                iter_done += 1
                if float(iter_done)/niters >= perc:
                    str_perc = str(perc * 100.0) + "%"
                    perc += 0.1
                    print str_perc,

            self.database.index("transformations_%s_%s" % (self.ion_mode, self.name_peak_list), "transformations_%s_%s" % (self.ion_mode, self.name_peak_list), ["mz_a", "transformation", "type"])
            self.database.save()
            
            print
            print "------------------------------------------------------------------"
            print "Finished"
            print


        self.total_peaks = len(self.peak_list.intensities)
        self.peaks_done, perc = 0, 0.0
        print
        print "Assigning reactant pairs - %s" % (self.fn_peak_list)
        print "------------------------------------------------------------------"

        self.database.execute("DROP TABLE IF EXISTS SPS_TEMP_%s_%s;" % (self.ion_mode, self.name_peak_list))
        self.database.execute("""
            CREATE TEMPORARY TABLE SPS_TEMP_%s_%s (
            mz decimal(12,7) DEFAULT NULL,
            intensity decimal(18,8) DEFAULT NULL,
            formula text DEFAULT NULL,
            ion char(20) DEFAULT NULL,
            name text DEFAULT NULL,
            exact_mass decimal(12,7) DEFAULT NULL,
            mass decimal(12,7) DEFAULT NULL,
            ppm_error decimal(12,2) DEFAULT NULL,
            KEGG_COMPOUND text DEFAULT NULL,
            primary key  (KEGG_COMPOUND, ion, mz, intensity, name)
            );""" % (self.ion_mode, self.name_peak_list))
        
        for mz in copy.deepcopy(self.peak_list.intensities):
            for ion in self.ions.lib:
                
                self.calc.tolerance_peak(mz, self.ppm)

                where = "WHERE exact_mass >= %s AND exact_mass <= %s and KEGG_COMPOUND NOT NULL" % (self.calc.min_tol-self.ions.lib[ion], self.calc.max_tol-self.ions.lib[ion])
                columns_insert = ["mz", "intensity", "formula", "ion", "name", "exact_mass", "mass", "ppm_error", "KEGG_COMPOUND"]
                columns_select = [mz, self.peak_list.intensities[mz], "formula", '"%s"' % (ion), "name", "exact_mass", "exact_mass + %s" % (self.ions.lib[ion]), "round((%s-(exact_mass+%s))/((exact_mass+%s)*0.000001),2)" % (mz, self.ions.lib[ion], self.ions.lib[ion]), "KEGG_COMPOUND"]
                self.database.insert_select("SPS_TEMP_%s_%s" % (self.ion_mode, self.name_peak_list), columns_insert, columns_select, "compounds", where)

            self.peaks_done += 1
            if float(self.peaks_done)/self.total_peaks >= perc:
                str_perc = str(perc * 100.0) + "%"
                perc += 0.1
                print str_perc,

        self.database.execute("DROP TABLE IF EXISTS reactant_pairs_%s_%s;" % (self.ion_mode, self.name_peak_list))
        self.database.execute("""
        CREATE TABLE reactant_pairs_%s_%s AS SELECT
            SPS1.mz as mz_a, SPS2.mz as mz_b,
            SPS2.mz - SPS1.mz as mz_diff,
            SPS2.mass - SPS1.mass as mass_diff,
            (SPS2.mz - (SPS2.mass - SPS2.exact_mass)) - (SPS1.mz - (SPS1.mass - SPS1.exact_mass)) as mz_diff_neutral,
            rpairs.mass_pair_diff as mass_diff_netural,           
            SPS1.ion as ion_a, SPS2.ion as ion_b,
            rpairs.rpairid,
            rpairs.KEGG_compound_a, rpairs.KEGG_compound_b,
            rpairs.transformation,
            rpairs.type
        FROM 
            %srpairs as rpairs
            JOIN SPS_TEMP_%s_%s as SPS1 ON rpairs.KEGG_compound_a = SPS1.KEGG_COMPOUND
            JOIN SPS_TEMP_%s_%s as SPS2 ON rpairs.KEGG_compound_b = SPS2.KEGG_COMPOUND
            WHERE SPS1.KEGG_COMPOUND IS NOT SPS2.KEGG_COMPOUND;
        """ % (self.ion_mode, self.name_peak_list, self.subset_MIDB_table, self.ion_mode, self.name_peak_list, self.ion_mode, self.name_peak_list))
        self.database.save()
        
        self.database.execute("DROP TABLE IF EXISTS TM_%s_%s;" % (self.ion_mode, self.name_peak_list))
        self.database.execute("""CREATE TABLE TM_%s_%s AS SELECT
        * FROM SPS_TEMP_%s_%s
            WHERE KEGG_COMPOUND IN 
            (SELECT DISTINCT KEGG_compound_a from reactant_pairs_%s_%s);
        """ % (self.ion_mode, self.name_peak_list, self.ion_mode, self.name_peak_list, self.ion_mode, self.name_peak_list))
        self.database.save()
        
        print
        print "------------------------------------------------------------------"
        print "Finished"
        print


    def Maps(self, strategy):
        
        if strategy == "SPS":
            type = 0
            
        elif strategy == "TM":
            type = self.transformation_type
        
        try:
            self.database.create_table("maps_coverage") 
        except:
            pass

        #print strategy
        #print self.subset_MIDB_KEGG_organisms
        
        if len(self.subset_MIDB_KEGG_organisms) > 0 and self.subset_MIDB_KEGG_organisms != ["*"] and len(self.subset_MIDB_KEGG_pathways) == 0:

            for org in self.subset_MIDB_KEGG_organisms:
                records = self.database.select("map_compound", ["DISTINCT mapid"], "where map_compound.mapid like '""" + str(org) + "%'")
                for record in records:
                    if str(record[0]) not in self.subset_MIDB_KEGG_pathways:
                        self.subset_MIDB_KEGG_pathways.append(str(record[0]))
        
        elif self.subset_MIDB_KEGG_organisms == ["*"] and len(self.subset_MIDB_KEGG_pathways) == 0:

            records = self.database.select("map_compound", ["DISTINCT mapid"], "where map_compound.mapid like 'map%' OR (map_compound.mapid like 'ko%' and length(map_compound.mapid) = 7)""")
            for record in records:
                if str(record[0]) not in self.subset_MIDB_KEGG_pathways:
                    self.subset_MIDB_KEGG_pathways.append(str(record[0]))
        else:
            # Pathway selection KEGG already exists!
            pass           
        
        sql_str_compound = """ AND map_compound.mapid in ('%s')""" % ("','".join(map(str,self.subset_MIDB_KEGG_pathways)))
        sql_str_formula = """ AND map_formula.mapid in ('%s')""" % ("','".join(map(str,self.subset_MIDB_KEGG_pathways)))
        
        self.database.create_table("maps_%s"  % (strategy), self.ion_mode, self.name_peak_list)
        self.database.execute_complex("maps_%s"  % (strategy), (sql_str_compound, sql_str_formula), self.name_peak_list, self.ion_mode)
        
        maps_cpds = self.database.select("map_cpd_%s_%s_%s" % (strategy, self.ion_mode, self.name_peak_list), ["DISTINCT mapid"])

        for m in maps_cpds:

            str_compounds = ""
            
            compoundids = self.database.select("map_cpd_%s_%s_%s" % (strategy, self.ion_mode, self.name_peak_list), ["DISTINCT KEGG_COMPOUND"], "WHERE MAPID = '%s'" % (m[0]))
            subtotal_fs = self.database.count("(SELECT DISTINCT formula FROM map_f_%s_%s_%s WHERE MAPID = '%s')" % (strategy, self.ion_mode, self.name_peak_list, m[0]), "formula")

            subtotal_cpds = 0
            for compoundid in compoundids:
                str_compounds += "+" + compoundid[0]
                subtotal_cpds += 1
                
            url = "http://www.genome.jp/dbget-bin/show_pathway?%s%s" % (m[0], str_compounds)
            url_fs = "http://www.genome.jp/dbget-bin/show_pathway?%s" % (m[0])
            total_cpds = self.database.count("(SELECT DISTINCT KEGG_COMPOUND FROM map_compound WHERE MAPID = '%s')" % (m[0]), "KEGG_COMPOUND")
            total_fs = self.database.count("(SELECT DISTINCT formula FROM map_formula WHERE MAPID = '%s')" % (m[0]), "formula")

            self.database.execute("""INSERT INTO maps_coverage_%s_%s_%s
            SELECT '%s', mapname, '%s', '%s', '%s', '%s', 'C'
            FROM map_names WHERE mapid like '%s'
            ;""" % (strategy, self.ion_mode, self.name_peak_list, m[0], url, subtotal_cpds, total_cpds, round(float(subtotal_cpds)/total_cpds,2), "%"+str(m[0])[-5:]))

            self.database.cursor.execute("""INSERT INTO maps_coverage_%s_%s_%s
            SELECT '%s', mapname, '%s', '%s', '%s', '%s', 'F'
            FROM map_names WHERE mapid like '%s'
            ;""" % (strategy, self.ion_mode, self.name_peak_list, m[0], url_fs, subtotal_fs, total_fs, round(float(subtotal_fs)/total_fs,2), "%"+str(m[0])[-5:]))

        self.database.save()
        #self.database.close() 
            


    def _EFSjobs(self, mass, ppm, atoms, ions):
        calc = Calculate()
        min_max = calc.tolerance_peak(mass, ppm)
        atom_names = []
        atoms_min_limits = []
        atoms_max_limits = []
        mass_atoms = []
        final_out = {}
        for atom in atoms.limits:
            if max(atoms.limits[atom]) != 0:
                atoms_min_limits.append(min(atoms.limits[atom]))
                atoms_max_limits.append(max(atoms.limits[atom]))
                atom_names.append(atom)
                mass_atoms.append(atoms.lib[atom])
        index = -1
        
        for ion in ions.lib:
            atom_names.append(ion)     
            index += 1
            atoms_min_limits_inp = []
            atoms_max_limits_inp = []
            mass_atoms_inp = []

            temp = len(ions.lib)*[0]
            atoms_min_limits_inp.extend(atoms_min_limits)
            atoms_max_limits_inp.extend(atoms_max_limits)
            temp[index] = 1
            atoms_min_limits_inp.extend(temp)
            atoms_max_limits_inp.extend(temp)

            temp = len(ions.lib)*[0]
            mass_atoms_inp.extend(mass_atoms)
            temp[index] = ions.lib[ion]
            mass_atoms_inp.extend(temp)

            n = len(mass_atoms_inp)
            counts = [0] * n
            
            out = EFScpp.calculate(mass, 0.0, min(mass_atoms_inp), min_max[0], min_max[1], counts, atom_names, mass_atoms_inp, atoms_min_limits_inp, atoms_max_limits_inp, n-1, [])
            final_out[ion] = out
        return final_out



    def _EFS(self, mz_values, rules):

        self.total_peaks = len(mz_values)
        self.peaks_done, perc = 0, 0.0

        for mz in mz_values: # selection of mz values
            job = self._EFSjobs(mz, self.ppm, self.atoms, self.ions)
            for ion_name in job:
                for assignment in job[ion_name]:

                    peak_inp = assignment[0]
                    theo_mass = assignment[3]
                    
                    atom_names = assignment[1]
                    atom_counts = assignment[2]
                    
                    count = {}
                    
                    for i in range(len(atom_names)):
                        count[atom_names[i]] = atom_counts[i]
    
                    f = Formula("", self.atoms, self.ions, self.isotopes)
                    emp_formula = f.print_it(count)
    
                    f = Formula(emp_formula, self.atoms, self.ions, self.isotopes)
                    parent_formula = f.parent()

                    sumRules = sum([f.lewis_senior()["lewis"], f.lewis_senior()["senior"], f.HNOPS()["HC"], f.HNOPS()["NOPSC"]])

                    if (rules == True and sumRules == 4) or (rules == False):
                    
                        self.database.insert("EFS_%s_%s" % (self.ion_mode, self.name_peak_list),
                            {"mz":peak_inp, "intensity":self.peak_list.intensity(peak_inp), "ion":ion_name, 
                             "ppm_error":(peak_inp-theo_mass)/(theo_mass*0.000001), "mass":theo_mass, 
                             "exact_mass":theo_mass-self.ions.lib[ion_name],
                             "formula":emp_formula, 
                             "lewis":f.lewis_senior()["lewis"], "senior":f.lewis_senior()["senior"],
                             "HC":f.HNOPS()["HC"], "NOPSC":f.HNOPS()["NOPSC"], 
                             "rules":sumRules})
                    else:
                        pass
                    
            self.peaks_done += 1
            if float(self.peaks_done)/self.total_peaks >= perc:
                str_perc = str(perc * 100.0) + "%"
                perc += 0.1
                self.database.save()
                print str_perc,

    def _EFSpp(self, mz_values, rules, ncpus=2):

        ppservers = ()
        job_server = pp.Server(ncpus, ppservers=ppservers)

        #def create_subsets(l, n):
        #    return [l[i:i+n] for i in range(0, len(l), n)]

        print "Starting Parallel Python with", job_server.get_ncpus(), "workers"
        print "------------------------------------------------------------------"
        
        #subsets = create_subsets(mz_values, 100)

        length = 100
        subsets = [mz_values[i:i+length] for i in range(0, len(mz_values), length)]

        jobs = []

        self.total_peaks = len(mz_values)
        self.peaks_done, perc = 0, 0.0

        for subset in subsets:
            subset_jobs = []
            for mz in subset: # selection of mz values
                subset_jobs.append(job_server.submit(self._EFSjobs, (mz, self.ppm, self.atoms, self.ions,), (Calculate,), ("from MI_Pack import EFScpp",)))
            jobs.append(subset_jobs)

        for subset in jobs:
            for job in subset:
                job = job()

                for ion_name in job:

                    for assignment in job[ion_name]:
    
                        peak_inp = assignment[0]
                        theo_mass = assignment[3]
                        
                        atom_names = assignment[1]
                        atom_counts = assignment[2]
                        
                        count = {}
                        
                        for i in range(len(atom_names)):
                            count[atom_names[i]] = atom_counts[i]
        
                        f = Formula("", self.atoms, self.ions, self.isotopes)
                        emp_formula = f.print_it(count)
        
                        f = Formula(emp_formula, self.atoms, self.ions, self.isotopes)
                        parent_formula = f.parent()

                        sumRules = sum([f.lewis_senior()["lewis"], f.lewis_senior()["senior"], f.HNOPS()["HC"], f.HNOPS()["NOPSC"]])

                        if (rules == True and sumRules == 4) or (rules == False):
                        
                            self.database.insert("EFS_%s_%s" % (self.ion_mode, self.name_peak_list),
                                {"mz":peak_inp, "intensity":self.peak_list.intensity(peak_inp), "ion":ion_name, 
                                 "ppm_error":(peak_inp-theo_mass)/(theo_mass*0.000001), "mass":theo_mass, 
                                 "exact_mass":theo_mass-self.ions.lib[ion_name],
                                 "formula":emp_formula,
                                 "lewis":f.lewis_senior()["lewis"], "senior":f.lewis_senior()["senior"],
                                 "HC":f.HNOPS()["HC"], "NOPSC":f.HNOPS()["NOPSC"], 
                                 "rules":sumRules})
                        else:
                            pass
                self.peaks_done += 1

                if float(self.peaks_done)/self.total_peaks >= perc:
                    str_perc = str(perc * 100.0) + "%"
                    perc += 0.1
                    self.database.save()
                    print str_perc,


    def _EFSdb(self, mz_values, rules):

        try:
            conn = MySQLdb.connect('127.0.0.1', 'root', '', 'EFS');
            curMySQLdb = conn.cursor()

            ranges = []
            curMySQLdb.execute("show tables;")
            for table in curMySQLdb.fetchall():
                prefixes = table[0].split("__")
                ranges.append([float(prefixes[1].replace("_", ".")), float(prefixes[2].replace("_", "."))])
            
        except:
            print "Can't connect to MySQL server"
            return

        efs_out = []
        ef_id = 1
        calc = Calculate()

        self.total_peaks = len(mz_values)
        self.peaks_done, perc = 0, 0.0
    
        for mz in mz_values:
            
            for ion_name in self.ions.lib:

                min_max = calc.tolerance_peak(mz, self.ppm)

                minTolTemp, maxTolTemp = min_max[0] - self.ions.lib[ion_name], min_max[1] - self.ions.lib[ion_name]
                if minTolTemp >= 0 and maxTolTemp >= 0 and minTolTemp <= maxTolTemp:
                    tables = []
                    
                    for StartEnd in ranges:
                        table = "EFS__%s__%s" % (str(StartEnd[0]).replace(".", "_"), str(StartEnd[1]).replace(".", "_"))
                        if minTolTemp >= StartEnd[0] and minTolTemp <= StartEnd[1] and table not in tables:
                            tables.append(table)
                        if maxTolTemp >= StartEnd[0] and maxTolTemp <= StartEnd[1] and table not in tables:
                            tables.append(table)

                    if len(tables) > 0:
                        names = ["C","H","N","O","P","S"]
                        temp = []
                        for table in tables:
                            
                            if rules == True:
                                temp.append("""select * from %s where lewis = 1 and senior = 1 and HC = 1 and NOPSC = 1 and ExactMass >= %s
                                            and ExactMass <= %s %s""" % (table, minTolTemp, maxTolTemp, ""))
                                #temp.append("""select C,H,N,O,P,S,DoubleBondEquivalents,lewis,senior,HC,NOPSC,ExactMass+%s from %s where lewis = 1 and senior = 1 and HC = 1 and NOPSC = 1 and ExactMass >= %s
                                #            and ExactMass <= %s %s""" % (IonsLib[ion], table, minTolTemp, maxTolTemp, sql_filter))
                            else:
                                temp.append("select * from %s where ExactMass >= %s and ExactMass <= %s" % (table, minTolTemp, maxTolTemp))                        
                                #temp.append("select C,H,N,O,P,S,DoubleBondEquivalents,lewis,senior,HC,NOPSC,ExactMass+%s from %s where ExactMass >= %s and ExactMass <= %s" % (IonsLib[ion], table, minTolTemp, maxTolTemp))

                        #print " UNION ".join(map(str,temp))
                        curMySQLdb.execute(" UNION ".join(map(str,temp)))# + " order by ExactMass")
                        efs = curMySQLdb.fetchall()
                        for ef in efs:
                            #print {"DBE":ef[6], "LEWIS":ef[7], "SENIOR":ef[8], "HC":ef[9], "NOPSC":ef[10]}
                            #efs_out.append({"mass" : float(ef[-1])+IonsLib[ion], "atoms" : dict(zip(names, ef)), "ion" : ion, "DBE":ef[6], "LEWIS":ef[7], "SENIOR":ef[8], "HC":ef[9], "NOPSC":ef[10]})
                            ef_id += 1

                            theo_mass = float(ef[-1])+self.ions.lib[ion_name]

                            f = Formula("", self.atoms, self.ions, self.isotopes)
                            emp_formula = f.print_it(dict(zip(names, ef)))

                            self.database.insert("EFS_%s_%s" % (self.ion_mode, self.name_peak_list),
                            {"mz":mz, "intensity":self.peak_list.intensity(mz), "ion":ion_name, 
                             "ppm_error":(mz-float(theo_mass))/(float(theo_mass)*0.000001), "mass":theo_mass, 
                             "exact_mass":float(ef[-1]),
                             "formula":emp_formula,
                             "lewis":int(ef[7]), "senior":int(ef[8]),
                             "HC":int(ef[9]), "NOPSC":int(ef[10]), 
                             "rules":sum([int(ef[7]), int(ef[8]), int(ef[9]), int(ef[10])])})
                        
            self.peaks_done += 1
            if float(self.peaks_done)/self.total_peaks >= perc:
                str_perc = str(perc * 100.0) + "%"
                perc += 0.1
                self.database.save()
                print str_perc,



    def EFS(self, db=False, rules=False, ncpus=2, add_records=False):
  
        self.database.execute("""SELECT name FROM sqlite_master WHERE name='EFS_%s_%s'""" % (self.ion_mode, self.name_peak_list));

        if add_records == True and self.database.fetchone() != None:
            
            self.database.execute("""DROP INDEX IF EXISTS idx_mz_formula_%s_%s_EFS""" % (self.ion_mode, self.name_peak_list))
            self.database.execute("""DROP INDEX IF EXISTS idx_rules_mz_formula_%s_%s_EFS""" % (self.ion_mode, self.name_peak_list))
            
            self.database.execute("""select distinct mz from EFS_%s_%s order by mz""" % (self.ion_mode, self.name_peak_list))
            mz_values = []
            mzs_in_db = self.database.fetchall()
            for mz in self.peak_list.intensities:
                a = False
                for mz_in_db in mzs_in_db:
                    if mz == mz_in_db[0]:
                        a = True
                        break

                if a == False:
                    mz_values.append(mz)
            
        else:

            # Add additional columns for uncommon atoms
            #sql_str = ""
            #for atom in self.atoms.lib:
            #    if atom not in ["C", "H", "N", "O", "P", "S"]:
            #        sql_str += "'%s' INTEGER DEFAULT NULL, " % (atom)

            self.database.create_table("formulae", self.ion_mode, self.name_peak_list)   
            mz_values = self.peak_list.intensities.keys()

        print
        print "Empirical formula search - %s" % (self.fn_peak_list)
        print "------------------------------------------------------------------"
        
        if db == True:
            self._EFSdb(mz_values, rules)
        elif _has_pp == 0 or ncpus <= 1:
            self._EFS(mz_values, rules)
        elif _has_pp == 1 or ncpus > 1:
            self._EFSpp(mz_values, rules, ncpus)
        else:
            print "Check parameters"

        self.database.save()
        self.database.index("mz_formula_%s_%s_EFS" % (self.ion_mode, self.name_peak_list),
            "EFS_%s_%s" % (self.ion_mode, self.name_peak_list),["mz", "formula"])
                
        self.database.index("rules_mz_formula_%s_%s_EFS" % (self.ion_mode, self.name_peak_list),
            "EFS_%s_%s" % (self.ion_mode, self.name_peak_list),["rules", "mz", "formula"])

        self.database.save()
        #self.database.close()  
        print
        print "------------------------------------------------------------------"


    def PPS(self, library, min_peaks=2, fn_SNR=None, offset=-1.0):

        if fn_SNR != None and fn_SNR != "":
            self.SNR_values = PeakList(fn_SNR, True)
            
        def Assign(sort_ids, lib_masses, pattern_type, abundance):
            print
            print "Assigning peak patterns - %s - (%s)" % (self.fn_peak_list, ", ".join(map(str,abundance.keys())))
            print "------------------------------------------------------------------"

            if "Isotopes" in str(library.__class__):
                isotopes_count = self.database.count("PPS_%s_%s" % (self.ion_mode, self.name_peak_list), "label_a",
                "where type == 'isotopes' and ((label_a = '%s' and label_b = '%s') or (label_b = '%s' and label_a = '%s'))" % (abundance.keys()[0], abundance.keys()[1], abundance.keys()[0], abundance.keys()[1]))
                
                if isotopes_count > 0:
                    self.database.delete("PPS_%s_%s" % (self.ion_mode, self.name_peak_list),
                    "where type == 'isotopes' and ((label_a = '%s' and label_b = '%s') or (label_b = '%s' and label_a = '%s'))" % (abundance.keys()[0], abundance.keys()[1], abundance.keys()[0], abundance.keys()[1]))
                    self.database.save()
                    
            check_ions_count = 0
            if "Ions" in str(library.__class__):
                for ion in library.lib:
                    ions_count = self.database.count("PPS_%s_%s" % (self.ion_mode, self.name_peak_list), "label_a",
                    "where type == 'ions' and (label_a = '%s' or label_b = '%s')" % (ion, ion)) 
                    if ions_count > 0:
                        check_ions_count += 1
                
                if check_ions_count == len(library.lib) or (len(library.lib) > check_ions_count and check_ions_count >= 2):
                    for ion in library.lib:
                        self.database.delete("PPS_%s_%s" % (self.ion_mode, self.name_peak_list),
                        "where type == 'ions' and (label_a = '%s' or label_b = '%s')" % (ion, ion))
                        self.database.save()
            
            pattern_id = self.database.select("PPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["max(id)"])
            if pattern_id[0][0] == None:
                pattern_id = 0
            else:
                pattern_id = int(pattern_id[0][0])
            
            self.total_peaks =  len(sort_ids)
            self.peaks_done, perc = 0, 0.0
            copy_sort_ids = deepcopy(sort_ids)
            
            c = 0
            for peak_id in sort_ids:
                c += 1
                self.peaks_done += 1
                if peak_id in copy_sort_ids: # check if not already assigned - avoiding double records
                    
                    peak_error = self.calc.abs_error(lib_masses[peak_id][1], self.ppm)
                    min_tol = float(lib_masses[peak_id][0]) - peak_error
                    max_tol = float(lib_masses[peak_id][0]) + peak_error
    
                    pattern = []
                    ions_pattern = []
                    for comp_peak_id in copy_sort_ids:
                        if lib_masses[comp_peak_id][0] >= min_tol and lib_masses[comp_peak_id][0] <= max_tol:
                            pattern.append(comp_peak_id)
                            ions_pattern.append(lib_masses[comp_peak_id][2])
    
                    if len(pattern) >= min_peaks: # check number of peaks in pattern
                        inserted = 0
                        pattern_id += 1
                        for i in pattern:
                            copy_sort_ids.remove(i) # remove peak id to avoid duplications
                            for j in pattern:

                                if lib_masses[i][1] - lib_masses[j][1] > 0.0:

                                    rules_out = 1
                                    n_atoms = "NULL"
                                    mean_atoms = "NULL"
                                    sd_atoms = "NULL"
                                    SNR_mean_a = "NULL"
                                    SNR_mean_b = "NULL"

                                    if lib_masses[j][3] > lib_masses[i][3] and abundance.values()[0] < abundance.values()[1] and pattern_type == "isotopes":
                                        rules_out = 0

                                    if lib_masses[j][3] < lib_masses[i][3] and abundance.values()[0] > abundance.values()[1] and pattern_type == "isotopes":
                                        rules_out = 0
                                        
                                    if pattern_type == "isotopes" and rules_out != 0 and "(13C)" in ions_pattern: # check inensities isotopes if isotope pattern
                                        max_intensity = (-0.0213 * lib_masses[i][1]) + 100.24
                                        min_intensity = (-0.0768 * lib_masses[i][1]) + 98.693
                                        check_inten = 100 - (lib_masses[i][3]/(lib_masses[i][3]+lib_masses[j][3]))*100
                                    
                                        if check_inten >= min_intensity and check_inten <= max_intensity:
                                            rules_out = 1
                                        else:
                                            rules_out = 0

                                        n_atoms = round((lib_masses[i][3]*abundance[lib_masses[j][2]])/(lib_masses[j][3]*abundance[lib_masses[i][2]]),1)
                                        
                                    if pattern_type == "isotopes" and rules_out != 0 and "(13C)" not in ions_pattern: # check inensities isotopes if isotope pattern
                                        if abundance.values()[0] > abundance.values()[1]:
                                            n_atoms = round((lib_masses[i][3]*abundance[lib_masses[j][2]])/(lib_masses[j][3]*abundance[lib_masses[i][2]]),1)
                                        else: # e.g. Lithium
                                            n_atoms = round((lib_masses[j][3]*abundance[lib_masses[i][2]])/(lib_masses[i][3]*abundance[lib_masses[j][2]]),1)
                                        
                                    if pattern_type == "isotopes" and rules_out != 0 and fn_SNR != None and fn_SNR != "":
                                        ########################################
                                        # SNR VALUES
                                        ########################################
                                        la = self.SNR_values.snrs[lib_masses[i][1]]
                                        lb = self.SNR_values.snrs[lib_masses[j][1]]

                                        temp = []
                                        for snr_value in range(len(la)):
                                            if la[snr_value] != 0.0 and lb[snr_value] != 0.0: 
                                                #print la[snr_value], lb[snr_value]
                                                #print la[snr_value], lb[snr_value], (la[snr_value]*abundance[lib_masses[j][2]])/(lb[snr_value]*abundance[lib_masses[i][2]])
                                                n_atoms_SNR = round((la[snr_value]*abundance[lib_masses[j][2]])/(lb[snr_value]*abundance[lib_masses[i][2]]),1)
                                                temp.append(n_atoms_SNR)

                                        if len(temp) > 1:
                                            mean_atoms, sd_atoms = Calculate().mean_stdv(temp)
                                            mean_atoms_incl_offset = mean_atoms + offset

                                            SNR_mean_a, SNR_sd_a = Calculate().mean_stdv(la)
                                            SNR_mean_b, SNR_sd_b = Calculate().mean_stdv(lb)
                                            #print mean_atoms_incl_offset, mean_atoms, sd_atoms, sd_atoms*3
                                            #print mean_atoms - sd_atoms*3, mean_atoms + sd_atoms*3, rules_out, check_inten, "=>", min_intensity, check_inten, "<=", max_intensity, check_inten >= min_intensity and check_inten <= max_intensity, ions_pattern, len(set(ions_pattern)), len(ions_pattern)
                                            #raw_input()

                                        ########################################
                                        # SNR VALUES
                                        ########################################

                                    self.database.insert("PPS_%s_%s" %(self.ion_mode, self.name_peak_list),
                                        {"id":pattern_id, "mz_pair_diff":lib_masses[i][1] - lib_masses[j][1],
                                        "mz_a":lib_masses[j][1], "mz_b":lib_masses[i][1],
                                        "intensity_a":lib_masses[j][3], "intensity_b":lib_masses[i][3],
                                        "SNR_mean_a":SNR_mean_a, "SNR_mean_b":SNR_mean_b,
                                        "label_a":lib_masses[j][2], "label_b":lib_masses[i][2],
                                        "type":pattern_type, "num_peaks":len(pattern),
                                        "n_atoms": n_atoms, "mean_atoms": mean_atoms, "sd_atoms": sd_atoms, "rules":rules_out}) 
                                    inserted = 1
                        if inserted == 0:
                            pattern_id -= 1
                if float(self.peaks_done)/self.total_peaks >= perc:
                    print str(perc*100)+"%",
                    perc += 0.1
                    self.database.save()
            
            self.database.save()
            print

        table_flag = self.database.select("sqlite_master", ["name"], "WHERE type='table' and name = 'PPS_%s_%s'" % (self.ion_mode, self.name_peak_list))
        if len(table_flag) == 0:
            self.database.create_table("peak_patterns", self.ion_mode, self.name_peak_list)       

        if "Ions" in str(library.__class__):
            pattern_type = "ions"
            peak_id = 0
            lib_masses = {}
            sort_ids = []
            abundance = OrderedDict()

            for ion in library.lib:
                abundance[ion] = 0.0
                for mz in self.peak_list.intensities:
                                
                    peak_min_ion = mz - library.lib[ion]
                    
                    if peak_min_ion > 0:                   
                        peak_id += 1
                        lib_masses[peak_id] = [peak_min_ion, mz, ion, self.peak_list.intensity(mz), library.lib[ion]]
    
                        #### sort peak list from high to low based on original peak
                        
                        peak_error = self.calc.abs_error(mz, self.ppm)
                        min_tol = peak_min_ion - peak_error
                        max_tol = peak_min_ion + peak_error
    
                        sort_ids.append([peak_id, min_tol+max_tol])
                   
            sort_ids.sort(lambda x, y: cmp(y[1],x[1]))
    
            for x in range(len(sort_ids)):
                sort_ids[x] = sort_ids[x][0]
            
            Assign(sort_ids, lib_masses, pattern_type, abundance)
    
            ############################################################

        if "Isotopes" in str(library.__class__):
            pattern_type = "isotopes"
            for ion in library.lib:
                peak_id = 0
                lib_masses = {}
                sort_ids = []
                
                for mz in self.peak_list.intensities:
                                
                    peak_min_ion = mz - ion[2]
                    
                    if peak_min_ion > 0:
                        peak_id += 1
                        lib_masses[peak_id] = [mz, mz, ion[0], self.peak_list.intensity(mz), 0.0]
                        peak_error = self.calc.abs_error(mz, self.ppm)
                        min_tol = peak_min_ion - peak_error
                        max_tol = peak_min_ion + peak_error
                        sort_ids.append([peak_id, min_tol+max_tol])
                   
                        peak_id += 1
                        lib_masses[peak_id] = [peak_min_ion, mz, ion[1], self.peak_list.intensity(mz), ion[2]]
                        peak_error = self.calc.abs_error(mz, self.ppm)
                        min_tol = peak_min_ion - peak_error
                        max_tol = peak_min_ion + peak_error
                        sort_ids.append([peak_id, min_tol+max_tol])   

                sort_ids.sort(lambda x, y: cmp(y[1],x[1]))
    
                for x in range(len(sort_ids)):
                    sort_ids[x] = sort_ids[x][0]
                
                abundance = OrderedDict([(ion[0],ion[3]), (ion[1],ion[4])])
                Assign(sort_ids, lib_masses, pattern_type, abundance)
            
        try:
            self.database.index("id_%s_%s_PPS" % (self.ion_mode, self.name_peak_list),
                "PPS_%s_%s" % (self.ion_mode, self.name_peak_list),["id"])

            self.database.index("mz_ab_%s_%s_PPS" % (self.ion_mode, self.name_peak_list),
                "PPS_%s_%s" % (self.ion_mode, self.name_peak_list),["mz_a", "mz_b"])

            self.database.index("many_%s_%s_PPS" % (self.ion_mode, self.name_peak_list),
                "PPS_%s_%s" % (self.ion_mode, self.name_peak_list),
                ["rules", "intensity_a", "intensity_b", "num_peaks", "mz_a", "mz_b"])
        except:
            pass
        
        #self.database.close()  
        print
        print "------------------------------------------------------------------"
        print "Finished"
        print

    def output(self, strategy, max_records=10000, fn_out = None):

        print
        print "Creating output - %s" % (self.fn_peak_list)
        print "------------------------------------------------------------------"

        db_names = []
        if strategy == "SPS":
            
            self.fn_output = os.path.join(self.path, "MI_Pack__output_SPS_%s_%s.txt" % (self.ion_mode, self.name_peak_list))

            ######### MIDB #########################
            # WHICH DATABASES ARE SELECTED
            ########################################         
            self.database.execute("PRAGMA table_info(SPS_%s_%s)" % (self.ion_mode, self.name_peak_list))
            col_info = self.database.fetchall()
            for entry in col_info[8:]:
                db_names.append(str(entry[1]))
            ########################################
            
        elif strategy == "TM":
            
            self.fn_output = os.path.join(self.path, "MI_Pack__output_TM_%s_%s.txt" % (self.ion_mode, self.name_peak_list))
            
            ######### MIDB #########################
            # WHICH DATABASES ARE SELECTED
            ########################################              
            db_names = ["KEGG_COMPOUND"] # TM is KEGG specific
            ########################################
                        
        elif strategy == "EFS":
            
            self.fn_output = os.path.join(self.path, "MI_Pack__output_EFS_%s_%s.txt" % (self.ion_mode, self.name_peak_list))
        
        else:
            print "Wrong input, please use SPS, TM or EFS to create output"
            sys.exit()

        if fn_out != None: # GALAXY
            self.fn_output = fn_out # GALAXY
            
        self.output = open(self.fn_output, "w")

        # check if tables exist in database
        tables_selected = self.database.select("sqlite_master", ["name"], "WHERE type='table'")
        tables = []
        for t in tables_selected:
            tables.append(str(t["name"]))

        EFS, PPS, records_isotopes = 0, 0, []

        if "SPS_%s_%s" % (self.ion_mode, self.name_peak_list) not in tables and strategy == "SPS":
            print "SPS data not available"
            return

        elif "TM_%s_%s" % (self.ion_mode, self.name_peak_list) not in tables and strategy == "TM":
            print "TM data not available"
            return

        elif "EFS_%s_%s" % (self.ion_mode, self.name_peak_list) not in tables and strategy == "EFS":
            print "EFS data not available"
            return

        if "PPS_%s_%s" % (self.ion_mode, self.name_peak_list) in tables:
            records_isotopes = self.database.select("PPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct label_a", "label_b"], "where type = 'isotopes' and rules == 1")
            PPS = 1
        if "EFS_%s_%s" % (self.ion_mode, self.name_peak_list) in tables:
            EFS = 1
        # END - check if tables exist in database
        

        # Read input (e.g. m/z, Intensity, Fold change, p-value and PC loading)
        stats = open(self.fn_peak_list, "r")
        stats.readline()
        lines = stats.readlines()
        dict_stats = OrderedDict()
        for line in lines:
            line = line.replace("\n","").split("\t")
            if len(line) > 2:
                dict_stats[float(line[0]),float(line[1]), line[2]] = line[1:]
            else:
                dict_stats[float(line[0])] = line[1:]
        # END - Read input
        
        
        # Salt clusters [NaCl] & [Na(37Cl)] (Temporary Solution)
        salts = []
        for i in range(0,6):
            for ii in range(0,6):
                if i >= ii or i == 0:
                    if i == 1:
                        temp_i = ""
                    else:
                        temp_i = i
                    if ii == 1:
                        temp_ii = ""
                    else:
                        temp_ii = ii
                    if ii >= 1 and i >= 1: 
                        salts.append("[NaCl]%s[Na(37Cl)]%s" % (temp_i, temp_ii))
                    if i >= 1 and ii == 0:
                        salts.append("[NaCl]%s" % (temp_i))
                    if ii >= 1 and i == 0:
                        salts.append("[Na(37Cl)]%s" % (temp_ii))

        salts_sql = """ and formula in ('%s')""" % ("', '".join(map(str,salts)))
        # END - Salt clusters [NaCl] & [Na(37Cl)] (Temporary Solution)
        
        temp_parent_formulae = [] # to count the number of daughter formulae (including parent formula)
        self.rows = []

        self.row_default = OrderedDict([("input",""),
            ("no_emp_formula", 0),
            ("no_metabolite_names", 0),
            ("EFS", 0),
            ("strategy", ""),
            ("emp_formula_parent", ""),
            ("no_emp_formulae_related_to_parent", 0),
            ("emp_formula_peak", ""),
            ("ion_form", ""),
            ("exact_mass", 0.0),
            ("exact_adduct_mass", 0.0),
            ("mass_error", 0.0)])

        for row in dict_stats:

            if type(row) == float:
                mz = row
                stats_out = "%s\t%s" % (mz, "\t".join(map(str,dict_stats[mz]))) # stats data peaklist            
            else:
                mz = row[0]
                stats_out = "%s\t%s" % (mz, "\t".join(map(str,dict_stats[mz, row[1], row[2]]))) # stats data peaklist

            self.row = copy.copy(self.row_default)
            self.row["input"] = stats_out

            # COLLECTING isotopes & adducts
            sql_like_isotopes = []
            if len(sql_like_isotopes) == 0:
                for ip in self.isotopes.lib:
                    sql_like = "and formula not like '%" + ip[1] + "%'" # DO NOT PRINT ISOTOPES > NO ISOTOPE PATTERN FOUND!!
                    if sql_like not in sql_like_isotopes:
                        sql_like_isotopes.append(str(sql_like)) 

            # SELECT ISOTOPES AND PRINT HEADERS FOR EACH ISOTOPE PAIR            
            if len(records_isotopes) > 0 or "PPS_%s_%s" % (self.ion_mode, self.name_peak_list) in tables:
                ip_out = OrderedDict()
                                
                for ip in records_isotopes:
                    ip_out["%s-%s" % (ip[0], ip[1])] = []
                    
                    isotope_records_selected = self.database.select("PPS_%s_%s" % (self.ion_mode, self.name_peak_list), ['mz_a', 'n_atoms', 'id'],
                    "where (mz_a = %s or mz_b = %s) and rules = 1 and label_a = '%s' and label_b = '%s'" % (mz, mz, ip[0], ip[1]))

                    if len(isotope_records_selected) > 0:
                        for assignment in isotope_records_selected:
                            ip_out["%s-%s" % (ip[0], ip[1])].append([assignment[1], assignment[2]])
                            
                            if assignment[0] != mz:
                                sql_like = "and formula not like '%" + ip[1] + "%'"
                                if sql_like in sql_like_isotopes:
                                    sql_like_isotopes.remove(str(sql_like))
                            
                isotope_info_header = len(ip_out) * 2 * [""]
                isotope_info = len(ip_out) * 2 * [""]

                index = -1
                for pair in ip_out:
                    index += 1
                    
                    #### header
                    isotope_info_header[index] = pair
                    isotope_info_header[index+(len(ip_out))] = pair + "_ID"
                    #### header
    
                    #### calculations
                    for a in ip_out[pair]:
                        if isotope_info[index] != "":
                            isotope_info[index] += "," + str(a[0]) + " "
                            isotope_info[index+(len(ip_out))] += "," + str(a[1]) + " "
                        else:
                            isotope_info[index] += str(a[0])
                            isotope_info[index+(len(ip_out))] += str(a[1])                        
                    #### calculations
                isotope_info_header = "\t".join(map(str,isotope_info_header)) + "\t"
                self.row["isotopes"] = "\t".join(map(str,isotope_info))
                self.row_default["isotopes"] = "\t".join(map(str,isotope_info))
            else:
                isotope_info_header = ""

            self.row["databases"] = ""

            EFS_flag, MIDB_flag = 0, 0
            records_f = []

            if EFS == 1:# empirical formulae table

                extra_records_f = self.database.select("EFS_%s_%s" % (self.ion_mode, self.name_peak_list), ["formula", "ion", "exact_mass", "mass", "ppm_error"],
                "where mz= %s %s" % (mz, salts_sql))
                
                records_f = self.database.select("EFS_%s_%s" % (self.ion_mode, self.name_peak_list), ["formula", "ion", "exact_mass", "mass", "ppm_error"],
                "where rules = 4 and mz= %s %s order by exact_mass" % (mz, " ".join(map(str,sql_like_isotopes))))         
                
                for record in extra_records_f:
                    records_f.append(record)
                
                if len(records_f) > 0:
                    EFS_flag = 1
            
            if len(records_f) > 0 and len(records_f) <= max_records:

                formula_to_check = []
                for record_f in records_f:
    
                    f = Formula(str(record_f["formula"]), self.atoms, self.ions, self.isotopes)
                    parent_formula = f.parent()
                    
                    if strategy == "SPS" or strategy == "TM":
                        MIDB_flag = 0
                        
                        temp_names = OrderedDict()
                        temp_names_count = OrderedDict()
                        for db_name in db_names:
                            if db_name not in temp_names:
                                temp_names[db_name] = []
                                temp_names_count[db_name] = 0
                        
                        formula_to_check.append(str(record_f["formula"])) # RULES SEE BELOW

                        records_f_MIDB = []

                        if strategy == "SPS":                            
                            records_f_MIDB = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT name"],
                            "WHERE formula = '%s'" % (record_f["formula"]))

                        elif strategy == "TM":
                            records_f_MIDB = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT name"],
                            "WHERE formula = '%s'" % (record_f["formula"])) #  and type = 1

                        if len(records_f_MIDB) > 0:
                            MIDB_flag = 1
                            for db_name in db_names:
                                for record_MIDB in records_f_MIDB:
                                    if strategy == "SPS":

                                        check_db_name = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (record_MIDB["name"], db_name))
                                        if len(check_db_name) > 0:
                                            temp_names[db_name].append(str(record_MIDB["name"]))
                                        else:
                                            temp_names[db_name].append("NA")
                                            
                                    elif strategy == "TM":
                                        check_db_name = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (record_MIDB["name"], db_name))   #  and type = 1                            
                                        if len(check_db_name) > 0:
                                            temp_names[db_name].append(str(record_MIDB["name"]))
                                        else:
                                            temp_names[db_name].append("NA")
                                temp_names_count[db_name] = "%s(%s)" % (len(temp_names[db_name]), len(temp_names[db_name]) - temp_names[db_name].count("NA"))
                    
                        self.row["no_emp_formula"] = len(records_f)
                        self.row["EFS"] = EFS_flag
                        self.row["strategy"] = MIDB_flag
                        self.row["emp_formula_parent"] = parent_formula
                        self.row["no_emp_formulae_related_to_parent"] = 0
                        self.row["emp_formula_peak"] = record_f["formula"]
                        self.row["ion_form"] = record_f["ion"]
                        self.row["exact_mass"] = record_f['exact_mass']
                        self.row["exact_adduct_mass"] = record_f["mass"]
                        self.row["mass_error"] = round(float(record_f["ppm_error"]),2)
                        
                            
                        if len(temp_names.values()[0]) > 0 and len(temp_names.values()[0])<= max_records:
                            self.row["no_metabolite_names"] = len(records_f_MIDB)
                            self.row["databases"] = "\t".join(map(str,temp_names.values()))
                        elif len(temp_names.values()[0]) > max_records:
                            self.row["no_metabolite_names"] = len(records_f_MIDB)
                            self.row["databases"] = "\t".join(map(str,temp_names_count.values()))
                        elif len(temp_names.values()[0]) == 0:
                            self.row["no_metabolite_names"] = 0
                            self.row["databases"] = "\t".join(map(str,[0]*len(db_names)))

                        #print self.row.values(), 1
                        
                        self.rows.append(self.row) ####################################################
                        #print 1, self.row_default
                        self.row = copy.copy(self.row_default)
                        self.row["input"] = stats_out
                        
                    else:
                        self.row["no_emp_formula"] = len(records_f)
                        self.row["no_metabolite_names"] = 0
                        self.row["EFS"] = EFS_flag
                        self.row["strategy"] = MIDB_flag
                        self.row["emp_formula_parent"] = parent_formula
                        self.row["no_emp_formulae_related_to_parent"] = 0
                        self.row["emp_formula_peak"] = record_f["formula"]
                        self.row["ion_form"] = record_f["ion"]
                        self.row["exact_mass"] = record_f['exact_mass']
                        self.row["exact_adduct_mass"] = record_f["mass"]
                        self.row["mass_error"] = round(float(record_f["ppm_error"]),2)
                        self.row["databases"] = "\t".join(map(str,[0]*len(db_names)))

                        #print self.row.values(), 2
                        
                        self.rows.append(self.row) ####################################################
                        #print 2, self.row_default
                        self.row = copy.copy(self.row_default)
                        self.row["input"] = stats_out

                    temp_parent_formulae.append(parent_formula)
                
                #####################################################
                # IF FORMULA IN MIDB THAT DOES NOT PASS THE RULES!!!#
                #####################################################
                if strategy == "SPS" or strategy == "TM":
                    
                    formulae_MIDB = []
                    #### SELECT FORMULAE THAT DID NOT PASS THE RULES
                    if strategy == "SPS":
                        formulae_MIDB = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula"],
                        "WHERE mz = %s and formula not in ('%s')" % (str(mz), "', '".join(map(str,formula_to_check))))
                    elif strategy == "TM":
                        formulae_MIDB = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula"],
                        "WHERE mz = %s and formula not in ('%s')" % (str(mz), "', '".join(map(str,formula_to_check)))) # and type = 1
                    #### SELECT FORMULAE THAT DID NOT PASS THE RULES
                    
                    # correct counts empirical formulae
                    for row in self.rows:
                        if row["input"] == self.rows[-1]["input"]:
                            row["no_emp_formula"] = row["no_emp_formula"] + len(formulae_MIDB)
                    # correct counts empirical formulae
                    
                    
                    for formula_MIDB in formulae_MIDB:
                        if strategy == "SPS":
                            records_f_MIDB = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula, exact_mass, ion, mass, ppm_error, name"],
                            "WHERE formula = '%s' and formula not in ('%s')" % (formula_MIDB["formula"], "', '".join(map(str,formula_to_check))))

                        elif strategy == "TM":
                            records_f_MIDB = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula, exact_mass, ion, mass, ppm_error, name"],
                            "WHERE formula = '%s' and formula not in ('%s')" % (formula_MIDB["formula"], "', '".join(map(str,formula_to_check)))) #  and type = 1
                    
                        temp_names = OrderedDict()
                        temp_names_count = OrderedDict()
                        for db_name in db_names:
                            if db_name not in temp_names:
                                temp_names[db_name] = []
                                temp_names_count[db_name] = 0

                        if len(records_f_MIDB) > 0:
                            MIDB_flag = 1
                            
                            for db_name in db_names:
                                for record_MIDB in records_f_MIDB:
                                    if strategy == "SPS":
                                        check_db_name = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (record_MIDB["name"], db_name))
                                        
                                        
                                        if len(check_db_name) > 0:
                                            temp_names[db_name].append(str(record_MIDB["name"]))
                                        else:
                                            temp_names[db_name].append("NA")
                                            
                                    elif strategy == "TM":
                                        check_db_name = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (record_MIDB["name"], db_name)) # and type = 1                          
                                        if len(check_db_name) > 0:
                                            temp_names[db_name].append(record_MIDB["name"])
                                        else:
                                            temp_names[db_name].append("NA")
                                temp_names_count[db_name] = "%s(%s)" % (len(temp_names[db_name]), len(temp_names[db_name]) - temp_names[db_name].count("NA"))

                        
                        self.row["no_emp_formula"] = self.rows[-1]["no_emp_formula"]
                        self.row["no_metabolite_names"] = len(records_f_MIDB)
                        self.row["EFS"] = 0
                        self.row["strategy"] = MIDB_flag
                        
                        if len(temp_names.values()[0]) > 0 and len(temp_names.values()[0])<= max_records:
                            temp_parent_formulae.append(record_MIDB["formula"])

                            self.row["emp_formula_parent"] = record_MIDB["formula"]
                            self.row["no_emp_formulae_related_to_parent"] = 0
                            self.row["emp_formula_peak"] = record_MIDB["formula"]
                            self.row["ion_form"] = record_MIDB["ion"]
                            self.row["exact_mass"] = record_MIDB["exact_mass"]
                            self.row["exact_adduct_mass"] = record_MIDB["mass"]
                            self.row["mass_error"] = round(float(record_MIDB["ppm_error"]),2)
                            self.row["databases"] = "\t".join(map(str,temp_names.values()))
                            
                            #print self.row.values(), 3

                            self.rows.append(self.row) ####################################################
                            #print 3, self.row_default
                            self.row = copy.copy(self.row_default)
                            self.row["input"] = stats_out
                            
                        elif len(temp_names.values()[0]) > max_records:
                            temp_parent_formulae.append(record_MIDB["formula"]) 

                            self.row["emp_formula_parent"] = ""
                            self.row["no_emp_formulae_related_to_parent"] = ""
                            self.row["emp_formula_peak"] = ""
                            self.row["ion_form"] = ""
                            self.row["exact_mass"] = ""
                            self.row["exact_adduct_mass"] = ""
                            self.row["mass_error"] = ""
                            self.row["databases"] = "\t".join(map(str,temp_names_count.values()))

                            #print self.row.values(), 4
                        

                            self.rows.append(self.row) ####################################################
                            #print 4, self.row_default
                            self.row = copy.copy(self.row_default)
                            self.row["input"] = stats_out
                                                
                        elif len(temp_names.values()[0]) == 0:
                            print "#########################"
                            print "len(temp_names.values()[0]) == 0:"
                            raw_input()

                ###########################################################
                # END - IF FORMULA IN MIDB THAT DOES NOT PASS THE RULES!!!#
                ###########################################################

            ################################################################
            # END - IF EFS IS NOT APPLIED OR TO MANY NO. OF EF ASSIGNMENTS ARE FOUND!!!#
            ################################################################
            else:
                ######## PARENT FORMULAE ########################################################
                if len(records_f) != 0:
                    add_f_records = False
                    for record_f in records_f:
                        f = Formula(str(record_f["formula"]), self.atoms, self.ions, self.isotopes)
                        parent_formula = f.parent()
                        temp_parent_formulae.append(parent_formula)
                else:
                    add_f_records = True 
                ######## PARENT FORMULAE ########################################################
                
                if strategy == "SPS" or strategy == "TM":
                    records_f_MIDB = []
                    #### SELECT ALL FORMULAE ASSIGNMENTS ####
                    if strategy == "SPS":
                        records_f_MIDB = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula, ion, mass, ppm_error, exact_mass"],
                        "WHERE mz= '%s'" % (mz))
                    elif strategy == "TM":
                        records_f_MIDB = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["DISTINCT formula, ion, mass, ppm_error, exact_mass"],
                        "WHERE mz= '%s'" % (mz)) # AND type = 1
                    
                    if len(records_f_MIDB) > 0:
                        MIDB_flag = 1

                        for record_MIDB in records_f_MIDB:

                            if add_f_records == True:                           
                                temp_parent_formulae.append(record_MIDB["formula"])

                            #### CHECK TOTAL NO. OF NAME ASSIGNMENTS ####
                            if strategy == "SPS":
                                records_names_MIDB = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct name"],
                                "WHERE formula = '%s'" % (str(record_MIDB["formula"])))
                            elif strategy == "TM":
                                records_names_MIDB = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct name"],
                                "WHERE formula = '%s'" % (str(record_MIDB["formula"]))) #  AND type = 1
    
                            temp_names_count = OrderedDict()
                            temp_names = OrderedDict()                        
                            for db_name in db_names:
                                if db_name not in temp_names:
                                    temp_names_count[db_name] = 0
                                    temp_names[db_name] = []
                            
                                for record_MIDB_name in records_names_MIDB:
                                    #### CHECK DB NAME ####
                                    if strategy == "SPS":
                                        
                                        check_db_name = self.database.select("SPS_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (str(record_MIDB_name["name"]), db_name))
                                    elif strategy == "TM":
                                        check_db_name = self.database.select("TM_%s_%s" % (self.ion_mode, self.name_peak_list), ["distinct %s" % (db_name)],
                                        """WHERE name = "%s" AND %s is not NULL""" % (str(record_MIDB_name["name"]), db_name)) #  and type = 1
                                    
                                    if len(check_db_name) > 0:
                                        temp_names[db_name].append(str(record_MIDB_name["name"]))
                                    else:
                                        temp_names[db_name].append("NA")                        
                                    
                                temp_names_count[db_name] = "%s(%s)" % (len(temp_names[db_name]), len(temp_names[db_name]) - temp_names[db_name].count("NA"))

                            ############## SPS ONLY ############################ 
                            if len(records_f) == 0:
                                records_f = records_f_MIDB
                            ############## SPS ONLY ############################

                            self.row["no_emp_formula"] = len(records_f)
                            self.row["no_metabolite_names"] = len(records_names_MIDB)
                            self.row["EFS"] = EFS_flag
                            self.row["strategy"] = MIDB_flag
                            
                            if len(temp_names_count.values()[0]) > 0 and len(temp_names.values()[0]) <= max_records: # CHECK NO. OF NAME ASSIGNMENTS
                                self.row["emp_formula_parent"] = record_MIDB["formula"]
                                self.row["no_emp_formulae_related_to_parent"] = 0
                                self.row["emp_formula_peak"] = record_MIDB["formula"]
                                self.row["ion_form"] = record_MIDB["ion"]
                                self.row["exact_mass"] = record_MIDB['exact_mass']
                                self.row["exact_adduct_mass"] = record_MIDB["mass"]
                                self.row["mass_error"] = round(float(record_MIDB["ppm_error"]),2)
                                self.row["databases"] = "\t".join(map(str,temp_names.values()))

                                #print self.row.values(), 5
                        

                                self.rows.append(self.row) ####################################################
                                #print 5, self.row_default
                                self.row = copy.copy(self.row_default)
                                self.row["input"] = stats_out
                                
                            elif len(temp_names_count.values()[0]) > max_records:
                                self.row["emp_formula_parent"] = ""
                                self.row["no_emp_formulae_related_to_parent"] = 0
                                self.row["emp_formula_peak"] = len(records_f)
                                self.row["ion_form"] = ""
                                self.row["exact_mass"] = ""
                                self.row["exact_adduct_mass"] = ""
                                self.row["mass_error"] = ""
                                self.row["databases"] = "\t".join(map(str,temp_names_count.values()))

                                #print self.row.values(), 6
                        
                                self.rows.append(self.row) ####################################################
                                #print 6, self.row_default
                                self.row = copy.copy(self.row_default)
                                self.row["input"] = stats_out
                                
                            elif len(temp_names_count.values()[0]) == 0:
                                self.row["emp_formula_parent"] = ""
                                self.row["no_emp_formulae_related_to_parent"] = 0
                                self.row["emp_formula_peak"] = len(records_f)
                                self.row["ion_form"] = ""
                                self.row["exact_mass"] = ""
                                self.row["exact_adduct_mass"] = ""
                                self.row["mass_error"] = ""
                                self.row["databases"] = "\t".join(map(str,[0]*len(db_names)))
                                
                                #print self.row.values(), 7
                            
                                self.rows.append(self.row) ####################################################
                                #print 7, self.row_default
                                self.row = copy.copy(self.row_default)
                                self.row["input"] = stats_out
                                
                            else:
                                print "###############"
                                raw_input
                            
                            
                    
                    elif len(records_f_MIDB) == 0:
                        self.row["no_emp_formula"] = len(records_f)
                        self.row["no_metabolite_names"] = 0
                        self.row["EFS"] = EFS_flag
                        self.row["strategy"] = MIDB_flag
                        self.row["emp_formula_parent"] = ""
                        self.row["no_emp_formulae_related_to_parent"] = 0
                        self.row["emp_formula_peak"] = len(records_f)
                        self.row["ion_form"] = ""
                        self.row["exact_mass"] = ""
                        self.row["exact_adduct_mass"] = ""
                        self.row["mass_error"] = ""
                        self.row["databases"] = "\t".join(map(str,[0]*len(db_names)))
                        
                        #print self.row.values(), 8
                        
                        self.rows.append(self.row) ####################################################
                        #print 8, self.row_default
                        self.row = copy.copy(self.row_default)
                        self.row["input"] = stats_out
                        
                    else:
                        print "###############"
                        raw_input
                    
                
                elif strategy == "EFS":
                    self.row["no_emp_formula"] = len(records_f)
                    self.row["no_metabolite_names"] = 0
                    self.row["EFS"] = EFS_flag
                    self.row["strategy"] = MIDB_flag
                    self.row["emp_formula_parent"] = ""
                    self.row["no_emp_formulae_related_to_parent"] = ""
                    self.row["emp_formula_peak"] = len(records_f)
                    self.row["ion_form"] = ""
                    self.row["exact_mass"] = ""
                    self.row["exact_adduct_mass"] = ""
                    self.row["mass_error"] = ""
                    self.row["databases"] = "\t".join(map(str,[0]*len(db_names)))

                    #print self.row.values(), 9
                    
                    self.rows.append(self.row) ####################################################
                    #print 9, self.row_default
                    self.row = copy.copy(self.row_default)
                    self.row["input"] = stats_out
                    
        if strategy == "EFS":
            self.output.write("%s\tNo. empirical formulae\tNo. metabolite names\tEFS\tSPS/TM\tEmpirical formula (parent)\tNo. emperical formulae (parent related)\tEmpirical formula (peak)\tIon form\tTheoretical mass (neutral) (Da)\tTheoretical mass (Da)\tMass error (ppm)\t%s%s\n" % (self.peak_list.header.replace("\n",""), isotope_info_header, "\t".join(map(str,db_names))))
        else:
            self.output.write("%s\tNo. empirical formulae\tNo. metabolite names\tEFS\t%s\tEmpirical formula (parent)\tNo. emperical formulae (parent related)\tEmpirical formula (peak)\tIon form\tTheoretical mass (neutral) (Da)\tTheoretical mass (Da)\tMass error (ppm)\t%s%s\n" % (self.peak_list.header.replace("\n",""), strategy, isotope_info_header, "\t".join(map(str,db_names))))
        
        for row in self.rows:
            row["no_emp_formulae_related_to_parent"] = temp_parent_formulae.count(row["emp_formula_parent"])
            self.output.write("%s\n" % ("\t".join(map(str,row.values()))))

        self.output.close()
        print "%s" % (self.fn_output)
        print "------------------------------------------------------------------"
        print "Finished"
        print


def CombineOutputFiles(fn_POS_out, fn_NEG_out):
        
    inp_pos = open(fn_POS_out, "r")
    inp_neg = open(fn_NEG_out, "r")
    
    header_pos = inp_pos.readline()
    header_neg = inp_neg.readline()
    
    if header_pos == header_neg:
        
        if "SPS" in fn_POS_out:
            fn_out = "MI_Pack__output_SPS_POS_NEG_combined.txt"
        elif "TM" in fn_POS_out:
            fn_out = "MI_Pack__output_TM_POS_NEG_combined.txt"
        elif "EFS" in fn_POS_out:
            fn_out = "MI_Pack__output_EFS_POS_NEG_combined.txt"
        else:
            fn_out = "MI_Pack__output_POS_NEG_combined.txt"
        try:
            column = header_pos.split("\t").index("Empirical formula (parent)")
        except:
            return "Please, check format files."
        
        all_formulae = []
        rows = []

        for open_f in [inp_pos, inp_neg]:
            
            lines = open_f.readlines()
                
            for line in lines:
                                
                emp_formula = line.split("\t")[column]
                all_formulae.append(emp_formula)
        
            for line in lines:                
                line = line.split("\t")                
                emp_formula = line[column] 
                line[column+1] = str(all_formulae.count(emp_formula))
                line[0] = float(line[0])
                rows.append(line)

        rows = sorted(rows, key=operator.itemgetter(0, column))
        out = open(os.path.join(fn_POS_out.replace(os.path.basename(fn_POS_out), ""), fn_out), "w")
        out.write(header_pos)
        for row in rows:
            out.write("\t".join(map(str,row)))
        out.close()
        
        return
    else:
        return "Format of both files is different"

