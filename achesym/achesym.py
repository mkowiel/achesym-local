# -*- coding: utf-8 -*-
#
# Author: Marcin Kowiel, 2013-2019
# Email : mkowiel@ump.edu.pl
# Last modified: 25-09-2019 version 1.0.11
#
# 1.0.11 (13.02.2020):
#   fix bug with atom Mismatch
#   pep8 fixes
#   unification of test
#
# 1.0.10 (25.09.2019):
#   remove usage of mmtbx.utils.write_pdb_file
#   increase MAX_COMBINATIONS 4 times
#
# 1.0.9 (03.06.2018):
#   fix problem with memory error (for server with 2GB RAM)
#   validate input missing CRYST and SCALE
#
# 1.0.8 (08-10-2017):
#   fix normalizer for space group I41 (#80)
#   add __slots__ to GridPoint class
#   few performance improvements
#   add Exception if chain contact is not found


import math
import sys
import os
import datetime
import gc

from itertools import product, izip, combinations
from operator import attrgetter
from collections import OrderedDict, defaultdict

import iotbx.pdb
from iotbx import cif
from cctbx import sgtbx
from cctbx import xray
from cctbx.array_family import flex
from cctbx.sgtbx import rt_mx
from cctbx import adptbx
from scitbx.math import principal_axes_of_inertia
from scitbx import matrix
from mmtbx.tls import tools

ACHESYM_VERSION = '1.0.11'

DEG_AS_RAD = math.pi / 180.0
# 2000000 should fit 2GB of memory so 2000000*4 should fit into 8GB limit
MAX_COMBINATIONS = 2000000 * 4
EPS = 10e-15


class Normalizer(object):
    __slots__ = []

    @classmethod
    def factory(cls, space_group_number):  # pragma: no mccabe
        if space_group_number == 1:
            return NormalizerSG1()
        elif space_group_number in (3, 4):
            return NormalizerSG3()
        elif space_group_number == 5:
            return NormalizerSG5()
        elif space_group_number in (16, 17, 18, 19):
            return NormalizerSG16()
        elif space_group_number in (20, 21):
            return NormalizerSG20()
        elif space_group_number == 22:
            return NormalizerSG22()
        elif space_group_number in (23, 24):
            return NormalizerSG23()
        elif space_group_number in (75, 76, 77, 78):
            return NormalizerSG75()
        elif space_group_number == 79:
            return NormalizerSG79()
        elif space_group_number == 80:
            return NormalizerSG80()
        elif space_group_number in (89, 90, 91, 92, 93, 94, 95, 96):
            return NormalizerSG89()
        elif space_group_number == 97:
            return NormalizerSG97()
        elif space_group_number == 98:
            return NormalizerSG98()
        elif space_group_number in (143, 144, 145):
            return NormalizerSG143()
        elif space_group_number == 146:
            return NormalizerSG146()
        elif space_group_number in (149, 151, 153):
            return NormalizerSG149()
        elif space_group_number in (150, 152, 154):
            return NormalizerSG150()
        elif space_group_number == 155:
            return NormalizerSG155()
        elif space_group_number in (168, 169, 170, 171, 172, 173):
            return NormalizerSG168()
        elif space_group_number in (177, 178, 179, 180, 181, 182):
            return NormalizerSG177()
        elif space_group_number == 195:
            return NormalizerSG195()
        elif space_group_number == 196:
            return NormalizerSG196()
        elif space_group_number == 197:
            return NormalizerSG197()
        elif space_group_number == 198:
            return NormalizerSG198()
        elif space_group_number == 199:
            return NormalizerSG198()
        elif space_group_number in (207, 208):
            return NormalizerSG207()
        elif space_group_number == 209:
            return NormalizerSG209()
        elif space_group_number == 210:
            return NormalizerSG210()
        elif space_group_number == 211:
            return NormalizerSG211()
        elif space_group_number in (212, 213):
            return NormalizerSG212()
        elif space_group_number == 214:
            return NormalizerSG214()
        assert 0, "Unsupportet space group: No" + space_group_number

    def coset(self):
        return ("x, y, z",)

    def continuous_traslation(self):
        """ for example if z may be arbitrary we set it to 0.0 """
        return (0.0, 0.0, 0.0)

    def is_in_anticheshire_cell(self, x, y, z):
        """ x, y, z shoule be >= 0.0 and < 1.0 """
        return 0.0 <= x < 1.0 and 0.0 <= y < 1.0 and 0.0 <= z < 1.0

    def anticheshire_post_translation(self, x, y, z):
        """ for cubic space group it may be necessary to translate the point """
        return (0.0, 0.0, 0.0)


class NormalizerSG1(Normalizer):

    def continuous_traslation(self):
        return (1.0, 1.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return x == 0.0 and y == 0.0 and z == 0.0
        return x == 0.5 and y == 0.5 and z == 0.5


class NormalizerSG3(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y,z", "x, y, z+1/2", "x+1/2,y,z+1/2")

    def continuous_traslation(self):
        return (0.0, 1.0, 0.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x < 0.25 and y == 0.0 and 0.0 <= z < 0.5
        return 0.0 <= x < 0.25 and y == 0.5 and 0.0 <= z < 0.5


class NormalizerSG5(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2")

    def continuous_traslation(self):
        return (0.0, 1.0, 0.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x < 0.25 and y == 0.0 and 0.0 <= z < 0.5
        return 0.0 <= x < 0.25 and y == 0.5 and 0.0 <= z < 0.5


class NormalizerSG16(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y,z", "x,y+1/2,z", "x+1/2,y+1/2,z", "x, y, z+1/2", "x+1/2,y,z+1/2", "x,y+1/2,z+1/2",
                "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x < 0.25 and 0.0 <= y < 0.25 and 0.0 <= z < 0.5


class NormalizerSG20(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y,z", "x, y, z+1/2", "x+1/2,y,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x < 0.25 and 0.0 <= y < 0.25 and 0.0 <= z < 0.5


class NormalizerSG22(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/4,y+1/4,z+1/4", "x+1/2,y+1/2,z+1/2", "x+3/4,y+3/4,z+3/4")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x < 0.25 and 0.0 <= y < 0.25 and 0.0 <= z < 0.25


class NormalizerSG23(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y,z", "x,y+1/2,z", "x+1/2,y+1/2,z")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x < 0.25 and 0.0 <= y < 0.25 and 0.0 <= z < 0.5


class NormalizerSG75(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z", "y+1/2,x+1/2,-z", "y,x,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 0.25 and x-1 <= y < 0.5-x and z == 0.0
        return 0.0 <= x <= 0.25 and x <= y < 0.5 - x and z == 0.5


class NormalizerSG79(Normalizer):

    def coset(self):
        return ("x, y, z", "y,x,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 0.25 and x <= y <= 0.5-x and z == 0.0
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and z == 0.5


class NormalizerSG80(Normalizer):

    def coset(self):
        return ("x, y, z", "y,x,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 0.25 and x <= y <= 0.5-x and z == 0.0
        # return 0.0 <= x <= 0.25 and x <= y <= 0.5-x and z == 0.5
        return 0.0 <= x <= 0.25 and 0.5 - x <= y <= 0.5 + x and z == 0.5


class NormalizerSG89(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z", "x, y, z+1/2", "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and 0.0 <= z < 0.5


class NormalizerSG97(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and 0.0 <= z < 0.5


class NormalizerSG98(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and 0.5 - x <= y < 0.5 + x and 0.0 <= z < 0.5


class NormalizerSG143(Normalizer):

    def coset(self):
        return ("x, y, z", "x+2/3,y+1/3,z", "x+1/3,y+2/3,z", "-x,-y,z", "-x-2/3,-y-1/3,z", "-x-1/3,-y-2/3,z",
                "y,x,-z", "y+1/3,x+2/3,-z", "y+2/3,x+1/3,-z", "-y,-x,-z", "-y-1/3,-x-2/3,-z", "-y-2/3,-x-1/3,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 1.0/3.0 and 0.0 <= y <= x*0.5 and z == 0.0
        return 0.0 <= x <= 1.0 / 3.0 and 0.0 <= y <= x * 0.5 and z == 0.5


class NormalizerSG146(Normalizer):

    def coset(self):
        return ("x, y, z", "y,x,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 1.0/3.0 and 0.0 <= y <= x and z == 0.0
        return 0.0 <= x <= 1.0 / 3.0 and 0.0 <= y <= x and z == 0.5


class NormalizerSG149(Normalizer):

    def coset(self):
        return ("x, y, z", "x+2/3,y+1/3,z", "x+1/3,y+2/3,z", "x, y, z+1/2", "x+2/3,y+1/3,z+1/2", "x+1/3,y+2/3,z+1/2",
                "-x-2/3,-y-1/3,z", "-x-1/3,-y-2/3,z", "-x,-y,z+1/2", "-x-2/3,-y-1/3,z+1/2", "-x,-y,z",
                "-x-1/3,-y-2/3,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 1.0 / 3.0 and 0.0 <= y <= 0.5 * x and 0 <= z < 0.5


class NormalizerSG150(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2", "-x,-y,z+1/2", "-x,-y,z")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 2.0 / 3.0 and 0.0 <= y <= 0.5 * x and 2.0 * x - 1.0 <= y and 0 <= z < 0.5


class NormalizerSG155(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 1.0 / 3.0 and 0.0 <= y < 1.0 / 3.0 and 0 <= z <= 0.25


class NormalizerSG168(Normalizer):

    def coset(self):
        return ("x, y, z", "y,x,-z")

    def continuous_traslation(self):
        return (0.0, 0.0, 1.0)

    def is_in_anticheshire_cell(self, x, y, z):
        # return 0.0 <= x <= 2.0/3.0 and 0.0 <= y <= 0.5*x and 2.0*x-1.0 <= y and z == 0.0
        return 0.0 <= x <= 2.0 / 3.0 and 0.0 <= y <= 0.5 * x and 2.0 * x - 1.0 <= y and z == 0.5


class NormalizerSG177(Normalizer):

    def coset(self):
        return ("x, y, z", "x, y, z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 2.0 / 3.0 and 0.0 <= y <= 0.5 * x and 2.0 * x - 1.0 <= y and 0 <= z < 0.5


class NormalizerSG195(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2", "y,x,-z", "y+1/2,x+1/2,-z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and x <= z <= 0.5 - x


class NormalizerSG196(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/4,y+1/4,z+1/4", "x+1/2,y+1/2,z+1/2", "x+3/4,y+3/4,z+3/4",
                "y,x,-z", "y+1/4,x+1/4,-z-1/4", "y+1/2,x+1/2,-z-1/2", "y+3/4,x+3/4,-z-3/4")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.125 and x <= y <= 0.25 - x and x <= z <= 0.25 - x


class NormalizerSG197(Normalizer):

    def coset(self):
        return ("x, y, z", "y,x,-z")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and x <= z <= 0.5 - x


class NormalizerSG198(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2", "-y+1/4,-x+1/4,-z+1/4", "-y+3/4,-x+3/4,-z+3/4")

    def is_in_anticheshire_cell(self, x, y, z):
        x = x - 1.0 if x > 0.5 else x
        y = y - 1.0 if y > 0.5 else y
        z = z - 1.0 if z > 0.5 else z
        # return -0.125 <= y <= 0.125 and y < z <= 0.25+y and x <= 0.125 and y-z-0.125 <= x <= z
        return -0.375 < x < 0.125 and -0.125 < y < 0.125 and max(x, y, y - x - 0.125) < z < y + 0.25

    def anticheshire_post_translation(self, x, y, z):
        a = -1.0 if x > 0.5 else 0.0
        b = -1.0 if y > 0.5 else 0.0
        c = -1.0 if z > 0.5 else 0.0
        return (a, b, c)


class NormalizerSG199(Normalizer):

    def coset(self):
        return ("x, y, z", "-y+1/4,-x+1/4,-z+1/4")

    def is_in_anticheshire_cell(self, x, y, z):
        x = x - 1.0 if x > 0.5 else x
        y = y - 1.0 if y > 0.5 else y
        z = z - 1.0 if z > 0.5 else z
        # return -0.125 <= y <= 0.125 and y < z <= 0.25+y and x <= 0.125 and y-z-0.125 <= x <= z
        return -0.375 < x < 0.125 and -0.125 < y < 0.125 and max(x, y, y - x - 0.125) < z < y + 0.25

    def anticheshire_post_translation(self, x, y, z):
        a = -1.0 if x > 0.5 else 0.0
        b = -1.0 if y > 0.5 else 0.0
        c = -1.0 if z > 0.5 else 0.0
        return (a, b, c)


class NormalizerSG207(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and x <= z <= 0.5 - x


class NormalizerSG209(Normalizer):
    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.5 and x <= y <= 0.5 - x and 0 <= z <= 0.25


class NormalizerSG210(Normalizer):
    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and 0 <= y <= 0.25 and x <= z <= 0.25 - x and y <= z <= 0.25 - y


class NormalizerSG211(Normalizer):

    def is_in_anticheshire_cell(self, x, y, z):
        return 0.0 <= x <= 0.25 and x <= y <= 0.5 - x and x <= z <= 0.5 - x


class NormalizerSG212(Normalizer):

    def coset(self):
        return ("x, y, z", "x+1/2,y+1/2,z+1/2")

    def is_in_anticheshire_cell(self, x, y, z):
        x = x - 1.0 if x > 0.5 else x
        y = y - 1.0 if y > 0.5 else y
        z = z - 1.0 if z > 0.5 else z
        # return -0.125 <= y <= 0.125 and y < z <= 0.25+y and x <= 0.125 and y-z-0.125 <= x <= z
        return -0.375 < x < 0.125 and -0.125 < y < 0.125 and max(x, y, y - x - 0.125) < z < y + 0.25

    def anticheshire_post_translation(self, x, y, z):
        a = -1.0 if x > 0.5 else 0.0
        b = -1.0 if y > 0.5 else 0.0
        c = -1.0 if z > 0.5 else 0.0
        return (a, b, c)


class NormalizerSG214(Normalizer):

    def is_in_anticheshire_cell(self, x, y, z):
        x = x - 1.0 if x > 0.5 else x
        y = y - 1.0 if y > 0.5 else y
        z = z - 1.0 if z > 0.5 else z
        # return -0.125 <= y <= 0.125 and y < z <= 0.25+y and x <= 0.125 and y-z-0.125 <= x <= z
        return -0.375 < x < 0.125 and -0.125 < y < 0.125 and max(x, y, y - x - 0.125) < z < y + 0.25

    def anticheshire_post_translation(self, x, y, z):
        a = -1.0 if x > 0.5 else 0.0
        b = -1.0 if y > 0.5 else 0.0
        c = -1.0 if z > 0.5 else 0.0
        return (a, b, c)


# =================
def fix_zero(x):
    return 0.0 if abs(x) < EPS else x


def fix_zero_3(vec3):
    return (fix_zero(vec3[0]), fix_zero(vec3[1]), fix_zero(vec3[2]))


def translate_in_fractional(xray_structure, t):
    sites = xray_structure.sites_frac()
    sites = sites + t
    xray_structure.set_sites_frac(sites)


def fround(x, digits):
    pos = 10 ** int(math.floor(math.log10(x) - digits + 1))
    return round(x / pos) * pos


def find_best_tanslation(p1, p2, unit_cell, translations=(0.0, -1.0, 1.0)):
    """ p1 and p2 should be within 0.0 x 1.0 box"""
    # assert(p1[0] >= 0.0 and p1[0] <= 1.0)
    # assert(p1[1] >= 0.0 and p1[1] <= 1.0)
    # assert(p1[2] >= 0.0 and p1[2] <= 1.0)

    assert (p2[0] >= 0.0 and p2[0] <= 1.0)
    assert (p2[1] >= 0.0 and p2[1] <= 1.0)
    assert (p2[2] >= 0.0 and p2[2] <= 1.0)

    d = unit_cell.orthogonalize((p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]))
    min_len_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
    min_xyz = (0.0, 0.0, 0.0)
    for x in translations:
        for y in translations:
            for z in translations:
                d = unit_cell.orthogonalize((p2[0] - p1[0] + x, p2[1] - p1[1] + y, p2[2] - p1[2] + z))
                len_sq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2]
                if (len_sq < min_len_sq):
                    min_len_sq = len_sq
                    min_xyz = (x, y, z)
    return min_xyz


class TLSTranformator(object):
    def __init__(self, pdb):
        self.pdb_filename = pdb
        self.pdb_inp = iotbx.pdb.input(file_name=self.pdb_filename)
        self.pdb_hierarchy = self.pdb_inp.construct_hierarchy()

    def read_tls(self):
        remark_3_records = self.pdb_inp.extract_remark_iii_records(iii=3)
        input_tls = tools.tls_from_pdb_inp(remark_3_records, self.pdb_hierarchy)
        selection_strings = [item.selection_string for item in input_tls.tls_params]
        tlsos = [tools.tlso(item.t, item.l, item.s, origin=item.origin) for item in input_tls.tls_params]
        return tlsos, selection_strings

    @classmethod
    def get_mat_form_vec(cls, xyz):
        px, py, pz = xyz
        return matrix.sqr((0.0, pz, -py,
                           -pz, 0.0, px,
                           py, -px, 0.0))

    @classmethod
    def get_TLSSt(cls, tlso):
        T = matrix.sym(sym_mat3=tlso.t)
        L = DEG_AS_RAD * DEG_AS_RAD * matrix.sym(sym_mat3=tlso.l)
        S = DEG_AS_RAD * matrix.sqr(tlso.s)
        S = S + matrix.sym(sym_mat3=(0.0, 0.0, -S[0] - S[4], 0.0, 0.0, 0.0))
        St = S.transpose()
        return (T, L, S, St)

    @classmethod
    def get_RRtPPt(cls, rot, tra):
        R = matrix.sqr(rot)
        Rt = R.transpose()
        P = cls.get_mat_form_vec(tra)
        Pt = P.transpose()
        return (R, Rt, P, Pt)

    @classmethod
    def get_RMRt(cls, rot, M, symetric):
        R = matrix.sqr(rot)
        Rt = R.transpose()
        if symetric:
            M = matrix.sym(sym_mat3=M)
            return matrix.sym(sym_mat3=matrix.sqr(R * M * Rt).as_sym_mat3())
        M = matrix.sqr(M)
        return matrix.sqr(R * M * Rt)

    @classmethod
    def transform_t(cls, tlso, rot, tra):
        T, L_rad, S_rad, St_rad = cls.get_TLSSt(tlso)
        R, Rt, P, Pt = cls.get_RRtPPt(rot, tra)
        # R(T + (PLP' + PS + S'P' ))R'
        return matrix.sym(
            sym_mat3=matrix.sqr(R * (T + (P * L_rad * Pt) + (P * S_rad) + (St_rad * Pt)) * Rt).as_sym_mat3())

    @classmethod
    def transform_l(cls, tlso, rot, tra):
        L = DEG_AS_RAD * DEG_AS_RAD * matrix.sym(sym_mat3=tlso.l)
        R = matrix.sqr(rot)
        Rt = R.transpose()
        # RLR'
        to_deg = 1. / DEG_AS_RAD
        to_deg = to_deg * to_deg
        return to_deg * matrix.sym(sym_mat3=matrix.sqr(R * L * Rt).as_sym_mat3())

    @classmethod
    def transform_s(cls, tlso, rot, tra):
        T, L_rad, S_rad, St_rad = cls.get_TLSSt(tlso)
        R, Rt, P, Pt = cls.get_RRtPPt(rot, tra)
        # R(S + LP')R' = R(S - LP)R'
        to_deg = 1. / DEG_AS_RAD
        return to_deg * R * (S_rad + L_rad * Pt) * Rt

    @classmethod
    def transform_o(cls, tlso, rot, tra):
        xyz = flex.vec3_double(1, tlso.origin)
        xyz = rot * xyz + tra
        return xyz[0]

    @classmethod
    def transform_t_alt(cls, tlso, rot, tra):
        return cls.get_RMRt(rot, tlso.t, True)

    @classmethod
    def transform_l_alt(cls, tlso, rot, tra):
        return cls.get_RMRt(rot, tlso.l, True)

    @classmethod
    def transform_s_alt(cls, tlso, rot, tra):
        return cls.get_RMRt(rot, tlso.s, False)

    @classmethod
    def transform_o_alt(cls, tlso, rot, tra):
        return cls.transform_o(tlso, rot, tra)

    @classmethod
    def calc_u(cls, tlso, r):
        T, L_rad, S_rad, St_rad = cls.get_TLSSt(tlso)
        A = cls.get_mat_form_vec(r) - cls.get_mat_form_vec(tlso.origin)
        At = A.transpose()
        U = T + (A * L_rad * At) + (A * S_rad) + (St_rad * At)
        return matrix.sqr(U)  # matrix.sym(sym_mat3=matrix.sqr(U).as_sym_mat3())

    @classmethod
    def calc_U(cls, tlso, r, ueq):
        ucalc = cls.calc_u(tlso, r)
        trucalc = adptbx.b_as_u(ueq) - ucalc.trace() / 3.0
        u = matrix.sym(sym_mat3=[trucalc, trucalc, trucalc, 0, 0, 0])
        return ucalc + u

    @classmethod
    def is_unitary(cls, mat):
        R = matrix.sqr(mat)
        I_ = R * R.transpose()
        if abs(sum(I_ - matrix.sqr([1, 0, 0, 0, 1, 0, 0, 0, 1]))) < 0.001:
            return True
        return False

    def transform(self, chains, rotations_in_cart, translations_in_cart):
        tlsos, selection_strings = self.read_tls()
        # tools.show_tls(tlsos = tlsos)
        out_tlsos = []
        for _tlso, sel_str in zip(tlsos, selection_strings):
            tlso = _tlso
            for chain, rot, tra in zip(chains, rotations_in_cart, translations_in_cart):
                if self.is_unitary(rot) is False:
                    raise Exception("Rotaton matrix should be unitary matrix")
                if chain == 'ALL' or chain in sel_str:
                    t_mat = list(self.transform_t_alt(tlso, rot, tra))
                    t_mat = [t_mat[0], t_mat[4], t_mat[8], t_mat[1], t_mat[2], t_mat[5]]
                    l_mat = list(self.transform_l_alt(tlso, rot, tra))
                    l_mat = [l_mat[0], l_mat[4], l_mat[8], l_mat[1], l_mat[2], l_mat[5]]
                    s_mat = list(self.transform_s_alt(tlso, rot, tra))
                    origin = list(self.transform_o_alt(tlso, rot, tra))
                    tlso = tools.tlso(t_mat, l_mat, s_mat, origin=origin)
            out_tlsos.append(tlso)
        return out_tlsos, selection_strings


class GridPoint(object):
    __slots__ = ['chain_tag', 'sym_op', 'sym_op_id', 'tx', 'ty', 'tz', 'prev_tag']

    def __init__(self, chain_tag, sym_op, tx=0, ty=0, tz=0, prev_tag='', sym_op_id=None):
        self.chain_tag = chain_tag
        self.sym_op = sym_op
        self.sym_op_id = sym_op_id
        self.tx = tx
        self.ty = ty
        self.tz = tz
        self.prev_tag = prev_tag

    def __unicode__(self):
        # return u'%s: %s + %s,%s,%s' % (
        #   unicode(self.chain_tag), unicode(self.sym_op),
        #   int(self.tx), int(self.ty), int(self.tz)
        #   )
        sym = 1 + int(self.sym_op_id) if self.sym_op_id is not None else unicode(self.sym_op)
        return u'%s:%s%s%s%s' % (unicode(self.chain_tag), sym, 5 + int(self.tx), 5 + int(self.ty), 5 + int(self.tz))

    __repr__ = __unicode__

    def pretty_print(self):
        return u'%s[%s]: %s + %s,%s,%s' % (
            unicode(self.chain_tag), unicode(self.prev_tag), unicode(self.sym_op),
            int(self.tx), int(self.ty), int(self.tz)
        )

    def _tx_lt(self, other):
        return int(self.tx) < int(other.tx)

    def _ty_equal(self, other):
        return int(self.ty) < int(other.ty)

    def _tz_equal(self, other):
        return int(self.tz) < int(other.tz)

    def _t_lt(self, other):
        return (
            int(self.tx) < int(other.tx) or (
                int(self.tx) == int(other.tx) and (
                    int(self.ty) < int(other.ty) or (
                        int(self.ty) == int(other.ty) and int(self.tz) < int(other.tz)
                    )
                )
            )
        )

    def _t_equal(self, other):
        return (int(self.tx) == int(other.tx)) and int(self.ty) == int(other.ty) and int(self.tz) == int(other.tz)

    def __eq__(self, other):
        return self.chain_tag == other.chain_tag and str(self.sym_op) == str(other.sym_op) and self._t_equal(other)

    def __ne__(self, other):
        return not self.__eq__(other)

    def __lt__(self, other):
        return (
            (self.chain_tag < other.chain_tag) or (
                self.chain_tag == other.chain_tag and (
                    (str(self.sym_op) < str(other.sym_op)) or
                    (str(self.sym_op) == str(other.sym_op) and self._t_lt(other))
                )
            )
        )

    def __le__(self, other):
        return self.__lt__(other) or self.__eq__(other)

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __hash__(self):
        return hash(self.__repr__())


class GridCountBase(OrderedDict):

    def __init__(self, grid, log_stream):
        super(GridCountBase, self).__init__()
        self.log_stream = log_stream
        self.d = grid.grid_spacing
        group = defaultdict(int)
        for points in grid.itervalues():
            for tag in self.get_tags(points):
                group[tag] += 1
        for tag, count in self.sort_group(group):
            self[tag] = count

    def __str__(self):
        txt = ""
        for tag, count in self.iteritems():
            txt += '%s: ~%.1f [A^3]\n' % (', '.join((str(point) for point in tag)), count * self.d * self.d * self.d)
        return txt

    def get_tags(self, point_list):
        return [tuple(sorted(point_list))]

    def sort_group(self, group):
        a = sorted(group.iteritems())
        b = sorted(a, key=lambda x: x[1], reverse=True)
        return sorted(b, key=lambda x: len(x[0]), reverse=True)

    def reinit(self):
        group = {}
        for tag, count in self.iteritems():
            tags = self.get_tags(tag)
            for tag in tags:
                group[tag] = group.setdefault(tag, 0) + count
        self.clear()
        for tag, count in self.sort_group(group):
            self[tag] = count


class GridCountDetails(GridCountBase):

    def __init__(self, grid, log_stream=sys.stdout):
        super(GridCountDetails, self).__init__(grid, log_stream)

    def __str__(self):
        txt = ""
        for tag, count in self.iteritems():
            txt += '%s: ~%.1f [A^3]\n' % (', '.join((str(point) for point in tag)), count * self.d * self.d * self.d)
        return txt

    def get_tags(self, point_list):
        point_for_chain_id = defaultdict(list)
        for p in point_list:
            point_for_chain_id[p.chain_tag].append(p)

        unique_chain_point_list = []
        for chain_id in sorted(point_for_chain_id.iterkeys()):
            unique_chain_point_list.append(sorted(point_for_chain_id[chain_id]))
        return self.subtract_first_t(list(product(*unique_chain_point_list)))

    def subtract_first_t(self, tags):
        new_tags = []
        for tag in tags:
            new_tag = []
            if len(tag) > 0:
                first_tag = tag[0]
                tx = first_tag.tx
                ty = first_tag.ty
                tz = first_tag.tz
                for point in tag:
                    new_tag.append(GridPoint(
                        point.chain_tag,
                        point.sym_op,
                        point.tx - tx,
                        point.ty - ty,
                        point.tz - tz,
                        point.prev_tag,
                        point.sym_op_id
                    ))
            new_tags.append(tuple(new_tag))
        return new_tags

    def collapse(self, chain_ids, function='sum', found_tags=None):  # pragma: no mccabe
        aggregation_func = sum
        if function == 'max' or function == 'Max':
            aggregation_func = max
        translations = {}
        out_translations = {}
        p_first = chain_ids[0]
        for p in found_tags:
            out_translations[p.chain_tag] = (-p.tx, -p.ty, -p.tz)
        for p in chain_ids:
            translations[(p.chain_tag, p.sym_op)] = (p_first.tx - p.tx, p_first.ty - p.ty, p_first.tz - p.tz)
        chain_ids = [GridPoint(p.chain_tag, p.sym_op, 0, 0, 0, '', p.sym_op_id) for p in chain_ids]
        chain_tag_sorted = sorted(list(set(chain_ids)))
        ids_set = set([p.chain_tag for p in chain_tag_sorted])
        new_chain_tag = ''.join([p.chain_tag for p in chain_tag_sorted])
        new_data = {}
        for tag, count in self.iteritems():
            to_translate = {}
            for point in tag:
                if GridPoint(point.chain_tag, point.sym_op, 0, 0, 0, '', point.sym_op_id) in chain_tag_sorted:
                    dtx = translations[(point.chain_tag, point.sym_op)][0] + point.tx
                    dty = translations[(point.chain_tag, point.sym_op)][1] + point.ty
                    dtz = translations[(point.chain_tag, point.sym_op)][2] + point.tz
                    if (dtx, dty, dtz) not in to_translate:
                        to_translate[(dtx, dty, dtz)] = point
            if len(to_translate) == 0:
                to_translate[(0, 0, 0)] = None
            for dt, dt_point in to_translate.iteritems():
                new_tag = []
                new_data_part = {}
                for point in tag:
                    if (
                        dt_point is not None and
                        point.chain_tag == dt_point.chain_tag and
                        point.sym_op == dt_point.sym_op
                    ):
                        reference = point.chain_tag if point.prev_tag == '' else point.prev_tag
                        tr_tag_x, tr_tag_y, tr_tag_z = out_translations[reference]
                        prev = found_tags[0].chain_tag
                        new_tag.append(GridPoint(new_chain_tag, point.sym_op, point.tx + tr_tag_x, point.ty + tr_tag_y,
                                                 point.tz + tr_tag_z, prev, point.sym_op_id))
                    elif (point.chain_tag in ids_set and point.prev_tag != ''):
                        tr_tag_x, tr_tag_y, tr_tag_z = out_translations[point.prev_tag]
                        prev = found_tags[0].chain_tag
                        new_tag.append(GridPoint(new_chain_tag, point.sym_op, point.tx + tr_tag_x, point.ty + tr_tag_y,
                                                 point.tz + tr_tag_z, prev, point.sym_op_id))
                    elif (point.chain_tag in ids_set):
                        pass
                    else:
                        new_tag.append(
                            GridPoint(point.chain_tag, point.sym_op, point.tx, point.ty, point.tz, point.prev_tag,
                                      point.sym_op_id))
                new_data_part[tuple(new_tag)] = count
            for temp_tag, temp_count in new_data_part.iteritems():
                new_data[temp_tag] = aggregation_func((new_data.setdefault(temp_tag, 0), temp_count))
        self.clear()
        for tag, count in new_data.iteritems():
            self[tag] = count
        return new_chain_tag

    def find_tags(self, tags):
        tags = set(tags)
        for key, count in self.iteritems():
            unique_tags = set((k.chain_tag for k in key))
            if len(key) == len(tags) and tags.issubset(unique_tags):
                return key
        # if not found try without len
        for key, count in self.iteritems():
            unique_tags = set((k.chain_tag for k in key))
            if tags.issubset(unique_tags):
                return key

    def find_any_tag_and_symmetry(self, tag_point_dict):
        points = set(tag_point_dict.values())
        points_tag = set((point.chain_tag for point in points))
        for key in self.keys():
            point_set = set(key)
            point_set_tag = set((point.chain_tag for point in point_set))
            if not point_set_tag.issubset(points_tag) and points.isdisjoint(point_set) is False:
                return key

    def get_symmetries(self, tag):
        out = OrderedDict()
        for point in tag:
            if point.chain_tag not in out:
                out[point.chain_tag] = point
        return out

    def get_step(self, tag, tag_found):
        tag_int = set(list(tag)).intersection(set(tag_found))
        tag_dif = set(list(tag)).difference(set(tag_found))
        steps = [sorted(list(tag_int))]
        for key in sorted(tag_dif, key=attrgetter('chain_tag')):
            new = list(steps[-1]) if len(steps) > 0 else []
            new.append(key)
            steps.append(list(new))
        result = [step[0:-1] for step in steps[1:]]
        return result

    def get_transformations(self, grid_count_unique):  # pragma: no mccabe
        chains = set()
        for tags in grid_count_unique.iterkeys():
            chains = chains.union(set(tags))
        # print 'Chains:', chains

        len_chain_last = len(chains)
        tags = OrderedDict()
        steps = []
        starting_tag = grid_count_unique.keys()[0]
        tag = self.find_tags(starting_tag)
        self.pop(tag)
        tag_dict = self.get_symmetries(tag)
        steps.extend(self.get_step(tag, tags.values()))
        tags.update(tag_dict)
        chains = chains.difference(set(tags.keys()))
        while len(chains) > 0 and len_chain_last != len(chains):
            tag = self.find_any_tag_and_symmetry(tags)
            self.pop(tag)
            tag_dict = self.get_symmetries(tag)
            # print tag, row, tag_dict
            steps.extend(self.get_step(tag, tags.values()))
            tags.update(tag_dict)
            len_chain_last = len(chains)
            chains = chains.difference(set(tags.keys()))
        tags_dict = OrderedDict()
        for key, val in tags.iteritems():
            tags_dict[key] = val.sym_op
        steps_out = []
        for step_list in steps:
            tmp = [step.chain_tag for step in step_list]
            steps_out.append(tmp)
        # print 'transformations:', tags_dict, steps_out
        return tags_dict, steps_out

    def find_tag_t(self, tag):
        tag = set(tag)
        for key, count in self.iteritems():
            unique_tag = set(key)
            if tag.issubset(unique_tag):
                return key

    def find_chain_tag_t(self, chain_tag):
        for tag in self.iterkeys():
            if chain_tag in [p.chain_tag for p in tag]:
                return tag
        return None

    def find_t(self, tag_list, chain_tag):
        if chain_tag is None or tag_list is None:
            return 0, 0, 0, None
        for p in tag_list:
            if p.chain_tag == chain_tag:
                return p.tx, p.ty, p.tz, p.prev_tag
        return 0, 0, 0, None

    def get_transformations_t(self, function='sum'):  # pragma: no mccabe
        tag_out = []
        chains = set()
        for tags in self.iterkeys():
            for point in tags:
                chains.add(point.chain_tag)

        starting_tag = self.keys()[0]
        tag_out.extend(list(starting_tag))
        # print >> self.log_stream, 'TAGS', tag_out
        # print >> self.log_stream, 'TABLE'
        # print >> self.log_stream, self
        while len(chains) > 0:
            tag = self.find_tag_t(starting_tag)
            self.pop(tag)
            tags_found = set([p.chain_tag for p in tag])
            chains = chains.difference(tags_found)
            collapsed_chain_tag = self.collapse(tag, function, tag_out)
            self.reinit()
            # print >> self.log_stream, 'TAGS', tag_out
            # print >> self.log_stream, 'TABLE'
            # print >> self.log_stream, self
            starting_tag = self.find_chain_tag_t(collapsed_chain_tag)
            if starting_tag is not None:
                tx, ty, tz, prev_tag = self.find_t(starting_tag, collapsed_chain_tag)
                tx2, ty2, tz2 = 0, 0, 0
                tag_out.extend([GridPoint(p.chain_tag, p.sym_op, p.tx + tx2 - tx, p.ty + ty2 - ty, p.tz + tz2 - tz,
                                          p.prev_tag, p.sym_op_id) for p in starting_tag if
                                p.chain_tag != collapsed_chain_tag])
            elif len(chains) > 0:
                raise Exception(
                    'Contact between chain(s) %s and other chains not found. '
                    'Try to increase radius in advanced options.' % str(collapsed_chain_tag)
                )

        return tag_out


class GridCountPairSimple(GridCountDetails):

    def __init__(self, grid, log_stream=sys.stdout):
        super(GridCountPairSimple, self).__init__(grid, log_stream)

    def get_tags(self, point_list):
        """
        return list with tags based on grid points in the point_list list
        """
        unique_tags = super(GridCountPairSimple, self).get_tags(point_list)
        new_tag = []
        for unique_tag in unique_tags:
            new_tag.extend(list(combinations(unique_tag, r=2)))
        return self.subtract_first_t(new_tag)


class GridCountPair(GridCountDetails):

    def __init__(self, grid, log_stream=sys.stdout):
        super(GridCountPair, self).__init__(grid, log_stream)

    def get_tags(self, point_list):
        """
        return list with tags based on grid points in the point_list list
        """
        unique_tags = super(GridCountPair, self).get_tags(point_list)
        new_tag = []
        for unique_tag in unique_tags:
            new_tag.extend(list(combinations(unique_tag, r=2)))
        return self.subtract_first_t(new_tag)

    def find_tag_t(self, tag):
        return super(GridCountPair, self).find_tag_t(self.keys()[0])

    def get_transformations_t(self, function='sum'):  # pragma: no mccabe
        tag_out = []
        assemblies = {}

        starting_tag = self.keys()[0]
        tag_out.extend(list(starting_tag))

        chains = set()
        for tags in self.iterkeys():
            for point in tags:
                chains.add(point.chain_tag)
        # print >> self.log_stream, 'TAGS', tag_out
        # print >> self.log_stream, 'TABLE'
        # print >> self.log_stream, self
        kk = 56
        # while len(assemblies[collapsed_chain_tag]) == len(chains):
        while len(self) > 0 and kk > 0:
            kk -= 1
            tag = self.find_tag_t(starting_tag)
            self.pop(tag)

            tag0 = tag[0].chain_tag
            tag1 = tag[1].chain_tag

            tag_out = list(assemblies[tag0]) if tag0 in assemblies else [tag[0]]
            if tag1 in assemblies:
                tag_out.extend([GridPoint(p.chain_tag, p.sym_op, p.tx + tag[1].tx, p.ty + tag[1].ty, p.tz + tag[1].tz,
                                          p.prev_tag, p.sym_op_id) for p in assemblies[tag1]])  # assemblies[tag1]
            else:
                tag_out.append(tag[1])
            collapsed_chain_tag = self.collapse(tag, function, tag_out)

            assemblies_cout = 0
            if tag0 in assemblies:
                assemblies_cout += 1
            if tag1 in assemblies:
                assemblies_cout += 1

            if assemblies_cout == 0:
                # powstaje nowy aglomerat
                assemblies[collapsed_chain_tag] = list(tag)
            elif assemblies_cout == 1:
                # dolacza element do istnejacego assembly
                assemblies[collapsed_chain_tag] = []
                if tag0 in assemblies:
                    assemblies[collapsed_chain_tag].extend(assemblies[tag[0].chain_tag])
                else:
                    assemblies[collapsed_chain_tag].extend([tag[0]])
                if tag1 in assemblies:
                    assemblies[collapsed_chain_tag].extend([GridPoint(p.chain_tag, p.sym_op, p.tx + tag[1].tx,
                                                                      p.ty + tag[1].ty, p.tz + tag[1].tz, p.prev_tag,
                                                                      p.sym_op_id) for p in assemblies[tag1]])
                else:
                    assemblies[collapsed_chain_tag].extend([tag[1]])
            else:
                assemblies[collapsed_chain_tag] = []
                assemblies[collapsed_chain_tag].extend(assemblies[tag0])
                assemblies[collapsed_chain_tag].extend([GridPoint(p.chain_tag, p.sym_op, p.tx + tag[1].tx,
                                                                  p.ty + tag[1].ty, p.tz + tag[1].tz, p.prev_tag,
                                                                  p.sym_op_id) for p in assemblies[tag1]])
                # złacz assemblies
            # print >> self.log_stream, 'TAGS', assemblies
            self.reinit()
            # print >> self.log_stream, 'TABLE'
            # print >> self.log_stream, self
            starting_tag = self.find_chain_tag_t(collapsed_chain_tag)

        return assemblies[collapsed_chain_tag]


class GridCountUnique(GridCountBase):

    def __init__(self, grid, log_stream):
        super(GridCountUnique, self).__init__(grid, log_stream)

    def get_tags(self, point_list):
        return [tuple(sorted(set((p.chain_tag for p in point_list))))]


class Grid(dict):

    def __init__(self, grid_spacing=0.7, log_stream=None):
        self.grid_spacing = grid_spacing
        self.log_stream = log_stream

    def grid_index(self, x):
        return int(math.floor(x / self.grid_spacing))

    def coord_from_grid_index(self, n):
        return n * self.grid_spacing

    def round_coord_to_grid(self, x):
        n = self.grid_index(x)
        return self.coord_from_grid_index(n)

    def cartesian_to_index(self, xyz):
        return tuple((self.grid_index(cart) for cart in xyz))

    def index_to_cartesian(self, nop):
        return tuple((self.coord_from_grid_index(index) for index in nop))

    def round_point_to_grid(self, xyz_c):
        return tuple((self.round_coord_to_grid(coord) for coord in xyz_c))

    def is_inside_unit_cell(self, grid_point, uc):
        xyz = uc.fractionalize(site_cart=grid_point)
        if (-EPS <= xyz[0] < 1.0 + EPS and -EPS <= xyz[1] < 1.0 + EPS and -EPS <= xyz[2] < 1.0 + EPS):
            return True
        return False

    def index_inside_unit_cell(self, grid_point, uc):
        xyz_f_in = uc.fractionalize(site_cart=grid_point)
        xyz_f = sgtbx.fractional_mod_positive(xyz_f_in)
        xyz_c = uc.orthogonalize(xyz_f)
        # make sure that index is not too low
        # add 1/4 of grid spacing
        # is save sine xyz_c is close to grid node
        shift = 0.25 * self.grid_spacing
        xyz_c = tuple((x + shift for x in xyz_c))
        nop = self.cartesian_to_index(xyz_c)
        # we could return here
        # but since the grid_spacing is not multiple of unit cell
        # there are problems with the border points
        # the index in fractional minus fractional(0.0, 0,0, 1.0)
        # may not be index point
        # and then rounded may fall outside of the unit cell region
        df = (xyz_f[0] - xyz_f_in[0], xyz_f[1] - xyz_f_in[1], xyz_f[2] - xyz_f_in[2])
        n_list = (0,) if int(math.floor(df[0] + EPS)) == 0 else (0, -1, 1)
        o_list = (0,) if int(math.floor(df[1] + EPS)) == 0 else (0, -1, 1)
        p_list = (0,) if int(math.floor(df[2] + EPS)) == 0 else (0, -1, 1)
        for n in n_list:
            for o in o_list:
                for p in p_list:
                    ind = (nop[0] + n, nop[1] + o, nop[2] + p)
                    xyz_c = self.index_to_cartesian(ind)
                    if self.is_inside_unit_cell(xyz_c, uc):
                        return ind
        for n in (0, -1, 1):
            for o in (0, -1, 1):
                for p in (0, -1, 1):
                    ind = (nop[0] + n, nop[1] + o, nop[2] + p)
                    xyz_c = self.index_to_cartesian(ind)
                    if self.is_inside_unit_cell(xyz_c, uc):
                        return ind
        raise Exception("index inside unit cell not found")

    # finds and initialize the grid indexes that
    # lie inside unit cell with empty list
    def create(self, uc):
        return None
        """
        n = []
        o = []
        p = []
        for x in (0.0, 1.0):
            for y in (0.0, 1.0):
                for z in (0.0, 1.0):
                    xyz_c = uc.orthogonalize((x, y, z))
                    i_n, i_o, i_p = self.cartesian_to_index(xyz_c)
                    n.append(i_n)
                    o.append(i_o)
                    p.append(i_p)

        X = max(n)+2-min(n)
        Y = max(o)+2-min(o)
        Z = max(p)+2-min(p)
        print 'X', X, 'Y', Y, 'Z', Z
        N =  X*Y*Z
        print 'ALL',
        if N > 14000000:
            raise Exception("Unit cell too big, maximal number of grid points exceeded")
        # +2 to make sure we don't ommit the last points
        for n_i in range(min(n), max(n)+2):
            for o_i in range(min(o), max(o)+2):
                for p_i in range(min(p), max(p)+2):
                    nop = (n_i, o_i, p_i)
                    xyz = self.index_to_cartesian(nop)
                    if self.is_inside_unit_cell(xyz, uc):
                        self.setdefault(nop, [])
        """

    def add(self, indexes, grid_point, uc):
        for ind in indexes:
            if ind not in self:
                xyz = self.index_to_cartesian(ind)
                if self.is_inside_unit_cell(xyz, uc):
                    self.setdefault(ind, [grid_point])
            elif (grid_point not in self[ind]):
                self[ind].append(grid_point)
        if len(self) > 8000000:
            raise Exception(
                "Memory error, maximal number of grid points exceeded, please increase grid size in advanced options")


class Packer(object):

    def __init__(self, in_pdb='in.pdb', chain_grouping=None, r=2.3, grid_spacing=0.7, function='sum', log_stream=None):
        self.in_pdb = in_pdb
        self.r = r
        self.grid_spacing = grid_spacing
        self.chain_grouping = chain_grouping
        self.function = function
        # open log
        if log_stream is None:
            self.log_stream = sys.stdout
        elif (type(log_stream) == str or type(log_stream) == unicode):
            self.log_stream = open(log_stream, 'w')
        else:
            self.log_stream = log_stream

    def is_in_sphere(self, p, center, r):
        x = p[0] - center[0]
        y = p[1] - center[1]
        z = p[2] - center[2]
        return x * x + y * y + z * z - r * r < 0.0

    # returns indexes that are on the grid that are inside sphere
    # with center in xyz_c point and self.r radius
    def grid_indexes(self, grid, xyz_c, uc):
        x, y, z = xyz_c

        x_min = x - self.r
        y_min = y - self.r
        z_min = z - self.r

        x_max = x + self.r + self.grid_spacing
        y_max = y + self.r + self.grid_spacing
        z_max = z + self.r + self.grid_spacing

        n_min = grid.cartesian_to_index((x_min, y_min, z_min))
        n_max = grid.cartesian_to_index((x_max, y_max, z_max))

        indexes = []
        n_range = range(n_min[0], n_max[0] + 1)
        o_range = range(n_min[1], n_max[1] + 1)
        p_range = range(n_min[2], n_max[2] + 1)

        for nop in product(n_range, o_range, p_range):
            p = grid.index_to_cartesian(nop)
            if self.is_in_sphere(p, xyz_c, self.r):
                if grid.is_inside_unit_cell(p, uc):
                    indexes.append(nop)
                else:
                    nop_tr = grid.index_inside_unit_cell(p, uc)
                    indexes.append(nop_tr)
        return indexes

    def find_indices(self, xyz_f, uc, mean_point_translation):
        xyz_f_m = sgtbx.fractional_mod_positive(xyz_f)
        trans = flex.double(xyz_f) - flex.double(xyz_f_m)
        trans = mean_point_translation - trans
        xyz_c = uc.orthogonalize(xyz_f_m)
        return trans, xyz_c

    def find_minimal_chain_distances(self, model, crystal_symmetry):  # pragma: no mccabe
        sps = crystal_symmetry.special_position_settings(min_distance_sym_equiv=0.5)
        min_dist = {}
        coords = {}
        for chain_1 in model.chains():
            chain_id_1 = str(chain_1.id).strip()
            if chain_id_1 not in coords:
                coords[chain_id_1] = flex.vec3_double()

            for atom_1 in chain_1.atoms():
                if (not atom_1.hetero and not atom_1.element_is_hydrogen()):
                    coords[chain_id_1].append(atom_1.xyz)

        for chain_id_1 in coords.keys():
            coords[chain_id_1] = crystal_symmetry.unit_cell().fractionalize(coords[chain_id_1])
            print >> self.log_stream, chain_id_1, coords[chain_id_1].size()

        for chain_id_1, xyz_1 in coords.items():
            for chain_id_2, xyz_2 in coords.items():
                code = "%s_%s" % (chain_id_1, chain_id_2)
                sym_code = "%s_%s" % (chain_id_2, chain_id_1)
                if chain_id_1 != chain_id_2 and sym_code not in min_dist:
                    min_dist[code] = float('inf')
                    for atom_1_xyz in xyz_1:
                        for atom_2_xyz in xyz_2:
                            reference = sps.sym_equiv_sites(atom_1_xyz)
                            info = sgtbx.min_sym_equiv_distance_info(
                                reference_sites=reference,
                                other=atom_2_xyz
                            )
                            dd = info.dist()
                            if min_dist[code] > dd:
                                min_dist[code] = dd
                    print >> self.log_stream, code, min_dist[code]
        return min_dist

    def test_distances(self, data_pdb, model):
        print >> self.log_stream, datetime.datetime.now()
        crystal_symmetry = data_pdb.crystal_symmetry_from_cryst1()
        distance = self.test_minimal_distances(model, crystal_symmetry)
        print >> self.log_stream, 'Minimal distance betwee chains %s' % distance
        print >> self.log_stream, datetime.datetime.now()

    def pack(self):  # pragma: no mccabe
        data_pdb = iotbx.pdb.input(file_name=self.in_pdb)
        pdb_hierarchy = data_pdb.construct_hierarchy()

        xray_structure = data_pdb.xray_structure_simple()
        sg_ops = xray_structure.space_group().all_ops()
        print >> self.log_stream, '######'

        model = pdb_hierarchy.models()[0]
        uc = xray_structure.unit_cell()
        yield None, None

        print >> self.log_stream, 'Generating grid points...'
        self.log_stream.flush()
        # generate grid indexes
        grid = Grid(self.grid_spacing, log_stream=self.log_stream)
        grid.create(uc)
        yield None, None

        print >> self.log_stream, 'Generating tags for grid points...'
        self.log_stream.flush()
        # for each space group operation
        mean_point_translation = {}
        point_dict = {}

        for op_id, op in enumerate(sg_ops):
            # for each chain
            group_coords = {}
            for chain in model.chains():
                chain_id = self.chain_grouping[str(chain.id).strip()]
                coords = group_coords.setdefault(chain_id, flex.vec3_double())
                for atom_xyz in (op(uc.fractionalize(atom.xyz)) for atom in chain.atoms() if
                                 (not atom.hetero and not atom.element_is_hydrogen())):
                    coords.append(atom_xyz)

            for chain_id, coords in group_coords.iteritems():
                if len(coords) > 0:
                    xyz = coords.mean()
                    mean_point_translation[(chain_id, str(op))] = flex.double(xyz) - flex.double(
                        sgtbx.fractional_mod_positive(xyz))
                else:
                    mean_point_translation[(chain_id, str(op))] = (0.0, 0.0, 0.0)

            for chain_str_id in self.chain_grouping.keys():
                # start = datetime.datetime.now()
                chain_id = self.chain_grouping[chain_str_id]
                coords = group_coords[chain_id]

                print >> self.log_stream, 'Adding tag %s, chain %s' % (str(op), chain_str_id)
                self.log_stream.flush()
                for xyz_f in coords:
                    m_point = mean_point_translation[(chain_id, str(op))]
                    trans, xyz_c = self.find_indices(xyz_f, uc, m_point)
                    indexes = self.grid_indexes(grid, xyz_c, uc)

                    cache_id = (chain_id, op_id, trans[0], trans[1], trans[2])
                    try:
                        p = point_dict[cache_id]
                    except Exception:
                        p = GridPoint(chain_id, op, trans[0], trans[1], trans[2], '', op_id)
                        point_dict[cache_id] = p
                    grid.add(indexes, p, uc)
                gc.collect()
                # print >> self.log_stream, '[%d grid points]' % (len(grid))
                # print >> self.log_stream, "took %.2f s" % ((datetime.datetime.now()-start).total_seconds())

        group_coords = None
        point_dict = None
        gc.collect()
        yield None, None

        print >> self.log_stream, '\nAnalyzing tags...\n'

        summed = GridCountUnique(grid, self.log_stream)
        gc.collect()
        print >> self.log_stream, "Common volume:"
        print >> self.log_stream, summed
        print >> self.log_stream, '###'
        yield None, None
        summed = None

        self.log_stream.flush()
        # start = datetime.datetime.now()
        details = GridCountPairSimple(grid, self.log_stream)
        gc.collect()
        # print >> self.log_stream, "took %.2f s" % ((datetime.datetime.now()-start).total_seconds())
        print >> self.log_stream, "Common volume:"
        print >> self.log_stream, details
        print >> self.log_stream, '###'
        yield None, None

        if len(details) > 0:
            transform = details.get_transformations_t(self.function)
        else:
            transform = None
        gc.collect()
        yield transform, mean_point_translation


class Chesym(object):

    def __init__(self, log_stream=None):
        self.log_stream = sys.stdout
        if (type(log_stream) == str or type(log_stream) == unicode):
            self.log_stream = open(log_stream, 'w')
        elif log_stream is not None:
            self.log_stream = log_stream

    def report_progress(self, percent, total):
        pass

    def report_progres_post_grouping(self, percent, total, test_grouping):
        if (test_grouping):
            self.report_progress(0.4 * percent + 0.3 * total, total)
        else:
            self.report_progress(0.7 * percent, total)

    # group like "A,B;C,D"
    # should produce 2 grapus made from chains AB and CD
    def get_group_dict(self, group, chains, letters=False):  # pragma: no mccabe
        group_dict = {}
        chain_id = []
        for chain in chains:
            tag = str(chain.id).strip()
            if tag not in chain_id:
                chain_id.append(tag)

        if (group is not None):
            for (i_part, part) in enumerate(str(group).split(';')):
                for chain_name in part.split(','):
                    tag = chain_name.strip()
                    if tag in chain_id:
                        group_dict[tag] = i_part + 1

        for tag in sorted(chain_id):
            if (tag not in group_dict):
                group_dict[tag] = len(set(group_dict.values())) + 1

        print >> self.log_stream, "Defining group names"
        if letters is True:
            group_dict_renamed = {}
            for code in group_dict.itervalues():
                new_code = ''
                for key in group_dict.iterkeys():
                    if group_dict[key] == code:
                        new_code = new_code + key
                for key in group_dict.iterkeys():
                    if group_dict[key] == code:
                        group_dict_renamed[key] = new_code

            for key in sorted(group_dict_renamed, key=group_dict_renamed.get):
                print >> self.log_stream, "%s belongs to group: %s" % (key, group_dict_renamed[key])
            return group_dict_renamed
        else:
            for key in sorted(group_dict, key=group_dict.get):
                print >> self.log_stream, "%s belongs to group: %s" % (key, group_dict[key])
            return group_dict

    def mean_xyz_without_hetatom(self, pdb_model):
        atoms = flex.vec3_double()
        for chain in pdb_model.chains():
            for atom in chain.atoms():
                if (not atom.hetero):
                    atoms.append(atom.xyz)
        return atoms.mean()

    def get_group_coord(self, chains_coord, tags):
        coord = flex.vec3_double()
        for tag in tags:
            coord.extend(chains_coord[tag])
        return coord

    def get_translation(self, group_coord, chain_coord, unit_cell):
        if len(group_coord) > 0:
            group_mean_cart = group_coord.mean()
            group_mean = unit_cell.fractionalize(group_mean_cart)
        else:
            group_mean = (0.5, 0.5, 0.5)
        chain_mean_cart = chain_coord.mean()
        chain_mean_frac = unit_cell.fractionalize(chain_mean_cart)
        chain_mean_inside = flex.double(sgtbx.fractional_mod_positive(chain_mean_frac))
        t1 = chain_mean_inside - flex.double(chain_mean_frac)
        t2 = find_best_tanslation(group_mean, chain_mean_inside, unit_cell, translations=(0.0, -1.0, 1.0, -2.0, 2.0))
        t_sum = (t1[0] + t2[0], t1[1] + t2[1], t1[2] + t2[2])
        return t_sum

    def packing_by_contact_analysis(self, chain_aliases, in_pdb, grid_spacing, radius, function):
        p = Packer(in_pdb, chain_grouping=chain_aliases, log_stream=self.log_stream, grid_spacing=grid_spacing,
                   r=radius, function=function)
        progress = 0.0
        chain_sym_op = None
        steps = None
        for chain_sym_op_tmp, steps_tmp in p.pack():
            progress += 4.0
            self.report_progress(progress, 100.0)
            chain_sym_op = chain_sym_op_tmp
            steps = steps_tmp
            self.log_stream.flush()
        return chain_sym_op, steps

    def packing_by_bounding_boxes(self, group_dict, model, xray_structure):  # pragma: no mccabe
        print >> self.log_stream, 'Searching the most compact chains placement'
        # find the most complex placement
        points = {}
        matching = {}
        translation = {}
        best_box_vol = None
        best_box_area = None
        best_box_lenght = None
        best_op = None
        best_translation = None
        sg_ops = xray_structure.space_group().all_ops()
        # sg_cob = [sgtbx.change_of_basis_op(op) for op in sg_ops]

        # prepare list of atoms in chains
        for chain in model.chains():
            point_flex_vec = flex.vec3_double()
            for atom in chain.atoms():
                # print >> self.log_stream, chain.id, atom.name, atom.xyz
                if (not atom.hetero):
                    point_flex_vec.append(atom.xyz)
            vec = points.setdefault(group_dict[chain.id], flex.vec3_double())
            vec.extend(xray_structure.unit_cell().fractionalize(point_flex_vec))
        k = len(points)
        print >> self.log_stream, 'Structure', model.id, 'has', k, 'chain groups'

        self.report_progress(1.0, 100.0)

        # for each chain generate symmetry (len(sg_ops)^k combinations)
        progress_len = len(sg_ops) ** k
        if progress_len > MAX_COMBINATIONS:
            raise Exception(
                "Too many combinations too test (%s). Reduce number of chains or group them together!" % progress_len)

        for op_i, op_list in enumerate(product(sg_ops, repeat=k)):
            translation_list = []
            points_trasposed_temp_list = []
            points_trasposed_mean_list = []
            # for each chain
            print >> self.log_stream, 'Testing:', '; '.join(
                ["Group %s: %s" % (key, str(op)) for (key, op) in izip(points.keys(), op_list)])
            for op, p_key in izip(op_list, points.keys()):
                point_flex_vec = points[p_key]
                points_trasposed_temp = flex.vec3_double()

                # apply symmetry operation
                # points are fractional
                for p in point_flex_vec:
                    points_trasposed_temp.append(op(p))

                mean = points_trasposed_temp.mean()
                tr_uc = tuple(flex.double(sgtbx.fractional_mod_positive(mean)) - flex.double(mean))
                transltion_to_uc = flex.vec3_double(len(points_trasposed_temp), tr_uc)

                points_trasposed_temp = points_trasposed_temp + transltion_to_uc
                mean = points_trasposed_temp.mean()

                tr = (0.0, 0.0, 0.0)
                if (len(points_trasposed_mean_list) > 0):
                    tr = find_best_tanslation(points_trasposed_mean_list[0], mean, xray_structure.unit_cell())
                    print >> self.log_stream, "    Group", p_key, "best translation ->", tr
                    vec_tr = flex.vec3_double(len(points_trasposed_temp), tr)
                    points_trasposed_temp = points_trasposed_temp + vec_tr
                    mean = points_trasposed_temp.mean()

                points_trasposed_mean_list.append(mean)
                points_trasposed_temp_list.append(points_trasposed_temp)
                translation_list.append((tr_uc[0] + tr[0], tr_uc[1] + tr[1], tr_uc[2] + tr[2]))
                gc.collect()

            points_trasposed = flex.vec3_double()
            # should be in cartesian
            for points_vec in points_trasposed_temp_list:
                points_trasposed.extend(xray_structure.unit_cell().orthogonalize(points_vec))

            # calculate PCA
            pca = principal_axes_of_inertia(points=points_trasposed)
            # calculate bounding box
            es = pca.eigensystem()
            eigen_vec = flex.vec3_double()
            for i in range(3):
                eigen_vec.append((es.vectors()[3 * i], es.vectors()[3 * i + 1], es.vectors()[3 * i + 2]))

            mini = [0.0] * len(eigen_vec)
            maxi = [0.0] * len(eigen_vec)
            for (i, vec) in enumerate(eigen_vec):
                p = flex.vec3_double([points_trasposed[0]])
                vec = flex.vec3_double([vec])
                u = p.dot(vec) / vec.dot(vec)
                mini[i] = u[0]
                maxi[i] = u[0]
                for p in points_trasposed:
                    p = flex.vec3_double([p])
                    u = p.dot(vec) / vec.dot(vec)
                    mini[i] = min(mini[i], u[0])
                    maxi[i] = max(maxi[i], u[0])
            box_dim = sorted(list(flex.double(maxi) - flex.double(mini)))
            box_vol = box_dim[0] * box_dim[1] * box_dim[2]
            box_area = 2.0 * (box_dim[0] * box_dim[1] + box_dim[0] * box_dim[2] + box_dim[1] * box_dim[2])
            box_lenght = 4.0 * (box_dim[0] + box_dim[1] + box_dim[2])

            # approx equal
            box_vol = fround(box_vol, 6)
            box_area = fround(box_area, 6)
            box_lenght = fround(box_lenght, 6)

            print >> self.log_stream, '       Box Volume %.5e Box Area %.5e Box Length %.5e' % (
                box_vol, box_area, box_lenght
            )

            # find minimal volume
            if (
                best_box_vol is None or
                box_vol < best_box_vol or
                (box_vol == best_box_vol and best_box_area < box_area) or
                (box_vol == best_box_vol and best_box_area == box_area and box_lenght < best_box_lenght)
            ):
                best_op = op_list
                best_box_vol = box_vol
                best_box_area = box_area
                best_box_lenght = box_lenght
                best_translation = translation_list

            self.report_progress(0.3 * float(op_i) / float(progress_len), 1.0)

        [
            translation.setdefault(key, best_tr) for (key, best_tr) in izip(points.keys(), best_translation)
        ]
        [
            matching.setdefault(key, op) for (key, op) in izip(points.keys(), best_op)
        ]
        print >> self.log_stream, '\n#\nThe best:', '; '.join(
            ["Group %s: %s" % (key, str(op)) for (key, op) in izip(points.keys(), best_op)]
        )
        print >> self.log_stream, '       Box Volume %.5e Box Area %.5e Box Lenght %.5e' % (
            best_box_vol, best_box_area, best_box_lenght
        )
        print >> self.log_stream, ''  # Applying symmetry to chains
        return translation, matching

    def print_tansformations(self, applied_chains, applied_rot_in_cart, applied_trans_in_cart, grid_spacing, radius,
                             function, out=sys.stdout):
        # print >> out, "         1         2         3         4         5         6         7         8"
        # print >> out, "12345678901234567890123456789012345678901234567890123456789012345678901234567890"
        print >> out, 'REMARK ACHESYM    Structure transformed with aCHESYM v%s, %s' % (
            ACHESYM_VERSION, datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        )
        print >> out, 'REMARK ACHESYM    "ACHESYM: an algorithm and program for standardized '
        print >> out, 'REMARK ACHESYM    placement of macromolecular models in the unit cell"'
        print >> out, 'REMARK ACHESYM    Kowiel, M., Jaskolski, M. & Dauter, Z. (2014). Acta Cryst. D70'
        print >> out, 'REMARK ACHESYM    '
        print >> out, 'REMARK ACHESYM    grid spacing: %.2f, radius: %.2f, function: %s' % (
            grid_spacing, radius, function
        )
        print >> out, 'REMARK ACHESYM    '
        i = 1
        for chain, rot, trans in zip(applied_chains, applied_rot_in_cart, applied_trans_in_cart):
            print >> out, 'REMARK ACHESYM1%3d %3s%10.6f%10.6f%10.6f %10.5f' % (
                i, chain, fix_zero(rot[0]), fix_zero(rot[1]), fix_zero(rot[2]), fix_zero(trans[0])
            )
            print >> out, 'REMARK ACHESYM2%3d %3s%10.6f%10.6f%10.6f %10.5f' % (
                i, chain, fix_zero(rot[3]), fix_zero(rot[4]), fix_zero(rot[5]), fix_zero(trans[1])
            )
            print >> out, 'REMARK ACHESYM3%3d %3s%10.6f%10.6f%10.6f %10.5f' % (
                i, chain, fix_zero(rot[6]), fix_zero(rot[7]), fix_zero(rot[8]), fix_zero(trans[2])
            )
            i += 1

    def apply_reindex(self, block, block_name, vector):
        if (block_name in block):
            for i, hkl in enumerate(vector):
                block[block_name][i] = hkl

    def shift_phase(self, block, block_name, t):
        if ('_refln_index_h' in block and '_refln_index_k' in block and '_refln_index_l' in block):
            block_h = block['_refln_index_h']
            block_k = block['_refln_index_k']
            block_l = block['_refln_index_l']
        elif ('_refln.index_h' in block and '_refln.index_k' in block and '_refln.index_l' in block):
            block_h = block['_refln.index_h']
            block_k = block['_refln.index_k']
            block_l = block['_refln.index_l']

        ta, tb, tc = t
        if (block_name in block):
            for i, (h, k, l) in enumerate(zip(block_h, block_k, block_l)):
                if block[block_name][i].strip() != '?':
                    block[block_name][i] = ("%.1f" % math.fmod(
                        float(block[block_name][i]) + 360.0 * (int(h) * ta + int(k) * tb + int(l) * tc), 360.0))

    def shift_phase_AB(self, block, block_name_A, block_name_B, t):
        if ('_refln_index_h' in block and '_refln_index_k' in block and '_refln_index_l' in block):
            block_h = block['_refln_index_h']
            block_k = block['_refln_index_k']
            block_l = block['_refln_index_l']
        elif ('_refln.index_h' in block and '_refln.index_k' in block and '_refln.index_l' in block):
            block_h = block['_refln.index_h']
            block_k = block['_refln.index_k']
            block_l = block['_refln.index_l']

        ta, tb, tc = t
        if (block_name_A in block and block_name_B in block):
            for i, (h, k, l) in enumerate(zip(block_h, block_k, block_l)):
                if block[block_name_A][i].strip() != '?' and block[block_name_B][i].strip() != '?':
                    phase = 2.0 * math.pi * (int(h) * ta + int(k) * tb + int(l) * tc)
                    cos = math.cos(phase)
                    sin = math.sin(phase)
                    phase_A = float(block[block_name_A][i])
                    phase_B = float(block[block_name_B][i])
                    block[block_name_A][i] = ("%.1f" % (phase_A * cos - phase_B * sin))
                    block[block_name_B][i] = ("%.1f" % (phase_A * sin + phase_B * cos))

    def apply_phase_shift(self, block, t):
        # x' = x Q +q
        # P = Q^-1
        #
        # h' = h P
        # h' = h Q^-1
        #
        # q = -P^-1 p
        # p = -P q
        # p = -Q^-1 q
        # phi'(h') = phi(h) - 2pi h p
        # phi'(h') = phi(h) - 2pi h -Q^-1 q
        # phi'(h') = phi(h) + 2pi h Q^-1 q
        # phi'(h') = phi(h) + 2pi h' q
        # phi'(h') = phi(h) + 360 h' q

        self.shift_phase(block, '_refln_phase_calc', t)
        self.shift_phase(block, '_refln.phase_calc', t)
        self.shift_phase(block, '_refln_phase_meas', t)
        self.shift_phase(block, '_refln.phase_meas', t)

        self.shift_phase_AB(block, '_refln_A_meas', '_refln_B_meas', t)
        self.shift_phase_AB(block, '_refln.A_meas', '_refln.B_meas', t)
        self.shift_phase_AB(block, '_refln_A_calc', '_refln_B_calc', t)
        self.shift_phase_AB(block, '_refln.A_calc', '_refln.B_calc', t)

    def copy_remarks(self, in_pdb, out):
        in_pdb = open(in_pdb, 'r')
        ignored_tags = (
            'ATOM',
            'HETATM',
            'ANISOU',
            # 'REMARK   3',
            'CRYST',
            'ORIGX',
            'SCALE',
            'TER',
            'END',
            'CONECT ',
            'MASTER ',
        )
        print >> out, 'REMARK ACHESYM    DATA FROM THE ORIGINAL FILE'
        print >> out, 'REMARK ACHESYM    MAY NEED MANUAL CORRECTIONS'
        for line in in_pdb:
            if not line.startswith(ignored_tags):
                print >> out, line,
        in_pdb.close()

    def copy_atoms(self, stringio_pdb, out):
        ignored_tags = (
            'CONECT',
            'END',
        )
        for line in stringio_pdb.getvalue().splitlines():
            if line.startswith(ignored_tags):
                print >> out, line

    def copy_remarks_after(self, in_pdb, out):
        in_pdb = open(in_pdb, 'r')
        after_tags = (
            'CONECT',
            'END',
        )
        for line in in_pdb:
            if line.startswith(after_tags):
                print >> out, line,
        in_pdb.close()

    def get_model_chains(self, model):
        chains = set()
        for chain in model.chains():
            chains.add(str(chain.id))
        chains = sorted(chains)
        return chains

    def combine_transformations(self, chains, sym_in_frac):
        chains_sorted = set()
        for chain in chains:
            chains_sorted.add(str(chain))
        chains_sorted = sorted(chains_sorted)
        out_sym_in_frac = []
        for chain in chains_sorted:
            rot_chain = rt_mx().r()
            tr_chain = (0, 0, 0)
            for op_chain, rot, tra in sym_in_frac:
                if chain == op_chain or op_chain == 'ALL':
                    rot_chain = rot.multiply(rot_chain)
                    if isinstance(tra, sgtbx.tr_vec):
                        tra = tra.as_double()
                    if isinstance(tr_chain, sgtbx.tr_vec):
                        tr_chain = tr_chain.as_double()
                    tr_chain = matrix.sqr(rot.as_double()) * matrix.rec(tr_chain, (3, 1)) + matrix.rec(tra, (3, 1))
            out_sym_in_frac.append((chain, rot_chain, tr_chain.elems))
        return out_sym_in_frac

    def orthogonalize(self, sym_applied_in_fractional, unit_cell):
        applied_chains = []
        applied_rot_in_cart = []
        applied_trans_in_cart = []
        for chain, rot, trans in sym_applied_in_fractional:
            r = unit_cell.matrix_cart(rot)
            if isinstance(trans, sgtbx.tr_vec):
                trans = trans.as_double()
            t = unit_cell.orthogonalize(trans)
            applied_chains.append(chain)
            applied_rot_in_cart.append(r)
            applied_trans_in_cart.append(t)
        return applied_chains, applied_rot_in_cart, applied_trans_in_cart

    def chesym(self, in_pdb='in.pdb', out_pdb='out.pdb', test_compactness=False,
               in_cif='', out_cif='out.cif', group=None, analyze_contacts=False,
               grid_spacing=0.7, radius=2.3, function='sum'):  # pragma: no mccabe

        # validate cell
        with open(in_pdb, 'r') as pdb_file:
            txt = pdb_file.read()
            if 'CRYST' not in txt:
                raise Exception('Missing unit cell informaton, CRYST record is required')
            if 'SCALE' not in txt:
                raise Exception('Missing fractional transcormation matrix, SCALE records are required')
            txt = ''

        data_pdb = iotbx.pdb.input(file_name=in_pdb)
        pdb_hierarchy = data_pdb.construct_hierarchy()

        # the xray_structure_simple mix up order
        # xray_structure = data_pdb.xray_structure_simple()
        xray_structure = pdb_hierarchy.extract_xray_structure(crystal_symmetry=data_pdb.crystal_symmetry())

        sg = xray_structure.space_group()
        sg_type = sg.type()
        sg_number = sg_type.number()
        sg_ops = sg.all_ops()
        print >> self.log_stream, "Space Group No: ", sg_number
        print >> self.log_stream, "Space Group operations: "
        for od_id, op in enumerate(sg_ops):
            print >> self.log_stream, "%s: %s" % (od_id + 1, op)
        self.log_stream.flush()

        normalizer = Normalizer.factory(sg_number)

        self.report_progress(0.0, 100.0)

        assert len(pdb_hierarchy.models()) == 1
        model = pdb_hierarchy.models()[0]
        chains_size = model.chains_size()  # number of chains
        test_grouping = chains_size > 1 and test_compactness is True
        ########################
        mass = self.mean_xyz_without_hetatom(model)
        fmass = xray_structure.unit_cell().fractionalize(mass)
        print >> self.log_stream, "\nINPUT Mean in Cartesian coordinates: (%.3f, %.3f, %.3f)" % mass
        print >> self.log_stream, "INPUT Mean in fractional coordinates: (%.3f, %.3f, %.3f)" % fmass
        self.log_stream.flush()
        sym_applied_in_fractional = []
        sym_applied_in_fractional_tmp = OrderedDict()
        sym_applied_to_print = OrderedDict()
        ########################
        if chains_size > 1 and analyze_contacts is True:
            print >> self.log_stream, '\nANALYZE CONTACTS\n'
            self.log_stream.flush()

            chain_aliases = self.get_group_dict(group, model.chains(), letters=True)
            if len(chain_aliases) > 1:
                grid_points, mean_trans = self.packing_by_contact_analysis(chain_aliases, in_pdb, grid_spacing, radius,
                                                                           function)
                if grid_points is not None:
                    unit_cell = xray_structure.unit_cell()

                    for point in grid_points:
                        for chain in model.chains():
                            if chain_aliases[str(chain.id).strip()] == point.chain_tag:
                                op = point.sym_op
                                # op_r = op.r().as_double()
                                op_in_car = matrix.sqr(unit_cell.orthogonalization_matrix())
                                op_in_car = op_in_car * matrix.sqr(op.r().as_double())
                                op_in_car = op_in_car * matrix.sqr(unit_cell.orthogonalization_matrix()).inverse()
                                op_in_car_t = op_in_car.transpose()
                                op_in_car_elems = op_in_car.elems
                                t0 = mean_trans[(point.chain_tag, str(op))]
                                t0 = (-t0[0], -t0[1], -t0[2])
                                # print >> self.log_stream, "MM %s %s %s %s" % (point.chain_tag, t0[0], t0[1], t0[2])
                                op_t0_in_car = unit_cell.orthogonalize(t0)
                                t1 = op.t().as_double()
                                # print >> self.log_stream, "OO %s %s %s %s" % (point.chain_tag, t1[0], t1[1], t1[2])
                                op_t1_in_car = unit_cell.orthogonalize(t1)
                                t2 = (point.tx, point.ty, point.tz)
                                # print >> self.log_stream, "DD %s %s %s %s" % (point.chain_tag, t2[0], t2[1], t2[2])
                                op_t2_in_car = unit_cell.orthogonalize(t2)
                                t_sum = (point.tx + t0[0], point.ty + t0[1], point.tz + t0[2])
                                xyz = chain.atoms().extract_xyz()
                                chain.atoms().set_xyz(
                                    (op_in_car_elems * xyz) + op_t0_in_car + op_t1_in_car + op_t2_in_car)
                                for atom in chain.atoms():
                                    if atom.uij_is_defined():
                                        new_uij = op_in_car * matrix.sym(sym_mat3=atom.uij) * op_in_car_t
                                        atom.set_uij(new_uij.as_sym_mat3())
                                # may be added two times
                                sym_applied_in_fractional_tmp[str(chain.id).strip() + '@1'] = (
                                    str(chain.id).strip(), op.r(), op.t()
                                )
                                sym_applied_in_fractional_tmp[str(chain.id).strip() + '@2'] = (
                                    str(chain.id).strip(), rt_mx().r(), t_sum
                                )
                                sym_applied_to_print[str(chain.id).strip()] = (
                                    point.sym_op.as_xyz(), t_sum[0], t_sum[1], t_sum[2]
                                )

                    print >> self.log_stream, 'Proposed transformation:'
                    for chain_tag, value in sym_applied_to_print.iteritems():
                        sym_op_xyz, ptx, pty, ptz = value
                        print >> self.log_stream, "chain %s:" % (chain_tag)
                        print >> self.log_stream, "  R: %s" % (sym_op_xyz)
                        print >> self.log_stream, "  T: %d,%d,%d" % (ptx, pty, ptz)
                    sym_applied_in_fractional.extend(sym_applied_in_fractional_tmp.values())

            xray_structure = pdb_hierarchy.extract_xray_structure(crystal_symmetry=data_pdb.crystal_symmetry())
            self.report_progress(30.0, 100.0)
            self.log_stream.flush()
        ########################

        elif (test_compactness and test_grouping):
            group_dict = self.get_group_dict(group, model.chains())
            translation, matching = self.packing_by_bounding_boxes(group_dict, model, xray_structure)

            for chain in model.chains():
                for atom in chain.atoms():
                    xyz_in_frac = xray_structure.unit_cell().fractionalize(atom.xyz)
                    xyz_in_frac = matching[group_dict[chain.id]](xyz_in_frac)
                    xyz_in_frac = (xyz_in_frac[0] + translation[group_dict[chain.id]][0],
                                   xyz_in_frac[1] + translation[group_dict[chain.id]][1],
                                   xyz_in_frac[2] + translation[group_dict[chain.id]][2])
                    xyz_in_cart = xray_structure.unit_cell().orthogonalize(xyz_in_frac)
                    atom.set_xyz(xyz_in_cart)
            xray_structure = pdb_hierarchy.extract_xray_structure(crystal_symmetry=data_pdb.crystal_symmetry())
            self.report_progress(30.0, 100.0)

        ########################

        gc.collect()
        print >> self.log_stream, '#'
        mass = self.mean_xyz_without_hetatom(pdb_hierarchy.models()[0])
        fmass = xray_structure.unit_cell().fractionalize(mass)
        print >> self.log_stream, "\nMean in Cartesian coordinates: (%.3f, %.3f, %.3f)" % mass
        print >> self.log_stream, "Mean in fractional coordinates: (%.3f, %.3f, %.3f)" % fmass

        # print >> self.log_stream, "Initial translation T1: (%.1f, %.1f, %.1f)" % (
        #   fix_zero(-math.floor(fmass[0])), fix_zero(-math.floor(fmass[1])), fix_zero(-math.floor(fmass[2]))
        # )

        t1 = (-math.floor(fmass[0]), -math.floor(fmass[1]), -math.floor(fmass[2]))
        fmass = (fmass[0] - math.floor(fmass[0]), fmass[1] - math.floor(fmass[1]), fmass[2] - math.floor(fmass[2]))

        xray_structure_cpy = xray_structure.deep_copy_scatterers()
        xray_structure_cpy.erase_scatterers()
        xray_structure_cpy.add_scatterer(xray.scatterer(
            label="C",
            site=sgtbx.fractional_mod_positive(fmass),
            u=0.2))

        # print >> self.log_stream, "... Mean in fractional coordinates: (%.3f, %.3f, %.3f)" % (
        #   fix_zero_3(xray_structure_cpy.sites_frac()[0])
        # )

        (x, y, z) = normalizer.continuous_traslation()
        t2 = (0.0, 0.0, 0.0)
        if (x > 0.0 or y > 0.0 or z > 0.0):
            t2 = ((-fmass[0] + 0.5) * x, (-fmass[1] + 0.5) * y, (-fmass[2] + 0.5) * z)
            # print >> self.log_stream, "Applying translation T2: (%.3f, %.3f, %.3f)" % fix_zero_3(t2)
            translate_in_fractional(xray_structure_cpy, t2)

        t1 = (t1[0] + t2[0], t1[1] + t2[1], t1[2] + t2[2])
        print >> self.log_stream, "Initial translation T1: (%.3f, %.3f, %.3f)" % (
            fix_zero(t1[0]), fix_zero(t1[1]), fix_zero(t1[2])
        )
        print >> self.log_stream, "... Mean in fractional coordinates: (%.3f, %.3f, %.3f)" % fix_zero_3(
            xray_structure_cpy.sites_frac()[0]
        )

        print >> self.log_stream, "#"

        best_t1 = t1
        best_m1 = sgtbx.change_of_basis_op('x, y, z')
        best_m2 = sgtbx.change_of_basis_op('x, y, z')
        best_t3 = (0.0, 0.0, 0.0)
        found = False

        gc.collect()
        progress_len = float(len(normalizer.coset()) * len(sg_ops))
        progress_i = 0
        for op in normalizer.coset():
            cob_op = sgtbx.change_of_basis_op(op)

            for sym in sg_ops:
                sym = sgtbx.change_of_basis_op(sym)
                print >> self.log_stream, "\n  Testing: Normalizer operation N: %s, Symmetry operation R: %s" % (
                    cob_op.as_xyz(), sym.as_xyz()
                )
                mult_op = cob_op.new_denominators(sym)
                xs = xray_structure_cpy.change_basis(mult_op)
                xs = xs.change_basis(sym)
                pos = xs.sites_frac()[0]
                (x, y, z) = normalizer.continuous_traslation()
                t3 = (-math.floor(pos[0]) - x * (pos[0] - math.floor(pos[0]) - 0.5),
                      -math.floor(pos[1]) - y * (pos[1] - math.floor(pos[1]) - 0.5),
                      -math.floor(pos[2]) - z * (pos[2] - math.floor(pos[2]) - 0.5))
                print >> self.log_stream, "  Translation after N and R operations, T2: (%.3f, %.3f, %.3f)" % fix_zero_3(
                    t3)

                (x, y, z) = (pos[0] + t3[0], pos[1] + t3[1], pos[2] + t3[2])
                print >> self.log_stream, "  ... Mean in fractional coordinates: (%.3f, %.3f, %.3f)" % (x, y, z)
                if (normalizer.is_in_anticheshire_cell(x, y, z)):
                    (ta, tb, tc) = normalizer.anticheshire_post_translation(x, y, z)
                    t3 = (t3[0] + ta, t3[1] + tb, t3[2] + tc)
                    print >> self.log_stream, "\n  Success: "
                    print >> self.log_stream, "    T1: (%.3f, %.3f, %.3f)" % fix_zero_3(t1)
                    print >> self.log_stream, "    N: xyz: %s      hkl: %s" % (
                        cob_op.as_xyz(), cob_op.c_inv().r().as_hkl()
                    )
                    print >> self.log_stream, "    R: xyz: %s" % (sym.as_xyz())
                    print >> self.log_stream, "    T2: (%.3f, %.3f, %.3f)" % fix_zero_3(t3)
                    print >> self.log_stream, "    N_hkl:     ", mult_op.c_inv().r().as_hkl(False, 'hkl', ',')
                    print >> self.log_stream, "    Det(N_hkl):", mult_op.c().r().determinant()
                    best_m1 = cob_op
                    best_m2 = sym
                    best_t3 = t3
                    found = True
                progress_i += 1
                self.report_progres_post_grouping(0.9 * float(progress_i) / progress_len, 1.0, test_grouping)

        print >> self.log_stream, '####'
        self.log_stream.flush()

        if found is False:
            raise Exception("Transformation not found")

        print >> self.log_stream, 'saving %s' % os.path.basename(out_pdb)

        unit_cell = xray_structure.unit_cell()
        translate_in_fractional(xray_structure, best_t1)
        best_m1 = best_m1.new_denominators(best_m2)
        xray_structure = xray_structure.change_basis(best_m1)
        xray_structure = xray_structure.change_basis(best_m2)
        sites = xray_structure.sites_frac()
        sites = sites + best_t3
        xray_structure.set_sites_frac(sites)

        sym_applied_in_fractional.append(("ALL", rt_mx().r(), best_t1))
        sym_applied_in_fractional.append(("ALL", best_m1.c().r(), best_m1.c().t()))
        sym_applied_in_fractional.append(("ALL", best_m2.c().r(), best_m2.c().t()))
        sym_applied_in_fractional.append(("ALL", rt_mx().r(), best_t3))

        applied_chains, applied_rot_in_cart, applied_trans_in_cart = self.orthogonalize(sym_applied_in_fractional,
                                                                                        unit_cell)

        transformator = TLSTranformator(in_pdb)
        out_tlsos, selection_strings = transformator.transform(applied_chains, applied_rot_in_cart,
                                                               applied_trans_in_cart)

        chains = self.get_model_chains(model)
        sym_combined_in_frac = self.combine_transformations(chains, sym_applied_in_fractional)
        combined_chains, combined_rot_in_cart, combined_trans_in_cart = self.orthogonalize(sym_combined_in_frac,
                                                                                           unit_cell)

        out = open(out_pdb, 'w')
        self.print_tansformations(combined_chains, combined_rot_in_cart, combined_trans_in_cart, grid_spacing, radius,
                                  function, out=out)

        if len(out_tlsos) > 0:
            tools.remark_3_tls(out_tlsos, selection_strings, out=out)

        self.copy_remarks(in_pdb=in_pdb, out=out)

        pdb_hierarchy.adopt_xray_structure(xray_structure, assert_identical_id_str=True)
        pdb_str = pdb_hierarchy.as_pdb_string(
            append_end=False,
            crystal_symmetry=xray_structure.crystal_symmetry(),
            atoms_reset_serial_first_value=None,
            interleaved_conf=1,
        )
        out.write(pdb_str)

        self.copy_remarks_after(in_pdb=in_pdb, out=out)
        out.close()

        gc.collect()
        print >> self.log_stream, '####'
        self.report_progres_post_grouping(95.0, 100.0, test_grouping)

        if (in_cif != ''):
            cif_reader = cif.reader(file_path=in_cif)
            cif_model = cif_reader.model()

            print >> self.log_stream, 'saving %s' % os.path.basename(out_cif)
            for key, block in cif_model.items():

                vec_h = None
                vec_k = None
                vec_l = None
                if ('_refln_index_h' in block):
                    vec_h = block['_refln_index_h']
                if ('_refln_index_k' in block):
                    vec_k = block['_refln_index_k']
                if ('_refln_index_l' in block):
                    vec_l = block['_refln_index_l']

                if ('_refln.index_h' in block):
                    vec_h = block['_refln.index_h']
                if ('_refln.index_k' in block):
                    vec_k = block['_refln.index_k']
                if ('_refln.index_l' in block):
                    vec_l = block['_refln.index_l']

                if vec_h is None or vec_k is None or vec_l is None:
                    raise Exception(
                        'hkl data missing. Program reads coordinates only in .pdb format and '
                        'structure factors in .cif format'
                    )
                if len(vec_h) != len(vec_k) or len(vec_h) != len(vec_l):
                    raise Exception('Program failed to load all hkl')
                out_vec_h = list()
                out_vec_k = list()
                out_vec_l = list()
                assert len(vec_h) == len(vec_k)
                assert len(vec_h) == len(vec_l)

                for (h_, k_, l_) in zip(vec_h, vec_k, vec_l):
                    h_ = int(h_)
                    k_ = int(k_)
                    l_ = int(l_)

                    new_hkl = best_m1.apply((h_, k_, l_))

                    out_vec_h.append(str(new_hkl[0]))
                    out_vec_k.append(str(new_hkl[1]))
                    out_vec_l.append(str(new_hkl[2]))

                self.apply_phase_shift(block, best_t1)

                self.apply_reindex(block, '_refln_index_h', out_vec_h)
                self.apply_reindex(block, '_refln_index_k', out_vec_k)
                self.apply_reindex(block, '_refln_index_l', out_vec_l)

                self.apply_reindex(block, '_refln.index_h', out_vec_h)
                self.apply_reindex(block, '_refln.index_k', out_vec_k)
                self.apply_reindex(block, '_refln.index_l', out_vec_l)

                t = best_m1.c().t().as_double()
                t = (t[0] + best_t3[0], t[1] + best_t3[1], t[2] + best_t3[2])
                self.apply_phase_shift(block, t)

                out = open(out_cif, 'w')
                out.write(str(cif_model))
                out.close()

        print >> self.log_stream, '####'
        self.report_progres_post_grouping(100.0, 100.0, test_grouping)


def main():  # pragma: no mccabe
    usage = """
Usage: achesym.py input_file.pdb output_file.pdb [structure_factors.cif output_structure_factors.cif]
                 [-analyze] [-pack] [-group "A,B;C,D"]

-analyze: Find the most compact assembly (with "common volume" analysis)
-pack: Find the most compact assembly (with bounding boxes) - deprecated
-group "A,B;C,D": force chain grouping. In the example chains A and B will be
                  treat like one object and chains C and D as a second object
                  group separator ";"
                  chain separator ","

Program requires cctbx, and can be run by cctbx aware python interpreter
For example:
#> cctbx.python achesym.py -analyze in.pdb out.pdb
    """

    if '-h' in sys.argv or '-help' in sys.argv:
        sys.stderr.write(usage)
        sys.exit(1)

    analyze_contacts = False
    if '-analyze' in sys.argv:
        analyze_contacts = True
        sys.argv.remove('-analyze')
    if '-a' in sys.argv:
        analyze_contacts = True
        sys.argv.remove('-a')

    test_compactness = False
    if '-pack' in sys.argv:
        test_compactness = True
        sys.argv.remove('-pack')
    if '-p' in sys.argv:
        test_compactness = True
        sys.argv.remove('-p')

    group = None
    if '-group' in sys.argv:
        i_group = sys.argv.index('-group')
        sys.argv.pop(i_group)
        group = sys.argv.pop(i_group)
    if '-g' in sys.argv:
        i_group = sys.argv.index('-g')
        sys.argv.pop(i_group)
        group = sys.argv.pop(i_group)

    if (len(sys.argv) < 2):
        sys.stderr.write(usage)
        sys.exit(1)

    if (not os.path.exists(sys.argv[1])):
        sys.stderr.write('ERROR: File %s was not found!' % sys.argv[1])
        sys.exit(1)

    if (not str(sys.argv[1]).endswith(".pdb")):
        sys.stderr.write('ERROR: File %s should be .pdb file!' % sys.argv[1])
        sys.exit(1)
    in_pdb = sys.argv[1]

    if (len(sys.argv) > 2):
        out_pdb = sys.argv[2]
    else:
        out_pdb = 'out_%s' % in_pdb

    if (len(sys.argv) > 3):
        if (not os.path.exists(sys.argv[3])):
            sys.stderr.write('ERROR: File %s was not found!' % sys.argv[3])
            sys.exit(1)

        if (not str(sys.argv[3]).endswith(".cif")):
            sys.stderr.write('ERROR: File %s should be .cif file!' % sys.argv[3])
            sys.exit(1)

        in_cif = sys.argv[3]
    else:
        in_cif = ''

    if (len(sys.argv) > 4):
        out_cif = sys.argv[4]
    else:
        out_cif = 'out_%s' % in_cif

    chesym = Chesym()
    chesym.chesym(
        in_pdb=in_pdb,
        out_pdb=out_pdb,
        test_compactness=test_compactness,
        in_cif=in_cif,
        out_cif=out_cif,
        group=group,
        analyze_contacts=analyze_contacts,
        radius=2.3,
        grid_spacing=0.7,
        function='sum',
    )


if __name__ == '__main__':
    main()
