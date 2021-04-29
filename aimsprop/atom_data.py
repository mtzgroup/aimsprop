from . import units

"""
    Contains three dictionaries:
        atom_symbol_table: {atomic number (int) : atom symbol (str) - all caps}
        mass_table: {atom symbol (str) - all caps : mass (float) - amu }
        extended_mass_table: {atom symbol (str) - all caps : [(isotope mass (float) - amu, isotope composition (float))}

    Masses in dictionaries are automatically converted to a.u.

    Citation:
        NIST Standard Reference Database 144
        Online: September 1999 | Last update: January 2015
        URL: https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ele=&ascii=ascii&isotype=all
        Values obtained: September 15, 2017

"""

atom_symbol_table = {
    1: "H",
    2: "HE",
    3: "LI",
    4: "BE",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "NE",
    11: "NA",
    12: "MG",
    13: "AL",
    14: "SI",
    15: "P",
    16: "S",
    17: "CL",
    18: "AR",
    19: "K",
    20: "CA",
    21: "SC",
    22: "TI",
    23: "V",
    24: "CR",
    25: "MN",
    26: "FE",
    27: "CO",
    28: "NI",
    29: "CU",
    30: "ZN",
    31: "GA",
    32: "GE",
    33: "AS",
    34: "SE",
    35: "BR",
    36: "KR",
    37: "RB",
    38: "SR",
    39: "Y",
    40: "ZR",
    41: "NB",
    42: "MO",
    43: "TC",
    44: "RU",
    45: "RH",
    46: "PD",
    47: "AG",
    48: "CD",
    49: "IN",
    50: "SN",
    51: "SB",
    52: "TE",
    53: "I",
    54: "XE",
    55: "CS",
    56: "BA",
    57: "LA",
    58: "CE",
    59: "PR",
    60: "ND",
    61: "PM",
    62: "SM",
    63: "EU",
    64: "GD",
    65: "TB",
    66: "DY",
    67: "HO",
    68: "ER",
    69: "TM",
    70: "YB",
    71: "LU",
    72: "HF",
    73: "TA",
    74: "W",
    75: "RE",
    76: "OS",
    77: "IR",
    78: "PT",
    79: "AU",
    80: "HG",
    81: "TL",
    82: "PB",
    83: "BI",
    84: "PO",
    85: "AT",
    86: "RN",
    87: "FR",
    88: "RA",
    89: "AC",
    90: "TH",
    91: "PA",
    92: "U",
    93: "NP",
    94: "PU",
    95: "AM",
    96: "CM",
    97: "BK",
    98: "CF",
    99: "ES",
    100: "FM",
    101: "MD",
    102: "NO",
    103: "LR",
    104: "RF",
    105: "DB",
    106: "SG",
    107: "BH",
    108: "HS",
    109: "MT",
    110: "DS",
    111: "RG",
    112: "CN",
    113: "NH",
    114: "FL",
    115: "MC",
    116: "LV",
    117: "TS",
    118: "OG",
}

mass_table = {
    "H": 1.00782503223,
    "HE": 3.0160293201,
    "LI": 6.0151228874,
    "BE": 9.012183065,
    "B": 10.01293695,
    "C": 12.0,
    "N": 14.0030740044,
    "O": 15.9949146196,
    "F": 18.9984031627,
    "NE": 19.9924401762,
    "NA": 22.989769282,
    "MG": 23.985041697,
    "AL": 26.98153853,
    "SI": 27.9769265347,
    "P": 30.9737619984,
    "S": 31.9720711744,
    "CL": 34.968852682,
    "AR": 35.967545105,
    "K": 38.9637064864,
    "CA": 39.962590863,
    "SC": 44.95590828,
    "TI": 45.95262772,
    "V": 49.94715601,
    "CR": 49.94604183,
    "MN": 54.93804391,
    "FE": 53.93960899,
    "CO": 58.93319429,
    "NI": 57.93534241,
    "CU": 62.92959772,
    "ZN": 63.92914201,
    "GA": 68.9255735,
    "GE": 69.92424875,
    "AS": 74.92159457,
    "SE": 73.922475934,
    "BR": 78.9183376,
    "KR": 77.92036494,
    "RB": 84.9117897379,
    "SR": 83.9134191,
    "Y": 88.9058403,
    "ZR": 89.9046977,
    "NB": 92.906373,
    "MO": 91.90680796,
    "TC": 96.9063667,
    "RU": 95.90759025,
    "RH": 102.905498,
    "PD": 101.9056022,
    "AG": 106.9050916,
    "CD": 105.9064599,
    "IN": 112.90406184,
    "SN": 111.90482387,
    "SB": 120.903812,
    "TE": 119.9040593,
    "I": 126.9044719,
    "XE": 123.905892,
    "CS": 132.905451961,
    "BA": 129.9063207,
    "LA": 137.9071149,
    "CE": 135.90712921,
    "PR": 140.9076576,
    "ND": 141.907729,
    "PM": 144.9127559,
    "SM": 143.9120065,
    "EU": 150.9198578,
    "GD": 151.9197995,
    "TB": 158.9253547,
    "DY": 155.9242847,
    "HO": 164.9303288,
    "ER": 161.9287884,
    "TM": 168.9342179,
    "YB": 167.9338896,
    "LU": 174.9407752,
    "HF": 173.9400461,
    "TA": 179.9474648,
    "W": 179.9467108,
    "RE": 184.9529545,
    "OS": 183.9524885,
    "IR": 190.9605893,
    "PT": 189.9599297,
    "AU": 196.96656879,
    "HG": 195.9658326,
    "TL": 202.9723446,
    "PB": 203.973044,
    "BI": 208.9803991,
    "PO": 208.9824308,
    "AT": 209.9871479,
    "RN": 210.9906011,
    "FR": 223.019736,
    "RA": 223.0185023,
    "AC": 227.0277523,
    "TH": 230.0331341,
    "PA": 231.0358842,
    "U": 233.0396355,
    "NP": 236.04657,
    "PU": 238.0495601,
    "AM": 241.0568293,
    "CM": 243.0613893,
    "BK": 247.0703073,
    "CF": 249.0748539,
    "ES": 252.08298,
    "FM": 257.0951061,
    "MD": 258.0984315,
    "NO": 259.10103,
    "LR": 262.10961,
    "RF": 267.12179,
    "DB": 268.12567,
    "SG": 271.13393,
    "BH": 272.13826,
    "HS": 270.13429,
    "MT": 276.15159,
    "DS": 281.16451,
    "RG": 280.16514,
    "CN": 285.17712,
    "NH": 284.17873,
    "FL": 289.19042,
    "MC": 288.19274,
    "LV": 293.20449,
    "TS": 292.20746,
    "OG": 294.21392,
}

extended_mass_table = {
    "H": [
        (1.00782503223, 0.999885),
        (2.01410177812, 0.000115),
        (3.0160492779, 0.0),
    ],
    "HE": [
        (4.00260325413, 0.99999866),
        (3.0160293201, 1.34e-06),
    ],
    "LI": [
        (7.0160034366, 0.9241),
        (6.0151228874, 0.0759),
    ],
    "BE": [
        (9.012183065, 1.0),
    ],
    "B": [
        (11.00930536, 0.801),
        (10.01293695, 0.199),
    ],
    "C": [
        (12.0, 0.9893),
        (13.0033548351, 0.0107),
        (14.0032419884, 0.0),
    ],
    "N": [
        (14.0030740044, 0.99636),
        (15.0001088989, 0.00364),
    ],
    "O": [
        (15.9949146196, 0.99757),
        (17.9991596129, 0.00205),
        (16.9991317565, 0.00038),
    ],
    "F": [
        (18.9984031627, 1.0),
    ],
    "NE": [
        (19.9924401762, 0.9048),
        (21.991385114, 0.0925),
        (20.993846685, 0.0027),
    ],
    "NA": [
        (22.989769282, 1.0),
    ],
    "MG": [
        (23.985041697, 0.7899),
        (25.982592968, 0.1101),
        (24.985836976, 0.1),
    ],
    "AL": [
        (26.98153853, 1.0),
    ],
    "SI": [
        (27.9769265347, 0.92223),
        (28.9764946649, 0.04685),
        (29.973770136, 0.03092),
    ],
    "P": [
        (30.9737619984, 1.0),
    ],
    "S": [
        (31.9720711744, 0.9499),
        (33.967867004, 0.0425),
        (32.9714589098, 0.0075),
        (35.96708071, 0.0001),
    ],
    "CL": [
        (34.968852682, 0.7576),
        (36.965902602, 0.2424),
    ],
    "AR": [
        (39.9623831237, 0.996035),
        (35.967545105, 0.003336),
        (37.96273211, 0.000629),
    ],
    "K": [
        (38.9637064864, 0.932581),
        (40.9618252579, 0.067302),
        (39.963998166, 0.000117),
    ],
    "CA": [
        (39.962590863, 0.96941),
        (43.95548156, 0.02086),
        (41.95861783, 0.00647),
        (47.95252276, 0.00187),
        (42.95876644, 0.00135),
        (45.953689, 4e-05),
    ],
    "SC": [
        (44.95590828, 1.0),
    ],
    "TI": [
        (47.94794198, 0.7372),
        (45.95262772, 0.0825),
        (46.95175879, 0.0744),
        (48.94786568, 0.0541),
        (49.94478689, 0.0518),
    ],
    "V": [
        (50.94395704, 0.9975),
        (49.94715601, 0.0025),
    ],
    "CR": [
        (51.94050623, 0.83789),
        (52.94064815, 0.09501),
        (49.94604183, 0.04345),
        (53.93887916, 0.02365),
    ],
    "MN": [
        (54.93804391, 1.0),
    ],
    "FE": [
        (55.93493633, 0.91754),
        (53.93960899, 0.05845),
        (56.93539284, 0.02119),
        (57.93327443, 0.00282),
    ],
    "CO": [
        (58.93319429, 1.0),
    ],
    "NI": [
        (57.93534241, 0.68077),
        (59.93078588, 0.26223),
        (61.92834537, 0.036346),
        (60.93105557, 0.011399),
        (63.92796682, 0.009255),
    ],
    "CU": [
        (62.92959772, 0.6915),
        (64.9277897, 0.3085),
    ],
    "ZN": [
        (63.92914201, 0.4917),
        (65.92603381, 0.2773),
        (67.92484455, 0.1845),
        (66.92712775, 0.0404),
        (69.9253192, 0.0061),
    ],
    "GA": [
        (68.9255735, 0.60108),
        (70.92470258, 0.39892),
    ],
    "GE": [
        (73.921177761, 0.365),
        (71.922075826, 0.2745),
        (69.92424875, 0.2057),
        (72.923458956, 0.0775),
        (75.921402726, 0.0773),
    ],
    "AS": [
        (74.92159457, 1.0),
    ],
    "SE": [
        (79.9165218, 0.4961),
        (77.91730928, 0.2377),
        (75.919213704, 0.0937),
        (81.9166995, 0.0873),
        (76.919914154, 0.0763),
        (73.922475934, 0.0089),
    ],
    "BR": [
        (78.9183376, 0.5069),
        (80.9162897, 0.4931),
    ],
    "KR": [
        (83.9114977282, 0.56987),
        (85.9106106269, 0.17279),
        (81.91348273, 0.11593),
        (82.91412716, 0.115),
        (79.91637808, 0.02286),
        (77.92036494, 0.00355),
    ],
    "RB": [
        (84.9117897379, 0.7217),
        (86.909180531, 0.2783),
    ],
    "SR": [
        (87.9056125, 0.8258),
        (85.9092606, 0.0986),
        (86.9088775, 0.07),
        (83.9134191, 0.0056),
    ],
    "Y": [
        (88.9058403, 1.0),
    ],
    "ZR": [
        (89.9046977, 0.5145),
        (93.9063108, 0.1738),
        (91.9050347, 0.1715),
        (90.9056396, 0.1122),
        (95.9082714, 0.028),
    ],
    "NB": [
        (92.906373, 1.0),
    ],
    "MO": [
        (97.90540482, 0.2439),
        (95.90467612, 0.1667),
        (94.90583877, 0.1584),
        (91.90680796, 0.1453),
        (99.9074718, 0.0982),
        (96.90601812, 0.096),
        (93.9050849, 0.0915),
    ],
    "TC": [
        (96.9063667, 1.0),
        (98.9062508, 0.0),
        (97.9072124, 0.0),
    ],
    "RU": [
        (101.9043441, 0.3155),
        (103.9054275, 0.1862),
        (100.9055769, 0.1706),
        (98.9059341, 0.1276),
        (99.9042143, 0.126),
        (95.90759025, 0.0554),
        (97.9052868, 0.0187),
    ],
    "RH": [
        (102.905498, 1.0),
    ],
    "PD": [
        (105.9034804, 0.2733),
        (107.9038916, 0.2646),
        (104.9050796, 0.2233),
        (109.9051722, 0.1172),
        (103.9040305, 0.1114),
        (101.9056022, 0.0102),
    ],
    "AG": [
        (106.9050916, 0.51839),
        (108.9047553, 0.48161),
    ],
    "CD": [
        (113.90336509, 0.2873),
        (111.90276287, 0.2413),
        (110.90418287, 0.128),
        (109.90300661, 0.1249),
        (112.90440813, 0.1222),
        (115.90476315, 0.0749),
        (105.9064599, 0.0125),
        (107.9041834, 0.0089),
    ],
    "IN": [
        (114.903878776, 0.9571),
        (112.90406184, 0.0429),
    ],
    "SN": [
        (119.90220163, 0.3258),
        (117.90160657, 0.2422),
        (115.9017428, 0.1454),
        (118.90331117, 0.0859),
        (116.90295398, 0.0768),
        (123.9052766, 0.0579),
        (121.9034438, 0.0463),
        (111.90482387, 0.0097),
        (113.9027827, 0.0066),
        (114.903344699, 0.0034),
    ],
    "SB": [
        (120.903812, 0.5721),
        (122.9042132, 0.4279),
    ],
    "TE": [
        (129.906222748, 0.3408),
        (127.90446128, 0.3174),
        (125.9033109, 0.1884),
        (124.9044299, 0.0707),
        (123.9028171, 0.0474),
        (121.9030435, 0.0255),
        (122.9042698, 0.0089),
        (119.9040593, 0.0009),
    ],
    "I": [
        (126.9044719, 1.0),
    ],
    "XE": [
        (131.904155086, 0.269086),
        (128.904780861, 0.264006),
        (130.90508406, 0.212324),
        (133.90539466, 0.104357),
        (135.907214484, 0.088573),
        (129.903509349, 0.04071),
        (127.903531, 0.019102),
        (123.905892, 0.000952),
        (125.9042983, 0.00089),
    ],
    "CS": [
        (132.905451961, 1.0),
    ],
    "BA": [
        (137.905247, 0.71698),
        (136.90582714, 0.11232),
        (135.90457573, 0.07854),
        (134.90568838, 0.06592),
        (133.90450818, 0.02417),
        (129.9063207, 0.00106),
        (131.9050611, 0.00101),
    ],
    "LA": [
        (138.9063563, 0.9991119),
        (137.9071149, 0.0008881),
    ],
    "CE": [
        (139.9054431, 0.8845),
        (141.9092504, 0.11114),
        (137.905991, 0.00251),
        (135.90712921, 0.00185),
    ],
    "PR": [
        (140.9076576, 1.0),
    ],
    "ND": [
        (141.907729, 0.27152),
        (143.910093, 0.23798),
        (145.9131226, 0.17189),
        (142.90982, 0.12174),
        (144.9125793, 0.08293),
        (147.9168993, 0.05756),
        (149.9209022, 0.05638),
    ],
    "PM": [
        (144.9127559, 1.0),
        (146.915145, 0.0),
    ],
    "SM": [
        (151.9197397, 0.2675),
        (153.9222169, 0.2275),
        (146.9149044, 0.1499),
        (148.9171921, 0.1382),
        (147.9148292, 0.1124),
        (149.9172829, 0.0738),
        (143.9120065, 0.0307),
    ],
    "EU": [
        (152.921238, 0.5219),
        (150.9198578, 0.4781),
    ],
    "GD": [
        (157.9241123, 0.2484),
        (159.9270624, 0.2186),
        (155.9221312, 0.2047),
        (156.9239686, 0.1565),
        (154.9226305, 0.148),
        (153.9208741, 0.0218),
        (151.9197995, 0.002),
    ],
    "TB": [
        (158.9253547, 1.0),
    ],
    "DY": [
        (163.9291819, 0.2826),
        (161.9268056, 0.25475),
        (162.9287383, 0.24896),
        (160.9269405, 0.18889),
        (159.9252046, 0.02329),
        (157.9244159, 0.00095),
        (155.9242847, 0.00056),
    ],
    "HO": [
        (164.9303288, 1.0),
    ],
    "ER": [
        (165.9302995, 0.33503),
        (167.9323767, 0.26978),
        (166.9320546, 0.22869),
        (169.9354702, 0.1491),
        (163.9292088, 0.01601),
        (161.9287884, 0.00139),
    ],
    "TM": [
        (168.9342179, 1.0),
    ],
    "YB": [
        (173.9388664, 0.32026),
        (171.9363859, 0.2168),
        (172.9382151, 0.16103),
        (170.9363302, 0.1409),
        (175.9425764, 0.12996),
        (169.9347664, 0.02982),
        (167.9338896, 0.00123),
    ],
    "LU": [
        (174.9407752, 0.97401),
        (175.9426897, 0.02599),
    ],
    "HF": [
        (179.946557, 0.3508),
        (177.9437058, 0.2728),
        (176.9432277, 0.186),
        (178.9458232, 0.1362),
        (175.9414076, 0.0526),
        (173.9400461, 0.0016),
    ],
    "TA": [
        (180.9479958, 0.9998799),
        (179.9474648, 0.0001201),
    ],
    "W": [
        (183.95093092, 0.3064),
        (185.9543628, 0.2843),
        (181.94820394, 0.265),
        (182.95022275, 0.1431),
        (179.9467108, 0.0012),
    ],
    "RE": [
        (186.9557501, 0.626),
        (184.9529545, 0.374),
    ],
    "OS": [
        (191.961477, 0.4078),
        (189.9584437, 0.2626),
        (188.9581442, 0.1615),
        (187.9558352, 0.1324),
        (186.9557474, 0.0196),
        (185.953835, 0.0159),
        (183.9524885, 0.0002),
    ],
    "IR": [
        (192.9629216, 0.627),
        (190.9605893, 0.373),
    ],
    "PT": [
        (194.9647917, 0.3378),
        (193.9626809, 0.3286),
        (195.96495209, 0.2521),
        (197.9678949, 0.07356),
        (191.9610387, 0.00782),
        (189.9599297, 0.00012),
    ],
    "AU": [
        (196.96656879, 1.0),
    ],
    "HG": [
        (201.9706434, 0.2986),
        (199.96832659, 0.231),
        (198.96828064, 0.1687),
        (200.97030284, 0.1318),
        (197.9667686, 0.0997),
        (203.97349398, 0.0687),
        (195.9658326, 0.0015),
    ],
    "TL": [
        (204.9744278, 0.7048),
        (202.9723446, 0.2952),
    ],
    "PB": [
        (207.9766525, 0.524),
        (205.9744657, 0.241),
        (206.9758973, 0.221),
        (203.973044, 0.014),
    ],
    "BI": [
        (208.9803991, 1.0),
    ],
    "PO": [
        (208.9824308, 1.0),
        (209.9828741, 0.0),
    ],
    "AT": [
        (209.9871479, 1.0),
        (210.9874966, 0.0),
    ],
    "RN": [
        (210.9906011, 1.0),
        (222.0175782, 0.0),
        (220.0113941, 0.0),
    ],
    "FR": [
        (223.019736, 1.0),
    ],
    "RA": [
        (223.0185023, 1.0),
        (228.0310707, 0.0),
        (226.0254103, 0.0),
        (224.020212, 0.0),
    ],
    "AC": [
        (227.0277523, 1.0),
    ],
    "TH": [
        (230.0331341, 1.0),
        (232.0380558, 0.0),
    ],
    "PA": [
        (231.0358842, 1.0),
    ],
    "U": [
        (235.0439301, 0.992742),
        (234.0409523, 0.007204),
        (233.0396355, 5.4e-05),
        (238.0507884, 0.0),
        (236.0455682, 0.0),
    ],
    "NP": [
        (236.04657, 1.0),
        (237.0481736, 0.0),
    ],
    "PU": [
        (238.0495601, 1.0),
        (244.0642053, 0.0),
        (242.0587428, 0.0),
        (241.0568517, 0.0),
        (240.0538138, 0.0),
        (239.0521636, 0.0),
    ],
    "AM": [
        (241.0568293, 1.0),
        (243.0613813, 0.0),
    ],
    "CM": [
        (243.0613893, 1.0),
        (248.0723499, 0.0),
        (247.0703541, 0.0),
        (246.0672238, 0.0),
        (245.0654915, 0.0),
        (244.0627528, 0.0),
    ],
    "BK": [
        (247.0703073, 1.0),
        (249.0749877, 0.0),
    ],
    "CF": [
        (249.0748539, 1.0),
        (252.0816272, 0.0),
        (251.0795886, 0.0),
        (250.0764062, 0.0),
    ],
    "ES": [
        (252.08298, 1.0),
    ],
    "FM": [
        (257.0951061, 1.0),
    ],
    "MD": [
        (258.0984315, 1.0),
        (260.10365, 0.0),
    ],
    "NO": [
        (259.10103, 1.0),
    ],
    "LR": [
        (262.10961, 1.0),
    ],
    "RF": [
        (267.12179, 1.0),
    ],
    "DB": [
        (268.12567, 1.0),
    ],
    "SG": [
        (271.13393, 1.0),
    ],
    "BH": [
        (272.13826, 1.0),
    ],
    "HS": [
        (270.13429, 1.0),
    ],
    "MT": [
        (276.15159, 1.0),
    ],
    "DS": [
        (281.16451, 1.0),
    ],
    "RG": [
        (280.16514, 1.0),
    ],
    "CN": [
        (285.17712, 1.0),
    ],
    "NH": [
        (284.17873, 1.0),
    ],
    "FL": [
        (289.19042, 1.0),
    ],
    "MC": [
        (288.19274, 1.0),
    ],
    "LV": [
        (293.20449, 1.0),
    ],
    "TS": [
        (292.20746, 1.0),
    ],
}


# Convert masses in mass tables to a.u.
mass_table = {k: v * units.units["au_per_amu"] for k, v in list(mass_table.items())}

for k, v in list(extended_mass_table.items()):
    new_v = []
    for val in v:
        val2 = val[0] * units.units["au_per_amu"]
        new_v.append((val2, val[1]))
    extended_mass_table[k] = new_v


def from_Ns_to_widths(Ns):
    """
    Colton is a sweet boi.
    From the code of FMS
    case(‘H’,‘h’)      = 4.7d0
    case(‘D’,‘d’)      = 6.3d0 # TODO what is the number of this?
    case(‘T’,‘t’)      = 6.3d0 # TODO what is the number of this?
    case(‘C’,‘c’)      = 22.7d0
    case(‘N’,‘n’)      = 19.0d0
    case(‘O’,‘o’)      = 12.2d0
    case(‘S’,‘s’)      = 16.7d0
    case(‘F’,‘f’)      = 8.5d0
    """
    d = {
        1: 4.7,
        6: 22.7,
        7: 9.0,
        8: 12.2,
        9: 8.5,
        16: 16.7,
    }
    return [d[x] for x in Ns]
