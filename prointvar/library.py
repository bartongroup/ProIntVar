#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""

This is where dictionaries and lists defining global-like
parameters or variables are defined. I will try to added
references and/or dates to every entry as a commented header.

There are several entries used nowhere, but here for future
reference or simply deprecated.

Fábio Madeira, 2015+

"""

_mmcif_types = {
    'group_PDB': str,
    'id': int,
    'type_symbol': str,
    'label_atom_id': str,
    'label_alt_id': str,
    'label_comp_id': str,
    'label_asym_id': str,
    'label_entity_id': str,
    'label_seq_id': str,
    'new_asym_id': str,
    'new_seq_id': str,
    'pdbx_PDB_ins_code': str,
    'Cartn_x': float,
    'Cartn_y': float,
    'Cartn_z': float,
    'occupancy': float,
    'B_iso_or_equiv': float,
    'Cartn_x_esd': float,
    'Cartn_y_esd': float,
    'Cartn_z_esd': float,
    'occupancy_esd': float,
    'B_iso_or_equiv_esd': float,
    'pdbx_formal_charge': int,
    'auth_seq_id': str,
    'auth_comp_id': str,
    'auth_asym_id': str,
    'auth_atom_id': str,
    'pdbx_PDB_model_num': str,
    'pdbe_label_seq_id': str,
    'orig_label_asym_id': str,
    'orig_auth_asym_id': str,
    'auth_seq_id_full': str,
    'label_seq_id_full': str,
    'contact_indexes': str,
}

_dssp_types = {
    'LINE': int,
    'RES': str,
    'CHAIN': str,
    'CHAIN_FULL': str,
    'AA': str,
    'SS': str,
    'SS_CLASS': str,
    'STRUCTURE': str,
    'BP1': str,
    'BP2': str,
    'BP2_CHAIN': str,
    'ACC': int,
    'RSA': float,
    'RSA_class': str,
    'NH_O_1': int,
    'NH_O_1_nrg': float,
    'O_HN_1': int,
    'O_HN_1_nrg': float,
    'NH_O_2': int,
    'NH_O_2_nrg': float,
    'O_HN_2': int,
    'O_HN_2_nrg': float,
    'TCO': float,
    'KAPPA': float,
    'ALPHA': float,
    'PHI': float,
    'PSI': float,
    'X-CA': float,
    'Y-CA': float,
    'Z-CA': float
}

_sifts_types = {
    'PDB_regionId': int,
    'PDB_regionStart': int,
    'PDB_regionEnd': int,
    'PDB_regionResNum': str,
    'PDB_dbVersion': str,
    'PDB_dbAccessionId': str,
    'PDB_dbResNum': str,
    'PDB_dbResName': str,
    'PDB_dbChainId': str,
    'PDB_Annotation': str,
    'PDB_entityId': str,
    'PDB_codeSecondaryStructure': str,
    'PDB_nameSecondaryStructure': str,
    'UniProt_regionId': int,
    'UniProt_regionStart': int,
    'UniProt_regionEnd': int,
    'UniProt_regionResNum': str,
    'UniProt_dbVersion': str,
    'UniProt_dbAccessionId': str,
    'UniProt_dbResNum': str,
    'UniProt_dbResName': str,
    'CATH_regionId': int,
    'CATH_regionStart': int,
    'CATH_regionEnd': int,
    'CATH_regionResNum': str,
    'CATH_dbVersion': str,
    'CATH_dbAccessionId': str,
    'SCOP_regionId': int,
    'SCOP_regionStart': int,
    'SCOP_regionEnd': int,
    'SCOP_regionResNum': str,
    'SCOP_dbVersion': str,
    'SCOP_dbAccessionId': str,
    'Pfam_regionId': int,
    'Pfam_regionStart': int,
    'Pfam_regionEnd': int,
    'Pfam_regionResNum': str,
    'Pfam_dbVersion': str,
    'Pfam_dbAccessionId': str
}

_hbplus_types = {
    "CHAIN_D": str,
    "RES_D": str,
    "INSCODE_D": str,
    "ATOM_D": str,
    "CHAIN_A": str,
    "RES_A": str,
    "INSCODE_A": str,
    "ATOM_A": str,
    "DIST_DA": float,
    "CATEG_DA": str,
    "NUM_AAS": int,
    "DIST_CA-CA": float,
    "ANGLE_DHA": float,
    "DIST_H-A": float,
    "ANGLE_H-A-AA": float,
    "ANGLE_D-A-AA": float,
    "ID": int
}

_arpeggio_types = {
    "ENTRY_A": str,
    "CHAIN_A": str,
    "RES_A": str,
    "INSCODE_A": str,
    "RES_FULL_A": str,
    "ATOM_A": str,
    "ENTRY_B": str,
    "CHAIN_B": str,
    "RES_B": str,
    "INSCODE_B": str,
    "RES_FULL_B": str,
    "ATOM_B": str,
    "STERIC_CLASH": int,
    "COVALENT": int,
    "VDW_CLASH": int,
    "VDW_INTER": int,
    "PROXIMAL": int,
    "HYDROGEN": int,
    "WEAK_HYDROGEN": int,
    "HALOGEN": int,
    "IONIC": int,
    "METAL_COMPLEX": int,
    "AROMATIC": int,
    "HYDROPHOBIC": int,
    "CARBONYL": int,
    "POLAR": int,
    "WEAK_POLAR": int,
    "DIST": float,
    "VDW_DIST": float,
    "ENTITIES": str
}

_probe_types = {}

_stamp_types = {
    "Domain1": str,
    "Domain2": str,
    "Fits": int,
    "Sc": float,
    "RMS": float,
    "A_Len": int,
    "B_Len": int,
    "Align_Len": int,
    "N_Fit": int,
    "N_Equiv": int,
    "N_SS_Equiv": int,
    "PID": float,
    "SS_PID": float,
    "Pm": float,
}

_uni_ens_var_types = {
    'begin': int,
    'end': int,
    'polyphenScore': float,
    'siftScore': float,
}

_dtypes_convert = {
    int: 'int64',
    float: 'float64',
    str: 'object'
}

mmcif_types = {k: _dtypes_convert[v] for k, v in _mmcif_types.items()}
dssp_types = {k: _dtypes_convert[v] for k, v in _dssp_types.items()}
sifts_types = {k: _dtypes_convert[v] for k, v in _sifts_types.items()}
hbplus_types = {k: _dtypes_convert[v] for k, v in _hbplus_types.items()}
arpeggio_types = {k: _dtypes_convert[v] for k, v in _arpeggio_types.items()}
probe_types = {k: _dtypes_convert[v] for k, v in _probe_types.items()}
stamp_types = {k: _dtypes_convert[v] for k, v in _stamp_types.items()}
uni_ens_var_types = {k: _dtypes_convert[v] for k, v in _uni_ens_var_types.items()}

arpeggio_col_renames = {
    "STERIC_CLASH": "Steric-Clash",
    "COVALENT": "Covalent-Bond",
    "VDW_CLASH": "VDW-Clash",
    "VDW_INTER": "VDW-Bond",
    "PROXIMAL": "VDW-Proximal",
    "HYDROGEN": "Hydrogen-Bond",
    "WEAK_HYDROGEN": "Weak-Hydrogen-Bond",
    "HALOGEN": "Halogen-Bond",
    "IONIC": "Ionic-Bond",
    "METAL_COMPLEX": "Metal-Complex",
    "AROMATIC": "Aromatic-Bond",
    "HYDROPHOBIC": "Hydrophobic-Bond",
    "CARBONYL": "Carbonyl-Bond",
    "POLAR": "Polar-Bond",
    "WEAK_POLAR": "Weak-Polar-Bond",
    "Amide-Amide": "Amide-Amide",
    "Amide-Aromatic": "Amide-Aromatic",
    "Aromatic-Aromatic": "Aromatic-Aromatic",
    "Atom-Ring": "Atom-Ring",
}

# working species in Ensembl Variants, as of November 2014
# based on ftp://ftp.ensembl.org/pub/release-77/variation/vcf/
ensembl_species = [
    "bos_taurus",
    "canis_familiaris",
    "danio_rerio",
    "drosophila_melanogaster",
    "equus_caballus",
    "felis_catus",
    "gallus_gallus",
    "homo_sapiens",
    "macaca_mulatta",
    "meleagris_gallopavo",
    "monodelphis_domestica",
    "mus_musculus",
    "nomascus_leucogenys",
    "ornithorhynchus_anatinus",
    "ovis_aries",
    "pan_troglodytes",
    "pongo_abelii",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
    "sus_scrofa",
    "taeniopygia_guttata",
    "tetraodon_nigroviridis"
]

ensembl_species_proteome = [
    "bos_taurus",
    "canis_familiaris",
    "danio_rerio",
    "drosophila_melanogaster",
    # "equus_caballus",
    # "felis_catus",
    "gallus_gallus",
    "homo_sapiens",
    # "macaca_mulatta",
    # "meleagris_gallopavo",
    # "monodelphis_domestica",
    "mus_musculus",
    # "nomascus_leucogenys",
    # "ornithorhynchus_anatinus",
    # "ovis_aries",
    # "pan_troglodytes",
    # "pongo_abelii",
    "rattus_norvegicus",
    "saccharomyces_cerevisiae",
    "sus_scrofa",
    # "taeniopygia_guttata",
    # "tetraodon_nigroviridis"]
]

# variant types and impact
# http://www.ensembl.org/info/genome/variation/predicted_data.html
ensembl_variant_types = {
    'transcript_ablation': 'HIGH',
    'splice_acceptor_variant': 'HIGH',
    'splice_donor_variant': 'HIGH',
    'stop_gained': 'HIGH',
    'frameshift_variant': 'HIGH',
    'stop_lost': 'HIGH',
    'start_lost': 'HIGH',
    'transcript_amplification': 'HIGH',
    'inframe_insertion': 'MODERATE',
    'inframe_deletion': 'MODERATE',
    'missense_variant': 'MODERATE',
    'protein_altering_variant': 'MODERATE',
    'splice_region_variant': 'LOW',
    'incomplete_terminal_codon_variant': 'LOW',
    'stop_retained_variant': 'LOW',
    'synonymous_variant': 'LOW',
    'coding_sequence_variant': 'MODIFIER',
    'mature_miRNA_variant': 'MODIFIER',
    '5_prime_UTR_variant': 'MODIFIER',
    '3_prime_UTR_variant': 'MODIFIER',
    'non_coding_transcript_exon_variant': 'MODIFIER',
    'intron_variant': 'MODIFIER',
    'NMD_transcript_variant': 'MODIFIER',
    'non_coding_transcript_variant': 'MODIFIER',
    'upstream_gene_variant': 'MODIFIER',
    'downstream_gene_variant': 'MODIFIER',
    'TFBS_ablation': 'MODIFIER',
    'TFBS_amplification': 'MODIFIER',
    'TF_binding_site_variant': 'MODIFIER',
    'regulatory_region_ablation': 'MODERATE',
    'regulatory_region_amplification': 'MODIFIER',
    'feature_elongation': 'MODIFIER',
    'regulatory_region_variant': 'MODIFIER',
    'feature_truncation': 'MODIFIER',
    'intergenic_variant': 'MODIFIER'
}

# updating terms in ensembl output so that they match
# UniProt Proteins API counterparts
update_ensembl_to_uniprot = {
    'minor_allele_frequency': 'frequency',
    'start': 'begin',
    'end': 'end',
    'sift': 'siftScore',
    'polyphen': 'polyphenScore',
    'type': 'consequenceType',
    'id': 'xrefs_id'
}

# Default valid protein residues in single-letter alphabet
aa_codes_1to3_common = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP',
    'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
    'H': 'HIS', 'K': 'LYS', 'I': 'ILE',
    'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL',
    'Y': 'TYR', 'W': 'TRP'
}

# Asx	B	Aspartic acid or Asparagine
# Glx	Z	Glutamic acid or Glutamine
# Xaa	X	Any amino acid
# Xle	J	Leucine or Isoleucine
aa_codes_1to3_extended = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP',
    'E': 'GLU', 'F': 'PHE', 'G': 'GLY',
    'H': 'HIS', 'K': 'LYS', 'I': 'ILE',
    'L': 'LEU', 'M': 'MET', 'N': 'ASN',
    'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL',
    'Y': 'TYR', 'W': 'TRP',
    'X': 'LNT', 'B': 'ASX', 'Z': 'GLX',
    'J': 'XLE', 'U': 'SEC', 'O': 'PYL',
    '-': '---',
}

aa_codes_3to1_common = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D',
    'GLU': 'E', 'PHE': 'F', 'GLY': 'G',
    'HIS': 'H', 'LYS': 'K', 'ILE': 'I',
    'LEU': 'L', 'MET': 'M', 'ASN': 'N',
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V',
    'TYR': 'Y', 'TRP': 'W'
}

# obtained from biopython at
# https://github.com/biopython/biopython/blob/master/Bio/Data/SCOPData.py
aa_codes_3to1_extended = {
    "00C": "C", "01W": "X", "02K": "A", "03Y": "C", "07O": "C",
    "08P": "C", "0A0": "D", "0A1": "Y", "0A2": "K", "0A8": "C",
    "0AA": "V", "0AB": "V", "0AC": "G", "0AD": "G", "0AF": "W",
    "0AG": "L", "0AH": "S", "0AK": "D", "0AM": "A", "0AP": "C",
    "0AU": "U", "0AV": "A", "0AZ": "P", "0BN": "F", "0C ": "C",
    "0CS": "A", "0DC": "C", "0DG": "G", "0DT": "T", "0FL": "A",
    "0G ": "G", "0NC": "A", "0SP": "A", "0U ": "U", "0YG": "YG",
    "10C": "C", "125": "U", "126": "U", "127": "U", "128": "N",
    "12A": "A", "143": "C", "175": "ASG", "193": "X", "1AP": "A",
    "1MA": "A", "1MG": "G", "1PA": "F", "1PI": "A", "1PR": "N",
    "1SC": "C", "1TQ": "W", "1TY": "Y", "1X6": "S", "200": "F",
    "23F": "F", "23S": "X", "26B": "T", "2AD": "X", "2AG": "A",
    "2AO": "X", "2AR": "A", "2AS": "X", "2AT": "T", "2AU": "U",
    "2BD": "I", "2BT": "T", "2BU": "A", "2CO": "C", "2DA": "A",
    "2DF": "N", "2DM": "N", "2DO": "X", "2DT": "T", "2EG": "G",
    "2FE": "N", "2FI": "N", "2FM": "M", "2GT": "T", "2HF": "H",
    "2LU": "L", "2MA": "A", "2MG": "G", "2ML": "L", "2MR": "R",
    "2MT": "P", "2MU": "U", "2NT": "T", "2OM": "U", "2OT": "T",
    "2PI": "X", "2PR": "G", "2SA": "N", "2SI": "X", "2ST": "T",
    "2TL": "T", "2TY": "Y", "2VA": "V", "2XA": "C", "32S": "X",
    "32T": "X", "3AH": "H", "3AR": "X", "3CF": "F", "3DA": "A",
    "3DR": "N", "3GA": "A", "3MD": "D", "3ME": "U", "3NF": "Y",
    "3QN": "K", "3TY": "X", "3XH": "G", "4AC": "N", "4BF": "Y",
    "4CF": "F", "4CY": "M", "4DP": "W", "4F3": "GYG", "4FB": "P",
    "4FW": "W", "4HT": "W", "4IN": "W", "4MF": "N", "4MM": "X",
    "4OC": "C", "4PC": "C", "4PD": "C", "4PE": "C", "4PH": "F",
    "4SC": "C", "4SU": "U", "4TA": "N", "4U7": "A", "56A": "H",
    "5AA": "A", "5AB": "A", "5AT": "T", "5BU": "U", "5CG": "G",
    "5CM": "C", "5CS": "C", "5FA": "A", "5FC": "C", "5FU": "U",
    "5HP": "E", "5HT": "T", "5HU": "U", "5IC": "C", "5IT": "T",
    "5IU": "U", "5MC": "C", "5MD": "N", "5MU": "U", "5NC": "C",
    "5PC": "C", "5PY": "T", "5SE": "U", "5ZA": "TWG", "64T": "T",
    "6CL": "K", "6CT": "T", "6CW": "W", "6HA": "A", "6HC": "C",
    "6HG": "G", "6HN": "K", "6HT": "T", "6IA": "A", "6MA": "A",
    "6MC": "A", "6MI": "N", "6MT": "A", "6MZ": "N", "6OG": "G",
    "70U": "U", "7DA": "A", "7GU": "G", "7JA": "I", "7MG": "G",
    "8AN": "A", "8FG": "G", "8MG": "G", "8OG": "G", "9NE": "E",
    "9NF": "F", "9NR": "R", "9NV": "V", "A  ": "A", "A1P": "N",
    "A23": "A", "A2L": "A", "A2M": "A", "A34": "A", "A35": "A",
    "A38": "A", "A39": "A", "A3A": "A", "A3P": "A", "A40": "A",
    "A43": "A", "A44": "A", "A47": "A", "A5L": "A", "A5M": "C",
    "A5N": "N", "A5O": "A", "A66": "X", "AA3": "A", "AA4": "A",
    "AAR": "R", "AB7": "X", "ABA": "A", "ABR": "A", "ABS": "A",
    "ABT": "N", "ACB": "D", "ACL": "R", "AD2": "A", "ADD": "X",
    "ADX": "N", "AEA": "X", "AEI": "D", "AET": "A", "AFA": "N",
    "AFF": "N", "AFG": "G", "AGM": "R", "AGT": "C", "AHB": "N",
    "AHH": "X", "AHO": "A", "AHP": "A", "AHS": "X", "AHT": "X",
    "AIB": "A", "AKL": "D", "AKZ": "D", "ALA": "A", "ALC": "A",
    "ALM": "A", "ALN": "A", "ALO": "T", "ALQ": "X", "ALS": "A",
    "ALT": "A", "ALV": "A", "ALY": "K", "AN8": "A", "AP7": "A",
    "APE": "X", "APH": "A", "API": "K", "APK": "K", "APM": "X",
    "APP": "X", "AR2": "R", "AR4": "E", "AR7": "R", "ARG": "R",
    "ARM": "R", "ARO": "R", "ARV": "X", "AS ": "A", "AS2": "D",
    "AS9": "X", "ASA": "D", "ASB": "D", "ASI": "D", "ASK": "D",
    "ASL": "D", "ASM": "X", "ASN": "N", "ASP": "D", "ASQ": "D",
    "ASU": "N", "ASX": "B", "ATD": "T", "ATL": "T", "ATM": "T",
    "AVC": "A", "AVN": "X", "AYA": "A", "AYG": "AYG", "AZK": "K",
    "AZS": "S", "AZY": "Y", "B1F": "F", "B1P": "N", "B2A": "A",
    "B2F": "F", "B2I": "I", "B2V": "V", "B3A": "A", "B3D": "D",
    "B3E": "E", "B3K": "K", "B3L": "X", "B3M": "X", "B3Q": "X",
    "B3S": "S", "B3T": "X", "B3U": "H", "B3X": "N", "B3Y": "Y",
    "BB6": "C", "BB7": "C", "BB8": "F", "BB9": "C", "BBC": "C",
    "BCS": "C", "BE2": "X", "BFD": "D", "BG1": "S", "BGM": "G",
    "BH2": "D", "BHD": "D", "BIF": "F", "BIL": "X", "BIU": "I",
    "BJH": "X", "BLE": "L", "BLY": "K", "BMP": "N", "BMT": "T",
    "BNN": "F", "BNO": "X", "BOE": "T", "BOR": "R", "BPE": "C",
    "BRU": "U", "BSE": "S", "BT5": "N", "BTA": "L", "BTC": "C",
    "BTR": "W", "BUC": "C", "BUG": "V", "BVP": "U", "BZG": "N",
    "C  ": "C", "C12": "TYG", "C1X": "K", "C25": "C", "C2L": "C",
    "C2S": "C", "C31": "C", "C32": "C", "C34": "C", "C36": "C",
    "C37": "C", "C38": "C", "C3Y": "C", "C42": "C", "C43": "C",
    "C45": "C", "C46": "C", "C49": "C", "C4R": "C", "C4S": "C",
    "C5C": "C", "C66": "X", "C6C": "C", "C99": "TFG", "CAF": "C",
    "CAL": "X", "CAR": "C", "CAS": "C", "CAV": "X", "CAY": "C",
    "CB2": "C", "CBR": "C", "CBV": "C", "CCC": "C", "CCL": "K",
    "CCS": "C", "CCY": "CYG", "CDE": "X", "CDV": "X", "CDW": "C",
    "CEA": "C", "CFL": "C", "CFY": "FCYG", "CG1": "G", "CGA": "E",
    "CGU": "E", "CH ": "C", "CH6": "MYG", "CH7": "KYG", "CHF": "X",
    "CHG": "X", "CHP": "G", "CHS": "X", "CIR": "R", "CJO": "GYG",
    "CLE": "L", "CLG": "K", "CLH": "K", "CLV": "AFG", "CM0": "N",
    "CME": "C", "CMH": "C", "CML": "C", "CMR": "C", "CMT": "C",
    "CNU": "U", "CP1": "C", "CPC": "X", "CPI": "X", "CQR": "GYG",
    "CR0": "TLG", "CR2": "GYG", "CR5": "G", "CR7": "KYG", "CR8": "HYG",
    "CRF": "TWG", "CRG": "THG", "CRK": "MYG", "CRO": "GYG", "CRQ": "QYG",
    "CRU": "EYG", "CRW": "ASG", "CRX": "ASG", "CS0": "C", "CS1": "C",
    "CS3": "C", "CS4": "C", "CS8": "N", "CSA": "C", "CSB": "C",
    "CSD": "C", "CSE": "C", "CSF": "C", "CSH": "SHG", "CSI": "G",
    "CSJ": "C", "CSL": "C", "CSO": "C", "CSP": "C", "CSR": "C",
    "CSS": "C", "CSU": "C", "CSW": "C", "CSX": "C", "CSY": "SYG",
    "CSZ": "C", "CTE": "W", "CTG": "T", "CTH": "T", "CUC": "X",
    "CWR": "S", "CXM": "M", "CY0": "C", "CY1": "C", "CY3": "C",
    "CY4": "C", "CYA": "C", "CYD": "C", "CYF": "C", "CYG": "C",
    "CYJ": "X", "CYM": "C", "CYQ": "C", "CYR": "C", "CYS": "C",
    "CZ2": "C", "CZO": "GYG", "CZZ": "C", "D11": "T", "D1P": "N",
    "D3 ": "N", "D33": "N", "D3P": "G", "D3T": "T", "D4M": "T",
    "D4P": "X", "DA ": "A", "DA2": "X", "DAB": "A", "DAH": "F",
    "DAL": "A", "DAR": "R", "DAS": "D", "DBB": "T", "DBM": "N",
    "DBS": "S", "DBU": "T", "DBY": "Y", "DBZ": "A", "DC ": "C",
    "DC2": "C", "DCG": "G", "DCI": "X", "DCL": "X", "DCT": "C",
    "DCY": "C", "DDE": "H", "DDG": "G", "DDN": "U", "DDX": "N",
    "DFC": "C", "DFG": "G", "DFI": "X", "DFO": "X", "DFT": "N",
    "DG ": "G", "DGH": "G", "DGI": "G", "DGL": "E", "DGN": "Q",
    "DHA": "S", "DHI": "H", "DHL": "X", "DHN": "V", "DHP": "X",
    "DHU": "U", "DHV": "V", "DI ": "I", "DIL": "I", "DIR": "R",
    "DIV": "V", "DLE": "L", "DLS": "K", "DLY": "K", "DM0": "K",
    "DMH": "N", "DMK": "D", "DMT": "X", "DN ": "N", "DNE": "L",
    "DNG": "L", "DNL": "K", "DNM": "L", "DNP": "A", "DNR": "C",
    "DNS": "K", "DOA": "X", "DOC": "C", "DOH": "D", "DON": "L",
    "DPB": "T", "DPH": "F", "DPL": "P", "DPP": "A", "DPQ": "Y",
    "DPR": "P", "DPY": "N", "DRM": "U", "DRP": "N", "DRT": "T",
    "DRZ": "N", "DSE": "S", "DSG": "N", "DSN": "S", "DSP": "D",
    "DT ": "T", "DTH": "T", "DTR": "W", "DTY": "Y", "DU ": "U",
    "DVA": "V", "DXD": "N", "DXN": "N", "DYG": "DYG", "DYS": "C",
    "DZM": "A", "E  ": "A", "E1X": "A", "ECC": "Q", "EDA": "A",
    "EFC": "C", "EHP": "F", "EIT": "T", "ENP": "N", "ESB": "Y",
    "ESC": "M", "EXB": "X", "EXY": "L", "EY5": "N", "EYS": "X",
    "F2F": "F", "FA2": "A", "FA5": "N", "FAG": "N", "FAI": "N",
    "FB5": "A", "FB6": "A", "FCL": "F", "FFD": "N", "FGA": "E",
    "FGL": "G", "FGP": "S", "FHL": "X", "FHO": "K", "FHU": "U",
    "FLA": "A", "FLE": "L", "FLT": "Y", "FME": "M", "FMG": "G",
    "FMU": "N", "FOE": "C", "FOX": "G", "FP9": "P", "FPA": "F",
    "FRD": "X", "FT6": "W", "FTR": "W", "FTY": "Y", "FVA": "V",
    "FZN": "K", "G  ": "G", "G25": "G", "G2L": "G", "G2S": "G",
    "G31": "G", "G32": "G", "G33": "G", "G36": "G", "G38": "G",
    "G42": "G", "G46": "G", "G47": "G", "G48": "G", "G49": "G",
    "G4P": "N", "G7M": "G", "GAO": "G", "GAU": "E", "GCK": "C",
    "GCM": "X", "GDP": "G", "GDR": "G", "GFL": "G", "GGL": "E",
    "GH3": "G", "GHG": "Q", "GHP": "G", "GL3": "G", "GLH": "Q",
    "GLJ": "E", "GLK": "E", "GLM": "X", "GLN": "Q", "GLQ": "E",
    "GLU": "E", "GLX": "Z", "GLY": "G", "GLZ": "G", "GMA": "E",
    "GMS": "G", "GMU": "U", "GN7": "G", "GND": "X", "GNE": "N",
    "GOM": "G", "GPL": "K", "GS ": "G", "GSC": "G", "GSR": "G",
    "GSS": "G", "GSU": "E", "GT9": "C", "GTP": "G", "GVL": "X",
    "GYC": "CYG", "GYS": "SYG", "H2U": "U", "H5M": "P", "HAC": "A",
    "HAR": "R", "HBN": "H", "HCS": "X", "HDP": "U", "HEU": "U",
    "HFA": "X", "HGL": "X", "HHI": "H", "HHK": "AK", "HIA": "H",
    "HIC": "H", "HIP": "H", "HIQ": "H", "HIS": "H", "HL2": "L",
    "HLU": "L", "HMR": "R", "HOL": "N", "HPC": "F", "HPE": "F",
    "HPH": "F", "HPQ": "F", "HQA": "A", "HRG": "R", "HRP": "W",
    "HS8": "H", "HS9": "H", "HSE": "S", "HSL": "S", "HSO": "H",
    "HTI": "C", "HTN": "N", "HTR": "W", "HV5": "A", "HVA": "V",
    "HY3": "P", "HYP": "P", "HZP": "P", "I  ": "I", "I2M": "I",
    "I58": "K", "I5C": "C", "IAM": "A", "IAR": "R", "IAS": "D",
    "IC ": "C", "IEL": "K", "IEY": "HYG", "IG ": "G", "IGL": "G",
    "IGU": "G", "IIC": "SHG", "IIL": "I", "ILE": "I", "ILG": "E",
    "ILX": "I", "IMC": "C", "IML": "I", "IOY": "F", "IPG": "G",
    "IPN": "N", "IRN": "N", "IT1": "K", "IU ": "U", "IYR": "Y",
    "IYT": "T", "IZO": "M", "JJJ": "C", "JJK": "C", "JJL": "C",
    "JW5": "N", "K1R": "C", "KAG": "G", "KCX": "K", "KGC": "K",
    "KNB": "A", "KOR": "M", "KPI": "K", "KST": "K", "KYQ": "K",
    "L2A": "X", "LA2": "K", "LAA": "D", "LAL": "A", "LBY": "K",
    "LC ": "C", "LCA": "A", "LCC": "N", "LCG": "G", "LCH": "N",
    "LCK": "K", "LCX": "K", "LDH": "K", "LED": "L", "LEF": "L",
    "LEH": "L", "LEI": "V", "LEM": "L", "LEN": "L", "LET": "X",
    "LEU": "L", "LEX": "L", "LG ": "G", "LGP": "G", "LHC": "X",
    "LHU": "U", "LKC": "N", "LLP": "K", "LLY": "K", "LME": "E",
    "LMF": "K", "LMQ": "Q", "LMS": "N", "LP6": "K", "LPD": "P",
    "LPG": "G", "LPL": "X", "LPS": "S", "LSO": "X", "LTA": "X",
    "LTR": "W", "LVG": "G", "LVN": "V", "LYF": "K", "LYK": "K",
    "LYM": "K", "LYN": "K", "LYR": "K", "LYS": "K", "LYX": "K",
    "LYZ": "K", "M0H": "C", "M1G": "G", "M2G": "G", "M2L": "K",
    "M2S": "M", "M30": "G", "M3L": "K", "M5M": "C", "MA ": "A",
    "MA6": "A", "MA7": "A", "MAA": "A", "MAD": "A", "MAI": "R",
    "MBQ": "Y", "MBZ": "N", "MC1": "S", "MCG": "X", "MCL": "K",
    "MCS": "C", "MCY": "C", "MD3": "C", "MD6": "G", "MDH": "X",
    "MDO": "ASG", "MDR": "N", "MEA": "F", "MED": "M", "MEG": "E",
    "MEN": "N", "MEP": "U", "MEQ": "Q", "MET": "M", "MEU": "G",
    "MF3": "X", "MFC": "GYG", "MG1": "G", "MGG": "R", "MGN": "Q",
    "MGQ": "A", "MGV": "G", "MGY": "G", "MHL": "L", "MHO": "M",
    "MHS": "H", "MIA": "A", "MIS": "S", "MK8": "L", "ML3": "K",
    "MLE": "L", "MLL": "L", "MLY": "K", "MLZ": "K", "MME": "M",
    "MMO": "R", "MMT": "T", "MND": "N", "MNL": "L", "MNU": "U",
    "MNV": "V", "MOD": "X", "MP8": "P", "MPH": "X", "MPJ": "X",
    "MPQ": "G", "MRG": "G", "MSA": "G", "MSE": "M", "MSL": "M",
    "MSO": "M", "MSP": "X", "MT2": "M", "MTR": "T", "MTU": "A",
    "MTY": "Y", "MVA": "V", "N  ": "N", "N10": "S", "N2C": "X",
    "N5I": "N", "N5M": "C", "N6G": "G", "N7P": "P", "NA8": "A",
    "NAL": "A", "NAM": "A", "NB8": "N", "NBQ": "Y", "NC1": "S",
    "NCB": "A", "NCX": "N", "NCY": "X", "NDF": "F", "NDN": "U",
    "NEM": "H", "NEP": "H", "NF2": "N", "NFA": "F", "NHL": "E",
    "NIT": "X", "NIY": "Y", "NLE": "L", "NLN": "L", "NLO": "L",
    "NLP": "L", "NLQ": "Q", "NMC": "G", "NMM": "R", "NMS": "T",
    "NMT": "T", "NNH": "R", "NP3": "N", "NPH": "C", "NPI": "A",
    "NRP": "LYG", "NRQ": "MYG", "NSK": "X", "NTY": "Y", "NVA": "V",
    "NYC": "TWG", "NYG": "NYG", "NYM": "N", "NYS": "C", "NZH": "H",
    "O12": "X", "O2C": "N", "O2G": "G", "OAD": "N", "OAS": "S",
    "OBF": "X", "OBS": "X", "OCS": "C", "OCY": "C", "ODP": "N",
    "OHI": "H", "OHS": "D", "OIC": "X", "OIP": "I", "OLE": "X",
    "OLT": "T", "OLZ": "S", "OMC": "C", "OMG": "G", "OMT": "M",
    "OMU": "U", "ONE": "U", "ONH": "A", "ONL": "X", "OPR": "R",
    "ORN": "A", "ORQ": "R", "OSE": "S", "OTB": "X", "OTH": "T",
    "OTY": "Y", "OXX": "D", "P  ": "G", "P1L": "C", "P1P": "N",
    "P2T": "T", "P2U": "U", "P2Y": "P", "P5P": "A", "PAQ": "Y",
    "PAS": "D", "PAT": "W", "PAU": "A", "PBB": "C", "PBF": "F",
    "PBT": "N", "PCA": "E", "PCC": "P", "PCE": "X", "PCS": "F",
    "PDL": "X", "PDU": "U", "PEC": "C", "PF5": "F", "PFF": "F",
    "PFX": "X", "PG1": "S", "PG7": "G", "PG9": "G", "PGL": "X",
    "PGN": "G", "PGP": "G", "PGY": "G", "PHA": "F", "PHD": "D",
    "PHE": "F", "PHI": "F", "PHL": "F", "PHM": "F", "PIA": "AYG",
    "PIV": "X", "PLE": "L", "PM3": "F", "PMT": "C", "POM": "P",
    "PPN": "F", "PPU": "A", "PPW": "G", "PQ1": "N", "PR3": "C",
    "PR5": "A", "PR9": "P", "PRN": "A", "PRO": "P", "PRS": "P",
    "PSA": "F", "PSH": "H", "PST": "T", "PSU": "U", "PSW": "C",
    "PTA": "X", "PTH": "Y", "PTM": "Y", "PTR": "Y", "PU ": "A",
    "PUY": "N", "PVH": "H", "PVL": "X", "PYA": "A", "PYO": "U",
    "PYX": "C", "PYY": "N", "QLG": "QLG", "QMM": "Q", "QPA": "C",
    "QPH": "F", "QUO": "G", "R  ": "A", "R1A": "C", "R4K": "W",
    "RC7": "HYG", "RE0": "W", "RE3": "W", "RIA": "A", "RMP": "A",
    "RON": "X", "RT ": "T", "RTP": "N", "S1H": "S", "S2C": "C",
    "S2D": "A", "S2M": "T", "S2P": "A", "S4A": "A", "S4C": "C",
    "S4G": "G", "S4U": "U", "S6G": "G", "SAC": "S", "SAH": "C",
    "SAR": "G", "SBL": "S", "SC ": "C", "SCH": "C", "SCS": "C",
    "SCY": "C", "SD2": "X", "SDG": "G", "SDP": "S", "SEB": "S",
    "SEC": "A", "SEG": "A", "SEL": "S", "SEM": "S", "SEN": "S",
    "SEP": "S", "SER": "S", "SET": "S", "SGB": "S", "SHC": "C",
    "SHP": "G", "SHR": "K", "SIB": "C", "SIC": "DC", "SLA": "P",
    "SLR": "P", "SLZ": "K", "SMC": "C", "SME": "M", "SMF": "F",
    "SMP": "A", "SMT": "T", "SNC": "C", "SNN": "N", "SOC": "C",
    "SOS": "N", "SOY": "S", "SPT": "T", "SRA": "A", "SSU": "U",
    "STY": "Y", "SUB": "X", "SUI": "DG", "SUN": "S", "SUR": "U",
    "SVA": "S", "SVV": "S", "SVW": "S", "SVX": "S", "SVY": "S",
    "SVZ": "X", "SWG": "SWG", "SYS": "C", "T  ": "T", "T11": "F",
    "T23": "T", "T2S": "T", "T2T": "N", "T31": "U", "T32": "T",
    "T36": "T", "T37": "T", "T38": "T", "T39": "T", "T3P": "T",
    "T41": "T", "T48": "T", "T49": "T", "T4S": "T", "T5O": "U",
    "T5S": "T", "T66": "X", "T6A": "A", "TA3": "T", "TA4": "X",
    "TAF": "T", "TAL": "N", "TAV": "D", "TBG": "V", "TBM": "T",
    "TC1": "C", "TCP": "T", "TCQ": "Y", "TCR": "W", "TCY": "A",
    "TDD": "L", "TDY": "T", "TFE": "T", "TFO": "A", "TFQ": "F",
    "TFT": "T", "TGP": "G", "TH6": "T", "THC": "T", "THO": "X",
    "THR": "T", "THX": "N", "THZ": "R", "TIH": "A", "TLB": "N",
    "TLC": "T", "TLN": "U", "TMB": "T", "TMD": "T", "TNB": "C",
    "TNR": "S", "TOX": "W", "TP1": "T", "TPC": "C", "TPG": "G",
    "TPH": "X", "TPL": "W", "TPO": "T", "TPQ": "Y", "TQI": "W",
    "TQQ": "W", "TRF": "W", "TRG": "K", "TRN": "W", "TRO": "W",
    "TRP": "W", "TRQ": "W", "TRW": "W", "TRX": "W", "TS ": "N",
    "TST": "X", "TT ": "N", "TTD": "T", "TTI": "U", "TTM": "T",
    "TTQ": "W", "TTS": "Y", "TY1": "Y", "TY2": "Y", "TY3": "Y",
    "TY5": "Y", "TYB": "Y", "TYI": "Y", "TYJ": "Y", "TYN": "Y",
    "TYO": "Y", "TYQ": "Y", "TYR": "Y", "TYS": "Y", "TYT": "Y",
    "TYU": "N", "TYW": "Y", "TYX": "X", "TYY": "Y", "TZB": "X",
    "TZO": "X", "U  ": "U", "U25": "U", "U2L": "U", "U2N": "U",
    "U2P": "U", "U31": "U", "U33": "U", "U34": "U", "U36": "U",
    "U37": "U", "U8U": "U", "UAR": "U", "UCL": "U", "UD5": "U",
    "UDP": "N", "UFP": "N", "UFR": "U", "UFT": "U", "UMA": "A",
    "UMP": "U", "UMS": "U", "UN1": "X", "UN2": "X", "UNK": "X",
    "UR3": "U", "URD": "U", "US1": "U", "US2": "U", "US3": "T",
    "US5": "U", "USM": "U", "VAD": "V", "VAF": "V", "VAL": "V",
    "VB1": "K", "VDL": "X", "VLL": "X", "VLM": "X", "VMS": "X",
    "VOL": "X", "WCR": "GYG", "X  ": "G", "X2W": "E", "X4A": "N",
    "X9Q": "AFG", "XAD": "A", "XAE": "N", "XAL": "A", "XAR": "N",
    "XCL": "C", "XCN": "C", "XCP": "X", "XCR": "C", "XCS": "N",
    "XCT": "C", "XCY": "C", "XGA": "N", "XGL": "G", "XGR": "G",
    "XGU": "G", "XPR": "P", "XSN": "N", "XTH": "T", "XTL": "T",
    "XTR": "T", "XTS": "G", "XTY": "N", "XUA": "A", "XUG": "G",
    "XX1": "K", "XXY": "THG", "XYG": "DYG", "Y  ": "A", "YCM": "C",
    "YG ": "G", "YOF": "Y", "YRR": "N", "YYG": "G", "Z  ": "C",
    "Z01": "A", "ZAD": "A", "ZAL": "A", "ZBC": "C", "ZBU": "U",
    "ZCL": "F", "ZCY": "C", "ZDU": "U", "ZFB": "X", "ZGU": "G",
    "ZHP": "N", "ZTH": "T", "ZU0": "T", "ZZJ": "A", 'XXX': 'X',
    '---': 'X', '***': 'X', '-': 'X', 'PYL': 'O', 'HSD': 'X',
    'HSE': 'S', 'HSP': 'X', 'XAA': 'X', 'MPR': 'X', 'FRD': 'X',
    'GLM': 'X', 'PPH': 'X', 'IVA': 'X', 'LOV': 'X', 'STA': 'X',
    'ETA': 'X', 'CYH': 'X',
}

aa_default_atoms = {
    'ALA': ['N', 'O', 'CA', 'C', 'CB'],
    'CYS': ['N', 'O', 'CA', 'C', 'SG', 'CB'],
    'ASP': ['N', 'O', 'CA', 'C', 'CG', 'OD2', 'OD1', 'CB'],
    'GLU': ['N', 'O', 'CA', 'C', 'CD', 'CG', 'OE1', 'CB', 'OE2'],
    'PHE': ['N', 'O', 'CA', 'C', 'CE1', 'CD1', 'CZ', 'CG', 'CB', 'CE2', 'CD2'],
    'GLY': ['N', 'O', 'CA', 'C'],
    'HIS': ['N', 'O', 'CA', 'C', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'LYS': ['N', 'O', 'CA', 'C', 'NZ', 'CD', 'CE', 'CG', 'CB'],
    'ILE': ['N', 'O', 'CA', 'C', 'CG2', 'CD1', 'CB', 'CG1'],
    'LEU': ['N', 'O', 'CA', 'C', 'CD2', 'CD1', 'CG', 'CB'],
    'MET': ['N', 'O', 'CA', 'C', 'SD', 'CG', 'CB', 'CE'],
    'ASN': ['N', 'O', 'CA', 'C', 'ND2', 'OD1', 'CB', 'CG'],
    'PRO': ['N', 'O', 'CA', 'C', 'CD', 'CG', 'CB'],
    'GLN': ['N', 'O', 'CA', 'C', 'CD', 'NE2', 'OE1', 'CB', 'CG'],
    'ARG': ['N', 'O', 'CA', 'C', 'CZ', 'CD', 'NE', 'CG', 'CB'],
    'SER': ['N', 'O', 'CA', 'C', 'CB', 'OG'],
    'THR': ['N', 'O', 'CA', 'C', 'OG1', 'CG2', 'CB'],
    'VAL': ['N', 'O', 'CA', 'C', 'CB', 'CG2', 'CG1'],
    'TYR': ['N', 'O', 'CA', 'C', 'CE1', 'CD1', 'CZ', 'CG', 'CB', 'CE2', 'CD2'],
    'TRP': ['N', 'O', 'CA', 'C', 'NE1', 'CZ3', 'CD1', 'CE3', 'CG', 'CB', 'CZ2', 'CE2', 'CD2'],
}

# obtained from biopython at
# https://github.com/biopython/biopython/blob/master/Bio/PDB/DSSP.py
# Miller max acc: Miller et al. 1987 http://dx.doi.org/10.1016/0022-2836(87)90038-6
# Wilke: Tien et al. 2013 http://dx.doi.org/10.1371/journal.pone.0080635
# Sander: Sander & Rost 1994 http://dx.doi.org/10.1002/prot.340200303
ASA_Miller = {
    'ALA': 113.0, 'ARG': 241.0, 'ASN': 158.0, 'ASP': 151.0,
    'CYS': 140.0, 'GLN': 189.0, 'GLU': 183.0, 'GLY': 85.0,
    'HIS': 194.0, 'ILE': 182.0, 'LEU': 180.0, 'LYS': 211.0,
    'MET': 204.0, 'PHE': 218.0, 'PRO': 143.0, 'SER': 122.0,
    'THR': 146.0, 'TRP': 259.0, 'TYR': 229.0, 'VAL': 160.0
}

ASA_Wilke = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0,
    'CYS': 167.0, 'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0,
    'HIS': 224.0, 'ILE': 197.0, 'LEU': 201.0, 'LYS': 236.0,
    'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 'SER': 155.0,
    'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
}

ASA_Sander = {
    'ALA': 106.0, 'ARG': 248.0, 'ASN': 157.0, 'ASP': 163.0,
    'CYS': 135.0, 'GLN': 198.0, 'GLU': 194.0, 'GLY': 84.0,
    'HIS': 184.0, 'ILE': 169.0, 'LEU': 164.0, 'LYS': 205.0,
    'MET': 188.0, 'PHE': 197.0, 'PRO': 136.0, 'SER': 130.0,
    'THR': 142.0, 'TRP': 227.0, 'TYR': 222.0, 'VAL': 142.0
}

nucleotides_common = [
    'A', 'C', 'G', 'I', 'U', 'DA', 'DC', 'DG', 'DI', 'DT', 'DU', 'N'
]

metals_common = [
    'LI', 'BE', 'NA', 'MG', 'AL', 'K', 'CA', 'SC', 'TI', 'V', 'CR', 'MN',
    'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO',
    'TC', 'RU', 'RH', 'PD', 'AG', 'CD', 'IN', 'SN', 'CS', 'BA', 'LA', 'CE',
    'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM', 'YB',
    'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL', 'PB',
    'BI', 'PO', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP', 'PU', 'AM', 'CM',
    'BK', 'CF',
     # confirm the following ones
    'LR', 'ES', 'NO', 'GE', 'RF', 'DB', 'SB', 'FM', 'MD', 'SG'
]

# Protein sequence alignments: a strategy for the hierarchical
# analysis of residue conservation. Livingstone & Barton (1993).
# hydrophobic,polar,small,proline,tiny,aliphatic,aromatic,positive,negative,charged
# Among the most problematic types of substitutions are those that involve:
# disruption of hydrogen-bond networks;
# disruption of disulphide contacts;
# change in folding/conformation dynamics;
# insertion of steric clashes;
# insertion of void space;
# change in ionization state by mutating to opposite charges;
# change in hydrophocity;
# mutate from or to Pro (also cis-Pro);
# mutate from Gly to any other amino acid;

aa_physicochemical = {
    'A': ['hydrophobic', 'small', 'tiny'],
    'C': ['hydrophobic', 'small'],
    'E': ['polar', 'negative', 'charged'],
    'D': ['polar', 'small', 'negative', 'charged'],
    'G': ['hydrophobic', 'small', 'tiny'],
    'F': ['hydrophobic', 'aromatic'],
    'I': ['hydrophobic', 'aliphatic'],
    'H': ['hydrophobic', 'polar', 'aromatic', 'positive', 'charged'],
    'K': ['hydrophobic', 'polar', 'aromatic', 'positive', 'charged'],
    'M': ['hydrophobic'],
    'L': ['hydrophobic', 'aliphatic'],
    'N': ['polar', 'small'],
    'Q': ['polar'],
    'P': ['small', 'proline'],
    'S': ['polar', 'small', 'tiny'],
    'R': ['polar', 'positive', 'charged'],
    'T': ['hydrophobic', 'polar', 'small'],
    'W': ['hydrophobic', 'polar', 'aromatic'],
    'V': ['hydrophobic', 'small', 'aliphatic'],
    'Y': ['hydrophobic', 'polar']
}

# David, A. & Sternberg, M. J. The contribution of missense mutations in core and Rim residues
# of protein-protein interfaces to human disease. J. Mol. Biol. 1–34 (2015).
# doi:10.1016/j.jmb.2015.07.004
aa_polar = ['Q', 'N', 'H', 'S', 'T', 'Y', 'C', 'M', 'W']
aa_hydrophobic = ['A', 'I', 'L', 'F', 'V']
aa_proline = ['P']
aa_glycine = ['G']
aa_pos_charge = ['R', 'K']
aa_neg_charge = ['D', 'E']
aa_physicochemical_minimal = {
    'N': 'polar',
    'Q': 'polar',
    'S': 'polar',
    'T': 'polar',
    'W': 'polar',
    'Y': 'polar',
    'H': 'polar',
    'M': 'polar',
    'C': 'polar',
    'A': 'hydrophobic',
    'F': 'hydrophobic',
    'I': 'hydrophobic',
    'L': 'hydrophobic',
    'V': 'hydrophobic',
    'R': 'positive_charged',
    'K': 'positive_charged',
    'E': 'negative_charged',
    'D': 'negative_charged',
    'G': 'glycine',
    'P': 'proline'
}

sequence_similarity = {
    'A': ['A', 'G', 'V', 'L', 'I'],
    'C': ['C', 'M'],
    'E': ['D', 'E', 'N', 'Q'],
    'D': ['D', 'E', 'N', 'Q'],
    'G': ['A', 'G', 'V', 'L', 'I'],
    'F': ['F', 'Y', 'W'],
    'I': ['A', 'G', 'V', 'L', 'I'],
    'H': ['K', 'R', 'H'],
    'K': ['K', 'R', 'H'],
    'M': ['C', 'M'],
    'L': ['A', 'G', 'V', 'L', 'I'],
    'N': ['D', 'E', 'N', 'Q'],
    'Q': ['D', 'E', 'N', 'Q'],
    'P': ['P'],
    'S': ['S', 'T'],
    'R': ['K', 'R', 'H'],
    'T': ['S', 'T'],
    'W': ['F', 'Y', 'W'],
    'V': ['A', 'G', 'V', 'L', 'I'],
    'Y': ['F', 'Y', 'W']
}
