# In this file, lowercase s at the start of a string denotes a SMARTS string
# (i.e., a full SMARTS pattern or just a fragment of a SMARTS pattern)

# ####################################################
# Common SMARTS patterns.
# ####################################################
s_carboxyl = "[CX3](=O)[OX2H1]"
s_carbonyl = "[$([CX3]=[OX1]),$([CX3+]-[OX1-])]"
s_halogen = "[F,Cl,Br,I]"
s_halide = f"*{s_halogen}"
s_acyl_halide = f"[CX3](=[OX1]){s_halogen}"
# ####################################################

# SMARTS string with ID 1. It can be used to match atoms without certain adjacent ring
# groups defined in the string.
s_1 = f"&!$([*R]*{s_carboxyl})&!$([*R]{s_carbonyl})&!$([*R]{s_halide})&!$([*R]*{s_acyl_halide})"

# The primary SMARTS string used for ROMP reactions.
s_romp_rxn_smarts = f"[CH1R{s_1}:0]=[CH1R{s_1}:1]>>([#0]=[C:0].[#0]=[C:1])"
romp_rxn_smarts = s_romp_rxn_smarts  # alias
