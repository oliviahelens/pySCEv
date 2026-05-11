"""
Curated mouse gene-set panel for pancreatic endocrinogenesis.

Hand-picked from well-established markers in the pancreas / islet-biology
literature, used as a fallback when external pathway libraries
(MSigDB Hallmark via Enrichr, Reactome, PROGENy, CollecTRI) are unavailable.

Symbols are UPPERCASE to match the rank-metric uppercasing in
validation/pancreas_endocrinogenesis_enrichment/run_enrichment.py:_signed_rank.
That's the same human-MSigDB convention; preserved here so the same script
can be re-run against MSigDB without re-mangling symbols.

Sets are scoped to pancreas biology. They are NOT a substitute for a full
MSigDB Hallmark scan -- they bias toward signals we'd expect to see based on
the operational read of the metric (low entropy = coordinated motion -> early
/ proliferative; high entropy = scattered -> terminal / secreting). The whole
point is to test that hypothesis, so the panel includes confirmatory sets
(proliferation, hormone biosynthesis) and contrastive ones (UPR, hypoxia,
EMT) so the result is informative either way.
"""

from __future__ import annotations

PANCREAS_PANEL: dict[str, list[str]] = {
    # Confirmatory: expected UP in low pySCEv if hypothesis holds
    "PROLIFERATION_CELL_CYCLE": [
        "MKI67", "TOP2A", "CCNB1", "CCNB2", "CDK1", "CCNA2",
        "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "PCNA",
        "BIRC5", "AURKA", "AURKB", "CENPA", "CENPE", "CENPF",
        "KIF11", "KIF20A", "KIF23", "RRM1", "RRM2", "TYMS",
        "FOXM1", "TUBA1B", "TUBB5", "H2AFZ", "HMGB2", "STMN1",
    ],
    "ENDOCRINE_PROGENITOR_TFS": [
        "NEUROG3", "NEUROD1", "PAX4", "PAX6", "INSM1",
        "FOXA2", "FOXA3", "PDX1", "RFX6", "MNX1", "HNF1B",
        "SOX9", "NKX2-2", "NKX6-1",
    ],
    # Confirmatory: expected UP in high pySCEv if hypothesis holds
    "PANCREATIC_HORMONES": [
        "INS1", "INS2", "GCG", "SST", "PYY", "IAPP", "GHRL",
        "CHGA", "CHGB", "SCG2", "SCG3", "SCG5", "NPY", "PPY",
    ],
    "MATURE_BETA": [
        "INS1", "INS2", "IAPP", "GLP1R", "GIPR", "MAFA",
        "NKX6-1", "PDX1", "UCN3", "SLC2A2", "PCSK1", "PCSK2",
        "ERO1B", "G6PC2",
    ],
    "MATURE_ALPHA": [
        "GCG", "MAFB", "ARX", "IRX1", "IRX2", "TTR",
        "GC", "MEIS2", "POU3F4", "FEV",
    ],
    "SECRETORY_MACHINERY": [
        "SEC61A1", "SEC61B", "SEC61G", "SSR1", "SSR2", "SSR3", "SSR4",
        "SRP9", "SRP14", "SRP54", "SRP68", "SRP72",
        "ERLEC1", "EDEM1", "SEC23A", "SEC24A", "SEC24B",
        "SYT4", "SYT13", "SNAP25", "STX1A", "VAMP2",
    ],
    # Contrastive: would refine interpretation in either direction
    "UNFOLDED_PROTEIN_RESPONSE": [
        "ATF4", "ATF6", "ATF3", "DDIT3", "XBP1",
        "HSPA5", "HSP90B1", "ERN1", "EIF2AK3", "HERPUD1",
        "PDIA3", "PDIA4", "PDIA6", "CALR", "CANX",
    ],
    "HYPOXIA_GLYCOLYSIS": [
        "HIF1A", "VEGFA", "SLC2A1", "PGK1", "LDHA",
        "BNIP3", "PDK1", "HK1", "HK2", "ENO1", "ALDOA",
        "GAPDH", "PFKM", "TPI1", "PFKP",
    ],
    "OXPHOS": [
        "NDUFA1", "NDUFA2", "NDUFB1", "NDUFB2", "NDUFS1",
        "ATP5A1", "ATP5B", "ATP5C1", "ATP5F1", "ATP5O",
        "UQCR11", "UQCRC1", "UQCRC2", "COX4I1", "COX5A",
        "COX6A1", "COX7A1", "SDHA", "SDHB", "SDHC", "SDHD",
    ],
    "APOPTOSIS_INTRINSIC": [
        "BAX", "BAK1", "BCL2", "BCL2L1", "BID",
        "CASP3", "CASP6", "CASP7", "CASP8", "CASP9",
        "CYCS", "APAF1", "FAS", "FASLG", "TP53", "MDM2",
    ],
    "NOTCH_SIGNALING": [
        "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4",
        "HES1", "HEY1", "HEY2", "HEYL",
        "DLL1", "DLL3", "DLL4", "JAG1", "JAG2",
        "RBPJ", "MAML1", "MAML2",
    ],
    "WNT_SIGNALING": [
        "WNT3A", "WNT5A", "WNT5B", "WNT7A", "WNT7B",
        "CTNNB1", "AXIN1", "AXIN2", "LEF1", "TCF7", "TCF7L2",
        "FZD1", "FZD2", "FZD7", "DKK1", "DKK3",
        "APC", "GSK3B",
    ],
    "EMT_MIGRATION": [
        "VIM", "SNAI1", "SNAI2", "TWIST1", "TWIST2",
        "ZEB1", "ZEB2", "CDH2", "CDH11", "FN1",
        "MMP2", "MMP9", "MMP14", "ACTA2", "TAGLN", "S100A4",
    ],
    "DUCTAL_MARKERS": [
        "KRT19", "KRT7", "KRT8", "KRT18",
        "SOX9", "HNF1B", "CFTR", "SPP1", "MUC1", "MUC6",
        "ANXA2", "ANXA4", "ANXA5", "CLU",
    ],
    "EXOCRINE_ACINAR": [
        "CPA1", "CPA2", "PRSS1", "PRSS2", "PRSS3",
        "CTRB1", "CTRC", "CTRL", "AMY1", "AMY2A",
        "ELA1", "ELA2A", "ELA3A", "ELA3B", "RBPJL",
        "PNLIP", "PLA2G1B",
    ],
}
