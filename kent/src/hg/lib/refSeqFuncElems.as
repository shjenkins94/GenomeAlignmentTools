table refSeqFuncElems
"Bed 9+ representation of RefSeq functional elements."
    (
    string chrom;          "Reference sequence chromosome or scaffold"
    uint   chromStart;     "Start position in chromosome"
    uint   chromEnd;       "End position in chromosome"
    string name;           "type of element"
    uint   score;          "unused; placeholder for BED format"
    char[1] strand;        "+ for forward strand, - for reverse strand"
    uint   thickStart;     "Start position in chromosome"
    uint   thickEnd;       "End position in chromosome"
    uint reserved;         "Used as itemRgb: color based on type of element"
    string soTerm;         "Sequence ontology (SO) term"
    lstring note;          "A note describing the element"
    lstring geneIds;       "Entrez Gene ID of element and/or associated gene(s)"
    lstring pubMedIds;     "PubMed ID of associated publication(s)"
    lstring experiment;    "Experimental evidence"
    lstring function;      "Predicted function"
    lstring _mouseOver;    "Mouse over label"
    )