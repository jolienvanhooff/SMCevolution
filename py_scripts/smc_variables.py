#!/usr/bin/env python3

# Module used to import various variables - such as dictionaries - often used in SMCevolution project

complex_codes = {
    "Cond": "condensin",
    "CondI":"condensin I",
    "CondII":"condensin II",
    "SMC56": "SMC5/6",
    "Coh": "cohesin"
}

complex_members = {
    "condensin":["SMC2","SMC4"],
    "condensin I":["CAPH","CAPG","CAPD2"],
    "condensin II":["CAPH2", "CAPG2", "CAPD3"],
    "SMC5/6":["SMC5", "SMC6", "Nse4", "Nse1", "Nse3", "Nse2", "Nse5", "Nse6"],
    "cohesin":["SMC1", "SMC3", "Scc1", "Rec8", "Scc3", "PDS5", "NIPBL", "MAU2", "Separase", "WAPL", "Eco1", "Securin", "Sororin", "CTCF", "DYAD", "Haspin", "Shugoshin"]
}

# 0: fill, 1: stroke + text
complex_colours = {
    "condensin":("#FF8233","#710000"),
    "condensin I":("#FFC300","#AA7A38"),
    "condensin II":("#C70039","#C70039"),
    "cohesin":("#3362FF","#155289"),
    "SMC5/6":("#A633FF","#883268")
}