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
    "cohesin":["SMC1", "SMC3", "Scc1", "Rec8", "Scc3", "PDS5", "NIPBL", "MAU2", "WAPL", "Eco1", "Securin", "Sororin", "Haspin", "Shugoshin", "Separase", "CTCF"]
}

# 0: fill, 1: stroke + text
complex_colours = {
    "condensin":("#ddcc77","#AA7A38"),
    "condensin I":("#999933","#335A30"),
    "condensin II":("#117733","#024217"),
    "cohesin":("#4477aa","#155289"),
    "SMC5/6":("#aa3377","#7B255A")
}