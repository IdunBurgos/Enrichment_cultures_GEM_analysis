def translation_dicts():
    
    compounds_dict={'Acetic acid  (mg/L)':"ac", 
      'Ethanol (mg/L)':"etoh", 
      'Propionic acid (mg/L)':"ppa",
      'Butyric acid (mg/L)':"but",
      'Lactic acid (mg/L)':"lac__L",
      'Formic acid (mg/L)':"for",
      'Valeric acid (mg/L)':"pta",
      '1-Propanol (mg/L)':"ppoh",
      'Isovaleric acid (mg/L)':"ival",
      'Caproic acid (mg/L)':"hxa", 
      'Isobutyric acid (mg/L)':"ibt",
      'Levulinic acid':"4oxptn"}
    
    compounds_dict = {key.replace(" ",""):value for key,value in compounds_dict.items()}

    

    source_dict = {"Marshland soil":"M",
    "Compost and Digestate":"CD",
    "Cow manure":"CM"}

    substrate_dict = {"Xylan":"X",
    "Avicel":"A",
    "PASC":"P"}


    gas_sheet_dict = {"Marshland soil":"Marshland soil gas",
    "Compost and Digestate":"Compost and Digestate Gas",
    "Cow manure":"Cow manure Gas"}

    community_dict = {('Compost_Digestate', 'Avicel'):"CD_A",
                       ('Compost_Digestate', 'PASC'):"CD_P",
                       ('Compost_Digestate', 'Xylan'):"CD_X",
                       ('Cow_Manure', 'Avicel'):"CM_A",
                       ('Cow_Manure', 'PASC'):"CM_P", 
                       ('Cow_Manure', 'Xylan'):"CM_X", 
                       ('Marshland', 'PASC'):"M_P", 
                       ('Marshland', 'Xylan'):"M_X"}
    
    
    return compounds_dict, source_dict,substrate_dict, gas_sheet_dict, community_dict
