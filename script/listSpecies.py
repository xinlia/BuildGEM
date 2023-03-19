import libsbml

r = libsbml.SBMLReader()
doc = r.readSBML(
    ".\\Metabolic_model_Methylocystis_bryophila.xml"
)
m = doc.getModel()

for n in range(0, m.getNumSpecies()):
    s = m.getSpecies(n)
    # print Species id
    print(
        "Species ",
        s.getId(),
        ": ",
    )
    print(
        s.getName(),
        " ",
    )
    #print (s.getChemicalFormula() , " ",)
    #fbc:charge="-8" fbc:chemicalFormula
    print("\n\n")
