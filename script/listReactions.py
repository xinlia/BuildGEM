# conda install <libsbml>
# pip install python-libsbml
import libsbml

r = libsbml.SBMLReader()
'''
doc = r.readSBML(
    ".\\5GB1_base_fermentation.xml"
)
'''
doc = r.readSBML(
    ".\\Metabolic_model_Methylocystis_bryophila.xml"
)

m = doc.getModel()

for n in range(0, m.getNumReactions()):
    r = m.getReaction(n)
    # print reaction id
    print(
        #"Reaction ",
        r.getId(),
        ": ",
        end=""
    )

    # look at reactants
    numReactants = r.getNumReactants()
    if (numReactants > 1):
        s = m.getSpecies(r.getReactant(0).getSpecies())
        print(
            r.getReactant(0).getStoichiometry(),
            " ",
            s.getName(),
            " ",
            end=""
        )
        for k in range(1, numReactants):
            # get species referred to by the reaction
            s = m.getSpecies(r.getReactant(k).getSpecies())
            print(
                "+ ",
                r.getReactant(k).getStoichiometry(),
                " ",
                s.getName(),
                " ",
                end=""
            )
    elif (numReactants == 1):
        # get species referred to by the reaction
        s = m.getSpecies(r.getReactant(0).getSpecies())
        print(
            r.getReactant(0).getStoichiometry(),
            " ",
            s.getName(),
            " ",
            end=""
        )

    if (r.getReversible() == True):
        print("<=> ", end="")
    else:
        print("=> ", end="")

# look at products
    numProducts = r.getNumProducts()
    if (numProducts > 1):
        s = m.getSpecies(r.getProduct(0).getSpecies())
        print(
            r.getProduct(0).getStoichiometry(),
            " ",
            s.getName(),
            " ",
            end=""
        )
        for k in range(1, numProducts):
            # get species referred to by the reaction
            s = m.getSpecies(r.getProduct(k).getSpecies())
            print(
                "+ ",
                r.getProduct(k).getStoichiometry(),
                " ",
                s.getName(),
                " ",
                end=""
            )
    elif (numProducts == 1):
        # get species referred to by the reaction
        s = m.getSpecies(r.getProduct(0).getSpecies())
        print(
            r.getProduct(0).getStoichiometry(),
            " ",
            s.getName(),
            " ",
            end=""
        )
    print("\n")
