from LammpsUnifiedCleaner import file_unifier

DIR = '/home/matt/Documents/AddAtoms/'

file_unifier(DIR, 'system.in.settings', ['pre_reaction.data', 'post_reaction.data'])

if createAtoms is not None:
    # Extend renumberedAtomDict with placeholder values for createAtoms
    # This get prevents key errors later when reducing the molecule file
    for createAtom in createAtoms:
        postRenumberedAtomDict[createAtom] = createAtom