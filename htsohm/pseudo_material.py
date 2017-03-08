class PseudoMaterial:

    def __init__(self, uuid):
        self.uuid = uuid

    def __repr__(self):
        return ('%s(' % self.__class__.__name__ +
            'uuid=%r, ' % self.uuid +
            'lattice_constants=%r, ' % self.lattice_constants +
            'atom_types=%r, ' % self.atom_types +
            'atom_sites=%r)' % self.atom_sites
        )
