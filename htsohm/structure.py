class Structure:
    """Class for storing pseudomaterial structural data.
    Attributes:
        uuid (str) : Version 4 UUID used to identify pseudomaterial record in
            `materials` database.
        run_id (str) : identification string distinguishing runs.
        lattice_constants (dict) : crystal lattice parameters:
                {
                    'a' : float,
                    'b' : float,
                    'c' : float
                }
        atom_types (list : dict) : Lennard-Jones parameters and partial charge:
                [
                    {
                        "chemical-id"  : str,
                        "epsilon"      : float,
                        "sigma"        : float
                    }
                ]
        atom_sites (list : dict) : atom-site locations as fractions:
                [
                    {
                        "chemical-id"  : str,
                        "x-frac"       : float,
                        "y-frac"       : float,
                        "z-frac"       : float,
                        "charge"       : float
                    }
                ]
    """

    def __init__(self, uuid):
        """Instantiates PseudoMaterial object with real uuid and null values for
        all other attributes.
        
        Args:
            uuid (str) : Version 4 UUID identifying pseudomaterial record in
                `materials` database.
        Returns:
            None
        """
        self.uuid = uuid
        self.run_id = None
        self.lattice_constants = {"a" : None, "b" : None, "c" : None}
        self.atom_types= []
        self.atom_sites = []

    def __repr__(self):
        return ('{0.__class__.__name__!s}('
                '{0.uuid!r}, '
                '{0.run_id!r}, '
                '{0.lattice_constants!r}, '
                '{0.atom_types!r}, '
                '{0.atom_sites!r})').format(self)

    def volume(self):
        return (
                self.lattice_constants["a"] *
                self.lattice_constants["b"] *
                self.lattice_constants["c"]
            )

    def number_density(self):
        return len(self.atom_sites) / self.volume()
