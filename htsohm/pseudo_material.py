class PseudoMaterial:
    """Class for storing pseudomaterial structural data.

    Attributes:
        uuid (str) : Version 4 UUID used to identify pseudomaterial record in
            `materials` database.
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
                        "charge"       : float,
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
                        "z-frac"       : float
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
        self.lattice_constants = {
                "a" : None,
                "b" : None, 
                "c" : None
                }
        self.atom_types= [
                {
                    "chemical-id"  : None,
                    "charge"       : None,
                    "epsilon"      : None,
                    "sigma"        : None
                    }
                ]
        self.atom_sites = [
                {
                    "chemical-id"  : None,
                    "x-frac"       : None,
                    "y-frac"       : None,
                    "z-frac"       : None
                    }
                ]

    def __repr__(self):
        return ('{0.__class__.__name__!s}('
                '{0.uuid!r}, '
                '{0.lattice_constants!r}, '
                '{0.atom_types!r}, '
                '{0.atom_sites!r})').format(self)
