"""
/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the New BSD License, (the "License").
******************************************************************************/
"""

import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.core.surface import SlabGenerator
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def cjson_to_structure(cjson):
    numbers = cjson["atoms"]["elements"]["number"]
    uc = cjson["unitCell"]
    lattice = Lattice.from_parameters(
        uc["a"], uc["b"], uc["c"],
        uc["alpha"], uc["beta"], uc["gamma"]
    )

    frac = np.array(cjson["atoms"]["coords"]["3dFractional"]).reshape(-1, 3)
    # Filter out periodic image atoms (fractional coord ~1.0 == duplicate of 0.0)
    keep = ~np.any(np.abs(frac - 1.0) < 1e-4, axis=1)
    species = [numbers[i] for i, k in enumerate(keep) if k]
    return Structure(lattice, species, frac[keep], coords_are_cartesian=False)


def slab_to_cjson(slab):
    lat = slab.lattice
    cart = slab.cart_coords
    frac = slab.frac_coords
    numbers = [site.specie.Z for site in slab]
    return {
        "chemicalJson": 1,
        "atoms": {
            "elements": {"number": numbers},
            "coords": {
                "3d": cart.flatten().tolist(),
                "3dFractional": frac.flatten().tolist(),
            },
        },
        "unitCell": {
            "a": lat.a, "b": lat.b, "c": lat.c,
            "alpha": lat.alpha, "beta": lat.beta, "gamma": lat.gamma,
            "cellVectors": lat.matrix.flatten().tolist(),
        },
    }


def run(avo_input):
    cjson = avo_input.get("cjson", {})
    opts = avo_input.get("options", {})

    if "unitCell" not in cjson:
        return {"message": "No unit cell found. Please open a periodic structure first."}

    structure = cjson_to_structure(cjson)

    # SlabGenerator expects the conventional standard cell for correct Miller indices
    conv = SpacegroupAnalyzer(structure).get_conventional_standard_structure()

    h = int(opts.get("h", 1))
    k = int(opts.get("k", 1))
    l = int(opts.get("l", 1))
    layers = int(opts.get("layers", 4))
    vacuum = float(opts.get("vacuum", 10.0))
    term_idx = int(opts.get("termination", 0))

    gen = SlabGenerator(conv, (h, k, l), layers, vacuum,
                        primitive=False, center_slab=True)
    slabs = gen.get_slabs()

    if not slabs:
        return {"message": f"No slabs generated for Miller index ({h} {k} {l})."}

    term_idx = min(term_idx, len(slabs) - 1)
    slab = slabs[term_idx]

    return {
        "cjson": slab_to_cjson(slab),
        "message": f"Termination {term_idx + 1} of {len(slabs)} "
                   f"for ({h} {k} {l}) surface.",
    }
