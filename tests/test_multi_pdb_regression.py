import pytest

from compute_distance import get_biopython_structure


@pytest.mark.parametrize(
    "pdb_path",
    [
        "pdb1a1m.ent",
        "pdb1fat.ent",
        "esmfold_predicted.pdb",
    ],
)
def test_pdb_can_be_parsed(pdb_path):
    structure = get_biopython_structure("test", pdb_path)
    assert structure is not None