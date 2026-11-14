"""Check an installed wheel from outside the source tree, without extra packages."""

import importlib.metadata
import json
import platform

import nucleide
from nucleide import alara, r2s

distribution_version = importlib.metadata.version("nucleide")
assert nucleide.__version__ == distribution_version

regrouped = alara.regroup_flux([0, 2, 5], [4, 9], [0, 1, 3, 5], model="flat_energy")
assert regrouped["flux"] == [2.0, 5.0, 6.0]
assert regrouped["source_total"] == regrouped["target_total"] == 13.0

tags = r2s.tag_zone_strength(
    [120.0, 40.0], [0, 0, 1], split=True, voxel_volumes_cm3=[1.0, 3.0, 2.0]
)
assert tags["source_strength"] == [30.0, 90.0, 40.0]
assert tags["total"] == 160.0

print(
    json.dumps(
        {
            "version": distribution_version,
            "python": platform.python_version(),
            "libc": platform.libc_ver(),
            "package": nucleide.__file__,
            "regrouping_and_volume_split": "passed",
        }
    )
)
