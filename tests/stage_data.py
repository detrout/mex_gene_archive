###
# Stage test data
import bz2
from contextlib import contextmanager
from pathlib import Path
import tarfile
from tempfile import TemporaryDirectory


test_dir = Path(__file__).parent
test_data_dir = test_dir / "testdata"

@contextmanager
def srx5908538():
    try:
        test_archive = test_data_dir / "SRX5908538.tar.bz2"
        analysis_dir = TemporaryDirectory(suffix="_srx5908538")
        with bz2.open(test_archive, "rb") as decompressed:
            tarball = tarfile.TarFile(fileobj=decompressed, mode="r")
            tarball.extractall(path=analysis_dir.name, numeric_owner=True)
            yield Path(analysis_dir.name)
    finally:
        analysis_dir.cleanup()


if __name__ == "__main__":
    import time
    import os
    t0 = time.monotonic()
    with srx5908538() as analysis_dir:
        tnow = time.monotonic()
        print("Time: {:.3f} seconds".format(tnow-t0))
        t0 = tnow
        print(type(analysis_dir), analysis_dir)
        print(os.listdir(analysis_dir))
