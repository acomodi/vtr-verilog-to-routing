#!/usr/bin/env python3
"""
    Script to download the SymbiFlow Series-7 architectures
"""

import sys
import os
import argparse
import math
import textwrap
import fnmatch
import tempfile
import shutil
import subprocess
from urllib import request

GCS_URL = (
    "https://storage.googleapis.com/symbiflow-arch-defs-gha/latest?authuser=0"
)

SYMBIFLOW_URL_MIRRORS = {"google": GCS_URL}


class ExtractionError(Exception):
    """
    Extraction error exception class
    """


def parse_args():
    """
    Parses and returns script's arguments
    """

    description = textwrap.dedent(
        """
            Download and extract a symbiflow benchmark release into a
            VTR-style directory structure.

            If a previous matching symbiflow release tar.gz file is found
            does nothing (unless --force is specified).
        """
    )
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description=description
    )

    parser.add_argument(
        "--vtr_flow_dir",
        required=True,
        help="The 'vtr_flow' directory under the VTR tree. "
        "If specified this will extract the symbiflow release, "
        "placing benchmarks under vtr_flow/benchmarks/symbiflow ",
    )
    parser.add_argument(
        "--force",
        default=False,
        action="store_true",
        help="Run extraction step even if directores etc. already exist",
    )

    parser.add_argument("--mirror", default="google", choices=["google"], help="Download mirror")

    parser.add_argument(
        "--upgrade_archs",
        action="store_true",
        default=True,
        help="Try to upgrade included architecture files (using the upgrade_archs.py)",
    )

    return parser.parse_args()


def main():
    """
    Main function
    """

    args = parse_args()

    try:
        tar_xz_url = SYMBIFLOW_URL_MIRRORS[args.mirror]

        tar_xz_filename = "symbiflow.tar.xz"

        print("Downloading {}".format(tar_xz_url))
        download_url(tar_xz_filename, tar_xz_url)

        print("Extracting {}".format(tar_xz_filename))
        extract_to_vtr_flow_dir(args, tar_xz_filename)

    except ExtractionError as error:
        print("Failed to extract data: ", error)
        sys.exit(1)

    sys.exit(0)


def download_url(filename, url):
    """
    Downloads the symbiflow release
    """
    latest_package_url = request.urlopen(url).read().decode("utf-8")
    print("Downloading latest package:\n{}".format(latest_package_url))
    request.urlretrieve(latest_package_url, filename, reporthook=download_progress_callback)


def download_progress_callback(block_num, block_size, expected_size):
    """
    Callback for urllib.urlretrieve which prints a dot for every percent of a file downloaded
    """
    total_blocks = int(math.ceil(expected_size / block_size))
    progress_increment = int(math.ceil(total_blocks / 100))

    if block_num % progress_increment == 0:
        sys.stdout.write(".")
        sys.stdout.flush()
    if block_num * block_size >= expected_size:
        print("")


def extract_to_vtr_flow_dir(args, tar_xz_filename):
    """
    Extracts the 'benchmarks' directory of the symbiflow release
    into its corresponding vtr directory
    """

    # Reference directories
    arch_dir = os.path.join(args.vtr_flow_dir, "arch")
    symbiflow_arch_extract_dir = os.path.join(arch_dir, "symbiflow")

    arch_upgrade_script = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "upgrade_arch.py"
    )

    if not args.force:
        # Check that all expected directories exist
        expected_dirs = [
            args.vtr_flow_dir,
            symbiflow_arch_extract_dir,
        ]
        for directory in expected_dirs:
            if not os.path.isdir(directory):
                raise ExtractionError("{} should be a directory".format(directory))

    # Create a temporary working directory
    tmpdir = tempfile.mkdtemp(suffix="download_symbiflow", dir=".")

    # Extract matching files into the temporary directory
    subprocess.call(
        "tar -C {} -xf {} share/symbiflow/arch/xc7a50t_test".format(tmpdir, tar_xz_filename),
        shell=True,
    )

    # Move the extracted files to the relevant directories, SDC files first (since we
    # need to look up the BLIF name to make it match)
    for dirpath, _, filenames in os.walk(tmpdir):
        for filename in filenames:
            src_file_path = os.path.join(dirpath, filename)
            dst_file_path = None

            if fnmatch.fnmatch(src_file_path, "*/xc7a50t_test/arch.timing.xml"):
                dst_file_path = os.path.join(symbiflow_arch_extract_dir, "arch.timing.xml")

                subprocess.call("{} {}".format(arch_upgrade_script, src_file_path), shell=True)

            elif fnmatch.fnmatch(src_file_path, "*/xc7a50t_test/*.bin"):
                dst_file_path = os.path.join(symbiflow_arch_extract_dir, filename)

            if dst_file_path:
                shutil.move(src_file_path, dst_file_path)

    shutil.rmtree(tmpdir)

    print("Done")


if __name__ == "__main__":
    main()
