"""Utilities for tests."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import os
from absl import flags


def test_dir(relative_path):
    """Gets the path to a test_dataset file in genomics at relative path.

    Args:
      relative_path: a directory path relative to base directory of this module.
    Returns:
      The absolute path to a test_dataset file.
    """

    return os.path.join(flags.FLAGS.test_srcdir,
                        os.path.split(os.path.abspath(__file__))[0],
                        relative_path)
