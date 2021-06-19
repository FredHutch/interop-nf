import os
from pathlib import Path

tests_dir = os.path.dirname(__file__)


def get_resource_path(path: str):
    return str(Path(os.path.join(tests_dir, path)))


miseq_demo_path = get_resource_path('data/MiSeq/')
